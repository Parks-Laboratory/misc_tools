#------------------------------------------------------------------------
# Title: Upload Clinical Traits into SQL Server
# Date: June 6 2016
# Author: Jenny Nguyen
# Email: jnnguyen2@wisc.edu
#------------------------------------------------------------------------

# CONTENTS:


# NOTES:
#------------------------------------------------------------------------

library(data.table)
library(stringr)
library(dplyr)
library(RODBC)


# Chow Clinical Traits
#------------------------------------------------------------------------

# load file
chow <- fread("file:///E:/MICROARRAY DATA/Clinical_Traits/clincal_traits_chow.csv", colClasses = "character")

# change names to lower case
setnames(chow, str_to_lower(colnames(chow)))

# replace all spaces with "_"
spaces <- colnames(chow) %>% str_subset("\\s")
replace <- str_replace_all(spaces, " ", "_")
setnames(chow, spaces, replace)

# change all "NULL" to NA
fix_null <- function(x) ifelse(x == "NULL", NA, x) 
chow <- mutate_each(chow, funs(fix_null))

# convert values to numeric
chow <- mutate_each(chow, funs(as.numeric), -number, -mouseid, -jackson_name, -sex)

# rename
setnames(chow, c("mouseid", "jackson_name"), c("mouse_id", "strain"))


# Highfat Clinical Traits
#------------------------------------------------------------------------

# load file
highfat <- fread("file:///E:/MICROARRAY DATA/Clinical_Traits/clinical_traits_highfat.csv", colClasses = "character")

# change names to lower case 
setnames(highfat, str_to_lower(colnames(highfat)))

# replace all spaces with "_"
spaces <- colnames(highfat) %>% str_subset("\\s")
replace <- str_replace_all(spaces, " ", "_")
setnames(highfat, spaces, replace)

# change all % 
perc <- colnames(highfat) %>% str_subset("%")
p_replace <- str_replace_all(perc, "%", "percent")
setnames(highfat, perc, p_replace)

# change all -
dash <- colnames(highfat) %>% str_subset("/")
ds_replace <- str_replace_all(dash, "/", "_")
setnames(highfat, dash, ds_replace)

# change all /
dash <- colnames(highfat) %>% str_subset("-")
ds_replace <- str_replace_all(dash, "-", "")
setnames(highfat, dash, ds_replace)

# change all dots at end of name to nothing
dots <- colnames(highfat) %>% str_subset(".*\\.")
d_replace <- str_replace(dots, "\\.$", "")
setnames(highfat, dots, d_replace)

# convert values to numeric
highfat <- mutate_each(highfat, funs(as.numeric), -mouse_id, -strain, -sex)

# fix for mouse id
highfat <- mutate(highfat, mouse_id = mouse_id %>% str_extract("\\d+") %>% str_pad(width = 4, side = "left", pad = "0"))


# Aggregate Data
#------------------------------------------------------------------------

fix_nan <- function(x) ifelse(is.nan(x), NA, x)

chow_agg <- chow %>% 
  group_by(strain) %>% 
  summarise_each(funs(mean(., na.rm = TRUE)), -number, -mouse_id, -sex) %>% 
  mutate_each(funs(fix_nan))

highfat_agg <- highfat %>% 
  group_by(strain, sex) %>% 
  summarise_each(funs(mean(., na.rm = TRUE)), -mouse_id) 

highfat_agg_m <- highfat_agg %>% 
  subset(sex == "Male") %>% 
  dplyr::select(-sex) %>% 
  mutate_each(funs(fix_nan))

highfat_agg_f <- highfat_agg %>% 
  subset(sex == "Female") %>% 
  dplyr::select(-sex) %>% 
  mutate_each(funs(fix_nan))


# Upload to SQL
#------------------------------------------------------------------------

db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={SQL Server}')
# sqlSave(channel = db, dat = chow, tablename = "clinical_traits_chow", rownames = FALSE)
# sqlSave(channel = db, dat = highfat, tablename = "clinical_traits_highfat", rownames = FALSE)
# sqlSave(channel = db, dat = chow_agg, tablename = "avg_clinical_traits_chow_male", rownames = FALSE)
# sqlSave(channel = db, dat = highfat_agg_m, tablename = "avg_clinical_traits_highfat_male", rownames = FALSE)
# sqlSave(channel = db, dat = highfat_agg_f, tablename = "avg_clinical_traits_highfat_female", rownames = FALSE)
odbcClose(db)


# Find Overlap Columns
#------------------------------------------------------------------------

# load data from data base (in case column names change in transition)
db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={SQL Server}')
chow <- sqlQuery(db, "select * from dbo.clinical_traits_chow")
highfat <- sqlQuery(db, "select * from dbo.clinical_traits_highfat")

# save column names to a csv file to find similarities
x <- colnames(chow)[5:ncol(chow)]
y <- colnames(highfat)[4:ncol(highfat)]
x <- c(x, rep("", ncol(highfat) - ncol(chow) + 1))
data.frame(chow = x, highfat = "", notes = "", highfat_original = y) %>% 
  write.csv("file:///E:/MICROARRAY DATA/Clinical_Traits/clinical_traits_in_common.csv", row.names = FALSE)

# load file of overlap clinical traits
traits <- fread("file:///E:/MICROARRAY DATA/Clinical_Traits/clinical_traits_in_common.csv")
traits <- traits %>% 
  dplyr::select(chow, highfat, common_name, notes) %>% 
  subset(chow != "")

# upload clinical traits to SQL Server
sqlSave(channel = db, dat = traits, tablename = "clinical_traits_in_common", rownames = FALSE)
odbcClose(db)


