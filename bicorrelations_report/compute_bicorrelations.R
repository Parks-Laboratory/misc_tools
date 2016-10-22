#------------------------------------------------------------------------
# Title: Compute Gene-Gene Bicorrelations (Expression) 
# Date: May 17 2016
# Author: Jenny Nguyen
# Email: jnnguyen2@wisc.edu
#------------------------------------------------------------------------

# CONTENTS:
# Part 1: Download Data from Database
# Part 2: Run Bicor Computations on Condor
# Part 3: Export Data as CSV to Upload to SQL


# NOTES:
  # More than 1 probeset id for gene_symbol
  # Need to find correlation by probeid and then switch to gene_symbol
  # IE may have same gene_symbol correlated to itself
#------------------------------------------------------------------------


# Input Command Line Args
#------------------------------------------------------------------------

library(optparse)

# set arguments
option_list <- list(
  make_option(c("--condor"), default = FALSE, help = "Whether code is running on Condor [default %default]"), 
  make_option(c("--sql"), default = FALSE, help = "Whether code is running to upload to SQL [default %default]"), 
  make_option(c("--option"), default = "bicor", help = "Whether extract bicor or p"),
  make_option(c("--i"), default = 0, help = "What element of data to run (1-4)")
)

# process arguments
opt <- parse_args(OptionParser(option_list = option_list))

# set input parameters
condor <- opt$condor
sql <- opt$sql
option <- opt$option
i <- as.numeric(opt$i) + 1


# Load Libraries
#------------------------------------------------------------------------
library(reshape2)
library(stringr)
library(dplyr)
library(data.table)
library(RODBC)
library(WGCNA)


# Part 1: Download Data from Database
#------------------------------------------------------------------------

if(!condor){
  
  setwd("E:/MICROARRAY DATA")
  
  # connect to SQL Server
  db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={SQL Server}')
  
  # obtain the names of table/views from SQL Server
  data_sources <- sqlTables(db) %>% 
    subset(TABLE_SCHEM == "dbo") %>% 
    dplyr::select(TABLE_CAT, TABLE_SCHEM, TABLE_NAME) %>% 
    subset(str_detect(TABLE_NAME, "avg_expression|expression_differences"))
  
}


if(!condor & !sql){
  
  # Pull Data from Database
  #------------------------------------------------------------------------
  
  # generate queries to extract from SQL server
  queries <- paste("select * from", paste(data_sources$TABLE_SCHEM, data_sources$TABLE_NAME, sep = "."))
  
  # extract data from SQL server
  data <- lapply(queries, function(q) sqlQuery(db, q))
  
  # close the database
  odbcClose(db)
  
  
  # Save data to run on HTCondor
  #------------------------------------------------------------------------
  
  save(data, file = "bicor_data.Rdata")
  
}


# Part 2: Run Bicor Computations on Condor
#------------------------------------------------------------------------

if(condor){
  
  # Load data to run on HTCondor
  #------------------------------------------------------------------------
  
  load("bicor_data.Rdata")
  data <- data.table(data[[i]])
  
  
  # Compute Bicorrelations
  #------------------------------------------------------------------------
  
  # convert probes wide and remove strain info
  b <- data %>% 
    dcast(strain ~ probe_id, value.var = "mean_rma_expression") %>% 
    dplyr::select(-strain)
  
  # compute bicorrelations between probes
  b <- bicorAndPvalue(b, alternative = "two.sided")
  
  # extract bicorrelations and p-values
  # convert from cor matrix format to data.frame format
  b <- b[[option]] %>% 
    melt() %>% 
    subset(Var1 != Var2)
  colnames(b) <- c("probeid1", "probeid2", "value")
  
  # save data
  bicor_data <- list(b)
  names(bicor_data) <- option
  
  # Save data to upload to SQL server
  #------------------------------------------------------------------------
  
  save(bicor_data, file = paste0("bicor_data_out_", option, "_", i, ".Rdata"))
  
}


# Part 3: Export Data as CSV to Upload to SQL
#------------------------------------------------------------------------

if(sql){
  
  # Export Bicorrelations to SQL server
  #------------------------------------------------------------------------
  
  # close the database
  odbcClose(db)
  
  # generate table_names
  table_names <- str_replace(data_sources$TABLE_NAME, "_avg_expression_by_strain", "")
  
  
  # this doesn't even really work - have to run this one by one
  # size of files too big for disk space
  Sys.time()
  print("Starting the Loop")
  
  for(j in 1:8){
  
    print(paste0("iteration = ", j))
    
    # save the correlations
    print( paste0("Loading bicor data ", Sys.time()) )
    load( url(paste0("http://pages.stat.wisc.edu/~jnguyen/", "bicor_data_out_bicor_", j, ".Rdata")) )
    print(dim(bicor_data[[1]]))

    data.table::fwrite(bicor_data[[1]], file.path = paste0("E:/MICROARRAY DATA/bicor_data", "_out_data_bicor_", table_names[j], ".csv"))
    rm(bicor_data)

    print( Sys.time() )
    
    # save the p-value
    print( paste0("Loading pvalue data ", Sys.time()) )
    load( url(paste0("http://pages.stat.wisc.edu/~jnguyen/", "bicor_data_out_p_", j, ".Rdata")) )
    print(dim(bicor_data[[1]]))

    data.table::fwrite(bicor_data[[1]], file.path = paste0("E:/MICROARRAY DATA/bicor_data", "_out_data_bicor_pval_", table_names[j], ".csv"))
    rm(bicor_data)
    print( Sys.time() )
  }
  
  print("Exiting the Loop")
  Sys.time()
  
  
}
