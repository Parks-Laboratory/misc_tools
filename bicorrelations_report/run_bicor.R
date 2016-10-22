
# Markdown Report to Query Experimental Details from Lab Database

# Function:
  # 1) Correlations: intersects top genes that are highly correlated (bicorrelation) to gene of interest across genders and diet within a tissue set

# Parameters:
  # database: (string) name of the database to query data from
  # table_prefix: (string) the prefix of the tables to extract data from (refer to database); tables should have the column names gene_symbol1, gene_symbol2, bicor, pvalue
  # gene_name: (string) gene of interest
  # select_top: (integer) select top n genes that are highly correlated to gene of interest (ranked by pvalue); top genes from intersected for each tissue group
  
# Returns: PDF report

#######################################################################
## SET PARAMETERS #####################################################
#######################################################################

# set parameters
gene <- "sbno1"
database <- "external_labs" 
# current options for database are "external_labs" and "HMDP"
table_prefix <- "Chick_Munger_bicorrelations_all_rna"
# current options are 
# "Chick_Munger_bicorrelations_all_protein" or "Chick_Munger_bicorrelations_all_rna" for external labs database
# "bicorrelations_all" for HMDP database

#######################################################################

# create output folder if needed
setwd("E:/MICROARRAY DATA/bicorrelations")
if(!file.exists(file.path("outputs", stringr::str_trim(paste(database, table_prefix))))) dir.create(file.path("outputs", stringr::str_trim(paste(database, table_prefix))))

# run program
rmarkdown::render(
  # file to run the function
  input = "E:/MICROARRAY DATA/bicorrelations/intersect_bicorrelations.Rmd",
  
  # parameters to pass to function
  params = list(
    database = database, 
    table_prefix = table_prefix,
    gene_symbol = gene,
    select_top = 1000
  ),
  
  # location to save the output file
  output_file = file.path("outputs", stringr::str_trim(paste(database, table_prefix)), paste0(gene, "_output.pdf"))
)