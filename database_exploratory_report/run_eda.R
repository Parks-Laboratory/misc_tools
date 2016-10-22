
# Markdown Report to Query Experimental Details from Lab Database

# Function:
# 1) Correlations: intersects top genes that are highly correlated (bicorrelation) to gene of interest across genders and diet within a tissue set; also intersects across mouse and human if adipose is specified
# 2) Transcript Abundance: plots the RNA Expression of gene of interest in different tissues for different diets and genders
# 3) eQTL: obtains SNPS that are highly correlated (pvalue < 1e-05) with gene of interest
      # Plots SNP location by chromosome
      # Prints name of and location of select number of SNPs 

# Parameters:
  # gene_name: (string) gene of interest
  # select_top: (integer) select top n genes that are highly correlated to gene of interest (ranked by pvalue); top genes from intersected for each tissue group
  # tissues: (vector) tissues of interest
  # opt_output_all_cor_tables: (boolean) whether to output csv's of genes that are highly correlated with gene of interest (ordered by pvalue) in all tissues

# Returns: PDF report


# gene to query
gene <- "ypel5"

rmarkdown::render(
  # file to run the function
  input = "M:/R/database_to_pdf/db_to_pdf_output.Rmd",
  
  # parameters to pass to function
  params = list(
    gene_name = gene,
    tissues = c("Adipose", "Liver"),
    select_top = 500,
    opt_output_all_cor_tables = FALSE
  ),
  
  # location to save the output file
  output_file = paste0("M:/R/database_to_pdf/output/", gene, "_output.pdf")
)
