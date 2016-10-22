#------------------------------------------------------------------------
# Title: Wrapper Function to Plot Expressions of Given Gene by Strain and Diet
#------------------------------------------------------------------------

# Note:
  # all tables should have a column named "gene_symbol"
  # all tables should have a column with the word "expression" somewhere in the name; there should only be one expression column
  # if strain information is provided, tables should have a column with the word "strain" somewhere in the name

# Parameters: 
  # gene_symbol, the name of the gene to query
  # database, the name of the database to query from
  # table_prefix, the prefix of the tables to extract data from (refer to database); tables should all have the same column names
  # use_tissue, tissues to plot (adipose, liver, etc)
  # plot_boxplots, whether to do boxplots; if TRUE makes boxplots, if FALSE makes scatter plots
  # group_by_strain, option to place strain in x-axis (only relevant if plot_boxplots is FALSE); otherwise orders every diet by expression
# Returns: list of expression plots


#######################################################################
## SET PARAMETERS #####################################################
#######################################################################

# set parameters
gene <- "sun1"
database <- "HMDP" 
# current options for database are "external_labs" and "HMDP"
table_prefix <- "avg_expression_by_strain"
# current options are 
  # "Chick_Munger_full_protein_expression" or "Chick_Munger_full_rna_expression" for external labs database
  # "avg_expression_by_strain" for HMDP database
tissues <- c("liver")  #adipose
group_by_strain <- FALSE
plot_boxplots <- FALSE

#######################################################################

# get plot data
setwd("E:/MICROARRAY DATA/expression_plots/")
source('E:/MICROARRAY DATA/expression_plots/expression_plots.R')

p <- expression_plots(
  gene_symbol = gene, 
  database = database, 
  table_prefix = table_prefix,
  use_tissue = tissues,
  plot_boxplots <- plot_boxplots,
  group_by_strain <- group_by_strain 
)

# print pdf
opt_title <- ifelse(plot_boxplots, "boxplots_", "")
if(!file.exists(file.path("outputs", table_prefix))) dir.create(file.path("outputs", table_prefix))

pdf( file.path("outputs", table_prefix, paste0("transcript_abundance_", opt_title, gene, ".pdf") ), height = 6, width = 10)
plyr::l_ply(p, print)
dev.off()
