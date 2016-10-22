#------------------------------------------------------------------------
# Title: Function to Plot Expressions of Given Gene by Strain and Diet
# Date: May 16 2016
# Author: Jenny Nguyen
# Email: jnnguyen2@wisc.edu
#------------------------------------------------------------------------

# CONTENTS:

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


# NOTES:
#
# Code to run function:
# source('E:/MICROARRAY DATA/expression_plots/expression_plots.R')
# expression_plots(gene_symbol = "agpat5")
#------------------------------------------------------------------------

# load libraries
library(reshape2)
library(plyr)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
library(RODBC)
library(ggplot2)
library(jn.general)
theme_set(theme_bw())

# Define function
expression_plots <- function(gene_symbol, 
                             database = "HMDP", table_prefix = 'avg_expression_by_strain', 
                             use_tissue = c("adipose", "liver"), 
                             plot_boxplots = FALSE, group_by_strain = TRUE
                             ){
  
  # connect to SQL Server
  db <- odbcDriverConnect( paste0('SERVER=PARKSLAB;DATABASE=', database, ';Trusted_Connection=Yes;DRIVER={SQL Server}') )
  
  # generate the search value for tissues
  tissue_search <- paste(use_tissue, collapse = "|")
  
  # obtain the names of table/views from SQL Server: tables that contain expression averaged by strain
  data_sources <- sqlTables(db) %>% 
    subset(TABLE_SCHEM == "dbo") %>% 
    dplyr::select(TABLE_CAT, TABLE_SCHEM, TABLE_NAME) %>% 
    subset(str_detect(TABLE_NAME, table_prefix)) %>% 
    dplyr::mutate(tissue = str_extract(TABLE_NAME, tissue_search)) %>% 
    subset(tissue %in% tissue_search)
  
  # loop through tissues and generate plot: expression vs strain by diet
  plots <- lapply( unique(data_sources$tissue), function(in_tissue){
    
    # subset names of SQL tables to the specified tissue
    sub_data_sources <- subset(data_sources, tissue == in_tissue)
    
    # generate SQL query to extract data for specified gene_symbol
    queries <- paste0("select * from ", sub_data_sources$TABLE_SCHEM, ".", sub_data_sources$TABLE_NAME, " where gene_symbol = '", gene_symbol, "'")
    
    # extract data from the SQL server
    data <- lapply(queries, function(q) sqlQuery(db, q))
    names(data) <- sub_data_sources$TABLE_NAME %>% 
      str_replace(table_prefix, "") %>% 
      str_replace(in_tissue, "") %>% 
      str_replace("(^_*)|(_$)", "")
    
    # pre-processing to obtain only the relevant columns
    data <- lapply(data, dplyr::select, matches("strain"), matches("expression"))
    strain_col <- unique( laply(data, function(x) str_subset(colnames(x), "strain")) )
    expression_col <- unique( laply(data, function(x) str_subset(colnames(x), "expression")) )
    
    # expression plots or box plots
    if(plot_boxplots){
      
      # add the name of the list as an identifier column
      data <- jn.general::rename_list(data)
      data <- lapply(data, function(x) tidyr::separate(x, file_id, c("diet", "sex"), sep = "[_]"))
      
      # make diet/sex columns; combine data to make plot data
      plot_data <- rbindlist(data)
      setnames(plot_data, expression_col, "expression")
      plot_data <- dplyr::select(plot_data, diet, sex, expression)
      
      # pre-processing for missing groups
      summ <- plot_data %>% 
        group_by(diet, sex) %>% 
        dplyr::summarise(n = n())
      
      if( nrow(summ) %% 2 != 0 ){
        
        # find the missing group
        non_dup_diet <- with(summ, diet[ !(duplicated(diet) | duplicated(diet, fromLast = TRUE)) ])
        non_dup_sex <- with(summ, sex[ !(duplicated(sex) | duplicated(sex, fromLast = TRUE)) ])
        
        # loop through potential 
        for(diet in non_dup_diet){
          for(sex in non_dup_sex){
            
            if( nrow(subset(summ, diet == diet, sex == sex)) ) plot_data <- rbind(plot_data, data.frame(diet = diet, sex = sex, expression = 0))
            
          }
        }
      }
      
      # make box plots
      g <- ggplot(data = plot_data, aes(sex, expression)) + 
        geom_boxplot(alpha = 0.8) +
        facet_grid(~diet) +
        labs(x = "Sex", y = "Gene Expression") +
        ggtitle(paste("Transcript Abundance:", gene_symbol, str_to_title(in_tissue), "Tissue")) +
        theme(legend.position = "bottom")
      
    } else{
      
      # if strain information available, merge on strains; otherwise do a plot each ordered
      if( group_by_strain & length(strain_col) == 1 ){
        
        # average strains fix for multiple probe_ids that match to the same gene
        duplicated_strains <- any(laply( data, function(x) any(duplicated(x[[strain_col]])) ))
        if(duplicated_strains){
          data <- lapply(data, function(x){
            
            y <- copy( data.table(x) )
            setnames(y, expression_col, "expr.x")
            
            y <- y %>% 
              group_by_(strain_col) %>% 
              dplyr::summarise(expr.x = mean(expr.x))
            
            setnames(y, "expr.x", expression_col)
            return(data.frame(y))
            
          })
        }
        
        # merge the chow and highfat data for plotting
        merged_data <- jn.general::merge_mult(data, by = strain_col, all = TRUE, suffixes = names(data))
        
        # reformat data for plotting: wide to long with descriptors
        plot_data <- dplyr::select(merged_data, strain, matches(expression_col))
        
        # reformat names
        old_names <- str_subset( colnames(plot_data), "expression")
        new_names <- str_replace( old_names, "expression_", "expression.")
        setnames(plot_data, old_names, new_names)
        
        # reformat data for plotting
        baseline_col <- str_subset(colnames(plot_data), "chow.*male|male.*chow")
        plot_data <- plot_data %>% 
          arrange_(baseline_col) %>% 
          mutate(order = 1:nrow(.)) %>% 
          melt(id.vars = c("strain", "order")) %>% 
          tidyr::separate(variable, c("expression", "diet"), sep = "[.]") 
        
      } else if( !group_by_strain | length(strain_col) == 0 ){
        
        # arrange data
        data <- lapply(data, function(x) arrange_(x, expression_col))
        
        # add extra space to the expression data as necessary
        need_to_add <- max( laply(data, nrow) )
        merged_data <- lapply(data, function(x) c(x[,expression_col], rep(NA, need_to_add - nrow(x)))) %>% 
          as.data.frame %>% 
          mutate(order = 1:nrow(.))
        
        # reformat for plotting
        plot_data <- melt(merged_data, id.vars = "order", variable.name = "diet")
        
      } else{
        
        error("All tables must have the same column names")
        
      }
      
      # additional pre-processing
      name_x_axis <- ifelse( group_by_strain & length(strain_col) == 1, str_to_title(strain_col), "")
      plot_data <- tidyr::separate(plot_data, diet, c("diet", "sex"), sep = "_")
      
      # generate plot
      g <- ggplot(data = plot_data, aes(order, value, color = diet)) + 
        geom_point(alpha = 0.8) +
        facet_grid(~sex) +
        scale_x_continuous(breaks = NULL) +
        labs(x = name_x_axis, y = "Gene Expression", color = "Type") +
        ggtitle(paste("Transcript Abundance:", gene_symbol, str_to_title(in_tissue), "Tissue")) +
        theme(legend.position = "bottom")
      
    }

    # return plots
    return(g)
  })
  names(plots) <- unique(data_sources$tissue)

  # close the database
  odbcClose(db)
  
  # return the plots in a list
  return(plots)
}



