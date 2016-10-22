#------------------------------------------------------------------------
# Title: Find Cholesterol Module Intersections Across Diets
# Date: June 2 2016
# Author: Jenny Nguyen
# Email: jnnguyen2@wisc.edu
#------------------------------------------------------------------------

# CONTENTS:

# NOTES:
#------------------------------------------------------------------------



library(reshape2)
library(stringr)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(WGCNA)
library(xlsx)
library(jn.general)


tmp_merge_mult <- function (..., by, all = FALSE, all.x = NULL, suffixes = letters) {
  l <- list(...)
  if (length(l) < 2) stop("At least 2 data frames must be supplied for merge")
  if (missing(by)) stop("Missing by argument")
  if (length(suffixes) < length(l)) stop("The number of suffixes must be equal than the number of data frames")
  if (!((length(by) == 1) | (length(l) - length(by) == 1))) stop("Length of by argument must be either 1 or 1 less than the total data frames")
  
  suffix <- paste0("_", suffixes)
  
  merged <- data.table::data.table(l[[1]])
  
  if (length(by) == 1) by <- rep(by, length(l) - 1)
  
  if(length(l) > 2){
    for ( i in 2:(length(l)-1) ) {
      if (is.null(all.x)) {
        merged <- merge(merged, l[[i]], by = by[[i - 1]], all = all, suffix = c(suffix[i - 1], ""))
      }
      else {
        merged <- merge(merged, l[[i]], by = by[[i - 1]], all.x = all.x, suffix = c(suffix[i - 1], ""))
      }
    }
    if(is.null(all.x)){
      merged <- merge(merged, l[[length(l)]], by = by[[length(l)-1]], all = all, suffix = c(suffix[length(l)-1], suffix[length(l)]))
    }  else{
      merged <- merge(merged, l[[length(l)]], by = by[[length(l)-1]], all.x = all.x, suffix = c(suffix[length(l)-1], suffix[length(l)]))
    }
  } else{
    if(is.null(all.x)){
      merged <- merge(l[[1]], l[[2]], by = by[[1]], all = all, suffix = c(suffix[1], suffix[2]))
    } else{
      merged <- merge(l[[1]], l[[2]], by = by[[1]], all.x = all.x, suffix = c(suffix[1], suffix[2]))
    }
  }

  return(merged)
}

cleanup <- function(x){
  x %>% 
    dplyr::select(gene_symbol) %>% 
    arrange(gene_symbol) %>% 
    distinct %>% 
    return
}

append_to_csv <- function(data, name){
  
  write.xlsx(data, file = file_name, sheetName = name, append = TRUE, row.names = FALSE)
}

# Intersections: Liver HMDP
#------------------------------------------------------------------------

# colors corresponding to cholesterol
  # chow male: black
  # highfat male: lightyellow
  # highfat female: greenyellow
  # reactive: black

chow_male <- fread("file:///E:/MICROARRAY DATA/INDIVIDUAL_NETWORKS/network_outputs/liver_chow_male_network.csv") %>% subset(clusterL == "black")
highfat_male <- fread("file:///E:/MICROARRAY DATA/INDIVIDUAL_NETWORKS/network_outputs/liver_highfat_male_network.csv") %>% subset(clusterL == "lightyellow")
highfat_female <- fread("file:///E:/MICROARRAY DATA/INDIVIDUAL_NETWORKS/network_outputs/liver_highfat_female_network.csv") %>% subset(clusterL == "greenyellow")
reactive <- fread("file:///E:/MICROARRAY DATA/REACTIVE_HDMP_NETWORKS/code_outputs/0.2_with_0.1_threshold/liver_network.csv") %>% subset(clusterL == "black")

intersecting_male <- merge_mult(highfat_male, chow_male, by = "probeid", suffixes = c("hf.m", "chow.m")) %>% 
  # clean up matching genes
  dplyr::select(probeid, matches("gene_symbol")) %>% 
  subset(gene_symbol_hf.m == gene_symbol_chow.m) %>% 
  mutate(gene_symbol = gene_symbol_chow.m) %>% 
  cleanup

intersecting_hf <- merge_mult(highfat_male, highfat_female, by = "probeid", suffixes = c("hf.m", "hf.f")) %>% 
  # clean up matching genes
  dplyr::select(probeid, matches("gene_symbol")) %>% 
  subset(gene_symbol_hf.m == gene_symbol_hf.f) %>% 
  mutate(gene_symbol = gene_symbol_hf.m) %>% 
  cleanup

# intersect genes
intersecting_genes <- merge_mult(chow_male, highfat_male, highfat_female, by = "probeid", suffixes = c("chow.m", "hf.m", "hf.f")) %>% 
  # clean up matching genes
  dplyr::select(probeid, matches("gene_symbol")) %>% 
  subset(gene_symbol_chow.m == gene_symbol_hf.m & gene_symbol_hf.m == gene_symbol_hf.f) %>% 
  mutate(gene_symbol = gene_symbol_chow.m) %>% 
  cleanup

intersecting_genes_w_reactive <- merge_mult(chow_male, highfat_male, highfat_female, reactive, by = "probeid", suffixes = c("chow.m", "hf.m", "hf.f", "reactive")) %>% 
  # clean up to matching genes
  dplyr::select(probeid, matches("gene_symbol")) %>% 
  subset(gene_symbol_chow.m == gene_symbol_hf.m & gene_symbol_hf.m == gene_symbol_hf.f & gene_symbol_hf.f == gene_symbol_reactive) %>% 
  mutate(gene_symbol = gene_symbol_chow.m) %>% 
  cleanup


# write intersecting genes file
file_name <- "E:/MICROARRAY DATA/cholesterol_intersections/cholesterol_intersection_HMDP_liver.xlsx"

write.xlsx(intersecting_genes_w_reactive, file = file_name, sheetName = "intersect_individual_reactive", append = FALSE, row.names = FALSE)
append_to_csv(intersecting_genes, name = "intersect_all")
append_to_csv(intersecting_hf, name = "intersect_highfat")
append_to_csv(intersecting_male, name = "intersect_male")
append_to_csv(chow_male, name = "chow.m_module_all")
append_to_csv(highfat_male, name = "hf.m_module_all")
append_to_csv(highfat_female, name = "hf.f_module_all")
append_to_csv(reactive, name = "reactive_module_all")


# Intersections: CM RNA data
#------------------------------------------------------------------------

# colors corresponding to cholesterol
  # chow male: darkorange2
  # chow female: darkorange2
  # highfat male: cyan
  # highfat female: purple

chow_male <- fread("file:///E:/MICROARRAY DATA/ChickMungeretal2016/network_outputs_rna/chow_M_network.csv") %>% subset(clusterL == "darkorange2")
chow_female <- fread("file:///E:/MICROARRAY DATA/ChickMungeretal2016/network_outputs_rna/chow_F_network.csv") %>% subset(clusterL == "darkorange2")
hf_male <- fread("file:///E:/MICROARRAY DATA/ChickMungeretal2016/network_outputs_rna/HF_M_network.csv") %>% subset(clusterL == "cyan")
hf_female <- fread("file:///E:/MICROARRAY DATA/ChickMungeretal2016/network_outputs_rna/HF_F_network.csv") %>% subset(clusterL == "purple")


CM_all <- merge_mult(chow_male, chow_female, hf_male, hf_female, by = "probeid", suffixes = c("chow.m", "chow.f", "hf.m", "hf.f")) %>% 
  # clean up to matching genes
  dplyr::select(probeid, matches("gene_symbol")) %>% 
  subset(gene_symbol_chow.m == gene_symbol_chow.f & gene_symbol_chow.f == gene_symbol_hf.m & gene_symbol_hf.m == gene_symbol_hf.f) %>% 
  mutate(gene_symbol = gene_symbol_chow.m) %>% 
  cleanup

CM_HF <- merge_mult(hf_male, hf_female, by = "probeid", suffixes = c("hf.m", "hf.f")) %>% 
  # clean up to matching genes
  dplyr::select(probeid, matches("gene_symbol")) %>% 
  subset(gene_symbol_hf.m == gene_symbol_hf.f) %>% 
  mutate(gene_symbol = gene_symbol_hf.m) %>% 
  cleanup

CM_C <- merge_mult(chow_male, chow_female, by = "probeid", suffixes = c("chow.m", "chow.f")) %>% 
  # clean up to matching genes
  dplyr::select(probeid, matches("gene_symbol")) %>% 
  subset(gene_symbol_chow.m == gene_symbol_chow.f) %>% 
  mutate(gene_symbol = gene_symbol_chow.m) %>% 
  cleanup

CM_male_all <- merge_mult(chow_male, hf_male, by = "probeid", suffixes = c("chow.m", "hf.m")) %>% 
  # clean up to matching genes
  dplyr::select(probeid, matches("gene_symbol")) %>% 
  subset(gene_symbol_chow.m == gene_symbol_hf.m) %>% 
  mutate(gene_symbol = gene_symbol_chow.m) %>% 
  cleanup

CM_female_all <- merge_mult(chow_female, hf_female, by = "probeid", suffixes = c("chow.f", "hf.f")) %>% 
  # clean up to matching genes
  dplyr::select(probeid, matches("gene_symbol")) %>% 
  subset(gene_symbol_chow.f == gene_symbol_hf.f) %>% 
  mutate(gene_symbol = gene_symbol_chow.f) %>% 
  cleanup


# write intersecting genes file
file_name <- "E:/MICROARRAY DATA/cholesterol_intersections/cholesterol_intersection_CM.xlsx"

write.xlsx(CM_all, file = file_name, sheetName = "intersect_all", append = FALSE, row.names = FALSE)
append_to_csv(CM_C, name = "intersect_chow")
append_to_csv(CM_HF, name = "intersect_highfat")
append_to_csv(CM_female_all, name = "intersect_female")
append_to_csv(CM_male_all, name = "intersect_male")
append_to_csv(chow_male, name = "chow_male_module_all")
append_to_csv(chow_female, name = "chow_female_module_all")
append_to_csv(hf_male, name = "highfat_male_module_all")
append_to_csv(hf_female, name = "highfat_female_module_all")


# Intersections: CM protein data
#------------------------------------------------------------------------

# colors corresponding to cholesterol
  # chow male: white
  # chow female: NONE
  # highfat male: darkred
  # highfat female: saddlebrown

chow_male <- fread("file:///E:/MICROARRAY DATA/ChickMungeretal2016/network_outputs_protein/chow_M_network.csv") %>% subset(clusterL == "white")
# chow_female <- fread("file:///E:/MICROARRAY DATA/ChickMungeretal2016/network_outputs_protein/chow_F_network.csv") %>% subset(clusterL == "NONE")
hf_male <- fread("file:///E:/MICROARRAY DATA/ChickMungeretal2016/network_outputs_protein/HF_M_network.csv") %>% subset(clusterL == "darkred")
hf_female <- fread("file:///E:/MICROARRAY DATA/ChickMungeretal2016/network_outputs_protein/HF_F_network.csv") %>% subset(clusterL == "saddlebrown")


CM_all <- merge_mult(chow_male, hf_male, hf_female, by = "probeid", suffixes = c("chow.m", "hf.m", "hf.f")) %>% 
  # clean up to matching genes
  dplyr::select(probeid, matches("gene_symbol")) %>% 
  subset(gene_symbol_chow.m == gene_symbol_hf.m & gene_symbol_hf.m == gene_symbol_hf.f) %>% 
  mutate(gene_symbol = gene_symbol_chow.m) %>% 
  cleanup

CM_HF <- merge_mult(hf_male, hf_female, by = "probeid", suffixes = c("hf.m", "hf.f")) %>% 
  # clean up to matching genes
  dplyr::select(probeid, matches("gene_symbol")) %>% 
  subset(gene_symbol_hf.m == gene_symbol_hf.f) %>% 
  mutate(gene_symbol = gene_symbol_hf.m) %>% 
  cleanup

CM_male_all <- merge_mult(chow_male, hf_male, by = "probeid", suffixes = c("chow.m", "hf.m")) %>% 
  # clean up to matching genes
  dplyr::select(probeid, matches("gene_symbol")) %>% 
  subset(gene_symbol_chow.m == gene_symbol_hf.m) %>% 
  mutate(gene_symbol = gene_symbol_chow.m) %>% 
  cleanup


# write intersecting genes file
file_name <- "E:/MICROARRAY DATA/cholesterol_intersections/cholesterol_intersection_CM_protein.xlsx"

write.xlsx(CM_all, file = file_name, sheetName = "intersect_all", append = FALSE, row.names = FALSE)
append_to_csv(CM_HF, name = "intersect_highfat")
append_to_csv(CM_male_all, name = "intersect_male")
append_to_csv(chow_male, name = "chow_male_module_all")
append_to_csv(data.frame(NOTE = "chow female does not have a cholesterol module"), name = "chow_female_module_all")
append_to_csv(hf_male, name = "highfat_male_module_all")
append_to_csv(hf_female, name = "highfat_female_module_all")
