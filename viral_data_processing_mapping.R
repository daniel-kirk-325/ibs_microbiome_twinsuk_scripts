## Author: 
  # Panayiotis Louca 

## Purpose of script: 
  # Create mapping file to facilitate downstream viral analyses 

## Date Created: 
  # 08 June 2026 

## Notes: 
  #  

## Clear environment 
  rm(list = ls()) 

## Set seed 
  set.seed(1234)

## Set functions: 
      
## load up packages: 

  ### core 
  library(tidyverse)
  library(readxl)

# -------------------------------------------------------------------------- # 
  
# ----------------- # 
#   Import data  ---- 
# ----------------- # 
  
  df <- read_delim("/Users/panayiotislouca/Documents/KCL_files/Data/TwinsUK/Microbiome/Phages/Database_independent_approach_revised/tpm abundance table/merged_tpm_b1_fixed_ids_lastline.tsv")
  
# -------------------------------------------------------------------------- # 
  
# ---------------------------------------- # 
#   Create mapping file of phage names  ---- 
# ---------------------------------------- # 
  
  df_map <- df %>%
    rename(phage_ID = ID) %>%
    select(phage_ID) %>%
    mutate(sudo_phage_ID = paste0("phage_", row_number()))

# -------------------------------------------------------------------------- # 
  
# ------------------------------------------- # 
#   Add in taxonomy and replication cycle  ---- 
# ------------------------------------------- # 
  
        ######   Taxonomy  ---- 
  df_taxonomy <- read_delim("/Users/panayiotislouca/Documents/KCL_files/Data/TwinsUK/Microbiome/Phages/Database_independent_approach_revised/viral contigs metadata/merged_genomad_taxonomy_final_contigs.tsv")

  df_map <- df_map %>%
    left_join(df_taxonomy, 
              by = c("phage_ID" = "seq_name")) %>%
    rename(taxonomy_lineage = lineage)
    
        ######  Replication cycle ---- 
  df_rep_cycle <- read_delim("/Users/panayiotislouca/Documents/KCL_files/Data/TwinsUK/Microbiome/Phages/Database_independent_approach_revised/viral contigs metadata/merged_bacphlip_life_cycle.tsv")
  
  df_map <- df_map %>%
    left_join(df_rep_cycle, 
              by = c("phage_ID" = "ID"))
  
# -------------------------------------------------------------------------- # 

# ------------------ # 
#   Save dataset  ---- 
# ------------------ # 
  
  write.csv(df_map,
            "/Users/panayiotislouca/Documents/KCL_files/Data/TwinsUK/Microbiome/Phages/Database_independent_approach_revised/viral contigs metadata/phage_contig_map_with_taxonomy_and_life_cycle.csv",
            row.names = FALSE)
  