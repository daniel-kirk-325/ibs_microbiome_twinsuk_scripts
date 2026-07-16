  ## Author: 
    # Panayiotis Louca 

    # Modified by: Dan Kirk
  
  ## Purpose of script: 
    # Merges the viral data processed in different batches
    # Calculates alpha diversity 
    # Saves merged dataset
    # Processes merged dataset to remove repeated GSID entries, leaving one unique entry per sample
  
  ## Date Created: 
    # 08 June 2026 
  
  ## Notes: 
    # Original name: Combine_TPM_phage_batches.R 
  
  ## Clear environment 
    rm(list = ls()) 
  
  ## Set seed 
    set.seed(1234)
  
  ## Set functions: 
        
  ## load up packages: 
  
    ### core 
    library(tidyverse)
    library(readxl)
    
    # vegan 
    library(vegan)
  
# -------------------------------------------------------------------------- # 
    
  # ------------------------- # 
  #   Import TPM datasets  ---- 
  # ------------------------- # 
   
# Import all normalised phage batches 
    
    # base path 
    base_path <- "base_path"
    
    # Get all files matching the pattern 
    file_pattern <- "merged_tpm_b[0-9]+_fixed_ids_lastline\\.tsv$"
    file_list <- list.files(path = base_path, pattern = file_pattern, full.names = TRUE)
    
    # Import each file as a separate dataframe 
    for (file in file_list) {
      batch_id <- sub(".*merged_tpm_(b[0-9]+)_fixed_ids_lastline\\.tsv$", "\\1", basename(file))
      df_name <- paste0("df_", batch_id)
      assign(df_name, read.delim(file, sep = "\t", header = TRUE), envir = .GlobalEnv)
    }   

# -------------------------------------------------------------------------- # 

# ------------------------------- # 
#   Import phage mapping file  ---- 
# ------------------------------- # 
    
    df_map <- read.csv("phage_contig_map_with_taxonomy_and_life_cycle.csv")
    
# -------------------------------------------------------------------------- # 

# ------------------------ # 
#   Process each batch  ---- 
# ------------------------ # 
    
        ######   BATCH 1  ---- 
    
  head(df_b1)
  
  df_b1t <- df_b1 %>%
    t() %>% # transpose data 
    as.data.frame() %>%
    `colnames<-`(.[1, ]) %>% # use first row as column names 
    slice(-1) %>% # drop first row 
    rownames_to_column(var = "sample_id")
  
  # Clean sample ID to extract GSID 
  df_b1t <- df_b1t %>%
    mutate(gsid = stringr::str_extract(sample_id, "(?<=_)\\d+\\.\\d{3}\\.\\d{3}(?=_)"),
           gsid = str_replace_all(gsid, "\\.", "_")) %>%
    relocate(gsid, .after = sample_id)
  
  unique(df_b1t$sample_id)
  unique(df_b1t$gsid)
  
  # replace phage names with sudo phage names from mapping 
  length(unique(grep("nr_contigs_", names(df_b1t)))) # 116,740 
  
  phage_idx <- match(df_map$phage_ID, names(df_b1t))
  
  names(df_b1t)[phage_idx[!is.na(phage_idx)]] <- df_map$sudo_phage_ID[!is.na(phage_idx)]
  
  # check phage class  
  table(sapply(df_b1t[ ,grep("^phage_", names(df_b1t))], class))
  
  # clean up phage columns 
  phage_idx <- grep("^phage_", names(df_b1t))
  
  df_b1t[phage_idx] <- lapply(df_b1t[phage_idx], function(x) as.numeric(trimws(x)))
  
  # Check NA's 
  na_counts <- sapply(df_b1t[phage_idx], function(x) sum(is.na(x)))
  sort(na_counts, decreasing = TRUE)[1:10]
  
  # check Row sums 
  summary(rowSums(df_b1t[ ,grep("^phage_", names(df_b1t))])) 
  
  # Extract TPM columns (those containing "phage_") 
  tpm_columns <- grep("^phage_", names(df_b1t))
  tpm_data <- df_b1t[, tpm_columns]
  
  # Convert to relative abundances (each row sums to 100) but those with rowSum == 0 remain 0 
  rs <- rowSums(tpm_data, na.rm = TRUE)
  
  relative_abundance <- tpm_data               # start with original values
  relative_abundance[rs > 0, ] <- sweep(
    tpm_data[rs > 0, ],
    1,
    rs[rs > 0],
    "/"
  ) * 100
  
  relative_abundance[rs == 0, ] <- 0
  
  # Verify conversion worked - all rows should sum to 100 or 0 
  summary(rowSums(relative_abundance)) 
  
  # Replace the TPM columns in original dataframe 
  df_b1t[, tpm_columns] <- relative_abundance

# -------------------------------------------------------------------------- #  
  
  ######   BATCH 2  ---- 
  
  head(df_b2)
  
  df_b2t <- df_b2 %>%
    t() %>% # transpose data 
    as.data.frame() %>%
    `colnames<-`(.[1, ]) %>% # use first row as column names 
    slice(-1) %>% # drop first row 
    rownames_to_column(var = "sample_id")
  
  # Clean sample ID to extract GSID 
  df_b2t <- df_b2t %>%
    mutate(gsid = stringr::str_extract(sample_id, "(?<=_)\\d+\\.\\d{3}\\.\\d{3}(?=_)"),
           gsid = str_replace_all(gsid, "\\.", "_")) %>%
    relocate(gsid, .after = sample_id)
  
  unique(df_b2t$sample_id)
  unique(df_b2t$gsid)
  
  # replace phage names with sudo phage names from mapping 
  length(unique(grep("nr_contigs_", names(df_b2t)))) # 116,740 
  
  phage_idx <- match(df_map$phage_ID, names(df_b2t))
  
  names(df_b2t)[phage_idx[!is.na(phage_idx)]] <- df_map$sudo_phage_ID[!is.na(phage_idx)]
  
  # check phage class  
  table(sapply(df_b2t[ ,grep("^phage_", names(df_b2t))], class))
  
  # clean up phage columns 
  phage_idx <- grep("^phage_", names(df_b2t))
  
  df_b2t[phage_idx] <- lapply(df_b2t[phage_idx], function(x) as.numeric(trimws(x)))
  
  # Check NA's 
  na_counts <- sapply(df_b2t[phage_idx], function(x) sum(is.na(x)))
  sort(na_counts, decreasing = TRUE)[1:10]
  
  # check Row sums 
  summary(rowSums(df_b2t[ ,grep("^phage_", names(df_b2t))])) 
  
  # Extract TPM columns (those containing "phage_") 
  tpm_columns <- grep("^phage_", names(df_b2t))
  tpm_data <- df_b2t[, tpm_columns]
  
  # Convert to relative abundances (each row sums to 100) but those with rowSum == 0 remain 0 
  rs <- rowSums(tpm_data, na.rm = TRUE)
  
  relative_abundance <- tpm_data               # start with original values
  relative_abundance[rs > 0, ] <- sweep(
    tpm_data[rs > 0, ],
    1,
    rs[rs > 0],
    "/"
  ) * 100
  
  relative_abundance[rs == 0, ] <- 0
  
  # Verify conversion worked - all rows should sum to 100 or 0 
  summary(rowSums(relative_abundance)) 
  
  # Replace the TPM columns in original dataframe 
  df_b2t[, tpm_columns] <- relative_abundance
  
# -------------------------------------------------------------------------- #  
  
  ######   BATCH 3  ---- 
  
  head(df_b3)
  
  df_b3t <- df_b3 %>%
    t() %>% # transpose data 
    as.data.frame() %>%
    `colnames<-`(.[1, ]) %>% # use first row as column names 
    slice(-1) %>% # drop first row 
    rownames_to_column(var = "sample_id")
  
  # Clean sample ID to extract GSID 
  df_b3t <- df_b3t %>%
    mutate(gsid = stringr::str_extract(sample_id, "(?<=_)\\d+\\.\\d{3}\\.\\d{3}(?=_)"),
           gsid = str_replace_all(gsid, "\\.", "_")) %>%
    relocate(gsid, .after = sample_id)
  
  unique(df_b3t$sample_id)
  unique(df_b3t$gsid)
  
  # replace phage names with sudo phage names from mapping 
  length(unique(grep("nr_contigs_", names(df_b3t)))) # 116,740 
  
  phage_idx <- match(df_map$phage_ID, names(df_b3t))
  
  names(df_b3t)[phage_idx[!is.na(phage_idx)]] <- df_map$sudo_phage_ID[!is.na(phage_idx)]
  
  # check phage class  
  table(sapply(df_b3t[ ,grep("^phage_", names(df_b3t))], class))
  
  # clean up phage columns 
  phage_idx <- grep("^phage_", names(df_b3t))
  
  df_b3t[phage_idx] <- lapply(df_b3t[phage_idx], function(x) as.numeric(trimws(x)))
  
  # Check NA's 
  na_counts <- sapply(df_b3t[phage_idx], function(x) sum(is.na(x)))
  sort(na_counts, decreasing = TRUE)[1:10]
  
  # check Row sums 
  summary(rowSums(df_b3t[ ,grep("^phage_", names(df_b3t))])) 
  
  # Extract TPM columns (those containing "phage_") 
  tpm_columns <- grep("^phage_", names(df_b3t))
  tpm_data <- df_b3t[, tpm_columns]
  
  # Convert to relative abundances (each row sums to 100) but those with rowSum == 0 remain 0 
  rs <- rowSums(tpm_data, na.rm = TRUE)
  
  relative_abundance <- tpm_data               # start with original values
  relative_abundance[rs > 0, ] <- sweep(
    tpm_data[rs > 0, ],
    1,
    rs[rs > 0],
    "/"
  ) * 100
  
  relative_abundance[rs == 0, ] <- 0
  
  # Verify conversion worked - all rows should sum to 100 or 0 
  summary(rowSums(relative_abundance)) 
  
  # Replace the TPM columns in original dataframe 
  df_b3t[, tpm_columns] <- relative_abundance
  
# -------------------------------------------------------------------------- #  
  
  ######   BATCH 4  ---- 
  
  head(df_b4)
  
  df_b4t <- df_b4 %>%
    t() %>% # transpose data 
    as.data.frame() %>%
    `colnames<-`(.[1, ]) %>% # use first row as column names 
    slice(-1) %>% # drop first row 
    rownames_to_column(var = "sample_id")
  
  # Clean sample ID to extract GSID 
  df_b4t <- df_b4t %>%
    mutate(gsid = stringr::str_extract(sample_id, "(?<=_)\\d+\\.\\d{3}\\.\\d{3}\\w*(?=_)"),
           gsid = str_replace_all(gsid, "\\.", "_")) %>%
    relocate(gsid, .after = sample_id)
  
  unique(df_b4t$sample_id)
  unique(df_b4t$gsid)
  
  # replace phage names with sudo phage names from mapping 
  length(unique(grep("nr_contigs_", names(df_b4t)))) # 116,740 
  
  phage_idx <- match(df_map$phage_ID, names(df_b4t))
  
  names(df_b4t)[phage_idx[!is.na(phage_idx)]] <- df_map$sudo_phage_ID[!is.na(phage_idx)]
  
  # check phage class  
  table(sapply(df_b4t[ ,grep("^phage_", names(df_b4t))], class))
  
  # clean up phage columns 
  phage_idx <- grep("^phage_", names(df_b4t))
  
  df_b4t[phage_idx] <- lapply(df_b4t[phage_idx], function(x) as.numeric(trimws(x)))
  
  # Check NA's 
  na_counts <- sapply(df_b4t[phage_idx], function(x) sum(is.na(x)))
  sort(na_counts, decreasing = TRUE)[1:10]
  
  # check Row sums 
  summary(rowSums(df_b4t[ ,grep("^phage_", names(df_b4t))])) 
  
  # Extract TPM columns (those containing "phage_") 
  tpm_columns <- grep("^phage_", names(df_b4t))
  tpm_data <- df_b4t[, tpm_columns]
  
  # Convert to relative abundances (each row sums to 100) but those with rowSum == 0 remain 0 
  rs <- rowSums(tpm_data, na.rm = TRUE)
  
  relative_abundance <- tpm_data               # start with original values
  relative_abundance[rs > 0, ] <- sweep(
    tpm_data[rs > 0, ],
    1,
    rs[rs > 0],
    "/"
  ) * 100
  
  relative_abundance[rs == 0, ] <- 0
  
  # Verify conversion worked - all rows should sum to 100 or 0 
  summary(rowSums(relative_abundance))
  
  # Replace the TPM columns in original dataframe 
  df_b4t[, tpm_columns] <- relative_abundance
  
# -------------------------------------------------------------------------- #  
  
  ######   BATCH 5  ---- 
  
  head(df_b5)
  
  df_b5t <- df_b5 %>%
    t() %>% # transpose data 
    as.data.frame() %>%
    `colnames<-`(.[1, ]) %>% # use first row as column names 
    slice(-1) %>% # drop first row 
    rownames_to_column(var = "sample_id")
  
  # Clean sample ID to extract GSID 
  df_b5t <- df_b5t %>%
    mutate(gsid = stringr::str_extract(sample_id, "(?<=_)\\d+\\.\\d{3}\\.\\d{3}(?=_)"),
           gsid = str_replace_all(gsid, "\\.", "_")) %>%
    relocate(gsid, .after = sample_id)
  
  unique(df_b5t$sample_id)
  unique(df_b5t$gsid)
  
  # replace phage names with sudo phage names from mapping 
  length(unique(grep("nr_contigs_", names(df_b5t)))) # 116,740 
  
  phage_idx <- match(df_map$phage_ID, names(df_b5t))
  
  names(df_b5t)[phage_idx[!is.na(phage_idx)]] <- df_map$sudo_phage_ID[!is.na(phage_idx)]
  
  # check phage class  
  table(sapply(df_b5t[ ,grep("^phage_", names(df_b5t))], class))
  
  # clean up phage columns 
  phage_idx <- grep("^phage_", names(df_b5t))
  
  df_b5t[phage_idx] <- lapply(df_b5t[phage_idx], function(x) as.numeric(trimws(x)))
  
  # Check NA's 
  na_counts <- sapply(df_b5t[phage_idx], function(x) sum(is.na(x)))
  sort(na_counts, decreasing = TRUE)[1:10]
  
  # check Row sums 
  summary(rowSums(df_b5t[ ,grep("^phage_", names(df_b4t))])) 
  
  # Extract TPM columns (those containing "phage_") 
  tpm_columns <- grep("^phage_", names(df_b5t))
  tpm_data <- df_b5t[, tpm_columns]
  
  # Convert to relative abundances (each row sums to 100) but those with rowSum == 0 remain 0 
  rs <- rowSums(tpm_data, na.rm = TRUE)
  
  relative_abundance <- tpm_data               # start with original values
  relative_abundance[rs > 0, ] <- sweep(
    tpm_data[rs > 0, ],
    1,
    rs[rs > 0],
    "/"
  ) * 100
  
  relative_abundance[rs == 0, ] <- 0
  
  # Verify conversion worked - all rows should sum to 100 or 0 
  summary(rowSums(relative_abundance))
  
  # Replace the TPM columns in original dataframe 
  df_b5t[, tpm_columns] <- relative_abundance

# -------------------------------------------------------------------------- #  
  
  # -------------------- # 
  #   Merge datasets  ---- 
  # -------------------- # 
  
  # merge them together 
  head(df_b1t)[1:10]
  head(df_b2t)[1:10]
  head(df_b3t)[1:10]
  head(df_b4t)[1:10]
  head(df_b5t)[1:10]
  
  # Rbind the dataframes & add batch variable 
  df_all <- list(
    df_b1t %>% mutate(phage_batch = "1"),
    df_b2t %>% mutate(phage_batch = "2"),
    df_b3t %>% mutate(phage_batch = "3"),
    df_b4t %>% mutate(phage_batch = "4"),
    df_b5t %>% mutate(phage_batch = "5")
  ) %>%
    bind_rows() %>%
    relocate(phage_batch, .after = gsid)
  
  length(unique(df_all$sample_id)) # 5,124 
  length(unique(df_all$gsid)) # 2,301 
  
  phage_cols <- grep("^phage_[0-9]+", names(df_all))
  
  summary(rowSums(df_all %>% select(all_of(phage_cols))))
  
  sum(rowSums(df_all %>% select(all_of(phage_cols)), na.rm = TRUE) == 0) # 123 samples have no phages detected 

# -------------------------------------------------------------------------- #  
  
# ------------------------------------ # 
#   Create alpha diversity metrics  ---- 
# ------------------------------------ # 
  
  phage_cols <- grep("^phage_[0-9]+", names(df_all))
  
  df_all <- df_all %>%
    mutate(
      phage_alpha_div_shannon = vegan::diversity(df_all[ ,phage_cols], index = "shannon"),
      phage_alpha_div_taxonomic_richness = vegan::specnumber(df_all[ ,phage_cols])
    ) %>%
    relocate(
      phage_alpha_div_shannon,
      phage_alpha_div_taxonomic_richness,
      .after = phage_batch
    )
  
# -------------------------------------------------------------------------- # 

# --------------- # 
#   Save full file  ---- 
# --------------- # 
  
  write.csv(df_all,
            file.path("merged_tpm_COMBINED.csv"),
            row.names = FALSE)
 
  df_all <- read.csv("merged_tpm_COMBINED.csv")
  

  # --------------- # 
#   Save unique file  ---- 
# --------------- #  
  
  # Here, repeats for the same GSID in different lanes are dropped, leaving one sample per GSID
  #
  # In cases where a GSID has some runs with all zeros and others with non-zeros, a non-zero run is randomly selected
  #
  # Otherwise, one of the runs is randomly selected
  

  # Re-establish phage_cols vector
  phage_cols <- names(df_all)[grepl("^phage_[0-9]+$", names(df_all))]
  
  
  # Get sample IDs of those that are all zeros
  zero_sums <- df_all$sample_id[rowSums(df_all[phage_cols]) == 0]
  
  
  # Get GSIDs of these samples 
  zero_gsids <- unique(df_all$gsid[df_all$sample_id %in% zero_sums])
  
  
  # Save those GSIDs with only zeros for all runs 
  non_zero_gsids <- unique(df_all$gsid[!df_all$sample_id %in% zero_sums])
  all_zeros <- zero_gsids[!zero_gsids %in% non_zero_gsids]  # 31 samples had all zeros for all runs
  
  
  # Check this that all non-zeros are all 0 
  sum(df_all[df_all$gsid %in% all_zeros, phage_cols]) # Equals zero - expected
  
  
  # Check that all remaining gsids in zero_gsids have at least one non-zero row
  mixed_gsids <- zero_gsids[!zero_gsids %in% all_zeros]
  hits <- mixed_gsids[sapply(mixed_gsids, function(i) {all(rowSums(df_all[df_all$gsid == i, phage_cols]) == 0)})]
  if (length(hits) > 0) print(hits) else message("Check passed: no mixed_gsids are all zero")
  
  
  # Drop zero runs from df_all in those with a mix of zeros and non-zeros
  zero_sums_id_gsids <- df_all[df_all$sample_id %in% zero_sums, c('sample_id', 'gsid')]  # Get GSIDs of all samples with zeros
  zero_mixed_id_gsids <- zero_sums_id_gsids[!zero_sums_id_gsids$gsid %in% all_zeros, ]   # Drop GSIDs who are all zeros
  df_all_reduced <- df_all[!df_all$sample_id %in% zero_mixed_id_gsids$sample_id, ]       # Drop sample IDs that are all zeros in GSIDs with a mix of zero and non-zero entries
  
  rm(df_all)
  
  # Randomly select one sample from all remaining repeats
  set.seed(123)
  sampled_ids <- df_all_reduced[, c('sample_id', 'gsid')] |>
    group_by(gsid) |>
    slice_sample(n = 1) |>
    pull(sample_id)
  
  # Make unique df, where all GSIDs have only one sample
  df_unique <- df_all_reduced[df_all_reduced$sample_id %in% sampled_ids, ]
  
  rm(df_all_reduced)
  
  # Save dataframe
  write.csv(df_unique,
            file.path("merged_tpm_COMBINED_UNIQUE.csv"),
            row.names = FALSE)
  
  

  













