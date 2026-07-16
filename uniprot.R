library(dplyr)        # 1.1.4
library(tibble)
library(stringr)      # 1.5.2
require(lme4)         # 1.1.37
library(compositions) # 2.0.9
library(ggplot2)


# Functions ---------------------------------------------------------------


lmm_genes <- function(input_df, var_list, f_null, pathname) {
  
  # Check if the folder exists
  if (!file.exists(dirname(pathname))) {stop(paste("Error: Folder does not exist"))}
  
  null <- glmer(f_null, data = input_df, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), nAGQ = 40)
  
  results_df <- data.frame(EntryName = character(0), BETA = numeric(0), SE = numeric(0), OR = numeric(0), 
                           Lower = numeric(0), Upper = numeric(0), warnings = numeric(0), PVAL = numeric(0))
  
  count <- 0 
  for (v in var_list) {
    
    print(match(v, var_list)*100 / length(var_list))
    #count <- count + 1
    #print(count)
    
    fsplit <- strsplit(f_null, '~')
    f1 <- as.formula(paste(fsplit[[1]][1], '~', v, '+', fsplit[[1]][2]))
    
    full <- glmer(f1, data = input_df, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), nAGQ = 40)
    
    full_sum <- summary(full)
    
    warnings <- length(full_sum$optinfo$warnings) + length(full_sum$optinfo$conv$lme4$messages) + length(full_sum$fitMsgs)
    
    pos <- 2 # Variable is second in the output (including intercept)
    
    Vcov <- vcov(full, useScale = FALSE)
    BETA <- full@beta[pos]
    se1 <- sqrt(diag(Vcov))
    SE <- se1[pos]
    OR <- exp(BETA)
    LL <- exp(BETA - 1.96 * SE)
    UL <- exp(BETA + 1.96 * SE)
    PVAL <- anova(null, full)$`Pr(>Chisq)`[2]
    
    # Add the results to the results dataframe
    row <- data.frame(BETA = BETA, SE = SE, OR = OR, Lower = LL, Upper = UL, warnings = warnings, PVAL = PVAL)
    rownames(row) <- v
    results_df <- rbind(results_df, row)
  }
  
  results_df_ibs <- results_df 
  
  results_df_ibs <- round(results_df_ibs[, names(results_df_ibs) != 'PVAL'], 4)
  results_df_ibs['PVAL'] <- round(results_df$PVAL, 8)
  results_df_ibs <- results_df_ibs %>% arrange(PVAL)
  
  
  # Calculate q-values using the Benjamini-Hochberg procedure
  q_values <- p.adjust(results_df_ibs$PVAL, method = "fdr")
  bonf <- p.adjust(results_df_ibs$PVAL, method = "bonferroni")
  
  # Add q-values to the DataFrame
  results_df_ibs$q_value <- q_values
  results_df_ibs$bonf <- bonf
  
  results_df_ibs$warnings <- ifelse(results_df_ibs$warnings > 0, 1, 0)
  
  # Save the results to a CSV file
  write.csv(results_df_ibs, pathname)
  
}


# Data pre-processing -----------------------------------------------------

## Load data
df <- read.csv('ibs_microbiome_twins_drug_processed.csv')[, -1]

## Drop males
df <- df[df$SEX != 'M', ]

## Keep only those with phage data
phage_ids <- read.csv("viroprofiler_iids.csv")[, -1]
df <- df[df$iid %in% phage_ids, ]

## Add family ID
df$fid <- as.integer(substr(df$iid, 1, nchar(df$iid) - 1))

## Set batch as factor 
df$batch <- as.factor(df$batch)


# UniProt pathways  -------------------------------------------------------

df_gene_fam <- vroom::vroom("genefamily_uniprot_all.tsv") 

## Keep only those in our dataset

prefix_ab_cpm <- paste0(df$prefix, '_Abundance-CPM')
df_gene_fam <- cbind(df_gene_fam[, c("uniprotid", "EC number", "Protein names")], df_gene_fam[prefix_ab_cpm])

## Get minimum non-zero for later analyses 
sample_abundance <- unlist(df_gene_fam[grepl('Abun', names(df_gene_fam))])
one_tenth_min_non_zero <- min(sample_abundance[sample_abundance > 0])/10
rm(sample_abundance)

## Some entries have multiple ECs per uniprotid, meaning if we group by EC number, these will be left out. 
## Solution: Add new columns to df_gene_fam to denote each individual EC number 
## Then, for all unique EC numbers, do colsums for all the samples 

## Get all separate ECs as a df
ec_split <- setNames(strsplit(df_gene_fam$`EC number`, ";"), df_gene_fam$uniprotid)

## Get unique EC numbers and do sum the same columns if any of the newly added 
## columns contain the EC number

ec_to_df <- function(ec_split) {
  # Find the maximum length of the lists
  max_length <- max(sapply(ec_split, length))
  
  # Convert the list to a data frame, filling with NA for shorter vectors
  df <- do.call(rbind, lapply(ec_split, function(x) {
    length(x) <- max_length
    x
  }))
  
  # Convert to a data frame and set the row names
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  rownames(df) <- names(ec_split)
  
  # Rename columns
  colnames(df) <- paste0("C", seq_len(ncol(df)))
  
  df <- tibble::rownames_to_column(df, var = "uniprotid")
  rownames(df) <- NULL
  
  return(df)
}
ec_df <- ec_to_df(ec_split)

df_gene_fam <- merge(df_gene_fam, ec_df, by = 'uniprotid')

## Sum the sample columns of all rows containing each unique EC
unique_ecs <- unique(unlist(df_gene_fam[names(ec_df)[2:length(ec_df)]]))
unique_ecs <- unique_ecs[!is.na(unique_ecs)]
C_cols <- names(ec_df)[names(ec_df) != 'uniprotid']

rows <- list()
for (i in unique_ecs) {
  print(match(i, unique_ecs))
  # Get rows which contain ec
  one_ec_df <- df_gene_fam %>% filter(if_any(all_of(C_cols), ~ . == i))
  
  # Select the columns containing 'Abundance' and summarize
  summed_abundance <- one_ec_df %>% 
    select(contains('Abundance')) %>% 
    summarise_all(sum, na.rm = TRUE) %>% 
    unlist()
  
  # Store the summed_abundance as a named vector with the EC number as the name
  rows[[i]] <- summed_abundance
}

# Convert the list to a dataframe and add the EC_number column
combined_proteins <- bind_rows(rows, .id = "EC_number")
combined_proteins$EC_number <- unique_ecs

# Transpose df
combined_proteins <- as.data.frame(t(combined_proteins))

# Set protein names to column names
names(combined_proteins) <- combined_proteins[1, ]

# Remove first row (redundant)
combined_proteins <- combined_proteins[-1, ]

# Save rownames (i.e., prefix) first column and remove redundant information
# to add to dataframe below
prefix <- sub("_[^_]+$", "", paste0(rownames(combined_proteins)))


# Make enzyme concentrations numeric
combined_proteins <- as.data.frame(lapply(combined_proteins[, names(combined_proteins) != 'prefix'], as.numeric))


## Save low prevalence columns
to_drop <- names(combined_proteins)[(colSums(combined_proteins == 0)) > 0.9 * nrow(combined_proteins)]


## CLR transform uniprot dataset
combined_proteins <- as.data.frame(compositions::clr(combined_proteins + one_tenth_min_non_zero))


## Drop low prevalence columns
combined_proteins_filtered <- combined_proteins[, !names(combined_proteins) %in% to_drop]
rm(combined_proteins)

# Map existing names to pseudonames to suit model formula
original_enzyme_names <- names(combined_proteins_filtered)
enzyme_indices <- paste0('X', 1:length(names(combined_proteins_filtered)))
names(combined_proteins_filtered) <- enzyme_indices
mapping <- setNames(original_enzyme_names, enzyme_indices)


# Add IBS, batch, fid, age and BMI 
combined_proteins_filtered$prefix <- gsub('X', '', prefix)
input_df <- merge(df[, c('prefix', 'iid', 'dov', 'IBS_overall', 'batch', 'fid', 'age_Mb_sample', 'bmi_Mb_sample')], combined_proteins_filtered, by = 'prefix')
input_df[, c("age_scaled", "BMI_scaled")] <- scale(input_df[, c("age_Mb_sample", "bmi_Mb_sample")])

f_null <- paste("IBS_overall ~ age_scaled + BMI_scaled + batch + (1 | fid)")

### Run next 6 lines all together
path <- "uniprot.csv"
lmm_genes(input_df, enzyme_indices, f_null, path)
results_df <- read.csv(path)
results_df$EC_number <- mapping[results_df$X]
#write.csv(results_df, path)



# Drug sensitivity analysis -----------------------------------------------


top_twins_df <- read.csv('uniprot.csv')
top_twins <- top_twins_df$X[top_twins_df$q_value < 0.05]

drug_vars <- c("Proton.Pump.Inhibitors", "Anti_Depressants.anxiolytics", "Antispasmodics.antimotility", "Analgesics", 'Laxatives')

input_df <- merge(input_df, df[, c('prefix', drug_vars)], by = 'prefix')

f_null <- paste("IBS_overall ~ age_scaled + BMI_scaled + batch + (1 | fid) + ", paste0(drug_vars, collapse = '+'))

### Run next 6 lines all together
path <- "uniprot_drug_adjusted.csv"
lmm_genes(input_df, top_twins, f_null, path)
results_df <- read.csv(path)
results_df$EC_number <- mapping[results_df$X]
write.csv(results_df, path)



# Diet sensitivity analysis -----------------------------------------------

top_twins_df <- read.csv('uniprot.csv')
top_twins <- top_twins_df$X[top_twins_df$q_value < 0.05]

## Make copy of df since diet is only available in a subset of participants
df_diet <- input_df[, c('iid', 'fid', 'age_scaled', 'BMI_scaled', 'batch', top_twins, 'dov', 'IBS_overall')]

## Load diet data
food_groups <- read.csv("ffq_foodgroups_nutrients_energyadjvar_1_2_3_3b_4.csv")[-1]

## Estblish food group vars
food_group_vars <- c(
  "animal_fat_adj", "dairy_adj", "eggs_adj", "fish_adj", "fruit_juice_adj",
  "fruit_adj", "legumes_adj", "meat_adj", "misc_animal_adj", "nuts_adj",
  "potatoes_adj", "refined_grains_adj", "sugar_beverage_adj", "sweet_dessert_adj",
  "tea_coffee_adj", "veg_oil_adj", "vegetables_adj", "whole_grain_adj",
  "margarine_adj", "alcohol_adj"
) 

## Convert dates to dates in all datasets
food_groups$FFQ_Date <- as.Date(food_groups$FFQ_Date, format = "%Y-%m-%d")
df_diet$dov <- as.Date(df_diet$dov, format = "%Y-%m-%d")


## Keep observation closest to metagenomics date per sample
df_diet_unique <- df_diet %>%
  left_join(food_groups, by = c("iid" = "Participant_ID")) %>%
  mutate(date_diff = abs(FFQ_Date - dov)) %>%
  group_by(iid) %>%
  slice_min(order_by = date_diff, n = 1, with_ties = FALSE) %>%
  ungroup()

df_full <- df_diet_unique[complete.cases(df_diet_unique[food_group_vars]), ]

## Impute missing diet data with means 
df_diet_unique <- df_diet_unique %>%
  mutate(across(all_of(food_group_vars), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

df_diet_unique[, c(food_group_vars)] <- scale(df_diet_unique[, c(food_group_vars)])

f_null <- paste("IBS_overall ~ age_scaled + BMI_scaled + batch + (1 | fid) + ", paste0(food_group_vars, collapse = '+'))

### Run next 6 lines all together
path <- "uniprot_diet_adjusted.csv"
lmm_genes(df_diet_unique, top_twins, f_null, path)
results_df <- read.csv(path)
results_df$EC_number <- mapping[results_df$X]
write.csv(results_df, path)


# Validate in PREDICT -----------------------------------------------------

#### Get the proteins which are to be replicated
top_twins_df <- read.csv('uniprot.csv')
top_twins <- top_twins_df$EC_number[top_twins_df$q_value < 0.05]


## Get mapping file for converting SAMEA IDs to iids
mapping <- read.csv('P1_map.tsv', sep = '\t')

# Add iid column
mapping$iid <- gsub('predict', '', mapping$username)

# Add "_Abundance-CPM" to the sample_accession number to match format in predict_ko
mapping$sample_accession <- paste0(mapping$sample_accession, "_Abundance-CPM")


## Load predict bacterial dataset 
predict <- read.csv('ibs_microbiome_predict_drug_processed.csv')[, -1]
predict$fid <- as.integer(substr(predict$iid, 1, nchar(predict$iid) - 1))
predict <- predict[predict$sex != 'M', ]


## Drop those not in our sample from mapping data
mapping <- mapping[mapping$iid %in% predict$iid, ]


## Load PREDICT uniprot data 
predict_up <- vroom::vroom("genefamily_uniprot_all_predict1.tsv")


## Drop those not in our sample from predict uniprot
ids_to_keep <-mapping$sample_accession[mapping$sample_accession %in% names(predict_up)]
predict_up <- predict_up[, c('uniprotid', 'EC number', 'Protein names', ids_to_keep)]


## Get minimum non-zero for later analyses 
sample_abundance <- unlist(predict_up[grepl('Abun', names(predict_up))])
one_tenth_min_non_zero <- min(sample_abundance[sample_abundance > 0])/10
rm(sample_abundance)


########
#### Same procedure as in Twins for grouping together entries with the same EC number
#### Follow procedure for twins 
########


#### Get all separate ECs as a df
ec_split <- setNames(strsplit(predict_up$`EC number`, ";"), predict_up$uniprotid)

## Get unique EC numbers and sum the same columns if any of the newly added 
## columns contain the EC number

ec_to_df <- function(ec_split) {
  # Find the maximum length of the lists
  max_length <- max(sapply(ec_split, length))
  
  # Convert the list to a data frame, filling with NA for shorter vectors
  df <- do.call(rbind, lapply(ec_split, function(x) {
    length(x) <- max_length
    x
  }))
  
  # Convert to a data frame and set the row names
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  rownames(df) <- names(ec_split)
  
  # Rename columns
  colnames(df) <- paste0("C", seq_len(ncol(df)))
  
  df <- tibble::rownames_to_column(df, var = "uniprotid")
  rownames(df) <- NULL
  
  return(df)
}
ec_df <- ec_to_df(ec_split)

predict_up <- merge(predict_up, ec_df, by = 'uniprotid')

## Sum the sample columns of all rows containing each unique EC
unique_ecs <- unique(unlist(predict_up[names(ec_df)[2:length(ec_df)]]))
unique_ecs <- unique_ecs[!is.na(unique_ecs)]
C_cols <- names(ec_df)[names(ec_df) != 'uniprotid']

rows <- list()
for (i in unique_ecs) {
  print(match(i, unique_ecs))
  # Get rows which contain ec
  one_ec_df <- predict_up %>% filter(if_any(all_of(C_cols), ~ . == i))
  
  # Select the columns containing 'Abundance' and summarize
  summed_abundance <- one_ec_df %>% 
    select(contains('Abundance')) %>% 
    summarise_all(sum, na.rm = TRUE) %>% 
    unlist()
  
  # Store the summed_abundance as a named vector with the EC number as the name
  rows[[i]] <- summed_abundance
}

# Convert the list to a dataframe and add the EC_number column
combined_proteins <- bind_rows(rows, .id = "EC_number")
combined_proteins$EC_number <- unique_ecs


# Transpose df
combined_proteins <- as.data.frame(t(combined_proteins))

# Set protein names to column names
names(combined_proteins) <- combined_proteins[1, ]

# Remove first row (redundant)
combined_proteins <- combined_proteins[-1, ]


# Save rownames (i.e., accession ID) as first column and remove redundant information
# to add to dataframe below
sample_accession <- sub("_[^_]+$", "", paste0(rownames(combined_proteins)))


# Make enzyme concentrations numeric
combined_proteins <- as.data.frame(lapply(combined_proteins[, names(combined_proteins) != 'prefix'], as.numeric))


## Check prevalence of enzymes to replicate
round(100*(colSums(combined_proteins[top_twins] == 0) / nrow(combined_proteins)),2)

#X.6.3.5.. X2.7.1.45  X2.7.7.2 X1.1.1.81 
# 5.92      0.00      0.99      8.22


## CLR transform uniprot dataset
combined_proteins <- as.data.frame(compositions::clr(combined_proteins + one_tenth_min_non_zero))


# Map existing names to pseudonames to suit model formula
enzyme_indices <- paste0('X', 1:length(top_twins))
names(top_twins) <- enzyme_indices

for (i in seq_along(top_twins)) {
  names(combined_proteins)[names(combined_proteins) == top_twins[i]] <- names(top_twins)[i] 
}


# Add IBS, batch, fid, age and BMI 
combined_proteins$sample_accession <- paste0(sample_accession, '_Abundance-CPM')
combined_proteins <- merge(mapping[, c('iid', 'sample_accession')], combined_proteins[c('sample_accession', names(top_twins))], by = 'sample_accession')
input_df <- merge(predict[, c('iid', 'IBS_overall', 'fid', 'age', 'bmi')], combined_proteins, by = 'iid')
input_df[, c("age_scaled", "BMI_scaled")] <- scale(input_df[, c("age", "bmi")])


f_null <- paste("IBS_overall ~ age_scaled + BMI_scaled + (1 | fid)")

### Run next 6 lines all together
path <- "uniprot_predict.csv"
lmm_genes(input_df, names(top_twins), f_null, path)
results_df <- read.csv(path)
results_df$EC_number <- top_twins[results_df$X]
write.csv(results_df, path)



























