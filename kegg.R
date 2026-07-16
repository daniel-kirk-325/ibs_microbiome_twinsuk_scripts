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

  # Add q-values to the DataFrame
  results_df_ibs$q_value <- q_values

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


# KEGG pathways  -------------------------------------------------------

KEGG <- read.csv("strata_removed_KEGG_pathways-cpm_20211201.tsv", sep = '\t')

## Transpose KEGG
KEGG_t <- as.data.frame(t(KEGG))

# First row to col names & remove first row
names(KEGG_t) <- KEGG_t[1,]
KEGG_t <- KEGG_t[-1, ]

# Rownames to first col and reset rownames
KEGG_t <- cbind(rownames(KEGG_t), KEGG_t) ; rownames(KEGG_t) <- NULL
names(KEGG_t)[names(KEGG_t) == names(KEGG_t)[1]] <- 'prefix'

# Update prefix
KEGG_t$prefix <- gsub(".*X(.*?)_Ab.*", "\\1", KEGG_t$prefix)


# Keep only those in our dataset
KEGG_t <- KEGG_t[KEGG_t$prefix %in% df$prefix, ]

# Convert to numeric 
KEGG_names <- names(KEGG_t)[4:length(KEGG_t)]
KEGG_t[, 4:length(KEGG_t)] <- lapply(KEGG_t[, 4:length(KEGG_t)], as.numeric)


KEGG_to_ko <- read.csv("KEGG_modules_to_KO (1).tsv", sep = '\t')


## Kegg to ko -> Add abundance of KOs in each module
# Column names should be first column
KEGG_to_ko <- rbind(names(KEGG_to_ko), KEGG_to_ko); names(KEGG_to_ko) <- c('Module', 2:length(KEGG_to_ko))

# Convert dataframe to 1D vector
KEGG_to_ko_vector <- as.vector(t(KEGG_to_ko))
KEGG_to_ko_vector <- KEGG_to_ko_vector[KEGG_to_ko_vector != '']

# Get index locations of Module
Ms <- grep("^M", KEGG_to_ko_vector)

# Module contents are index of each M until index of next M
module_sums <- list()
for (i in seq_along(Ms)) {
  if (i != length(Ms)) {
    KOs <- KEGG_to_ko_vector[(Ms[i]+1) : (Ms[i+1]-1)]
  } else {
    KOs <- KEGG_to_ko_vector[(Ms[i]+1) : length(KEGG_to_ko_vector)]
  }
  KOs <- KOs[KOs %in% names(KEGG_t)]
  mod_sums <- rowSums(KEGG_t[KOs])
  module_sums[[KEGG_to_ko_vector[Ms[i]]]] <- c(mod_sums)
}

modules_df <- data.frame(module_sums)
modules_df_raw <- modules_df


## Save low prevalence modules
low_prevalence <- names(modules_df)[(colSums(modules_df == 0)) > 0.9 * nrow(modules_df)]

# Establish vector for running model below
modules <- names(modules_df)

# CLR transform modules
# Use min-non zero of original df, not the updated modules df
min_non_zero <- min(KEGG_t[, 4:length(KEGG_t)][KEGG_t[, 4:length(KEGG_t)] > 0])/10
modules_df[modules] <- compositions::clr(modules_df[modules] + min_non_zero)

# Bind to input df containing covs and IBS status
modules_df <- cbind(KEGG_t$prefix, modules_df)
modules_df <- modules_df %>% dplyr::rename(prefix = "KEGG_t$prefix")
modules_df <- merge(df, modules_df, by = 'prefix')

# Drop low prevalence and low RA
modules <- modules[!modules %in% low_prevalence]



#### Remove modules with identical (or near identical) information

sp <- cor(modules_df[modules], method = 'spearman', use = "pairwise.complete.obs")

# Convert correlation matrix to long format
sp_long <- as.data.frame(as.table(sp))
colnames(sp_long) <- c("Var1", "Var2", "Correlation")

# Remove self-correlations and duplicate pairs
sp_long_unique <- sp_long[as.character(sp_long$Var1) < as.character(sp_long$Var2), ]

# Sort by absolute correlation
sp_sorted <- sp_long_unique[order(-abs(sp_long_unique$Correlation)), ]

# Get extreme correlations
extrm <- sp_sorted[abs(sp_sorted$Correlation) > 0.999, ]

# Remove one from each pair
sum(extrm$Var1 %in% extrm$Var2) # Equals 0; means one column can be removed at random
modules <- modules[!modules %in% extrm$Var2]
 


#### Do univariate analyses
modules_df[, c("age_scaled", "BMI_scaled")] <- scale(modules_df[, c("age_Mb_sample", "bmi_Mb_sample")])
f_null <- paste("IBS_overall ~ age_scaled + BMI_scaled + batch + (1 | fid)")
#lmm_genes(modules_df, modules, f_null, 'kegg.csv')



# Drug sensitivity analysis -----------------------------------------------

top_twins_df <- read.csv('kegg.csv')
top_twins <- top_twins_df$X[top_twins_df$q_value < 0.05]

drug_vars <- c("Proton.Pump.Inhibitors", "Anti_Depressants.anxiolytics", "Antispasmodics.antimotility", "Analgesics", 'Laxatives')
f_null <- paste("IBS_overall ~ age_scaled + BMI_scaled + batch + (1 | fid) + ", paste0(drug_vars, collapse = '+'))
#lmm_genes(modules_df, top_twins, f_null, 'kegg_drug_adjusted.csv')




# Diet sensitivity analysis -----------------------------------------------

top_twins_df <- read.csv('kegg.csv')
top_twins <- top_twins_df$X[top_twins_df$q_value < 0.05]

## Make copy of df since diet is only available in a subset of participants
df_diet <- modules_df

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

## Scale age and BMI
df_diet_unique[, c(food_group_vars)] <- scale(df_diet_unique[, c(food_group_vars)])
f_null <- paste("IBS_overall ~ age_scaled + BMI_scaled + batch + (1 | fid) + ", paste0(food_group_vars, collapse = '+'))
#lmm_genes(df_diet_unique, top_twins, f_null, 'kegg_diet_adjusted.csv')




# Validate in PREDICT -----------------------------------------------------

## Get mapping file for converting SAMEA IDs to iids
mapping <- read.csv('P1_map.tsv', sep = '\t')

# Add iid column
mapping$iid <- gsub('predict', '', mapping$username)

# Add "_Abundance-CPM" to the sample_accession number to match format in predict_ko
mapping$sample_accession <- paste0(mapping$sample_accession, "_Abundance-CPM")


## Load predict dataset 
predict <- read.csv('ibs_microbiome_predict_drug_processed.csv')[, -1]
predict$fid <- as.integer(substr(predict$iid, 1, nchar(predict$iid) - 1))
predict <- predict[predict$sex != 'M', ]


## Load predict KEGG dataset
predict_ko <- vroom::vroom("PREDICT1_humann_genefamilies_KO_uniref50_cpm.tsv")

# Transpose df
predict_ko <- as.data.frame(t(predict_ko))
# Set first row to column names and remove first row
names(predict_ko) <- predict_ko[1, ]
predict_ko <- predict_ko[-1, ]
# Set rownames to first column 
predict_ko <- cbind(rownames(predict_ko), predict_ko) ; rownames(predict_ko) <- NULL
# Rename first column
names(predict_ko)[names(predict_ko) == names(predict_ko)[1]] <- 'sample_accession'


# Drop those not in our dataset
mapping <- mapping[mapping$iid %in% predict$iid, ]

# Add predict iid to predict_ko
predict_ko <- merge(mapping, predict_ko, by = 'sample_accession')


## Keep only community-level KEGG pathways 

ko_columns <- names(predict_ko %>% select(matches("^K\\d+")))
# Remove bacterial species information
ko_only <- ko_columns[!grepl('\\|', ko_columns)]


#### Merge KEGG pathways from the same module

KEGG_to_ko <- read.csv("KEGG_modules_to_KO (1).tsv", sep = '\t')

## Kegg to ko -> Add abundance of KOs in each module
# Column names should be first column
KEGG_to_ko <- rbind(names(KEGG_to_ko), KEGG_to_ko); names(KEGG_to_ko) <- c('Module', 2:length(KEGG_to_ko))

# Convert dataframe to 1D vector
KEGG_to_ko_vector <- as.vector(t(KEGG_to_ko))
KEGG_to_ko_vector <- KEGG_to_ko_vector[KEGG_to_ko_vector != '']

# Get index locations of Module
Ms <- grep("^M", KEGG_to_ko_vector)

# Convert relevant columns to numeric 
predict_ko[ko_only] <- lapply(predict_ko[ko_only], as.numeric)

  
# Module contents are index of each M until index of next M
module_sums <- list()
for (i in seq_along(Ms)) {
  if (i != length(Ms)) {
    KOs <- KEGG_to_ko_vector[(Ms[i]+1) : (Ms[i+1]-1)]
  } else {
    KOs <- KEGG_to_ko_vector[(Ms[i]+1) : length(KEGG_to_ko_vector)]
  }
  KOs <- KOs[KOs %in% ko_only]
  mod_sums <- rowSums(predict_ko[KOs])
  module_sums[[KEGG_to_ko_vector[Ms[i]]]] <- c(mod_sums)
}

modules_df <- data.frame(module_sums)
modules_df_raw <- modules_df

modules <- names(modules_df)

# CLR transform modules
# Use min-non zero of original df, not the updated modules df
min_non_zero <- min(predict_ko[, ko_only][predict_ko[, ko_only] > 0])/10
modules_df[modules] <- compositions::clr(modules_df[modules] + min_non_zero)

# Bind to input df containing covs and IBS status
modules_df <- cbind(predict_ko$iid, modules_df)
modules_df <- modules_df %>% rename(iid = "predict_ko$iid")
modules_df <- merge(predict[, c('iid', 'fid', 'IBS_overall', 'age', 'bmi')], modules_df, by = 'iid')


# Get Modules to replicate 
top_twins_df <- read.csv('kegg.csv')
top_twins <- top_twins_df$X[top_twins_df$q_value < 0.05]


#### Run LMM
modules_df[, c("age_scaled", "BMI_scaled")] <- scale(modules_df[, c("age", "bmi")])
f_null <- paste("IBS_overall ~ age_scaled + BMI_scaled + (1 | fid)")
#lmm_genes(modules_df, top_twins, f_null, 'kegg_predict.csv')























