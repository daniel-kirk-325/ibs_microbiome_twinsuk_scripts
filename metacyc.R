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


## Load functional data
df <- read.csv("ibs_microbiome_twins_drug_processed.csv")

## Drop males
df <- df[df$SEX != 'M', ]

## Keep only those with phage data
phage_ids <- read.csv("viroprofiler_iids.csv")[, -1]
df <- df[df$iid %in% phage_ids, ]

## Estbalish vector with genes   

## Add family ID
df$fid <- as.integer(substr(df$iid, 1, nchar(df$iid) - 1))

## Load metacyc data
metacyc_labels <- read.csv("MetaCyc_pathway_map_update.tsv", sep = '\t')
metacyc_data <- read.csv("strata_removed_pathabundance-cpm_20211201.tsv", sep = '\t')


# Pathways  -----------------------------------------------------------------

# Save pathways for later steps
pathways <- metacyc_data$Pathway[3:nrow(metacyc_data)]

metacyc_data <- as.data.frame(t(metacyc_data)) %>%
  rownames_to_column(var = "prefix") %>% 
  setNames(.[1, ]) %>%
  dplyr::slice(-1) %>% 
  dplyr::rename(prefix = Pathway)

metacyc_data$prefix <- substr(metacyc_data$prefix, 2, (nchar(metacyc_data$prefix) - nchar('_Abundance.CPM')))

input_df <- merge(df, metacyc_data, by = 'prefix')

input_df[pathways] <- lapply(input_df[pathways], as.numeric)
# Save copy of df before CLR scaling
input_df_original_scale <- input_df

## Note pathways present in >10% of samples
pathways_to_keep <- names(input_df[pathways])[colSums(input_df[pathways] > 0) >= 0.1*nrow(input_df)] 

## CLR-transform pathways
min_non_zero <- min(input_df[pathways][input_df[pathways] > 0])/10
input_df[pathways] <- compositions::clr(input_df[pathways] + min_non_zero)

## Keep only pathways present in >10%
pathways <- pathways_to_keep

## Recode nams for function
pathway_indices <- paste0('X', 1:length(pathways))
colnames(input_df)[colnames(input_df) %in% pathways] <- pathway_indices
mapping <- setNames(pathways, pathway_indices)


## Run LMM
input_df[, c("age_scaled", "BMI_scaled")] <- scale(input_df[, c("age_Mb_sample", "bmi_Mb_sample")])
f_null <- paste("IBS_overall ~ age_scaled + BMI_scaled + batch + (1 | fid)")

lmm_genes(input_df, pathway_indices, f_null, "metacyc.csv")
results_df <- read.csv("metacyc.csv")
results_df$X <- mapping[results_df$X]
results_df <- results_df %>% dplyr::rename(pathway = X)
write.csv(results_df, "metacyc.csv")




