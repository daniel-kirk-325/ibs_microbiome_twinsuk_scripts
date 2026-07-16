require(lme4)
library(dplyr)
library(compositions)
library(ggplot2)

## --------------------- --------------------- ---------------------
## Calculates correlation matrix between all filtered bacteria and VCs
## in IBS vs controls and then compares if these differ via Frobenius
## norm

## This reveals whether global correlations between phages and bacteria 
## differ in IBS
## --------------------- --------------------- ---------------------

## Load data
df <- read.csv("ibs_microbiome_twins_drug_processed.csv")[, -1]

## Drop males
df <- df[df$SEX != 'M', ]

## Add family ID
df$fid <- as.integer(substr(df$iid, 1, nchar(df$iid) - 1))

## Set batch as factor 
df$batch <- as.factor(df$batch)

## Establish bacteria 
bacteria <- names(df)[!(names(df) %in% c("ParticipantID", "ResponseDate", "ibs_status", "Origin", "IBS_overall", "prevalent_ibs", "incident_ibs", "age_incident_ibs", "prefix", 'iid', 'dov', 
                                         'SEX', 'age_Mb_sample', 'bmi_Mb_sample', 'YEAR_BIRTH', 'fid', 'sampling_diff', 'age', 'ibs_status_upd', 'batch', 
                                         "Proton.Pump.Inhibitors", "Anti_Depressants.anxiolytics", "Antispasmodics.antimotility", "Analgesics", 'Laxatives', 
                                         'subtype'))]


# Load phage data ------------------------------------------------------------------

## Add phage data
phages <- read.csv("merged_tpm_COMBINED_UNIQUE.csv")

phages_vigs <- names(phages)[grep("^phage_[0-9]+", names(phages))]
print("LENGTH OF PHAGES VIGS:")
print(length(phages_vigs))

## Add GSID key to dataframe
GSID_key <- read.csv("datasets/raw_microbiome_data/list of ids with metagenomes and date of visit.txt", sep = '\t')
df <- merge(GSID_key[, c('prefix', 'gsid')], df, by = 'prefix')

# Replace dashed with underscores in df$gsid to match format in phage data
df$gsid <- gsub('-', '_', df$gsid)

## Merge phages based on GSID
df <- merge(df, phages[!names(phages) %in% c("X", "prefix", "iid", "sex", "dov", "batch", 'Batch', "col_name")], by = 'gsid')



# Apply prevalence & relative abundance filters ---------------------------

## Identify low abundance bacteria 
low_prev <- names(df[bacteria][colSums(df[bacteria] == 0) > 0.9 * nrow(df)])

# CLR transformation of bacteria 
one_tenth_min_non_zero <- min(df[bacteria][df[bacteria] > 0])/10
df[bacteria] <- clr(df[bacteria] + one_tenth_min_non_zero)

# Drop low prevalent bacteria from bacteria vector
bacteria <- bacteria[!bacteria %in% low_prev]


## Identify low abundance phages
phage_prev <- colSums(df[phages_vigs] == 0)
to_drop <- names(phage_prev[phage_prev > 0.9 * nrow(df)])


## Identify phages with RA below threshold
ra_threshold <- 1e-3
phage_ra <- colMeans(df[phages_vigs[!phages_vigs %in% to_drop]])
low_ra <- names(phage_ra[phage_ra < ra_threshold])


## Combine low abundant and low relative abundance phages
to_drop <- unique(c(to_drop, low_ra))

# CLR tranformation of phages
min_nonzero <- min(unlist(df[phages_vigs])[unlist(df[phages_vigs]) > 0])
df[phages_vigs] <- clr(df[phages_vigs] + (min_nonzero / 10))


## Drop low prev phages from df
df <- df[, !(colnames(df) %in% to_drop)]

## Update phages_vigs to exclude dropped columns
phages_vigs <- setdiff(phages_vigs, to_drop)

print("Dim bacteria:")
print(length(bacteria))

print("Dim phages_vigs:")
print(length(phages_vigs))

ibs_df <- df[df$IBS_overall == 1, ]
hc_df <- df[df$IBS_overall == 0, ]



# Calculate observed correlation matrices ---------------------------------


## Correlations in IBS
ibs_cor <- cor(ibs_df[phages_vigs], ibs_df[bacteria], method = 'spearman')

#write.csv(ibs_cor, 'ibs_cor_matrix.csv')
print("ibs_cor done")

## Correlations in HC
hc_cor <- cor(hc_df[phages_vigs], hc_df[bacteria], method = 'spearman')

#write.csv(hc_cor, 'hc_cor_matrix.csv')
print("hc_cor done")

## Calculate global distance measure
obs_stat <- sum((hc_cor - ibs_cor)^2)
print(paste0('Obsered sum of differences: ', obs_stat))



# Shuffle labels and calculate permuted correlation matrices --------------

#### Maintain proportions of IBS and relatedness when shuffling

## Identify singlets
singlets <- names(table(df$fid)[table(df$fid) == 1])
singlets <- df[df$fid %in% singlets, ]


## Identify concordant and discordant twins
related = df[!df$fid %in% singlets$fid, ]
concordant <- related %>% group_by(fid) %>% filter(n_distinct(IBS_overall) == 1) %>% ungroup()
discordant <- related %>% group_by(fid) %>% filter(n_distinct(IBS_overall) > 1) %>% ungroup()


## Establish length of healthy cases in singlets and concordant twins
healthy_singlets <- table(singlets$IBS_overall)[1]
healthy_conc <- as.integer(table(concordant$IBS_overall)[1] / 2)

seeds <- c()
greater_than_observed <- c()
permuted_squared_differences <- c()
n_sims = 5000

for (seed in 1:n_sims) {
  set.seed(seed)
  print(seed)

  ## Permuate labels of singlets, discordant and concordant individually, 
  ## then merge for permutation tests
  
  # Establish S1 singlets
  S1_singlets <- sample(singlets$iid, size = healthy_singlets, replace = FALSE)
  
  # Establish S1 for concordant twins
  S1_conc <- sample(unique(concordant$fid), size = healthy_conc, replace = FALSE)

  # Randomly select one iid from each discordant pair to be in S1
  S1_disc <- discordant %>%
    group_by(fid) %>%
    slice_sample(n = 1, replace = FALSE) %>%
    pull(iid)
  
  S1 <- df[(df$iid %in% c(S1_singlets,S1_disc)) | (df$fid %in% S1_conc), ]
  S2 <- df[!df$iid %in% S1$iid, ]
  
  S1_cor <- cor(S1[phages_vigs], S1[bacteria], method = 'spearman')
  S2_cor <- cor(S2[phages_vigs], S2[bacteria], method = 'spearman')
  
  seeds <- c(seeds, seed)
  permuted_squared_differences <- c(permuted_squared_differences, (sum((S1_cor - S2_cor)^2)))
  greater_than_observed <- c(greater_than_observed, (sum((S1_cor - S2_cor)^2) > obs_stat))
}

stat_diffs <- data.frame(seed = seeds, permuted_squared_differences = permuted_squared_differences, greater_than_observed = greater_than_observed)


# Calculate p-value  
p_perm <- signif(sum(stat_diffs$greater_than_observed/nrow(stat_diffs)),2)

## p=0.002














