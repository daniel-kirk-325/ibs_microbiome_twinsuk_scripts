library(dplyr)      # 1.2.1
library(SpiecEasi)  # 1.1.2
library(NetCoMi)    # 1.2.0


# CODE EXPLANATION  -------------------------------------------------------

## Performs network analysis on the phage and bacterial counts data after applying
## filtering and adjusting for covariates.

## The script was run in CREATE HPC. Relevant outputs and the network itself are 
## then saved to and results are generated locally (see network_results_compilation.R) 

# Load bacterial data ---------------------------------------------------------------

#### Load bacterial relative abundance data 
bacterial_ra <- read.table("merged_metaphlan_bugs_abundance_table_20211201_version_4.tsv", sep = '\t')

print("Bacteria loaded")

# Make first row column names 
names(bacterial_ra) <- bacterial_ra[1, ]

# Keep species only
bacterial_ra <- bacterial_ra[(grepl("s__", bacterial_ra$clade_name)) & (!grepl("t__", bacterial_ra$clade_name)), ]

# Format taxa name so that only species remains 
bacterial_ra$clade_name <- sub(".*s__", "", bacterial_ra$clade_name)


#### Keep samples used in main analysis (determined by prefix)

# 1) load IDs with phage data
phage_ids <- read.csv("viroprofiler_iids.csv")[, -1]

# 2) load dataset used in bacterial analysis
bacteria_ids <- read.csv(".ibs_microbiome_twins_drug_processed.csv")[, -1]

# 3) Remove males 
bacteria_ids <- bacteria_ids[bacteria_ids$SEX != 'M', ]

# 4) Drop those without phage data
bacteria_ids <- bacteria_ids[bacteria_ids$iid %in% phage_ids, ]

# 5) Keep samples based on prefix
to_keep <- names(bacterial_ra)[names(bacterial_ra) %in% bacteria_ids$prefix]
bacterial_ra <- bacterial_ra[, c('clade_name', to_keep)]

# 6) Transpose dataframe so that bacteria are columns and samples are rows
bacterial_ra <- as.data.frame(t(bacterial_ra), stringsAsFactors = FALSE)
colnames(bacterial_ra) <- bacterial_ra[1, ]                              # use first row as column names
bacterial_ra <- bacterial_ra[-1, ]                                       # remove first row
bacterial_ra <- cbind(Prefix = rownames(bacterial_ra), bacterial_ra)     # move rownames into a column
rownames(bacterial_ra) <- NULL

print("Bacteria transposed")

#### Transform relative abundance to counts by multiplying each value by reads_QC per sample

## Load metadata containing QCed reads
reads_data <- read.csv("metadata_release_20211201_PL.csv")

## Add this to bacterial dataframe
bacterial_ra <- merge(reads_data[, c('Prefix', 'Mreads_QC', 'GSID')], bacterial_ra, by = 'Prefix')

## Multiply bacterial relative abundance by Mreads_QC per sample
bacteria <- names(bacterial_ra)[!names(bacterial_ra) %in% c('Prefix', 'Mreads_QC', 'GSID')]
bacterial_ra[bacteria] <- sapply(bacterial_ra[bacteria], as.numeric)
bacterial_ra[bacteria] <- bacterial_ra[bacteria] * bacterial_ra$Mreads_QC

## Rename df for clarity 
bacterial_counts <- bacterial_ra ; rm(bacterial_ra)

print("Bacteria converted to counts")

# Add phage data ----------------------------------------------------------

## Add phage data 
phages <- read.csv("merged_tpm_COMBINED_UNIQUE.csv")

print("Phages loaded")

## Replace _ with - for comparability with gsids in reads_data
phages$gsid <- gsub('_', '-', phages$gsid)

## Keep only samples used in main analysis
phages <- phages[phages$gsid %in% bacterial_counts$GSID, ]

## Establish vector containing phage columns
phages_vigs <- names(phages)[grepl('phage_[0-9]', names(phages))]

## Merge bacteria with phage data based on GSID
df <- merge(bacterial_counts, phages[,c('gsid', phages_vigs)], by.x = 'GSID', by.y = 'gsid')
df[phages_vigs] <- lapply(df[phages_vigs], as.numeric)

print("Phages merged with bacteria")

# Set rownames to prefix (used to establish cases and controls below)
rownames(df) <- df$Prefix 

# Sort df by row names
df <- df[order(rownames(df)), ]

## Multiply phage columns by Mreads_QC to obtain counts
df[phages_vigs] <- df[phages_vigs] * as.numeric(df$Mreads_QC)

print("Phage relative abundance converted to counts")

## Drop unnecessary columns to simplify downstream operations
df <- df[, -which(names(df) %in% c("GSID", 'Prefix', 'Mreads_QC'))]


# Data pre-processing ----------------------------------------------------------

## Establish cases and controls 
case_ids <- bacteria_ids$prefix[bacteria_ids$IBS_overall == 1]
control_ids <- bacteria_ids$prefix[bacteria_ids$IBS_overall == 0]

cases <- df[case_ids, ] 
controls <- df[control_ids, ]

## Save low prevalent (missing in >80% of samples) VCs and bacteria in both groups
low_prev_cases <- names(cases)[colSums(cases == 0)  > 0.8*nrow(cases)]
low_prev_controls <- names(controls)[colSums(controls == 0)  > 0.8*nrow(controls)]
low_prev <- unique(low_prev_cases, low_prev_controls)


#### Apply additional relative abundance threshold to phages

## Get GSID of cases and controls (since phage relative abundance data does not contain prefix)
case_gsids <- bacterial_counts$GSID[bacteria_ids$prefix %in% case_ids]
control_gsids <- bacterial_counts$GSID[bacteria_ids$prefix %in% control_ids]

# Set relative abundance threshold
ra_threshold <- 1e-2 # (0.01%)

# Calculate mean relative abundance of phages in cases and controls
case_ras <- colMeans(phages[phages$gsid %in% case_gsids, phages_vigs])         
control_ras <- colMeans(phages[phages$gsid %in% control_gsids, phages_vigs])   

# Get phages with a low relative abundance 
low_ra <- c(names(case_ras[case_ras < ra_threshold]), names(control_ras[control_ras < ra_threshold]))

# Add these to the low prevalence features, giving final list of features to drop from analysis
features_to_drop <- unique(c(low_prev, low_ra))
filtered_species <- names(df)[!names(df) %in% features_to_drop]

print("Species filtering done. Numner of remaining species:")
print(length(filtered_species))

# Adjust batch and covariates ---------------------------------------------

# 1) CLR-transform the counts data (with pseudocount)
df[] <- compositions::clr(df + 1)

print("Counts data CLR-transformed")

# 2) Regress against age, BMI and batch on features which passed low prevalence
#    and low abundance filtering

# 2.1)  Get batch, age and BMI and from data used in univariate analysis
metadata <- bacteria_ids[c('prefix', 'age_Mb_sample', 'bmi_Mb_sample', 'batch')]

# Sort by prefix to match order in df
metadata <- metadata[order(metadata$prefix), ]

# Check samples in metadata align with order in df
all(rownames(df) == metadata$prefix) # TRUE

# Ensure batch is treated as a factor
metadata$batch <- as.factor(metadata$batch)


# 2.2) Make empty matrix in which to add residuals
adjusted_df <- matrix(NA, nrow = nrow(df), ncol = length(filtered_species))
rownames(adjusted_df) <- rownames(df)
colnames(adjusted_df) <- filtered_species

for (j in filtered_species) {
  fit <- lm(df[, j] ~ age_Mb_sample + bmi_Mb_sample + batch, data = metadata)
  adjusted_df[, j] <- residuals(fit)
}

print("CLR-transformed counts data residualised")

# Create and compare networks ---------------------------------------------

# Re-establish cases and controls with adjusted datasets & dropping
# low prevlance species (and low relative abundance VCs)

cases <- adjusted_df[case_ids, ]
controls <- adjusted_df[control_ids, ]

print("Running network")

## Make combined network
se_mb_est  <- netConstruct(data = cases,
                           data2 = controls,
                           measure = "spieceasi",
                           measurePar = list(method = "mb",                       # ~ mb reported to have superior performance in Kurtz et al. 2015
                                             lambda.min.ratio = 1e-2,
                                             nlambda = 25,                        # Regularization,
					     pulsar.params = list(rep.num = 70)), # Stability 
                           normMethod = "none",  # Dealt with internally by spieceasi (Kurtz 2015)
                           zeroMethod = "none",  # Dealt with internally by spieceasi (Kurtz 2015)
                           sparsMethod = "none", # Dealt with internally by spieceasi (https://microbiome.github.io/OMA/docs/devel/pages/network_learning.html#sec-spieceasi-mb)
                           verbose = 3,
                           seed = 123)

print("Network done")

print("Check stability (should be near 0.05)")

print("Stability (cases):")
tryCatch({
  print(SpiecEasi::getStability(se_mb_est$measureOut1))
}, error = function(e) {
  message("Error computing stability (cases): ", conditionMessage(e))
})

print("Stability (controls):")
tryCatch({
  print(SpiecEasi::getStability(se_mb_est$measureOut2))
}, error = function(e) {
  message("Error computing stability (controls): ", conditionMessage(e))
})

#### Remove VC-VC correlations

# Get index locs of viruses
virus_idx <- which(grepl("phage", colnames(cases)))[1] : length(colnames(cases))

## Update adjacency matrices

adja1 <- se_mb_est$adjaMat1   # IBS
adja2 <- se_mb_est$adjaMat2   # Controls

# Zero out the virus-virus block in both networks and add updated matrices to network object
adja1[virus_idx, virus_idx] <- 0
adja2[virus_idx, virus_idx] <- 0

se_mb_est$adjaMat1 <- adja1
se_mb_est$adjaMat2 <- adja2

## Repeat for association matrices

asso1 <- se_mb_est$assoMat1 # IBS
asso2 <- se_mb_est$assoMat2 # Control

asso1[virus_idx, virus_idx] <- 0
asso2[virus_idx, virus_idx] <- 0

se_mb_est$assoMat1 <- asso1
se_mb_est$assoMat2 <- asso2


## Repeat for similarity matrices

if (!is.null(se_mb_est$dissMat1)) {
  se_mb_est$dissMat1[virus_idx, virus_idx] <- Inf   # Infinite distance = no connection
  se_mb_est$dissMat2[virus_idx, virus_idx] <- Inf
}

if (!is.null(se_mb_est$simMat1)) {
  se_mb_est$simMat1[virus_idx, virus_idx] <- 0
  se_mb_est$simMat2[virus_idx, virus_idx] <- 0
}


print("VC-VC connections zeroed-out across all network objects")

## Save adjacency and association matrices
write.csv(adja1, "results/original_data/adja1_ibs.csv", row.names = F)
write.csv(adja2, "results/original_data/adja2_ctrl.csv", row.names = F)
write.csv(asso1, "results/original_data/asso1_ibs.csv", row.names = F)
write.csv(asso2, "results/original_data/asso2_ctrl.csv", row.names = F)


# Clean the edge lists 
if (!is.null(se_mb_est$edgelist1)) {
  se_mb_est$edgelist1 <- se_mb_est$edgelist1[!(grepl('phage', se_mb_est$edgelist1$v1) & 
                                                 (grepl('phage', se_mb_est$edgelist1$v2))), ]
  
  se_mb_est$edgelist2 <- se_mb_est$edgelist2[!(grepl('phage', se_mb_est$edgelist2$v1) & 
                                                 (grepl('phage', se_mb_est$edgelist2$v2))), ]
}



#### Do network comparison on filtered networks

se_net <- netAnalyze(se_mb_est,
                     clustMethod = "cluster_fast_greedy",
                     hubPar = "eigenvector",
		     hubQuant = 0.98, # Hub quantile = 2%
                     normDeg = FALSE)

# Save net 
saveRDS(se_net, file = "results/original_data/se_net.rds")

net_summary <- summary(se_net, groupNames = c("IBS", "Controls"))
net_stats <- as.data.frame(net_summary$glob_probs_lcc)
net_stats <- cbind(Statistic = rownames(net_stats), net_stats)
rownames(net_stats) <- NULL

print("Net comparison done")

print(net_summary)


#### Add phage-specific statistics

## Get association matrices
assom_ibs <- as.data.frame(se_net$input$assoMat1)
assom_ctl <- as.data.frame(se_net$input$assoMat2)

## Convert to long format (upper triangle)
assom_ibs <- assom_ibs %>%
  tibble::rownames_to_column("Var1") %>%
  tidyr::pivot_longer(cols = -Var1, names_to = "Var2", values_to = "value") %>%
  dplyr::filter(Var1 < Var2)

assom_ctl <- assom_ctl %>%
  tibble::rownames_to_column("Var1") %>%
  tidyr::pivot_longer(cols = -Var1, names_to = "Var2", values_to = "value") %>%
  dplyr::filter(Var1 < Var2)


## Keep non-zero phage correlations only
non_zero_ibs <- assom_ibs[(grepl("phage", assom_ibs$Var1) | grepl("phage", assom_ibs$Var2)) & abs(assom_ibs$value) != 0, ]
non_zero_ctl <- assom_ctl[(grepl("phage", assom_ctl$Var1) | grepl("phage", assom_ctl$Var2)) & abs(assom_ctl$value) != 0, ]


## How many non-zero phage correlations in total in each network?
net_stats <- rbind(net_stats, c('N non-zero VC correlations', nrow(non_zero_ibs), nrow(non_zero_ctl)))


## Direction of associations
neg_ibs <- nrow(non_zero_ibs[non_zero_ibs$value < 0, ])
neg_ctl <- nrow(non_zero_ctl[non_zero_ctl$value < 0, ])

net_stats <- rbind(net_stats, c('N negative VC correlations', neg_ibs, neg_ctl))


## Save network statistics
write.csv(net_stats, "results/original_data/network_statistics.csv")

## Save summary output
writeLines(
  capture.output(summary(se_net, groupNames = c("IBS", "Controls"))),
  "results/original_data/network_summary.txt"
)


