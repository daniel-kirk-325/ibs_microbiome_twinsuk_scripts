library(dplyr)      # 1.2.1
library(SpiecEasi)  # 1.1.2
library(NetCoMi)    # 1.2.0


# CODE EXPLANATION  -------------------------------------------------------

## Performs network analysis on the phage and bacterial counts data after applying
## filtering, adjusting for covariates and shuffling IBS labels.

## The script was run in as an array job in CREATE HPC (see b_perm.sh). Relevant outputs and the network itself are 
## then saved to and results are generated locally (see network_results_compilation.R and 
## net_figure.R) 


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
bacteria_ids <- read.csv("ibs_microbiome_twins_drug_processed.csv")[, -1]

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

print("Species filtering done")

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

# Permute labels ----------------------------------------------------------

# Permute labels on the group level and ensure that proportions of labels in 
# concordant, discordant and singlets remain constant

I <- as.numeric(Sys.getenv("I")) # Load array number from bash script (bash.sh)
set.seed(I)

## 1) Establish concordancy lists

# Get singlets (fids with only one twin)
singlets <- bacteria_ids %>%
  group_by(fid) %>%
  filter(n() == 1) %>%
  ungroup()

# Get family members concordant for IBS status 
concordant <- bacteria_ids %>%
  filter(!iid %in% singlets$iid) %>% 
  group_by(fid) %>%
  filter(n_distinct(IBS_overall) == 1) %>%
  ungroup()

# Get family members discordant for IBS status 
disconcordant <- bacteria_ids %>%
  filter(!iid %in% singlets$iid) %>% 
  group_by(fid) %>%
  filter(n_distinct(IBS_overall) > 1) %>%
  ungroup()


## 2) Establish length of healthy cases in singlets and concordant twins
healthy_singlets <- sum(singlets$IBS_overall == 0)
healthy_conc <- sum(concordant$IBS_overall == 0) / 2   # each concordant pair contributes 2 rows


## 3) Create permuted datasets with same length and number of concordant 
## & discordant twins as original dataset

# Establish S1 for singlets
S1_singlets <- singlets %>%
  slice_sample(n = healthy_singlets, replace = FALSE) %>%
  pull(prefix)

# Establish S1 for concordant twins 
family_ids <- concordant %>%
  distinct(fid) %>%
  slice_sample(n = healthy_conc, replace = FALSE) %>%
  pull(fid)

S1_conc <- concordant %>%
  filter(fid %in% family_ids) %>%
  pull(prefix)

# Randomly select one iid from each discordant pair to be in S1
S1_disc <- disconcordant %>%
  group_by(fid) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  pull(prefix)

# Combine all selected iids for S1
S1 <- c(S1_singlets, S1_conc, S1_disc)
S2 <- setdiff(bacteria_ids$prefix, S1)  

## 4) Use shuffled data to run analyses

S1_mat <- adjusted_df[S1, ]
S2_mat <- adjusted_df[S2, ]

print("Datasets permuted")

print("Running network")

## 5) Make combined network using shuffled data

se_mb_est  <- netConstruct(data = S2_mat,
                           data2 = S1_mat,
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

#### Remove VC-VC correlations

# Get index locs of viruses
virus_idx <- which(grepl("phage", colnames(S1_mat)))[1] : length(colnames(S1_mat))

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

# Clean the edge lists 
if (!is.null(se_mb_est$edgelist1)) {
  se_mb_est$edgelist1 <- se_mb_est$edgelist1[!(grepl('phage', se_mb_est$edgelist1$v1) & 
                                                 (grepl('phage', se_mb_est$edgelist1$v2))), ]
  
  se_mb_est$edgelist2 <- se_mb_est$edgelist2[!(grepl('phage', se_mb_est$edgelist2$v1) & 
                                                 (grepl('phage', se_mb_est$edgelist2$v2))), ]
}


#### Do network comparison on filtered networks
print("Comparing networks")
se_net <- netAnalyze(se_mb_est,
                     clustMethod = "cluster_fast_greedy",
                     hubPar = "eigenvector",
		     hubQuant = 0.98, # Hub quantile = 2%
                     normDeg = FALSE, 
                     gcmHeat = FALSE)

net_summary <- summary(se_net, groupNames = c("S2", "S1"))
net_stats <- as.data.frame(net_summary$glob_probs_lcc)
net_stats <- cbind(Statistic = rownames(net_stats), net_stats)
rownames(net_stats) <- NULL

print("Network statistics obtained")


#### Get phage-specific statistics

## Get adjacency matrices
assom_S2 <- as.data.frame(se_net$input$assoMat1)
assom_S1 <- as.data.frame(se_net$input$assoMat2)

## Convert to long format (upper triangle)
assom_S2 <- assom_S2 %>%
  tibble::rownames_to_column("Var1") %>%
  tidyr::pivot_longer(cols = -Var1, names_to = "Var2", values_to = "value") %>%
  dplyr::filter(Var1 < Var2)

assom_S1 <- assom_S1 %>%
  tibble::rownames_to_column("Var1") %>%
  tidyr::pivot_longer(cols = -Var1, names_to = "Var2", values_to = "value") %>%
  dplyr::filter(Var1 < Var2)

## Keep non-zero phage correlations only
non_zero_S2 <- assom_S2[(grepl("phage", assom_S2$Var1) | grepl("phage", assom_S2$Var2)) & abs(assom_S2$value) != 0, ]
non_zero_S1 <- assom_S1[(grepl("phage", assom_S1$Var1) | grepl("phage", assom_S1$Var2)) & abs(assom_S1$value) != 0, ]


## How many non-zero phage correlations in total in each network?
net_stats <- rbind(net_stats, c('N non-zero VC correlations', nrow(non_zero_S2), nrow(non_zero_S1)))

## Direction of associations
neg_S2 <- nrow(non_zero_S2[non_zero_S2$value < 0, ])
neg_S1 <- nrow(non_zero_S1[non_zero_S1$value < 0, ])

net_stats <- rbind(net_stats, c('N negative VC correlations', neg_S2, neg_S1))


## Save network statistics
write.csv(net_stats, paste0("results/permutations/permuted_netstats", I, ".csv"))


## Save hubs in each permutation (top 2%)
thresh <- quantile(se_net$centralities$eigenv1, 0.98)
hubs1 <- names(se_net$centralities$eigenv1)[se_net$centralities$eigenv1 >= thresh]

thresh <- quantile(se_net$centralities$eigenv2, 0.98)
hubs2 <- names(se_net$centralities$eigenv2)[se_net$centralities$eigenv2 >= thresh]

write.csv(hubs1, paste0("results/permutation_hubs/S1_hubs_perm", I, ".csv"))
write.csv(hubs2, paste0("results/permutation_hubs/S2_hubs_perm", I, ".csv"))



