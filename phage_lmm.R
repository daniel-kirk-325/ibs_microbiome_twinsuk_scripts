require(lme4) 		# 1.1.37
library(dplyr) 		# 1.1.4
library(readxl) 
library(compositions) 	# 2.0.9
library(stringdist)
library(ggplot2)


# Functions ---------------------------------------------------------------

lmm_phage <- function(input_df, var_list, f_null) {
  
  null <- glmer(f_null, data = input_df, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), nAGQ = 40)
  
  results_df <- data.frame(EntryName = character(0), BETA = numeric(0), SE = numeric(0), OR = numeric(0), 
                           Lower = numeric(0), Upper = numeric(0), warnings = numeric(0), PVAL = numeric(0))
  
  count <- 0 
  for (v in var_list) {
    
    # progress_indices <- c(floor(length(var_list) * 0.25), floor(length(var_list) * 0.5), floor(length(var_list) * 0.75))
    # progress_indicators <- var_list[progress_indices]
    # 
    # if (v %in% progress_indicators) {
    #   print(match(v, var_list) * 100 / length(var_list))
    # }
    
    print(match(v, var_list)*100 / length(var_list))
    
    
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
  
  results_df_ibs$warnings <- ifelse(results_df_ibs$warnings > 0, 1, 0)
  
  results_df_ibs <- cbind(rownames(results_df_ibs), results_df_ibs) ; rownames(results_df_ibs) <- NULL
  
  colnames(results_df_ibs)[1] <- "Virus"
  
  return(results_df_ibs) 
  
}

# Data pre-processing -----------------------------------------------------

## Load data
df <- read.csv('ibs_microbiome_twins_drug_processed.csv')[, -1]

## Convert or drop sex
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

drug_vars <- c("Proton.Pump.Inhibitors", "Anti_Depressants.anxiolytics", "Antispasmodics.antimotility", "Analgesics", 'Laxatives')


## CLR transformation of bacteria - required for analysis after phages 
one_tenth_min_non_zero <- min(df[bacteria][df[bacteria] > 0])/10
df[bacteria] <- clr(df[bacteria] + one_tenth_min_non_zero)



# Load phage data ------------------------------------------------------------------

## Add phage data
phages <- read.csv("merged_tpm_COMBINED_UNIQUE.csv")
print(dim(phages))

phages_vigs <- names(phages)[grep("^phage_[0-9]+", names(phages))]
print("LENGTH OF PHAGES VIGS:")
print(length(phages_vigs))


## Add GSID key to dataframe
GSID_key <- read.csv("list of ids with metagenomes and date of visit.txt", sep = '\t')
GSID_key <- GSID_key[GSID_key$prefix %in% df$prefix, ]
df <- merge(GSID_key[, c('prefix', 'gsid')], df, by = 'prefix')

# Replace dashed with underscores in df$gsid to match format in phage data
df$gsid <- gsub('-', '_', df$gsid)

## Merge phages based on GSID
df <- merge(df, phages[!names(phages) %in% c("X", "prefix", "iid", "sex", "dov", "batch", 'Batch', "col_name")], by = 'gsid')

## Check number of all zero phage relative abundance rows after merging
sum(rowSums(df[df$IBS_overall == 0, phages_vigs]) == 0) # 2
sum(rowSums(df[df$IBS_overall == 1, phages_vigs]) == 0) # 0



# Replication cycle comparison --------------------------------------------

####
# Sum the relative abundance of all virulent phages and see if this differs
# in IBS compared to healthy controls
####

## Load taxonomy data
taxonomy <- read.csv("phage_contig_map_with_taxonomy_and_life_cycle.csv")

## Get phages based on replication cycle
virulent_phages <- taxonomy$sudo_phage_ID[taxonomy$life_cycle == 'Virulent']
temp_phages <- taxonomy$sudo_phage_ID[taxonomy$life_cycle == 'Temperate']
equal_phages <- taxonomy$sudo_phage_ID[taxonomy$life_cycle == 'EQUAL']


## Isolate subsets and sum relative abundances
sum_lytic <- rowSums(df[, virulent_phages])
sum_temp <- rowSums(df[, temp_phages])
sum_equal <- rowSums(df[, equal_phages])


## Make df containing only replication cycle sums and metadata
rep_df <- data.frame('sum_lytic' = sum_lytic, 'sum_temp' = sum_temp, 'sum_equal' = sum_equal)
rep_df[, c('age_Mb_sample', 'bmi_Mb_sample', 'batch', 'IBS_overall', 'fid')] <- df[, c('age_Mb_sample', 'bmi_Mb_sample', 'batch', 'IBS_overall', 'fid')]
rep_df[, c('age_scaled', 'BMI_scaled')] <- scale(rep_df[, c('age_Mb_sample', 'bmi_Mb_sample')])


## CLR-transform the summed columns
rep_cols <- c('sum_lytic', 'sum_temp', 'sum_equal')
rep_df[, rep_cols] <- compositions::clr(rep_df[, rep_cols])


## Compare differences
f_null <- paste("IBS_overall ~ age_scaled + BMI_scaled + batch + (1 | fid)")
null <- glmer(f_null, data = rep_df, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), nAGQ = 40)

f_full <- paste("IBS_overall ~ age_scaled + BMI_scaled + batch + (1 | fid) + sum_temp")
full <- glmer(f_full, data = rep_df, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), nAGQ = 40)

anova(null, full) # p = 0.8668

#### Plot data to verify 
plot_df <- data.frame('sum_lytic' = sum_lytic, 'sum_temp' = sum_temp, 'sum_equal' = sum_equal, 'IBS' = df$IBS_overall)

plot_long <- reshape2::melt(
  plot_df,
  id.vars = "IBS",
  variable.name = "replication_cycle",
  value.name = "value"
)

# Clean labels
plot_long$replication_cycle <- factor(plot_long$replication_cycle,
                                      levels = c("sum_lytic", "sum_temp", "sum_equal"),
                                      labels = c("Virulent", "Temperate", "Equal"))

plot_long$IBS <- ifelse(plot_long$IBS == 0, 'Healthy Control', 'IBS')

ggplot(plot_long, aes(x = replication_cycle, y = value, fill = IBS)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = 21,        
    outlier.fill = NA,         
    outlier.color = "black",   
    outlier.alpha = 0.3,       
    outlier.size = 2,         
    outlier.stroke = 0.8
  ) +
  scale_fill_manual(
    values = c("Healthy Control" = "#00BFC4", 
               "IBS"  = "#F8766D") 
  ) +
  labs(
    title = 'Relative Abundance of Viral Contigs by \nReplication Cycle in IBS vs Healthy Controls',
    x = "",
    y = "Relative Abundance",
    fill = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


# LMM ---------------------------------------------------------------------

# Use LMM to identify IBS-associated phages 

# Identify column names where more than 90% of values are 0
phage_prev <- colSums(df[phages_vigs] == 0)
to_drop <- names(phage_prev[phage_prev > 0.9 * nrow(df)])


# Identify phages with RA below threshold
ra_threshold <- 1e-3 # (0.001%)
phage_ra <- colMeans(df[phages_vigs[!phages_vigs %in% to_drop]])
low_ra <- names(phage_ra[phage_ra < ra_threshold])

to_drop <- unique(c(to_drop, low_ra))


## CLR tranformation of phages
min_nonzero <- min(unlist(df[phages_vigs])[unlist(df[phages_vigs]) > 0])
df[phages_vigs] <- compositions::clr(df[phages_vigs] + (min_nonzero / 10))


## Drop low prev phages from df
df <- df[, !(colnames(df) %in% to_drop)]


## Update phages_vigs to exclude dropped columns
phages_vigs <- setdiff(phages_vigs, to_drop)


input_df <- df

## Scale variables
input_df[, c('age_scaled', 'BMI_scaled')] <- scale(input_df[, c('age_Mb_sample', 'bmi_Mb_sample')])

## Establish formula
f_null <- paste("IBS_overall ~ age_scaled + BMI_scaled + batch + (1 | fid)")

## Run LMM function
final_results <- lmm_phage(input_df, phages_vigs, f_null)

## Sort by p-values
final_results <- final_results %>% arrange(PVAL)

## Add q-values
final_results$q_value <- p.adjust(final_results$PVAL, method = 'fdr')

write.csv(final_results, "phage_lmm.csv")


# Adjust for predicted host ------------------------------------------------

#### Load IBS-associated phages
top_phages_df <- read.csv("phage_lmm.csv")
top_phages <- top_phages_df$Virus[top_phages_df$q_value < 0.05]


#### Get predicted host of IBS-associated phages
all_names <- read.csv("all_Host_prediction_to_genome.csv")

# Add pseudo ID to naming dataset
viral_naming_map <- read.csv("phage_contig_map_with_taxonomy_and_life_cycle.csv")
viral_naming_map <- viral_naming_map %>% rename(Virus = phage_ID)
all_names <- merge(viral_naming_map[, c('Virus', 'sudo_phage_ID')], all_names, by = 'Virus')

## Merge predicted host if genus is the same (e.g., Blautia_A --> Blautia)

# Retain genus only 
all_names$host_genus <- sub(".*g__(.*);s__.*", "\\1", all_names$Host.taxonomy)

# Replace underscore plus capital letter with blank space
all_names$host_genus_merged <- all_names$host_genus
all_names[grepl('_', all_names$host_genus), 'host_genus_merged'] <- gsub("_(?!.*_).*$", "", all_names[grepl('_[A-Z]', all_names$host_genus), 'host_genus'], perl = TRUE)

# Keep only IBS-associated phages
to_keep <- viral_naming_map$Virus[viral_naming_map$sudo_phage_ID %in% top_phages]
all_names <- all_names[all_names$Virus %in% to_keep, ]



#### Get genus-level bacterial data
bacterial_mg = read.csv('merged_metaphlan_bugs_abundance_table_20211201_version_4.tsv', sep = '\t')
genus_level <- bacterial_mg[!grepl("s__", bacterial_mg$clade_name) & grepl("g__", bacterial_mg$clade_name), ]

# Drop samples not in our dataset
genus_level <- genus_level[, c('clade_name', paste0('X', df$prefix))]

## Keep only genus name
genus_level$clade_name <- sub(".*g__", "", genus_level$clade_name)


# Transpose genus-level dataset
genus_level <- as.data.frame(t(genus_level))
genus_level <- tibble::rownames_to_column(genus_level, var = "prefix")
genus_level$prefix <- gsub('X', '', genus_level$prefix)
names(genus_level) <- genus_level[1, ]
genus_level <- genus_level[-1, ]
names(genus_level)[1] <- 'prefix'


#### CLR-transform genus-level bacterial relative abundance 
genus_level[-1] <- lapply(genus_level[-1], as.numeric)

# Identify low prev genera 
low_prev <- names(genus_level[-1][colSums(genus_level[-1] == 0) > 0.9 * nrow(genus_level)])

min_nonzero <- 1e-6 # Use same as species analysis 
genus_level[-1] <- compositions::clr(genus_level[-1] + one_tenth_min_non_zero)

# Add CLR-transformed genus-level bacterial relative abundance to df
df <-merge(df, genus_level, by = 'prefix')


#### Keep only genera in our dataset
unique(all_names$host_genus_merged) # 25 genera
genus_present <- unique(all_names$host_genus_merged)[unique(all_names$host_genus_merged) %in% names(genus_level)] # 12
genus_absent <- unique(all_names$host_genus_merged)[!unique(all_names$host_genus_merged) %in% names(genus_level)] # 13

print(genus_absent)
# [1] "Hominicoprocola"      "CAG-103"              "Vescimonas"           "UBA11524"             "Acetatifactor"        "Ventricola"          
# [7] "Alitiscatomonas"      "Faecivicinus"         "SFHK01"               "CAG-170"              "Hominenteromicrobium" "Thomasclavelia"   
# [13] "Agathobacter" 

# Retain IBS-associated phages which infect hosts present in bacterial data
host_in_twins <- all_names[all_names$host_genus %in% genus_present, ]
known_host <- top_phages[top_phages %in% host_in_twins$sudo_phage_ID] # length = 18



#### Repeat LLM while adjusting for species relative abundance

input_df <- df[c('fid', 'age_Mb_sample', 'bmi_Mb_sample', 'batch', known_host, genus_present, 'IBS_overall')]
input_df[, c('age_scaled', 'BMI_scaled')] <- scale(input_df[, c('age_Mb_sample', 'bmi_Mb_sample')])


results_df <- data.frame(EntryName = character(0), BETA = numeric(0), SE = numeric(0), OR = numeric(0), 
                         Lower = numeric(0), Upper = numeric(0), warnings = numeric(0), PVAL = numeric(0))


for (phage in known_host) {
  
  print(match(phage, known_host) * 100 / length(unique(known_host)))
  
  # Identify predicted hosts of phage
  predicted_hosts <- unique(host_in_twins$host_genus_merged[host_in_twins$sudo_phage_ID == phage])
  
  # Add predicted hosts to formula
  f_null <- paste("IBS_overall ~ age_scaled + BMI_scaled + batch + (1 | fid) + ", paste0(predicted_hosts, collapse = '+'))
  
  print(phage)
  print(f_null)
  cat("\n")
  
  # Run null model
  null <- glmer(f_null, data = input_df, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), nAGQ = 40)
  
  # Run full model
  fsplit <- strsplit(f_null, '~')
  f1 <- as.formula(paste(fsplit[[1]][1], '~', phage, '+', fsplit[[1]][2]))
  full <- glmer(f1, data = input_df, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), nAGQ = 40)
  full_sum <- summary(full)
  
  # Get model fit statistics
  warnings <- length(full_sum$optinfo$warnings) + length(full_sum$optinfo$conv$lme4$messages) + length(full_sum$fitMsgs)
  
  pos <- 2 # Variable is second in the output (including intercept)
  
  Vcov <- vcov(full, useScale = FALSE)
  BETA <- full@beta[pos]
  se1 <- sqrt(diag(Vcov))
  SE <- se1[pos]
  OR <- exp(BETA)
  LL <- exp(BETA - 1.96 * SE)
  UL <- exp(BETA + 1.96 * SE)
  PVAL <- anova(null, full)$`Pr(>Chisq)`[pos]
  
  # Get string of predicted hosts
  Bacteria <- paste0(predicted_hosts, collapse = ',')
  
  # Add the results to the results dataframe
  row <- data.frame(Bacteria = Bacteria, BETA = BETA, SE = SE, OR = OR, Lower = LL, Upper = UL, warnings = warnings, PVAL = PVAL)
  rownames(row) <- phage
  results_df <- rbind(results_df, row)
}

# Formmat results df
results_df[c('BETA', 'SE', 'OR', 'Lower', 'Upper')] <- round(results_df[c('BETA', 'SE', 'OR', 'Lower', 'Upper')], 4)
results_df['PVAL'] <- round(results_df$PVAL, 8)
results_df <- results_df %>% arrange(PVAL)

write.csv(results_df, "phage_lmm_host_adjusted.csv")



# Correlation between IBS-associated phages and bacteria ------------

## IBS-associated bacteria
top_bacteria_df <- read.csv("bacteria_analysis.csv")
top_bacteria <- top_bacteria_df$X[top_bacteria_df$q_value < 0.05]

## IBS-associated phages
top_phages_df <- read.csv("phage_lmm.csv")
top_phages <- top_phages_df$Virus[top_phages_df$q_value < 0.05]


## Calculate correlations between IBS-associated bacteria and IBS-associated phages
cor_res <- data.frame(bacterium = character(), phage = character(), rho = numeric(), p_value = numeric(), stringsAsFactors = FALSE)
for (bact in top_bacteria) {
  for (phage in top_phages) {
    
    test_res <- cor.test(df[[bact]], df[[phage]], method = "spearman")
    
    cor_res <- rbind(cor_res,
                     data.frame(
                       bacterium = bact,
                       phage     = phage,
                       rho       = as.numeric(test_res$estimate),   
                       p_value   = test_res$p.value
                     )
    )
  }
}

# Sort by p-values & add FDR
cor_res <- cor_res %>% arrange(p_value)
cor_res$q_value <- p.adjust(cor_res$p_value, method = 'fdr')

#write.csv(cor_res, "ibs_assoc_phage_bacteria_corr.csv")

# Number of significant correlations:
sum(cor_res$q_value < 0.05) # 1708


## Keep bacteria which replicated in PREDICT
predict <- read.csv("bacteria_analysis_predict.csv")
predict <- predict$X[predict$PVAL < 0.05]
cor_res <- cor_res[cor_res$bacterium %in% predict, ]

## Keep only bacteria and VCs with at least one significant assoc
cor_res$q_value <- p.adjust(cor_res$p_value, method = 'fdr')

# Format phage pseudo IDs
cor_res$phage_format <- gsub('_', ' ', cor_res$phage)

#### Group phages based on hosts they infect 

## Add genus data to cor_res
host_info <- all_names[c('Virus', 'host_genus_merged')] %>% distinct()

## Add pseudo phage ID to host_info
host_info$sudo_phage_ID <- viral_naming_map$sudo_phage_ID[match(host_info$Virus, viral_naming_map$Virus)]

## Add genus data to cor_res
cor_res_genus <- merge(cor_res, host_info[, c('sudo_phage_ID', "host_genus_merged")], by.x = 'phage', by.y = 'sudo_phage_ID')

## Genera which only appear once --> Other
appears_once <- names(table(cor_res_genus$host_genus_merged)[table(cor_res_genus$host_genus_merged) == 24])
cor_res_genus[cor_res_genus$host_genus_merged %in% appears_once, 'host_genus_merged'] <- 'Other'

## Add family and order to unknown hosts
cor_res_genus[cor_res_genus$host_genus_merged == 'UBA11524', 'host_genus_merged'] <- '(Clostridia) UBA11524'

## Sort according to most prevalent hosts
ordering <- sort(table(cor_res$host_genus_merged))


cor_res_genus <- cor_res_genus %>%
  mutate(
    # Replace underscores with spaces using base R gsub
    bacterium = gsub("_", " ", bacterium) %>%
      trimws(),                    # equivalent of str_trim()
    
    # Order host genera by frequency, force "Other" to the very end
    host_genus_merged = forcats::fct_infreq(host_genus_merged) %>%
      forcats::fct_relevel("Other", after = Inf)
  ) %>%
  droplevels()


ggplot(cor_res_genus, aes(x = bacterium, y = phage_format, fill = rho)) +
  geom_tile(color = "black", size = 0.3) +
  scale_fill_gradient2(
    low = "#2C7BB6", mid = "white", high = "#D7191C",
    midpoint = 0, limits = c(-1, 1),
    name = "Spearman's\nRho",
    na.value = "grey90"
  ) +
  labs(
    title = "Correlation between IBS-associated \nViral Contigs & IBS-associated Bacteria",
    x = "",
    y = "Phage Pseudo ID Grouped by Predicted Host Genus"
  ) +
  
  facet_grid(host_genus_merged ~ .,
             scales = "free_y", space = "free_y", switch = "y") +
  
  theme_minimal(base_size = 16) +
  
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    
    # X-axis labels: clean names + italics
    axis.text.x = element_text(
      angle = 60, hjust = 1, vjust = 1,
      size = 15, face = "italic"
    ),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    
    axis.text.y = element_text(size = 8, face = "italic"),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(color = "black", size = 0.5),    
    # Legend
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.key.size = unit(1.3, "cm"),
    
    # Facet strips (host genera on the left)
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 0, size = 12, face = "bold", hjust = 1,
      margin = margin(r = 10)
    ),
    strip.background = element_rect(fill = NA, colour = NA),
    
    panel.grid = element_blank(),
    panel.spacing.y = unit(1.2, "lines"),
    plot.margin = margin(l = 45, r = 15, t = 15, b = 15)
  )
path_png <- "ibs_assoc_phage_bacteria_corr.png"
#ggsave(path_png, width = 12, height = 19, dpi = 600, bg = "white")



# Bacteria-phage correlation in IBS-phages with unkn. hosts ---------------

####
# See which bacterial species the IBS-associated phages infecting UBA11524,
# CAG-170, CAG-103 and SFHK01 correlate most strongly with
####

## Load IBS-associated phages
top_phages_df <- read.csv("phage_lmm.csv")
top_phages <- top_phages_df$Virus[top_phages_df$q_value < 0.05]


## Load all tested bacteria (10% filtered)
filtered_bacteria_df <- read.csv("bacteria_analysis.csv")
filtered_bacteria <- filtered_bacteria_df$X


## Establish unknown hosts
unknown_hosts <- c('UBA11524', 'CAG-170', 'CAG-103', 'SFHK01')


## Get IBS-associated phages infecting these hosts
unknown_hosts_phages <- unique(all_names$sudo_phage_ID[all_names$host_genus_merged %in% unknown_hosts])
bacteria_phage_cor <- cor(df[unknown_hosts_phages], df[filtered_bacteria], method = 'spearman') 

# Convert matrix to long format
long_cor <- reshape2::melt(bacteria_phage_cor) ; names(long_cor) <- c('sudo_phage_ID', 'bacteria', 'rho')

## Add host information to long_cor
long_cor <- merge(long_cor, unique(all_names[, c('sudo_phage_ID', 'host_genus_merged')]), by = 'sudo_phage_ID') ; long_cor <- long_cor[, c(1:2, 4, 3)]

# Drop hosts that are not relevant to the current analysis (some of these phages infect multiple hosts)
long_cor <- long_cor[long_cor$host_genus_merged %in% unknown_hosts, ]


#### Add bacterial family information to long_cor 

## Load data
bacteria_taxonomy <- read.csv("microbiome_renamed.csv")[-1] # Load dataset containing formatted names for plotting

## Get species and corresponding family
bacteria_taxonomy <- data.frame(
  clade_name = bacteria_taxonomy$clade_name,
  bacteria = gsub(".*s__", '\\1', bacteria_taxonomy$clade_name),                   # Match name in long_cor for merging below
  bacteria_family = gsub(".*f__(.*)\\|g__.*", '\\1', bacteria_taxonomy$clade_name) # Format name as "bacteria_family" for clarity to distinguish from host family after merging below 
)

## Add this information to long_cor
long_cor <- merge(long_cor, bacteria_taxonomy[, c('bacteria', 'bacteria_family')], by = 'bacteria') ; long_cor <- long_cor[, c(1, 5, 2:4)]

## Get top correlations for phages infecting each unknown predicted host
long_cor %>% group_by(host_genus_merged, sudo_phage_ID) %>% arrange(desc(abs(rho)), by_group = T) %>% slice(1:30) %>% group_split()














