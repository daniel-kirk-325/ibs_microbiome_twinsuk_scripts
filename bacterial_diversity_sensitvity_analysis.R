require(lme4)         # 1.1.37
library(dplyr)        # 1.1.4
library(readxl) 
library(compositions) # 2.0.9
library(ggplot2)
library(ggtext)
library(gridExtra)
library(grid)
library(patchwork)
library(vegan)        # 2.7.1


## Load data
df <- read.csv('ibs_microbiome_twins_drug_processed.csv')[, -1]

## Convert or drop sex
df <- df[df$SEX != 'M', ]

## Keep only those with phage data
phage_ids <- read.csv("viroprofiler_iids.csv")[, -1]
df <- df[df$iid %in% phage_ids, ]

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


#### CLR transform compositional data following pseduocount (1/10th lowest value)
low_prev <- names(df[bacteria][colSums(df[bacteria] == 0) > 0.9 * nrow(df)])


## CLR transformation
one_tenth_min_non_zero <- min(df[bacteria][df[bacteria] > 0])/10
clr_df <- df
clr_df[bacteria] <- compositions::clr(clr_df[bacteria] + one_tenth_min_non_zero)


# Create diet dataset -----------------------------------------------------

## Make copy of df since diet is only available in a subset of participants
df_diet <- clr_df

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


# Drug adjustment ---------------------------------------------------------

# Alpha diversity  --------------------------------------------------------

alpha_diversity <- read.csv("datasets/raw_microbiome_data/merged_alpha_diversity_GS.tsv", sep = '\t') %>% 
  dplyr::rename(prefix = Prefix) # Rename column for merging

# Prepare and scale
diversity_df <- merge(alpha_diversity, df[, c('prefix', 'batch', 'fid', 'age_Mb_sample', 'bmi_Mb_sample', 'IBS_overall', 'SEX', drug_vars)], by = 'prefix')
div_vars <- div_vars <- c("observed", "diversity_shannon")

diversity_df[, c('age_scaled', 'BMI_scaled', 'observed', 'diversity_shannon')] <- 
  scale(diversity_df[, c('age_Mb_sample', 'bmi_Mb_sample', 'observed', 'diversity_shannon')])

# Run models and collect p-values
plot_labs <- c('Observed', 'Shannon')
pvals <- c()
names(plot_labs) <- div_vars

## Establish null model
f_null <- paste0("IBS_overall ~ age_scaled + BMI_scaled + batch + (1 | fid) + ", paste0(drug_vars, collapse = '+'))
null <- glmer(f_null, data = diversity_df, family = binomial,
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)), nAGQ = 40)

for (var in div_vars[1:2]) {
  f_drug <- as.formula(paste0("IBS_overall ~ ", var, " + age_scaled + BMI_scaled + batch + (1|fid) + ", paste0(drug_vars, collapse = '+')))
  full <- glmer(f_drug, data = diversity_df, family = binomial,
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)), nAGQ = 40)
  print(summary(full))
  p <- anova(null, full)$`Pr(>Chisq)`[2]
  print(anova(null, full))
  pvals[var] <- ifelse(p < 0.001, "<0.001", paste0("= ", signif(p, 2)))
}


# Beta diversity  ----------------------------------------------------------

input_df <- clr_df[bacteria]

## Add IBS and covs
input_df[, c('IBS', 'iid', 'age_Mb_sample', 'bmi_Mb_sample', 'fid', 'batch', drug_vars)] <- df[, c('IBS_overall', 'iid', 'age_Mb_sample', 'bmi_Mb_sample', 'fid', 'batch', drug_vars)]

## Make IBS factor (easier for plotting)
input_df$IBS <- as.factor(input_df$IBS)
input_df$IBS <- ifelse(input_df$IBS == '0', 'Control', 'IBS')

## Scale age and BMI
input_df[, c('age_scaled', 'BMI_scaled')] <- scale(input_df[, c('age_Mb_sample', 'bmi_Mb_sample')])

aitchison_dist <- dist(input_df[bacteria], method = "euclidean")

# Perform PCA on distance matrix
pca_result <- prcomp(input_df[bacteria])
pca_df <- as.data.frame(pca_result$x)
pca_sum <- summary(pca_result)

pca_df$IBS <- input_df$IBS # Add group information to the data frame

# Calculate the percentage of variance explained
explained_variance <- round(100*pca_sum$importance[2, ],2)
pc1_var <- paste0('PC1 (', round(explained_variance[1],2), '%)') 
pc2_var <- paste0('PC2 (', round(explained_variance[2],2), '%)') 

# Define the formula for PERMANOVA
adonis_f <- as.formula(paste0("aitchison_dist ~ age_scaled + BMI_scaled + batch + ", paste0(drug_vars, collapse = "+"), " + IBS"))
set.seed(0)
permanova_result <- adonis2(adonis_f, data = input_df, permutations = how(blocks = input_df$fid, nperm = 999), by = "margin")

IBS_explained_var <- paste0('IBS R-squared: ', round(100*permanova_result['IBS', 'R2'],2), '%')
IBS_p <- paste0('P-value: ', permanova_result['IBS', 'Pr(>F)'], ' ')
text_box_content <- paste(IBS_explained_var, IBS_p, sep = "\n")


bd_drug_adjusted <- ggplot(pca_df, aes(x = PC1, y = PC2, color = IBS)) +
  geom_point() + ggtitle('Gut Bacterial Beta Diversity in Twins UK\n Aitchison distance (Drug-Adjusted)') +
  stat_ellipse(aes(group = IBS), level = 0.95, geom = "polygon", alpha = 0, size = 1) +
  labs(x = pc1_var, y = pc2_var, color = "IBS Status") +
  scale_color_manual(values = c("#00AFBB", "#F8766D")) +  
  theme_minimal() + 
  theme(plot.margin = unit(c(0.5, 2, 1, 0.05), "lines"), 
        plot.title = element_text(hjust = 0.5, size = 18)) + 
  coord_fixed() + 
  annotate(
    "label", x = Inf, y = Inf, label = text_box_content, hjust = 1.1, vjust = 1.1, 
    fill = "grey90", color = "black", size = 4, label.size = 0.3
  )
bd_drug_adjusted



# Diet adjustment ---------------------------------------------------------

# Alpha diversity  --------------------------------------------------------

alpha_diversity <- read.csv("merged_alpha_diversity_GS.tsv", sep = '\t') %>% 
  dplyr::rename(prefix = Prefix) # Rename column for merging

# Prepare and scale
diversity_df <- merge(alpha_diversity, df_diet_unique[, c('prefix', 'batch', 'fid', 'age_Mb_sample', 'bmi_Mb_sample', 'IBS_overall', 'SEX', drug_vars, food_group_vars)], by = 'prefix')
div_vars <- names(alpha_diversity)[names(alpha_diversity) != 'prefix']

diversity_df[, c('age_scaled', 'BMI_scaled', 'observed', 'diversity_shannon', 'diversity_gini_simpson', 'dominance_gini', food_group_vars)] <- 
  scale(diversity_df[, c('age_Mb_sample', 'bmi_Mb_sample', 'observed', 'diversity_shannon', 'diversity_gini_simpson', 'dominance_gini', food_group_vars)])

# Run models and collect p-values
plot_labs <- c('Observed', 'Shannon', 'Gini Simpson', 'Gini Dominance')
pvals <- c()
names(plot_labs) <- div_vars

## Establish null model
f_null <- as.formula(paste0("IBS_overall ~ age_scaled + BMI_scaled + batch + (1|fid) + ", paste0(food_group_vars, collapse = '+')))

null <- glmer(f_null, data = diversity_df, family = binomial,
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)), nAGQ = 40)

for (var in div_vars) {
  f_diet <- as.formula(paste0("IBS_overall ~ ", var, " + age_scaled + BMI_scaled + batch + (1|fid) + ", paste0(food_group_vars, collapse = '+')))
  full <- glmer(f_diet, data = diversity_df, family = binomial,
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)), nAGQ = 40)
  print(summary(full))
  p <- anova(null, full)$`Pr(>Chisq)`[2]
  print(anova(null, full))
  pvals[var] <- ifelse(p < 0.001, "<0.001", paste0("= ", signif(p, 2)))
}


# Beta diversity  ----------------------------------------------------------

# Ordination using Aitchson distance 

input_df <- df_diet_unique[c(bacteria, food_group_vars)]

## Add IBS and covs
input_df[, c('IBS', 'age_Mb_sample', 'bmi_Mb_sample', 'fid', 'batch')] <- df_diet_unique[, c('IBS_overall', 'age_Mb_sample', 'bmi_Mb_sample', 'fid', 'batch')]

## Make IBS factor (easier for plotting)
input_df$IBS <- as.factor(input_df$IBS)
input_df$IBS <- ifelse(input_df$IBS == '0', 'Control', 'IBS')

## Scale age and BMI
input_df[, c('age_scaled', 'BMI_scaled')] <- scale(input_df[, c('age_Mb_sample', 'bmi_Mb_sample')])

aitchison_dist <- dist(input_df[bacteria], method = "euclidean")

# Perform PCA on distance matrix
pca_result <- prcomp(input_df[bacteria])
pca_df <- as.data.frame(pca_result$x)
pca_sum <- summary(pca_result)

pca_df$IBS <- input_df$IBS # Add group information to the data frame

# Calculate the percentage of variance explained
explained_variance <- round(100*pca_sum$importance[2, ],2)
pc1_var <- paste0('PC1 (', round(explained_variance[1],2), '%)') 
pc2_var <- paste0('PC2 (', round(explained_variance[2],2), '%)') 


# Define the formula for PERMANOVA
adonis_f <- as.formula(paste0("aitchison_dist ~ age_scaled + BMI_scaled + batch + ", paste0(food_group_vars, collapse = "+"), " + IBS"))
set.seed(0)
permanova_result <- adonis2(adonis_f, data = input_df, permutations = how(blocks = input_df$fid, nperm = 999), by = "margin")

IBS_explained_var <- paste0('IBS R-squared: ', round(100*permanova_result['IBS', 'R2'],2), '%')
IBS_p <- paste0('P-value: ', permanova_result['IBS', 'Pr(>F)'], ' ')
text_box_content <- paste(IBS_explained_var, IBS_p, sep = "\n")


bd_diet_adjusted <- ggplot(pca_df, aes(x = PC1, y = PC2, color = IBS)) +
  geom_point() + ggtitle('Aitchison distance (Diet-Adjusted)') +
  stat_ellipse(aes(group = IBS), level = 0.95, geom = "polygon", alpha = 0, size = 1) +
  labs(x = pc1_var, y = pc2_var, color = "IBS Status") +
  scale_color_manual(values = c("#00AFBB", "#F8766D")) +  
  theme_minimal() + 
  theme(plot.margin = unit(c(0.5, 1, 1, 0.05), "lines"), 
        plot.title = element_text(hjust = 0.5, size = 18)) + 
  coord_fixed() + 
  annotate(
    "label", x = Inf, y = Inf, label = text_box_content, hjust = 1.1, vjust = 1.1, 
    fill = "grey90", color = "black", size = 4, label.size = 0.3
  )
bd_diet_adjusted


# Put individual figures on the same plot  --------------------------------

combined_plot_v <- bd_drug_adjusted / bd_diet_adjusted

path <- "beta_diversity_adjusted.png"
ggsave(path, combined_plot_v, width = 8, height = 12, dpi = 300)

















