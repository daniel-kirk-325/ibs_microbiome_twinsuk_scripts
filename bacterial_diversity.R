require(lme4)         # 1.1.37
library(dplyr)        # 1.1.4
library(readxl) 
library(compositions) # 2.0.9
library(ggplot2)
library(ggtext)
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


## CLR transformation
one_tenth_min_non_zero <- min(df[bacteria][df[bacteria] > 0])/10
clr_df <- df
clr_df[bacteria] <- compositions::clr(clr_df[bacteria] + one_tenth_min_non_zero)


# Alpha diversity  --------------------------------------------------------

## Load previously calculated alpha diversity
alpha_diversity <- read.csv("merged_alpha_diversity_GS.tsv", sep = '\t') %>% 
  dplyr::rename(prefix = Prefix) # Rename column for merging

# Prepare and scale
diversity_df <- merge(alpha_diversity, df[, c('prefix', 'batch', 'fid', 'age_Mb_sample', 'bmi_Mb_sample', 'IBS_overall', 'SEX', drug_vars)], by = 'prefix')
div_vars <- c('observed', 'diversity_shannon')
diversity_df[, c('age_scaled', 'BMI_scaled', 'observed', 'diversity_shannon')] <- scale(diversity_df[, c('age_Mb_sample', 'bmi_Mb_sample', 'observed', 'diversity_shannon')])

# Run models and collect p-values
plot_labs <- c('Observed', 'Shannon')
OR_p_labs <- c()
names(plot_labs) <- div_vars

## Establish null model
f_null <- "IBS_overall ~ age_scaled + BMI_scaled + batch + (1 | fid)"
null <- glmer(f_null, data = diversity_df, family = binomial,
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)), nAGQ = 40)

for (var in div_vars) {
  f1 <- as.formula(paste0("IBS_overall ~ ", var, " + age_scaled + BMI_scaled + batch + (1|fid)"))
  full <- glmer(f1, data = diversity_df, family = binomial,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)), nAGQ = 40)
  
  full_sum <- summary(full)
  print(full_sum)
  
  ## Get OR[95%CI] of metrics
  beta <- full_sum$coefficients[2, 1]
  SE <- full_sum$coefficients[2, 2]
  OR <- round(exp(beta),2)
  Lower <- round(exp(beta - (1.96*SE)),2)
  Upper <- round(exp(beta + (1.96*SE)),2)
  
  OR_lab <- paste0(OR, ' (', Lower, '-', Upper, ')')
  
  ## Get p-value from anova
  p <- anova(null, full)$`Pr(>Chisq)`[2]
  print(anova(null, full))
  pval <- ifelse(p < 0.001, "<0.001", paste0("= ", signif(p, 2)))
  p_lab <- paste0("\n(p ", pval, ")")
  
  OR_p_labs[var] <- paste0(OR_lab, p_lab)
}

# Reshape data for ggplot
long_df <- reshape2::melt(diversity_df[, c("IBS_overall", div_vars)], id.vars = "IBS_overall")
long_df$plot_lab <- factor(long_df$variable, levels = div_vars, labels = plot_labs)
long_df$IBS_label <- factor(long_df$IBS_overall, levels = c(0, 1), labels = c("Control", "IBS"))

# Add p-value labels for facet titles
plot_lab_with_p <- paste0(plot_labs, '\n', OR_p_labs, '\n')
names(plot_lab_with_p) <- plot_labs
long_df$plot_lab <- factor(long_df$plot_lab, levels = plot_labs, labels = plot_lab_with_p)


# Final plot
alpha_diversity_plot <- ggplot(long_df, aes(x = plot_lab, y = value, fill = IBS_label)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6) +
  scale_fill_manual(values = c("Control" = "#00BFC4", "IBS" = "#F8766D")) +
  labs(x = "Diversity Metric, OR (95% CI)", y = "Standardised Diversity Metric (Z-score)", fill = "") +
  ggtitle("Gut Bacterial Alpha Diversity in TwinsUK") +
  theme_bw() + 
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = .5),
    plot.title = element_text(hjust = 0.5),
    plot.margin = unit(c(2, 2, 2, 2), "lines")
  )
alpha_diversity_plot



# Beta Diversity  ----------------------------------------------------------

# • Aitchson distance --------------------------------------

input_df <- clr_df[bacteria]

## Add IBS and covs
input_df[, c('IBS', 'iid', 'age_Mb_sample', 'bmi_Mb_sample', 'fid', 'batch', drug_vars)] <- df[, c('IBS_overall', 'iid', 'age_Mb_sample', 'bmi_Mb_sample', 'fid', 'batch', drug_vars)]

## Make IBS factor (easier for plotting)
input_df$IBS <- as.factor(input_df$IBS)
input_df$IBS <- ifelse(input_df$IBS == '0', 'Control', 'IBS')

idf_reduced <- input_df # For analysis without related individuals

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
adonis_f <- as.formula(paste0("aitchison_dist ~ age_scaled + BMI_scaled + batch + IBS"))
set.seed(0)
permanova_result <- adonis2(adonis_f, data = input_df, permutations = how(blocks = input_df$fid, nperm = 999), by = "margin")

IBS_explained_var <- paste0('IBS R-squared: ', round(100*permanova_result['IBS', 'R2'],2), '%')
IBS_p <- paste0('P-value: ', permanova_result['IBS', 'Pr(>F)'], ' ')
text_box_content <- paste(IBS_explained_var, IBS_p, sep = "\n")


EU_aitch <- ggplot(pca_df, aes(x = PC1, y = PC2, color = IBS)) +
  geom_point() + ggtitle('Gut Bacterial Beta Diversity in\n TwinsUK (Aitchison distance)') +
  stat_ellipse(aes(group = IBS), level = 0.95, geom = "polygon", alpha = 0, size = 1) +
  labs(x = pc1_var, y = pc2_var, color = "IBS Status") +
  scale_color_manual(values = c("#00AFBB", "#F8766D")) + 
  theme_bw() + 
  theme_minimal() + 
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), 
        plot.title = element_text(hjust = 0.5, size = 18)) + 
  #coord_fixed() + 
  annotate(
    "label", x = Inf, y = Inf, label = text_box_content, hjust = 1.1, vjust = 1.1, 
    fill = "grey90", color = "black", size = 4, label.size = 0.3
  )
p <- EU_aitch
p

disper <- betadisper(aitchison_dist, input_df$IBS)
disp_df <- data.frame(
  IBS_Status = disper$group,
  Dispersion = disper$distances
)
set.seed(0)
perm_test <- permutest(disper, permutations = 999)
p_value <- perm_test$tab[1, "Pr(>F)"]

boxplot(Dispersion ~ IBS_Status, data = disp_df,
        xlab = "IBS Status", ylab = "Dispersion",
        col = c("#00AFBB", "#F8766D"), 
        main = paste0("Aitchison Distance\np = ", signif(p_value, digits = 3)))




# • Non-related subsampling ------------------------------


#### Repeat analysis on non-related individuals to account for inability to 
#### adjust with random effects

pcoa_plots <- c()
par(mfrow = c(2,2))
for (i in c(1:8)) {
  
  set.seed(i)
  IBS_df <- idf_reduced[idf_reduced$IBS == 'IBS', ] %>% group_by(fid) %>% sample_n(1) %>% ungroup()
  
  # 2) Drop fids in the IBS df above
  idf_reduced_min_IBS_dfids <- idf_reduced[!(idf_reduced$fid %in% IBS_df$fid), ]
  
  # 3) Randomly drop one fid of the concordant controls
  control_df <- idf_reduced_min_IBS_dfids %>% group_by(fid) %>% sample_n(1) %>% ungroup()
  
  # 4) Merge df and do permanova
  red_df <- rbind(IBS_df, control_df)
  
  # 5) Process data and perform PERMANOVA
  
  ## Scale age and BMI
  red_df[, c('age_scaled', 'BMI_scaled')] <- scale(red_df[, c('age_Mb_sample', 'bmi_Mb_sample')])
  
  beta_d <- vegdist(red_df[bacteria], method = "euclidean")
  
  pcoa_results <- cmdscale(beta_d, eig = TRUE, k = 2) # k = 2 for 2D ordination
  pcoa_coords <- as.data.frame(pcoa_results$points)
  colnames(pcoa_coords) <- c("pcoa1", "pcoa2")
  pcoa_coords$IBS <- red_df$IBS # Add group information to the data frame
  
  # Calculate the percentage of variance explained
  explained_variance <- (pcoa_results$eig / sum(pcoa_results$eig)) * 100
  pc1_var <- paste0('PCoa1 (', round(explained_variance[1],2), '%)') 
  pc2_var <- paste0('PCoa2 (', round(explained_variance[2],2), '%)') 
  
  # Define the formula for PERMANOVA
  adonis_f <- as.formula(paste0("beta_d ~ age_scaled + BMI_scaled + batch + IBS"))
  set.seed(1)
  permanova_result <- adonis2(adonis_f, data = red_df, by = "margin")
  
  IBS_explained_var <- paste0('IBS R-squared: ', round(100*permanova_result['IBS', 'R2'],2), '%')
  IBS_p <- paste0('P-value: ', permanova_result['IBS', 'Pr(>F)'], ' ')
  text_box_content <- paste(IBS_explained_var, IBS_p, sep = "\n")
  
  p <- ggplot(pcoa_coords, aes(x = pcoa1, y = pcoa2, color = IBS)) +
    geom_point() + ggtitle('', ) +
    stat_ellipse(aes(group = IBS), level = 0.95, geom = "polygon", alpha = 0, size = 1) +
    labs(x = pc1_var, y = pc2_var, color = "IBS Status") +
    scale_color_manual(values = c("#00AFBB", "#F8766D")) +  
    theme_minimal() + 
    theme(plot.margin = unit(c(0.5, 2, 1, 0.05), "lines"), 
          plot.title = element_text(hjust = 0.5)) + 
    coord_fixed() + 
    annotate(
      "label", x = Inf, y = Inf, label = text_box_content, hjust = 1.1, vjust = 1.1, 
      fill = "grey90", color = "black", size = 4, label.size = 0.3
    )
  pcoa_plots[[i]] <- p
  
  disper <- betadisper(beta_d, red_df$IBS)
  disp_df <- data.frame(
    IBS_Status = disper$group,
    Dispersion = disper$distances
  )
  perm_test <- permutest(disper, permutations = 999)
  p_value <- perm_test$tab[1, "Pr(>F)"]
  boxplot(Dispersion ~ IBS_Status, data = disp_df,
          xlab = "IBS Status", ylab = "Dispersion",
          col = c("#00AFBB", "#F8766D"), 
          main = paste0("Test for multivariate homogeneity of \ngroups dispersions: p = ", signif(p_value, digits = 3)))
  
}
par(mfrow = c(1,1))
grid.arrange(grobs = pcoa_plots[1:4], nrow = 2, ncol = 2)
grid.arrange(grobs = pcoa_plots[5:8], nrow = 2, ncol = 2)



# • BC on relative abundance ------------------------------------------------

input_df <- df[bacteria]/100 # Scale relative abundance to between 0 and 1

## Add IBS and covs
input_df[, c('IBS', 'iid', 'age_Mb_sample', 'bmi_Mb_sample', 'fid', 'batch', drug_vars)] <- df[, c('IBS_overall', 'iid', 'age_Mb_sample', 'bmi_Mb_sample', 'fid', 'batch', drug_vars)]

## Make IBS factor (easier for plotting)
input_df$IBS <- as.factor(input_df$IBS)
input_df$IBS <- ifelse(input_df$IBS == '0', 'Control', 'IBS')
idf_reduced <- input_df # For analysis without related individuals


## Scale age and BMI
input_df[, c('age_scaled', 'BMI_scaled')] <- scale(input_df[, c('age_Mb_sample', 'bmi_Mb_sample')])

beta_d <- vegdist(input_df[bacteria], method = "bray")

pcoa_results <- cmdscale(beta_d, eig = TRUE, k = 2) # k = 2 for 2D ordination
pcoa_coords <- as.data.frame(pcoa_results$points)
colnames(pcoa_coords) <- c("pcoa1", "pcoa2")
pcoa_coords$IBS <- input_df$IBS # Add group information to the data frame


# Calculate the percentage of variance explained
explained_variance <- (pcoa_results$eig / sum(pcoa_results$eig)) * 100
pc1_var <- paste0('PCoa1 (', round(explained_variance[1],2), '%)') 
pc2_var <- paste0('PCoa2 (', round(explained_variance[2],2), '%)') 


# Define the formula for PERMANOVA
adonis_f <- as.formula(paste0("beta_d ~ age_scaled + BMI_scaled + batch + IBS"))
set.seed(0)
permanova_result <- adonis2(adonis_f, data = input_df, permutations = permutations <- how(blocks = input_df$fid, nperm = 999), by = 'margin')

IBS_explained_var <- paste0('IBS R-squared: ', round(100*permanova_result['IBS', 'R2'],2), '%')
IBS_p <- paste0('P-value: ', permanova_result['IBS', 'Pr(>F)'], ' ')
text_box_content <- paste(IBS_explained_var, IBS_p, sep = "\n")


BC_ra <- ggplot(pcoa_coords, aes(x = pcoa1, y = pcoa2, color = IBS)) +
  geom_point() + ggtitle('Bray-Curtis (Relative Abundance)') +
  stat_ellipse(aes(group = IBS), level = 0.95, geom = "polygon", alpha = 0, size = 1) +
  labs(x = pc1_var, y = pc2_var, color = "IBS Status") +
  scale_color_manual(values = c("#00AFBB", "#F8766D")) +  
  theme_minimal() + 
  theme(plot.margin = unit(c(0.5, 2, 1, 0.05), "lines"), 
        plot.title = element_text(hjust = 0.5)) + 
  coord_fixed() + 
  annotate(
    "label", x = Inf, y = Inf, label = text_box_content, hjust = 1.1, vjust = 1.1, 
    fill = "grey90", color = "black", size = 4, label.size = 0.3
  )
BC_ra


disper <- betadisper(beta_d, input_df$IBS)
disp_df <- data.frame(
  IBS_Status = disper$group,
  Dispersion = disper$distances
)

set.seed(0)
perm_test <- permutest(disper, permutations = 9999)
p_value <- perm_test$tab[1, "Pr(>F)"]

boxplot(Dispersion ~ IBS_Status, data = disp_df,
        xlab = "IBS Status", ylab = "Dispersion",
        col = c("#00AFBB", "#F8766D"), 
        main = paste0("Bray-Curtis (Relative Abundance)\np = ", signif(p_value, digits = 3)))



# • BC on Arcsin transformed data -------------------------

input_df <- asin(sqrt(df[bacteria]/100)) # Scale relative abundance to between 0 and 1 in order to use arcsine

## Add IBS and covs
input_df[, c('IBS', 'iid', 'age_Mb_sample', 'bmi_Mb_sample', 'fid', 'batch', drug_vars)] <- df[, c('IBS_overall', 'iid', 'age_Mb_sample', 'bmi_Mb_sample', 'fid', 'batch', drug_vars)]

## Make IBS factor (easier for plotting)
input_df$IBS <- as.factor(input_df$IBS)
input_df$IBS <- ifelse(input_df$IBS == '0', 'Control', 'IBS')
idf_reduced <- input_df # For analysis without related individuals


## Scale age and BMI
input_df[, c('age_scaled', 'BMI_scaled')] <- scale(input_df[, c('age_Mb_sample', 'bmi_Mb_sample')])

beta_d <- vegdist(input_df[bacteria], method = "bray")

pcoa_results <- cmdscale(beta_d, eig = TRUE, k = 2) # k = 2 for 2D ordination
pcoa_coords <- as.data.frame(pcoa_results$points)
colnames(pcoa_coords) <- c("pcoa1", "pcoa2")
pcoa_coords$IBS <- input_df$IBS # Add group information to the data frame


# Calculate the percentage of variance explained
explained_variance <- (pcoa_results$eig / sum(pcoa_results$eig)) * 100
pc1_var <- paste0('PCoa1 (', round(explained_variance[1],2), '%)') 
pc2_var <- paste0('PCoa2 (', round(explained_variance[2],2), '%)') 


# Define the formula for PERMANOVA
adonis_f <- as.formula(paste0("beta_d ~ age_scaled + BMI_scaled + batch + IBS"))
set.seed(0)
permanova_result <- adonis2(adonis_f, data = input_df, permutations = permutations <- how(blocks = input_df$fid, nperm = 999), by = 'margin')

IBS_explained_var <- paste0('IBS R-squared: ', round(100*permanova_result['IBS', 'R2'],2), '%')
IBS_p <- paste0('P-value: ', permanova_result['IBS', 'Pr(>F)'], ' ')
text_box_content <- paste(IBS_explained_var, IBS_p, sep = "\n")


BC_as <- ggplot(pcoa_coords, aes(x = pcoa1, y = pcoa2, color = IBS)) +
  geom_point() + ggtitle('Bray-Curtis (Arcsine)') +
  stat_ellipse(aes(group = IBS), level = 0.95, geom = "polygon", alpha = 0, size = 1) +
  labs(x = pc1_var, y = pc2_var, color = "IBS Status") +
  scale_color_manual(values = c("#00AFBB", "#F8766D")) +  
  theme_minimal() + 
  theme(plot.margin = unit(c(0.5, 2, 1, 0.05), "lines"), 
        plot.title = element_text(hjust = 0.5)) + 
  coord_fixed() + 
  annotate(
    "label", x = Inf, y = Inf, label = text_box_content, hjust = 1.1, vjust = 1.1, 
    fill = "grey90", color = "black", size = 4, label.size = 0.3
  )
BC_as


disper <- betadisper(beta_d, input_df$IBS)
disp_df <- data.frame(
  IBS_Status = disper$group,
  Dispersion = disper$distances
)

set.seed(0)
perm_test <- permutest(disper, permutations = 9999)
p_value <- perm_test$tab[1, "Pr(>F)"]

boxplot(Dispersion ~ IBS_Status, data = disp_df,
        xlab = "IBS Status", ylab = "Dispersion",
        col = c("#00AFBB", "#F8766D"), 
        main = paste0("Bray-Curtis (Arcsine)\np = ", signif(p_value, digits = 3)))



# • Diversity comparison plot -----------------------------------------------

## Diversity comparison plot
comb_plot <- BC_ra / BC_as
path <- "beta_diversity_comparison.png"
ggsave(path, comb_plot, width = 8, height = 12, dpi = 600)




# Replication in PREDICT  -------------------------------------------------

### Data pre-processing

predict <- read.csv('ibs_microbiome_predict_drug_processed.csv')[, -1]
predict$fid <- as.integer(substr(predict$iid, 1, nchar(predict$iid) - 1))
predict <- predict[predict$sex != 'M', ]

## If adjusting for drugs:
predict$drug_user <- 0
predict$drug_user <- ifelse(rowSums(predict[drug_vars] == 1) > 0, 1, predict$drug_user)



# • Alpha diversity -------------------------------------------------------

predict_diversity <- read.csv("P1_alpha.tsv", sep = '\t') %>% 
  dplyr::rename(Observed = Richness)

diversity_df <- merge(predict_diversity, predict[, c('iid', 'fid', 'age', 'bmi', 'IBS_overall', 'drug_user')], by = 'iid')

diversity_df[, c('age_scaled', 'BMI_scaled', 'Shannon', 'Observed')] <- scale(diversity_df[, c('age', 'bmi', 'Shannon', 'Observed')])


# Run models and collect p-values
plot_labs <- c('Observed', 'Shannon')
OR_p_labs <- c()
names(plot_labs) <- c('Observed', 'Shannon')

f_null <- "IBS_overall ~ age_scaled + BMI_scaled + (1 | fid)"
null <- glmer(f_null, data = diversity_df, family = binomial,
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)), nAGQ = 40)

for (var in plot_labs) {
  f <- as.formula(paste0("IBS_overall ~ ", var, " + age_scaled + BMI_scaled + (1|fid)"))
  full <- glmer(f, data = diversity_df, family = binomial,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)), nAGQ = 40)
  
  full_sum <- summary(full)
  print(full_sum)
  
  ## Get OR[95%CI] of metrics
  beta <- full_sum$coefficients[2, 1]
  SE <- full_sum$coefficients[2, 2]
  OR <- round(exp(beta),2)
  Lower <- round(exp(beta - (1.96*SE)),2)
  Upper <- round(exp(beta + (1.96*SE)),2)
  
  OR_lab <- paste0(OR, ' (', Lower, '-', Upper, ')')
  
  ## Get p-value from anova
  p <- anova(null, full)$`Pr(>Chisq)`[2]
  print(anova(null, full))
  pval <- ifelse(p < 0.001, "<0.001", paste0("= ", signif(p, 2)))
  p_lab <- paste0("\n(p ", pval, ")")
  
  OR_p_labs[var] <- paste0(OR_lab, p_lab)

}

# Reshape data for ggplot
long_df <- reshape2::melt(diversity_df[, c("IBS_overall", plot_labs)], id.vars = "IBS_overall")
long_df$plot_lab <- factor(long_df$variable, levels = plot_labs, labels = plot_labs)
long_df$IBS_label <- factor(long_df$IBS_overall, levels = c(0, 1), labels = c("Control", "IBS"))

# Add p-value labels for facet titles
plot_lab_with_p <- paste0(plot_labs, '\n', OR_p_labs, '\n')
names(plot_lab_with_p) <- plot_labs
long_df$plot_lab <- factor(long_df$plot_lab, levels = plot_labs, labels = plot_lab_with_p)

# Final plot
predict_alpha <- ggplot(long_df, aes(x = plot_lab, y = value, fill = IBS_label)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6) +
  scale_fill_manual(values = c("Control" = "#00BFC4", "IBS" = "#F8766D")) +
  labs(x = "Diversity Metric, OR (95% CI)", y = "Standardised Diversity Metric (Z-score)", fill = "") +
  ggtitle("Gut Bacterial Alpha Diversity in PREDICT-1") +
  theme_bw() + 
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = .5),
    plot.title = element_text(hjust = 0.5)
  )
#path <- "for_manuscript/diversity/diversity_figures/predict_alpha_diversity.png"
#ggsave(path, width = 10, height = 6, dpi = 600, bg = "white")




# • Beta diversity --------------------------------------------------------

predict_bacteria <- names(predict)[!names(predict) %in% c('iid', 'sex', 'bmi', 'age', 'dov', 'fid', 'IBS_overall', 'drug_user', drug_vars)]


# •• Aitchson distance --------------------------------------

## CLR transformation
one_tenth_min_non_zero <- min(predict[predict_bacteria][predict[predict_bacteria] > 0])/10
clr_df <- predict
clr_df[predict_bacteria] <- compositions::clr(clr_df[predict_bacteria] + one_tenth_min_non_zero)

input_df <- clr_df[predict_bacteria]

## Add IBS and covs
input_df[, c('IBS', 'iid', 'age', 'bmi', 'fid', 'drug_user')] <- predict[, c('IBS_overall', 'iid', 'age', 'bmi', 'fid', 'drug_user')]

## Make IBS factor (easier for plotting)
input_df$IBS <- as.factor(input_df$IBS)
input_df$IBS <- ifelse(input_df$IBS == '0', 'Control', 'IBS')


## Scale age and BMI
input_df[, c('age_scaled', 'BMI_scaled')] <- scale(input_df[, c('age', 'bmi')])


aitchison_dist <- dist(input_df[predict_bacteria], method = "euclidean")

# Perform PCA on distance matrix
pca_result <- prcomp(input_df[predict_bacteria])
pca_df <- as.data.frame(pca_result$x)
pca_sum <- summary(pca_result)

pca_df$IBS <- input_df$IBS # Add group information to the data frame

# Calculate the percentage of variance explained
explained_variance <- round(100*pca_sum$importance[2, ],2)
pc1_var <- paste0('PC1 (', round(explained_variance[1],2), '%)') 
pc2_var <- paste0('PC2 (', round(explained_variance[2],2), '%)') 

# Define the formula for PERMANOVA
adonis_f <- as.formula(paste0("aitchison_dist ~ age_scaled + BMI_scaled + IBS"))
set.seed(0)
permanova_result <- adonis2(adonis_f, data = input_df, permutations = permutations <- how(blocks = input_df$fid, nperm = 999), by = 'margin')


IBS_explained_var <- paste0('IBS R-squared: ', round(100*permanova_result['IBS', 'R2'],2), '%')
IBS_p <- paste0('P-value: ', permanova_result['IBS', 'Pr(>F)'], ' ')
text_box_content <- paste(IBS_explained_var, IBS_p, sep = "\n")


predict_EU_aitch <- ggplot(pca_df, aes(x = PC1, y = PC2, color = IBS)) +
  geom_point() + ggtitle('Gut Bacterial Beta Diversity in\n PREDICT-1 (Aitchison distance)') +
  stat_ellipse(aes(group = IBS), level = 0.95, geom = "polygon", alpha = 0, size = 1) +
  labs(x = pc1_var, y = pc2_var, color = "IBS Status") +
  scale_color_manual(values = c("#00AFBB", "#F8766D")) +  
  theme_minimal() + 
  theme(plot.margin = unit(c(0.5, 2, 1, 0.5), "lines"), 
        plot.title = element_text(hjust = 0.5, size = 18)) + 
  #coord_fixed() + 
  annotate(
    "label", x = Inf, y = Inf, label = text_box_content, hjust = 1.1, vjust = 1.1, 
    fill = "grey90", color = "black", size = 4, label.size = 0.3
  )
predict_EU_aitch

disper <- betadisper(aitchison_dist, input_df$IBS)
disp_df <- data.frame(
  IBS_Status = disper$group,
  Dispersion = disper$distances
)
set.seed(0)
perm_test <- permutest(disper, permutations = 999)
p_value <- perm_test$tab[1, "Pr(>F)"]

boxplot(Dispersion ~ IBS_Status, data = disp_df,
        xlab = "IBS Status", ylab = "Dispersion",
        col = c("#00AFBB", "#F8766D"), 
        main = paste0("Aitchison Distance\np = ", signif(p_value, digits = 3)))
dev.off()



















