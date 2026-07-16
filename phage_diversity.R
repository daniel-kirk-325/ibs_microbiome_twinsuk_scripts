library(dplyr)
require(lme4) # 1.1.37
library(vegan) # 2.7.1
library(compositions) # 2.0.9
library(ggplot2)
library(ggtext)
library(gridExtra)
library(grid)
library(patchwork)

# Data pre-processing -----------------------------------------------------

## Load data
df <- read.csv('ibs_microbiome_twins_drug_processed.csv')[, -1]

## Convert or drop sex
df <- df[df$SEX != 'M', ]

## Add family ID
df$fid <- as.integer(substr(df$iid, 1, nchar(df$iid) - 1))

## Set batch as factor 
df$batch <- as.factor(df$batch)


drug_vars <- c("Proton.Pump.Inhibitors", "Anti_Depressants.anxiolytics", "Antispasmodics.antimotility", "Analgesics", 'Laxatives')


# Phages ------------------------------------------------------------------

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



# Alpha diversity  --------------------------------------------------------

div_vars <- c('phage_alpha_div_taxonomic_richness', 'phage_alpha_div_shannon')

## Add metadata
diversity_df <- df[, c('fid', 'age_Mb_sample', 'bmi_Mb_sample', 'batch', 'IBS_overall', div_vars)]

## Scale variables
diversity_df[, c('age_scaled', 'BMI_scaled', div_vars)] <- scale(diversity_df[, c('age_Mb_sample', 'bmi_Mb_sample', div_vars)])

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
  labs(x = "OR (95% CI)", y = "Standardised Diversity Metric (Z-score)", fill = "") +
  ggtitle("Alpha diversity of Viral Contigs\n in IBS and Healthy Controls") +
  theme_bw() + 
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.title = element_text(hjust = 0.5, size = 15)
  )
path <- "virome_alpha_diversity.png"
#ggsave(path, alpha_diversity_plot, width = 8, height = 7.7, dpi = 600, bg = "white")


# Beta diversity by Aitchson distance -------------------------------------

## CLR tranformation of phages
min_nonzero <- min(unlist(df[phages_vigs])[unlist(df[phages_vigs]) > 0])
df[phages_vigs] <- clr(df[phages_vigs] + (min_nonzero / 10))

## Make IBS factor (easier for plotting)
df$IBS <- as.factor(df$IBS)
df$IBS <- ifelse(df$IBS == '0', 'Control', 'IBS')

## Scale age and BMI
df[, c('age_scaled', 'BMI_scaled')] <- scale(df[, c('age_Mb_sample', 'bmi_Mb_sample')])

aitchison_dist <- dist(df[phages_vigs], method = "euclidean")

# Perform PCA on distance matrix
pca_result <- prcomp(df[phages_vigs])
pca_df <- as.data.frame(pca_result$x)
pca_sum <- summary(pca_result)

pca_df$IBS <- df$IBS # Add group information to the data frame

# Calculate the percentage of variance explained
explained_variance <- round(100*pca_sum$importance[2, ],2)
pc1_var <- paste0('PC1 (', round(explained_variance[1],2), '%)') 
pc2_var <- paste0('PC2 (', round(explained_variance[2],2), '%)') 


# Define the formula for PERMANOVA
adonis_f <- as.formula(paste0("aitchison_dist ~ age_scaled + BMI_scaled + batch + IBS"))
set.seed(0)
permanova_result <- adonis2(adonis_f, data = df, permutations = how(blocks = df$fid, nperm = 999), by = "margin")

IBS_explained_var <- paste0('IBS R-squared: ', round(100*permanova_result['IBS', 'R2'],2), '%')
IBS_p <- paste0('P-value: ', permanova_result['IBS', 'Pr(>F)'], ' ')
text_box_content <- paste(IBS_explained_var, IBS_p, sep = "\n")


beta_diversity <- ggplot(pca_df, aes(x = PC1, y = PC2, color = IBS)) +
  geom_point() + 
  ggtitle('Viral Contigs Beta Diversity\n (Aitchison distance)') +
  stat_ellipse(aes(group = IBS), level = 0.95, geom = "polygon", alpha = 0, size = 1) +
  labs(x = pc1_var, y = pc2_var, color = "") +
  scale_color_manual(values = c("#00AFBB", "#F8766D")) +  
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off") +
  theme_minimal() + 
  theme(
    #plot.margin = margin(0.5,0.5,0.5,0.5),
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.title = element_text(hjust = 0.5, size = 15)
  ) + 
  annotate(
    "label",
    x = max(pca_df$PC1) * 1.02,
    y = max(pca_df$PC2) * 1.02,
    label = text_box_content,
    hjust = 1, vjust = 1,
    fill = "grey90", color = "black",
    size = 4, label.size = 0.3
  )
path <- "virome_beta_diversity.png"
ggsave(path, beta_diversity, width = 8, height = 7.7, dpi = 600, bg = "white")


disper <- betadisper(aitchison_dist, df$IBS)
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
        main = paste0("Aitchison Distance p = ", signif(p_value, digits = 3)))




