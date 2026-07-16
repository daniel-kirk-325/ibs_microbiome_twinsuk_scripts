require(lme4)         # 1.1.37
library(dplyr)        # 1.1.4
library(readxl) 
library(DHARMa)       # 0.4.7
library(compositions) # 2.0.9


# Functions ---------------------------------------------------------------

lmm_Mb <- function(input_df, var_list, f_null, cohort, pathname) {
  
  # Check if the folder exists
  if (!file.exists(dirname(pathname))) {stop(paste("Error: Folder does not exist"))}
  
  null <- glmer(f_null, data = input_df, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), nAGQ = 40)
  
  results_df <- data.frame(EntryName = character(0), BETA = numeric(0), SE = numeric(0), OR = numeric(0), 
                           Lower = numeric(0), Upper = numeric(0), warnings = numeric(0), PVAL = numeric(0))
  
  count <- 0 
  for (v in var_list) {
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
  
  results_df_ibs[, !(names(results_df_ibs) == 'PVAL')] <- round(results_df_ibs[, !(names(results_df_ibs) == 'PVAL')], 4)
  results_df_ibs[, names(results_df_ibs) == 'PVAL'] <- round(results_df_ibs[, names(results_df_ibs) == 'PVAL'], 8)
  results_df_ibs <- results_df_ibs %>% arrange(PVAL)
  
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
df[bacteria] <- compositions::clr(df[bacteria] + one_tenth_min_non_zero)

bacteria <- bacteria[!bacteria %in% low_prev]



# Species-level analysis  ----------------------------------------------------------------

input_df <- df

## Scale age and BMI
input_df[, c('age_scaled', 'BMI_scaled')] <- scale(input_df[, c('age_Mb_sample', 'bmi_Mb_sample')])

f_null <- paste("IBS_overall ~ age_scaled + BMI_scaled + batch + (1 | fid)")

#lmm_Mb(input_df, bacteria, f_null, 'twins', "bacteria_analysis.csv")



# Check model assumptions and distributions -------------------------------

## Check GMM assumptions using dharma

top_twins_df <- read.csv("bacteria_analysis.csv")
top_twins <- top_twins_df$X[top_twins_df$q_value < 0.05]


pdf("bacterial_dharma.pdf")
count <- 0
for (m in top_twins) {
  count <- count + 1
  print(100 * (count / length(top_twins)))
  
  f1 <- as.formula(paste0("IBS_overall ~ age_scaled + BMI_scaled + batch +", m, " + (1 | fid)"))
  model <- glmer(f1, data = input_df, family = binomial, 
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)), nAGQ = 40)
  

  # Diagnostics
  simulated_residuals <- simulateResiduals(model, n = 1000)
  
  par(mfrow = c(1, 2))
  plotQQunif(simulated_residuals, testOutliers = FALSE)
  plotResiduals(simulated_residuals)
  mtext(m, side = 3, outer = TRUE, line = -1, cex = 1.5)
  
  par(mfrow = c(1, 1))
  plotResiduals(simulated_residuals, form = input_df[[m]], main = '', xlab = m)
  title(main = m, line = 2.5)
}
dev.off()



# Drug Sensitivity analysis --------------------------------------

f_drugs <- paste("IBS_overall ~ age_scaled + BMI_scaled + batch + (1 | fid) + ", paste0(drug_vars, collapse = '+'))

## Get microbes to replicate
top_twins_df <- read.csv("bacteria_analysis.csv")
top_twins <- top_twins_df$X[top_twins_df$q_value < 0.05]

lmm_Mb(input_df, top_twins, f_drugs, 'twins', "bacteria_analysis_drug_adjusted.csv")



# Diet sensitivity analysis -----------------------------------------------

## Make copy of df since diet is only available in a subset of participants
df_diet <- df

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

hist(as.numeric(df_full$FFQ_Date - df_full$dov), 
     main = 'Time difference between closest FFQ and \nstool sampling date per participant', 
     xlab = 'Days', 
     breaks = 15, 
     col = 'grey90')

## Get number of participants within X amount of years
df_full$df_full <- abs(df_full$FFQ_Date - df_full$dov)
sum(df_full$date_diff < 365.25*1) #78
sum(df_full$date_diff < 365.25*2) #145
sum(df_full$date_diff < 365.25*3) #469
sum(df_full$date_diff < 365.25*4) #578
sum(df_full$date_diff < 365.25*5) #639

mean(df_full$date_diff) # Time difference of 1085.202 days (2.97 years)
sd(df_full$date_diff) # 2.57 years

## Impute missing diet data with means 
df_diet_unique <- df_diet_unique %>%
  mutate(across(all_of(food_group_vars), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

## Scale age and BMI
df_diet_unique[, c('age_scaled', 'BMI_scaled', food_group_vars)] <- scale(df_diet_unique[, c('age_Mb_sample', 'bmi_Mb_sample', food_group_vars)])

f_diet <- paste("IBS_overall ~ age_scaled + BMI_scaled + batch + (1 | fid) + ", paste0(food_group_vars, collapse = '+')) 

lmm_Mb(df_diet_unique, top_twins, f_diet, 'twins', "bacteria_analysis_diet_adjusted.csv")


# BSS adjustment ----------------------------------------------------------

bss <- read.csv("BSS_14days_twins.csv")
bss <- bss %>% rename(BSS_score = Response)

bss_input <- merge(bss[, c('iid', 'BSS_score')], df)

## ! Combine BSS scores 6&7 due to low number of 7 scores
bss_input[bss_input$BSS_score == 7, 'BSS_score'] <- 6 
bss_input$BSS_score <- as.factor(bss_input$BSS_score )

## Scale age and BMI
bss_input[, c('age_scaled', 'BMI_scaled')] <- scale(bss_input[, c('age_Mb_sample', 'bmi_Mb_sample')])

f_null <- paste("IBS_overall ~ age_scaled + BMI_scaled + batch + (1 | fid) + BSS_score")

lmm_Mb(bss_input, top_twins, f_null, 'twins', "bacteria_analysis_BSS_adjusted.csv")


# Replication in PREDICT --------------------------------------------------

top_twins_df <- read.csv("bacteria_analysis.csv")
top_twins <- top_twins_df$X[top_twins_df$q_value < 0.05]


predict <- read.csv('ibs_microbiome_predict_drug_processed.csv')[, -1]
predict$fid <- as.integer(substr(predict$iid, 1, nchar(predict$iid) - 1))
predict <- predict[predict$sex != 'M', ]


#### Perform replication

## Calculate low prevalence bacteria
predict_bacteria_all <- names(predict)[!names(predict) %in% c('iid', 'age', 'sex', 'bmi', 'dov', 'IBS_overall', drug_vars, 'fid')]
low_prev <- names(predict[predict_bacteria_all][colSums(predict[predict_bacteria_all] == 0) > 0.9 * nrow(predict)])
predict_bacteria <- predict_bacteria_all[!predict_bacteria_all %in% low_prev]


## Calculate mean and median relative abundances in final dataset
mean_ras_predict <- colMeans(predict[predict_bacteria_all])
med_ras_predict <- apply(predict[predict_bacteria_all], 2, median)
mean_ra_df_predict <- data.frame(Mean_RA = mean_ras_predict, Median_RA = med_ras_predict)
# write.csv(mean_ra_df_predict, "datasets/processed_microbiome_data/predict/relative_abundances_predict.csv")


## CLR transformation
one_tenth_min_non_zero <- min(predict[predict_bacteria_all][predict[predict_bacteria_all] > 0])/10
predict[predict_bacteria_all] <- clr(predict[predict_bacteria_all] + one_tenth_min_non_zero)



## Drop low prevalence bacteria from top_twins
top_twins <- top_twins[top_twins %in% predict_bacteria] # 41/45 bacteria were prevalent >= 10% in predict


input_df <- predict

## Scale age and BMI
input_df[, c('age_scaled', 'BMI_scaled')] <- scale(input_df[, c('age', 'bmi')])


f_null <- paste("IBS_overall ~ age_scaled + BMI_scaled + (1 | fid)")
#lmm_Mb(input_df, top_twins, f_null, 'predict', 'bacteria_analysis_predict')




