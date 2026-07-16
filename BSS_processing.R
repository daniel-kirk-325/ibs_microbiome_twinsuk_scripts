library(readxl) # 1.4.5
library(dplyr)  # 1.1.4


# CODE EXPLANATION  -------------------------------------------------------

## Gets the closest BSS for participants in with metagenomic and retains only
## those within 14 days of stool sample deposition date

## Only data from Twins is eventually retained, as BSS for PREDICT-1 is 
## available in a separate file and is only used in descriptive statistics
## (but not for analyses)

# Data processing ---------------------------------------------------------

bss1 <- read_xlsx("BSS.xlsx", sheet = 'Phenobase')

bss1 <- bss1[!grepl('last 3 months', bss1$`Phenotype (Column Heading)`), ]

bss2 <- read_xlsx("BSS.xlsx", sheet = 'MSGQv2_Data File_22.08.2023')
bss3 <- read_xlsx("BSS.xlsx", sheet = 'MSGQv3.1_Data File_19.09.2023')

bss1_ <- as.data.frame(bss1[, c('ParticipantID', 'Response Date', 'Response')])
bss2_ <- as.data.frame(bss2[, c('ParticipantID', 'VISIT_DATE', 'Q12')])[-1, ]
bss3_ <- as.data.frame(bss3[, c('ParticipantID', 'VISIT_DATE', 'Q30')])[-1, ]

names(bss2_) <- c('ParticipantID', 'Response Date', 'Response')
names(bss3_) <- c('ParticipantID', 'Response Date', 'Response')

bss_df <- rbind(bss1_, bss2_)
bss_df <- rbind(bss_df, bss3_)

## Drop duplicated entries 

bss_df <- distinct(bss_df) %>% arrange('ParticipantID', 'Response Date')
# Drop unknown responses
bss_df <- bss_df[!bss_df$Response %in% c(999905, 999906, 999911), ]

## See if there are any conflicting entries for BSS submitted on same date
bss_df %>% group_by(ParticipantID, `Response Date`) %>% filter(n_distinct(`Response Date`) > 1) # 0 rows


### Find closest entries 
twins_df <- read.csv("ibs_microbiome_twins_drug_processed.csv")[, -1]
twins_df <- twins_df[twins_df$SEX != 'M', ]

# Keep only those with phage data
phage_ids <- read.csv("viroprofiler_iids.csv")[, -1]
twins_df <- twins_df[twins_df$iid %in% phage_ids, ]

twins <- twins_df[, c('iid', 'dov')]
twins$Cohort <- 'Twins'

predict_df <- read.csv("ibs_microbiome_predict_drug_processed.csv")[, -1]
predict_df <- predict_df[predict_df$sex != 'M', ]
predict <- predict_df[, c('iid', 'dov')]
predict$Cohort <- 'Predict'

data <- rbind(twins, predict) %>% arrange(iid)

## Keep only those IDs that have both BSS data and MG data
common_ids <- intersect(data$iid, bss_df$ParticipantID) # N = 690
print(paste0(nrow(data) - length(common_ids), ' IDs do not have BSS')) # 315


##### Keep closest replicates in BSS df

## Format date
data$dov <- as.Date(data$dov, format = '%Y-%m-%d') 
bss_df$`Response Date` <- as.Date(bss_df$`Response Date`, format = '%Y-%m-%d') 
bss_df$ParticipantID <- as.integer(bss_df$ParticipantID)

## Add index column to pull out index locations
data['index'] <- 1:nrow(data)
bss_df['index'] <- 1:nrow(bss_df)
bss_df_ixs <- c() # to store index locations

for (i in common_ids) {
  dov <- data$dov[data$iid == i]
  bdf <- bss_df[bss_df$ParticipantID == i, ]
  min_diff <- min(abs(bdf$`Response Date` - dov))
  min_ix <- which(abs((bdf$`Response Date` - dov)) == min_diff)
  bdf_ix <- bdf$index[min_ix]
  bss_df_ixs <- c(bss_df_ixs, bdf_ix)
}
bss_df1 <- bss_df[bss_df$index %in% bss_df_ixs, ] %>% arrange(ParticipantID)

bs_data <- merge(data, bss_df1, by.x = 'iid', by.y = 'ParticipantID')
bs_data <- bs_data %>% select(-c('index.x', 'index.y')) # Index columns no longer required

## Calculate days between MS sample and BSS response date
bs_data$diffs <- bs_data$`Response Date` - bs_data$dov 

## Keep twins only
bs_twins <- bs_data[bs_data$Cohort == 'Twins', ]

sum(abs(bs_twins$diffs) < 28) # 539
sum(abs(bs_twins$diffs) < 14) # 509
sum(abs(bs_twins$diffs) < 10) # 481
sum(abs(bs_twins$diffs) < 7) # 451

## Keep those with BSS within two weeks before or after MG dov
bs_twins <- bs_twins[abs(bs_twins$diffs) < 14, ] # N = 509 within 2 weeks

#write.csv(bs_twins, "BSS_14days_twins.csv")

