library(readxl)  # 1.4.5
library(stringr) # 1.6.0
library(dplyr)   # 1.1.4


# CODE EXPLANATION --------------------------------------------------------

## Purpose

# Takes the pre-processed PREDICT-1 metagenomic dataset and applies processing
# based on medication use

## Data pre-processing

# 1) Removes missing values to be re-added later
# 2) Keep only the nearest time point between drug questionniare and IBS&metagenome sampling

## Drug pre-processing 

# 3) Gets all unique classes and drugs
# 4) Assigns each drug its appropriate class 

## Get drug counts 

# 5) Count participants using each class of drugs
# 6) Remove antibiotic users
# 7) Keep only those drugs that are adjusted for in the main analysis


# DATA PRE-PROCESSING  -----------------------------------------------------

## Load PREDICT-1 metagenomic data from file predict_metagenomic_processing.py
df <- read.csv("ibs_microbiome_predict.csv")[, -1]
df$dov <- as.Date(df$dov)

## Load drug data
drugs <- read_excel("TwinsUK_Prescribed_Medication_V1.0_2021.xlsx")
drugs <- drugs[drugs$ParticipantID %in% df$iid, ]

## Exclude missing values 
missing <- df[!df$iid %in% drugs$ParticipantID, ]
completecases <- df[df$iid %in% drugs$ParticipantID, ]

## Get nearest time point between drugs and df
drugs['index'] <- 1:nrow(drugs)
ixs <- c()
date_diffs <- c()
for (i in completecases$iid){
  dt <- as.Date(completecases$dov[completecases$iid == i])
  dd <- drugs[drugs$ParticipantID == i, ]
  date_diff <- abs(as.Date(dd$Date) - dt)
  date_diffs <- c(date_diffs, min(date_diff))
  ix <- which.min(date_diff)
  ixs <- c(ixs, dd$index[ix])
}

drugs_fil <- drugs[drugs$index %in% ixs, ]

# See differences in dates
plot(table(round((date_diffs/365.25),4)))

n <- 365*5
length(date_diffs[date_diffs < n])
length(date_diffs[date_diffs < n]) / nrow(completecases)

## Sort dfs to make binding easier 

drugs_fil <- arrange(drugs_fil, ParticipantID)
completecases <- arrange(completecases, iid)


# DRUG PRE-PROCESSING -----------------------------------------------------

## Get each class and drug
all_classes <- unique(unlist(strsplit(as.character(drugs_fil$Class), ";"), recursive = TRUE))
all_drugs1 <- unique(unlist(strsplit(as.character(drugs_fil$Substance), ";"), recursive = TRUE))
all_drugs2 <- unique(unlist(strsplit(as.character(drugs_fil$Product), ";"), recursive = TRUE))
all_drugs <- unique(c(all_drugs1, all_drugs2))

## Make variable that stores the class of each drug

# For substance And for product 
keys <- unlist(strsplit(as.character(drugs_fil$Class), ";"), recursive = TRUE)
items_sub <- unlist(strsplit(as.character(drugs_fil$Substance), ";"), recursive = TRUE)
items_prod <- unlist(strsplit(as.character(drugs_fil$Product), ";"), recursive = TRUE)


# GET DRUG COUNTS  -----------------------------------------------------

## List of drug classes relevant to IBS (excluding antibiotics)

ibs_classes <- drug_classes <- c(
  "Proton Pump Inhibitors",
  "Tricyclic & Related Antidepressant Drugs",
  "Antispasmod.&Other Drgs Alt.Gut Motility",
  "Antimotility Drugs",
  "Selective Serotonin Re-Uptake Inhibitors",
  "Other Drugs Used In Constipation",
  "Non-Opioid Analgesics And Compound Prep",
  "Opioid Analgesics",
  "Drugs Affecting Biliary Composition&Flow",
  "Drugs Used In Nausea And Vertigo", 
  'Bulk-Forming Laxatives', 
  'Soothing Haemorrhoidal Preparations', 
  'Stimulant Laxatives', 
  'Osmotic Laxatives', 
  'Anxiolytics',
  "Thyroid Hormones",
  "Oestrogens And Hrt",
  "Male Sex Hormones And Antagonists",
  "Parenteral Progestogen-only Contracep",
  "Oral Progestogen-only Contraceptives",
  "Replacement Therapy",
  "Combined Hormonal Contraceptives",
  "Progestogens & Progesterone Receptor Mod",
  "Progestogens"
)

class_counts_df <- data.frame(matrix(0, nrow = nrow(completecases), ncol = length(ibs_classes)))
colnames(class_counts_df) <- ibs_classes
class_counts_df$ParticipantID <- completecases$iid

for (C in ibs_classes) {
  ids_to_one <- drugs_fil$ParticipantID[grepl(C, drugs_fil$Class, fixed = T)]
  class_counts_df[class_counts_df$ParticipantID %in% ids_to_one, C] <- 1
}


#### Merge similar drug classes 

merge_cols <- function(data, col1_name, col2_name, new_col_name) {
  # Merge the two columns and replace values greater than 0 with 1
  merged_col <- pmax(data[[col1_name]], data[[col2_name]], na.rm = TRUE)
  merged_col <- ifelse(merged_col > 0, 1, 0)
  
  # Drop the original columns
  data <- data[, !(names(data) %in% c(col1_name, col2_name))]
  
  # Add the merged column to the dataframe
  data[[new_col_name]] <- merged_col
  
  return(data)
}

# Merge SSRIs anxiolytics and other antidepressant drugs
class_counts_df <- merge_cols(class_counts_df, "Tricyclic & Related Antidepressant Drugs", "Selective Serotonin Re-Uptake Inhibitors", "Anti_Depressants")
class_counts_df <- merge_cols(class_counts_df, "Anti_Depressants", "Anxiolytics", "Anti_Depressants&anxiolytics")

# Merge Antispasmod.&Other Drgs Alt.Gut Motility & Other Drugs Used In Constipation
class_counts_df <- merge_cols(class_counts_df, "Stimulant Laxatives", "Osmotic Laxatives", "Laxatives")
class_counts_df <- merge_cols(class_counts_df, "Laxatives", "Other Drugs Used In Constipation", "Anti_Constipation")
class_counts_df <- merge_cols(class_counts_df, "Anti_Constipation", "Bulk-Forming Laxatives", "Laxatives")

# Merge Antispasmod.&Other Drgs Alt.Gut Motility and Antimotility Drugs
class_counts_df <- merge_cols(class_counts_df, "Antispasmod.&Other Drgs Alt.Gut Motility", "Antimotility Drugs", "Antispasmodics&antimotility")

# Merge Analgesics 
class_counts_df <- merge_cols(class_counts_df, "Opioid Analgesics", "Non-Opioid Analgesics And Compound Prep", "Analgesics")

# Merge progestogens and progesterones 
class_counts_df <- merge_cols(class_counts_df, "Progestogens & Progesterone Receptor Mod", "Progestogens", "Progestogens")
class_counts_df <- merge_cols(class_counts_df, "Progestogens", "Oral Progestogen-only Contraceptives", "Progestogens")
class_counts_df <- merge_cols(class_counts_df, "Progestogens", "Parenteral Progestogen-only Contracep", "Progestogens")

# Merge Oestrogens And Hrt and Replacement Therapy
class_counts_df <- merge_cols(class_counts_df, "Oestrogens And Hrt", "Replacement Therapy", "HRT")


colSums(class_counts_df)


## Remove antibiotic users

abs <- c(
  'Cephalosporins',
  'Metronidazole',
  'Macrolides',
  'Tetracycline',
  'Quinolones',
  'Metronidazole, Tinidazole & Ornidazole',
  'Broad-Spectrum Penicillins',
  'Benzylpenicillin&Phenoxymethylpenicillin',
  'Penicillinase-Resistant Penicillins', 
  'Some Other Antibacterials')

abs_iid <- c()
for (i in abs) {
  abs_iid <- c(abs_iid, drugs_fil$ParticipantID[grepl(i, drugs_fil$Class)])
}

ab_df <- drugs_fil[drugs_fil$ParticipantID %in% unique(abs_iid), ]

abs_prod <- unlist(strsplit(as.character(ab_df$Product), ";"))
abs_sub <- unlist(strsplit(as.character(ab_df$Substance), ";"))
abs_class <- unlist(strsplit(as.character(ab_df$Class), ";"))

ixs <- which(abs_class %in% abs)

abs_prod <- abs_prod[ixs]
abs_sub <- abs_sub[ixs]
abs_class <- abs_class[ixs]

unique_abs <- unique(c(abs_prod, abs_sub))

## Get the participants using antibiotics 

ids_using_abs <- c()

for (ab in unique_abs) {
  ids_using_abs <- c(ids_using_abs, drugs_fil$ParticipantID[grepl(ab, drugs_fil$Substance, fixed = T) | grepl(ab, drugs_fil$Product, fixed = T)])
}

ids_using_abs <- unique(ids_using_abs)

## Recode these IDs to ones in class_counts_df

class_counts_df$antibiotics <- 0
class_counts_df[class_counts_df$ParticipantID %in% ids_using_abs, 'antibiotics'] <- 1

## Remove those using antibiotics
class_counts_df<- class_counts_df[class_counts_df$antibiotics == 0, ]

## Merge with complete cases

# Keep only those drugs that are adjusted for in the main analysis 
keep_cols <- c("Proton Pump Inhibitors", "Anti_Depressants&anxiolytics", "Antispasmodics&antimotility", "Analgesics", 'Laxatives')

completecases_plus_drug <- merge(completecases, class_counts_df[, c('ParticipantID', keep_cols)], by.x = 'iid', by.y = 'ParticipantID')
# write.csv(completecases_plus_drug, "ibs_microbiome_predict_drug_processed.csv")










