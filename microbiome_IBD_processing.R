library(dplyr)  # 1.1.4

# CODE EXPLANATION --------------------------------------------------------

# Uses questionniare data to identify individuals suspected of having IBD to 
# be removed from the main analysis


# DATA PROCESSING ---------------------------------------------------------

## Get all IDs with microbiome data
twins <- read.csv("merged_metaphlan_bugs_abundance_table_20211201_version_4.tsv", sep = '\t')
predict <- read.csv("PREDICT1_metaphlan-4.beta.2_vJan21_CHOCOPhlAnSGB_202103.tsv", sep = '\t')

microbiome_ids <- unique(c(sub("^.(.*?)_.*", "\\1", paste0(names(twins)[-1])), sub(".*predict", "", paste0(names(predict)[-1]))))

## Load IBD responses
ibd <- read.csv("datasets/IBD/survey_responses.csv")

## Drop useless columns
ibd <- ibd[!grepl("What is/are the reason", ibd$Phenotype..Column.Heading.), ]

length(unique(ibd$ParticipantID))
#10756

ibd <- ibd[ibd$ParticipantID %in% microbiome_ids, ]

length(unique(ibd$ParticipantID))
#2439

ibd$Response.Date <- as.Date(ibd$Response.Date, format = "%d/%m/%Y")

table(ibd$Response)

'''
     0      1      4     10     13     15     17     18     22     23     24     25     26     27     29     30     31     32     33 
 10241    289      1      2      1      1      3      1      2      2      5      9      4      3      1     10      1      1      1 
    34     35     36     37     39     40     41     42     44     45     46     47     48     49     50     51     52     53     54 
     2      5      2      2      4     11      3      4      2     15      2      3      2      2      7      1      3      3      1 
    55     57     59     60     61     62     63     65     66     67     68     70     74     76     80 999902 999905 999906 999911 
     8      2      2      8      1      4      1      1      2      2      1      4      1      1      1     35      1  12440  38535
'''

## Drop the age at diagnosis questions
ibd1 <- ibd[!grepl('age', ibd$Phenotype..Column.Heading., ignore.case = T), ]

table(ibd1$Response)
distinct(ibd1[, c('Response', 'Response.Code.Description')])

"""
     0      1 999902 999905 999906 999911 
 10241    289     35      1   8285  22473 
  
  Response                                       Response.Code.Description
1   999911             Question skipped due to the questionnaire branching
2   999906                                               No data available
3        1                                                             Yes
4        0                                                              No
5   999911                                                  No Description
6   999906                                                  No Description
7   999902                                                      Don't know
8   999905 Unable to determine (unclear answered given by the participant)
"""

# Get total number of unique one responses
ibd_positive_N <- length(unique(ibd1$ParticipantID[ibd1$Response == 1]))
# 64 individuals report having IBD in at least one time point 

## Investigate these cases more closely
ibd_positive_IDs <- ibd1$ParticipantID[ibd1$Response == 1]
ibd_positive <- ibd[ibd$ParticipantID %in% ibd_positive_IDs, ] %>% arrange(ParticipantID, Response.Date)

ibd_positive_predict_ids <- ibd_positive$ParticipantID[grepl('PREDICT', ibd_positive$Survey.Version)]
ibd_positive_predict <- ibd_positive[ibd_positive$ParticipantID %in% ibd_positive_predict_ids, ]

## Drop suspected mistaken IDs that are also in our dataset
to_drop <- c(9412, 95611, 23222, 73991)
ibd_positive <- ibd_positive[!ibd_positive$ParticipantID %in% to_drop, ]

#write.csv(unique(ibd_positive$ParticipantID), "microiome_IBD_iids.csv")

















