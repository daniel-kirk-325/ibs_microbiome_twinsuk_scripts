import pandas as pd # 2.1.4
import sys
import os
import numpy as np # 1.26.4


'''

Explanation:

    Important note: The viral data used in this file is NOT the that was used in the Gut manuscript.
    
    This is old data that was since revised for the manuscript. 
    
    However, the iids did not change. 
    
    Hence, when retaining only those iids with viral data, I continued using the viroprofiler_iids.csv file generated 
    at the end of this script.
    
    Since the same iids underwent viral profiling with both techniques, the same results should be obtained
    if this list of iids was also generated using the R script (viral_data_processing.R)
    
    However, I am placing this script here and the associated data purely for reproducibility purposes. 

'''

## Data file for mapping GSIDs to metadata
mb_meta_data = pd.read_csv(r'datasets/raw_microbiome_data/list of ids with metagenomes and date of visit.txt', sep = '\t')
mb_meta_data = mb_meta_data.rename(columns = {'gsid':'GSID'})

## Viral TPM normalised data
b1 = pd.read_csv(r"data_deposition/gut_paper/gut_paper_datasets/metagenome_data/phage_data/old_phage_data/tpm_merged_filtered_b1.tsv", delimiter='\t')
b2 = pd.read_csv(r"data_deposition/gut_paper/gut_paper_datasets/metagenome_data/phage_data/old_phage_data/tpm_merged_filtered_b2.tsv", delimiter='\t')
b3 = pd.read_csv(r"data_deposition/gut_paper/gut_paper_datasets/metagenome_data/phage_data/old_phage_data/tpm_merged_filtered_b3.tsv", delimiter='\t')
b4 = pd.read_csv(r"data_deposition/gut_paper/gut_paper_datasets/metagenome_data/phage_data/old_phage_data/tpm_merged_filtered_b4.tsv", delimiter='\t')
b5 = pd.read_csv(r"data_deposition/gut_paper/gut_paper_datasets/metagenome_data/phage_data/old_phage_data/tpm_merged_filtered_b5.tsv", delimiter='\t')
b6 = pd.read_csv(r"data_deposition/gut_paper/gut_paper_datasets/metagenome_data/phage_data/old_phage_data/tpm_merged_filtered_b6.tsv", delimiter='\t')


## Remove columns which sum to 0 
dfs = [b1, b2, b3, b4, b5, b6] 
names = ["b1", "b2", "b3", "b4", "b5", "b6"] 
def drop_zero_sum_cols(df, name): 
    data_cols = df.iloc[:, 1:] 
    zero_sum_cols = data_cols.columns[data_cols.sum() == 0] 
    print(f"{name}: removed {len(zero_sum_cols)} columns") 
    return df.drop(columns=zero_sum_cols) 
b1, b2, b3, b4, b5, b6 = [drop_zero_sum_cols(df, name) for df, name in zip(dfs, names)]


all_batches = pd.concat([pd.Series(b1.columns[1:]), pd.Series(b2.columns[1:]), pd.Series(b3.columns[1:]), pd.Series(b4.columns[1:]), pd.Series(b5.columns[1:]), pd.Series(b6.columns[1:])], ignore_index=True)

## Obtain GSID from column name for each participant
col_name, gsids = [], []
for j in all_batches:
    col_name.append(j)
    ix_start = j.find('_103667') + 1
    ix_end = ix_start + len('103667-004-XXX')
    gsids.append(j[ix_start:ix_end])
    
batch_nums = np.repeat(1, len(b1.columns[1:])).tolist() + np.repeat(2, len(b2.columns[1:])).tolist() + np.repeat(3, len(b3.columns[1:])).tolist() + np.repeat(4, len(b4.columns[1:])).tolist() + np.repeat(5, len(b5.columns[1:])).tolist() + np.repeat(6, len(b6.columns[1:])).tolist()
gsids_df = pd.DataFrame({'Batch':batch_nums, 'GSID':gsids, 'col_name':col_name})

len(gsids_df) # 5528
len(gsids_df.GSID.unique()) # 2480

merged_df = pd.merge(mb_meta_data, gsids_df, on = 'GSID')
len(merged_df) # 4685
len(merged_df.GSID.unique()) # 2068
len(merged_df.iid.unique()) # 1678
len(merged_df.drop_duplicates(['GSID', 'batch'])) # 2068
len(merged_df.drop_duplicates(['GSID', 'dov'])) # 2068

## Duplicate GSIDs were the same sample run in different lanes. 
## They do not differ by characteristic (confirmed by Sherry) and 
## one can be selected at random

merged_df = merged_df.drop_duplicates('GSID')

## Add phage data
phage_cols = merged_df.col_name
phage_samples = pd.concat([b1['ID'], 
                           b1.loc[:, b1.columns.isin(phage_cols)], 
                           b2.loc[:, b2.columns.isin(phage_cols)],
                           b3.loc[:, b3.columns.isin(phage_cols)],
                           b4.loc[:, b4.columns.isin(phage_cols)],
                           b5.loc[:, b5.columns.isin(phage_cols)]], axis = 1).T



# Set column names to first row and drop first row
phage_samples.columns = phage_samples.iloc[0].values
phage_samples = phage_samples.iloc[1:, :]

# Reset index and rename to merge with meta data
phage_samples = phage_samples.reset_index().rename(columns={'index':'col_name'})

phage_df_raw = pd.merge(merged_df, phage_samples, on = 'col_name')

# Drop column with name asterisk, as suggested in the readme
phage_df_raw = phage_df_raw.drop('*', axis = 1)

print(phage_df_raw.shape)
print(phage_df_raw.prefix.nunique())
print(phage_df_raw.GSID.nunique())
print(phage_df_raw.col_name.nunique())

# Save IDs only
phage_df_raw['iid'].to_csv(r"data_deposition/gut_paper/gut_paper_datasets/metagenome_data/phage_data/viroprofiler_iids.csv")





































