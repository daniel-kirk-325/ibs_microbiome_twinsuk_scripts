import os
import pandas as pd # 2.1.4
import numpy as np  # 1.26.4
import sys

'''
1) Load microbiome data
2) Process metadata
3) Add IBS data
4) Keep only observations with IBS data <1 year after metabolomics/microbiome 
    sampling and <3 years before
5) Keep true controls only 
'''

#%% Load functions

def least_diff(data1, datecol1, idcol1, data2, datecol2, idcol2, diff_thresh, date_time = False):    
    '''
    Description:
        
        Looks through all combinations of each ID and retuns the index combinations
        with the smallest difference
        
        Notes: 
            
            Convert datecol to date.time before using the function
            
            Remember to adjust diff_threshold according to the format of your
            date column, whether full date, year, or age.
            
            The function will only return the minimum difference for a given combination 
            of dates for two datasets. If you would like all the differences below
            diff_thresh, modify the last two "if" clases in PartB to select
            all those lower than diff_thresh, not just the minimum one
        
    ----------
    
    Parameters:
        - data1 & data2 are the dataset
        
        -datecol1 & datecol2 are the names of the columns that denote sampling date
        in each of the datasets 
        
        -idcol1 & idcol2 are the names of the ID columns in each dataset
                
        -diff_thresh is the threshold by which combinations will be retained. 
        Pairs with a difference about this threshold will not be returned from 
        the function
        
        - date_time should be True if date column is in datetime format
    ----------
    Returns:
        
        - Two lists that contain the indices of the datapoints  with the smallest 
        difference in time per ID, all of which are less than diff_thresh
    -------

    '''
    
   # ========================================================================== 
   # A) Indices that appear only once in each dataset 
   # ========================================================================== 
   
    idcol1_counts = data1[idcol1].value_counts()
    idcol2_counts = data2[idcol2].value_counts()
       
    mask1 = data1[idcol1].isin(idcol1_counts[idcol1_counts == 1].index)
    mask2 = data2[idcol2].isin(idcol2_counts[idcol2_counts == 1].index)
       
    unq1 = data1[mask1]
    unq2 = data2[mask2]
       
    commonID = list(set(unq1[idcol1]).intersection(set(unq2[idcol2])))
       
    unq1, unq2 = unq1[unq1[idcol1].isin(commonID)], unq2[unq2[idcol2].isin(commonID)]
              
    ixa1, ixa2, diff = [], [], []
    for i in unq1[idcol1]:
        time_diff = abs(unq1[datecol1][unq1[idcol1] == i].reset_index(drop=True) - unq2[datecol2][unq2[idcol2] == i].reset_index(drop=True))
        if date_time:
            time_diff = time_diff.dt.days.values[0]
            if time_diff < diff_thresh:
                ixa1.append(unq1[unq1[idcol1] == i].index.values[0])
                ixa2.append(unq2[unq2[idcol2] == i].index.values[0])
                diff.append(time_diff)
        else:
            if (time_diff < diff_thresh).values[0]:
                ixa1.append(unq1[unq1[idcol1] == i].index.values[0])
                ixa2.append(unq2[unq2[idcol2] == i].index.values[0])
                diff.append(time_diff)
        

       
    # ========================================================================== 
    # B) Indices that appear multiple times 
    # ========================================================================== 
       
    nunq1 = data1[~data1[idcol1].isin(commonID)]
    nunq2 = data2[~data2[idcol2].isin(commonID)]
    nunq1 = nunq1[nunq1[idcol1].isin(nunq2[idcol2])]
    nunq2 = nunq2[nunq2[idcol2].isin(nunq1[idcol1])]                  
    
    ixb1, ixb2, diffs = [], [], []
    for i in nunq1[idcol1].unique():
        df1 = nunq1[nunq1[idcol1] == i]
        df2 = nunq2[nunq2[idcol2] == i]   
        ib1, ib2, ds = [], [], []
        for d1 in df1.index:
            for d2 in df2.index:
                ib1.append(d1)
                ib2.append(d2)
                d = abs(df1.loc[d1][datecol1] - df2.loc[d2][datecol2])
                ds.append(d)
        mindiff = min(ds)
        if date_time:
            mindiff = mindiff.days       
            
            if mindiff < diff_thresh:
                ixx = ds.index(min(ds))
                ixb1.append(ib1[ixx])
                ixb2.append(ib2[ixx])
                diffs.append(mindiff)
        else:
            if (mindiff < diff_thresh):
                ixx = ds.index(min(ds))
                ixb1.append(ib1[ixx])
                ixb2.append(ib2[ixx])
                diffs.append(mindiff)
        
    return ixa1+ixb1, ixa2+ixb2, diff+diffs


#%% 0) Load IBS data from file IBS_data_processing.py

## Load IBS metadata 
ibs = pd.read_csv(r"IBS_dataset.csv").iloc[:, 1:]
ibs.ResponseDate = pd.to_datetime(ibs.ResponseDate)
## Note IBS incident cases
inc_cases = ibs[ibs.incident_ibs == 1]['ParticipantID']


#%% 1) Load microbiome data

mb = pd.read_csv(r'merged_metaphlan_bugs_abundance_table_20211201_version_4.tsv', sep = '\t')
dates = pd.read_csv(r'list of ids with metagenomes and date of visit.txt', sep = '\t')
pred = pd.read_csv(r'PREDICT1_metaphlan-4.beta.2_vJan21_CHOCOPhlAnSGB_202103.tsv', sep = '\t')

## Only keep bacteria annotated at species level 
spec_level = mb.loc[(mb['clade_name'].str.contains("s__")) & (~mb['clade_name'].str.contains("t__"))]
spec_name = spec_level['clade_name'].str.split('s__').str[-1]

mb = pd.concat([spec_name, spec_level.iloc[:, 1:]], axis = 1).rename(columns = {'clade_name':'spec_name'})
## Drop controls 
mb = mb.drop(['904591_T1a', '904591_T1b','904591_T1c','904591_T1d'], axis = 1)
dates = dates[~dates.prefix.isin(['904591_T1a', '904591_T1b','904591_T1c','904591_T1d'])]

## Transpose dataset, so that samples are rows and bacteria are columns
mbt = mb.T.reset_index().rename(columns = {'index':'prefix'})
print("Number of observations: ", len(mbt)-1)

sum(mbt.prefix.values[1:] != dates.prefix.values)

mbt['dov'] = ['dov'] + dates.dov.to_list()
mbt['iid'] = ['iid'] + dates.iid.to_list()
mbt['fid'] = ['fid'] + dates.iid.astype(str).str[:-1].to_list()
mbt['batch'] = ['batch'] + dates.batch.to_list()

## Update colnames
mbt = mbt.rename(columns = dict(zip(mbt.columns[1:], mbt.iloc[0, 1:])))
mbt = mbt[1:]

## Remove those in predict
pred_fids = pred.columns[1:].str[len('predict'):-1]
mbt = mbt[~mbt.fid.isin(pred_fids)]

#### Keep most recent date for each participant with multiple readings
mbt.dov = pd.to_datetime(mbt.dov, format = '%d/%m/%Y')

multiplicated_ids = mbt.iid.value_counts()[mbt.iid.value_counts() > 1].index

mbt = mbt.loc[mbt.groupby('iid')['dov'].idxmax()]


#%% 2) Add meta data

## Fill in missing ages & sex
meta = pd.read_excel(r"bmi_2020.xlsx", sheet_name = 'TwinDetails')
meta = meta[~meta.YEAR_BIRTH.isna()]

mbt = pd.merge(mbt, meta[['STUDY_NO', 'YEAR_BIRTH', 'SEX']], left_on='iid', right_on='STUDY_NO', how='inner').drop('STUDY_NO', axis = 1)
mbt['YEAR_BIRTH'] = pd.to_datetime(mbt['YEAR_BIRTH'].astype(int).astype(str) + '-07-01')
mbt['age_Mb_sample'] = round((mbt.dov - mbt.YEAR_BIRTH)  / pd.Timedelta(days=365.25), 2)

## Fill in missing BMI
bmidf = pd.read_excel(r"bmi_2020.xlsx", sheet_name = 'bmi')

data1ix, data2ix, diffs = least_diff(mbt, 'dov', 'iid', bmidf, 'Dates', 'TwinsID', diff_thresh=365*5, date_time=True)
mbt, bmi_data = mbt.loc[data1ix], bmidf.loc[data2ix]

## 16 did not have  BMI data within 5 years 
mbt['bmi_Mb_sample'] =  bmidf.loc[data2ix, 'BMI'].values
print("Len mbt_samples = {}".format(len(mbt)))

#%% 3) Add IBS data

overlap = list(set(ibs.ParticipantID).intersection(set(mbt.iid)))
print(f"{len(overlap)} participants have IBS data")

ibs = ibs[ibs.ParticipantID.isin(overlap)]
mbt = mbt[mbt.iid.isin(overlap)]

#%% 4) Keep only observations with IBS data <1 year after metabolomics/microbiome 
     # sampling and <3 years before

## Get IBS status at most recent response date before metab sampling
ibs_status_last_date_response = {}
ibs_status_last_date_date = {}

## Get index locations of IBS responses before each metabs reading per participant
pre_metabs_ixs, most_recent_observation_ix = [], []

## Track participants who have metabolomics data before IBS response date
## and get the number of days of this difference
metabs_before_ibs, date_diffs = [], []

for i in mbt.iid:
    date = pd.to_datetime(mbt[mbt.iid == i]['dov'].values[0]) + np.timedelta64(365, 'D')
    df = ibs.loc[(ibs.ParticipantID == i)]
    
    if np.any(df.ResponseDate <= date):
        df_pre_metabs = df.loc[df.ResponseDate <= date]
        
        # Get most recent response
        ibs_status_last_date_response[i] = df_pre_metabs['ibs_status_upd'].values[-1]
        
        # Get most recent date
        ibs_status_last_date_date[i] = df_pre_metabs['ResponseDate'].values[-1]
        
        # Get ix locations 
        pre_metabs_ixs.extend(df_pre_metabs.index)
        
        # Get ix location of most recent observation
        most_recent_observation_ix.append(df_pre_metabs.index[-1])

    else:
        # Otherwise there is no IBS response before most recent metabs reading 
        # so note the difference in days 
        metabs_before_ibs.append(i)
        date_diffs.append((df.ResponseDate.min() - date).days)
        

print(f"{len(metabs_before_ibs)} participants did not have IBS data before most recent sampling date")

mbt_pre_metabs = mbt[mbt.iid.isin(metabs_before_ibs)]

## Drop participants without IBS information prior to sampling
mbt_mrd = pd.merge(mbt, ibs.loc[most_recent_observation_ix][['ParticipantID', 'ResponseDate', 'ibs_status_upd']].rename(columns = {'ParticipantID':'iid'}), on = 'iid')

#### Get ID, age at metabs visit, IBSRomIII, IBS overall of the participants 
#### in mbt_mrd
table = mbt[mbt.iid.isin(ibs_status_last_date_date.keys())][['iid', 'age_Mb_sample']]

# Make df containing only observations where IBS came before metabolomics 
# sampling per participant
pre_metabs_df = ibs.loc[pre_metabs_ixs]
pre_metabs_df = pre_metabs_df[pre_metabs_df.ParticipantID.isin(mbt_mrd.iid)] 


## Get all Rome III responses
pmd_romeIII = pre_metabs_df[(pre_metabs_df.Origin == 'Rome_III_PH_codes') | (pre_metabs_df.Origin == 'Rome_III_Q18')]


## If any Rome III response before metabs sampling date == 1, participant 
## is given response code 1
rome_responses = []
for i in pmd_romeIII.ParticipantID.unique():
    date = mbt[mbt.iid == i]['dov'].values[0]
    
    if np.any(pmd_romeIII[(pmd_romeIII.ParticipantID == i) & (pmd_romeIII.ResponseDate <= date)]['ibs_status'] == 1):
        rome_responses.append(1)
    
    else:
        rome_responses.append(0)

rome_responses = pd.DataFrame({'iid':pmd_romeIII.ParticipantID.unique(),'IBS_RomeIII': rome_responses})

table = pd.merge(table, rome_responses, on = 'iid', how = 'outer')

ibs_status_last_date_response_series = pd.Series(ibs_status_last_date_response)
table['IBS_overall'] = ibs_status_last_date_response_series.values

table.IBS_overall.value_counts(True)

mbt_mrd['sampling_diff'] = abs((mbt_mrd['dov'] - mbt_mrd['ResponseDate']).dt.days / 365.25)
print(f"{len(mbt_mrd['sampling_diff'][mbt_mrd['sampling_diff'] > 3])} participants had a sampling difference >3 years and were removed")
# Drop participants with sampling diff > 3
mbt_mrd = mbt_mrd[mbt_mrd.sampling_diff < 3]

table = table[table.iid.isin(mbt_mrd.iid)]

#%% 5) Controls - Only include those that have answered a Rome III questionnaire 
#   and do NOT report to have IBS

controls = table[(table.IBS_RomeIII == 0) & (table.IBS_overall == 0)]
cases = table[(table.IBS_RomeIII == 1)]
cc = pd.concat([controls, cases], axis = 0)

mrd_cc = mbt_mrd[mbt_mrd.iid.isin(list(controls.iid) + list(cases.iid))]

mrd_cc = pd.merge(mrd_cc, cc[['iid', 'IBS_overall']], on = 'iid')

# Drop those who develop IBS in the future
inc_cases = ibs[ibs.incident_ibs == 1]['ParticipantID'].unique()
mrd_cc = mrd_cc[~((mrd_cc.iid.isin(inc_cases)) & (mrd_cc.IBS_overall == 0))]

# Remove IBD cases
ibd = pd.read_csv(r"microiome_IBD_iids.csv")
mrd_cc = mrd_cc[~mrd_cc.iid.isin(ibd.x)]

#mrd_cc.to_csv(r"ibs_microbiome_twins.csv")

