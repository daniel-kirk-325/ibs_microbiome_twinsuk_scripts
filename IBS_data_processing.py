# Modules & functions 
import pandas as pd # 2.1.4
import numpy as np  # 1.26.4
import os

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


#%% Import data 

# Year of birth
yob = pd.read_excel(r'ibs_data\ibs_data_request\DanielKirk_18.08.2023.xlsx', sheet_name = 'TwinDetails')
yob = yob.rename(columns = {'STUDY_NO':'ParticipantID'})

# Ever been told by doctor had IBS 
df1a = pd.read_excel(r'ibs_data\ibs_data_request\DanielKirk_18.08.2023.xlsx', sheet_name = 'Q11A')
# Has a doctor EVER diagnosed or treated you for any of the following conditions? \ Irritable Bowel Syndrome 
# PH0000536 - PH0000537
df1b = pd.read_excel(r'ibs_data\ibs_data_request\DanielKirk_18.08.2023.xlsx', sheet_name = 'Q22')
df1c = pd.read_excel(r'ibs_data\ibs_data_request\DanielKirk_18.08.2023.xlsx', sheet_name = 'Q23_Q34')
# Q29 - Has a doctor EVER told you that you have/had ANY of the conditions listed below in questions 2-14?
df1d = pd.read_excel(r'ibs_data\ibs_data_request\DanielKirk_18.08.2023.xlsx', sheet_name = 'Q29')
# Q36 - Has a doctor ever told you that you have/had ANY of the conditions listed below in questions 2-14?
df1e = pd.read_excel(r'ibs_data\ibs_data_request\DanielKirk_18.08.2023.xlsx', sheet_name = 'Q36')

 
#Please tick the most appropriate box.\Have you ever been diagnosed with:\a. Irritable bowel syndrome (IBS)
df3a = pd.read_excel(r'ibs_data\ibs_data_request\DanielKirk_2_01.09.2023.xlsx', sheet_name = 'Q18')
# Has a doctor ever diagnosed or treated you for any of the following conditions? \ Irritable Bowel Syndrome 
df3b = pd.read_excel(r'ibs_data\ibs_data_request\DanielKirk_2_01.09.2023.xlsx', sheet_name = 'Q21_mod')
# At what age did you develop IBS?
df3c = pd.read_excel(r'ibs_data\ibs_data_request\DanielKirk_2_01.09.2023.xlsx', sheet_name = 'Q22')
# Are you currently suffering from any of the following gastrointestinal illnesses?/ Irritable Bowel Syndrome
df3d = pd.read_excel(r'ibs_data\ibs_data_request\DanielKirk_2_01.09.2023.xlsx', sheet_name = 'Q27')

# Questions for Rome 3 
rc = pd.read_excel(r'ibs_data\ibs_data_request\DanielKirk_27.10.2023.xlsx', sheet_name = 'Data')
Q18 = pd.read_excel(r"ibs_data\ibs_data_request\DanielKirk_05.06.2024.xlsx", sheet_name = 'Q18').iloc[1:, :]
Q18_124 = pd.read_excel(r"ibs_data\ibs_data_request\DanielKirk_21.05.2024.xlsx", sheet_name = 'Q18').iloc[1:, :][['ParticipantID', 'Q18_124']]
Q18 = pd.merge(Q18, Q18_124, on = 'ParticipantID', how = 'outer')
Q18_115_118 = pd.read_excel(r"ibs_data\ibs_data_request\DanielKirk_12.06.2024.xlsx", sheet_name = 'Q18').iloc[1:, :]
Q18 = pd.merge(Q18_115_118, Q18, on = 'ParticipantID', how = 'outer')
# Get date for Q18
Q18_response_date = pd.read_csv(r"ibs_data\ibs_data_request\Q18_Autumn2009Questionnaire_v1_Date.csv")
Q18 = pd.merge(Q18, Q18_response_date, on = 'ParticipantID', how = 'outer')
Q18.ResponseDate = pd.to_datetime(Q18.ResponseDate)
    
#%% Cleaning 

# NOTE for df3b: *Where no response date is available, the questionnaire was paper based and 
# completed between November 2010 and April 2011

'''
Datasets containing information about whether participants have ever been told they have IBS
are combined. 

This information is used later before analysis to determine true controls (i.e., those who are negatve
according to Rome III criteria and never self-report IBS)

'''

## Print colnames 

df_list = [df1a, df1b, df1c, df1d, df1e, df3b, df3c, df3d]

for i, df in enumerate(df_list):
    print(i)
    print(df.columns, '\n')


## Remove question descriptions 

df1a = df1a.iloc[1:, :]
df1b = df1b.iloc[1:, :]
df1d = df1d.iloc[1:, :]
df1e = df1e.iloc[1:, :]
df3a = df3a.iloc[1:, :]
df3c = df3c.iloc[1:, :]
df3c.PH0001248 = df3c.PH0001248.astype(int)

## Extract relevant information as columns in remaining dataframes 

df1c_rel = df1c[['PhenotypeID', 'ParticipantID', 'ResponseDate', 'ResponseCode', 'ResponseDescription']]
df3d_rel = df3d[['PhenotypeID', 'ParticipantID', 'ResponseDate', 'ResponseCode', 'ResponseDescription']]

## Change "Study_No" to "ParticipantID" for unity 

df1b = df1b.rename(columns = {'Study_No' : 'ParticipantID'}) 
df1c_rel = df1c_rel.rename(columns = {'Study_No' : 'ParticipantID'}) 


#%% Investigating PH0000536, PH0000537, PH0000538

'''
df1b & df1c = Has a doctor EVER diagnosed or treated you for any of the following conditions? \ Irritable Bowel Syndrome 
df1d & df1e = Has a doctor EVER diagnosed or treated you for any of the following conditions? \ Irritable Bowel Syndrome 
'''

## Vertically combine df1d and df1e and keep only the relevant cols
df1de = pd.concat([df1d, df1e])[['ParticipantID', 'ResponseDate', 'PH0000536', 'PH0000537', 'PH0000538']]
df1de.iloc[:, -3:] = df1de.iloc[:, -3:].astype(int)
df1de['Origin'] = '1de'

## Separate df1c into two dataframes to capture IBS status and age of diagnosis 
df1c_p36 = df1c_rel[df1c_rel.PhenotypeID == 'PH0000536'][['ParticipantID', 'ResponseDate', 'ResponseCode']].rename(columns = {'ResponseCode':'PH0000536'})
df1c_p37 = df1c_rel[df1c_rel.PhenotypeID == 'PH0000537'][['ParticipantID', 'ResponseDate', 'ResponseCode']].rename(columns = {'ResponseCode':'PH0000537'})
df1c_p36['Origin'] = '1c'

df3a_p36 = df3a.rename(columns = {'Q18_143':'PH0000536', 'ResponseYear':'ResponseDate'})
df3a_p36['Origin'] = '3a'
## This was autunm questionniare 
df3a_p36.ResponseDate = pd.to_datetime('2009-11-01')

## Prepare df3d_rel for vertical concatenation below 
df3d_p36 = df3d_rel[['ParticipantID', 'ResponseDate', 'ResponseCode']].rename(columns = {'ResponseCode':'PH0000536'})
df3d_p36['Origin'] = '3d'


df3b.ResponseDate[df3b.ResponseDate.isnull()] = '2011-02-01'
df3b.ResponseDate = pd.to_datetime(df3b.ResponseDate)
df3b = df3b.rename(columns = {'Q21_93': 'PH0000536'})
df3b['Origin'] = '3b'


## Check for presence of IBS over time using PH0000536
ph36 = pd.concat([df1c_p36, df1de[['ParticipantID', 'ResponseDate', 'PH0000536', 'Origin']]])
ph36 = pd.concat([ph36, df3b])
ph36 = pd.concat([ph36, df3a_p36])
ph36 = pd.concat([ph36, df3d_p36])

ph36.PH0000536 = pd.to_numeric(ph36.PH0000536)
ph36 = ph36.sort_values(['ParticipantID', 'ResponseDate'])
ph36.PH0000536.unique() 

# 999906 = Don't know
# 999902 = Don't know 
# 999911 = Question skipped legitimately
# 999905 = Unable to determine (unclear answered given by the participant)

## 999911 can be converted to 0, the rest should be removed 

ph36.PH0000536[ph36.PH0000536 == 999911] = 0
drop_codes = [999902, 999906, 999905]
ph36 = ph36[~ph36.PH0000536.isin(drop_codes)]
ph36 = ph36.drop_duplicates(['ParticipantID', 'ResponseDate', 'PH0000536'])
ph36 = ph36.reset_index(drop = True)
# Keep a copy of original, unmodified ph36 df 
ph36_copy = ph36

# Make df with only those with all 0s or all 1s
all_zero_or_one_ids = ph36_copy.groupby('ParticipantID')['PH0000536'].apply(lambda x: np.all(x == 1) | np.all(x == 0)).index[ph36_copy.groupby('ParticipantID')['PH0000536'].apply(lambda x: np.all(x == 1) | np.all(x == 0))].to_list()
all_zero_or_one = ph36_copy[ph36_copy.ParticipantID.isin(all_zero_or_one_ids)]

#### Process age information 

# How many of these give at least one age of diagnosis 
df1de_age = df1de[df1de.PH0000537 < 1000][['ParticipantID', 'ResponseDate',  'PH0000537']].rename(columns = {'PH0000537':'ageD'})
df1a_age = df1a[df1a.Q11A_331.notna()][['ParticipantID', 'ResponseDate',  'Q11A_331']].rename(columns = {'Q11A_331':'ageD'})
df3c_age = df3c[df3c.PH0001248 < 1000]
df3c_age.columns = df1a_age.columns

## Combine dfs 
agel = pd.concat([pd.concat([df1de_age, df1a_age]), df3c_age]).reset_index(drop = True)
agel.ageD = agel.ageD.astype(int)

## When multiple ages are given, retain the lowest age 
lowest_ages = agel.groupby('ParticipantID')['ageD'].idxmin()
agel = agel.loc[lowest_ages]

# 1679 ages available

# Remove time from datecol 
agel['ResponseDate'] = agel['ResponseDate'].dt.date


#%% Conflicting IBS status 

'''

Some participants give different answers over time when asked if they have even been diagnosed with IBS.

For example:
    
    2009 - Yes
    2011 - No
    2013 - Yes
    2014 - No 
    etc.

Descriptives for various IBS statuses are provided before and after accounting for those with controversial IBS 
status. Since some of these counts will change after removing those with unascertainable IBS status or those with 
an assumed mistaken entry, the latter set of descriptives should be considered final.

'''
              
### Dataframe to summarise results & facilitate maniupulating the data

def ibs_stats(ph36, ibs_col):
    
    '''
    Makes a dataframe containing information about how many positive 
    reponses a participant has and whether they have 0s after 
    previously reports 1s in the response column (which would imply
    conflicting IBS status)
    '''
    
    number_instances, pc_ones, zero_after_1 = [], [], []
    
    for i in ph36.ParticipantID.unique():
        df = ph36[ph36.ParticipantID == i]
        number_instances.append(len(df))
        pc_ones.append(100 * df[ibs_col].sum() / len(df))
        
        zao = 'False'
        
        if df[ibs_col].sum() > 0:
            first_one = np.where(df[ibs_col] == 1)[0][0]   
            
            if np.any(df[ibs_col].values[first_one:] == 0):
                zao = 'True'
                
        zero_after_1.append(zao)
    
    col_names = ['iid', 'num_inst', 'percentage_ones', '0_after_1']
    ones_df = pd.DataFrame({
        col_names[0]: ph36.ParticipantID.unique(),
        col_names[1]: number_instances,
        col_names[2]: pc_ones,
        col_names[3]: zero_after_1
    })
    
    return ones_df

ones_df = ibs_stats(ph36, 'PH0000536')

## Get IDs that have a mixture of 0s and 1s in their responses
mixed = ones_df[(ones_df.percentage_ones > 0) & (ones_df.percentage_ones < 100)]
# len = 990

# Remove those with "0_after_1" == False because these are genuine (it means 
# they developed IBS during study)
mixed = mixed[mixed['0_after_1'] == 'True']
# len = 696


# Remove those for whom an entry may be a mistake (e.g., 1 or 0 only appears once in 4 instances)
def mistake_locs(mixed):
    locs = []
    for i in mixed.index:
        if mixed.loc[i]['num_inst'] >= 4:
            # Note: 0.05 added or subtracting to avoid errors due to rounding
            if (mixed.loc[i]['percentage_ones'] > (100*(1/mixed.loc[i]['num_inst']))*(mixed.loc[i]['num_inst']-1)-0.05) |  (mixed.loc[i]['percentage_ones'] < 0.5+(100*(1/mixed.loc[i]['num_inst']))):
                locs.append(i)
    return locs
    
                              
mistake = mixed.loc[mistake_locs(mixed)] # 287
controv = mixed.drop(mistake_locs(mixed)) # 409


#### Use age of diagnosis to resolve some conflicts 
## Calculate year of diagnosis based on age at which they report being diagnosed 
# Merge YOB with agel df
agel = pd.merge(agel, yob[['ParticipantID', 'YearOfBirth']], on = 'ParticipantID')
agel['YearDiagnosis'] = agel.YearOfBirth + agel.ageD

# Keep only controversial cases 
agel = agel[agel.ParticipantID.isin(controv.iid)] # Age of diagnosis available for 298 participants

# Now that we know year of diagnosis for these individuals, we can impute their 
# values for IBS diagnosis 
for i in agel.ParticipantID:
    date_cut = pd.to_datetime(agel.loc[agel.ParticipantID == i, 'YearDiagnosis'])
    mask_greater = ph36['ParticipantID'] == i
    ph36.loc[(mask_greater) & (ph36['ResponseDate'] > date_cut.iloc[0]), 'PH0000536'] = 1
    ph36.loc[(mask_greater) & (ph36['ResponseDate'] < date_cut.iloc[0]), 'PH0000536'] = 0

# Now redo key statistics:
ones_df_new = ibs_stats(ph36, 'PH0000536')

## Get IDs that have a mixture of 0s and 1s in their responses
mixed = ones_df_new[(ones_df_new.percentage_ones > 0) & (ones_df_new.percentage_ones < 100)]
# len = 692

# Remove "0_after_1" (i.e., developed IBS)
mixed = mixed[mixed['0_after_1'] == 'True']
# len = 398
 
mistake = mixed.loc[mistake_locs(mixed)] #287
controv = mixed.drop(mistake_locs(mixed)) #111
## IDs in controv will be dropped after accounting for Rome III criteria

## Impute those in mistake with 1
for i in mistake.iid:
    if int(mistake[mistake.iid == i]['percentage_ones'].values[0]) < 50:
        ph36.loc[ph36.ParticipantID == i, 'PH0000536'] = 0
    else:
        ph36.loc[ph36.ParticipantID == i, 'PH0000536'] = 1


#%% Rome 3

'''
Determine IBS status based on Rome 3 criteria 

Recurrent abdominal pain or discomfort (defined as an uncomfortable sensation 
not described as pain) for at least 3 days/month in the last 3 months, 
associated with two or more of the following:

    1. Improvement with defecation
    2. Onset associated with a change in the frequency of stool
    3. Onset associated with a change in the form (appearance) of stool



PH0001087	In the last 3 months, how often did you have discomfort or pain anywhere in your abdomen?
PH0001090	How often did this discomfort or pain get better or stop after you had a bowel movement?
PH0001091   When this discomfort or pain started, did you have more frequent bowel movements?
PH0001092	When this discomfort or pain started, did you have less frequent bowel movements?
PH0001093	When this discomfort or pain started, were your stools (bowel movements) looser?
PH0001094	When this discomfort or pain started, how often did you have harder stools?

'''

rc[rc.PhenotypeID == 'PH0001093'][['ResponseCode', 'ResponseDescription']].drop_duplicates()
rc.ResponseCode.value_counts(dropna = False)

# Change skipped legitimately and genuinely missing 
rc.loc[rc.ResponseCode == 999911, 'ResponseCode'] = 0
rc.loc[rc.ResponseCode == 999906, 'ResponseCode'] = np.nan

## PH0001087 - In the last 3 months, how often did you have discomfort or pain anywhere in your abdomen?
rc.loc[(rc.PhenotypeID == 'PH0001087') & (rc.ResponseCode < 4), 'ResponseCode'] = 0
rc.loc[(rc.PhenotypeID == 'PH0001087') & (rc.ResponseCode > 3), 'ResponseCode'] = 1

## PH0001090 - Improvement with defecation 
rc.loc[(rc.PhenotypeID == 'PH0001090') & (rc.ResponseCode < 2), 'ResponseCode'] = 0
rc.loc[(rc.PhenotypeID == 'PH0001090') & (rc.ResponseCode > 1), 'ResponseCode'] = 1

## Change in bowel movement frequency 
# More freq
rc.loc[(rc.PhenotypeID == 'PH0001091') & (rc.ResponseCode < 2), 'ResponseCode'] = 0
rc.loc[(rc.PhenotypeID == 'PH0001091') & (rc.ResponseCode > 1), 'ResponseCode'] = 1
# Less freq
rc.loc[(rc.PhenotypeID == 'PH0001092') & (rc.ResponseCode < 2), 'ResponseCode'] = 0
rc.loc[(rc.PhenotypeID == 'PH0001092') & (rc.ResponseCode > 1), 'ResponseCode'] = 1

## Change in bowel movements  consistency 
# Looser 
rc.loc[(rc.PhenotypeID == 'PH0001093') & (rc.ResponseCode < 2), 'ResponseCode'] = 0
rc.loc[(rc.PhenotypeID == 'PH0001093') & (rc.ResponseCode > 1), 'ResponseCode'] = 1
# Harder 
rc.loc[(rc.PhenotypeID == 'PH0001094') & (rc.ResponseCode < 2), 'ResponseCode'] = 0
rc.loc[(rc.PhenotypeID == 'PH0001094') & (rc.ResponseCode > 1), 'ResponseCode'] = 1

## Create column to determine if someone has IBS based on Rome 3

ibs_rc_df = rc.loc[(rc.PhenotypeID == 'PH0001087')][['ParticipantID', 'ResponseDate', 'ResponseCode']]
ibs_rc_df = ibs_rc_df.rename(columns = {'ResponseCode': 'PH0001087'})

ibs_rc_df = pd.merge(ibs_rc_df, rc.loc[(rc.PhenotypeID == 'PH0001090')].rename(columns={'ResponseCode':'PH0001090'})[['ParticipantID', 'PH0001090']], on = 'ParticipantID')
ibs_rc_df = pd.merge(ibs_rc_df, rc.loc[(rc.PhenotypeID == 'PH0001091')].rename(columns={'ResponseCode':'PH0001091'})[['ParticipantID', 'PH0001091']], on = 'ParticipantID')
ibs_rc_df = pd.merge(ibs_rc_df, rc.loc[(rc.PhenotypeID == 'PH0001092')].rename(columns={'ResponseCode':'PH0001092'})[['ParticipantID', 'PH0001092']], on = 'ParticipantID')
ibs_rc_df = pd.merge(ibs_rc_df, rc.loc[(rc.PhenotypeID == 'PH0001093')].rename(columns={'ResponseCode':'PH0001093'})[['ParticipantID', 'PH0001093']], on = 'ParticipantID')
ibs_rc_df = pd.merge(ibs_rc_df, rc.loc[(rc.PhenotypeID == 'PH0001094')].rename(columns={'ResponseCode':'PH0001094'})[['ParticipantID', 'PH0001094']], on = 'ParticipantID')

ibs_rc_df['change_freq'] = np.nan
# Change freq = 1 if either more freq, less freq, or both; else 0
ibs_rc_df.loc[(ibs_rc_df.PH0001091 == 1) | (ibs_rc_df.PH0001092 == 1), 'change_freq'] = 1
ibs_rc_df.loc[(ibs_rc_df.PH0001091 == 0) & (ibs_rc_df.PH0001092 == 0), 'change_freq'] = 0

ibs_rc_df['change_mvmnt'] = np.nan
# Change freq = 1 if either more freq, less freq, or both; else 0
ibs_rc_df.loc[(ibs_rc_df.PH0001093 == 1) | (ibs_rc_df.PH0001094 == 1), 'change_mvmnt'] = 1
ibs_rc_df.loc[(ibs_rc_df.PH0001093 == 0) & (ibs_rc_df.PH0001094 == 0), 'change_mvmnt'] = 0

## Final column for IBS by Rome 3
ibs_rc_df['ibs_rc'] = np.nan
ibs_rc_df.loc[ibs_rc_df.PH0001087 == 0, 'ibs_rc'] = 0

ibs_rc_df.loc[(ibs_rc_df['PH0001087'] == 1) & (ibs_rc_df[['PH0001090', 'change_freq', 'change_mvmnt']].sum(axis = 1) < 2), 'ibs_rc'] = 0
ibs_rc_df.loc[(ibs_rc_df['PH0001087'] == 1) & (ibs_rc_df[['PH0001090', 'change_freq', 'change_mvmnt']].sum(axis = 1) > 1), 'ibs_rc'] = 1

#### Establish subtypes 
ibs_rc_df['subtype'] = 'NA'

ibs_rc_df.loc[((ibs_rc_df.ibs_rc == 1) & ((ibs_rc_df.PH0001091 == 1) | (ibs_rc_df.PH0001093 == 1))), 'subtype'] = 'D'
ibs_rc_df.loc[((ibs_rc_df.ibs_rc == 1) & ((ibs_rc_df.PH0001092 == 1) | (ibs_rc_df.PH0001094 == 1))), 'subtype'] = 'C'
ibs_rc_df.loc[((ibs_rc_df.ibs_rc == 1) & ((((ibs_rc_df.PH0001091 == 1) | (ibs_rc_df.PH0001093 == 1))) & (((ibs_rc_df.PH0001092 == 1) | (ibs_rc_df.PH0001094 == 1))))), 'subtype'] = 'M'

### Add Rome III data from Q18 

'''

Q18_115	0	Never
Q18_115	1	Less than one day a month
Q18_115	2	One day a month
Q18_115	3	Two to three days a month
Q18_115	4	One day a week
Q18_115	5	More than one day a week
Q18_115	6	Every day


0	Never/rarely
1	Sometimes
2	Often
3	Most of the time
4	Always


0	Never/rarely
1	About 25% of the time
2	About 50% of the time
3	About 75% of the time
4	Always, 100% of the time

Anything above 0 is considered positive

Approach: positive response to any question relating to constipation (120, 122, 124)
or diarrhoea (119, 121, 123) is coded as having that subtype. 
'''

# Drop those with NA for abdominal pain
Q18 = Q18[~Q18.Q18_115.isnull()]

# abdominal pain < 3 days per month cannot have IBS
Q18_no_ibs = Q18[Q18.Q18_115 < 4]

Q18_pain = Q18[Q18.Q18_115 > 3]

Q18_pain.loc[Q18_pain.Q18_118 < 2, 'Q18_118'] = 0
Q18_pain.loc[Q18_pain.Q18_118 > 1, 'Q18_118'] = 1

Q18_pain.loc[Q18_pain.Q18_119 < 2, 'Q18_119'] = 0
Q18_pain.loc[Q18_pain.Q18_119 > 1, 'Q18_119'] = 1

Q18_pain.loc[Q18_pain.Q18_120 < 2, 'Q18_120'] = 0
Q18_pain.loc[Q18_pain.Q18_120 > 1, 'Q18_120'] = 1

Q18_pain.loc[Q18_pain.Q18_121 < 2, 'Q18_121'] = 0
Q18_pain.loc[Q18_pain.Q18_121 > 1, 'Q18_121'] = 1

Q18_pain.loc[Q18_pain.Q18_122 < 2, 'Q18_122'] = 0
Q18_pain.loc[Q18_pain.Q18_122 > 1, 'Q18_122'] = 1


Q18_pain['change_freq'] = 0
Q18_pain['change_mvmnt'] = 0

Q18_pain.loc[(Q18_pain.Q18_119 == 1) | (Q18_pain.Q18_120 == 1), 'change_freq'] = 1
Q18_pain.loc[(Q18_pain.Q18_121 == 1) | (Q18_pain.Q18_122 == 1), 'change_mvmnt'] = 1

Q18_pain['meets_rome'] = 0 
Q18_pain.loc[Q18_pain[['Q18_118', 'change_freq', 'change_mvmnt']].sum(axis = 1) > 1, 'meets_rome'] = 1


#### Establish subtypes 
Q18_pain['subtype'] = 'NA'

Q18_pain.loc[((Q18_pain.meets_rome == 1) & ((Q18_pain.Q18_119 == 1) | (Q18_pain.Q18_121 == 1))), 'subtype'] = 'D'
Q18_pain.loc[((Q18_pain.meets_rome == 1) & ((Q18_pain.Q18_120 == 1) | (Q18_pain.Q18_122 == 1))), 'subtype'] = 'C'
Q18_pain.loc[((Q18_pain.meets_rome == 1) & ((((Q18_pain.Q18_119 == 1) | (Q18_pain.Q18_121 == 1))) & (((Q18_pain.Q18_120 == 1) | (Q18_pain.Q18_122 == 1))))), 'subtype'] = 'M'

## Merge those with subtypes with those without IBS
Q18_subtypes = Q18_pain[['ParticipantID', 'ResponseDate', 'subtype']].rename(columns = {'ResponseYear':'ResponseDate'})
Q18_no_ibs_proc = Q18_no_ibs[['ParticipantID', 'ResponseDate']]
Q18_no_ibs_proc['subtype'] = 'NA'
Q18_merged = pd.concat([Q18_subtypes, Q18_no_ibs_proc], axis = 0)

## Merge two dataframes with subtype info
Q18_merged['Origin'] = 'Rome_III_Q18'
ibs_rc_df['Origin'] = 'Rome_III_PH_codes'
ibs_rc_df = ibs_rc_df[['ParticipantID', 'ResponseDate', 'subtype', 'Origin']]

subtypes = pd.concat([Q18_merged, ibs_rc_df], axis = 0)
subtypes = subtypes.sort_values(['ParticipantID', 'ResponseDate'])
subtypes = subtypes.reset_index(drop = True).rename(columns = {'ParticipantID':'iid'})

### Reclassify those who have multiple subtypes as having mixed 
# First, keep only those with multiple readings over time
to_reclassify = subtypes[subtypes.iid.isin(subtypes.iid.value_counts()[subtypes.iid.value_counts() > 1].index)]

# Then, remove those who do not have IBS at an earlier time point but then develop 
# it or those without IBS at both time points
to_drop = to_reclassify.groupby('iid').filter(lambda x: x['subtype'].values[0] == 'NA')['iid'].unique()
to_reclassify = to_reclassify[~to_reclassify.iid.isin(to_drop)]

# Next, if second observation is NA, convert to be same as first observation
for i in to_reclassify[to_reclassify.subtype == 'NA']['iid']:
    to_reclassify.loc[(to_reclassify.iid == i) & (to_reclassify.subtype == 'NA'), 'subtype'] = to_reclassify.loc[(to_reclassify.iid == i) & (to_reclassify.subtype != 'NA'), 'subtype'].values[0]

# Now, get any that have different subtypes at each date and relabel as mixed 
reclassify = to_reclassify.groupby('iid')['subtype'].apply(lambda x: len(x.unique()) > 1).index[to_reclassify.groupby('iid')['subtype'].apply(lambda x: len(x.unique()) > 1)]
to_reclassify.loc[to_reclassify.iid.isin(reclassify), 'subtype'] = 'M'

# Now recode in original subtypes df
for i in to_reclassify[['iid', 'subtype']].drop_duplicates()['iid']:
    subtypes.loc[subtypes.iid == i, 'subtype'] = to_reclassify[to_reclassify.iid == i]['subtype'].values[0]
    
#subtypes.to_csv(r"C:\Users\k2365231\OneDrive - King's College London\Year 2\metabolites_ibs\metabolite_codes\ibs_metabolite_datasets\processed\rome_subtypes.csv")

#%%  Return to overall IBS processing after processing subtypes 

meets_rome = Q18_pain[Q18_pain.meets_rome == 1]['ParticipantID']
Q18_pain = Q18_pain.dropna()

Q18_no_ibs = np.append(Q18_no_ibs.ParticipantID, Q18_pain[Q18_pain.meets_rome == 0]['ParticipantID'].values)

Q18_rc = pd.DataFrame({'ParticipantID':np.append(Q18_no_ibs, meets_rome), 
                       'ibs_status':np.append(np.repeat(0, len(Q18_no_ibs)), np.repeat(1, len(meets_rome))), 
                       'Origin':'Rome_III_Q18'})

Q18_rc = pd.merge(Q18_rc, Q18[['ResponseDate', 'ParticipantID']], on = 'ParticipantID')

## Convert subtype to IBS status in ibs_rc_df
ibs_rc_df.loc[ibs_rc_df.subtype != 'NA', 'subtype'] = 1
ibs_rc_df.loc[ibs_rc_df.subtype == 'NA', 'subtype'] = 0
ibs_rc_df, ph36 = ibs_rc_df.rename(columns = {'subtype':'ibs_status'}), ph36.rename(columns = {'PH0000536':'ibs_status'})

## Merge Q18_rc with ibs_rc_df
ibs_rc_df = pd.concat([ibs_rc_df[['ParticipantID','ResponseDate', 'ibs_status', 'Origin']], Q18_rc], axis = 0).reset_index(drop=True)
ph36_rc = pd.concat([ph36, ibs_rc_df], axis = 0).sort_values(['ParticipantID', 'ResponseDate', 'Origin']).reset_index(drop = True)
ph36_rc.ResponseDate = pd.to_datetime(ph36_rc.ResponseDate)


#### Append this information to IBS 
# First, drop controversial cases who do not have positive Rome III hits 
to_keep = ph36_rc[ph36_rc.ParticipantID.isin(controv.iid) & ((ph36_rc.Origin == 'Rome_III_PH_codes') | (ph36_rc.Origin == 'Rome_III_Q18')) & (ph36_rc.ibs_status == 1)]['ParticipantID']
controv = controv[~controv.iid.isin(to_keep)]
ph36_rc = ph36_rc[~ph36_rc.ParticipantID.isin(controv.iid)]

### Account for those who state not having IBS but qualify based on RomeIII 
ph36_rc['ibs_status_upd'] = ph36_rc['ibs_status']

'''
If someone has been diagnosed with IBS, Rome 3 criteria is ignored. 

If someone reports not being diagnosed but meets Rome 3 at a particular time point, 
from that point onwards they are coded as 1

This section of the code will use ibs_status_upd to recode those who report IBS 
at an earlier time point as having it at all future time points.

This will be later used to determine true controls from cases
'''

ph36_rc_upd = ph36_rc.reset_index(drop = True)

ph36_without_controv = ph36[~ph36.ParticipantID.isin(controv.iid)]

all_ones = ph36_without_controv.groupby('ParticipantID')['ibs_status'].apply(lambda x: np.all(x == 1)).index[ph36_without_controv.groupby('ParticipantID')['ibs_status'].apply(lambda x: np.all(x == 1))].to_list()
all_zeros = ph36_without_controv.groupby('ParticipantID')['ibs_status'].apply(lambda x: np.all(x == 0)).index[ph36_without_controv.groupby('ParticipantID')['ibs_status'].apply(lambda x: np.all(x == 0))].to_list()

## Convert those with all ones for diagnosis regardless of Rome_III
ph36_rc_upd.loc[ph36_rc_upd.ParticipantID.isin(all_ones), 'ibs_status_upd'] = 1

to_change = ph36_rc_upd.ParticipantID[~ph36_rc_upd.ParticipantID.isin(all_ones)].unique()

for i in to_change:
    if np.any(ph36_rc_upd[ph36_rc_upd.ParticipantID == i]['ibs_status'] > 0):
        # Get first date they report having IBS
        date_first_one = ph36_rc_upd.loc[(ph36_rc_upd.ibs_status == 1) & (ph36_rc_upd.ParticipantID == i), 'ResponseDate'].values[0]
        # Update all dates following this 
        ph36_rc_upd.loc[(ph36_rc_upd.ParticipantID == i) & (ph36_rc_upd.ResponseDate >= date_first_one), 'ibs_status_upd'] = 1


def check_zero_after_one(group):
    zero_after_one = []
    for idx, row in group.iterrows():
        ibs_status = row['ibs_status_upd']
        if ibs_status == 1:
            # Check if there's a 0 after the current row
            if 0 in group.loc[idx:, 'ibs_status_upd'].values:
                zero_after_one.append(1)
            else:
                zero_after_one.append(0)
        else:
            zero_after_one.append(0)
    return pd.Series(zero_after_one, index=group.index)

# Apply the function to each group of participantIDs
result = ph36_rc_upd.groupby('ParticipantID').apply(check_zero_after_one)
(result > 0).sum() # Should be 0


#%% Prev/inc file 

'''
Adds columns denoting prevalence, incidence and age of diagnosis:
    
    Prevalence - Whether an individual has IBS (according to ibs_status_upd) at a given time point
    Incidence - Whether an individual develops IBS at a later time point
    Age - Adds date on which they first have IBS (in case in incidence) or most recent date (in case of controls)
'''

prev_inc_age = ph36_rc_upd

all_ones = ph36_rc_upd.groupby('ParticipantID')['ibs_status_upd'].apply(lambda x: np.all(x == 1)).index[ph36_rc_upd.groupby('ParticipantID')['ibs_status_upd'].apply(lambda x: np.all(x == 1))].to_list()
all_zeros = ph36_rc_upd.groupby('ParticipantID')['ibs_status_upd'].apply(lambda x: np.all(x == 0)).index[ph36_rc_upd.groupby('ParticipantID')['ibs_status_upd'].apply(lambda x: np.all(x == 0))].to_list()

dev_ibs = ph36_rc_upd[~ph36_rc_upd.ParticipantID.isin(all_ones+all_zeros)]

def first_instance_of_one(df, participant_col, date_col, value_col):

    first_instance_dict = {}

    # Group by participant
    grouped = df.groupby(participant_col)

    # Iterate over each group
    for participant, group in grouped:
        # Find the index of the first occurrence of 1
        first_instance_index = group[value_col].eq(1).idxmax()
        # Get the date corresponding to the first instance
        first_instance_date = group.loc[first_instance_index, date_col]
        # Store participant ID and date in the dictionary
        first_instance_dict[participant] = first_instance_date

    return first_instance_dict

first_ones = first_instance_of_one(dev_ibs, 'ParticipantID', 'ResponseDate', 'ibs_status_upd')

prev_inc_age['prevalent_ibs'] = 0
prev_inc_age.loc[prev_inc_age.ParticipantID.isin(all_ones), 'prevalent_ibs'] = 1

prev_inc_age['incident_ibs'] = 0
prev_inc_age.loc[prev_inc_age.ParticipantID.isin(all_ones), 'incident_ibs'] = 'N/A'


for participant, date in first_ones.items():
    # Filter rows for the participant and date
    mask = (prev_inc_age.ParticipantID == participant) & (prev_inc_age.ResponseDate >= date)
    # Update values where conditions are met
    prev_inc_age.loc[mask, 'prevalent_ibs'] = 1
    
    prev_inc_age.loc[prev_inc_age.ParticipantID == participant, 'incident_ibs'] = 1
    

## Add age

prev_inc_age['age_incident_ibs'] = 'N/A'

most_recent_dates = ph36_rc_upd[ph36_rc_upd.ParticipantID.isin(all_zeros)].groupby('ParticipantID')['ResponseDate'].max()

## Add most recent date for those who never report IBS
for i in most_recent_dates.index:
    prev_inc_age.loc[prev_inc_age.ParticipantID == i, 'age_incident_ibs'] = most_recent_dates[i]

## Add date of incidence for those who develop IBS during study 
for i in first_ones.keys():
    prev_inc_age.loc[prev_inc_age.ParticipantID == i, 'age_incident_ibs'] = first_ones[i]

## Checks
prev_inc_age[(prev_inc_age.incident_ibs == 1) & (prev_inc_age.prevalent_ibs == 0)][['ParticipantID', 'ResponseDate', 'prevalent_ibs', 'incident_ibs', 'age_incident_ibs']].reset_index(drop = True)[:50]


# Remove time string from datetime
prev_inc_age['age_incident_ibs'] = pd.to_datetime(prev_inc_age['age_incident_ibs'], errors='coerce')
prev_inc_age['age_incident_ibs'] = prev_inc_age['age_incident_ibs'].dt.date.fillna('N/A')

#prev_inc_age.to_csv(r"IBS_dataset.csv")

