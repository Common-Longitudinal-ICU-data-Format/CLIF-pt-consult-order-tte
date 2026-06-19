#!/usr/bin/env python
# coding: utf-8

# # MIMIC-CLIF Early PT Consults Comparative Analysis
# ## Step 3: Used the Data Gathered for some Calculations
# 
# - Applies concensus criteria using the hourly dataframe.
# - Calculates outcomes
# - Converts date_time values to hours or days as needed.
# - Clustering of categorical values, assuming appropriate CLIF definitions.
# - Creates a summary TABLE 1.
# - Creates a graph.

# ## Setup

# In[1]:


### Import
#Import packages, config file and load CLIF orchestrator.
import pandas as pd
import pyarrow
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
# Force white backgrounds regardless of marimo/system dark theme
matplotlib.rcParams.update({
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'savefig.facecolor': 'white',
    'axes.edgecolor': 'black',
    'axes.labelcolor': 'black',
    'xtick.color': 'black',
    'ytick.color': 'black',
    'text.color': 'black',
})
import os
import sys
import shutil
from datetime import datetime, timedelta
import json
import warnings
warnings.filterwarnings('ignore')
import clifpy
import logging

#Our own helper function
import pthelperfunctions as helper

#file paths
work_dir = os.path.abspath('..')
output_folder = os.path.join(work_dir,'output')

#Config
with open(os.path.join(work_dir,'config','config.json'), 'r') as file:
    config = json.load(file)

#Load block data CLIF-Eligibility-for-mobilization output
block_df = helper.load_data('output_folder','block_df_2_end',folder='intermediate')

#Load Time Bin Object
time_bin = helper.time_bins(in_name='time_bin_step_2')

#Load Hourly DF Object
hourly = helper.hourly_blocks(in_name='hourly_df_two')


# In[2]:


_logger = logging.getLogger('clif_01')
_logger.setLevel(logging.INFO)
_logger.handlers.clear()

_log_dir = os.path.join(output_folder,'logs',f'{config['site_name']}_03_calculations_log.txt')
_fh = logging.FileHandler(_log_dir, mode='w')
_fh.setFormatter(logging.Formatter('%(asctime)s | %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
_logger.addHandler(_fh)

_ch = logging.StreamHandler()
_ch.setFormatter(logging.Formatter('%(message)s'))
_logger.addHandler(_ch)

def log(*args, **kwargs):
    _msg = ' '.join(str(a) for a in args)
    _logger.info(_msg)

log(f"=== CLIF Pipeline 03: Calculations ===")
log(f"Site: {config['site_name']}")


# ## Hourly Mobilization Analysis
# 
# Uses hourly data frame built earlier along with algorithm by Kaveri. Original code seen here: [CLIF-eligibility-for-mobilization](https://github.com/Common-Longitudinal-ICU-data-Format/CLIF-eligibility-for-mobilization/blob/main/code/02_mobilization_analysis.py)

# In[3]:


def compute_consensus_flags(df):
    # Derive helper columns
    df['recorded_hour'] = df['window_start_dttm'].dt.hour
    df['is_weekday'] = df['window_start_dttm'].dt.weekday < 5
    
    # --- RED flags ---
    df['red_resp_spo2_flag'] = ((df['spo2_min'] < 90) | df['spo2_min'].isna()).astype(int)
    df['red_map_flag'] = ((df['map_mean'] < 65) | df['map_mean'].isna()).astype(int)
    df['red_high_support_flag'] = ((df['ne_calc_last'] > 0.3) | (df['ne_calc_max'] > 0.3)).astype(int)
    df['red_hypertensive_flag'] = (
        (((df['sbp_max'] > 200) | (df['map_mean'] > 110)) &
        (df['red_med_flag'] == 1))
    ).astype(int)
    df['red_pulse_high_flag'] = (df['heart_rate_max'] > 150).astype(int)
    df['red_pulse_low_flag'] = ((df['heart_rate_min'] < 40) | df['heart_rate_min'].isna()).astype(int)

    # --- YELLOW flags ---
    df['yellow_resp_spo2_flag'] = ((df['spo2_min'] >= 90) | df['spo2_min'].isna()).astype(int)
    df['yellow_fio2_flag'] = (df['fio2_set_min'] > 0.6).astype(int)
    df['yellow_resp_rate_flag'] = (df['respiratory_rate_max'] > 30).astype(int)
    df['yellow_peep_flag'] = (df['peep_set_min'] > 10).astype(int)
    df['yellow_map_flag'] = ((df['map_mean'] >= 65) & (df['ne_calc_last'].between(0.1, 0.3))).astype(int)
    df['yellow_pulse_flag'] = (df['heart_rate_min'].between(120, 150)).astype(int)
    df['yellow_lactate_flag'] = (df['lactate_max'] > 4).astype(int)

    # --- GREEN flags ---
    df['green_resp_spo2_flag'] = ((df['spo2_min'] >= 90) | df['spo2_min'].isna()).astype(int)
    df['green_resp_rate_flag'] = ((df['respiratory_rate_max'] <= 30) | df['respiratory_rate_max'].isna()).astype(int)
    df['green_fio2_flag'] = ((df['fio2_set_min'] <= 0.6) | df['fio2_set_min'].isna()).astype(int)
    df['green_peep_flag'] = ((df['peep_set_min'] <= 10) | df['peep_set_min'].isna()).astype(int)
    df['green_map_flag'] = (((df['map_mean'] >= 65) & (df['ne_calc_last'] < 0.1)) | df['ne_calc_last'].isna()).astype(int)
    df['green_pulse_flag'] = ((df['heart_rate_min'] < 120) | df['heart_rate_min'].isna()).astype(int)
    df['green_lactate_flag'] = ((df['lactate_max'] <= 4) | df['lactate_max'].isna()).astype(int)
    df['green_hr_flag'] = ((df['heart_rate_min'] > 40) | df['heart_rate_min'].isna()).astype(int)

    # --- Composite flags (shared conditions) ---
    _base = (
        (df['tracheostomy_max'] == 0) & (df['paralytics_flag'] == 0) &
        (df['time_from_vent'] > 4)
    )
    _daytime = _base & (df['recorded_hour'] >= 8) & (df['recorded_hour'] < 17)
    _weekday = _daytime & (df['is_weekday'] == True)

    df['any_red'] = (
        (df['red_resp_spo2_flag'] | df['red_map_flag'] | df['red_high_support_flag'] |
         df['red_hypertensive_flag'] | df['red_pulse_high_flag'] | df['red_pulse_low_flag']) &
        _base
    ).astype(int)

    df['no_red'] = (
        ~(df['red_resp_spo2_flag'] | df['red_map_flag'] | df['red_high_support_flag'] |
          df['red_hypertensive_flag'] | df['red_pulse_high_flag'] | df['red_pulse_low_flag']) &
        _daytime
    ).astype(int)

    df['any_yellow'] = (
        (df['yellow_resp_spo2_flag'] | df['yellow_fio2_flag'] | df['yellow_resp_rate_flag'] |
         df['yellow_peep_flag'] | df['yellow_map_flag'] | df['yellow_pulse_flag'] |
         df['yellow_lactate_flag']) &
        _daytime
    ).astype(int)

    df['any_green'] = (
        (df['green_resp_spo2_flag'] | df['green_resp_rate_flag'] | df['green_fio2_flag'] |
         df['green_peep_flag'] | df['green_map_flag'] | df['green_pulse_flag'] |
         df['green_lactate_flag'] | df['green_hr_flag']) &
        _daytime
    ).astype(int)

    df['all_green'] = (
        df['green_resp_spo2_flag'] & df['green_resp_rate_flag'] & df['green_fio2_flag'] &
        df['green_peep_flag'] & df['green_map_flag'] & df['green_pulse_flag'] &
        df['green_lactate_flag'] & df['green_hr_flag'] & _daytime
    ).astype(int)

    df['all_green_all_hours'] = (
        df['green_resp_spo2_flag'] & df['green_resp_rate_flag'] & df['green_fio2_flag'] &
        df['green_peep_flag'] & df['green_map_flag'] & df['green_pulse_flag'] &
        df['green_lactate_flag'] & df['green_hr_flag'] & _base
    ).astype(int)

    df['all_green_weekday'] = (
        df['green_resp_spo2_flag'] & df['green_resp_rate_flag'] & df['green_fio2_flag'] &
        df['green_peep_flag'] & df['green_map_flag'] & df['green_pulse_flag'] &
        df['green_lactate_flag'] & df['green_hr_flag'] & _weekday
    ).astype(int)

    df['all_green_no_red'] = (
        df['green_resp_spo2_flag'] & df['green_resp_rate_flag'] & df['green_fio2_flag'] &
        df['green_peep_flag'] & df['green_map_flag'] & df['green_pulse_flag'] &
        df['green_lactate_flag'] & df['green_hr_flag'] & (df['any_red'] == 0) & _daytime
    ).astype(int)

    df['all_green_no_red_all_hours'] = (
        df['green_resp_spo2_flag'] & df['green_resp_rate_flag'] & df['green_fio2_flag'] &
        df['green_peep_flag'] & df['green_map_flag'] & df['green_pulse_flag'] &
        df['green_lactate_flag'] & df['green_hr_flag'] & (df['any_red'] == 0) & _base
    ).astype(int)

    df['all_green_no_red_weekday'] = (
        df['green_resp_spo2_flag'] & df['green_resp_rate_flag'] & df['green_fio2_flag'] &
        df['green_peep_flag'] & df['green_map_flag'] & df['green_pulse_flag'] &
        df['green_lactate_flag'] & df['green_hr_flag'] & (df['any_red'] == 0) & _weekday
    ).astype(int)

    df['all_green_no_red_yellow'] = (
        df['green_resp_spo2_flag'] & df['green_resp_rate_flag'] & df['green_fio2_flag'] &
        df['green_peep_flag'] & df['green_map_flag'] & df['green_pulse_flag'] &
        df['green_lactate_flag'] & df['green_hr_flag'] &
        (df['any_red'] == 0) & (df['any_yellow'] == 0) & _daytime
    ).astype(int)

    df['any_yellow_or_green_no_red'] = (
        (df['yellow_resp_spo2_flag'] | df['yellow_fio2_flag'] | df['yellow_resp_rate_flag'] |
         df['yellow_peep_flag'] | df['yellow_map_flag'] | df['yellow_pulse_flag'] |
         df['yellow_lactate_flag'] | df['green_resp_spo2_flag'] | df['green_resp_rate_flag'] |
         df['green_fio2_flag'] | df['green_peep_flag'] | df['green_map_flag'] |
         df['green_pulse_flag'] | df['green_lactate_flag'] | df['green_hr_flag']) &
        (df['any_red'] == 0) & _daytime
    ).astype(int)

    df['any_yellow_or_green_no_red_weekday'] = (
        (df['yellow_resp_spo2_flag'] | df['yellow_fio2_flag'] | df['yellow_resp_rate_flag'] |
         df['yellow_peep_flag'] | df['yellow_map_flag'] | df['yellow_pulse_flag'] |
         df['yellow_lactate_flag'] | df['green_resp_spo2_flag'] | df['green_resp_rate_flag'] |
         df['green_fio2_flag'] | df['green_peep_flag'] | df['green_map_flag'] |
         df['green_pulse_flag'] | df['green_lactate_flag'] | df['green_hr_flag']) &
        (df['any_red'] == 0) & _weekday
    ).astype(int)

    df['any_yellow_or_green_no_red_all_hours'] = (
        (df['yellow_resp_spo2_flag'] | df['yellow_fio2_flag'] | df['yellow_resp_rate_flag'] |
         df['yellow_peep_flag'] | df['yellow_map_flag'] | df['yellow_pulse_flag'] |
         df['yellow_lactate_flag'] | df['green_resp_spo2_flag'] | df['green_resp_rate_flag'] |
         df['green_fio2_flag'] | df['green_peep_flag'] | df['green_map_flag'] |
         df['green_pulse_flag'] | df['green_lactate_flag'] | df['green_hr_flag']) &
        (df['any_red'] == 0) & _base
    ).astype(int)

    df['green_resp_flag'] = (
        df['green_resp_spo2_flag'] & df['green_resp_rate_flag'] &
        df['green_fio2_flag'] & df['green_peep_flag'] & _daytime
    ).astype(int)

    df['green_cardio_flag'] = (
        df['green_map_flag'] & df['green_pulse_flag'] &
        df['green_lactate_flag'] & df['green_hr_flag'] & _daytime
    ).astype(int)

    df['yellow_resp_flag'] = (
        (df['yellow_resp_spo2_flag'] | df['yellow_fio2_flag'] | df['yellow_resp_rate_flag'] |
         df['yellow_peep_flag'] | df['green_resp_spo2_flag'] | df['green_resp_rate_flag'] |
         df['green_fio2_flag'] | df['green_peep_flag']) &
        (df['any_red'] == 0) & _daytime
    ).astype(int)

    df['yellow_cardio_flag'] = (
        (df['yellow_map_flag'] | df['yellow_pulse_flag'] | df['yellow_lactate_flag'] |
         df['green_map_flag'] | df['green_pulse_flag'] | df['green_lactate_flag'] | df['green_hr_flag']) &
        (df['any_red'] == 0) & _daytime
    ).astype(int)

    df['yellow_all_green'] = (df['all_green_no_red'] & (df['any_yellow'] == 0)).astype(int)
    df['yellow_not_all_green'] = (df['any_yellow_or_green_no_red'] & (df['all_green_no_red'] == 0)).astype(int)

    return df

hourly.df = compute_consensus_flags(hourly.df)


# In[4]:


hourly.save(suffix='_w_mob')
log('Mobilization calculations completed and summary saved.')
hourly.df['time_diff'] = hourly.df['time_from_vent']
hourly.df['time_bin'] = time_bin.classify_time_bin(hourly.df['time_diff'])


# ## Time to mobilization
# Use mobilization data to get a few variables.

# In[5]:


yellow_df = hourly.df.rename(columns={'any_yellow_or_green_no_red_all_hours':'yellow'}).copy()
yellow_df = yellow_df[['encounter_block','time_diff','time_bin','yellow']]

#First Time to eligibility
#Group and get the first hour
_mask = yellow_df['yellow'] == 1
grouped_yellow_df = (
    yellow_df[_mask] #Filter to only positive values
    .groupby('encounter_block')['time_diff']
    .min()
    .reset_index()
)
grouped_yellow_df.rename(columns={'time_diff': 'yellow_time_eligibility'}, inplace=True)
block_df = pd.merge(
    block_df,
    grouped_yellow_df[['encounter_block','yellow_time_eligibility']],
    on='encounter_block',
    how='left'
)
block_df['yellow_0_72h'] = (block_df['yellow_time_eligibility'] <= 72).astype(bool)
log('Calculated time to Yellow Readiness for Mobilization for at least 1 hour.')

#Total hours in first 24 hours
block_df = block_df.merge(
    helper.aggregate_by_time(
        yellow_df,
        'yellow',
        agg_func='sum'),
    on='encounter_block',
    how='left')
log('Calculated hours of Yellow Readiness for Mobilization in first 24-hours.')

#Eligibility for all 4 hours of each time_bin
time_bin.gather_time_bins(yellow_df[['encounter_block','time_bin','yellow']], 'yellow', agg_func='all')


#ELIGIBILITY FIRST 2-HOURS CONSECUTIVE
#Two consecutive hours (note this relies on the hourly_df being sorted which we do above.
yellow_df['yellow_2'] = (yellow_df['yellow'].shift(periods=1, fill_value=0) == 1) & yellow_df['yellow']
yellow_df = yellow_df[yellow_df['yellow_2'] & (yellow_df['time_diff'] > 1)] #The >1 is to prevent accidental counting of the previous encounter block.
#Group and get the first hour
_mask = yellow_df['yellow_2'] == 1
grouped_yellow_df = (
    yellow_df[_mask]
    .groupby('encounter_block')['time_diff']
    .min()
    .reset_index()
)
grouped_yellow_df.rename(columns={'time_diff': 'yellow_time_eligibility_2h'}, inplace=True)

block_df = pd.merge(
    block_df,
    grouped_yellow_df[['encounter_block','yellow_time_eligibility_2h']],
    on='encounter_block',
    how='left'
)
block_df['yellow_2h_0_72h'] = (block_df['yellow_time_eligibility_2h'] <= 72).astype(bool)
log('Calculated time to Yellow Readiness for Mobilization for at least 2 hour.')

del grouped_yellow_df, yellow_df


# ## Oversedation
# Based on 'coma' which was defined by RASS < -2 in the second notebook.

# In[6]:


coma_df = hourly.df[['encounter_block','time_diff','coma']].copy()

#Sum of hours in a coma in the first 24 hours
block_df = block_df.merge(
    helper.aggregate_by_time(
        coma_df,
        'coma',
        agg_func='sum'),
    on='encounter_block',
    how='left')
del coma_df
log('Calculated hours of oversedation.')


# ## Pressor Data

# In[7]:


#Pressor indicator
pressor_df = hourly.df[['encounter_block','time_diff','time_bin','ne_calc_max']].copy()
pressor_df['pressor'] = pressor_df['ne_calc_max'] > 0

#For 24 hour block data
block_df = block_df.merge(
    helper.aggregate_by_time(
        pressor_df[['encounter_block','time_diff','pressor']],
        'pressor',
        agg_func='flag'),
    on='encounter_block',
    how='left')
log('Calculated pressor use flag in the first 24-hours.')

#For time bins
time_bin.gather_time_bins(pressor_df[['encounter_block','time_bin','pressor']], 'pressor', agg_func='flag', fill_with=0)
del pressor_df
log('Calculated pressor use flag for time_bins.')


# ## Paralytics Data

# In[8]:


#Paralytics indicator
para_df = hourly.df[['encounter_block','time_diff','time_bin','paralytics_flag']].copy()
para_df.rename(columns={'paralytics_flag':'paralytics'}, inplace=True)

#For 24 hour block data
block_df = block_df.merge(
    helper.aggregate_by_time(
        para_df[['encounter_block','time_diff','paralytics']],
        'paralytics',
        agg_func='sum'),
    on='encounter_block',
    how='left')

#Convert to boolean
block_df['paralytics_0_24h_>3h'] = block_df['paralytics_0_24h_sum'] > 3
log('Calculated paralytics use flag for >4 horus in first 24-hours.')

#For time bins
time_bin.gather_time_bins(para_df, 'paralytics', agg_func='flag')
log('Calculated paralytics use flag for any amount of time in time_bins.')
del para_df


# ## Ventilator Data

# In[9]:


###VENT FREE DAYS
#If dead, 0
#Otherwise uses the last hour of IMV on the hourly data_frame
vent_df = hourly.df[['encounter_block','time_from_vent','time_diff','time_bin','hourly_on_vent']]
vent_df['hourly_on_vent'] = vent_df['hourly_on_vent'].astype(bool)
#Keep only values within 28-days and on-vent
vent_df = vent_df[(vent_df['time_from_vent'] <= 28*24) & vent_df['hourly_on_vent'] ]
#Get the MAX hour and merge it into DF
last_vent_df = (
    vent_df.groupby('encounter_block')['time_from_vent']
    .max()
    .reset_index()
)
last_vent_df['time_from_vent'] = last_vent_df['time_from_vent'].astype(int)
last_vent_df.rename(columns={'time_from_vent':'last_hour_on_vent'}, inplace=True)
block_df = pd.merge(
    block_df,
    last_vent_df,
    on='encounter_block',
    how='left'
)

#Get an 1 for patients alive at 28-days.
block_df['alive28'] = block_df['death_dttm'].isna() | ((block_df['death_dttm'] - block_df['block_vent_start_dttm']).dt.total_seconds() >= (28*24*60*60))
block_df['alive28'] = block_df['alive28'].astype(int)

#Calcute VFD
block_df['vent_free_days'] = block_df['alive28'] * (28 - block_df['last_hour_on_vent']/24)
block_df = block_df.drop(columns=['alive28','last_hour_on_vent'])

#REINTUBATIONS
vent_df = hourly.df[['encounter_block','time_from_vent','hourly_on_vent']] #Note it is already sorted in the proper order

#Intubation count for hospitalization only, not including other hospitalizations.
def count_intubations(series):
    """Counts the number of re-intubations as 0->1 transitions."""
    return ((series.shift(periods=1, fill_value=0) == 1) & (series == 0)).sum()

intubation_count_df = (
    vent_df
    .groupby('encounter_block')['hourly_on_vent']
    .apply(count_intubations)
    .reset_index()
    .rename(columns={'hourly_on_vent': 'intubation_count'})
)
block_df = pd.merge(
    block_df,
    intubation_count_df,
    on='encounter_block',
    how='left'
)
block_df['reintubation'] = (block_df['intubation_count'] > 1)

#VENT FLAG for BINS
hourly.df['vent'] = hourly.df['hourly_on_vent']
time_bin.gather_time_bins(hourly.df[['encounter_block','time_bin','vent']], 'vent', agg_func='flag')


# In[10]:


del intubation_count_df
del last_vent_df
del vent_df


# In[11]:


#SAVING POINT
path = os.path.join(output_folder, 'intermediate',"block_df_3_calculations.parquet")
block_df.to_parquet(path)
del path


# ## Close Time Bins Data Set

# In[12]:


#Censor out dead data
time_bin.remove_based_on_censor('death', keep_first=True)
#Save (which will save the data as well as a summary of it)
time_bin.save(suffix='_3_end')
#Save an additional version as a CSV for R.
path = os.path.join(output_folder, 'intermediate',"time_bins_3_end.csv")
time_bin.df.to_csv(path)
del path


# ## Date Time Calculations

# In[13]:


#Change relevant DTTM values to hours/days
block_df['imv_to_discharge_days'] = (block_df['discharge_dttm'] - block_df['block_vent_start_dttm']).dt.total_seconds()/(24*3600)
block_df['imv_to_end_hours'] = (block_df['block_vent_end_dttm'] - block_df['block_vent_start_dttm']).dt.total_seconds()/(3600)
block_df['adm_to_imv_hours'] = (block_df['block_vent_start_dttm'] - block_df['admission_dttm']).dt.total_seconds()/3600
block_df['imv_to_death_days'] = (block_df['death_dttm'] - block_df['block_vent_start_dttm']).dt.total_seconds()/(24*3600)
block_df['adm_to_discharge_days'] = (block_df['discharge_dttm'] - block_df['admission_dttm']).dt.total_seconds()/(24*3600)
block_df['icu_to_imv_hours'] = (block_df['block_vent_start_dttm'] - block_df['icu_in_dttm']).dt.total_seconds()/(3600) #Positive if in ICU first before IMV.
block_df['Time_first_PT'] = (block_df['pt_post_imv_dttm'] - block_df['block_vent_start_dttm']).dt.total_seconds()/3600
block_df['Time_last_PT'] = (block_df['pt_pre_imv_dttm'] - block_df['block_vent_start_dttm']).dt.total_seconds()/3600

#Add in a dichotomized outcomes variables
block_df['pt_ever'] = block_df['pt_post_imv_dttm'].notna()
block_df['pt_post48_IMV'] = block_df['Time_first_PT'].notna() & (block_df['Time_first_PT'] <= 48)
block_df['pt_pre24_IMV'] = block_df['Time_last_PT'].notna() & (block_df['Time_last_PT'] >= -24)
block_df['yellow_ever'] = block_df['yellow_time_eligibility_2h'].notna()
block_df['yellow_post48_IMV'] = block_df['yellow_ever'] & (block_df['yellow_time_eligibility_2h'] <= 48)
block_df['extubated_at_pt'] = block_df['imv_to_end_hours'] <= block_df['Time_first_PT']
block_df['is_dead'] = block_df['death_dttm'].notna()
block_df['pt_between_ICU_IMV'] = block_df['Time_first_PT'] < (-1*block_df['icu_to_imv_hours'])

# Add Hospital mortality: TRUE if Death_dttm < discharge_dttm or (discharge category is hospice or dead) 
block_df["is_dead_hosp"] = (
    (block_df["death_dttm"] <= block_df["discharge_dttm"]) |
    (block_df["discharge_category"].str.lower().isin(["hospice", "expired"]))
)
#48 hour mortality (in grace period)
block_df['is_dead_2'] = (block_df['death_dttm'] - block_df['block_vent_start_dttm']).dt.total_seconds() <= (2*24*60*60)
#30-day mortality
block_df['is_dead_30'] = (block_df['death_dttm'] - block_df['block_vent_start_dttm']).dt.total_seconds() <= (30*24*60*60)
#365-day mortality
block_df['is_dead_365'] = (block_df['death_dttm'] - block_df['block_vent_start_dttm']).dt.total_seconds() <= (365*24*60*60)


# ## Clustering of Categorical Data
# Individualy and manually cluster all the categorical variables we want to cluster

# ### Language

# In[14]:


log("LANGUAGE PRE:")
log(block_df['language_category'].value_counts(dropna=False))
keep = {"English","Spanish"}
missing = {'Unknown or NA'}
set_mask = block_df['language_category'].notna() & (~block_df['language_category'].isin(missing))
block_df["language_category"] = np.where(block_df["language_category"].isin(keep), block_df["language_category"], "Other")
block_df["language_category"] = np.where(set_mask, block_df['language_category'], None)
log("LANGUAGE POST:")
log(block_df['language_category'].value_counts(dropna=False)) #log results


# ### Race

# In[15]:


log("RACE PRE:")
log(block_df['race_category'].value_counts(dropna=False))
keep = {"White", "Black or African American"}
set_mask = block_df['race_category'].notna() & (~block_df['race_category'].eq("Unknown"))
block_df["race_category"] = np.where(block_df["race_category"].isin(keep), block_df["race_category"], "Other")
block_df["race_category"] = np.where(set_mask, block_df['race_category'], None)
log("RACE POST:")
log(block_df['race_category'].value_counts(dropna=False)) #log results


# ### Ethnicity

# In[16]:


#This just converts "Unknown" to None for better missingness tracking.
set_mask = block_df['ethnicity_category'].notna() & (~block_df['ethnicity_category'].eq("Unknown"))
block_df["ethnicity_category"] = np.where(set_mask, block_df['ethnicity_category'], None)


# ### ICU Type

# In[17]:


log("ICU TYPE PRE:")
log(block_df['ICU_type'].value_counts(dropna=False))
mapping = {
    "general_icu": "Medical ICU",
    "medical_icu": "Medical ICU",
    "cardiac_icu": "Cardiac ICU",
    "cardiothoracic_surgical_icu": "Cardiac ICU",
    "mixed_cardiothoracic_icu": "Cardiac ICU",
    "cvicu_icu":"Cardiac ICU",
    "surgical_icu": "Surgical ICU",
    "burn_icu": "Other",
    "neurosurgical_icu":"Other",
    "neuro_icu":"Other",
    "mixed_neuro_icu":"Other"
}
block_df['ICU_type'] = block_df['ICU_type'].map(mapping)
block_df['ICU_type'] = np.where(block_df['ICU_type'].notna(), block_df['ICU_type'], None)
log("ICU TYPE POST:")
log(block_df['ICU_type'].value_counts(dropna=False))


# ### Admission Category

# In[18]:


log("ADMISSION PRE:")
log(block_df['admission_type_category'].value_counts(dropna=False))
mapping = {
    "ed": "Emergency Department",
    "facility":"Other",
    "osh": "Transfer",
    "direct": "Other",
    "elective": "Other",
    "other": "Other"
}
block_df['admission_type_category'] = block_df['admission_type_category'].map(mapping)
log("ADMISSION POST:")
log(block_df['admission_type_category'].value_counts(dropna=False))


# ### Discharge Category

# In[19]:


log("DISCHARGE PRE:")
log(block_df['discharge_category'].value_counts(dropna=False))
mapping = {
    "Home": "Home",
    "Group Home":"Home",
    "Against Medical Advice (AMA)": "Home",
    "Assisted Living": "Home",
    "Hospice": "Hospice",
    "Expired": "Expired",
    "Skilled Nursing Facility (SNF)": "Rehabilitation",
    "Acute Inpatient Rehab Facility": "Rehabilitation",
    "Psychiatric Hospital": "Other",
    "Acute Care Hospital": "Other",
    "Long Term Care Hospital (LTACH)": "Other",
    "Other": "Other",
    "Chemical Dependency":"Other",
    "Shelter":"Home",
    "Jail":"Home"
}
block_df['discharge_category'] = block_df['discharge_category'].map(mapping)
log("DISCHARGE POST:")
log(block_df['discharge_category'].value_counts(dropna=False))


# In[20]:


#Column check point
helper.missing_summary(block_df,f_name='block_df_3_end')


# ## Remove obersvations with prior PT order
# Remove any `encounter_block` from both `block_df` and `time_bin.df` where `pt_pre24_IMV` == `True`.

# In[21]:


#Exclusion criteria
log(f"To be excluded based on PT 24 hours prior to IMV: {sum(block_df['pt_pre24_IMV'])}")
block_df = block_df[~block_df['pt_pre24_IMV']]
time_bin.df = time_bin.df[time_bin.df['encounter_block'].isin(block_df['encounter_block'])]


# In[22]:


block_df[block_df['discharge_category'].isna() & block_df['discharge_dttm'].isna()].head()


# ## Organize Columns and Summarize

# In[23]:


import scipy.stats as stats

column_order = pd.read_csv(os.path.join("..","config","column_def.csv"))
my_cols = column_order['name'].tolist()
column_order = column_order.set_index('name')

block_df = block_df[my_cols]
path = os.path.join(output_folder,'intermediate',"block_for_stats.parquet")
block_df.to_parquet(path)
path = os.path.join(output_folder, 'intermediate',"block_for_stats.csv")
block_df.to_csv(path)
del path

#Convert outcome to a categorical column
early_col = "pt_post48_IMV"
block_df["early_PT"] = np.where(block_df[early_col], "early_PT", "no_early_PT")
n_total = block_df["encounter_block"].count()
n_early = block_df[early_col].sum()
n_not = n_total - n_early

#SMD calculator function
def calculate_smd(group1, group2):
    # Calculate means
    mean1 = np.mean(group1)
    mean2 = np.mean(group2)
    
    # Calculate variances
    var1 = np.var(group1, ddof=1)
    var2 = np.var(group2, ddof=1)
    
    # Calculate sample sizes
    n1 = len(group1)
    n2 = len(group2)
    
    # Calculate Pooled Standard Deviation
    # s_p = sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    pooled_sd = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    
    # Calculate Cohen's d
    smd = (mean1 - mean2) / pooled_sd
    return smd

file_path = os.path.join(output_folder,'final',"table1.csv")
if os.path.exists(file_path):
    os.remove(file_path)

with open(file_path, mode="w") as file:
    file.write(f",,Overall,Early PT, No Early PT, P-value/SMD, Missing")
    for col in my_cols:
        lab = column_order.loc[col,"description"]
        if col == 'encounter_block':
            file.write(f"\nN,,{n_total},{n_early},{n_not},")
        elif block_df[col].dtype == 'object' or pd.api.types.is_string_dtype(block_df[col]):
            file.write(f"\n{lab}")
            cats = block_df[col].dropna().unique()
            p_df = pd.pivot_table(block_df[["encounter_block",col,"early_PT"]], index=col, columns="early_PT", values="encounter_block", aggfunc='count')
            chi2_stat, p_value, dof, expected = stats.chi2_contingency(p_df.to_numpy())
            for cc in cats:
                if cc: #Because there is a NAN category
                    cc_all = round(100*p_df.loc[cc,].sum()/block_df[col].count(),1)
                    cc_early = round(100*p_df.loc[cc,'early_PT'].sum()/p_df['early_PT'].sum(),1)
                    cc_not = round(100*p_df.loc[cc,'no_early_PT'].sum()/p_df['no_early_PT'].sum(),1)
                    file.write(f"\n,{cc},{cc_all}%,{cc_early}%,{cc_not}%")
            file.write(f", {round(p_value,5)}")
        elif block_df[col].dtype == "bool" or block_df[col].dtype == "boolean":
            if block_df[col].sum(skipna=True) > 0:#Pivot table does not work if there are no true values
                sub_df = block_df[block_df[col].notna()]
                sub_df['flag'] = np.where(sub_df[col],"TRUE", "FALSE")
                p_df = pd.pivot_table(sub_df[["encounter_block","flag","early_PT"]], index="flag", columns="early_PT", values="encounter_block", aggfunc='count')
                chi2_stat, p_value, dof, expected = stats.chi2_contingency(p_df.to_numpy())
                cc_all = round(100*p_df.loc["TRUE",].sum()/sub_df[col].count(),1)
                cc_early = round(100*p_df.loc["TRUE",'early_PT'].sum()/p_df['early_PT'].sum(),1)
                cc_not = round(100*p_df.loc["TRUE",'no_early_PT'].sum()/p_df['no_early_PT'].sum(),1)
                file.write(f"\n{lab},,{cc_all}%,{cc_early}%,{cc_not}%,{round(p_value,5)}")
            else:
                file.write(f"\n{lab},,0.00%,0.00%,0.00%,N/A")
        elif pd.api.types.is_numeric_dtype(block_df[col]):
            cc_all = block_df[col].dropna()
            cc_early = block_df.loc[block_df[early_col], col]
            cc_early = cc_early.dropna()
            cc_not = block_df.loc[~block_df[early_col], col]
            cc_not = cc_not.dropna()
            SMD = calculate_smd(cc_early.tolist(), cc_not.tolist())
            file.write(f"\n{lab} (Med & IQR),,{round(cc_all.median(),2)}  ({round(cc_all.quantile(0.25),2)} - {round(cc_all.quantile(0.75),2)})")
            file.write(f",{round(cc_early.median(),2)}  ({round(cc_early.quantile(0.25),2)} - {round(cc_early.quantile(0.75),2)})")
            file.write(f",{round(cc_not.median(),2)}  ({round(cc_not.quantile(0.25),2)} - {round(cc_not.quantile(0.75),2)})")
            file.write(f",{round(SMD,5)}")
        else:
            file.write(f"\n{lab},ERROR,,,,,")
        #Missing data column
        mis_pct = round(100*sum(block_df[col].isna())/n_total,2)
        file.write(f",{mis_pct}%")
            


# ## CIF Graph

# In[24]:


import matplotlib.pyplot as plt

#Create a new lists
yellow_list = block_df['yellow_time_eligibility_2h'].sort_values().dropna()
pt_list = block_df['Time_first_PT'].sort_values().dropna()

# Plot CIF
plt.figure(figsize=(8, 6))
plt.step(yellow_list, np.arange(yellow_list.size), color='orange')
plt.step(pt_list, np.arange(pt_list.size), color='blue')
plt.xlabel("Hours from IMV Initiation")
plt.xlim(0, 72)
plt.ylabel("Encounters")
plt.title("Physioliogic Readiness versus PT Consult Order CIF")
plt.legend()

#Save and show
path = os.path.join(output_folder, 'final','graphs',"CIF_Yellow_v_PT.png")
plt.savefig(path)
plt.show()
plt.close()


# ## Time Bin Summary Graphs

# In[25]:


pt_time_bin_df = time_bin.df.groupby('bin_start')['pt_order'].agg('sum').reset_index()
pt_time_bin_df.sort_values(by='bin_start', inplace=True)

plt.figure(figsize=(8, 6))
plt.bar(pt_time_bin_df['bin_start'].tolist(), pt_time_bin_df['pt_order'].tolist(), color='blue')
plt.xlabel("Hours from IMV Initiation")
plt.xticks(np.arange(0, 48, 4))
plt.ylabel("Encounters")
plt.title("Early PT Consult Order Prevalence Over Time")

#Save and show
path = os.path.join(output_folder, 'final','graphs',"Time_bin_PT.png")
plt.savefig(path)
plt.show()
plt.close()


# ## Merging for Stats
# Create a merged block_df and time_bin.df to be used for stats.

# In[26]:


column_order = column_order.reset_index()
mask_cols = (column_order['name'] == 'encounter_block') | (column_order['covariate'] == 1) | (column_order['outcome'] == 1) | column_order['other'].notna()
stats_cols = column_order[mask_cols]
stats_df = block_df[stats_cols['name'].tolist()].copy()
stats_df = stats_df.merge(time_bin.df, on='encounter_block', how='inner').reset_index()
log(f"Stats data set contains {block_df['encounter_block'].nunique()} encounter_blocks.")
log(f"Stats data set contains {stats_df['encounter_block'].nunique()} encounter_blocks.")
path = os.path.join(output_folder, 'intermediate',"block_and_time_bins_for_stats.csv")
stats_df.to_csv(path)

