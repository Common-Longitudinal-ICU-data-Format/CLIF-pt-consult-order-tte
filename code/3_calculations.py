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
block_df = helper.load_data('output_folder','block_df_2_aggregated',folder='intermediate')

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
    df['red_resp_spo2_flag'] = ((df['spo2_min'] < 90) | df['spo2_min'].isna()).astype("Int64")
    df['red_map_flag'] = ((df['map_mean'] < 65) | df['map_mean'].isna()).astype("Int64")
    df['red_high_support_flag'] = ((df['ne_calc_last'] > 0.3) | (df['ne_calc_max'] > 0.3)).astype("Int64")
    df['red_hypertensive_flag'] = (
        (((df['sbp_max'] > 200) | (df['map_mean'] > 110)) &
        (df['red_med_flag'] == 1))
    ).astype("Int64")
    df['red_pulse_high_flag'] = (df['heart_rate_max'] > 150).astype("Int64")
    df['red_pulse_low_flag'] = ((df['heart_rate_min'] < 40) | df['heart_rate_min'].isna()).astype("Int64")

    # --- YELLOW flags ---
    df['yellow_resp_spo2_flag'] = ((df['spo2_min'] >= 90) | df['spo2_min'].isna()).astype("Int64")
    df['yellow_fio2_flag'] = (df['fio2_set_min'] > 0.6).astype("Int64")
    df['yellow_resp_rate_flag'] = (df['respiratory_rate_max'] > 30).astype("Int64")
    df['yellow_peep_flag'] = (df['peep_set_min'] > 10).astype("Int64")
    df['yellow_map_flag'] = ((df['map_mean'] >= 65) & (df['ne_calc_last'].between(0.1, 0.3))).astype("Int64")
    df['yellow_pulse_flag'] = (df['heart_rate_min'].between(120, 150)).astype("Int64")
    df['yellow_lactate_flag'] = (df['lactate_max'] > 4).astype("Int64")

    # --- GREEN flags ---
    df['green_resp_spo2_flag'] = ((df['spo2_min'] >= 90) | df['spo2_min'].isna()).astype("Int64")
    df['green_resp_rate_flag'] = ((df['respiratory_rate_max'] <= 30) | df['respiratory_rate_max'].isna()).astype("Int64")
    df['green_fio2_flag'] = ((df['fio2_set_min'] <= 0.6) | df['fio2_set_min'].isna()).astype("Int64")
    df['green_peep_flag'] = ((df['peep_set_min'] <= 10) | df['peep_set_min'].isna()).astype("Int64")
    df['green_map_flag'] = (((df['map_mean'] >= 65) & (df['ne_calc_last'] < 0.1)) | df['ne_calc_last'].isna()).astype("Int64")
    df['green_pulse_flag'] = ((df['heart_rate_min'] < 120) | df['heart_rate_min'].isna()).astype("Int64")
    df['green_lactate_flag'] = ((df['lactate_max'] <= 4) | df['lactate_max'].isna()).astype("Int64")
    df['green_hr_flag'] = ((df['heart_rate_min'] > 40) | df['heart_rate_min'].isna()).astype("Int64")

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
    ).astype("Int64")

    df['no_red'] = (
        ~(df['red_resp_spo2_flag'] | df['red_map_flag'] | df['red_high_support_flag'] |
          df['red_hypertensive_flag'] | df['red_pulse_high_flag'] | df['red_pulse_low_flag']) &
        _daytime
    ).astype("Int64")

    df['any_yellow'] = (
        (df['yellow_resp_spo2_flag'] | df['yellow_fio2_flag'] | df['yellow_resp_rate_flag'] |
         df['yellow_peep_flag'] | df['yellow_map_flag'] | df['yellow_pulse_flag'] |
         df['yellow_lactate_flag']) &
        _daytime
    ).astype("Int64")

    df['any_green'] = (
        (df['green_resp_spo2_flag'] | df['green_resp_rate_flag'] | df['green_fio2_flag'] |
         df['green_peep_flag'] | df['green_map_flag'] | df['green_pulse_flag'] |
         df['green_lactate_flag'] | df['green_hr_flag']) &
        _daytime
    ).astype("Int64")

    df['all_green'] = (
        df['green_resp_spo2_flag'] & df['green_resp_rate_flag'] & df['green_fio2_flag'] &
        df['green_peep_flag'] & df['green_map_flag'] & df['green_pulse_flag'] &
        df['green_lactate_flag'] & df['green_hr_flag'] & _daytime
    ).astype("Int64")

    df['all_green_all_hours'] = (
        df['green_resp_spo2_flag'] & df['green_resp_rate_flag'] & df['green_fio2_flag'] &
        df['green_peep_flag'] & df['green_map_flag'] & df['green_pulse_flag'] &
        df['green_lactate_flag'] & df['green_hr_flag'] & _base
    ).astype("Int64")

    df['all_green_weekday'] = (
        df['green_resp_spo2_flag'] & df['green_resp_rate_flag'] & df['green_fio2_flag'] &
        df['green_peep_flag'] & df['green_map_flag'] & df['green_pulse_flag'] &
        df['green_lactate_flag'] & df['green_hr_flag'] & _weekday
    ).astype("Int64")

    df['all_green_no_red'] = (
        df['green_resp_spo2_flag'] & df['green_resp_rate_flag'] & df['green_fio2_flag'] &
        df['green_peep_flag'] & df['green_map_flag'] & df['green_pulse_flag'] &
        df['green_lactate_flag'] & df['green_hr_flag'] & (df['any_red'] == 0) & _daytime
    ).astype("Int64")

    df['all_green_no_red_all_hours'] = (
        df['green_resp_spo2_flag'] & df['green_resp_rate_flag'] & df['green_fio2_flag'] &
        df['green_peep_flag'] & df['green_map_flag'] & df['green_pulse_flag'] &
        df['green_lactate_flag'] & df['green_hr_flag'] & (df['any_red'] == 0) & _base
    ).astype("Int64")

    df['all_green_no_red_weekday'] = (
        df['green_resp_spo2_flag'] & df['green_resp_rate_flag'] & df['green_fio2_flag'] &
        df['green_peep_flag'] & df['green_map_flag'] & df['green_pulse_flag'] &
        df['green_lactate_flag'] & df['green_hr_flag'] & (df['any_red'] == 0) & _weekday
    ).astype("Int64")

    df['all_green_no_red_yellow'] = (
        df['green_resp_spo2_flag'] & df['green_resp_rate_flag'] & df['green_fio2_flag'] &
        df['green_peep_flag'] & df['green_map_flag'] & df['green_pulse_flag'] &
        df['green_lactate_flag'] & df['green_hr_flag'] &
        (df['any_red'] == 0) & (df['any_yellow'] == 0) & _daytime
    ).astype("Int64")

    df['any_yellow_or_green_no_red'] = (
        (df['yellow_resp_spo2_flag'] | df['yellow_fio2_flag'] | df['yellow_resp_rate_flag'] |
         df['yellow_peep_flag'] | df['yellow_map_flag'] | df['yellow_pulse_flag'] |
         df['yellow_lactate_flag'] | df['green_resp_spo2_flag'] | df['green_resp_rate_flag'] |
         df['green_fio2_flag'] | df['green_peep_flag'] | df['green_map_flag'] |
         df['green_pulse_flag'] | df['green_lactate_flag'] | df['green_hr_flag']) &
        (df['any_red'] == 0) & _daytime
    ).astype("Int64")

    df['any_yellow_or_green_no_red_weekday'] = (
        (df['yellow_resp_spo2_flag'] | df['yellow_fio2_flag'] | df['yellow_resp_rate_flag'] |
         df['yellow_peep_flag'] | df['yellow_map_flag'] | df['yellow_pulse_flag'] |
         df['yellow_lactate_flag'] | df['green_resp_spo2_flag'] | df['green_resp_rate_flag'] |
         df['green_fio2_flag'] | df['green_peep_flag'] | df['green_map_flag'] |
         df['green_pulse_flag'] | df['green_lactate_flag'] | df['green_hr_flag']) &
        (df['any_red'] == 0) & _weekday
    ).astype("Int64")

    df['any_yellow_or_green_no_red_all_hours'] = (
        (df['yellow_resp_spo2_flag'] | df['yellow_fio2_flag'] | df['yellow_resp_rate_flag'] |
         df['yellow_peep_flag'] | df['yellow_map_flag'] | df['yellow_pulse_flag'] |
         df['yellow_lactate_flag'] | df['green_resp_spo2_flag'] | df['green_resp_rate_flag'] |
         df['green_fio2_flag'] | df['green_peep_flag'] | df['green_map_flag'] |
         df['green_pulse_flag'] | df['green_lactate_flag'] | df['green_hr_flag']) &
        (df['any_red'] == 0) & _base
    ).astype("Int64")

    df['green_resp_flag'] = (
        df['green_resp_spo2_flag'] & df['green_resp_rate_flag'] &
        df['green_fio2_flag'] & df['green_peep_flag'] & _daytime
    ).astype("Int64")

    df['green_cardio_flag'] = (
        df['green_map_flag'] & df['green_pulse_flag'] &
        df['green_lactate_flag'] & df['green_hr_flag'] & _daytime
    ).astype("Int64")

    df['yellow_resp_flag'] = (
        (df['yellow_resp_spo2_flag'] | df['yellow_fio2_flag'] | df['yellow_resp_rate_flag'] |
         df['yellow_peep_flag'] | df['green_resp_spo2_flag'] | df['green_resp_rate_flag'] |
         df['green_fio2_flag'] | df['green_peep_flag']) &
        (df['any_red'] == 0) & _daytime
    ).astype("Int64")

    df['yellow_cardio_flag'] = (
        (df['yellow_map_flag'] | df['yellow_pulse_flag'] | df['yellow_lactate_flag'] |
         df['green_map_flag'] | df['green_pulse_flag'] | df['green_lactate_flag'] | df['green_hr_flag']) &
        (df['any_red'] == 0) & _daytime
    ).astype("Int64")

    df['yellow_all_green'] = (df['all_green_no_red'] & (df['any_yellow'] == 0)).astype("Int64")
    df['yellow_not_all_green'] = (df['any_yellow_or_green_no_red'] & (df['all_green_no_red'] == 0)).astype("Int64")

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
last_vent_df['time_from_vent'] = last_vent_df['time_from_vent'].astype("Int64")
last_vent_df.rename(columns={'time_from_vent':'last_hour_on_vent'}, inplace=True)
block_df = pd.merge(
    block_df,
    last_vent_df,
    on='encounter_block',
    how='left'
)

#Get an 1 for patients alive at 28-days.
block_df['alive28'] = block_df['death_dttm'].isna() | ((block_df['death_dttm'] - block_df['block_vent_start_dttm']).dt.total_seconds() >= (28*24*60*60))
block_df['alive28'] = block_df['alive28'].astype("Int64")

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


# ## Remove obersvations with prior PT order
# Remove any `encounter_block` from both `block_df` and `time_bin.df` where `pt_pre24_IMV` == `True`.

# In[20]:


#Exclusion criteria
log(f"To be excluded based on PT 24 hours prior to IMV: {sum(block_df['pt_pre24_IMV'])}")
block_df = block_df[~block_df['pt_pre24_IMV']]
time_bin.df = time_bin.df[time_bin.df['encounter_block'].isin(block_df['encounter_block'])]


# ## Save

# In[21]:


#Save
path = os.path.join(output_folder,'intermediate',"block_df_3_end.parquet")
block_df.to_parquet(path)
#Missing summary
helper.missing_summary(block_df,f_name='block_df_3_end')

