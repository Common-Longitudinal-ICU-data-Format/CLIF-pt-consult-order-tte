#!/usr/bin/env python
# coding: utf-8

# # MIMIC-CLIF Early PT Consults Comparative Analysis
# # Step 2: Data Gathering
# - Load clinical data from CLIF tables.
# - Add data to the block data frame, hourly blocks and create time bins.
# - This is mostly collecting the data and aggregating it in preparation for data analysis which happens separately.

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

#To use MIMIC which we should not need once all data is in CLIF format
use_mimic = 'mimic' in config['site_name'].lower()

#Load ClifOrchestrator
co = clifpy.ClifOrchestrator(
    data_directory=config['clif_folder'],
    filetype=config['file_type'],
    timezone=config['time_zone'],
    output_directory=output_folder
)


# In[2]:


#Create Logger
_logger = logging.getLogger('clif_01')
_logger.setLevel(logging.INFO)
_logger.handlers.clear()

_log_dir = os.path.join(output_folder,'logs',f'{config['site_name']}_02_data_gathering_log.txt')
_fh = logging.FileHandler(_log_dir, mode='w')
_fh.setFormatter(logging.Formatter('%(asctime)s | %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
_logger.addHandler(_fh)

_ch = logging.StreamHandler()
_ch.setFormatter(logging.Formatter('%(message)s'))
_logger.addHandler(_ch)

def log(*args, **kwargs):
    _msg = ' '.join(str(a) for a in args)
    _logger.info(_msg)

log(f"=== CLIF Pipeline 02: Data Gathering ===")
log(f"Site: {config['site_name']}")


# ## Load Prior Data

# In[3]:


#Load block data
block_df = helper.load_data('output_folder','block_df_1_end',folder='intermediate')
#Load Encounter Mapping
enc_map = helper.load_data('output_folder','encounter_mapping',folder='intermediate')
#Water Resp Support Table
rs_waterfall = helper.load_data('output_folder','respiratory_support_waterfall',folder='intermediate')
rs_waterfall = rs_waterfall[rs_waterfall['encounter_block'].isin(enc_map['encounter_block'])]


# In[4]:


#Initiliaze CLIF Orchestrator
#Filtered by hospitalization_id
#Apply stitching
hosp_list = enc_map['hospitalization_id'].tolist()
co.initialize(
    tables=['hospitalization','vitals'],
    filters={'hospitalization':{'hospitalization_id':hosp_list},
            'vitals':{'hospitalization_id':hosp_list}}
)
co.hospitalization.df = co.hospitalization.df.merge(enc_map[['hospitalization_id','encounter_block']], on='hospitalization_id', how='left')
co.encounter_mapping = enc_map

#Load Respiratory support but replace it with waterfall created in first data set.
co.load_table('respiratory_support')
co.respiratory_support.df = rs_waterfall


# ### Add weight_kg

# In[5]:


#Load the Vitals CLIF tables
clifpy.utils.apply_outlier_handling(co.vitals)

#Remove the ones not in cohort, add encounter_block and vent start time as reference, remove bad time values.
weight_df = co.vitals.df[co.vitals.df['vital_category'] ==  'weight_kg']
weight_df = weight_df.merge(enc_map, on='hospitalization_id', how='right').reset_index()
mweight_df = weight_df[weight_df['recorded_dttm'].notna()]

#Calculate time windows from block_vent_start_dttm
weight_df['time_diff'] = (weight_df['recorded_dttm'] - weight_df['block_vent_start_dttm']).dt.total_seconds()

# Define whether measurement is before or after vent_start_time
weight_df['before_vent_start'] = (weight_df['time_diff'] <= 0).astype("Int64")

# Calculate absolute time difference
weight_df['abs_time_diff'] = weight_df['time_diff'].abs()

# Sort data to prioritize measurements before vent start and closest in time
weight_df = weight_df.sort_values(['encounter_block', 'vital_category', 'before_vent_start', 'abs_time_diff'], 
                                    ascending=[True, True, False, True])

# Drop duplicates to keep the closest measurement for each vital_category per encounter block
weight_df = weight_df.drop_duplicates(subset=['encounter_block'], keep='first')
weight_df.rename(columns={'vital_value':'weight_kg'},inplace=True)

#Add back to block
block_df = pd.merge(
    block_df,
    weight_df[['encounter_block','weight_kg']],
    on='encounter_block',
    how='left'
)

del weight_df

log(f"Missing Weight_kg at block level: {sum(block_df['weight_kg'].isna())}")


# ## Columns and Filters

# In[6]:


rst_interest = [
    'device_category',
    'tracheostomy',
    'fio2_set',
    'lpm_set',
    'resp_rate_set',
    'peep_set',
    'resp_rate_obs'
]

vitals_required_columns = [
    'hospitalization_id',
    'recorded_dttm',
    'vital_category',
    'vital_value'
]
vitals_of_interest = ['heart_rate', 'respiratory_rate', 'sbp', 'dbp', 'map', 'spo2']

labs_required_columns = [
    'hospitalization_id',
    'lab_collect_dttm',
    'lab_result_dttm',
    'lab_category',
    'lab_value',
    'lab_value_numeric'
]
labs_of_interest = [
    'creatinine',
    'lactate',
    'platelet_count',
    'po2_arterial',
    'bilirubin_total'
]

meds_required_columns = [
    'hospitalization_id',
    'admin_dttm',
    'med_name',
    'med_category',
    'med_dose',
    'med_dose_unit'
]
meds_of_interest = [
    'norepinephrine', 'epinephrine', 'phenylephrine', 'vasopressin',
    'dopamine', 'angiotensin', 'nicardipine', 'nitroprusside',
    'clevidipine', 'cisatracurium', 'vecuronium', 'rocuronium',
    'metaraminol','dobutamine'
]

pat_ass_required_columns = [
    'hospitalization_id',
    'recorded_dttm',
    'assessment_category',
    'numerical_value',
    'categorical_value',
]

pat_ass_of_interest = [
    'braden_mobility',
    'RASS',
    'cam_total',
    'gcs_total',
]


# ## Load and Process Continous Medications

# In[7]:


#Load with filter
_meds_filters = {
    'hospitalization_id': hosp_list,
    'med_category': meds_of_interest
}
co.load_table('medication_admin_continuous', columns=meds_required_columns, filters=_meds_filters)
#Apply Stitching
co.medication_admin_continuous.df = co.medication_admin_continuous.df.merge(enc_map[['hospitalization_id','encounter_block']], on='hospitalization_id', how='left')

#Remove Outliers
clifpy.utils.apply_outlier_handling(co.medication_admin_continuous)

#Add weight coloumn
co.medication_admin_continuous.df = co.medication_admin_continuous.df.merge(block_df[['encounter_block','weight_kg']], on='encounter_block', how='left')

#Convert Units
_med_units_preferred = {
    'norepinephrine':"mcg/kg/min",
    'epinephrine':"mcg/kg/min",
    'phenylephrine':"mcg/kg/min",
    'vasopressin':"mcg/kg/min",
    'dopamine':"mcg/kg/min",
    'angiotensin':"mcg/kg/min",
    'metaraminol':"mcg/kg/min",
    'dobutamine':"mcg/kg/min"
}

from clifpy.utils.unit_converter import convert_dose_units_by_med_category
_converted_meds, _summary = convert_dose_units_by_med_category(
    co.medication_admin_continuous.df,
    preferred_units = _med_units_preferred,
    override = True
)
co.medication_admin_continuous.df = _converted_meds
co.medication_admin_continuous.df['med_dose'] = co.medication_admin_continuous.df['med_dose_converted']
co.medication_admin_continuous.df['med_dose_unit'] = co.medication_admin_continuous.df['med_dose_unit_converted']


# In[8]:


# Create summary tables
_summary_meds = co.medication_admin_continuous.df.groupby('med_category').agg(
    total_N=('med_category', 'size'),
    min=('med_dose', 'min'),
    max=('med_dose', 'max'),
    first_quantile=('med_dose', lambda x: x.quantile(0.25)),
    second_quantile=('med_dose', lambda x: x.quantile(0.5)),
    third_quantile=('med_dose', lambda x: x.quantile(0.75)),
    missing_values=('med_dose', lambda x: x.isna().sum())
).reset_index()
_summary_meds.to_csv(f'{output_folder}/final/summary_meds_by_category.csv', index=False)

_summary_meds_cat_dose = co.medication_admin_continuous.df.groupby(['med_category', 'med_dose_unit']).agg(
    total_N=('med_category', 'size'),
    min=('med_dose', 'min'),
    max=('med_dose', 'max'),
    first_quantile=('med_dose', lambda x: x.quantile(0.25)),
    second_quantile=('med_dose', lambda x: x.quantile(0.5)),
    third_quantile=('med_dose', lambda x: x.quantile(0.75)),
    missing_values=('med_dose', lambda x: x.isna().sum())
).reset_index()
_summary_meds_cat_dose.to_csv(f'{output_folder}/final/summary_meds_by_category_dose_units.csv', index=False)

# Diagnostic: Check which groups have all NaN values
log("Groups with all NaN med_dose values:")
for (_med_category, _med_dose_unit), _group in co.medication_admin_continuous.df.groupby(['med_category', 'med_dose_unit']):
    if _group['med_dose'].isna().all():
        log(f"  {_med_category} - {_med_dose_unit}: {len(_group)} rows, all NaN")

del _summary_meds, _summary_meds_cat_dose

#Filter non-usable data
_meds_filtered = co.medication_admin_continuous.df[co.medication_admin_continuous.df['med_dose'].notna()].copy()
def has_per_hour_or_min(unit):
    if pd.isnull(unit):
        return False
    unit = unit.lower()
    return '/hr' in unit or '/min' in unit
_meds_filtered = _meds_filtered[_meds_filtered['med_dose_unit'].apply(has_per_hour_or_min)]


# ### Pressor Dose Conversion

# In[9]:


# ── Norepinephrine equivalent calculation ──
_meds_list = [
    "norepinephrine", "epinephrine", "phenylephrine",
    "vasopressin", "dopamine","metaraminol",
    "angiotensin"
]
_pressor_mask = _meds_filtered['med_category'].isin(_meds_list)
_non_pressors = _meds_filtered[~_pressor_mask].copy()
_ne_df = _meds_filtered[_pressor_mask].copy()

for _med in _meds_list:
    if _med not in _ne_df['med_category'].unique():
        log(f"{_med} is not in the dataset.")
    else:
        log(f"{_med} is in the dataset.")

#Apply conversion factors
#Norepi to epi is 1:1.
_ne_df['med_dose'] = np.where(_ne_df['med_category'] == "phenylephrine", _ne_df['med_dose']/10, _ne_df['med_dose'])
_ne_df['med_dose'] = np.where(_ne_df['med_category'] == "dopamine", _ne_df['med_dose']/100, _ne_df['med_dose'])
_ne_df['med_dose'] = np.where(_ne_df['med_category'] == "metaraminol", _ne_df['med_dose']/8, _ne_df['med_dose'])
_ne_df['med_dose'] = np.where(_ne_df['med_category'] == "vasopressin", _ne_df['med_dose']*2.5, _ne_df['med_dose'])
_ne_df['med_dose'] = np.where(_ne_df['med_category'] == "angiotensin", _ne_df['med_dose']*10, _ne_df['med_dose'])
_ne_df['med_category'] = "norepinephrine"

co.medication_admin_continuous.df = pd.concat([_ne_df,_non_pressors])

log(f'Converted {_ne_df.shape[0]} pressors found to NE equivalents. Meds_df now shape: {co.medication_admin_continuous.df.shape}')
del _meds_filtered, _ne_df, _non_pressors


# ### Antihypertensives & Paralytics
# This is a bit unusual. We do not care about dose conversion because ultimately we just want to create a binary on/off flag. We are doing this inside the ClifOrchestrator data set because we want to take advantage of the optimizations of the *create_wide_dataset* function.

# In[10]:


# Convert IV continuous antihypertensives to all be the same drug.
_meds_list = [
    'nicardipine', 'nitroprusside', 'clevidipine'
]
_mask = co.medication_admin_continuous.df['med_category'].isin(_meds_list)
_sub_df = co.medication_admin_continuous.df[_mask]

for _med in _meds_list:
    if _med not in _sub_df['med_category'].unique():
        log(f"{_med} is not in the dataset.")
    else:
        log(f"{_med} is in the dataset.")

co.medication_admin_continuous.df['med_category'] = np.where(_mask, "nitroprusside", co.medication_admin_continuous.df['med_category'])
log(f"Found {sum(_mask)} antihypertensives and converted them all to nitroprusside.")

del _sub_df


# In[11]:


# Convert IV continuous paralytics to all be the same drug.
_meds_list = [
    'cisatracurium', 'vecuronium', 'rocuronium'
]
_mask = co.medication_admin_continuous.df['med_category'].isin(_meds_list)
_sub_df = co.medication_admin_continuous.df[_mask]

for _med in _meds_list:
    if _med not in _sub_df['med_category'].unique():
        log(f"{_med} is not in the dataset.")
    else:
        log(f"{_med} is in the dataset.")

co.medication_admin_continuous.df['med_category'] = np.where(_mask, 'rocuronium', co.medication_admin_continuous.df['med_category'])
log(f"Found {sum(_mask)} paralytics and converted them all to rocuronium.")

del _sub_df


# In[12]:


#Remove all the excess columns created
meds_required_columns.append('encounter_block')
co.medication_admin_continuous.df = co.medication_admin_continuous.df[meds_required_columns]


# ## Load Rest of CLIF Tables

# In[13]:


#Load with filter
_filters = {
    'hospitalization_id': hosp_list,
    'assessment_category': pat_ass_of_interest
}
co.load_table('patient_assessments', columns=pat_ass_required_columns, filters=_filters)
#Remove Outliers
clifpy.utils.apply_outlier_handling(co.patient_assessments)
log(f"Loaded patient_assessments table. Size: {co.patient_assessments.df.shape[0]}, hosp_ids: {co.patient_assessments.df['hospitalization_id'].nunique()}")

#Load with filter
_filters = {
    'hospitalization_id': hosp_list,
    'lab_category': labs_of_interest
}
co.load_table('labs', columns=labs_required_columns, filters=_filters)
#Remove Outliers
clifpy.utils.apply_outlier_handling(co.labs)
log(f"Loaded patient_assessments table. Size: {co.patient_assessments.df.shape[0]}, hosp_ids: {co.patient_assessments.df['hospitalization_id'].nunique()}")


# ## Create Wide Data Set

# In[14]:


_cohort = block_df[['encounter_block','patient_id','block_vent_start_dttm','block_last_vital_dttm']]
_cohort.rename(columns={'block_vent_start_dttm':'start_time','block_last_vital_dttm':'end_time'}, inplace=True)
_cohort['start_time'] = _cohort['start_time'] - pd.Timedelta(hours=24)

co.create_wide_dataset(
    tables_to_load = ['vitals','respiratory_support','labs','medication_admin_continuous','patient_assessments'],
    category_filters = {
        'vitals':vitals_of_interest,
        'respiratory_support':rst_interest,
        'labs':labs_of_interest,
        'medication_admin_continuous':meds_of_interest,
        'patient_assessments':pat_ass_of_interest
    },
    encounter_blocks = enc_map['encounter_block'].tolist(),
    cohort_df = _cohort
)
del _cohort
    


# ### Add MAP & RR where missing

# In[15]:


_add_map = co.wide_df['sbp'].notna() & co.wide_df['dbp'].notna() & co.wide_df['map'].isna()
co.wide_df['map'] = np.where(_add_map, 0.333*co.wide_df['sbp'] + 0.666*co.wide_df['dbp'], co.wide_df['map'])
log(f'Added {sum(_add_map)} missing MAP')

_add_rr = co.wide_df['resp_rate_obs'].notna() & co.wide_df['respiratory_rate'].isna()
co.wide_df['respiratory_rate'] = np.where(_add_rr, co.wide_df['resp_rate_obs'], co.wide_df['respiratory_rate'])
log(f'Added {sum(_add_rr)} respiratory rate rows from resp_supoprt to vitals.')


# ## SOFA Score
# Use the newly created wide_df

# In[16]:


#SOFA scores

#Cohort
sofa_input_df = enc_map.copy()
sofa_input_df['start_time'] = sofa_input_df['block_vent_start_dttm']
sofa_input_df['end_time'] = sofa_input_df['block_vent_start_dttm'] + pd.Timedelta(hours=24)
sofa_input_df.drop(columns=['block_vent_start_dttm'], inplace=True)

#Wide DF 
'''
It expects these med columns.
The conversion was already done prior to wide_df creation.
'''
sofa_wide = co.wide_df.copy()
_med_expected = ["norepinephrine","epinephrine", "dopamine","dobutamine"]
for col in _med_expected:
    if col in sofa_wide.columns:
        sofa_wide[f"{col}_mcg_kg_min"] =  sofa_wide[col]
    else:
        sofa_wide[f"{col}_mcg_kg_min"] = 0

sofa_df = co.compute_sofa_scores(
    wide_df = sofa_wide,
    cohort_df = sofa_input_df, # id, start_time, end_time  (local time)
    id_name="encounter_block",
    extremal_type = 'worst',
    fill_na_scores_with_zero=True,
    remove_outliers=False
)

'''
#THIS WAS ATTEMPTED TO ADDRESS MISSINGNESS BUT NOT USED
#Above fill_na=True was used instead.
#Re-assign missingness (the function assigns 0 to missing values but we want N/A).
sofa_missing_mask = (sofa_df['sofa_cv_97'].isna() |\
                     sofa_df['sofa_coag'].isna() |\
                     sofa_df['sofa_liver'].isna() |\
                     sofa_df['sofa_resp'].isna() |\
                     sofa_df['sofa_cns'].isna() |\
                     sofa_df['sofa_renal'].isna() )
#Keep missing values under the philosphy (no value = no monitoring = normal)
sofa_df['sofa_total'] = np.where(sofa_missing_mask, np.nan,sofa_df['sofa_total'])
log(f"Encounters with at least one missing value for SOFA: {sum(sofa_missing_mask)}")
'''
sofa_df.rename(columns={'sofa_total':"sofa_0_24h"},inplace=True)

#SOFA score merge to blocks
block_df = block_df.merge(
    sofa_df[['encounter_block','sofa_0_24h']],
    on='encounter_block',
    how='left'
)
log(f"Encounters SOFA score missing in final data: {sum(block_df['sofa_0_24h'].isna())}")

del sofa_input_df, sofa_df, sofa_wide


# ## Convert Wide to Hourly Data

# In[17]:


agg_plan = {
    'max':['rocuronium','nitroprusside','norepinephrine','tracheostomy','respiratory_rate','heart_rate','sbp','lactate'],
    'min':['spo2','fio2_set','heart_rate','peep_set'],
    'last':['norepinephrine','device_category','rocuronium','nitroprusside'],
    'mean':['map']
}

_temp_hourly_df = co.convert_wide_to_hourly(agg_plan,
                                            id_name='encounter_block',
                                            hourly_window=1,
                                            fill_gaps=True)


# In[18]:


log('Hourly_df created with columns:\n',_temp_hourly_df.dtypes)
log('Encounters in hourly_df:\n',_temp_hourly_df['encounter_block'].nunique())


# In[19]:


#Missing summary before filling anything in.
helper.missing_summary(_temp_hourly_df, f_name='hourly_df_2_clifpy_raw')

#Load Hourly_Blocks object
_temp_hourly_df = helper.convert_datetime_columns(_temp_hourly_df) #Currently clifpy implementation strips time zone but retains value.
_temp_hourly_df = _temp_hourly_df.merge(block_df[['encounter_block','block_vent_start_dttm']], on='encounter_block',how='left')
_temp_hourly_df['time_from_vent'] = np.ceil((_temp_hourly_df['window_end_dttm'] - _temp_hourly_df['block_vent_start_dttm']).dt.total_seconds()/3600)
_temp_hourly_df['time_from_vent'] = _temp_hourly_df['time_from_vent'].astype("Int64")
hourly = helper.hourly_blocks(in_df=_temp_hourly_df)


# ### Fix 'last' missingness
# PROBLEM:
# Conversion from wide to hourly as of now has a bug where 'last' takes the last value regardless of wether it is NA or not.
# This is generating missingness.
# 
# SOLUTION:
# Individualy fix those after creation of the hourly_df.
# 
# An issue on GitHub has been created for this, however as of April 30, 2026 it remains open. [clifpy issue 131](https://github.com/Common-Longitudinal-ICU-data-Format/clifpy/issues/131)

# In[20]:


for last_col in agg_plan['last']:
    new_name = f"{last_col}_last"
    log(f"Fixing missingness of {new_name}")
    last_df = co.wide_df[co.wide_df[last_col].notna()].copy()
    last_df = last_df.merge(block_df[['encounter_block','block_vent_start_dttm']], on='encounter_block',how='left')
    last_df['time_from_vent'] = np.ceil((last_df['event_time'] - last_df['block_vent_start_dttm']).dt.total_seconds()/3600)
    last_df['time_from_vent'] = last_df['time_from_vent'].astype("Int64")
    log(f"-- Sort step")
    last_df = last_df.sort_values(by=['encounter_block','event_time'])[['encounter_block','time_from_vent',last_col]]
    log(f"-- Agg step")
    last_df = last_df.groupby(['encounter_block','time_from_vent'])[last_col].agg('last').reset_index()
    last_df.rename(columns={last_col:new_name},inplace=True)
    hourly.df.drop(columns=[new_name],inplace=True)
    log(f"-- Merge step")
    hourly.df = hourly.df.merge(last_df, on=['encounter_block','time_from_vent'],how='left')


# In[21]:


co.wide_df.dtypes


# ### Forward and back fill

# In[22]:


#Forward Fill from last for max rows
hourly.hourly_fill('tracheostomy_max','ffill')
hourly.hourly_fill('tracheostomy_max',False)
hourly.df['tracheostomy_max'] = hourly.df['tracheostomy_max'].astype("Int64")
inter = list(set(agg_plan['max']) & set(agg_plan['last']))
for col in inter:
    log(f'Filling hourly for {col} _max and _last. Empty cells {sum(hourly.df[f'{col}_max'].isna()) + sum(hourly.df[f'{col}_last'].isna())}')
    hourly.hourly_fill(f'{col}_last','ffill') #Forward fill first
    hourly.df[f'{col}_max'] = np.where(hourly.df[f'{col}_max'].isna(), hourly.df[f'{col}_last'],hourly.df[f'{col}_max'])
    
#Forward Fill from last for min rows
inter = list(set(agg_plan['min']) & set(agg_plan['last']))
for col in inter:
    log(f'Filling hourly for {col} _min and _last. Empty cells {sum(hourly.df[f'{col}_min'].isna()) + sum(hourly.df[f'{col}_last'].isna())}')
    hourly.hourly_fill(f'{col}_last','ffill') #Forward fill first, note _last does not get bfilled.
    hourly.df[f'{col}_min'] = np.where(hourly.df[f'{col}_min'].isna(), hourly.df[f'{col}_last'],hourly.df[f'{col}_min'])
    hourly.hourly_fill(f'{col}_min','bffill') #To back fill anything left.
    log(f'Left {col}_last empty cells: {sum(hourly.df[f'{col}_last'].isna())}')

#Create vent maker and fill (1 if imv mentioned, 0 if anything, fill in otherwise.
hourly.hourly_fill('device_category_last','ffill')
hourly.df['hourly_on_vent'] = hourly.df['device_category_last'] == 'imv'
hourly.df['hourly_on_vent'] = hourly.df['hourly_on_vent'].astype("Int64")
hourly.hourly_fill('hourly_on_vent','ffill')

#Med flags specifically should get back filled with zeros.
hourly.hourly_fill(f'norepinephrine_last', 0)
hourly.hourly_fill(f'norepinephrine_max', 0)
hourly.hourly_fill(f'rocuronium_last', 0)
hourly.hourly_fill(f'rocuronium_max', 0)

#Other random fills
hourly.hourly_fill('respiratory_rate_max','bffill')
hourly.hourly_fill('spo2_min','bffill')
hourly.hourly_fill('fio2_set_min','bffill')
hourly.hourly_fill('sbp_max','bffill')
hourly.hourly_fill('peep_set_min','bffill')


# ### Create some binary flags and rename come columns

# In[23]:


_col_rename = {
    'norepinephrine_last':'ne_calc_last',
    'norepinephrine_max':'ne_calc_max',
    'nitroprusside_max':'red_med_flag',
    'rocuronium_max':'paralytics_flag'
}
hourly.df.rename(columns=_col_rename, inplace=True)
hourly.df['red_med_flag'] = hourly.df['red_med_flag'] > 0
hourly.df['red_med_flag'] = hourly.df['red_med_flag'].astype("Int64")
hourly.df['paralytics_flag'] = hourly.df['paralytics_flag'] > 0
hourly.df['paralytics_flag'] = hourly.df['paralytics_flag'].astype("Int64")


# ### RASS
# For some reason patient assessments do not seem to properly load into the wide data set so will add them the hourly manually.

# In[24]:


#Load assessments
co.patient_assessments.df = co.patient_assessments.df.merge(enc_map, on='hospitalization_id', how='right').reset_index()
co.patient_assessments.df['time_diff'] = co.patient_assessments.df['recorded_dttm'] - co.patient_assessments.df['block_vent_start_dttm']
rass_df = co.patient_assessments.df[co.patient_assessments.df['assessment_category'] == 'RASS'].copy()
rass_df.rename(columns={'numerical_value': 'RASS'}, inplace=True)

#Remove the ones not in cohort, add encounter_block and vent start time as reference
rass_df['time_from_vent'] = hourly.calc_time_from_vent(rass_df['time_diff'])

hourly.addto_blocks(rass_df,'RASS',agg_func='min', fill_with='bffill')


# In[25]:


#Define coma
hourly.df['coma'] = hourly.df['RASS_min'] < -2
hourly.df['coma'] = hourly.df['coma'].astype("Int64")
del rass_df


# ### Save Hourly Data

# In[26]:


#Remove negative time.
hourly.df = hourly.df[hourly.df['time_from_vent'] > 0]
#Summary and Save
hourly.save(suffix='_two')
log('Completed hourly.df and saved sumary')


# ## Time Bins and Other Aggregation

# In[27]:


#Create Time Bin Object
time_bin = helper.time_bins(in_eb = block_df[['encounter_block','block_vent_start_dttm']].copy().drop_duplicates())

#Add death event to time bins
death_df = block_df[['encounter_block','block_vent_start_dttm','death_dttm']].copy().drop_duplicates()
death_df['time_diff'] = death_df['death_dttm'] - death_df['block_vent_start_dttm']
time_bin.add_event(death_df[['encounter_block','time_diff']], 'death')
del death_df


# In[28]:


#Add PT consult
pt_df = block_df[['encounter_block','block_vent_start_dttm','pt_post_imv_dttm']].copy().drop_duplicates()
pt_df['time_diff'] = pt_df['pt_post_imv_dttm'] - pt_df['block_vent_start_dttm']
pt_df.sort_values(by=['encounter_block','time_diff'], inplace=True)
pt_df['time_bin'] = time_bin.classify_time_bin(pt_df['time_diff'])
pt_df = pt_df[pt_df['time_bin'].notna()][['encounter_block','time_bin','time_diff']]
pt_df.drop_duplicates(subset=['encounter_block'], keep='first', inplace=True) 
time_bin.add_event(pt_df[['encounter_block','time_diff']], 'pt_order')
#Now add the pt_now column to help with censoring models.
pt_df['pt_now'] = 1
time_bin.df = pd.merge(
    time_bin.df,
    pt_df,
    on=['encounter_block','time_bin'],
    how='left')
time_bin.bin_sort_fill('pt_now',0)
del pt_df


# In[29]:


#Create time columns for the wide data set
co.wide_df = co.wide_df.merge(block_df[['encounter_block','block_vent_start_dttm']].copy().drop_duplicates(), on='encounter_block',how='left')
co.wide_df['time_diff'] = (co.wide_df['event_time'] - co.wide_df['block_vent_start_dttm']).dt.total_seconds()/3600
co.wide_df['time_bin'] = time_bin.classify_time_bin(co.wide_df['time_diff'])

log(f'Finished creating time bins. Shape: {time_bin.df.shape}')


# ### Vitals
# 
# We have looked into pulling vitals data from the MIMIC ED database to see if it improves missingness of pre-intubations vitals.
# It added pre-intubation vital signs for about 10% of encounters but we still had 40% missingness so this portion of the code was removed.
# It can still be found in old versions of the code as needed.

# In[30]:


##HEART RATE
hr_df = co.wide_df[co.wide_df['heart_rate'].notna()]

#Default aggregation for the first 24 hours summary variable.
block_df = block_df.merge(helper.aggregate_by_time(hr_df, 'heart_rate'), on='encounter_block', how='left')
log(f"Missing heart rate at block level: {sum(block_df['heart_rate_0_24h_mean'].isna())}")

#Time Bin data
time_bin.gather_time_bins(hr_df,'heart_rate', fill_with='ffill')
log(f"Missing heart rate at time_bin level: {sum(time_bin.df['heart_rate_mean'].isna())}")

del hr_df


# In[31]:


##MAP
map_df = co.wide_df[co.wide_df['map'].notna()]

#Default aggregation for the first 24 hours
block_df = block_df.merge(helper.aggregate_by_time(map_df, 'map'), on='encounter_block', how='left')
log(f"Missing MAP at block level: {sum(block_df['map_0_24h_mean'].isna())}")

#Time Bin data
time_bin.gather_time_bins(map_df,'map', fill_with='ffill')
log(f"Missing MAP at time_bin level: {sum(time_bin.df['map_mean'].isna())}")

del map_df


# ### Respiratory Support Table

# In[32]:


#FiO2
fio2_df = co.wide_df[co.wide_df['fio2_set'].notna()]

#Default aggregation for the first 24 hours
block_df = block_df.merge(helper.aggregate_by_time(fio2_df, 'fio2_set'), on='encounter_block', how='left')
log(f"Missing FiO2 at block level: {sum(block_df['fio2_set_0_24h_mean'].isna())}")

#Time Bin data
time_bin.gather_time_bins(fio2_df,'fio2_set', fill_with='ffill')
log(f"Missing FiO2 at time_bin level: {sum(time_bin.df['fio2_set_mean'].isna())}")

del fio2_df


# In[33]:


#PEEP
peep_df = co.wide_df[co.wide_df['peep_set'].notna()]

#Default aggregation for the first 24 hours
block_df = block_df.merge(helper.aggregate_by_time(peep_df, 'peep_set'), on='encounter_block', how='left')
log(f"Missing FiO2 at block level: {sum(block_df['peep_set_0_24h_mean'].isna())}")

#Time Bin data
time_bin.gather_time_bins(peep_df,'peep_set', fill_with='ffill')
log(f"Missing FiO2 at time_bin level: {sum(time_bin.df['peep_set_mean'].isna())}")

del peep_df


# ### Patient Assessments

# In[34]:


co.patient_assessments.df['time_bin'] = time_bin.classify_time_bin(co.patient_assessments.df['time_diff'])


# In[35]:


#RASS
rass_df = co.patient_assessments.df[co.patient_assessments.df['assessment_category'] == 'RASS'].copy()
rass_df.rename(columns={'numerical_value': 'RASS'}, inplace=True)

#Default aggregation for the first 24 hours
block_df = block_df.merge(helper.aggregate_by_time(rass_df, 'RASS', agg_func='min'), on='encounter_block', how='left')
log(f"Missing RASS at block level: {sum(block_df['RASS_0_24h_min'].isna())}")

#Time Bin data
time_bin.gather_time_bins(rass_df,'RASS', agg_func='min', fill_with='ffill')
log(f"Missing RASS at time bin level: {sum(time_bin.df['RASS_min'].isna())}")

del rass_df


# In[36]:


#Braden Mobility Score

#Filter ass_df for braden_mob only
braden_df = co.patient_assessments.df[co.patient_assessments.df['assessment_category'] == 'braden_mobility'].copy()
braden_df.rename(columns={'numerical_value': 'braden_mobility'}, inplace=True)

#Default aggregation for the first 24 hours
block_df = block_df.merge(helper.aggregate_by_time(braden_df, 'braden_mobility',agg_func='max'), on='encounter_block', how='left')
log(f"Missing braden (first 24) at block level: {sum(block_df['braden_mobility_0_24h_max'].isna())}")

#Redefine time-diff to look at last braden within 24 hours of ICU discharge.
#This will give the max braden within 24 hours of ICU discharge.
braden_df = pd.merge(braden_df.copy(),block_df[['encounter_block','icu_out_dttm']], on='encounter_block', how='left')
braden_df['time_diff'] = braden_df['recorded_dttm'] - braden_df['icu_out_dttm']
block_df = pd.merge(
    block_df,
    helper.aggregate_by_time(braden_df,
                             'braden_mobility',
                             -24,0,
                             agg_func='max'),
    on='encounter_block',
    how='left'
)
block_df.rename(columns={'braden_mobility_-24_0h_max':'braden_mobility_last'},inplace=True)
log(f"Missing braden (ICU last) at block level: {sum(block_df['braden_mobility_last'].isna())}")

del braden_df


# In[37]:


#CAM-ICU
#Filter ass_df for braden_mob only
cam_df = co.patient_assessments.df[co.patient_assessments.df['assessment_category'] == 'cam_total'].copy()
cam_df['cam_icu'] = np.where(cam_df['categorical_value'].notna() & (cam_df['categorical_value'] == 'Positive'),1,0) #Convert categorical to binary

#Default aggregation for the first 24 hours
block_df = block_df.merge(helper.aggregate_by_time(cam_df, 'cam_icu', agg_func='flag'), on='encounter_block', how='left')
log(f"Missing CAM-ICU at block level: {sum(block_df['cam_icu_0_24h_flag'].isna())}")

del cam_df


# ## Elixhauser
# Using package comorbidipy with Quan et al mappings and Van Walraven weights
# 
# Outputs both unadjusted and age adjusted.
# 
# Note that this uses either ICD 9 or ICD 10 codes for any given encounter_block. There was a hard switch at some point so there should actually NO encounters with both ICD codes mixed.

# In[38]:


#If we need to convert admission year from MIMIC date to real dates.
if use_mimic:    
    log('Using MIMIC for admission year.')
    #load
    patient_mimic_df = helper.load_data("mimic","patients", folder='hosp', type='csv.gz')
    patient_mimic_df.rename(columns={'subject_id': 'patient_id'}, inplace=True)
    patient_mimic_df['patient_id'] = patient_mimic_df['patient_id'].astype(str)
    
    #filter
    patient_mimic_df = patient_mimic_df[patient_mimic_df['patient_id'].isin( block_df['patient_id'] )]
    print(f"Unique MIMIC patient id: {patient_mimic_df['patient_id'].nunique()}")
    print(f"Unique Block patient id: {block_df['patient_id'].nunique()}")
    
    #merge
    merged_patient_df = pd.merge(
        block_df[['encounter_block','patient_id','block_vent_start_dttm']],
        patient_mimic_df,
        on='patient_id',
        how='left'
    )

    #Admission Year (based on estimate from anchor year
    merged_patient_df["anchor_year_group"] = merged_patient_df["anchor_year_group"].str.slice(0, 4).astype("Int64") + 1 #Convert achor year group to an integer
    merged_patient_df["admission_year"] = merged_patient_df['block_vent_start_dttm'].dt.year.astype("Int64") - merged_patient_df["anchor_year"] + merged_patient_df["anchor_year_group"]
    
    block_df = block_df.merge(merged_patient_df[['admission_year','encounter_block']],on='encounter_block', how='left')
    
    del patient_mimic_df, merged_patient_df
else:
    block_df["admission_year"] = block_df['block_vent_start_dttm'].dt.year.astype("Int64")

log(f"Encounters with admission year missing in final data: {sum(block_df['admission_year'].isna())}")
print(f"Block Length: {len(block_df)}")
print(f"Unique Encounter Block: {block_df['encounter_block'].nunique()}")


# In[39]:


#Load with filter
_filters = {
    'hospitalization_id': hosp_list
}
co.load_table('hospital_diagnosis', filters=_filters)

#filter out primary diagnosis
diag_df = co.hospital_diagnosis.df[co.hospital_diagnosis.df['diagnosis_primary'] != 1] #Remove primary diagnosis

#Add EB
diag_df = diag_df.merge(enc_map, on='hospitalization_id', how='left')

#Merge with block data
diag_df = pd.merge(
        block_df[['encounter_block','admission_year','age']],
        diag_df[['encounter_block','diagnosis_code','diagnosis_code_format']],
        on='encounter_block',
        how='left'
    ).drop_duplicates().reset_index()

#Create separate data frames for each ICD version
diag_9_df = diag_df[diag_df['diagnosis_code_format'] == 'ICD9CM'].drop_duplicates()
diag_10_df = diag_df[diag_df['diagnosis_code_format'] == 'ICD10CM'].drop_duplicates()
diag_9_df = diag_9_df[['encounter_block','diagnosis_code','age']]
diag_10_df = diag_10_df[['encounter_block','diagnosis_code','age']]

#Check how many encnounterblocks are in both datasets and report it as an issue.
duplicate_icds_n = sum(diag_9_df['encounter_block'].isin(diag_10_df['encounter_block']))
log(f"Number of encnounters with comorbidities in both ICD codes: {duplicate_icds_n}")

#LAST SAVING POINT
path = os.path.join(output_folder,"intermediate", "diag_codes.parquet")
diag_df.to_parquet(path)
#del path, diag_df, diag_9_df, diag_10_df


# ## Save

# In[40]:


#SAVING POINT
path = os.path.join(output_folder,'intermediate', "block_df_2_aggregated.parquet")
block_df.to_parquet(path)
del path
#Censor out dead data
time_bin.remove_based_on_censor('death', keep_first=True)
#Save (which will save the data as well as a summary of it)
time_bin.save(suffix='_step_2')
del time_bin #To save memory

