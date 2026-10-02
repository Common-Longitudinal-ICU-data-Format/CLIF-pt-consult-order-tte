#!/usr/bin/env python
# coding: utf-8

# # MIMIC-CLIF Early PT Consults Comparative Analysis
# ## Step 1: Cohort Identification and Initial Data Gathering
# 
# - Uses CLIF dataset for plans of future roll out to CLIF
# - Cohort identification based on age, hours on vent and trach status.
# - Demographic data collection for cohort.
# - Hospitalization stitching.
# - Respiratory data restructuring.

# ## Setup

# In[ ]:


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

#MIMIC needs years are shifted so cannot be used for filteing.
use_mimic = 'mimic' in config['site_name'].lower()

#Admission year window, applied at every site so the pooled analysis covers one common era.
#Defaults are used when a site's config.json predates these keys. See config/README.md.
year_min = int(config.get('year_min', 2018))
year_max = int(config.get('year_max', 2024))


# In[ ]:


#output_folders
# NOTE: Uses print() not log() — logger is created AFTER this cell
# to ensure the log file points to the fresh output directory.
print("=== Output Folder Management ===")

# Create fresh output directory structure
print(f"Creating fresh output directory structure...")
os.makedirs(output_folder, exist_ok=True)
os.makedirs(os.path.join(output_folder,'logs'), exist_ok=True)
os.makedirs(os.path.join(output_folder,'final'), exist_ok=True)
os.makedirs(os.path.join(output_folder,'intermediate'), exist_ok=True)

# Create graphs subfolder
graphs_folder = os.path.join(output_folder,'final','graphs')
os.makedirs(graphs_folder, exist_ok=True)

print(f"Output directory structure ready:")
print(f"   {output_folder}/")
print(f"   +-- log/")
print(f"   +-- final/")
print(f"   |   +-- graphs/")
print(f"   +-- intermediate/")

#Create Logger
_logger = logging.getLogger('clif_01')
_logger.setLevel(logging.INFO)
_logger.handlers.clear()

_log_dir = os.path.join(output_folder,'logs',f'{config['site_name']}_01_cohort_log.txt')
_fh = logging.FileHandler(_log_dir, mode='w')
_fh.setFormatter(logging.Formatter('%(asctime)s | %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
_logger.addHandler(_fh)

_ch = logging.StreamHandler()
_ch.setFormatter(logging.Formatter('%(message)s'))
_logger.addHandler(_ch)

def log(*args, **kwargs):
    _msg = ' '.join(str(a) for a in args)
    _logger.info(_msg)

#Load ClifOrchestrator
co = clifpy.ClifOrchestrator(
    data_directory=config['clif_folder'],
    filetype=config['file_type'],
    timezone=config['time_zone'],
    output_directory=output_folder
)

log(f"=== CLIF Pipeline 01: Cohort Identification ===")
log(f"Site: {config['site_name']}")


# In[ ]:


#Load Clif Tables
co.initialize(
    tables=['patient', 'hospitalization', 'adt'],
)
log(f"Total Number of unique encounters in the hospitalization table: {co.hospitalization.df['hospitalization_id'].nunique()}")


# ## Cohort Identification 
# ### (A) Age Filter

# In[ ]:


log("\n=== STEP A: Filter by age ===\n")
_age_mask = (co.hospitalization.df['age_at_admission'] >= 18)
co.hospitalization.df = co.hospitalization.df[_age_mask]
del _age_mask

strobe_ab = {}
strobe_ab['A_after_age_filter'] = co.hospitalization.df['hospitalization_id'].nunique()
strobe_ab['A_unique_patients'] = co.hospitalization.df['patient_id'].nunique()
log(f"Number of unique hospitalizations after age filter: {strobe_ab['A_after_age_filter']}")
log(f"Number of unique patients after age filter: {strobe_ab['A_unique_patients']}")
log("\nMissing values in admission_dttm:", co.hospitalization.df['admission_dttm'].isna().sum())
log("Missing values in discharge_dttm:", co.hospitalization.df['discharge_dttm'].isna().sum())


# ### (B) Stitch Hospitalizations

# In[ ]:


log("\n=== STEP B: Stitch encounters ===\n")
_cohort_ids = co.hospitalization.df['hospitalization_id'].unique().tolist()
co.adt.df = co.adt.df[co.adt.df['hospitalization_id'].isin(_cohort_ids)]
del _cohort_ids
co.stitch_time_interval = 6
co.run_stitch_encounters()

strobe_ab['B_before_stitching'] = co.hospitalization.df['hospitalization_id'].nunique()
strobe_ab['B_after_stitching'] = co.hospitalization.df['encounter_block'].nunique()
strobe_ab['B_stitched_hosp_ids'] = strobe_ab['B_before_stitching'] - strobe_ab['B_after_stitching']
log(f"Number of unique hospitalizations before stitching: {strobe_ab['B_before_stitching']}")
log(f"Number of unique encounter blocks after stitching: {strobe_ab['B_after_stitching']}")
log(f"Number of linked hospitalization ids: {strobe_ab['B_before_stitching'] - strobe_ab['B_after_stitching']}")


# ### (B2) Year Filter

# In[ ]:


log("\n=== STEP B2: Admission year filter ===\n")

co.hospitalization.df['admission_year'] = co.hospitalization.df['admission_dttm'].dt.year.astype("Int64")
log(f"Mininum year before filtering: {co.hospitalization.df['admission_year'].min()}")
log(f"Maximum year before filtering: {co.hospitalization.df['admission_year'].max()}")

if use_mimic:
    ##No year filtering for MIMIC, include all.
    _eb_inyear = co.hospitalization.df['encounter_block'].drop_duplicates()
else:
    _eb_inyear = co.hospitalization.df[ co.hospitalization.df['admission_year'].between(year_min, year_max) ]['encounter_block'].drop_duplicates()

co.hospitalization.df = co.hospitalization.df[ co.hospitalization.df['encounter_block'].isin(_eb_inyear) ]

strobe_ab['B2_encounter_blocks_in_year'] = co.hospitalization.df['encounter_block'].nunique()
strobe_ab['B2_encounter_blocks_out_of_year'] = strobe_ab['B_after_stitching'] - strobe_ab['B2_encounter_blocks_in_year']
log(f"Encounter blocks in year: {strobe_ab['B2_encounter_blocks_in_year']}")
log(f"Encounter blocks out of year: {strobe_ab['B2_encounter_blocks_in_year']}")
log(f"Mininum year after filtering: {co.hospitalization.df['admission_year'].min()}")
log(f"Maximum year after filtering: {co.hospitalization.df['admission_year'].max()}")

del _eb_inyear


# ### (C) Identify Ventilator Usage
# - Filter first by any hospitalization that has a single IMV reportedd.
# - Create waterfall / hourly blocks of respiratory support data. See clifpy documentation for details of everything this entails.
# - Identify potential extubations based on `lpm_set` not being empty.
# - Impute missing FiO2 values

# In[ ]:


log("\n=== STEP C: Load & process respiratory support => Apply Waterfall & Identify IMV usage ===\n")

strobe_c = {}

co.initialize(tables=['respiratory_support'],
             filters={'respiratory_support':{'hospitalization_id':co.hospitalization.df['hospitalization_id'].tolist()}})
clifpy.utils.apply_outlier_handling(co.respiratory_support) #Remove outliers per CLIF standards
co.respiratory_support.df['device_category'] = co.respiratory_support.df['device_category'].str.lower()
_resp_support = co.respiratory_support.df.merge(co.encounter_mapping, on='hospitalization_id', how='left')
_resp_support = _resp_support.sort_values(['encounter_block', 'recorded_dttm'])

#Get only respiratory data for enconuters that have at least one IMV mention
_imv_mask = _resp_support['device_category'].str.contains("imv", case=False, na=False)
_imv_ids = _resp_support[_imv_mask][['encounter_block']].drop_duplicates()
_resp_support = _resp_support[
    _resp_support['encounter_block'].isin(_imv_ids['encounter_block'])
].reset_index(drop=True)
strobe_c['C_imv_encounter_blocks'] = _imv_ids['encounter_block'].nunique()
log(f"Number of encounter blocks with IMV: {strobe_c['C_imv_encounter_blocks']}")
del _imv_mask, _imv_ids

#Waterfall process
from clifpy.tables.respiratory_support import RespiratorySupport
_rs = RespiratorySupport(data=_resp_support)
rs_waterfall = _rs.waterfall(id_col="encounter_block", verbose=True, return_dataframe=True)
#Since we used EB to create waterfall the hosp_id's were not fowardfilled.
rs_waterfall['hospitalization_id'] = rs_waterfall.groupby('encounter_block')['hospitalization_id'].ffill()
log(f"Number of rows in respiratory support waterfall: {rs_waterfall.shape[0]}")

#Fixes an issue of the waterfall function using a different time-zone object although the time zone is still the same.
rs_waterfall = helper.convert_datetime_columns(rs_waterfall)


# In[ ]:


#Describe extend of missing device category problem
log(f"Number of rows missing device_category: {sum(rs_waterfall['device_category'].isna())}")
fix_mask = rs_waterfall['device_category'].isna() & (rs_waterfall['lpm_set'] > 0) & (rs_waterfall['lpm_set'].notna())
log(f"Of those rows, number with lpm_set: {sum(fix_mask)}")

#Identify extubations
def count_extubations(series:pd.Series):
    """
    input: clif_respiratory_support:device_category data series.
    output: int
    Counts the number of extubations as imv to non-imv transitions.
    """
    series_bool = series.dropna().str.contains("imv",case=False)
    ext_bool = (~series_bool) & (series_bool.shift(periods=1, fill_value=True))
    return sum(ext_bool)
def extubations_describe(rs_table:pd.DataFrame):
    """
    input: clif_respiratory_support:device_category data series.
    output: str
    Describe extubation counts.
    """
    extubation_df = rs_table.groupby('encounter_block')['device_category'].apply(count_extubations).reset_index()
    extubation_df.rename(columns={'device_category':'ext_n'},inplace=True)
    return f"""EXUBATION ANALYSIS:
            Total encounters: {extubation_df['encounter_block'].nunique()}
            Total extubations: {extubation_df['ext_n'].sum()}
            Encounters without extubations: {sum(extubation_df['ext_n'] < 1)}
            Encounters with > 1 extubations: {sum(extubation_df['ext_n'] > 1)}
            Encounters with > 2 extubations: {sum(extubation_df['ext_n'] > 2)}
            99th percentile: {extubation_df['ext_n'].quantile(.99):.2f}"""
log("PRE FIX EXTUBATION COUNTER")
log(extubations_describe(rs_waterfall))

#NOTE: The applied rules below are hiarcheal, meaning as one rule is applied these rows are removed from requiring a fix.

#Apply trach rule
trach_rule = fix_mask & (rs_waterfall['tracheostomy'] == 'True')
rs_waterfall.loc[trach_rule, 'device_category'] = 'trach collar'
fix_mask = fix_mask & ~trach_rule #Remove from list of needed to be fixed
log(f"Implied trach collar based on tracheostomy flag {sum(trach_rule)}")

#Apply NIPPV rule, assumes the lpm_set would be for bleed in.
nippv_rule = fix_mask & ((rs_waterfall['pressure_support_set'] > 0) | (rs_waterfall['peak_inspiratory_pressure_set'] > 0))
rs_waterfall.loc[nippv_rule, 'device_category'] = 'nippv'
fix_mask = fix_mask & ~nippv_rule #Remove from list of needed to be fixed
log(f"Implied NIPPV based on inspiratory pressure {sum(nippv_rule)}")

#Apply CPAP rule, assumes the lpm_set would be for bleed in.
cpap_rule = fix_mask & (rs_waterfall['peep_set'] > 0)
rs_waterfall.loc[cpap_rule, 'device_category'] = 'cpap'
fix_mask = fix_mask & ~cpap_rule #Remove from list of needed to be fixed
log(f"Implied CPAP based on peep pressure {sum(cpap_rule)}")

#Apply HFNC rule.
hfnc_rule = fix_mask & (rs_waterfall['lpm_set'] > 20)
rs_waterfall.loc[hfnc_rule, 'device_category'] = 'high flow nc'
fix_mask = fix_mask & ~hfnc_rule #Remove from list of needed to be fixed
log(f"Implied HFNC based on >20Lpm {sum(hfnc_rule)}")

#Apply Face Mask rule.
fm_rule = fix_mask & ((rs_waterfall['lpm_set'] > 10) | rs_waterfall['fio2_set'].notna())
rs_waterfall.loc[fm_rule, 'device_category'] = 'face mask'
fix_mask = fix_mask & ~fm_rule #Remove from list of needed to be fixed
log(f"Implied face mask based on >10Lpm or Fio2 set {sum(fm_rule)}")

#Apply LF rule. To all remainig (note implied that they do not meet any of the criteria above)
rs_waterfall.loc[fix_mask, 'device_category'] = 'nasa cannula'
log(f"Implied nasal cannula {sum(fix_mask)}")

log("POST FIX EXTUBATION COUNTER")
log(extubations_describe(rs_waterfall))


# In[ ]:


#Impute FIO2
'''
Function taken from https://github.com/Common-Longitudinal-ICU-data-Format/CLIF-eligibility-for-mobilization/blob/main/code/pyCLIF.py
'''
def impute_fio2_from_nasal_cannula_flow(df):
    """
    Impute missing FiO2 values for nasal cannula based on oxygen flow rate (LPM)
    using the standard clinical conversion table.
    
    Parameters:
    df (pd.DataFrame): DataFrame with device_category, lpm_set, and fio2_set columns
    
    Returns:
    pd.DataFrame: DataFrame with updated fio2_set values
    """
    # Create lookup table from the clinical standard (your image)
    fio2_lookup = {
        1: 0.24,   # 1 L/min → 24%
        2: 0.28,   # 2 L/min → 28%  
        3: 0.32,   # 3 L/min → 32%
        4: 0.36,   # 4 L/min → 36%
        5: 0.40,   # 5 L/min → 40%
        6: 0.44,   # 6 L/min → 44%
        7: 0.48,   # 7 L/min → 48%
        8: 0.52,   # 8 L/min → 52%
        9: 0.56,   # 9 L/min → 56%
        10: 0.60   # 10 L/min → 60%
    }
    
    # Create mask for rows that need imputation
    nasal_cannula_mask = (
        (df['device_category'] == 'nasal cannula') &
        (df['fio2_set'].isna()) &  # Missing FiO2
        (df['lpm_set'].notna()) &  # Have LPM value
        (df['lpm_set'] >= 1) &     # Within lookup range  
        (df['lpm_set'] <= 10) &    # Within lookup range
        (df['lpm_set'] == df['lpm_set'].round())  # Integer values only
    )
    
    # Apply the lookup for eligible rows
    df.loc[nasal_cannula_mask, 'fio2_set'] = df.loc[nasal_cannula_mask, 'lpm_set'].map(fio2_lookup)
    
    # Report what was imputed
    n_imputed = nasal_cannula_mask.sum()
    if n_imputed > 0:
        print(f"[OK] Imputed FiO2 for {n_imputed:,} nasal cannula rows using LPM lookup table")
        
        # Show breakdown by LPM
        imputed_breakdown = df[nasal_cannula_mask]['lpm_set'].value_counts().sort_index()
        print("   Breakdown by LPM:")
        for lpm, count in imputed_breakdown.items():
            fio2_pct = int(fio2_lookup[lpm] * 100)
            print(f"   {lpm}L/min → {fio2_pct}%: {count:,} rows")
    else:
        print(" No nasal cannula rows needed FiO2 imputation")
    
    return df

rs_waterfall = impute_fio2_from_nasal_cannula_flow(rs_waterfall)
#On vent flag
rs_waterfall['on_vent'] = np.where(rs_waterfall['device_category'].str.contains("imv", case=False, na=False), 1, 0)

log("Missing values in recorded_dttm:", rs_waterfall['recorded_dttm'].isna().sum())


# ### (D) Determine ventilation times (start/end) at encounter block level

# In[ ]:


log("\n=== STEP D: Determine ventilation times (start/end) at encounter block level ===\n")

strobe_d = {}

#IMV only data frame.
_resp_stitched_imv = rs_waterfall[rs_waterfall['on_vent'] == 1]

# at the block level
block_vent_times = _resp_stitched_imv.groupby('encounter_block', dropna=True).agg(
    block_vent_start_dttm=('recorded_dttm', 'min'),
    block_vent_end_dttm=('recorded_dttm', 'max')
).reset_index()

_block_same_vent = block_vent_times[block_vent_times['block_vent_start_dttm'] == block_vent_times['block_vent_end_dttm']].copy()
block_vent_times = block_vent_times[block_vent_times['block_vent_start_dttm'] != block_vent_times['block_vent_end_dttm']].copy()

strobe_d['D_blocks_with_valid_vent'] = block_vent_times['encounter_block'].nunique()
strobe_d['D_blocks_with_same_vent_start_end'] = _block_same_vent['encounter_block'].nunique()
log(f"Unique encounter blocks with valid IMV start/end: {strobe_d['D_blocks_with_valid_vent']}")


# In[ ]:


#Quick aside the start the block_df data frame and to filter our the CO data.
block_df = pd.merge(block_vent_times,
                    co.hospitalization.df[['encounter_block','patient_id']],
                    on='encounter_block',
                    how='left').drop_duplicates(subset='encounter_block').reset_index()

#Filter out
_eb_list = block_vent_times['encounter_block'].unique().tolist()
co.hospitalization.df = co.hospitalization.df[co.hospitalization.df['encounter_block'].isin(_eb_list)]
co.adt.df = co.adt.df[co.adt.df['encounter_block'].isin(_eb_list)]
co.encounter_mapping = co.encounter_mapping[co.encounter_mapping['encounter_block'].isin(_eb_list)]
co.patient.df = co.patient.df[co.patient.df['patient_id'].isin(co.hospitalization.df['patient_id'])]
rs_waterfall = rs_waterfall[rs_waterfall['encounter_block'].isin(_eb_list)]
    
log(f"01_cohort: ADULT and IMV for >0 hours : Block Length: {len(block_df)}, Encounter Blocks {block_df['encounter_block'].nunique()}")


# ### (E) Intubation episodes identification

# In[ ]:


log("\n=== STEP E: Intubation episodes identification ===\n")

# --- keep only rows with a known device category ---
_resp_df = rs_waterfall.loc[
    rs_waterfall['device_category'].notna()
    & (rs_waterfall['device_category'].astype(str).str.strip() != ''),
    ['encounter_block', 'recorded_dttm', 'device_category','tracheostomy','on_vent']
].copy()

_resp_df = _resp_df.sort_values(['encounter_block', 'recorded_dttm'], kind='stable').reset_index(drop=True)

# --- flag IMV vs. not, and label contiguous runs of the same state ---
prev_state = _resp_df.groupby('encounter_block')['on_vent'].shift()
_resp_df['run_id'] = (_resp_df['on_vent'] != prev_state).astype(bool).cumsum()   # first row of each encounter starts a new run

# --- collapse each run to its first/last timestamp ---
runs = (
    _resp_df.groupby('run_id', sort=True)
      .agg(encounter_block=('encounter_block', 'first'),
           is_imv=('on_vent', 'first'),
           run_start=('recorded_dttm', 'first'),
           run_last=('recorded_dttm', 'last'),
           tracheostomy=('tracheostomy', 'first'))
      .reset_index(drop=True)
)
log(f"Total respiratory support RUNS found: {runs.shape[0]}")
log(f"Unique encounter blocks in RUNS: {runs['encounter_block'].nunique()}")

# start of the next run in the same encounter = extubation time
runs['next_run_start'] = runs.groupby('encounter_block')['run_start'].shift(-1)

imv_episodes = runs[ runs['is_imv'] == 1 ].copy()
imv_episodes['vent_start_dttm'] = imv_episodes['run_start']
imv_episodes['vent_end_dttm'] = imv_episodes['next_run_start'].fillna(imv_episodes['run_last'])

imv_episodes = (
    imv_episodes[['encounter_block', 'vent_start_dttm', 'vent_end_dttm','tracheostomy']]
    .sort_values(['encounter_block', 'vent_start_dttm'])
    .drop_duplicates()
    .reset_index(drop=True)
)

log(f"Total intubation episodes found: {imv_episodes.shape[0]}")
log(f"Unique encounter blocks with >0 intubations episodes: {imv_episodes['encounter_block'].nunique()}")
    
#del _resp_df, runs


# In[ ]:


#Merge with the overall block vent initiation
imv_episodes = pd.merge(block_df[['encounter_block','block_vent_start_dttm']],
                        imv_episodes,
                        on='encounter_block',
                        how='left')

#Max out vent_end_dttm to block_vent_start_dttm + 72 hours
imv_episodes['72h'] = imv_episodes['block_vent_start_dttm'] + pd.Timedelta(hours=72)
imv_episodes = imv_episodes[ imv_episodes['vent_start_dttm'] <= imv_episodes['72h'] ]
imv_episodes['vent_end_dttm'] = imv_episodes[['72h','vent_end_dttm']].min(axis=1)

#Calculate the number of hours (by the filter above this would only include hours in the first 72 hours)
imv_episodes['vent_hours'] = (imv_episodes['vent_end_dttm'] - imv_episodes['vent_start_dttm']).dt.total_seconds()/3600
imv_hours_df = imv_episodes.groupby('encounter_block').agg(
    vent_hours = ('vent_hours','sum'),
    trach_at_start = ('tracheostomy','first')
).reset_index()

log(f"Unique encounter blocks with intubations data: {imv_hours_df['encounter_block'].nunique()}")

del imv_episodes


# ### (F) Exclusion Criteria

# In[ ]:


# Exclude blocks with imv for less than 4 hours
_blocks_under_4 = imv_hours_df[ imv_hours_df['vent_hours'] < 4 ]
imv_hours_df = imv_hours_df[ imv_hours_df['vent_hours'] >= 4 ]

strobe_excl = {}
strobe_excl['F_blocks_with_vent_4_or_more'] = imv_hours_df['encounter_block'].nunique()
strobe_excl['F_blocks_with_vent_less_than_4'] = len(_blocks_under_4)
log(f"Unique encounter blocks with IMV >=4 hours: {strobe_excl['F_blocks_with_vent_4_or_more']}")
log(f"Excluded {len(_blocks_under_4)} encounter blocks with <4 vent hours in first 72 hours of intubation.\n")

# Exclude blocks with trach at the time of intubation
_blocks_with_trach_at_intubation = imv_hours_df['trach_at_start'].astype(int)
imv_hours_df = imv_hours_df[ _blocks_with_trach_at_intubation == 0 ]
strobe_excl['F_final_blocks_with_trach_at_intubation'] = sum(_blocks_with_trach_at_intubation)
strobe_excl['F_final_blocks_without_trach_at_intubation'] = imv_hours_df['encounter_block'].nunique()
log(f"Blocks with trach at intubation: {strobe_excl['F_final_blocks_with_trach_at_intubation']}")
log(f"Cohort size in hourly blocks: {strobe_excl['F_final_blocks_without_trach_at_intubation']}")

del _blocks_under_4, _blocks_with_trach_at_intubation


# ### Save Cohort Sample
# - Apply filter as noted above
# - Save progress so far including encounter stitching and cohort sample.

# In[ ]:


#Filter out from final cohort
_eb_list = imv_hours_df['encounter_block'].unique().tolist()
block_df = block_df[block_df['encounter_block'].isin(_eb_list)]
co.hospitalization.df = co.hospitalization.df[co.hospitalization.df['encounter_block'].isin(_eb_list)]
co.adt.df = co.adt.df[co.adt.df['encounter_block'].isin(_eb_list)]
co.patient.df = co.patient.df[co.patient.df['patient_id'].isin(co.hospitalization.df['patient_id'])]
rs_waterfall = rs_waterfall[rs_waterfall['encounter_block'].isin(_eb_list)]

#Merge with encounter block data
enc_map = pd.merge(co.encounter_mapping,
                    block_df[['encounter_block','block_vent_start_dttm']],
                    on='encounter_block',
                    how='right')
enc_map = enc_map[enc_map['encounter_block'].isin(_eb_list)]

#Save progress
path = os.path.join(output_folder,'intermediate','block_df_1_creation.parquet')
block_df.to_parquet(path)
del path

#Save it for later use
path = os.path.join(output_folder,'intermediate','respiratory_support_waterfall.parquet')
rs_waterfall.to_parquet(path)
del path

log(f"01_cohort: FINAL COHORT: Block Length: {len(block_df)}, Encounter Blocks {block_df['encounter_block'].nunique()}")


# ## Start Data Collection
# 
# Collects the data for the cohort from data frames we already have loaded up.
# 
# - Add patient data
# - Add admission data
# - Add discharge data
# - Add ADT data

# ### (A) Patient Data

# In[ ]:


_columns_of_interest = ['patient_id','race_category','ethnicity_category','sex_category','death_dttm','language_category']
block_df = pd.merge(block_df,
                    co.patient.df[_columns_of_interest],
                    on='patient_id',
                    how='left')

for _col in _columns_of_interest:
    log(f"Blocks with {_col} missing: {sum(block_df[_col].isna())}")


# ### (B) Admission & Discharge Data
# 
# Unfortunately the CLIF admission_type_category mapping from MIMIC-CLIF does not include whether or not the patient came as an OSH transfer. So we will need to get data direct from MIMIC to determine OSH transfers if desired.

# In[ ]:


_hosp_df = co.hospitalization.df.copy()
_hosp_df.rename(columns = {'age_at_admission':'age'}, inplace=True)
#Sort to organize hospitalization by order of admission_dttm
_hosp_df = _hosp_df.sort_values(by = ['encounter_block','admission_dttm'])
#Aggregate each column by the following rules.
_hosp_df = _hosp_df.groupby('encounter_block').agg({'admission_dttm':'min',
                                                    'admission_type_category':'first',
                                                    'age':'min',
                                                    'discharge_dttm':'max',
                                                    'discharge_category':'last'}).reset_index()
#Merge back with block_data
block_df = pd.merge(block_df,
                    _hosp_df,
                    on='encounter_block',
                    how='left')

_columns_of_interest = ['encounter_block','admission_dttm','age','admission_type_category','discharge_dttm','discharge_category']
for _col in _columns_of_interest:
    log(f"Blocks with {_col} missing: {sum(block_df[_col].isna())}")


# ### (C) ADT Data
# - ICU in and out time.

# In[ ]:


#Merge with encounter block
merged_adt_df = pd.merge (
    block_df[['encounter_block','block_vent_start_dttm']],
    co.adt.df,
    on='encounter_block',
    how = 'inner'
)

#Remove location_categories not considered to be transfers out (ie. procedures and radiology)
units_not_transfer = ['radiology','procedural','dialysis']
merged_adt_df = merged_adt_df[ ~ merged_adt_df['location_category'].str.lower().isin(units_not_transfer) ]

#Keep only units that were more than an hour or an ICU
adt_mask = ((merged_adt_df["out_dttm"] - merged_adt_df["in_dttm"]).dt.total_seconds() > 3600) | merged_adt_df['location_category'].str.contains('icu',case=False)
merged_adt_df = merged_adt_df[adt_mask]

#Keep only data after vent and sort by in_dttm
merged_adt_df = merged_adt_df[merged_adt_df['out_dttm'] > merged_adt_df['block_vent_start_dttm']].copy()
merged_adt_df = merged_adt_df.sort_values(['encounter_block','in_dttm'], ascending=True)


# In[ ]:


#ICU Type, In Time
#remove anywhere the patient left before start of IMV
ICU_df = merged_adt_df[merged_adt_df['location_category'] == 'icu'].copy()
ICU_df = ICU_df.drop_duplicates(subset=['encounter_block'], keep='first')
ICU_df.rename(columns={'location_type':'ICU_type','in_dttm':'icu_in_dttm'}, inplace=True)
#Merge back to block
block_df = block_df.merge(
    ICU_df[['encounter_block','ICU_type','icu_in_dttm']],
    on='encounter_block',
    how='left'
)
block_df['ICU_type'] = block_df['ICU_type'].astype(str)

#ICU Out Time
#remove anywhere the patient left before start of IMV
ICU_df = merged_adt_df[merged_adt_df['location_category'] == 'icu'].copy()
ICU_df = ICU_df.drop_duplicates(subset=['encounter_block'], keep='last')
ICU_df.rename(columns={'out_dttm':'icu_out_dttm'}, inplace=True)
ICU_df = ICU_df[['encounter_block','icu_out_dttm']]
#Merge back to block
block_df = block_df.merge(
    ICU_df,
    on='encounter_block',
    how='left'
)
block_df['ICU_type'] = block_df['ICU_type'].astype(str)

#ICU LOS, including whole encounter block, not just post IMV
ICU_df = merged_adt_df[merged_adt_df['location_category'] == 'icu'].copy()
ICU_df["icu_los_days"] = (ICU_df["out_dttm"] - ICU_df["in_dttm"]).dt.total_seconds() / (24*3600)
icu_los_df = (
    ICU_df.groupby(["encounter_block"], as_index=False)
    .agg(icu_los_days=("icu_los_days", "sum"))
)
#Merge back to blocks
block_df = block_df.merge(
    icu_los_df[["encounter_block", "icu_los_days"]],
    on=["encounter_block"],
    how="left"
)

#First ICU discharge (based on first qualifying non ICU location) or if absent last ICU out.
transfer_df = merged_adt_df[merged_adt_df['location_category'] != 'icu'].copy()
transfer_df = transfer_df.drop_duplicates(subset=['encounter_block'], keep='first')
transfer_df.rename(columns={'in_dttm':'icu_first_out_dttm'}, inplace=True)
#Merge back to blocks
block_df = block_df.merge(
    transfer_df[["encounter_block", "icu_first_out_dttm"]],
    on=["encounter_block"],
    how="left"
)
#Take last ICU out if it is greater or if the transfer out time is blank
block_df["icu_first_out_dttm"] = np.fmax(block_df["icu_first_out_dttm"], block_df['icu_out_dttm'])

del merged_adt_df, ICU_df, icu_los_df, transfer_df

_columns_of_interest = ['icu_in_dttm','ICU_type','icu_first_out_dttm','icu_out_dttm','icu_los_days']
for _col in _columns_of_interest:
    log(f"Blocks with {_col} missing: {sum(block_df[_col].isna())}")
log('ICU types:')
log(block_df['ICU_type'].value_counts())


# ### (D) Vitals
# 
# First and last vital

# In[ ]:


#Load vitals
#NOTE: clifpy resolves filters/columns per table name (filters.get(table)), so both dicts must be
#keyed by 'vitals'. Keying by the column name instead silently loads the entire vitals table.
#Only hospitalization_id and recorded_dttm are used below (first/last vital per block).
co.initialize(tables=['vitals'],
    columns={'vitals': ['hospitalization_id', 'recorded_dttm']},
    filters={'vitals': {'hospitalization_id': enc_map['hospitalization_id'].unique().tolist()}}
)
#Merge with encounter maps
co.vitals.df = co.vitals.df.merge(enc_map, on='hospitalization_id', how='left')
_vital_times_df = co.vitals.df.groupby('encounter_block').agg(
    block_first_vital_dttm = ('recorded_dttm','min'),
    block_last_vital_dttm = ('recorded_dttm','max')).reset_index()

#Merge back with the full data
block_df = block_df.merge(
    _vital_times_df,
    on=["encounter_block"],
    how="left"
)

del _vital_times_df


# ### (E) Check conformity of time columns

# In[ ]:


# Define all chronological constraints as (earlier, later, description)
# Using your naming convention from the dataframe
chronological_checks = [
    # Admission must come before everything
    ("admission_dttm", "block_first_vital_dttm",  "first vital before admission"),
    ("admission_dttm", "icu_in_dttm",              "ICU in before admission"),
    ("admission_dttm", "block_vent_start_dttm",    "vent start before admission"),
    ("admission_dttm", "block_vent_end_dttm",      "vent end before admission [SEVERE]"),
    ("admission_dttm", "icu_out_dttm",             "ICU out before admission"),
    ("admission_dttm", "death_dttm",               "death before admission [SEVERE]"),
    ("admission_dttm", "discharge_dttm",           "discharge before admission"),
    ("admission_dttm", "block_last_vital_dttm",    "last vital before admission"),

    # Middle group must come before end group
    ("block_first_vital_dttm", "block_last_vital_dttm", "last vital before first vital"),
    ("block_vent_start_dttm",  "block_vent_end_dttm",   "vent end before vent start"),
    ("icu_in_dttm",            "icu_out_dttm",          "ICU out before ICU in"),
    ("icu_in_dttm",            "icu_first_out_dttm",    "ICU first out before ICU in"),

    # Start group must come before end group (cross-bracket)
    ("block_first_vital_dttm", "discharge_dttm",        "discharge before first vital [SEVERE]"),
    ("block_first_vital_dttm", "death_dttm",            "death before first vital [SEVERE]"),
    ("icu_in_dttm",            "discharge_dttm",        "discharge before ICU in"),
    ("icu_in_dttm",            "death_dttm",            "death before ICU in"),
    ("block_vent_start_dttm",  "discharge_dttm",        "discharge before vent start"),
    ("block_vent_start_dttm",  "death_dttm",            "death before vent start"),

    # Vent end / ICU out must come before discharge
    ("block_vent_end_dttm",    "discharge_dttm",        "discharge before vent end"),
    ("icu_out_dttm",           "discharge_dttm",        "discharge before ICU out"),
    ("icu_first_out_dttm",     "icu_out_dttm",          "ICU final out before ICU first out"),
]

log(f"==== TIME COMFORMITY CHECKS (before) =====")
for earlier, later, description in chronological_checks:
    # Only compare rows where both timestamps are non-null
    mask = block_df[earlier].notna() & block_df[later].notna()
    violations = (block_df.loc[mask, later] < block_df.loc[mask, earlier]).sum()
    log(f"Encounter blocks where {description}: {violations}")

old_block_df = block_df.copy()

#Attempt to fix some of the time violations by using first vital and last vital as proxy for admission and discharge
block_df["admission_dttm"] = block_df[["admission_dttm", "block_first_vital_dttm"]].min(axis=1)
block_df["discharge_dttm"] = block_df[["discharge_dttm", "block_last_vital_dttm"]].max(axis=1)
block_df["death_dttm"] = np.where(block_df["death_dttm"].notna(), block_df[["discharge_dttm", "block_last_vital_dttm","death_dttm"]].max(axis=1), np.nan)

log(f"==== TIME COMFORMITY CHECKS (after fix) =====")
severe_violation_mask = pd.Series(False, index=block_df.index)

for earlier, later, description in chronological_checks:
    mask = block_df[earlier].notna() & block_df[later].notna()
    violations = (block_df.loc[mask, later] < block_df.loc[mask, earlier]).sum()
    log(f"Encounter blocks where {description}: {violations}")
    
    if "[SEVERE]" in description:
        severe_violation_mask |= (mask & (block_df[later] < block_df[earlier]))

log(f"Encounter blocks removed due to SEVERE time violations: {sum(severe_violation_mask)}")
strobe_excl['X_blocks_with_bad_time_data'] = sum(severe_violation_mask)
strobe_excl['F_final_blocks_with_good_data'] = block_df['encounter_block'].nunique()
block_df = block_df[~severe_violation_mask].reset_index(drop=True)


# ### (F) PT Consult Orders

# In[ ]:


#load (loading from output since key_icu_orders is not a defined table in CLIFpy and we just created it in the prior script
_pt_df = helper.load_data("clif_folder","clif_key_icu_orders")

#Filter for PT orders only
_chart_mask = _pt_df["order_category"].isin(['pt_evaluation','pt_treat'])
_pt_df = _pt_df[_chart_mask]

#Merge
_pt_df = _pt_df.merge(enc_map, on='hospitalization_id', how='right').reset_index()
_pt_df['time_diff'] = _pt_df['order_dttm'] - _pt_df['block_vent_start_dttm']

#PT pre IMV
pt_pre_imv = _pt_df[_pt_df['time_diff'].dt.total_seconds() < 0]
pt_pre_imv = pt_pre_imv.groupby('encounter_block')['order_dttm'].agg('max').reset_index()
pt_pre_imv.rename(columns={'order_dttm':'pt_pre_imv_dttm'}, inplace=True)
block_df = pd.merge(
    block_df,
    pt_pre_imv,
    on='encounter_block',
    how='left'
)

#PT post IMV
pt_post_imv = _pt_df[_pt_df['time_diff'].dt.total_seconds() >= 0]
pt_post_imv = pt_post_imv.groupby('encounter_block')['order_dttm'].agg('min').reset_index()
pt_post_imv.rename(columns={'order_dttm':'pt_post_imv_dttm'}, inplace=True)
block_df = pd.merge(
    block_df,
    pt_post_imv,
    on='encounter_block',
    how='left'
)

'''
NOTE: Remove any encounter with a pt consult in the 24 hours preceding IMV initiation.
'''
block_df['pt_pre24_IMV'] = block_df['pt_pre_imv_dttm'].notna() & ((block_df['block_vent_start_dttm'] - block_df['pt_pre_imv_dttm'] ).dt.total_seconds() < 24*3600)
strobe_excl['X_blocks_with_pt_24h_prior'] = sum( block_df['pt_pre24_IMV'])
block_df = block_df[~block_df['pt_pre24_IMV']]

del _chart_mask, _pt_df, pt_pre_imv, pt_post_imv

print(f"Block Length: {len(block_df)}")
print(f"Unique Encounter Block: {block_df['encounter_block'].nunique()}")


# ## Save Data

# In[ ]:


#Filter out from final cohort
_eb_list = block_df['encounter_block'].unique().tolist()
enc_map = enc_map[enc_map['encounter_block'].isin(_eb_list)]

#Save encounter mappings
path = os.path.join(output_folder,'intermediate','encounter_mapping.parquet')
enc_map.to_parquet(path)
del path

#Save Block
path = os.path.join(output_folder,'intermediate','block_df_1_end.parquet')
block_df.to_parquet(path)
log(f"01_cohort: AFTER DATA COLLECTION: Block Length: {len(block_df)}, Encounter Blocks {block_df['encounter_block'].nunique()}")


# ## Flowchart

# In[ ]:


# Merge all strobe dicts
strobe_counts = {}
strobe_counts.update(strobe_ab)
strobe_counts.update(strobe_c)
strobe_counts.update(strobe_d)
strobe_counts.update(strobe_excl)

pd.DataFrame(list(strobe_counts.items()), columns=['Metric', 'Value']).to_csv(
    os.path.join(output_folder,'final','strobe_counts.csv'), index=False
)

log(strobe_counts)

_fig, _ax = plt.subplots(figsize=(10, 10))
_ax.axis('off')


_boxes = [
    {"text": f"All adult encounters after date filter\n(n = {strobe_counts['A_after_age_filter']})", "xy": (0.35, 0.9)},
    {"text": f"Linked Encounter Blocks\n(n = {strobe_counts['B_after_stitching']})", "xy": (0.35, 0.8)},
    {"text": f"Encounter Blocks in year range ({year_min} - {year_max})\n(n = {strobe_counts['B2_encounter_blocks_in_year']})", "xy": (0.35, 0.7)},
    {"text": f"Encounter blocks receiving IMV\n(n = {strobe_counts['C_imv_encounter_blocks']})", "xy": (0.35, 0.6)},
    {"text": f"Encounter blocks receiving IMV >= 4 hrs\n(n = {strobe_counts['F_blocks_with_vent_4_or_more']})", "xy": (0.35, 0.5)},
    {"text": f"Encounter blocks not on trach\n(n = {strobe_counts['F_final_blocks_without_trach_at_intubation']})", "xy": (0.35, 0.4)},
    {"text": f"Encounter blocks with valid data\n(n = {strobe_counts['F_final_blocks_with_good_data']})", "xy": (0.35, 0.3)},
    {"text": f"Encounter blocks cloned\n(n = {strobe_counts['F_final_blocks_with_good_data'] - strobe_counts['X_blocks_with_pt_24h_prior']})", "xy": (0.35, 0.2)},
]

_exclusions = [
    {"text": f"Linked hospitalizations\n(n = {strobe_counts['B_stitched_hosp_ids']})", "xy": (0.8, 0.85)},
    {"text": f"Excluded: Out of year range\n(n = {strobe_counts['B2_encounter_blocks_out_of_year']})", "xy":(0.8, 0.75)},
    {"text": f"Excluded: Encounters on vent for <4 hrs\n(n = {strobe_counts['D_blocks_with_same_vent_start_end'] + strobe_counts['F_blocks_with_vent_less_than_4']})", "xy": (0.8, 0.55)},
    {"text": f"Excluded: Encounters with Tracheostomy\n(n = {strobe_counts['F_final_blocks_with_trach_at_intubation']})", "xy": (0.8, 0.45)},
    {"text": f"Excluded: Encounters with bad date-time data (n = {strobe_counts['X_blocks_with_bad_time_data']})", "xy": (0.8, 0.35)},
    {"text": f"Excluded: Encounters with PT consult 24 hours\nprior to IMV (n = {strobe_counts['X_blocks_with_pt_24h_prior']})", "xy": (0.8, 0.25)},
]

# Draw main boxes and arrows
for _i, _box in enumerate(_boxes):
    _x, _y = _box["xy"]
    _ax.add_patch(Rectangle((_x - 0.25, _y - 0.05), 0.5, 0.1, edgecolor='black', facecolor='white'))
    _ax.text(_x, _y, _box["text"], ha='center', va='center', fontsize=10)
    if _i < len(_boxes) - 1:
        _ax.add_patch(FancyArrowPatch((_x, _y - 0.05), (_x, _y - 0.1), arrowstyle='->', mutation_scale=15))

# Draw exclusion boxes and connectors
for _excl in _exclusions:
    _x, _y = _excl["xy"]
    _ax.add_patch(Rectangle((_x - 0.20, _y - 0.04), 0.38, 0.08, edgecolor='black', facecolor='#f8d7da'))
    _ax.text(_x, _y, _excl["text"], ha='center', va='center', fontsize=9)

plt.tight_layout()
path = os.path.join(output_folder,'final','graphs',f'strobe_diagram_{config["site_name"]}.png')
plt.savefig(path)
plt.close(_fig)
log("Created STROBE diagram")

