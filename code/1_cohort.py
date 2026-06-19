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


# In[2]:


#output_folders
# NOTE: Uses print() not log() — logger is created AFTER this cell
# to ensure the log file points to the fresh output directory.
print("=== Output Folder Management ===")

'''
_output_old_folder = f'{output_folder}_old'
# Check if output folder exists
if os.path.exists(output_folder):
    print(f"Existing output folder found: {output_folder}")

    # If output_old already exists, remove it first
    if os.path.exists(_output_old_folder):
        print(f"Removing existing output_old folder...")
        shutil.rmtree(_output_old_folder)

    # Rename current output to output_old
    print(f"Renaming {output_folder} -> {_output_old_folder}")
    os.rename(output_folder, _output_old_folder)

    # Log what was backed up
    if os.path.exists(_output_old_folder):
        _backup_size = sum(
            os.path.getsize(os.path.join(_dirpath, _filename))
            for _dirpath, _dirnames, _filenames in os.walk(_output_old_folder)
            for _filename in _filenames
        ) / (1024 * 1024)
        print(f"Backup created: {_backup_size:.1f} MB")
'''

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


# In[3]:


#Load Clif Tables
co.initialize(
    tables=['patient', 'hospitalization', 'adt'],
)
log(f"Total Number of unique encounters in the hospitalization table: {co.hospitalization.df['hospitalization_id'].nunique()}")


# In[4]:


co.hospitalization.df.head()


# ## Cohort Identification 
# ### (A) Age Filter
# ### (B) Stitch Hospitalizations

# In[5]:


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


# ### (C) Identify Ventilator Usage
# - Filter first by any hospitalization that has a single IMV reportedd.
# - Create waterfall / hourly blocks of respiratory support data. See clifpy documentation for details of everything this entails.
# - Impute missing FiO2 values

# In[6]:


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
rs_waterfall['hospitalization_id'] = rs_waterfall['hospitalization_id'].ffill()
log(f"Number of rows in respiratory support waterfall: {rs_waterfall.shape[0]}")

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
#Save it for later use
path = os.path.join(output_folder,'intermediate','respiratory_support_waterfall.parquet')
rs_waterfall.to_parquet(path)
del path

log("Missing values in recorded_dttm:", rs_waterfall['recorded_dttm'].isna().sum())


# ### (D) Determine ventilation times (start/end) at encounter block level

# In[7]:


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


# In[8]:


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

log(f"01_cohort: ADULT and IMV for >0 hours : Block Length: {len(block_df)}, Encounter Blocks {block_df['encounter_block'].nunique()}")


# ### (E) Hourly sequence generation BLOCK level

# In[9]:


log("\n=== STEP E: Hourly sequence generation BLOCK level ===\n")

# 1) Load vitals
co.initialize(tables=['vitals'],
    filters={'hospitalization_id': co.encounter_mapping['hospitalization_id'].unique().tolist()}
)
#clifpy.utils.apply_outlier_handling(co.vitals) #Remove outliers per CLIF standards
co.vitals.df = co.vitals.df.merge(co.encounter_mapping, on='hospitalization_id', how='left')
#co.vitals.df = co.vitals.df.sort_values(['encounter_block', 'recorded_dttm'])


# Merge to get encounter_block on each vital
vitals_stitched = co.vitals.df.merge(block_vent_times, on='encounter_block', how='left')
# Group by block => find earliest & latest vital for that block
_vital_bounds_block = vitals_stitched.groupby('encounter_block', dropna=True)['recorded_dttm'].agg(['min', 'max']).reset_index()
_vital_bounds_block.columns = ['encounter_block', 'block_first_vital_dttm', 'block_last_vital_dttm']
block_df = block_df.merge(_vital_bounds_block, on='encounter_block', how='left')

# 2) Merge block_vent_times with vital_bounds_block
final_blocks = block_vent_times.merge(_vital_bounds_block, on='encounter_block', how='inner')

# 3) Check for bad blocks
_bad_block = final_blocks[final_blocks['block_last_vital_dttm'] < final_blocks['block_vent_start_dttm']]
final_blocks = final_blocks[final_blocks['block_last_vital_dttm'] >= final_blocks['block_vent_start_dttm']]
strobe_e = {}
strobe_e['E_blocks_with_vital_end_before_vent_start'] = _bad_block['encounter_block'].nunique()
if len(_bad_block) > 0:
    log("Warning: Some blocks have last vital < vent start (REMOVED):\n", strobe_e['E_blocks_with_vital_end_before_vent_start'])
else:
    log("There are no blocks that have last vital < vent start! Good job CLIF-ing")
#These are not removed since they may still be plausible.
_bad_block = final_blocks[final_blocks['block_first_vital_dttm'] > final_blocks['block_vent_start_dttm']]
strobe_e['E_blocks_with_vital_start_after_vent_start'] = _bad_block['encounter_block'].nunique()
if len(_bad_block) > 0:
    log("Warning: Some blocks have vent start < vital start(NOT REMOVED):\n", strobe_e['E_blocks_with_vital_start_after_vent_start'])
else:
    log("There are no blocks that have vent start < vital start! Good job CLIF-ing")

# 4) Generate the hourly sequence at block level
def _generate_hourly_sequence_block(_group):
    _blk = _group.name
    _start_time = _group['block_vent_start_dttm'].iloc[0]
    _end_time = _group['block_last_vital_dttm'].iloc[0]
    _hourly_timestamps = pd.date_range(start=_start_time, end=_end_time, freq='h')
    return pd.DataFrame({
        'encounter_block': _blk,
        'recorded_dttm': _hourly_timestamps
    })

with warnings.catch_warnings():
    warnings.simplefilter("ignore", DeprecationWarning)
    _hourly_seq_block = (
        final_blocks
        .groupby('encounter_block')
        .apply(_generate_hourly_sequence_block)
        .reset_index(drop=True)
    )
_hourly_seq_block = _hourly_seq_block.reset_index(drop=True)

_hourly_seq_block['recorded_date'] = _hourly_seq_block['recorded_dttm'].dt.date
_hourly_seq_block['recorded_hour'] = _hourly_seq_block['recorded_dttm'].dt.hour
_hourly_seq_block = _hourly_seq_block.drop(columns=['recorded_dttm'])
_hourly_seq_block = _hourly_seq_block.drop_duplicates(subset=['encounter_block', 'recorded_date', 'recorded_hour'])

# 6) Combine with respiratory support data
rs_waterfall['recorded_date'] = rs_waterfall['recorded_dttm'].dt.date
rs_waterfall['recorded_hour'] = rs_waterfall['recorded_dttm'].dt.hour
_hourly_vent_block = rs_waterfall.groupby(['encounter_block', 'recorded_date', 'recorded_hour']).agg(
    hourly_trach=('tracheostomy', 'max'),
    hourly_on_vent=('on_vent', 'max'),
).reset_index()

# Sanity check
_seq_blocks = set(_hourly_seq_block['encounter_block'].unique())
_vent_blocks = set(_hourly_vent_block['encounter_block'].unique())
_blocks_in_seq_not_vent = _seq_blocks - _vent_blocks
_blocks_in_vent_not_seq = _vent_blocks - _seq_blocks
log("Blocks in hourly_seq_block but not in hourly_vent_block:", len(_blocks_in_seq_not_vent))
if len(_blocks_in_seq_not_vent) > 0:
    log(sorted(list(_blocks_in_seq_not_vent)))
log("\nBlocks in hourly_vent_block but not in hourly_seq_block:", len(_blocks_in_vent_not_seq))

# Step 1: Reconstruct timestamps
_hourly_seq_block['recorded_dttm'] = pd.to_datetime(_hourly_seq_block['recorded_date']) + pd.to_timedelta(_hourly_seq_block['recorded_hour'], unit='h')
_hourly_vent_block['recorded_dttm'] = pd.to_datetime(_hourly_vent_block['recorded_date']) + pd.to_timedelta(_hourly_vent_block['recorded_hour'], unit='h')

# Step 2: Get max scaffold time per encounter
_max_times = (
    _hourly_seq_block.groupby('encounter_block')['recorded_dttm']
    .max().reset_index()
    .rename(columns={'recorded_dttm': 'max_seq_dttm'})
)

# Step 3: Identify extra vent rows beyond scaffold
_vent_plus_max = pd.merge(_hourly_vent_block, _max_times, on='encounter_block', how='left')
_extra_rows = _vent_plus_max[
    _vent_plus_max['recorded_dttm'] > _vent_plus_max['max_seq_dttm']
].copy()

# Step 4: Create gap-filler rows (O(1) dict lookup instead of O(N) scan)
_max_times_dict = dict(zip(_max_times['encounter_block'], pd.to_datetime(_max_times['max_seq_dttm'])))
_gap_rows = []
for _enc_id, _group in _extra_rows.groupby('encounter_block'):
    _max_time = _max_times_dict[_enc_id]
    _first_extra_time = _group['recorded_dttm'].min()

    if _first_extra_time <= _max_time + timedelta(hours=1):
        continue

    _gap_times = pd.date_range(
        start=_max_time + timedelta(hours=1),
        end=_first_extra_time - timedelta(hours=1),
        freq='H'
    )

    for _dt in _gap_times:
        _gap_rows.append({
            'encounter_block': _enc_id,
            'recorded_date': _dt.date(),
            'recorded_hour': _dt.hour,
            'recorded_dttm': _dt
        })

_gap_df = pd.DataFrame(_gap_rows)

# Step 5: Add all required columns to gap_df
_missing_cols = set(_hourly_vent_block.columns) - set(_gap_df.columns)
for _c in _missing_cols:
    _gap_df[_c] = np.nan
if len(_gap_df) > 0:
    _gap_df = _gap_df[_hourly_vent_block.columns]

# Step 6: Get scaffold rows with vent info via left join
_scaffold_df = pd.merge(
    _hourly_seq_block.drop(columns='recorded_dttm'),
    _hourly_vent_block.drop(columns='recorded_dttm'),
    on=['encounter_block', 'recorded_date', 'recorded_hour'],
    how='left'
)

_gap_df = _gap_df.drop(columns='recorded_dttm', errors='ignore')
_extra_rows = _extra_rows.drop(columns='recorded_dttm', errors='ignore')
_extra_rows = _extra_rows.drop(columns='max_seq_dttm', errors='ignore')

# Step 7: Combine all three
final_df_block_raw = pd.concat([_scaffold_df, _gap_df, _extra_rows], ignore_index=True)

# Step 8: Sort
final_df_block_raw = final_df_block_raw.sort_values(
    by=['encounter_block', 'recorded_date', 'recorded_hour']
).reset_index(drop=True)

# Step 9: Add time_from_vent
final_df_block_raw['time_from_vent'] = final_df_block_raw.groupby('encounter_block').cumcount()

_cols = ['encounter_block', 'recorded_date', 'recorded_hour', 'time_from_vent']
_cols += [col for col in final_df_block_raw.columns if col not in _cols]
final_df_block_raw = final_df_block_raw[_cols]

log("Final hourly block shape:", final_df_block_raw.shape)
log("Unique encounter_blocks:", final_df_block_raw['encounter_block'].nunique())


# ### (F) Exclusion Criteria

# In[10]:


# Count vent hours per block in first 72 hours
_first_72_hours = final_df_block_raw[(final_df_block_raw['time_from_vent'] >= 0) & (final_df_block_raw['time_from_vent'] < 72)].copy()
# Unbounded ffill is clinically appropriate here: IMV patients don't toggle on/off
# frequently, and gaps represent documentation holes, not extubation events.
_first_72_hours['hourly_on_vent'] = _first_72_hours.groupby('encounter_block')['hourly_on_vent'].ffill()
_first_72_hours['hourly_trach'] = _first_72_hours.groupby('encounter_block')['hourly_trach'].ffill()
_vent_hours_per_block = _first_72_hours.groupby('encounter_block')['hourly_on_vent'].sum()

# Exclude blocks with imv for less than 4 hours
_blocks_under_4 = _vent_hours_per_block[_vent_hours_per_block < 4].index
_final_df_block = final_df_block_raw[~final_df_block_raw['encounter_block'].isin(_blocks_under_4)].copy()

strobe_excl = {}
strobe_excl['F_blocks_with_vent_4_or_more'] = _final_df_block['encounter_block'].nunique()
strobe_excl['F_blocks_with_vent_less_than_4'] = len(_blocks_under_4)
log(f"Unique encounter blocks with IMV >=4 hours: {strobe_excl['F_blocks_with_vent_4_or_more']}")
log(f"Excluded {len(_blocks_under_4)} encounter blocks with <4 vent hours in first 72 hours of intubation.\n")

# Exclude blocks with trach at the time of intubation
_blocks_with_trach_at_intubation = _final_df_block[
    (_final_df_block['time_from_vent'] == 0) &
    (_final_df_block['hourly_trach'] == 1)
]['encounter_block'].unique()

log(f"Blocks with trach at intubation: {len(_blocks_with_trach_at_intubation)}")

final_df_block_clean = _final_df_block[
    ~_final_df_block['encounter_block'].isin(_blocks_with_trach_at_intubation)
].copy()
del _final_df_block, final_df_block_raw

strobe_excl['F_final_blocks_with_trach_at_intubation'] = len(_blocks_with_trach_at_intubation)
strobe_excl['F_final_blocks_without_trach_at_intubation'] = final_df_block_clean['encounter_block'].nunique()

log(f"Excluded {len(_blocks_with_trach_at_intubation)} blocks with trach at intubation")
log(f"Cohort size in hourly blocks: {strobe_excl['F_final_blocks_without_trach_at_intubation']}")


# ### Save Cohort Sample
# - Apply filter as noted above
# - Save progress so far including encounter stitching and cohort sample.

# In[11]:


#Filter out from final cohort
_eb_list = final_df_block_clean['encounter_block'].unique().tolist()
block_df = block_df[block_df['encounter_block'].isin(_eb_list)]
co.hospitalization.df = co.hospitalization.df[co.hospitalization.df['encounter_block'].isin(_eb_list)]
co.adt.df = co.adt.df[co.adt.df['encounter_block'].isin(_eb_list)]
co.patient.df = co.patient.df[co.patient.df['patient_id'].isin(co.hospitalization.df['patient_id'])]

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

# In[12]:


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

# In[13]:


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

# In[14]:


#Merge with encounter block
merged_adt_df = pd.merge (
    block_df[['encounter_block','block_vent_start_dttm']],
    co.adt.df,
    on='encounter_block',
    how = 'inner'
)
#keep ICUs only
merged_adt_df = merged_adt_df[merged_adt_df['location_category']== 'icu']

#Sort by in_dttm
merged_adt_df= merged_adt_df.sort_values(['encounter_block','in_dttm'], ascending=True)

#ICU Type, In Time and first out time.
#remove anywhere the patient left before start of IMV
ICU_df = merged_adt_df[merged_adt_df['out_dttm'] > merged_adt_df['block_vent_start_dttm']].copy()
ICU_df = ICU_df.drop_duplicates(subset=['encounter_block'], keep='first')
ICU_df.rename(columns={'location_type':'ICU_type','in_dttm':'icu_in_dttm','out_dttm':'icu_first_out_dttm'}, inplace=True)
#Merge back to block
block_df = block_df.merge(
    ICU_df[['encounter_block','ICU_type','icu_in_dttm','icu_first_out_dttm']],
    on='encounter_block',
    how='left'
)
block_df['ICU_type'] = block_df['ICU_type'].astype(str)

#ICU Out Time
#remove anywhere the patient left before start of IMV
ICU_df = merged_adt_df[merged_adt_df['out_dttm'] > merged_adt_df['block_vent_start_dttm']].copy()
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
merged_adt_df["icu_los_days"] = (merged_adt_df["out_dttm"] - merged_adt_df["in_dttm"]).dt.total_seconds() / (24*3600)
icu_los_df = (
    merged_adt_df.groupby(["encounter_block"], as_index=False)
    .agg(icu_los_days=("icu_los_days", "sum"))
)
#Merge back to blocks
block_df = block_df.merge(
    icu_los_df[["encounter_block", "icu_los_days"]],
    on=["encounter_block"],
    how="left"
)

del merged_adt_df, ICU_df, icu_los_df

_columns_of_interest = ['icu_in_dttm','ICU_type','icu_first_out_dttm','icu_out_dttm','icu_los_days']
for _col in _columns_of_interest:
    log(f"Blocks with {_col} missing: {sum(block_df[_col].isna())}")
log('ICU types:')
log(block_df['ICU_type'].value_counts())


# ### (D) Check conformity of time columns

# In[15]:


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


# ### (E) PT Consult Orders

# In[16]:


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
NOTE: We are not fully excluding these blocks from the analytics data set out of interest. We will keep them in the  block_df for now.
'''
block_df['pt_pre24_IMV'] = block_df['pt_pre_imv_dttm'].notna() & ((block_df['block_vent_start_dttm'] - block_df['pt_pre_imv_dttm'] ).dt.total_seconds() < 24*3600)
strobe_excl['X_blocks_with_pt_24h_prior'] = sum( block_df['pt_pre24_IMV'])

del _chart_mask, _pt_df, pt_pre_imv, pt_post_imv

print(f"Block Length: {len(block_df)}")
print(f"Unique Encounter Block: {block_df['encounter_block'].nunique()}")


# ## Save Data

# In[17]:


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

# In[18]:


# Merge all strobe dicts
strobe_counts = {}
strobe_counts.update(strobe_ab)
strobe_counts.update(strobe_c)
strobe_counts.update(strobe_d)
strobe_counts.update(strobe_e)
strobe_counts.update(strobe_excl)

pd.DataFrame(list(strobe_counts.items()), columns=['Metric', 'Value']).to_csv(
    os.path.join(output_folder,'final','strobe_counts.csv'), index=False
)

log(strobe_counts)

_fig, _ax = plt.subplots(figsize=(10, 10))
_ax.axis('off')


_boxes = [
    {"text": f"All adult encounters after date filter\n(n = {strobe_counts['A_after_age_filter']})", "xy": (0.5, 0.9)},
    {"text": f"Linked Encounter Blocks\n(n = {strobe_counts['B_after_stitching']})", "xy": (0.5, 0.75)},
    {"text": f"Encounter blocks receiving IMV\n(n = {strobe_counts['C_imv_encounter_blocks']})", "xy": (0.5, 0.6)},
    {"text": f"Encounter blocks receiving IMV >= 4 hrs\n(n = {strobe_counts['F_blocks_with_vent_4_or_more']})", "xy": (0.5, 0.45)},
    {"text": f"Encounter blocks not on trach\n(n = {strobe_counts['F_final_blocks_without_trach_at_intubation']})", "xy": (0.5, 0.3)},
    {"text": f"Encounter blocks with valid data\n(n = {strobe_counts['F_final_blocks_with_good_data']})", "xy": (0.5, 0.15)},
]

_exclusions = [
    {"text": f"Linked hospitalizations\n(n = {strobe_counts['B_stitched_hosp_ids']})", "xy": (0.8, 0.825)},
    {"text": f"Excluded: Encounters on vent for <4 hrs\n(n = {strobe_counts['D_blocks_with_same_vent_start_end'] + strobe_counts['E_blocks_with_vital_end_before_vent_start'] + strobe_counts['F_blocks_with_vent_less_than_4']})", "xy": (0.8, 0.525)},
    {"text": f"Excluded: Encounters with Tracheostomy\n(n = {strobe_counts['F_final_blocks_with_trach_at_intubation']})", "xy": (0.8, 0.375)},
    {"text": f"Excluded: Encounters with bad date-time data (n = {strobe_counts['X_blocks_with_bad_time_data']})", "xy": (0.8, 0.225)},
    {"text": f"Excluded: Encounters with PT consult 24 hours\nprior to IMV (n = {strobe_counts['X_blocks_with_pt_24h_prior']})", "xy": (0.8, 0.1)},
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

