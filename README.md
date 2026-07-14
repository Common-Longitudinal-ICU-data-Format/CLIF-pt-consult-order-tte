# Early PT Consult Order Target Trial Emulation Project

## CLIF VERSION
2.1.0


## Objective

- Determine the effect of early PT consults on ICU outcomes for mechanically ventilated adults. Early is defined as within 48 hours of invasive mechanical ventilation or ICU admission (if transferred already intubated).  
- Created using [MIMIC IV Database](https://physionet.org/content/mimiciv/3.1/) converted to [CLIF Format](https://github.com/Common-Longitudinal-ICU-data-Format/CLIF-MIMIC)
- With early mobilization critaria algorithm taken from [Eligibility for Mobilization Algorithm](https://github.com/Common-Longitudinal-ICU-data-Format/CLIF-eligibility-for-mobilization/tree/main)

## Required CLIF tables and fields

1. **patient**: `patient_id`, `race_category`, `ethnicity_category`, `sex_category`, `death_dttm`
2. **hospitalization**: `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`,`admission_category`,`discharge_category`, `age_at_admission`
3. **adt**:`hospitalization_id`, `in_dttm`, `out_dttm`, `location_category`, `location_type`
4. **vitals**: `hospitalization_id`, `recorded_dttm`, `vital_category`, `vital_value`
   - `vital_category` = 'heart_rate', 'resp_rate', 'sbp', 'dbp', 'map', 'spo2', 'weight_kg'
5. **labs**: `hospitalization_id`, `lab_result_dttm`, `lab_order_dttm`, `lab_category`, `lab_value`, `lab_value_numeric`
   - `lab_category` = 'lactate', 'creatinine', 'bilirubin_total', 'po2_arterial', 'platelet_count'
6. **medication_admin_continuous**: `hospitalization_id`, `admin_dttm`, `med_name`, `med_category`, `med_dose`, `med_dose_unit`, `med_group`
   - `med_category` = 'norepinephrine', 'epinephrine', 'phenylephrine', 'vasopressin','dopamine', 'angiotensin', 'nicardipine', 'nitroprusside','clevidipine','cisatracurium','vecuronium','rocuronium','metaraminol','dobutamine'
7. **respiratory_support**: `hospitalization_id`, `recorded_dttm`, `device_category`, `mode_category`, `tracheostomy`, `fio2_set`, `lpm_set`, `resp_rate_set`, `peep_set`, `resp_rate_obs`
8. **patient_assessments**: `hospitalization_id`, `recorded_dttm`, `assessment_category`, `numerical_value`, `categorical_value`
   - `assessment_category` = 'braden_mobility', 'RASS', 'cam_total', 'gcs_total',
10. **key_icu_orders**: `hospitalization_id`,'order_dttm', 'order_category'

## Cohort Identification

- Adults (age >= 18)
- On invasive mechanical ventiulation for at least 4 hours.
- Without a tracheostomy.
- Without a PT consult order in 24 hours prior to intubation.

## Detailed Instructions for running the project

### 1. Requirements
The project requires **Python 3.11+** with `uv` installed and **R 4.x**. The Jupyter notebooks are converted to just Python so Jupyter itself is not required. Uses `UV` and `renv`, respectively, for dependencies.

### 2. Download This Repository

### 3. Update Config:

Follow the instructions in [`config/README.md`](config/README.md) to set your site name, the path to your CLIF tables, the file type, and the time zone.

### 4. Run Pipeline
Run the entire pipeline using the commands.
```
chmod +x run_pipeline.sh   # make it executable (one time only)
source run_pipeline.sh
```
These scripts install the required Python and R dependencies.
### Pipeline steps

| Step | Script | Language | Description |
|------|--------|----------|-------------|
| 1 | `1_cohort.py` | Python | File organization, Cohort identification, STROBE diagram |
| 2 | `2_data_gathering.py` | Python | Gathers and aggregates data from multiple CLIF tables, creates "time_bin" and "hourly" data sets. |
| 3 | `3_calculations.py` | Python | Mobilization analysis, outcomes definitions, Table 1, setup for Rscript |
| 4 | `4_table_one.R` | R | Table 1, Graphs, setup for CCW |
| 5 | `5_ccw.R` | R | Clone censor weight, outcomes models, bootstrapping |

## Output

We want the output saved to `output/final` and `output/logs`.

## Authors

*Giulio C. Rottaro Castejon*

*Jinping Liang*

*Haidong Lu*

*Fan Li*

*Snigdha Jain*

Yale University