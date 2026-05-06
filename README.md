# PT Consult Order Target Trial Emulation
- Using [MIMIC IV Database](https://physionet.org/content/mimiciv/3.1/)
- Converted to [CLIF Format](https://github.com/Common-Longitudinal-ICU-data-Format/CLIF-MIMIC)
- With early mobilization critaria algorithm taken from [Eligibility for Mobilization Algorithm](https://github.com/Common-Longitudinal-ICU-data-Format/CLIF-eligibility-for-mobilization/tree/main)
- Should work with generalized CLIF data with the tables required below.

## Required CLIF Tables

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

## Running Pipeline

1. **Requires**: `Python` and `R`. The Jupyter notebooks are converted to just Python so Jupyter itself is not required. Uses Python environments for dependencies.

2. **Download This Repository**

3. **Update Config**:

- Open `config/config.json`
- Fill out the site_name, clif_folder path and time_zone. A MIMIC folder path will only be used if the site name string contains mimic (ie. MIMIC-CLIF)
- `time_bin_size` and `time_end` should be 12 and 48 respectively but they are included in the configuration to be customized if needed.

4. **Run Script**: Run the entire pipeline using the commands.
```
chmod +x run_pipeline.sh   # make it executable (one time only)
source run_pipeline.sh
```

## Output

We want the output saved to `output/final` and `output/logs`.

## Authors

*Snigdha Jain*

*Jinping Liang*

*Giulio C. Rottaro Castejon*

Yale University