# Buddy Test Report

| | |
|---|---|
| **Buddy site / institution** | *UCMC (tested against MIMIC-CLIF)* |
| **Tester** | *Kaveri Chhikara* |
| **Date** | *2026-07-27* |

## Environment

| | |
|---|---|
| **OS** | macOS 26.5.1 |
| **RAM** | 64 GB |
| **Python** | 3.12.11 (via `uv`) |
| **R** | 4.5.2 |

Python pipeline (steps 1-3, `uv`) plus an R pipeline (steps 4-5, `renv`). With the fixes below, all five steps now run green end to end on R 4.5.2.

## Checks

| # | Check | Result | Notes |
|:-:|-------|--------|-------|
| 1 | Environment reproduces (`uv sync` / `renv::restore`, nothing by hand) | Pass (fixed) | `uv sync` works. The shipped `renv.lock` (R 4.4.1) won't restore on R 4.5.2: `ragg` 1.3.2 fails to compile on current Apple clang and aborts all 116 packages. Regenerated `renv.lock` for R 4.5.2. |
| 2 | Configuration works from `config/README.md` alone; no hardcoding | Pass (fixed) | `config/README.md` was vague about the `mimic` key and used wrong key names (`timezone`/`filetype`). Updated to describe the raw MIMIC-IV path and the real keys (`time_zone`/`file_type`). |
| 3 | Required tables/fields match what the code reads (mCIDE-valid) | Pass w/ note | CLIF tables all present. Step 2 also reads a raw MIMIC table (`hosp/patients.csv.gz`) for `anchor_year`, gated to MIMIC sites, so still portable. |
| 4 | Runs end to end with no manual edits between steps | Pass (fixed) | Started with five crashes (missing MIMIC path, R CRAN mirror, `bal.tab` single-arm, Fine-Gray single-state, unassigned `p_mortality_curve`). All resolved; pipeline now completes. |
| 5 | Outputs in `output/final/` with right naming/type, no raw dumps | Pass | Produces `table1.csv`, `strobe_counts.csv`, weight/CIF/balance graphs, bootstrap summaries, and the Fine-Gray mortality curve. |
| 6 | **Data security**: no PHI, no raw data committed *(blocking)* | Pass (fixed) | Notebooks previously carried saved patient-level cell outputs. All five (`1_cohort`, `2_data_gathering`, `3_calculations`, `99_key_icu_orders`, `99_looking_for_mobilization`) now cleared to 0 output cells. Config paths are local and untracked. |
| 7 | Clinical sanity: aggregates plausible for the cohort | Pass | 32,124 IMV encounter blocks; STROBE counts, Table 1, IPCW weights, and clone sizes (N about 21,499 / E about 7,503) all plausible. |
| 8 | Documentation usable: could run from the README alone | Pass (fixed) | README said `source run_pipeline.sh`, which kills the terminal. Changed to `bash run_pipeline.sh`. |

## Overall verdict

**Verdict:** Pass with notes

The pipeline runs green end to end after the fixes below, and the data-security blocker is cleared (all notebook outputs stripped). Two durability items remain for the maintainer to address before wide distribution: the Fine-Gray survival-version pin, and the upstream `p_mortality_curve` bug we patched. Neither blocks a run.

## Changes made during the buddy test

| Change | File(s) |
|---|---|
| Set `clif_folder` to the local CLIF-MIMIC data and `mimic` to the raw MIMIC-IV 3.1 path (needed because `site_name` contains "mimic"). Local paths, not committed. | `config/config.json` |
| Pinned `survival` to 3.7-0 so Fine-Gray runs on R 4.5.2 (survival 3.8.x rejects the 2-level `dc_type` with "single state"; 3.7-0 accepts it). Verified 3.7-0 compiles on R 4.5.2. | `renv.lock` |
| Captured the mortality-curve plot into `p_mortality_curve` (it was built but never assigned, so the final `ggsave` errored `object 'p_mortality_curve' not found`). | `code/5_ccw.R` (~line 900) |
| Fixed a typo in the clone-N filter: `(PT_censor_N = 0) | (pt_now = 1)` should use `==`; the `=` version was a no-op. | `code/5_ccw.R` (line 164) |
| Regenerated `renv.lock` for R 4.5.2 (current CRAN binaries), which fixes the `ragg` compile failure and the `--vanilla` CRAN-mirror error in steps 4-5. | `renv.lock` |
| README run instruction `source` -> `bash run_pipeline.sh`, and documented the `mimic` raw-data dependency plus the correct config key names. | `README.md`, `config/README.md` |

The pandas 3.0 `.fillna(inplace=True)` no-op (which silently left `pt_order`/`pt_now` unfilled and emptied the clone-N arm, causing the `bal.tab` "single value" crash) was independently fixed upstream by the maintainer in `pthelperfunctions.py`. We reached the same fix; `main` already carries it.

### Blocking issues (all resolved during the buddy test)

1. **[FIXED] Patient-level outputs committed in the notebooks.** All five `.ipynb` files had saved cell
   outputs (16/30/16 in the pipeline notebooks plus the two `99_*` ones). Now cleared to 0 output cells.
   Recommend adding `nbstripout` (or a pre-commit hook) so outputs cannot be re-committed, and scrubbing the
   outputs from earlier git history since MIMIC is credentialed-access.

2. **[FIXED] Pipeline crashed at five points before completing.** In order:
   - **MIMIC raw path.** Step 2 (`2_data_gathering.py:933`) reads raw MIMIC `hosp/patients.csv.gz` for
     `anchor_year`. With `mimic` unset it built the relative path `hosp/patients.csv.gz` and hit
     `FileNotFoundError`. Fixed by setting `mimic`. Only MIMIC sites hit this branch.
   - **R CRAN mirror.** Step 4 failed with "trying to use CRAN without setting a mirror" because
     `renv::restore` never succeeded (check 1), so the scripts' `install.packages()` fallback ran under
     `Rscript --vanilla` with no repo. Fixed by regenerating the renv library for R 4.5.2.
   - **`bal.tab` single-arm.** The pandas-3.0 `.fillna` no-op left `PT_censor_N` all `1`, so the clone-N exit
     filter (`5_ccw.R:330`) dropped every N row and `bal.tab` errored "treatment must have at least two
     unique values." Fixed at the source in `pthelperfunctions.py`.
   - **Fine-Gray single-state.** `finegray(Surv(imv_to_discharge_days, dc_type) ~ ., etype="dead")`
     (`5_ccw.R:566`) with a 2-level `dc_type` errors on survival 3.8.x. Worked around by pinning survival
     3.7-0 (see non-blocking #1 for the durable fix).
   - **Unassigned plot.** The mortality-curve `ggplot` was built but not stored, so the final `ggsave`
     referenced a nonexistent `p_mortality_curve`. Fixed by assigning it.

### Non-blocking notes

1. **Fine-Gray is fragile on survival >= 3.8.** It only runs because we pinned survival 3.7-0. The durable
   fix is to code `dc_type` as three levels (censor / dead / discharge) so Fine-Gray works on any survival
   version. That changes the causal outcome coding, so it needs the author's sign-off.

2. **`p_mortality_curve` bug is in the committed code.** It fails for any site regardless of survival version
   (the `ggsave` references an object that was never assigned). Worth folding the one-line fix upstream.

3. **Watch for more pandas 3.0 breakage.** The move to `pandas>=3.0.3` (Copy-on-Write) is what turned
   `df[col].fillna(..., inplace=True)` into a silent no-op. Grep for other chained `inplace=True` writes, and
   either pin `pandas<3` with a committed `uv.lock` or audit for it.

4. **`renv.lock` portability.** The shipped lock pinned R 4.4.1 and versions that no longer compile on R
   4.5.x macOS (`ragg`). Regenerated here for 4.5.2. Decide whether to ship a lock that tracks current R or
   document the exact required R version.
