# Buddy Test Report

| | |
|---|---|
| **Buddy site / institution** | *UCMC (tested against MIMIC-CLIF)* |
| **Tester** | *Kaveri Chhikara* |
| **Date** | *2026-07-15* |

## Environment

| | |
|---|---|
| **OS** | macOS 26.5.1 |
| **RAM** | 64 GB |
| **Python** | 3.12.11 (via `uv`) |
| **R** | 4.5.2 |

Python pipeline (steps 1-3, `uv`) plus an R pipeline (steps 4-5, `renv`).

## Checks

| # | Check | Result | Notes |
|:-:|-------|--------|-------|
| 1 | Environment reproduces (`uv sync` / `renv::restore`, nothing by hand) | Pass (fixed) | `uv sync` works. `renv.lock` (R 4.4.1) won't restore on R 4.5.2: `ragg` 1.3.2 fails to compile on current Apple clang and aborts all 116 packages. Regenerated `renv.lock` for R 4.5.2. |
| 2 | Configuration works from `config/README.md` alone; no hardcoding | Pass | README doesn't mention that a MIMIC site also needs the raw MIMIC-IV path (`mimic` key). |
| 3 | Required tables/fields match what the code reads (mCIDE-valid) | Pass w/ note | CLIF tables all present. Step 2 also reads a raw MIMIC table (`hosp/patients.csv.gz`) for `anchor_year`. It's gated to MIMIC sites, so still portable, just undocumented. |
| 4 | Runs end to end with no manual edits between steps | Fail | Four crashes found (missing MIMIC path, R CRAN mirror, `bal.tab` single-arm, Fine-Gray single-state). Three fixed. Fine-Gray still fails (blocking #1). |
| 5 | Outputs in `output/final/` with right naming/type, no raw dumps | Pass w/ note | Steps 1-4 write `table1.csv`, `strobe_counts.csv`, and the weight/CIF graphs correctly. Step 5 dies before writing the outcome-model tables. |
| 6 | **Data security**: no PHI, no raw data committed *(blocking)* | Fail | Committed notebooks contain patient-level outputs (blocking #3). |
| 7 | Clinical sanity: aggregates plausible for the cohort | Pass |  |
| 8 | Documentation usable: could run from the README alone | Pass (fixed) | README says `source run_pipeline.sh`; that kills the terminal (non-blocking #1). Should be `bash run_pipeline.sh`. |

## Overall verdict

**Verdict:** Fail

Two blockers remain: patient-level outputs committed in notebooks (a data-security fail is always a Fail), and the pipeline doesn't finish end to end (Fine-Gray). Everything else was fixed during the test.

## Changes made during the buddy test

| Change | File(s) |
|---|---|
| Set `clif_folder` to the local CLIF-MIMIC data and `mimic` to the raw MIMIC-IV 3.1 path (needed because `site_name` contains "mimic") | `config/config.json` |
| Fixed `bin_sort_fill` / `hourly_fill`: `self.df[col].fillna(x, inplace=True)` is a silent no-op under pandas 3.0 Copy-on-Write, which left `pt_order`/`pt_now` as `1/NaN` instead of `0/1`. Now assigns the result back. | `code/pthelperfunctions.py` (lines 173, 391) |
| Fixed a typo in the clone-N filter: `(PT_censor_N = 0) | (pt_now = 1)` should use `==`. The `=` version was a no-op. | `code/5_ccw.R` (line 164) |
| Regenerated `renv.lock` for R 4.5.2 from current CRAN binaries. Fixes the `ragg` compile failure and the `--vanilla` CRAN-mirror error in steps 4-5. | `renv.lock` |
| Updated `pyproject.toml` to the correct Python version | `pyproject.toml` |

### Blocking issues (fix before distribution)

1. **[OPEN, needs a research decision] Step 5 Fine-Gray fails: "survival time has only a single state."**
   `5_ccw.R:566` calls `finegray(Surv(imv_to_discharge_days, dc_type) ~ ., etype = "dead")`, but `dc_type`
   (`5_ccw.R:146`) has only two levels (`alive`/`dead`). Fine-Gray is a competing-risks model and needs
   censoring plus at least two event states. Confirmed: a 2-level status throws this error, a 3-level one
   works. Since the time variable is time-to-discharge, the competing events should split death from
   discharge-alive (both are currently lumped into `alive` and treated as censoring). Not fixed here because
   it defines the causal estimand and needs author sign-off. A likely fix is coding `dc_type` as
   `censor` / `dead` / `discharge`. Everything upstream now runs, including `bal.tab` and the IPCW weighting.

2. **[FIXED] Pipeline crashed on config/dependency errors before finishing.** In order:
   - **MIMIC raw path.** Step 2 (`2_data_gathering.py:933`) reads raw MIMIC `hosp/patients.csv.gz` for
     `anchor_year` (to recover real admission years from MIMIC's shifted dates). With `mimic` unset it built
     the relative path `hosp/patients.csv.gz` and hit `FileNotFoundError`. Fixed by setting `mimic`. Only
     MIMIC sites hit this branch.
   - **R CRAN mirror.** Step 4 failed with "trying to use CRAN without setting a mirror" because
     `renv::restore` never succeeded (check 1), so the scripts' `install.packages()` fallback ran under
     `Rscript --vanilla` with no repo. Fixed by regenerating the renv library for R 4.5.2.
   - **`bal.tab` single-arm.** The pandas-3.0 `bin_sort_fill` no-op left `PT_censor_N` all `1`, so the
     clone-N exit filter (`5_ccw.R:330`, `filter(PT_censor_N == 0)`) dropped every N row. The N arm
     disappeared and `bal.tab` errored with "treatment must have at least two unique values." Fixed at the
     source in `pthelperfunctions.py`. Both arms are now present (N about 21,499, E about 7,503).

3. **[OPEN] Patient-level outputs committed to GitHub in the notebooks.** The `.ipynb` files carry saved
   patient-level cell outputs. Clear them and remove from history before distribution (MIMIC is
   credentialed-access). Either delete the notebooks (the `.py` files are what the pipeline runs) or run
   `nbstripout` and scrub git history.

### Non-blocking notes

1. **Run with `bash run_pipeline.sh`, not `source`.** The script sets `set -euo pipefail`, so under `source`
   any step failure takes down the interactive shell and the terminal restarts. Fix the README instruction.

2. **Watch for more pandas 3.0 breakage.** The re-pull moved to `pandas>=3.0.3`, where Copy-on-Write turns
   `df[col].fillna(..., inplace=True)` into a silent no-op (the cause of the `bal.tab` failure). Grep for
   other chained `inplace=True` writes, and either pin `pandas<3` with a committed `uv.lock` or audit for it.

3. **Document the MIMIC raw-data dependency in `config/README.md`.** A site whose `site_name` contains
   "mimic" must set `mimic` to the raw MIMIC-IV root (the folder with `hosp/`). Other sites can leave it
   blank. Right now you only find this out by hitting the crash.

4. **`renv.lock` isn't portable to newer R.** It pinned R 4.4.1 and versions that no longer compile on R
   4.5.x macOS (`ragg`). Regenerated here for 4.5.2. Decide whether to ship a lock that tracks current R or
   document the exact R version required.
