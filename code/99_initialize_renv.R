# =============================================================================
# 99_initialize_renv.R
# Create (or overwrite) renv.lock containing every package in
# `packages_to_install`.
#
# Run with:  Rscript 99_initialize_renv.R
# =============================================================================

project_dir <- "~/project_pi_sj692/shared/PT_consults"
setwd(project_dir)

# ---- Packages ---------------------------------------------------------------
packages_to_install <- c(
  "tidyverse", "pscl", "ggplot2", "dplyr", "openxlsx",
  "tibble", "cobalt", "this.path", "glue", "data.table",
  "scales", "arrow", "comorbidity", "mets", "survival","timereg"
)

# ---- Repos ------------------------------------------------------------------
# Posit Package Manager, RHEL 9 binaries. MUST be set before any install.
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/rhel9/latest"))

# PPM only serves precompiled binaries if R identifies its platform in the
# User-Agent header. Without this you silently get source packages (and the
# cmake / libuv build failures that come with them).
options(HTTPUserAgent = sprintf(
  "R/%s R (%s)",
  getRversion(),
  paste(getRversion(), R.version$platform, R.version$arch, R.version$os)
))

options(Ncpus = max(1L, parallel::detectCores() - 1L))

# ---- Bootstrap renv ---------------------------------------------------------
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# bare = TRUE   -> don't scan the project and auto-install discovered deps
# restart=FALSE -> never try to restart the session mid-script
if (!file.exists("renv/activate.R")) {
  renv::init(bare = TRUE, restart = FALSE)
} else {
  source("renv/activate.R")
}

# renv::activate() can reset repos from a previous lockfile; re-assert ours.
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/rhel9/latest"))

# ---- Install ----------------------------------------------------------------
# Installed one at a time so a single bad package name doesn't abort the run.
failed <- character(0)
for (pkg in packages_to_install) {
  message("\n>>> installing: ", pkg)
  ok <- tryCatch({
    renv::install(pkg, prompt = FALSE)
    TRUE
  }, error = function(e) {
    message("!!! FAILED: ", pkg, " -- ", conditionMessage(e))
    FALSE
  })
  if (!ok) failed <- c(failed, pkg)
}

# ---- Snapshot ---------------------------------------------------------------
# type = "all"    -> record everything in the project library, not just packages
#                    referenced by existing project code (default "implicit"
#                    would silently omit most of packages_to_install).
# prompt = FALSE  -> overwrite renv.lock without asking.
renv::snapshot(type = "implicit", prompt = FALSE)

# ---- Report -----------------------------------------------------------------
locked <- names(renv::lockfile_read()$Packages)
missing <- setdiff(packages_to_install, locked)

message("\n========== SUMMARY ==========")
message("Lockfile:  ", normalizePath("renv.lock", mustWork = FALSE))
message("Recorded:  ", length(locked), " packages")
if (length(failed))  message("Failed to install: ", paste(failed, collapse = ", "))
if (length(missing)) message("Requested but NOT in lockfile: ", paste(missing, collapse = ", "))
if (!length(failed) && !length(missing)) message("All requested packages recorded successfully.")

# Non-zero exit so run_pipeline.sh can catch problems with `set -e`.
if (length(failed) || length(missing)) quit(status = 1)
