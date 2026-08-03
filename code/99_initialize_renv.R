# Setup R environment using renv

#Change WD
setwd("~/project_pi_sj692/shared/PT_consults")

# Use Posit Package Manager binary repo (RHEL 9) so packages install as
# precompiled binaries instead of building from source. This avoids needing
# system dependencies (cmake, libuv) that are not available on this cluster
# node. This MUST be set before renv::init()/renv::install() below, since
# those steps trigger installs that use whatever repos is set at the time.
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/rhel9/latest"))

# Install renv if not already installed:
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

# Initialize renv for the project:
renv::init()

# Install required packages:
# ---- Packages ------------------------------------------------
packages_to_install <- c("tidyverse", "pscl", "ggplot2", "dplyr", "openxlsx",
                         "tibble", "cobalt", "this.path", "glue","data.table",
                         "scales","arrow","comorbidity","met")

renv::install("tidyverse", type = "binary")

renv::install(packages_to_install)
# Save the project's package state:
renv::snapshot()
