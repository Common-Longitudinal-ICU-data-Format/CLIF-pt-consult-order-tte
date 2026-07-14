# Setup R environment using renv

#Change WD
setwd("~/project_pi_sj692/shared/PT_consults")

# Install renv if not already installed:
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

# Initialize renv for the project:
renv::init()
# Install required packages:
# ---- Packages ----------------------------------------------------------------
packages_to_install <- c("tidyverse", "pscl", "ggplot2", "dplyr", "openxlsx",
              "tibble", "cobalt", "this.path", "glue","data.table",
              "survival","scales","arrow","comorbidity")

renv::install(packages_to_install)
# Save the project's package state:
renv::snapshot()