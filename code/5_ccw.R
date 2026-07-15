# =============================================================================
# 5_ccw.R
# Clone-Censor Weighting (CCW) implementation
#
#   - Stabilized versus unstabilized weights. (optional)
#   - For Clone E only, apply the "pt_now" (analogous to Webster et al. "recentstart")
#     logic when assigning interval weights (also optional):
#       * pt_now == 1  & uncensored  →  weight = 1 / P(uncensored)
#       * pt_now == 0  & uncensored  →  weight = 1   (no upweighting needed)
#       * censored                   →  weight = 0
#   - For Clone N: interval weight = 1 / P(uncensored) for all uncensored rows
#     (no pt_now logic; the paper's recentstart only applies to the treatment arm)
#   - Cumulative weight = product of all interval weights up to and including
#     the current time bin.
#   - Outcome analysis as follows:
#       * VFD - Zero Inflated Negative Binomials
#       * ICU LOS - Poisson
#       * Mortality (hosp, 30d, 1y) - Logistic
#       * Competing hosp mortality versus discharge alive - FineGrey
#   - Bootstrapping with non-parametric re-sampling with replacement.
# NEXT STEPS:
# Add ICU FG model, Add age subgroup analysis
# =============================================================================

# ---- Packages ----------------------------------------------------------------
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("renv")
library(this.path)
setwd(dirname(this.path()))
work_dir      <- normalizePath("..")
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::load(project = work_dir)

packages <- c("tidyverse", "pscl", "ggplot2", "dplyr", "openxlsx",
              "tibble","cobalt","glue","data.table","survival","scales")

installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])

library(tidyverse); library(pscl); library(ggplot2); library(dplyr); library(glue)
library(openxlsx); library(tibble); library(cobalt); library(this.path); library(data.table)
library(survival); library(scales)

# ---- Paths -------------------------------------------------------------------
setwd(dirname(this.path()))
work_dir      <- normalizePath("..")
output_folder <- file.path(work_dir, "output")

#----- Options -----------------------------------------------------------------
resample_N <- 10 #Effective bootstrapping resamples.
input_file_path <- file.path(output_folder, "intermediate",
                             "block_and_time_bins_for_stats.csv")
use_recent_start_logic <- FALSE
use_stabilized_weights <- FALSE
label <- "ALL" #Additional label to append at the end of output files.

# ---- Standardized output naming ------------------------------------------
# All graph and table filenames follow:
#   (name)_(trimmed|original)_(simple|MV)_(label).(ext)
# For outputs that aren't tied to a simple-vs-MV regression (e.g. diagnostic
# plots), pass model_type = NULL to drop that slot entirely.
make_filename <- function(name, trim_status = "NA", model_type = "NA", ext = "pdf") {
  parts <- c(name, trim_status, model_type, label)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  paste0(paste(parts, collapse = "_"), ".", ext)
}

# =============================================================================
# 1.  LOAD & PREP DATA
# =============================================================================
data <- read.csv(input_file_path)

# ---- Factorise -------------------------------
fac_vars <- c("sex_category", "race_category", "ethnicity_category",
              "language_category", "ICU_type")
data[fac_vars] <- lapply(data[fac_vars], function(x) {
  x[x == ""] <- NA
  factor(x)
})

# ---- Binaries --------------------------------
binary_vars <- c("pressor_flag", "paralytics_flag",
                 "is_dead", "is_dead_hosp", "is_dead_2", "is_dead_30", "is_dead_365",
                 "pt_order", "pt_now")
to_binary <- function(x, col_name) {
  if (is.logical(x)) {
    out <- as.integer(x)
  } else if (is.numeric(x)) {
    out <- x
  } else {
    x_chr <- trimws(as.character(x))
    x_chr[x_chr == ""] <- NA
    lookup <- c(
      "1" = 1, "0" = 0,
      "TRUE" = 1, "FALSE" = 0,
      "True" = 1, "False" = 0,
      "true" = 1, "false" = 0,
      "T" = 1, "F" = 0,
      "Y" = 1, "N" = 0,
      "YES" = 1, "NO" = 0,
      "Yes" = 1, "No" = 0
    )
    out <- unname(lookup[x_chr])
  }
}
for (col in binary_vars) {
  data[[col]] <- to_binary(data[[col]], col)
}

# ---- Integers -------------------------------
data$vent_free_days <- as.integer(data$vent_free_days)
data$icu_los_days   <- as.integer(data$icu_los_days)

# ---- Covariate lists & Ordering ----------------------------------------------
# Fixed (baseline) covariates — same as before
base_vars <- c("age", "sex_category", "race_category", "ethnicity_category",
               "weight_kg", "language_category","elixhauser","ICU_type")

# Time-varying covariates measured at each time_bin
tv_vars <- c("heart_rate_mean", "map_mean", "fio2_set_mean", "peep_set_mean",
             "pressor_flag", "paralytics_flag")

# Outcome Variables
out_vars <- c("vent_free_days","icu_los_days","is_dead_hosp","is_dead_30","is_dead_365",'imv_to_discharge_days')

#All covars, used to determine complete case analysis.
all_covars <- c(base_vars, tv_vars)

# Combined covariates & columns needed for complete analysis of the CCW.
all_vars <- c("encounter_block","time_bin","bin_start","bin_end","pt_order",
                "pt_now","pt_post48_IMV",base_vars, tv_vars, out_vars)

bin_df <- subset(data, select = all_vars)
bin_df <- bin_df[order(bin_df$encounter_block, bin_df$time_bin), ]

#Unique time bins
time_bins <- sort(unique(bin_df$time_bin))

# ---- Censor indicators -------------------------------------------------------
# Censor N: censored from the first PT bin onward (pt_order fills forward)
bin_df$PT_censor_N <- bin_df$pt_order

# Censor E: censored only at the final bin (bin_end == 48) if PT never occurred
bin_df$PT_censor_E <- ifelse((!bin_df$pt_post48_IMV) & bin_df$bin_end == 48, 1, 0)

# ---- Discharge Type for Fine Grey --------------------------------------------
bin_df$dc_type <- factor(
  case_when(
    bin_df$is_dead_hosp == 1 ~ "dead",
    TRUE                     ~ "alive"
  ),
  levels = c("alive", "dead")   # "alive" = censored/competing; "dead" = event
)

# =============================================================================
# 2.  COMPLETE-CASE FILTER
# =============================================================================
separate_complete_frames <- function(df,
                                     vars_needed = all_covars){
  vars_needed_N <- c("PT_censor_N", vars_needed)
  df_N <- df[complete.cases(df[, vars_needed_N]), ]
  df_N$clone <- factor("N", levels = c("N", "E"))
  df_N <- df_N[order(df_N$encounter_block, df_N$time_bin), ]
  #Keep uncensored row or rows where the censoring event happens.
  df_N <- df_N %>% filter( (PT_censor_N = 0) | (pt_now = 1))

  vars_needed_E <- c("PT_censor_E", vars_needed)
  df_E <- df[complete.cases(df[, vars_needed_E]), ]
  df_E <- df_E[order(df_E$encounter_block, df_E$time_bin), ]
  df_E$clone <- factor("E", levels = c("N", "E"))
  #All censoring occurs on the last time_bin so no need to filter.
  
  return(list(clones_N = df_N,clones_E = df_E))
  
}
# =============================================================================
# 3.  PER-TIME-BIN WEIGHTING  (Webster-Clark approach)
#
#  For each unique time_bin t:
#    a) Subset rows that are "active" at time t  (the row for that bin)
#    b) Fit:  PT_censor_* ~ base_vars + tv_vars   (binomial GLM, no time term)
#    c) P(uncensored | t) = 1 - predicted probability of censoring
#    d) Assign interval weight:
#         Clone N : 1 / P(uncensored)   for uncensored rows
#                   0                    for censored rows
#         Clone E : 1 / P(uncensored)   if pt_now == 1  AND uncensored or if not using pt_now logic.
#                   1                    if pt_now == 0  AND uncensored and using pt_now logic.
#                   0                    if censored
#    e) Collect interval weights, then compute cumulative product per encounter
#
# =============================================================================
fit_interval_weights <- function(clone_df,
                                 censor_col,   # "PT_censor_N" or "PT_censor_E"
                                 pt_now_logic = FALSE,
                                 stabilize = FALSE) {
  # pt_now_logic = TRUE  → Clone E weighting (use pt_now flag)
  # pt_now_logic = FALSE → Clone N weighting (always 1/p_uncens when uncensored)

  #Probability of censoring formula
  rhs_formula <- paste(c(base_vars, tv_vars), collapse = " + ")
  form        <- as.formula(paste(censor_col, "~", rhs_formula))
  
  # Stabilized weights formula
  rhs_stab <- paste(c(base_vars), collapse = " + ")
  form_stab <- as.formula(paste(censor_col, "~", rhs_stab))
  
  # We will collect one row per (encounter_block × time_bin) with its interval
  # weight.  Using a list then rbinding is memory-efficient.
  results_list <- vector("list", length(time_bins))

  for (i in seq_along(time_bins)) {
    tb <- time_bins[i]

    # Rows belonging to this time_bin
    bin_data <- clone_df[clone_df$time_bin == tb, ]

    # Need at least one censored AND one uncensored observation to fit the GLM.
    # If the outcome is perfectly separated (e.g., everyone uncensored at this
    # bin), skip the GLM and assign weights of 1 to uncensored, 0 to censored.
    n_cens   <- sum(bin_data[[censor_col]] == 1, na.rm = TRUE)
    n_uncens <- sum(bin_data[[censor_col]] == 0, na.rm = TRUE)

    if (n_cens == 0 || n_uncens == 0) {
      # No variation in censoring at this bin: P(uncensored) is trivially 0 or 1
      # Assign weight = 1 for uncensored, 0 for censored (no model needed)
      bin_data$p_cens   <- as.numeric(bin_data[[censor_col]])
      bin_data$p_uncens <- 1 - bin_data$p_cens
      bin_data$p_stab <- 1
    } else {
      # Fit the GLM for this time_bin
      # Response: P(censored at t)  — following Webster-Clark's formulation
      fit <- tryCatch(
        glm(form, data = bin_data, family = binomial(link = "logit")),
        error = function(e) {
          message(sprintf("GLM failed for time_bin %d (%s) probability. Defaulting to raw mean.",
                          tb, censor_col))
          NULL
        }
      )

      if (is.null(fit)) {
        # Fallback: use empirical proportion of censoring as the probability
        bin_data$p_cens   <- mean(bin_data[[censor_col]], na.rm = TRUE)
        bin_data$p_uncens <- 1 - bin_data$p_cens
      } else {
        bin_data$p_cens   <- predict(fit, newdata = bin_data, type = "response")
        bin_data$p_uncens <- 1 - bin_data$p_cens
      }
      #STABILIZATION STEP
      if (stabilize) {
        # Fit the GLM for this time_bin
        # Response: P(censored at t)  — using fix covariates only
        fit <- tryCatch(
          glm(form_stab, data = bin_data, family = binomial(link = "logit")),
          error = function(e) {
            message(sprintf("GLM failed for time_bin %d (%s) stabilization numerator. Defaulting to raw mean.",
                            tb, censor_col))
            NULL
          }
        )
        
        if (is.null(fit)) {
          # Fallback: use empirical proportion of censoring as the probability
          bin_data$p_stab   <- 1 - mean(bin_data[[censor_col]], na.rm = TRUE)
        } else {
          bin_data$p_stab   <- 1 - predict(fit, newdata = bin_data, type = "response")
        }
      } else {
        bin_data$p_stab <- 1
      }
    }

    # ---- Assign interval weight ----------------------------------------------
    if (!pt_now_logic) {
      # ---- Clone N ------------------------------------------------------------
      # Uncensored rows: weight = 1 / P(uncensored)
      # Censored rows:   weight = 0
      bin_data$interval_wt <- ifelse(
        bin_data[[censor_col]] == 0,
        bin_data$p_stab / bin_data$p_uncens,
        0
      )
    } else {
      # ---- Clone E (Webster-Clark recentstart logic) --------------------------
      # pt_now == 1 & uncensored: patient *just* started PT this bin.
      #   They could only plausibly remain in the study because they happened to
      #   start — upweight them by 1/P(uncensored).
      # pt_now == 0 & uncensored: patient did not start PT this bin.
      #   They are following the expected trajectory; weight = 1.
      # censored (regardless of pt_now): weight = 0.
      bin_data$interval_wt <- case_when(
        bin_data[[censor_col]] == 1              ~  0,                         # censored → 0
        bin_data[[censor_col]] == 0 & bin_data$pt_now == 1 ~ bin_data$p_stab / bin_data$p_uncens,  # just started PT → upweight
        bin_data[[censor_col]] == 0 & bin_data$pt_now == 0 ~ 1,               # not yet started, still in study → 1
        TRUE                                     ~  NA_real_
      )
    }

    results_list[[i]] <- bin_data[, c("encounter_block", "time_bin",
                                      "p_cens", "p_uncens", "p_stab", "interval_wt")]
  }

  # Combine all bins
  interval_weights <- do.call(rbind, results_list)
  interval_weights <- interval_weights[order(interval_weights$encounter_block,
                                             interval_weights$time_bin), ]

  # ---- Cumulative product of interval weights per encounter ------------------
  # This matches Webster-Clark's final multiplicative IPCW.
  # We use ave() with cumprod, which respects the ordering within each group.
  interval_weights$IPCW <- ave(
    interval_weights$interval_wt,
    interval_weights$encounter_block,
    FUN = cumprod
  )

  # Merge weights back onto the clone data frames for downstream use/diagnostics
  clone_df <- clone_df %>%
    left_join(interval_weights %>% dplyr::select(encounter_block, time_bin,
                                          p_cens, p_uncens, interval_wt, IPCW),
              by = c("encounter_block", "time_bin"))
  
  # Note we make this sample a global variable so we can review
  # results from outside the function.
  sample_df <<- clone_df
  
  #Make final weights for all clones by encounter block
  clone_df <- clone_df %>%
    group_by(encounter_block) %>%
    slice_tail(n = 1) %>%                        # last observed bin
    filter(!!sym(censor_col) == 0) %>%           # must be uncensored at exit
    ungroup()
  
  return(clone_df)
}

# =============================================================================
# 4.  COMBINE STEPS 2 & 3
# =============================================================================
all_vars_final <- c("encounter_block","pt_post48_IMV","clone","IPCW","dc_type", base_vars, out_vars)
clone_and_weight <- function(bin_data) {
  #Create clone data frames
  clone_frames <- separate_complete_frames(bin_data)
  
  #Give them weights
  clone_frames$clones_N <- fit_interval_weights(
    clone_df     = clone_frames$clones_N,
    censor_col   = "PT_censor_N",
    pt_now_logic = FALSE,
    stabilize = use_stabilized_weights
  )
  weights_N <<- sample_df #This is for analysis later of the original data set only.
  
  clone_frames$clones_E <- fit_interval_weights(
    clone_df     = clone_frames$clones_E,
    censor_col   = "PT_censor_E",
    pt_now_logic = use_recent_start_logic,
    stabilize = use_stabilized_weights
  )
  weights_E <<- sample_df #This is for analysis later of the original data set only.
  
  #Create one large analytic cohort
  out_df <- bind_rows(clone_frames$clones_N, clone_frames$clones_E)
  out_df <- subset(out_df, select = all_vars_final)
  
  # Trim weights at 1st / 99th percentile
  w_cut <- quantile(out_df$IPCW, probs = c(0.01, 0.99), na.rm = TRUE)
  out_df <- out_df %>%
    mutate(IPCW_trim = pmin(pmax(IPCW, w_cut[[1]]), w_cut[[2]]))
  
  out_df$eb_clone <- paste(out_df$encounter_block, out_df$clone, sep="_")
  return(out_df)
}

#Call the above function for the original data set (no bootstrapping)
original_analysis_df <- clone_and_weight(bin_df)

# =============================================================================
# 5.  IPCW Summaries  (diagnostic)
# It uses the weights from the original sample. Not for bootstrapping.
# =============================================================================

#Quick weight analytics
print(paste("Missing weights: ",sum(is.na(original_analysis_df$IPCW))))
print(paste("Infinite weights: ",sum(!is.finite(original_analysis_df$IPCW))))
print(paste("Zero weights: ",sum(original_analysis_df$IPCW == 0)))
print(paste("Weight summary (untrimmed): ", summary(original_analysis_df$IPCW)))
print(paste("Weight summary (trimmed): ", summary(original_analysis_df$IPCW_trim)))

original_analysis_df <- original_analysis_df %>%
  filter(!is.na(IPCW_trim), is.finite(IPCW_trim), IPCW_trim > 0)

# =============================================================================
# 6.  IPCW PLOTS  (diagnostic)
# It uses the weights from the original sample. Not for bootstrapping.
# =============================================================================

#Trajectory over time plots
ipcw_long <- bind_rows(
  weights_N %>% mutate(clone = "N"),
  weights_E %>% mutate(clone = "E")
) %>%
  filter(!is.na(IPCW), is.finite(IPCW), IPCW > 0) %>%
  mutate(
    clone      = factor(clone, levels = c("N", "E")),
    time_bin_f = factor(time_bin, levels = sort(unique(time_bin)))
  )

ipcw_cut  <- quantile(ipcw_long$IPCW, probs = c(0.01, 0.99), na.rm = TRUE)
ipcw_long <- ipcw_long %>%
  mutate(IPCW_trim = pmin(pmax(IPCW, ipcw_cut[[1]]), ipcw_cut[[2]]))

p_ipcw_time <- ggplot(ipcw_long, aes(x = time_bin_f, y = IPCW, fill = clone)) +
  geom_boxplot(outlier.alpha = 0.25, outlier.size = 0.8,
               position = position_dodge(width = 0.75)) +
  theme_bw() +
  labs(title = "Trajectory of Unstabilized IPCW Over Time",
       x = "Time bin", y = "Unstabilized IPCW", fill = "Clone")
ggsave(file.path(output_folder, "final", "graphs",
                 make_filename("weight_trajectory", "original", NULL)),
       plot = p_ipcw_time, width = 8, height = 5)

p_ipcw_time1 <- ggplot(ipcw_long, aes(x = time_bin_f, y = IPCW_trim, fill = clone)) +
  geom_boxplot(outlier.alpha = 0.25, outlier.size = 0.8,
               position = position_dodge(width = 0.75)) +
  theme_bw() +
  labs(title = "Trajectory of Trimmed Unstabilized IPCW Over Time",
       x = "Time bin", y = "Trimmed Unstabilized IPCW", fill = "Clone")
ggsave(file.path(output_folder, "final", "graphs",
                 make_filename("weight_trajectory", "trimmed", NULL)),
       plot = p_ipcw_time1, width = 8, height = 5)

#Final Weights Plot
g <- ggplot(original_analysis_df, aes(x = IPCW, fill = clone)) +
  geom_histogram(bins = 60, alpha = 0.5, position = "identity") +
  theme_bw() +
  labs(title = "Distribution of final weights by clone",
       x = "Final weight", y = "Count", fill = "Clone Group")
ggsave(file.path(output_folder, "final", "graphs",
                 make_filename("weight_distribution", "original", NULL)),
       plot = g, width = 7, height = 5)

g1 <- ggplot(original_analysis_df, aes(x = IPCW_trim, fill = clone)) +
  geom_histogram(bins = 60, alpha = 0.5, position = "identity") +
  theme_bw() +
  labs(title = "Distribution of trimmed final weights by clone",
       x = "Final weight", y = "Count", fill = "Clone Group")
ggsave(file.path(output_folder, "final", "graphs",
                 make_filename("weight_distribution", "trimmed", NULL)),
       plot = g1, width = 7, height = 5)


# =============================================================================
# 7.  BASELINE COVARIATE BALANCE CHECK
# It uses the weights from the original sample. Not for bootstrapping.
# =============================================================================
bal_ccw <- bal.tab(x = original_analysis_df[, base_vars], treat = original_analysis_df$clone,
                   weights = original_analysis_df$IPCW_trim, method = "weighting",
                   estimand = "ATE", s.d.denom = "pooled", un = TRUE)
print(bal_ccw)

p_balance <- love.plot(bal_ccw, stats = "mean.diffs", abs = TRUE,
                       thresholds = c(m = 0.1), var.order = "unadjusted",
                       stars = "raw", sample.names = c("Unweighted", "Weighted"),
                       title = "Baseline Covariate Balance Before and After IPCW")
print(p_balance)
ggsave(file.path(output_folder, "final", "graphs",
                 make_filename("balance_plot", "trimmed", NULL)),
       plot = p_balance, width = 8, height = 6)

# =============================================================================
# 8.  OUTCOME MODELS
# =============================================================================

#Create an outcomes data frame for bootstrapping
out_boot_df <- data.frame(
  iteration  = character(),
  type       = character(), #MV versus simple
  VFD_N      = numeric(),
  VFD_E      = numeric(),
  ICU_LOS_N  = numeric(),
  ICU_LOS_E  = numeric(),
  dead_hosp_N = numeric(),
  dead_hosp_E = numeric(),
  dead_30_N  = numeric(),
  dead_30_E  = numeric(),
  dead_365_N = numeric(),
  dead_365_E = numeric(),
  dead_FG_HR = numeric(),
  dead_FG_30_N = numeric(),
  dead_FG_30_E = numeric(),
  stringsAsFactors = FALSE
)
#Note that this will store both simple regression and MV regression as well
#as the original sample so the number of rows here will be larger than the
#number of bootstrap resamples.

#Function to calculate treatment effect with model + data
standardized_contrast <- function(fit, data, clone_var = "clone") {
  dE <- data; dN <- data
  dE[[clone_var]] <- factor("E", levels = levels(data[[clone_var]]))
  dN[[clone_var]] <- factor("N", levels = levels(data[[clone_var]]))
  pred_E <- predict(fit, newdata = dE, type = "response")
  pred_N <- predict(fit, newdata = dN, type = "response")
  pred <- tibble(mean_pred_E    = mean(pred_E, na.rm = TRUE),
         mean_pred_N    = mean(pred_N, na.rm = TRUE))
  
  return(pred)
}

#Function to calculate treatment effect with FG model + data
standardized_contrast_FG <- function(fit, data, time_point, clone_var = "clone") {
  dE <- data; dN <- data
  dE[[clone_var]] <- factor("E", levels = levels(data[[clone_var]]))
  dN[[clone_var]] <- factor("N", levels = levels(data[[clone_var]]))
  
  # Baseline cumulative hazard, computed once (independent of n subjects),
  # evaluated at time_point via a simple step-function lookup.
  H0 <- basehaz(fit, centered = FALSE)
  idx <- findInterval(time_point, H0$time)
  H0_t <- if (idx == 0) 0 else H0$hazard[idx]
  S0_t <- exp(-H0_t)
  
  # Per-subject relative risk exp(lp) at the same "zero" reference as H0,
  # so the two combine consistently. This is one vectorized calculation,
  # not a per-subject curve.
  risk_E <- predict(fit, newdata = dE, type = "risk", reference = "zero")
  risk_N <- predict(fit, newdata = dN, type = "risk", reference = "zero")
  
  # S_i(t) = 1 - S0(t) ^ exp(lp_i)  ->  P(dead at time_point) per subject
  pred_E <- 1 - S0_t ^ risk_E
  pred_N <- 1 - S0_t ^ risk_N
  
  tibble(
    frac_pred_E = mean(pred_E, na.rm = TRUE),
    frac_pred_N = mean(pred_N, na.rm = TRUE)
  )
}

model_outcomes <- function(sample_df, iteration_n, type_reg = "MV") {
  
  #Regression type sets the RHS of the formula
  if (type_reg == "simple") {
    mv_rhs <- 'clone'
  } else {
    mv_rhs <- paste(c("clone", base_vars), collapse = " + ")
  }
  
  #Models declared globally so they can be reviewed for the original sample.
  ##### VFD: ZINB #####
  fit_vfd       <<- zeroinfl(as.formula(paste("vent_free_days ~", mv_rhs, "| 1")),
                               data = sample_df, dist = "negbin", weights = IPCW_trim)
  ##### ICU LOS: Poisson #####
  fit_icu_los   <<- glm(as.formula(paste("icu_los_days ~", mv_rhs)),
                          data = sample_df, family = poisson(),   weights = IPCW_trim)
  ##### Hospital mortality: Binary #####
  fit_dead_hosp <<- glm(as.formula(paste("is_dead_hosp ~", mv_rhs)),
                          data = sample_df, family = binomial(),  weights = IPCW_trim)
  ##### 30-day mortality: Binary #####
  fit_dead_30   <<- glm(as.formula(paste("is_dead_30 ~", mv_rhs)),
                          data = sample_df, family = binomial(),  weights = IPCW_trim)
  ##### 1-year mortality: Binary #####
  fit_dead_365  <<- glm(as.formula(paste("is_dead_365 ~", mv_rhs)),
                          data = sample_df, family = binomial(),  weights = IPCW_trim)
  #### Hospital mortality: Fine-Grey (against discharge alive) ###
  data_dead_fg <<- finegray(
    Surv(imv_to_discharge_days, dc_type) ~.,
    data = sample_df,
    etype = "dead",
    id = sample_df$eb_clone,
    weights =  sample_df$IPCW_trim
  )
  fit_dead_fg <<- coxph(
    as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", mv_rhs)),
    data = data_dead_fg,
    id = eb_clone,
    weights = fgwt, #fgwt already incorporated IPCW based on line above.
    robust = TRUE
  )
  
  vfd_con       <- standardized_contrast(fit_vfd, sample_df)
  icu_con       <- standardized_contrast(fit_icu_los, sample_df)
  dead_hosp_con <- standardized_contrast(fit_dead_hosp, sample_df)
  dead_30_con   <- standardized_contrast(fit_dead_30, sample_df)
  dead_365_con  <- standardized_contrast(fit_dead_365, sample_df)
  dead_FG_30_con  <- standardized_contrast_FG(fit_dead_fg, data_dead_fg, 30)
  
  data.frame(
    iteration   = iteration_n,
    type        = type_reg,
    VFD_N       = vfd_con$mean_pred_N,
    VFD_E       = vfd_con$mean_pred_E,
    ICU_LOS_N   = icu_con$mean_pred_N,
    ICU_LOS_E   = icu_con$mean_pred_E,
    dead_hosp_N = dead_hosp_con$mean_pred_N,
    dead_hosp_E = dead_hosp_con$mean_pred_E,
    dead_30_N   = dead_30_con$mean_pred_N,
    dead_30_E   = dead_30_con$mean_pred_E,
    dead_365_N  = dead_365_con$mean_pred_N,
    dead_365_E  = dead_365_con$mean_pred_E,
    dead_FG_HR = exp(coef(fit_dead_fg)["cloneE"]),
    dead_FG_30_N = dead_FG_30_con$frac_pred_N,
    dead_FG_30_E = dead_FG_30_con$frac_pred_E
  )
}

#Create original models - Simple regressions
out_boot_df <- rbind(out_boot_df, model_outcomes(original_analysis_df,"original",type_reg = "simple"))
#Create original models - MV regressions
out_boot_df <- rbind(out_boot_df, model_outcomes(original_analysis_df,"original",type_reg = "MV"))

# =============================================================================
# 9. Marginal Fine-Gray mortality curve (0 to time_max) for clone E vs clone N
# =============================================================================

get_marginal_mortality_curve <- function(fit, data, times, clone_var = "clone") {
  dE <- data; dN <- data
  dE[[clone_var]] <- factor("E", levels = levels(data[[clone_var]]))
  dN[[clone_var]] <- factor("N", levels = levels(data[[clone_var]]))
  
  H0 <- basehaz(fit, centered = FALSE)
  H0_vec <- sapply(times, function(t) {
    idx <- findInterval(t, H0$time)
    if (idx == 0) 0 else H0$hazard[idx]
  })
  
  risk_E <- predict(fit, newdata = dE, type = "risk", reference = "zero")
  risk_N <- predict(fit, newdata = dN, type = "risk", reference = "zero")
  
  S_E <- outer(H0_vec, risk_E, function(h, r) exp(-h)^r)
  S_N <- outer(H0_vec, risk_N, function(h, r) exp(-h)^r)
  
  tibble(
    time   = times,
    mort_E = 1 - rowMeans(S_E, na.rm = TRUE),
    mort_N = 1 - rowMeans(S_N, na.rm = TRUE)
  )
}

time_grid <- 0:30

# ---- Point estimate curve from the ORIGINAL (non-bootstrapped) model ----
point_curve <- get_marginal_mortality_curve(fit_dead_fg, data_dead_fg, time_grid)

# =============================================================================
# 10.  RESULTS ORGANISATION
# Save the regression models for the original data set
# =============================================================================
extract_glm_table <- function(fit, model_name) {
  sm  <- summary(fit)$coefficients
  out <- as.data.frame(sm)
  out$term <- rownames(out)
  rownames(out) <- NULL
  p_col <- grep("Pr\\(", names(out), value = TRUE)
  out %>%
    transmute(model = model_name, component = "main", term = term,
              estimate = Estimate, se = `Std. Error`, p_value = .data[[p_col]])
}

extract_zeroinfl_table <- function(fit, model_name) {
  sm        <- summary(fit)
  count_tab <- as.data.frame(sm$coefficients$count)
  count_tab$term <- rownames(count_tab); rownames(count_tab) <- NULL
  zero_tab  <- as.data.frame(sm$coefficients$zero)
  zero_tab$term  <- rownames(zero_tab);  rownames(zero_tab)  <- NULL
  out_count <- count_tab %>%
    transmute(model = model_name, component = "count", term = term,
              estimate = Estimate, se = `Std. Error`, p_value = `Pr(>|z|)`)
  out_zero  <- zero_tab %>%
    transmute(model = model_name, component = "zero", term = term,
              estimate = Estimate, se = `Std. Error`, p_value = `Pr(>|z|)`)
  bind_rows(out_count, out_zero)
}

extract_finegray_table <- function(fit, model_name) {
  sm  <- summary(fit)$coefficients
  out <- as.data.frame(sm)
  out$term <- rownames(out)
  rownames(out) <- NULL
  p_col <- grep("Pr\\(", names(out), value = TRUE)
  out %>%
    transmute(model = model_name, component = "main", term = term,
              estimate = coef, hr = `exp(coef)`, se = `se(coef)`,
              p_value = .data[[p_col]])
}

#Create tabs
tab_vfd      <- extract_zeroinfl_table(fit_vfd,      "vent_free_days")
tab_icu_los  <- extract_glm_table(fit_icu_los,       "icu_los_days")
tab_dead_hosp<- extract_glm_table(fit_dead_hosp,     "is_dead_hosp")
tab_dead_30  <- extract_glm_table(fit_dead_30,       "is_dead_30")
tab_dead_365 <- extract_glm_table(fit_dead_365,      "is_dead_365")
tab_dead_fg <- extract_finegray_table(fit_dead_fg, "discharge_dead_alive")
#Save to excel file
wb <- createWorkbook()
addWorksheet(wb, "VFD");             writeData(wb, "VFD",              tab_vfd)
addWorksheet(wb, "ICU_LOS");         writeData(wb, "ICU_LOS",          tab_icu_los)
addWorksheet(wb, "Hosp");            writeData(wb, "Hosp",             tab_dead_hosp)
addWorksheet(wb, "30Day");           writeData(wb, "30Day",            tab_dead_30)
addWorksheet(wb, "1Year");           writeData(wb, "1Year",            tab_dead_365)
addWorksheet(wb, "FG");           writeData(wb, "FG",            tab_dead_fg)
addWorksheet(wb, "Predicted_Contrast"); writeData(wb, "Predicted_Contrast", out_boot_df[1,])
saveWorkbook(wb, file = file.path(output_folder, "final",
                                   make_filename("regression_results", "trimmed", "simple", "xlsx")),
             overwrite = TRUE)

#Create tabs
tab_vfd      <- extract_zeroinfl_table(fit_vfd,      "vent_free_days")
tab_icu_los  <- extract_glm_table(fit_icu_los,       "icu_los_days")
tab_dead_hosp<- extract_glm_table(fit_dead_hosp,     "is_dead_hosp")
tab_dead_30  <- extract_glm_table(fit_dead_30,       "is_dead_30")
tab_dead_365 <- extract_glm_table(fit_dead_365,      "is_dead_365")
tab_dead_fg <- extract_finegray_table(fit_dead_fg, "discharge_dead_alive")
#Save to excel file
wb <- createWorkbook()
addWorksheet(wb, "VFD");             writeData(wb, "VFD",              tab_vfd)
addWorksheet(wb, "ICU_LOS");         writeData(wb, "ICU_LOS",          tab_icu_los)
addWorksheet(wb, "Hosp");            writeData(wb, "Hosp",             tab_dead_hosp)
addWorksheet(wb, "30Day");           writeData(wb, "30Day",            tab_dead_30)
addWorksheet(wb, "1Year");           writeData(wb, "1Year",            tab_dead_365)
addWorksheet(wb, "FG");           writeData(wb, "FG",            tab_dead_fg)
addWorksheet(wb, "Predicted_Contrast"); writeData(wb, "Predicted_Contrast", out_boot_df[2,])
saveWorkbook(wb, file = file.path(output_folder, "final",
                                   make_filename("regression_results", "trimmed", "MV", "xlsx")),
             overwrite = TRUE)

# =============================================================================
# 11.  BOOTSTRAPPING
# =============================================================================
# ---- Bootstrap curves for confidence bands ----
curve_boot_E <- matrix(NA_real_, nrow = resample_N, ncol = length(time_grid))
curve_boot_N <- matrix(NA_real_, nrow = resample_N, ncol = length(time_grid))


#Simple Sampler Function
bootstrap_sample <- function(df, block_col = "encounter_block") {
  dt <- as.data.table(df)
  unique_blocks <- unique(dt[[block_col]])
  n_blocks <- length(unique_blocks)
  
  sampled_blocks <- sample(unique_blocks, size = n_blocks, replace = TRUE)
  
  lookup <- data.table(
    original = sampled_blocks,
    new_id   = seq_along(sampled_blocks)
  )
  
  bootstrap_df <- lookup[dt, on = c(original = block_col), allow.cartesian = TRUE, nomatch = 0]
  bootstrap_df[, (block_col) := new_id]
  bootstrap_df[, c("original", "new_id") := NULL]
  
  return(as.data.frame(bootstrap_df))
}

for (sample_i in 1:resample_N) {
  sample_df <- bootstrap_sample(bin_df)
  sample_df <- clone_and_weight(sample_df)
  #Create model (simple regressions)
  out_boot_df <- rbind(out_boot_df, model_outcomes(sample_df,sample_i,type_reg = "simple"))
  #Create model (simple regressions)
  out_boot_df <- rbind(out_boot_df, model_outcomes(sample_df,sample_i,type_reg = "MV"))
  print(paste("Completed resample ", sample_i))
  
  curve_b <- get_marginal_mortality_curve(fit_dead_fg, data_dead_fg, time_grid)
  curve_boot_E[sample_i, ] <- curve_b$mort_E
  curve_boot_N[sample_i, ] <- curve_b$mort_N
}

# =============================================================================
# 11.  BOOTSTRAPPING RESULTS
# =============================================================================

#Calculate differences and odd-ratios
out_boot_df$VFD_diff <- out_boot_df$VFD_E - out_boot_df$VFD_N
out_boot_df$VFD_OR <- out_boot_df$VFD_E / out_boot_df$VFD_N
out_boot_df$ICU_LOS_diff <- out_boot_df$ICU_LOS_E - out_boot_df$ICU_LOS_N
out_boot_df$ICU_LOS_OR <- out_boot_df$ICU_LOS_E / out_boot_df$ICU_LOS_N
out_boot_df$dead_hosp_diff <- out_boot_df$dead_hosp_E - out_boot_df$dead_hosp_N
out_boot_df$dead_hosp_OR <- out_boot_df$dead_hosp_E / out_boot_df$dead_hosp_N
out_boot_df$dead_30_diff <- out_boot_df$dead_30_E - out_boot_df$dead_30_N
out_boot_df$dead_30_OR <- out_boot_df$dead_30_E / out_boot_df$dead_30_N
out_boot_df$dead_365_diff <- out_boot_df$dead_365_E - out_boot_df$dead_365_N
out_boot_df$dead_365_OR <- out_boot_df$dead_365_E / out_boot_df$dead_365_N
out_boot_df$dead_FG_30_diff <- out_boot_df$dead_FG_30_E - out_boot_df$dead_FG_30_N
out_boot_df$dead_FG_30_OR <- out_boot_df$dead_FG_30_E / out_boot_df$dead_FG_30_N

#Divide into the simple and MV models
boots_simple <- filter(out_boot_df,type == "simple")
boots_MV <- filter(out_boot_df,type == "MV")

#Function to save results
save_summary_stats <- function(df, output_path) {
  # Sort columns alphabetically
  sorted_cols <- sort(names(df))
  
  # Compute summary stats for all numeric columns
  stats <- do.call(rbind, lapply(sorted_cols, function(col) {
    x <- df[[col]]
    
    if (is.numeric(x)) {
      x_clean <- x[!is.na(x)]
      data.frame(
        Column    = col,
        Min       = min(x_clean),
        P1        = quantile(x_clean, 0.01),
        P2_5      = quantile(x_clean, 0.025),
        P5        = quantile(x_clean, 0.05),
        P25       = quantile(x_clean, 0.25),
        Median    = median(x_clean),
        Mean      = mean(x_clean),
        SD        = sd(x_clean),
        P75       = quantile(x_clean, 0.75),
        P95       = quantile(x_clean, 0.95),
        P97_5     = quantile(x_clean, 0.975),
        P99       = quantile(x_clean, 0.99),
        Max       = max(x_clean),
        N         = length(x_clean),
        N_Missing = sum(is.na(x)),
        row.names = NULL
      )
    } else {
      NULL  # skip non-numeric columns
    }
  }))
  
  # Create workbook and styled sheet
  wb <- createWorkbook()
  addWorksheet(wb, "Summary Statistics")
  
  # Header style
  header_style <- createStyle(
    fontName    = "Arial",
    fontSize    = 11,
    fontColour  = "white",
    fgFill      = "#2F5496",
    halign      = "center",
    border      = "Bottom"
  )
  
  # Body style
  body_style <- createStyle(
    fontName  = "Arial",
    fontSize  = 10,
    numFmt    = "0.0000",
    halign    = "right",
    border    = "TopBottomLeftRight",
    borderColour = "#D9D9D9"
  )
  
  # Column name style
  col_style <- createStyle(
    fontName     = "Arial",
    fontSize     = 10,
    fontColour   = "#2F5496",
    halign       = "left"
  )
  
  # Write data
  writeData(wb, "Summary Statistics", stats, startRow = 1, startCol = 1, headerStyle = header_style)
  
  # Apply styles
  addStyle(wb, "Summary Statistics", body_style,
           rows = 2:(nrow(stats) + 1), cols = 2:ncol(stats), gridExpand = TRUE)
  addStyle(wb, "Summary Statistics", col_style,
           rows = 2:(nrow(stats) + 1), cols = 1)
  
  # Auto-fit column widths
  setColWidths(wb, "Summary Statistics", cols = 1:ncol(stats), widths = "auto")
  
  saveWorkbook(wb, output_path, overwrite = TRUE)
  message("Saved to: ", output_path)
}

# Save
boots_path <- file.path(output_folder, "final",
                        make_filename("bootstrap_summary", "trimmed", "simple", "xlsx"))
save_summary_stats(boots_simple, boots_path)
boots_path <- file.path(output_folder, "final",
                        make_filename("bootstrap_summary", "trimmed", "MV", "xlsx"))
save_summary_stats(boots_MV, boots_path)

# =============================================================================
# 11.  BOOTSTRAPPED SURVIVAL CURVE
# =============================================================================

# ---- Percentile CIs at each time point ----
ci_E <- apply(curve_boot_E, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci_N <- apply(curve_boot_N, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

plot_df <- tibble(
  time      = time_grid,
  mort_E    = point_curve$mort_E,
  mort_E_lo = ci_E[1, ],
  mort_E_hi = ci_E[2, ],
  mort_N    = point_curve$mort_N,
  mort_N_lo = ci_N[1, ],
  mort_N_hi = ci_N[2, ]
)

# ---- Plot ----
ggplot(plot_df, aes(x = time)) +
  geom_ribbon(aes(ymin = mort_E_lo, ymax = mort_E_hi, fill = "Clone E"), alpha = 0.2) +
  geom_ribbon(aes(ymin = mort_N_lo, ymax = mort_N_hi, fill = "Clone N"), alpha = 0.2) +
  geom_line(aes(y = mort_E, color = "Clone E"), linewidth = 1) +
  geom_line(aes(y = mort_N, color = "Clone N"), linewidth = 1) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Days since IMV", y = "Predicted cumulative mortality",
       color = "Clone", fill = "Clone",
       title = "Predicted 30-day mortality: Clone E vs Clone N") +
  theme_minimal()

ggsave(file.path(output_folder, "final", "graphs",
                 make_filename("mortality_curve", "95CI", "MV")),
       plot = p_mortality_curve, width = 8, height = 5)