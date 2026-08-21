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
#       * Competing hosp mortality versus discharge alive - FineGray
#           (with marginal effect incidence curves)
#       * Competing ICU mortality versus ICU discharge alive - FineGray
#           (with marginal effect incidence curves)
#   - Bootstrapping with non-parametric re-sampling with replacement.
# NEXT STEPS:
# Review
# =============================================================================

# ---- Packages ----------------------------------------------------------------
library(this.path)
library(tidyverse); library(pscl); library(ggplot2); library(dplyr); library(glue)
library(openxlsx); library(tibble); library(cobalt); library(this.path); library(data.table)
library(timereg); library(scales); library(arrow); library(survival)

# ---- Paths -------------------------------------------------------------------
work_dir      <- dirname(dirname(this.path()))
setwd(work_dir)
output_folder <- file.path(work_dir, "output")
if (!interactive()) pdf(NULL) #Remove automatic plots

#----- LOGGING -----------------------------------------------------------------
sink_log <- file(file.path(output_folder,"logs", "05_ccw_log.txt"), open="wt")
sink(sink_log, split=TRUE)
sink(sink_log, type = "message")
sessionInfo()

#----- Options -----------------------------------------------------------------
resample_N <- 50 #Effective bootstrapping resamples.
run_sub_group <- FALSE
input_file_path <- file.path(output_folder, "intermediate",
                             "block_and_time_bins_for_stats.parquet")
use_recent_start_logic <- FALSE
use_stabilized_weights <- FALSE
dc_horizon <- 90 #Time horizon censoring for Fine-Gray model.
icu_los_horizon <- 30

# ---- Standardized output naming ------------------------------------------
# All graph and table filenames follow:
#   (label)_(name)_(trimmed|original)_(simple|MV).(ext)
# For outputs that aren't tied to a simple-vs-MV regression (e.g. diagnostic
# plots), pass model_type = NULL to drop that slot entirely.
make_filename <- function(name, label_in, trim_status = NULL,
                          model_type = NULL, ext = "pdf")
  {
  parts <- c(label_in, name, trim_status, model_type)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  paste0(paste(parts, collapse = "_"), ".", ext)
}

# =============================================================================
# 1.  LOAD & PREP DATA
# =============================================================================
data <- read_parquet(input_file_path)

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
out_vars <- c("vent_free_days","icu_los_days","is_dead_hosp","is_dead_30",
              "is_dead_365","imv_to_discharge_days", "icu_los_first_days",
              "is_dead_icu")

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

# ---- Fine-Gray Variables --------------------------------------------
dc_fg_censor <- as.integer(bin_df$imv_to_discharge_days > dc_horizon)
#2 for alive, 1 for dead, 0 for censored.
bin_df$dc_fg_cause <-  (2 - bin_df$is_dead_hosp) * (1 - dc_fg_censor)
bin_df$dc_fg_time <-  pmin(bin_df$imv_to_discharge_days, dc_horizon)

icu_fg_censor <- as.integer(bin_df$icu_los_first_days > icu_los_horizon)
#1 for ICU DC alive, 2 for dead, 0 for censored.
bin_df$icu_fg_cause <-  (1 + bin_df$is_dead_icu) * (1- icu_fg_censor)
bin_df$icu_fg_time <-  pmin(bin_df$icu_los_first_days, icu_los_horizon)

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
  df_N <- df_N %>% filter( (PT_censor_N == 0) | (pt_now == 1))

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
clone_and_weight <- function(bin_data, log_summary = FALSE) {
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
  
  # Trim weights at 1st / 99th percentile
  w_cut <- quantile(out_df$IPCW, probs = c(0.01, 0.99), na.rm = TRUE)
  out_df <- out_df %>%
    mutate(IPCW_trim = pmin(pmax(IPCW, w_cut[[1]]), w_cut[[2]]))
  
  out_df$eb_clone <- paste(out_df$encounter_block, out_df$clone, sep="_")
  
  #Quick weight analytics
  if (log_summary) {
    print(paste("Missing weights: ",sum(is.na(out_df$IPCW))))
    print(paste("Infinite weights: ",sum(!is.finite(out_df$IPCW))))
    print(paste("Zero weights: ",sum(out_df$IPCW == 0)))
    print(paste("Weight summary (untrimmed): ", summary(out_df$IPCW)))
    print(paste("Weight summary (trimmed): ", summary(out_df$IPCW_trim)))
  }
  
  out_df <- out_df %>%
    filter(!is.na(IPCW_trim), is.finite(IPCW_trim), IPCW_trim > 0)
  
  return(out_df)
}

# =============================================================================
# 5.  IPCW Summaries  (diagnostic)
# - Plots of trajectory over time. (trimmed and untrimmed)
# - Plots distributions. (trimmed and untrimmed)
# - Plot covariate balance.
# =============================================================================

analyze_ipcw <- function(block_all, tv_e, tv_n, label_in) {
  #block_df must contain IPCW and IPCW_trim
  #tv_e and tv_n are time varying and must have time_bin and IPCW
  #label_in is string
  
  #Trajectory over time plots
  ipcw_long <- bind_rows(
    tv_n %>% mutate(clone = "N"),
    tv_e %>% mutate(clone = "E")
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
                   make_filename("weight_trajectory", label_in,
                                 trim_status = "original")),
         plot = p_ipcw_time, width = 8, height = 5)
  
  p_ipcw_time1 <- ggplot(ipcw_long, aes(x = time_bin_f, y = IPCW_trim, fill = clone)) +
    geom_boxplot(outlier.alpha = 0.25, outlier.size = 0.8,
                 position = position_dodge(width = 0.75)) +
    theme_bw() +
    labs(title = "Trajectory of Trimmed Unstabilized IPCW Over Time",
         x = "Time bin", y = "Trimmed Unstabilized IPCW", fill = "Clone")
  ggsave(file.path(output_folder, "final", "graphs",
                   make_filename("weight_trajectory", label_in,
                                 trim_status = "trimmed")),
         plot = p_ipcw_time1, width = 8, height = 5)

  #Final Weights Plot
  g <- ggplot(block_all, aes(x = IPCW, fill = clone)) +
    geom_histogram(bins = 60, alpha = 0.5, position = "identity") +
    theme_bw() +
    labs(title = "Distribution of final weights by clone",
         x = "Final weight", y = "Count", fill = "Clone Group")
  ggsave(file.path(output_folder, "final", "graphs",
                   make_filename("weight_distribution", label_in,
                                 trim_status = "original")),
         plot = g, width = 7, height = 5)
  
  g1 <- ggplot(block_all, aes(x = IPCW_trim, fill = clone)) +
    geom_histogram(bins = 60, alpha = 0.5, position = "identity") +
    theme_bw() +
    labs(title = "Distribution of trimmed final weights by clone",
         x = "Final weight", y = "Count", fill = "Clone Group")
  ggsave(file.path(output_folder, "final", "graphs",
                   make_filename("weight_distribution", label_in,
                                 trim_status = "trimmed")),
         plot = g1, width = 7, height = 5)
  
  #Covariate balance plot
  bal_ccw <- bal.tab(x = block_all[, base_vars], treat = block_all$clone,
                     weights = block_all$IPCW_trim, method = "weighting",
                     estimand = "ATE", s.d.denom = "pooled", un = TRUE)
  
  p_balance <- love.plot(bal_ccw, stats = "mean.diffs", abs = TRUE,
                         thresholds = c(m = 0.1), var.order = "unadjusted",
                         stars = "raw", sample.names = c("Unweighted", "Weighted"),
                         title = "Baseline Covariate Balance Before and After IPCW")
  print(p_balance)
  ggsave(file.path(output_folder, "final", "graphs",
                   make_filename("balance_plot", label_in,
                                 trim_status = "trimmed")),
         plot = p_balance, width = 8, height = 6)

}
# =============================================================================
# 6.  OUTCOME MODELS
# =============================================================================

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

# -----------------------------------------------------------------------------
# Fine-Gray helpers for timereg::comp.risk(..., model = "prop")
# -----------------------------------------------------------------------------
# The fitted model is
#     P(T <= t, cause = 1 | X, Z) = 1 - exp( -exp( X'A(t) + Z'gamma ) )
#   X : terms NOT wrapped in const() -> non-parametric, time-varying cumulative
#       coefficients A(t), stored in fit$cum (first column is time)
#   Z : terms wrapped in const()     -> constant log subdistribution hazard
#       ratios, stored in fit$gamma
# With mv_rhs_fg as built above, `clone` is in X and the baseline covariates
# are in Z.

#Per-subject predicted CIF, returned as a [subjects x times] matrix -- the same
#shape predict.cifreg()$cif used to return.
fg_predict_cif <- function(fit, newdata, times) {
  p <- predict(fit, newdata = as.data.frame(newdata), times = times,
               se = FALSE, uniform = FALSE, resample.iid = 0)
  cif <- p$P1
  if (is.null(dim(cif))) cif <- matrix(cif, nrow = nrow(newdata))

  #Before the first jump time every cumulative coefficient is 0, and the
  #cloglog link maps 0 to 1 - exp(-1) = 0.63 rather than 0.  The CIF is 0 by
  #definition there, so blank those columns out -- the plotting grids start at
  #t = 0, so without this every curve would start at ~0.6 instead of 0.
  t_first <- min(fit$cum[, 1], na.rm = TRUE)
  before  <- times < t_first
  if (any(before)) cif[, before] <- 0

  cif
}

#Log subdistribution hazard ratio for one term, as a length-1 named vector.
#  const() term -> its constant coefficient from fit$gamma
#  otherwise    -> the cumulative coefficient B(t) from fit$cum, at `time`
#                  (default: the last estimated time point)
fg_log_shr <- function(fit, term = "cloneE", time = NULL) {
  if (!is.null(fit$gamma)) {
    g <- as.numeric(fit$gamma)
    names(g) <- rownames(fit$gamma)
    if (term %in% names(g)) return(setNames(g[[term]], term))
  }
  if (term %in% colnames(fit$cum)) {
    tt  <- fit$cum[, 1]
    idx <- if (is.null(time)) length(tt) else max(1L, sum(tt <= time))
    return(setNames(fit$cum[idx, term], term))
  }
  setNames(NA_real_, term)
}

#Function to calculate treatment effect with FG model + data
# Adapted for timereg::comp.risk objects.  Same arguments and the same one-row
# tibble (frac_pred_E / frac_pred_N) as the mets::cifreg version; the actual
# prediction is delegated to get_marginal_curve() so there is a single
# standardisation implementation.
standardized_contrast_FG <- function(fit, data, time_point, clone_var = "clone") {
  curve <- get_marginal_curve(fit, data, time_point, clone_var = clone_var)
  tibble(
    frac_pred_E = curve$pred_E[1],
    frac_pred_N = curve$pred_N[1]
  )
}

model_outcomes <- function(sample_df, iteration_n, type_reg = "MV", trimmed_weights = FALSE) {
  
  #Regression type sets the RHS of the formula
  if (type_reg == "simple") {
    mv_rhs <- 'clone'
    mv_rhs_fg <- 'clone'
  } else {
    mv_rhs <- paste(c("clone", base_vars), collapse = " + ")
    mv_rhs_fg <- paste(c("clone",
                         paste0(paste0("const(",base_vars),")")),
                       collapse = " + ")
  }
  
  #Trimmed versus untrimmed weights
  if (trimmed_weights) {
    sample_df$IPCW <- sample_df$IPCW_trim
    trim_text <- "trimmed"
  } else {
    trim_text <- "original"
  }
  
  #Models declared globally so they can be reviewed for the original sample.
  ##### VFD: ZINB #####
  fit_vfd       <<- zeroinfl(as.formula(paste("vent_free_days ~", mv_rhs, "| 1")),
                               data = sample_df, dist = "negbin", weights = IPCW)
  vfd_con       <- standardized_contrast(fit_vfd, sample_df)
  
  ##### ICU LOS: Poisson #####
  fit_icu_los   <<- glm(as.formula(paste("icu_los_days ~", mv_rhs)),
                          data = sample_df, family = poisson(),   weights = IPCW)
  icu_con       <- standardized_contrast(fit_icu_los, sample_df)
  
  ##### Hospital mortality: Binary #####
  fit_dead_hosp <<- glm(as.formula(paste("is_dead_hosp ~", mv_rhs)),
                          data = sample_df, family = binomial(),  weights = IPCW)
  dead_hosp_con <- standardized_contrast(fit_dead_hosp, sample_df)
  
  ##### 30-day mortality: Binary #####
  fit_dead_30   <<- glm(as.formula(paste("is_dead_30 ~", mv_rhs)),
                          data = sample_df, family = binomial(),  weights = IPCW)
  dead_30_con   <- standardized_contrast(fit_dead_30, sample_df)
  
  ##### 1-year mortality: Binary #####
  fit_dead_365  <<- glm(as.formula(paste("is_dead_365 ~", mv_rhs)),
                          data = sample_df, family = binomial(),  weights = IPCW)
  dead_365_con  <- standardized_contrast(fit_dead_365, sample_df)
  
  #### Hospital mortality: Fine-Grey (against discharge alive) ###
  #timereg::comp.risk() takes `weights` as a plain vector and needs a plain
  #data.frame, so build that frame once and take the weights from it to
  #guarantee the rows line up.
  fg_data <- as.data.frame(sample_df)

  #model = "prop" is the Fine-Gray model on the cloglog scale:
  #  P(T <= t, cause = 1 | X, Z) = 1 - exp(-exp(X'A(t) + Z'gamma))
  #Event() defaults to cens.code = 0, which matches the coding above.
  #n.sim drives the resampling behind the non-parametric tests; set it to 0
  #inside the bootstrap if the run time becomes a problem.
  fit_dead_fg <<- comp.risk(
    as.formula(paste("Event(dc_fg_time,dc_fg_cause) ~", mv_rhs_fg)),
    data    = fg_data,
    cause   = 1,       #Per FG variables definitions above.
    n.sim   = 5000,
    model   = "prop",
    weights = fg_data$IPCW
    #Censoring is handled by a pooled Kaplan-Meier.  For the previous
    #cens.model = ~strata(clone) behaviour add:
    #, cens.model = "stratKM", cens.formula = ~clone
  )
  dead_FG_30_con  <- standardized_contrast_FG(fit_dead_fg, sample_df, 30)
  
  #### ICU LOS: Fine-Grey (against death) ###
  fit_icu_fg <<- comp.risk(
    as.formula(paste("Event(icu_fg_time,icu_fg_cause) ~", mv_rhs_fg)),
    data    = fg_data,
    cause   = 1,       #Per FG variables definitions above.
    n.sim   = 5000,
    model   = "prop",
    weights = fg_data$IPCW
  )
  icu_FG_10_con  <- standardized_contrast_FG(fit_icu_fg, sample_df, 10)
  
  #Organize outcome into a single row of a data frame
  output_df <- data.frame(
    iteration   = iteration_n,
    type        = type_reg,
    trim        = trim_text,
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
    #`clone` is a non-const() term, so there is no single coefficient for it:
    #fg_log_shr() falls back to the cumulative coefficient at the end of
    #follow-up.  (It returns the const() coefficient instead if clone is ever
    #wrapped in const(), so this line works either way.)
    dead_FG_HR = exp(fg_log_shr(fit_dead_fg, "cloneE")),
    dead_FG_30_N = dead_FG_30_con$frac_pred_N,
    dead_FG_30_E = dead_FG_30_con$frac_pred_E,
    icu_FG_10_N = icu_FG_10_con$frac_pred_N,
    icu_FG_10_E = icu_FG_10_con$frac_pred_E
  )
  
  #Calculate differences and odd-ratios
  output_df$VFD_diff <- output_df$VFD_E - output_df$VFD_N
  output_df$VFD_OR <- output_df$VFD_E / output_df$VFD_N
  output_df$ICU_LOS_diff <- output_df$ICU_LOS_E - output_df$ICU_LOS_N
  output_df$ICU_LOS_OR <- output_df$ICU_LOS_E / output_df$ICU_LOS_N
  output_df$dead_hosp_diff <- output_df$dead_hosp_E - output_df$dead_hosp_N
  output_df$dead_hosp_OR <- output_df$dead_hosp_E / output_df$dead_hosp_N
  output_df$dead_30_diff <- output_df$dead_30_E - output_df$dead_30_N
  output_df$dead_30_OR <- output_df$dead_30_E / output_df$dead_30_N
  output_df$dead_365_diff <- output_df$dead_365_E - output_df$dead_365_N
  output_df$dead_365_OR <- output_df$dead_365_E / output_df$dead_365_N
  output_df$dead_FG_30_diff <- output_df$dead_FG_30_E - output_df$dead_FG_30_N
  output_df$dead_FG_30_OR <- output_df$dead_FG_30_E / output_df$dead_FG_30_N
  output_df$icu_FG_10_diff <- output_df$icu_FG_10_E - output_df$icu_FG_10_N
  output_df$icu_FG_10_OR <- output_df$icu_FG_10_E / output_df$icu_FG_10_N
  
  return(output_df)
}

# =============================================================================
# 7. CIF Curves for clone E vs clone N
#- The marginal curves are for predictive models (ie. Fine-Gray)
#- The Aalen-Johansen curves are for weighted observed data.
# =============================================================================

time_grid_dc <- 0:dc_horizon
time_grid_icu <- 0:icu_los_horizon

get_marginal_curve <- function(fit, data, times, clone_var = "clone") {

  dE <- data; dN <- data
  dE[[clone_var]] <- factor("E", levels = levels(data[[clone_var]]))
  dN[[clone_var]] <- factor("N", levels = levels(data[[clone_var]]))

  # G-computation / marginal standardisation: predict every subject's CIF
  # under each clone assignment and average over the sample.  fg_predict_cif()
  # returns a [subjects x times] matrix, exactly like predict.cifreg()$cif did.
  cif_E <- fg_predict_cif(fit, dE, times)
  cif_N <- fg_predict_cif(fit, dN, times)

  tibble(
    time   = times,
    pred_E = as.numeric(colMeans(cif_E, na.rm = TRUE)),
    pred_N = as.numeric(colMeans(cif_N, na.rm = TRUE))
  )
}

get_aj_curve <- function(treat, time_to, cause, weights, times) {
  #Input the arrays separately for treatment, time_to (event) and cause.
  #Assumes cuase = 1 is the event of interest
  #times is the array to plot over time
  data <- data.frame(treat = treat, time = time_to, cause=cause, weight=weights)
  
  #Weighted survival fit
  fit_aj <- survfit(Surv(time, cause, type="mstate") ~ treat,
                     data = data, weights = weight, 
                     conf.type="none", robust=TRUE)
  
  # extend = TRUE carries the last estimate forward past the final event
  # time instead of dropping those rows
  s <- summary(fit_aj, times = times, extend = TRUE)
  
  est   <- s$pstate[, "1"]
  strat <- sub("^treat=", "", as.character(s$strata))
  parts <- split(est, strat)
  
  if (!all(lengths(parts) == length(times))) {
    stop("strata returned uneven lengths: ",
         paste(names(parts), lengths(parts), sep = "=", collapse = ", "))
  }
  
  tibble(
    time  = times,
    pred_E = parts[["E"]],
    pred_N = parts[["N"]]
  )
}
# =============================================================================
# 8.  RESULTS ORGANISATION
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

# Adapted for timereg::comp.risk objects.  summary() on a "comprisk" object only
# prints and returns NULL, so the table is assembled from the fit components:
#   * fit$gamma   -> const() terms: constant log subdistribution hazard ratios
#   * fit$cum     -> non-parametric terms: cumulative coefficients B(t)
#   * fit$obs./pval.* -> the resampling-based non-parametric tests (need n.sim > 0)
# Columns are the same as before (model, component, term, estimate, hr, se,
# p_value) with one extra trailing column, `time`, holding the time point a
# time-varying estimate refers to (NA for everything else).
extract_finegray_table <- function(fit, model_name) {
  rows <- list()

  # ---- Parametric const() terms -------------------------------------------
  # exp(estimate) is a constant subdistribution hazard ratio.
  if (!is.null(fit$gamma)) {
    est <- as.numeric(fit$gamma)
    V   <- if (!is.null(fit$robvar.gamma)) fit$robvar.gamma else fit$var.gamma
    se  <- sqrt(diag(as.matrix(V)))
    rows$const <- data.frame(
      model     = model_name,
      component = "main",
      term      = rownames(fit$gamma),
      estimate  = est,
      hr        = exp(est),
      se        = se,
      p_value   = 2 * pnorm(-abs(est / se)),
      time      = NA_real_,
      row.names = NULL, stringsAsFactors = FALSE
    )
  }

  # ---- Non-parametric (time-varying) terms ---------------------------------
  # These have no single coefficient; report the cumulative coefficient B(t)
  # at the last estimated time point, so exp(estimate) is the ratio of
  # cumulative subdistribution hazards over the whole follow-up window.
  cum <- as.matrix(fit$cum)
  if (ncol(cum) > 1) {
    i_last <- nrow(cum)
    B_last <- as.numeric(cum[i_last, -1])
    se_last <- if (!is.null(fit$var.cum)) {
      sqrt(as.numeric(as.matrix(fit$var.cum)[i_last, -1]))
    } else rep(NA_real_, length(B_last))
    rows$tvar <- data.frame(
      model     = model_name,
      component = "time_varying_cum",
      term      = colnames(cum)[-1],
      estimate  = B_last,
      hr        = exp(B_last),
      se        = se_last,
      p_value   = 2 * pnorm(-abs(B_last / se_last)),
      time      = cum[i_last, 1],
      row.names = NULL, stringsAsFactors = FALSE
    )
  }

  # ---- Non-parametric tests ------------------------------------------------
  # Supremum / Kolmogorov-Smirnov / Cramer-von Mises statistics from the
  # n.sim resamples.  The statistic goes in `estimate`; hr and se are NA.
  add_test <- function(comp, obs, pval) {
    if (is.null(obs) || length(obs) == 0) return(NULL)
    nm <- if (!is.null(dim(obs))) rownames(obs) else names(obs)
    if (is.null(nm)) nm <- paste0("term", seq_along(as.numeric(obs)))
    data.frame(
      model     = model_name,
      component = comp,
      term      = nm,
      estimate  = as.numeric(obs),
      hr        = NA_real_,
      se        = NA_real_,
      p_value   = if (is.null(pval)) NA_real_ else as.numeric(pval),
      time      = NA_real_,
      row.names = NULL, stringsAsFactors = FALSE
    )
  }
  # H0: B(t) = 0 over the whole time span (is the effect there at all?)
  rows$sup <- add_test("np_test_significance", fit$obs.testBeq0,    fit$pval.testBeq0)
  # H0: the effect is constant in time (i.e. a plain proportional FG effect)
  rows$ks  <- add_test("np_test_constant_KS",  fit$obs.testBeqC,    fit$pval.testBeqC)
  rows$cvm <- add_test("np_test_constant_CvM", fit$obs.testBeqC.is, fit$pval.testBeqC.is)

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

sumarize_regressions <- function (contrast_rows, file_name_in) {
  #Contrast rows includes the predicted contrast.
  #Note this replies in the regression models saved currently as global label.
  
  #Create tabs
  tab_vfd      <- extract_zeroinfl_table(fit_vfd,      "vent_free_days")
  tab_icu_los  <- extract_glm_table(fit_icu_los,       "icu_los_days")
  tab_dead_hosp<- extract_glm_table(fit_dead_hosp,     "is_dead_hosp")
  tab_dead_30  <- extract_glm_table(fit_dead_30,       "is_dead_30")
  tab_dead_365 <- extract_glm_table(fit_dead_365,      "is_dead_365")
  tab_dead_fg <- extract_finegray_table(fit_dead_fg, "discharge_dead_alive")
  tab_icu_fg <- extract_finegray_table(fit_icu_fg, "icu_dc_dead")
  
  #Save to excel file
  wb <- createWorkbook()
  addWorksheet(wb, "VFD");             writeData(wb, "VFD",              tab_vfd)
  addWorksheet(wb, "ICU_LOS");         writeData(wb, "ICU_LOS",          tab_icu_los)
  addWorksheet(wb, "Hosp");            writeData(wb, "Hosp",             tab_dead_hosp)
  addWorksheet(wb, "30Day");           writeData(wb, "30Day",            tab_dead_30)
  addWorksheet(wb, "1Year");           writeData(wb, "1Year",            tab_dead_365)
  addWorksheet(wb, "Hosp_MORT_FG");           writeData(wb, "Hosp_MORT_FG",            tab_dead_fg)
  addWorksheet(wb, "ICU_LOS_FG");           writeData(wb, "ICU_LOS_FG",            tab_icu_fg)
  addWorksheet(wb, "Predicted_Contrast"); writeData(wb, "Predicted_Contrast", contrast_rows)
  saveWorkbook(wb, file = file.path(output_folder, "final",file_name_in),
               overwrite = TRUE)
}
# =============================================================================
# 9.  BOOTSTRAPPING
# =============================================================================

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

bootstrap_all <- function(boot_in_df) {
  # ---- Bootstrap curves for confidence bands ----
  #FG curves
  curve_boot_dc_original <- array(NA_real_, dim = c(resample_N, length(time_grid_dc), 2),
                      dimnames = list(NULL, NULL, c("E", "N")))
  curve_boot_dc_trimmed <- array(NA_real_, dim = c(resample_N, length(time_grid_dc), 2),
                                  dimnames = list(NULL, NULL, c("E", "N")))
  
  curve_boot_icu_original <- array(NA_real_, dim = c(resample_N, length(time_grid_icu), 2),
                                  dimnames = list(NULL, NULL, c("E", "N")))
  curve_boot_icu_trimmed <- array(NA_real_, dim = c(resample_N, length(time_grid_icu), 2),
                                 dimnames = list(NULL, NULL, c("E", "N")))
  #AJ curves
  curve_boot_aj_dc_original <- array(NA_real_, dim = c(resample_N, length(time_grid_dc), 2),
                                  dimnames = list(NULL, NULL, c("E", "N")))
  curve_boot_aj_dc_trimmed <- array(NA_real_, dim = c(resample_N, length(time_grid_dc), 2),
                                 dimnames = list(NULL, NULL, c("E", "N")))
  
  curve_boot_aj_icu_original <- array(NA_real_, dim = c(resample_N, length(time_grid_icu), 2),
                                   dimnames = list(NULL, NULL, c("E", "N")))
  curve_boot_aj_icu_trimmed <- array(NA_real_, dim = c(resample_N, length(time_grid_icu), 2),
                                  dimnames = list(NULL, NULL, c("E", "N")))
  
  
  #Create an outcomes data frame for bootstrapping
  out_boot_df <- data.frame()

  for (sample_i in 1:resample_N) {
    sample_df <- bootstrap_sample(boot_in_df)
    sample_df <- clone_and_weight(sample_df)
    #Run 4x regression models (original/trimmed) and (MV/simple)
    ####### original - simple
    out_boot_df <- rbind(out_boot_df,
                         model_outcomes(sample_df,sample_i,
                                        type_reg = "simple",
                                        trimmed_weights = FALSE))
    ####### trimmed - simple
    out_boot_df <- rbind(out_boot_df,
                         model_outcomes(sample_df,sample_i,
                                        type_reg = "simple",
                                        trimmed_weights = TRUE))
    

    ####### original - MV (with curves)
    out_boot_df <- rbind(out_boot_df,
                         model_outcomes(sample_df,sample_i,
                                        type_reg = "MV",
                                        trimmed_weights = FALSE))
    
    curve_b <- get_marginal_curve(fit_dead_fg, sample_df, time_grid_dc)
    curve_boot_dc_original[sample_i,,"E"] <- curve_b$pred_E
    curve_boot_dc_original[sample_i,,"N"] <- curve_b$pred_N
    
    curve_b <- get_marginal_curve(fit_icu_fg, sample_df, time_grid_icu)
    curve_boot_icu_original[sample_i,,"E"] <- curve_b$pred_E
    curve_boot_icu_original[sample_i,,"N"] <- curve_b$pred_N
    
    ####### trimmed - MV (with curves)
    out_boot_df <- rbind(out_boot_df,
                         model_outcomes(sample_df,sample_i,
                                        type_reg = "MV",
                                        trimmed_weights = TRUE))
    
    curve_b <- get_marginal_curve(fit_dead_fg, sample_df, time_grid_dc)
    curve_boot_dc_trimmed[sample_i,,"E"] <- curve_b$pred_E
    curve_boot_dc_trimmed[sample_i,,"N"] <- curve_b$pred_N
    
    curve_b <- get_marginal_curve(fit_icu_fg, sample_df, time_grid_icu)
    curve_boot_icu_trimmed[sample_i,,"E"] <- curve_b$pred_E
    curve_boot_icu_trimmed[sample_i,,"N"] <- curve_b$pred_N
    
    #Plot AJ curves
    curve_b <- get_aj_curve(sample_df$clone,
                            sample_df$dc_fg_time,
                            sample_df$dc_fg_cause,
                            sample_df$IPCW,
                            time_grid_dc)
    curve_boot_aj_dc_original[sample_i,,"E"] <- curve_b$pred_E
    curve_boot_aj_dc_original[sample_i,,"N"] <- curve_b$pred_N
    
    curve_b <- get_aj_curve(sample_df$clone,
                            sample_df$icu_fg_time,
                            sample_df$icu_fg_cause,
                            sample_df$IPCW,
                            time_grid_icu)
    curve_boot_aj_icu_original[sample_i,,"E"] <- curve_b$pred_E
    curve_boot_aj_icu_original[sample_i,,"N"] <- curve_b$pred_N
    
    curve_b <- curve_b <- get_aj_curve(sample_df$clone,
                                       sample_df$dc_fg_time,
                                       sample_df$dc_fg_cause,
                                       sample_df$IPCW_trim,
                                       time_grid_dc)
    curve_boot_dc_trimmed[sample_i,,"E"] <- curve_b$pred_E
    curve_boot_dc_trimmed[sample_i,,"N"] <- curve_b$pred_N
    
    curve_b <- get_aj_curve(sample_df$clone,
                            sample_df$icu_fg_time,
                            sample_df$icu_fg_cause,
                            sample_df$IPCW_trim,
                            time_grid_icu)
    curve_boot_aj_icu_trimmed[sample_i,,"E"] <- curve_b$pred_E
    curve_boot_aj_icu_trimmed[sample_i,,"N"] <- curve_b$pred_N
    
    print(paste("Completed resample", sample_i))
  }
  
  return(list(outcome_df = out_boot_df,
              curve_dc_original = curve_boot_dc_original,
              curve_dc_trimmed = curve_boot_dc_trimmed,
              curve_icu_original = curve_boot_icu_original,
              curve_icu_trimmed = curve_boot_icu_trimmed,
              curve_aj_dc_original = curve_boot_aj_dc_original,
              curve_aj_dc_trimmed = curve_boot_aj_dc_trimmed,
              curve_aj_icu_original = curve_boot_aj_icu_original,
              curve_aj_icu_trimmed = curve_boot_aj_icu_trimmed))
}
# =============================================================================
# 10.  BOOTSTRAPPING RESULTS
# =============================================================================

#Function to summarize a bootstrapped set of summary statistics
create_boot_summary_table <- function(df) {
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
        P2_5      = quantile(x_clean, 0.025),
        P5        = quantile(x_clean, 0.05),
        P25       = quantile(x_clean, 0.25),
        Median    = median(x_clean),
        Mean      = mean(x_clean),
        SD        = sd(x_clean),
        P75       = quantile(x_clean, 0.75),
        P95       = quantile(x_clean, 0.95),
        P97_5     = quantile(x_clean, 0.975),
        Max       = max(x_clean),
        N         = length(x_clean),
        N_Missing = sum(is.na(x)),
        row.names = NULL
      )
    } else {
      NULL  # skip non-numeric columns
    }
  }))
  
  return(stats)
}

#Function to save results
save_boots_excel <- function(boot_df, label_in) {
  # Create workbook and styled sheet
  wb <- createWorkbook()
  n_col <- 15
  
  #Divide into the four options and run through it for them all
  sub_boot_df <- filter(boot_df, type == "simple" & trim=="original")
  sub_name <- "Simple_Original"
  addWorksheet(wb, sub_name)
  writeData(wb, sub_name, create_boot_summary_table(sub_boot_df),
            startRow = 1, startCol = 1)
  setColWidths(wb, sub_name, cols = 1:n_col, widths = "auto")
  
  sub_boot_df <- filter(boot_df, type == "simple" & trim=="trimmed")
  sub_name <- "Simple_Trimmed"
  addWorksheet(wb, sub_name)
  writeData(wb, sub_name, create_boot_summary_table(sub_boot_df),
            startRow = 1, startCol = 1)
  setColWidths(wb, sub_name, cols = 1:n_col, widths = "auto")
  
  sub_boot_df <- filter(boot_df, type == "MV" & trim=="original")
  sub_name <- "MV_Original"
  addWorksheet(wb, sub_name)
  writeData(wb, sub_name, create_boot_summary_table(sub_boot_df),
            startRow = 1, startCol = 1)
  setColWidths(wb, sub_name, cols = 1:n_col, widths = "auto")
  
  sub_boot_df <- filter(boot_df, type == "MV" & trim=="trimmed")
  sub_name <- "MV_Trimmed"
  addWorksheet(wb, sub_name)
  writeData(wb, sub_name, create_boot_summary_table(sub_boot_df),
            startRow = 1, startCol = 1)
  setColWidths(wb, sub_name, cols = 1:n_col, widths = "auto")
  
  output_path <- file.path(output_folder, "final",
                           make_filename("bootstrap_summary", label_in, ext = 'xlsx'))
  
  saveWorkbook(wb, output_path, overwrite = TRUE)
  message("Saved bootstrap summary to: ", output_path)
}

# =============================================================================
# 11.  BOOTSTRAPPED SURVIVAL CURVE
# =============================================================================
plot_marginal_curves <- function(prime_curve, boot_curve, time_grid_in,
                                 title_in = "", file_name = "fg_curve_unspecified")
{
  # ---- Percentile CIs at each time point ----
  ci_E <- apply(boot_curve[,,"E"], 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  ci_N <- apply(boot_curve[,,"N"], 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  
  plot_df <- tibble(
    time      = time_grid_in,
    pred_E    = prime_curve$pred_E,
    pred_E_lo = ci_E[1, ],
    pred_E_hi = ci_E[2, ],
    pred_N    = prime_curve$pred_N,
    pred_N_lo = ci_N[1, ],
    pred_N_hi = ci_N[2, ]
  )
  
  # ---- Plot ----
  p_mortality_curve <- ggplot(plot_df, aes(x = time)) +
    geom_ribbon(aes(ymin = pred_E_lo, ymax = pred_E_hi, fill = "Clone E"), alpha = 0.2) +
    geom_ribbon(aes(ymin = pred_N_lo, ymax = pred_N_hi, fill = "Clone N"), alpha = 0.2) +
    geom_line(aes(y = pred_E, color = "Clone E"), linewidth = 1) +
    geom_line(aes(y = pred_N, color = "Clone N"), linewidth = 1) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Days since IMV", y = "Predicted Cumulative Incidence",
         color = "Clone", fill = "Clone",
         title = title_in) +
    theme_minimal()
  
  ggsave(file.path(output_folder, "final", "graphs",file_name),
         plot = p_mortality_curve, width = 8, height = 5)
}
# =============================================================================
# 12.  FULL PIPELINE
# Run full pipeline on a data set (or subset)
# =============================================================================

run_pipeline <- function(pipe_in_df,label_in) {
  
  #Create and analyze weights for primary data
  print(paste0(label_in,": Creating Weights"))
  prime_df <- clone_and_weight(pipe_in_df, log_summary = TRUE)
  
  #Note weights_E and weights_N get save globaly by the function above
  print(paste0(label_in,": Analyzing Weights"))
  analyze_ipcw(prime_df, weights_E, weights_N, label_in)
  
  #Plot weighted mortality AJ curves (no regression)
  curve_prime_aj_dc_original <- get_aj_curve(prime_df$clone,
                                       prime_df$dc_fg_time,
                                       prime_df$dc_fg_cause,
                                       prime_df$IPCW,
                                       time_grid_dc)
  curve_prime_aj_dc_trimmed <- get_aj_curve(prime_df$clone,
                                       prime_df$dc_fg_time,
                                       prime_df$dc_fg_cause,
                                       prime_df$IPCW_trim,
                                       time_grid_dc)
  curve_prime_aj_icu_original <- get_aj_curve(prime_df$clone,
                                       prime_df$icu_fg_time,
                                       prime_df$icu_fg_cause,
                                       prime_df$IPCW,
                                       time_grid_icu)
  curve_prime_aj_icu_trimmed <- get_aj_curve(prime_df$clone,
                                      prime_df$icu_fg_time,
                                      prime_df$icu_fg_cause,
                                      prime_df$IPCW_trim,
                                      time_grid_icu)
  
  #Run the four outcomes types on the primary data.
  #For each outcome save the regression results and the curves (MV only)
  
  #-------original / simple
  print(paste0(label_in,": Primary Data Outcome Models"))
  outcome_now <- model_outcomes(prime_df, 0, type_reg = "simple",
                               trimmed_weights = FALSE)
  sumarize_regressions(outcome_now,
                       make_filename("regression_results",label_in,
                                     trim_status="original",
                                     model_type = "simple",
                                     ext="xlsx"))
  
  #-------original / MV
  outcome_now <- model_outcomes(prime_df, 0, type_reg = "MV",
                                trimmed_weights = FALSE)
  sumarize_regressions(outcome_now,
                       make_filename("regression_results",label_in,
                                     trim_status="original",
                                     model_type = "MV",
                                     ext="xlsx"))
  curve_prime_dc_original <- get_marginal_curve(fit_dead_fg, prime_df, time_grid_dc)
  curve_prime_icu_original <- get_marginal_curve(fit_icu_fg, prime_df, time_grid_icu)
  
  #-------trimmed / simple
  outcome_now <- model_outcomes(prime_df, 0, type_reg = "simple",
                                trimmed_weights = TRUE)
  sumarize_regressions(outcome_now,
                       make_filename("regression_results",label_in,
                                     trim_status="trimmed",
                                     model_type = "simple",
                                     ext="xlsx"))
  
  #-------trimmed / MV
  outcome_now <- model_outcomes(prime_df, 0, type_reg = "MV",
                                trimmed_weights = TRUE)
  sumarize_regressions(outcome_now,
                       make_filename("regression_results",label_in,
                                     trim_status="trimmed",
                                     model_type = "MV",
                                     ext="xlsx"))
  curve_prime_dc_trimmed <- get_marginal_curve(fit_dead_fg, prime_df, time_grid_dc)
  curve_prime_icu_trimmed <- get_marginal_curve(fit_icu_fg, prime_df, time_grid_icu)
  
  #Run the bootstrapping and save results
  print(paste0(label_in,": Running Bootstrapping for ",resample_N," re-samples."))
  boots_results <- bootstrap_all(pipe_in_df)
  save_boots_excel(boots_results$outcome_df, label_in)
  
  #Create FG predicted survival curves
  print(paste0(label_in,": Plotting FG Curves"))
  plot_marginal_curves(curve_prime_dc_original,
                       boots_results$curve_dc_original,
                       time_grid_dc,
                       title_in = "Predicted Hospital Mortality CIF (original)",
                       file_name = make_filename("mortality_curve",label_in,
                                                 trim_status="original",
                                                 model_type = "MV",
                                                 ext="pdf"))
  plot_marginal_curves(curve_prime_dc_trimmed,
                       boots_results$curve_dc_trimmed,
                       time_grid_dc,
                       title_in = "Predicted Hospital Mortality CIF (trimmed)",
                       file_name = make_filename("mortality_curve",label_in,
                                                 trim_status="trimmed",
                                                 model_type = "MV",
                                                 ext="pdf"))
  
  plot_marginal_curves(curve_prime_icu_original,
                       boots_results$curve_icu_original,
                       time_grid_icu,
                       title_in = "Predicted ICU LOS CIF (original)",
                       file_name = make_filename("icu_curve",label_in,
                                                 trim_status="original",
                                                 model_type = "MV",
                                                 ext="pdf"))
  plot_marginal_curves(curve_prime_icu_trimmed,
                       boots_results$curve_icu_trimmed,
                       time_grid_icu,
                       title_in = "Predicted ICU LOS CIF (trimmed)",
                       file_name = make_filename("icu_curve",label_in,
                                                 trim_status="trimmed",
                                                 model_type = "MV",
                                                 ext="pdf"))
  
  #Create AJ weighted survival curves
  print(paste0(label_in,": Plotting AJ Curves"))
  plot_marginal_curves(curve_prime_aj_dc_original,
                       boots_results$curve_aj_dc_original,
                       time_grid_dc,
                       title_in = "Weighted Hospital Mortality CIF (original)",
                       file_name = make_filename("mortality_curve",label_in,
                                                 trim_status="original",
                                                 model_type = "AJ",
                                                 ext="pdf"))
  plot_marginal_curves(curve_prime_aj_dc_trimmed,
                       boots_results$curve_aj_dc_trimmed,
                       time_grid_dc,
                       title_in = "Weighted Hospital Mortality CIF (trimmed)",
                       file_name = make_filename("mortality_curve",label_in,
                                                 trim_status="trimmed",
                                                 model_type = "AJ",
                                                 ext="pdf"))
  
  plot_marginal_curves(curve_prime_aj_icu_original,
                       boots_results$curve_aj_icu_original,
                       time_grid_icu,
                       title_in = "Weighted ICU LOS CIF (original)",
                       file_name = make_filename("icu_curve",label_in,
                                                 trim_status="original",
                                                 model_type = "AJ",
                                                 ext="pdf"))
  plot_marginal_curves(curve_prime_aj_icu_trimmed,
                       boots_results$curve_aj_icu_trimmed,
                       time_grid_icu,
                       title_in = "Weighted ICU LOS CIF (trimmed)",
                       file_name = make_filename("icu_curve",label_in,
                                                 trim_status="trimmed",
                                                 model_type = "AJ",
                                                 ext="pdf"))
  
  print(paste0(label_in,": DONE"))
}

run_pipeline(bin_df,"ALL")

# =============================================================================
# 13.  SUB GROUP ANALYSIS
# =============================================================================

if (run_sub_group) { 
  bin_65 <- bin_df %>% filter( (age >= 65) & (age < 75))
  run_pipeline(bin_65,"65")
  
  bin_75 <- bin_df %>% filter( (age >= 75) & (age < 85))
  run_pipeline(bin_75,"75")
  
  bin_85 <- bin_df %>% filter(age >= 85)
  run_pipeline(bin_85,"85")

}

# =============================================================================
# END
# =============================================================================
sink(type = "message")
sink()
close(sink_log)