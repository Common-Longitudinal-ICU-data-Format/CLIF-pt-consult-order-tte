# =============================================================================
# 4_ccw_MV_revised.R
# Clone-Censor Weighting (CCW) implementation
#
# KEY CHANGE FROM PRIOR VERSION:
# Previously we fit ONE logistic GLM per clone strategy, including time_bin_f
# as a factor covariate to handle all time points simultaneously.
#
# We now follow Webster-Clark et al. (2025, Pharmacoepidemiology & Drug Safety):
#   - Fit a SEPARATE GLM for each time_bin × clone combination
#     (i.e., one model per time point per strategy)
#   - Response is CENSORING (PT_censor_*), not uncensoring
#     (paper: Censor(t) ~ covariates(t); we get P(uncensored) = 1 - predicted)
#   - For Clone E only, apply the "pt_now" (analogous to paper's "recentstart")
#     logic when assigning interval weights:
#       * pt_now == 1  & uncensored  →  weight = 1 / P(uncensored)
#       * pt_now == 0  & uncensored  →  weight = 1   (no upweighting needed)
#       * censored                   →  weight = 0
#   - For Clone N: interval weight = 1 / P(uncensored) for all uncensored rows
#     (no pt_now logic; the paper's recentstart only applies to the treatment arm)
#   - Cumulative weight = product of all interval weights up to and including
#     the current time bin (matching paper's multiplicative IPCW)
#
# This yields up to N_time_bins × 2 GLMs rather than 2 GLMs.
# =============================================================================

# ---- Packages ----------------------------------------------------------------
packages <- c("tidyverse", "pscl", "ggplot2", "dplyr", "openxlsx",
              "tibble", "cobalt", "this.path", "glue")

installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])

library(tidyverse); library(pscl); library(ggplot2); library(dplyr); library(glue)
library(openxlsx); library(tibble); library(cobalt); library(this.path)

# ---- Paths -------------------------------------------------------------------
setwd(dirname(this.path()))
work_dir      <- normalizePath("..")
output_folder <- file.path(work_dir, "output")

# =============================================================================
# 1.  LOAD & PREP DATA
# =============================================================================
data <- read.csv(file.path(output_folder, "intermediate",
                           "block_and_time_bins_for_stats.csv"))
colnames(data)

# ---- Factorise -------------------------------
fac_vars <- c("sex_category", "race_category", "ethnicity_category",
              "language_category", "ICU_type")
data[fac_vars] <- lapply(data[fac_vars], function(x) {
  x[x == ""] <- NA
  factor(x)
})

# ---- Covariate lists ---------------------------------------------------------
# Fixed (baseline) covariates — same as before
base_vars <- c("age", "sex_category", "race_category", "ethnicity_category",
               "weight_kg", "language_category","elixhauser_age_adj","ICU_type")

# Time-varying covariates measured at each time_bin
tv_vars <- c("heart_rate_mean", "map_mean", "fio2_set_mean", "peep_set_mean",
             "pressor_flag", "paralytics_flag")

# Outcome Variables
out_vars <- c("vent_free_days","icu_los_days","is_dead_hosp","is_dead_30","is_dead_365")

#All covars, used to determine complete case analysis.
all_covars <- c(base_vars, tv_vars)

# Combined covariates & columns needed for complete analysis of the CCW.
all_vars <- c("encounter_block","time_bin","bin_start","bin_end","pt_order",
                "pt_now","pt_post48_IMV",base_vars, tv_vars)

df <- subset(data, select = all_vars)
df <- df[order(df$encounter_block, df$time_bin), ]

# Block level data frame for outcomes regression later
block_df <- subset(data, select = c("encounter_block", base_vars, out_vars))
block_df <- block_df[!duplicated(block_df$encounter_block), ]
                   
# ---- Censor indicators -------------------------------------------------------
# Clone N: censored from the first PT bin onward (pt_order fills forward)
df$PT_censor_N <- df$pt_order

# Clone E: censored only at the final bin (bin_end == 48) if PT never occurred
df$pt_post48_IMV = df$pt_post48_IMV == "True"
df$PT_censor_E <- ifelse((!df$pt_post48_IMV) & df$bin_end == 48, 1, 0)

# ---- Build clone-specific datasets -------------------------------------------
# Clone N: keep rows up to and including the first PT bin (or all rows if no PT)
keep_N <- df$PT_censor_N == 0 | df$pt_now == 1
df_N   <- df[keep_N, ]
d_c <- "PT_censor_E"
df_N <- df_N %>% select(-d_c)   # drop censor E column
d_c <- "PT_censor_N"
df_E   <- df %>% select(-d_c)

# =============================================================================
# 2.  COMPLETE-CASE FILTER
# =============================================================================
vars_needed_N <- c("PT_censor_N", all_covars)
df_N <- df_N[complete.cases(df_N[, vars_needed_N]), ]
df_N <- df_N[order(df_N$encounter_block, df_N$time_bin), ]

vars_needed_E <- c("PT_censor_E", all_covars)
df_E <- df_E[complete.cases(df_E[, vars_needed_E]), ]
df_E <- df_E[order(df_E$encounter_block, df_E$time_bin), ]

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
#         Clone E : 1 / P(uncensored)   if pt_now == 1  AND uncensored
#                   1                    if pt_now == 0  AND uncensored
#                   0                    if censored
#    e) Collect interval weights, then compute cumulative product per encounter
#
#  Why separate models per bin?
#    Each time bin may have a different risk set (some encounters have been
#    censored or have experienced the outcome in prior bins).  Fitting separate
#    models avoids the assumption that the relationship between covariates and
#    censoring is identical across all bins, which is the same rationale
#    Webster-Clark use for their 0–30 / 30–90 / >90 day window models.
# =============================================================================

fit_interval_weights <- function(clone_df,
                                 censor_col,   # "PT_censor_N" or "PT_censor_E"
                                 base_vars,
                                 tv_vars,
                                 pt_now_logic = FALSE,
                                 stabilize = FALSE) {
  # pt_now_logic = TRUE  → Clone E weighting (use pt_now flag)
  # pt_now_logic = FALSE → Clone N weighting (always 1/p_uncens when uncensored)

  rhs_formula <- paste(c(base_vars, tv_vars), collapse = " + ")
  form        <- as.formula(paste(censor_col, "~", rhs_formula))

  time_bins <- sort(unique(clone_df$time_bin))
  
  # Stabilized weights
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

  interval_weights
}

# ---- Run per-bin weighting for each clone ------------------------------------

message("Fitting per-time-bin GLMs for Clone N ...")
weights_N <- fit_interval_weights(
  clone_df     = df_N,
  censor_col   = "PT_censor_N",
  base_vars    = base_vars,
  tv_vars      = tv_vars,
  pt_now_logic = FALSE,
  stabilize = FALSE
)

message("Fitting per-time-bin GLMs for Clone E ...")
weights_E <- fit_interval_weights(
  clone_df     = df_E,
  censor_col   = "PT_censor_E",
  base_vars    = base_vars,
  tv_vars      = tv_vars,
  pt_now_logic = FALSE,
  stabilize = FALSE
)

# Merge weights back onto the clone data frames for downstream use/diagnostics
df_N <- df_N %>%
  left_join(weights_N %>% dplyr::select(encounter_block, time_bin,
                                         p_cens, p_uncens, interval_wt, IPCW),
            by = c("encounter_block", "time_bin"))

df_E <- df_E %>%
  left_join(weights_E %>% dplyr::select(encounter_block, time_bin,
                                         p_cens, p_uncens, interval_wt, IPCW),
            by = c("encounter_block", "time_bin"))

# =============================================================================
# 4.  FINAL WEIGHTS
#     Take the cumulative IPCW at the LAST observed time_bin for each encounter,
#     keeping only encounters that completed the strategy (uncensored at exit).
# =============================================================================
make_final_weight <- function(dat, clone_label, censor_col) {
  dat %>%
    group_by(encounter_block) %>%
    slice_tail(n = 1) %>%                        # last observed bin
    filter(!!sym(censor_col) == 0) %>%           # must be uncensored at exit
    ungroup() %>%
    dplyr::select(encounter_block, IPCW) %>%
    rename(weight = IPCW) %>%
    mutate(clone = clone_label)
}

w_E <- make_final_weight(df_E, "E", "PT_censor_E")
w_N <- make_final_weight(df_N, "N", "PT_censor_N")

# =============================================================================
# 5.  IPCW TRAJECTORY PLOT  (diagnostic — same structure as before)
# =============================================================================
ipcw_long <- bind_rows(
  df_N %>% mutate(clone = "N"),
  df_E %>% mutate(clone = "E")
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
ggsave(file.path(output_folder, "final", "graphs", "original_IPCW_trajectory.pdf"),
       plot = p_ipcw_time, width = 8, height = 5)

p_ipcw_time1 <- ggplot(ipcw_long, aes(x = time_bin_f, y = IPCW_trim, fill = clone)) +
  geom_boxplot(outlier.alpha = 0.25, outlier.size = 0.8,
               position = position_dodge(width = 0.75)) +
  theme_bw() +
  labs(title = "Trajectory of Trimmed Unstabilized IPCW Over Time",
       x = "Time bin", y = "Trimmed Unstabilized IPCW", fill = "Clone")
ggsave(file.path(output_folder, "final", "graphs", "trim_IPCW_trajectory.pdf"),
       plot = p_ipcw_time1, width = 8, height = 5)

# =============================================================================
# 6.  OUTCOME MODEL SETUP  (unchanged from original)
# =============================================================================

block_df$vent_free_days <- as.integer(block_df$vent_free_days)
block_df$icu_los_days   <- as.integer(block_df$icu_los_days)
block_df$is_dead_hosp   <- ifelse(block_df$is_dead_hosp  == "True", 1, 0)
block_df$is_dead_30     <- ifelse(block_df$is_dead_30    == "True", 1, 0)
block_df$is_dead_365    <- ifelse(block_df$is_dead_365   == "True", 1, 0)

# Join final weights onto the outcome data (one row per encounter × clone)
outcome_df <- bind_rows(
  block_df %>% mutate(clone = "E") %>% left_join(w_E, by = c("encounter_block", "clone")),
  block_df %>% mutate(clone = "N") %>% left_join(w_N, by = c("encounter_block", "clone"))
) %>%
  mutate(clone = factor(clone, levels = c("N", "E")))

summary(outcome_df$weight)
quantile(outcome_df$weight, c(0, .01, .05, .25, .5, .75, .95, .99, 1), na.rm = TRUE)
table(outcome_df$clone)
sum(is.na(outcome_df$weight))
sum(!is.finite(outcome_df$weight))

# Diagnose missing weights
missing_weights <- outcome_df[is.na(outcome_df$weight), ]
tapply(missing_weights$encounter_block, missing_weights$clone,
       function(x) length(unique(x)))

outcome_df <- outcome_df %>%
  filter(!is.na(weight), is.finite(weight), weight > 0)

# Trim weights at 1st / 99th percentile
w_cut <- quantile(outcome_df$weight, probs = c(0.01, 0.99), na.rm = TRUE)
outcome_df <- outcome_df %>%
  mutate(weight_trim = pmin(pmax(weight, w_cut[[1]]), w_cut[[2]]))

summary(outcome_df$weight_trim)
quantile(outcome_df$weight_trim, probs = c(0, .01, .05, .25, .5, .75, .95, .99, 1))

# Weight distribution plots
g <- ggplot(outcome_df, aes(x = weight, fill = clone)) +
  geom_histogram(bins = 60, alpha = 0.5, position = "identity") +
  theme_bw() +
  labs(title = "Distribution of final weights by clone",
       x = "Final weight", y = "Count", fill = "Clone Group")
ggsave(file.path(output_folder, "final", "graphs", "original_final_IPCW.pdf"),
       plot = g, width = 7, height = 5)

g1 <- ggplot(outcome_df, aes(x = weight_trim, fill = clone)) +
  geom_histogram(bins = 60, alpha = 0.5, position = "identity") +
  theme_bw() +
  labs(title = "Distribution of trimmed final weights by clone",
       x = "Final weight", y = "Count", fill = "Clone Group")
ggsave(file.path(output_folder, "final", "graphs", "trim_final_IPCW.pdf"),
       plot = g1, width = 7, height = 5)

# =============================================================================
# 7.  BASELINE COVARIATE BALANCE CHECK  (unchanged)
# =============================================================================
bal_ccw <- bal.tab(x = outcome_df[, base_vars], treat = outcome_df$clone,
                   weights = outcome_df$weight_trim, method = "weighting",
                   estimand = "ATE", s.d.denom = "pooled", un = TRUE)
print(bal_ccw)

p_balance <- love.plot(bal_ccw, stats = "mean.diffs", abs = TRUE,
                       thresholds = c(m = 0.1), var.order = "unadjusted",
                       stars = "raw", sample.names = c("Unweighted", "Weighted"),
                       title = "Baseline Covariate Balance Before and After IPCW")
print(p_balance)
ggsave(file.path(output_folder, "final", "graphs", "balance_plot_IPCW.pdf"),
       plot = p_balance, width = 8, height = 6)

# =============================================================================
# 8.  OUTCOME MODELS  (unchanged)
# =============================================================================

##### VFD: ZINB #####
fit_vfd <- zeroinfl(vent_free_days ~ clone | 1, data = outcome_df,
                    dist = "negbin", weights = weight_trim)
summary(fit_vfd)

##### ICU LOS: Poisson #####
fit_icu_los <- glm(icu_los_days ~ clone, data = outcome_df,
                   family = poisson(), weights = weight_trim)
summary(fit_icu_los)

##### Hospital mortality: Binary #####
fit_dead_hosp <- glm(is_dead_hosp ~ clone, data = outcome_df,
                     family = binomial(), weights = weight_trim)
summary(fit_dead_hosp)

##### 30-day mortality: Binary #####
fit_dead_30 <- glm(is_dead_30 ~ clone, data = outcome_df,
                   family = binomial(), weights = weight_trim)
summary(fit_dead_30)

##### 1-year mortality: Binary #####
fit_dead_365 <- glm(is_dead_365 ~ clone, data = outcome_df,
                    family = binomial(), weights = weight_trim)
summary(fit_dead_365)

# =============================================================================
# 9.  RESULTS ORGANISATION  (unchanged)
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

tab_vfd      <- extract_zeroinfl_table(fit_vfd,      "vent_free_days")
tab_icu_los  <- extract_glm_table(fit_icu_los,       "icu_los_days")
tab_dead_hosp<- extract_glm_table(fit_dead_hosp,     "is_dead_hosp")
tab_dead_30  <- extract_glm_table(fit_dead_30,       "is_dead_30")
tab_dead_365 <- extract_glm_table(fit_dead_365,      "is_dead_365")

standardized_contrast <- function(fit, data, outcome_name, clone_var = "clone") {
  dE <- data; dN <- data
  dE[[clone_var]] <- factor("E", levels = levels(data[[clone_var]]))
  dN[[clone_var]] <- factor("N", levels = levels(data[[clone_var]]))
  pred_E <- predict(fit, newdata = dE, type = "response")
  pred_N <- predict(fit, newdata = dN, type = "response")
  tibble(outcome        = outcome_name,
         mean_pred_E    = mean(pred_E, na.rm = TRUE),
         mean_pred_N    = mean(pred_N, na.rm = TRUE),
         diff_E_minus_N = mean(pred_E, na.rm = TRUE) - mean(pred_N, na.rm = TRUE),
         ratio_E_over_N = mean(pred_E, na.rm = TRUE) / mean(pred_N, na.rm = TRUE))
}

pred_vfd  <- standardized_contrast(fit_vfd,       outcome_df, "vent_free_days")
pred_icu  <- standardized_contrast(fit_icu_los,   outcome_df, "icu_los_days")
pred_hosp <- standardized_contrast(fit_dead_hosp, outcome_df, "is_dead_hosp")
pred_30   <- standardized_contrast(fit_dead_30,   outcome_df, "is_dead_30")
pred_365  <- standardized_contrast(fit_dead_365,  outcome_df, "is_dead_365")
pred_contrast_tab <- bind_rows(pred_vfd, pred_icu, pred_hosp, pred_30, pred_365)

wb <- createWorkbook()
addWorksheet(wb, "VFD");             writeData(wb, "VFD",              tab_vfd)
addWorksheet(wb, "ICU_LOS");         writeData(wb, "ICU_LOS",          tab_icu_los)
addWorksheet(wb, "Hosp");            writeData(wb, "Hosp",             tab_dead_hosp)
addWorksheet(wb, "30Day");           writeData(wb, "30Day",            tab_dead_30)
addWorksheet(wb, "1Year");           writeData(wb, "1Year",            tab_dead_365)
addWorksheet(wb, "Predicted_Contrast"); writeData(wb, "Predicted_Contrast", pred_contrast_tab)
saveWorkbook(wb, file = file.path(output_folder, "final", "ccw_IPCW_results.xlsx"),
             overwrite = TRUE)

# =============================================================================
# 10.  MULTIVARIATE OUTCOME MODELS  (unchanged)
# =============================================================================
mv_rhs <- paste(c("clone", base_vars), collapse = " + ")

fit_vfd_mv       <- zeroinfl(as.formula(paste("vent_free_days ~", mv_rhs, "| 1")),
                             data = outcome_df, dist = "negbin", weights = weight_trim)
fit_icu_los_mv   <- glm(as.formula(paste("icu_los_days ~", mv_rhs)),
                        data = outcome_df, family = poisson(),   weights = weight_trim)
fit_dead_hosp_mv <- glm(as.formula(paste("is_dead_hosp ~", mv_rhs)),
                        data = outcome_df, family = binomial(),  weights = weight_trim)
fit_dead_30_mv   <- glm(as.formula(paste("is_dead_30 ~", mv_rhs)),
                        data = outcome_df, family = binomial(),  weights = weight_trim)
fit_dead_365_mv  <- glm(as.formula(paste("is_dead_365 ~", mv_rhs)),
                        data = outcome_df, family = binomial(),  weights = weight_trim)

tab_vfd_mv       <- extract_zeroinfl_table(fit_vfd_mv,       "vent_free_days")
tab_icu_los_mv   <- extract_glm_table(fit_icu_los_mv,        "icu_los_days")
tab_dead_hosp_mv <- extract_glm_table(fit_dead_hosp_mv,      "is_dead_hosp")
tab_dead_30_mv   <- extract_glm_table(fit_dead_30_mv,        "is_dead_30")
tab_dead_365_mv  <- extract_glm_table(fit_dead_365_mv,       "is_dead_365")

pred_vfd_mv  <- standardized_contrast(fit_vfd_mv,       outcome_df, "vent_free_days")
pred_icu_mv  <- standardized_contrast(fit_icu_los_mv,   outcome_df, "icu_los_days")
pred_hosp_mv <- standardized_contrast(fit_dead_hosp_mv, outcome_df, "is_dead_hosp")
pred_30_mv   <- standardized_contrast(fit_dead_30_mv,   outcome_df, "is_dead_30")
pred_365_mv  <- standardized_contrast(fit_dead_365_mv,  outcome_df, "is_dead_365")
pred_contrast_tab_mv <- bind_rows(pred_vfd_mv, pred_icu_mv, pred_hosp_mv,
                                  pred_30_mv, pred_365_mv)

wb_mv <- createWorkbook()
addWorksheet(wb_mv, "VFD");             writeData(wb_mv, "VFD",              tab_vfd_mv)
addWorksheet(wb_mv, "ICU_LOS");         writeData(wb_mv, "ICU_LOS",          tab_icu_los_mv)
addWorksheet(wb_mv, "Hosp");            writeData(wb_mv, "Hosp",             tab_dead_hosp_mv)
addWorksheet(wb_mv, "30Day");           writeData(wb_mv, "30Day",            tab_dead_30_mv)
addWorksheet(wb_mv, "1Year");           writeData(wb_mv, "1Year",            tab_dead_365_mv)
addWorksheet(wb_mv, "Predicted_Contrast"); writeData(wb_mv, "Predicted_Contrast", pred_contrast_tab_mv)
saveWorkbook(wb_mv,
             file = file.path(output_folder, "final", "ccw_IPCW_results_multivariate.xlsx"),
             overwrite = TRUE)
