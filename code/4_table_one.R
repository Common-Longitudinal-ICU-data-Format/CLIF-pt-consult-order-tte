## ------------------------------------------------------------------
## Calculate Elixhouser
## Table 1
## CIF graph and PT per time bin graph
## ------------------------------------------------------------------

# ---- Packages ----------------------------------------------------------------
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("renv")
library(this.path)
setwd(dirname(this.path()))
work_dir      <- normalizePath("..")
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::load(project = work_dir)

packages <- c("arrow","dplyr","comorbidity")
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])

library(arrow)
library(dplyr)
library(comorbidity)

## ------------------------------------------------------------------
## 1. Read in the data
## ------------------------------------------------------------------

# ---- Paths -------------------------------------------------------------------
setwd(dirname(this.path()))
work_dir      <- normalizePath("..")
output_folder <- file.path(work_dir, "output")
block_file_path <- file.path(output_folder, "intermediate",
                             "block_df_3_end.parquet")
diag_file_path <- file.path(output_folder, "intermediate",
                             "diag_codes.parquet")

block_df <- read_parquet(block_file_path)          # one row per encounter_block
diag_df  <- read_parquet(diag_file_path)       # long format, one row per diagnosis code

## Sanity check what's in diagnosis_code_format -- this needs to map
## cleanly onto "9" (ICD-9-CM) vs "10" (ICD-10) for the comorbidity package
diag_df %>% count(diagnosis_code_format)

## ------------------------------------------------------------------
## 2. Clean codes
## ------------------------------------------------------------------
## The comorbidity package requires codes in upper case with NO
## punctuation (no decimal points), e.g. "E1140" not "E11.40".

diag_df <- diag_df %>%
  mutate(
    diagnosis_code_clean = toupper(gsub("[[:punct:][:space:]]", "", diagnosis_code))
  )

## Normalize the format flag -- EDIT this to match your actual values,
## e.g. if diagnosis_code_format is "ICD9"/"ICD10", "9"/"10", 9/10, etc.
diag_df <- diag_df %>%
  mutate(
    icd_version = case_when(
      grepl("ICD10CM", as.character(diagnosis_code_format)) ~ "10",
      grepl("ICD9CM",  as.character(diagnosis_code_format))  ~ "9",
      TRUE ~ NA_character_
    )
  )

if (any(is.na(diag_df$icd_version))) {
  warning(sprintf(
    "%d rows had an unrecognized diagnosis_code_format and were dropped -- check the case_when() mapping above.",
    sum(is.na(diag_df$icd_version))
  ))
}

## ------------------------------------------------------------------
## 3. Run comorbidity() separately for ICD-9 and ICD-10 subsets,
##    since the map argument is version-specific
## ------------------------------------------------------------------

diag9  <- diag_df %>% filter(icd_version == "9")
diag10 <- diag_df %>% filter(icd_version == "10")

elix9 <- if (nrow(diag9) > 0) {
  out <- comorbidity(
    x    = diag9,
    id   = "encounter_block",
    code = "diagnosis_code_clean",
    map  = "elixhauser_icd9_quan",
    assign0 = TRUE
  )
  # comorbidity() returns a data.frame with an extra "comorbidity" S3 class
  # (used internally for label-aware printing). That class has its own
  # `[.comorbidity` method which breaks dplyr's group_by()/across() column
  # resolution downstream -- strip it back to a plain data.frame/tibble.
  class(out) <- "data.frame"
  tibble::as_tibble(out)
} else NULL

elix10 <- if (nrow(diag10) > 0) {
  out <- comorbidity(
    x    = diag10,
    id   = "encounter_block",
    code = "diagnosis_code_clean",
    map  = "elixhauser_icd10_quan",
    assign0 = TRUE
  )
  class(out) <- "data.frame"
  tibble::as_tibble(out)
} else NULL

## Combine (both have the same 31 comorbidity flag columns + id)
elix_flags <- bind_rows(elix9, elix10)

## If a single encounter_block had both ICD-9 and ICD-10 coded diagnoses
## (e.g. spanning a coding transition), collapse to one row per block
## by taking the max of each flag (i.e. "present in either coding system")
elix_flags <- elix_flags %>%
  group_by(encounter_block) %>%
  summarise(across(everything(), ~ as.integer(any(. == 1))), .groups = "drop")

## ------------------------------------------------------------------
## 4. Score it
## ------------------------------------------------------------------
## weights = "vw"     -> van Walraven (2009) point system (most common single-number index)
## score() requires the "comorbidity" S3 class + a "map" attribute to know
## which weighting family applies.

elix_flags_for_scoring <- elix_flags
class(elix_flags_for_scoring) <- c("comorbidity", class(elix_flags_for_scoring))
attr(elix_flags_for_scoring, "map") <- "elixhauser_icd10_quan"

vw_score <- score(elix_flags_for_scoring, weights = "vw", assign0 = TRUE)

## Keep only encounter_block + a single summary score column -- not the
## 30 individual comorbidity flag columns
elix_score <- elix_flags %>%
  transmute(
    encounter_block,
    elixhauser = vw_score
  )

## ------------------------------------------------------------------
## 5. Join back onto block_df
## ------------------------------------------------------------------

block_df_scored <- block_df %>%
  left_join(elix_score, by = "encounter_block")

## Encounters with no diagnosis rows at all will have NA for elixhauser --
## they get assigned 0.
block_df_scored <- block_df_scored %>%
  mutate(elixhauser = tidyr::replace_na(elixhauser, 0))

## ------------------------------------------------------------------
## 6. Write out
## ------------------------------------------------------------------
out_file_path <- file.path(output_folder, "intermediate",
                             "block_df_4_end.parquet")
write_parquet(block_df_scored, out_file_path)

## ------------------------------------------------------------------
## 7. Organize Columns and Summarize
## ------------------------------------------------------------------
## Ported from 3_calculations.ipynb ("Organize Columns and Summarize"
## onward). Continues from block_df_scored (block_df + elixhauser score).

## ---- Column order / labels from config -------------------------------------
column_def_path <- file.path(work_dir, "config", "column_def.csv")
column_order    <- read.csv(column_def_path, stringsAsFactors = FALSE)
my_cols         <- column_order$name
rownames(column_order) <- column_order$name   # lookup by column name, like set_index('name')

block_df_final <- block_df_scored[, my_cols]

final_path <- file.path(output_folder, "intermediate", "block_df_final.parquet")
write_parquet(block_df_final, final_path)

## ---- Outcome flag ------------------------------------------------------------
early_col <- "pt_post48_IMV"
block_df_final$early_PT <- ifelse(block_df_final[[early_col]], "early_PT", "no_early_PT")

n_total <- sum(!is.na(block_df_final$encounter_block))
n_early <- sum(block_df_final[[early_col]], na.rm = TRUE)
n_not   <- n_total - n_early

## ---- SMD calculator ------------------------------------------------------------
calculate_smd <- function(group1, group2) {
  mean1 <- mean(group1)
  mean2 <- mean(group2)
  var1  <- var(group1)   # R's var() uses ddof = 1 by default, matching numpy's ddof=1
  var2  <- var(group2)
  n1    <- length(group1)
  n2    <- length(group2)
  pooled_sd <- sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
  smd <- (mean1 - mean2) / pooled_sd
  return(smd)
}

## ---- Build table1.csv ------------------------------------------------------------
table1_path <- file.path(output_folder, "final", "table1.csv")
if (file.exists(table1_path)) file.remove(table1_path)

con <- file(table1_path, open = "w")
cat(",,Overall,Early PT, No Early PT, P-value/SMD, Missing", file = con)

for (col in my_cols) {

  lab <- column_order[col, "description"]
  col_vals <- block_df_final[[col]]

  if (col == "encounter_block") {

    cat(sprintf("\nN,,%s,%s,%s,", n_total, n_early, n_not), file = con)

  } else if (is.character(col_vals) || is.factor(col_vals)) {

    cat(sprintf("\n%s", lab), file = con)
    col_chr <- as.character(col_vals)
    cats <- unique(col_chr[!is.na(col_chr)])
    tab  <- table(col_chr, block_df_final$early_PT, useNA = "no")
    p_value <- suppressWarnings(chisq.test(tab)$p.value)
    for (cc in cats) {
      if (nzchar(cc)) {   # skips the "" category, mirroring Python's `if cc:`
        cc_all   <- 100 * sum(tab[cc, ]) / sum(!is.na(col_chr))
        cc_early <- 100 * tab[cc, "early_PT"]    / sum(tab[, "early_PT"])
        cc_not   <- 100 * tab[cc, "no_early_PT"] / sum(tab[, "no_early_PT"])
        cat(sprintf("\n,%s,%.1f%%,%.1f%%,%.1f%%", cc, cc_all, cc_early, cc_not), file = con)
      }
    }
    cat(sprintf(", %.5f", p_value), file = con)

  } else if (is.logical(col_vals)) {

    if (sum(col_vals, na.rm = TRUE) > 0) {
      sub_df <- block_df_final[!is.na(col_vals), ]
      flag <- ifelse(sub_df[[col]], "TRUE", "FALSE")
      tab  <- table(flag, sub_df$early_PT)
      p_value <- suppressWarnings(chisq.test(tab)$p.value)
      cc_all   <- 100 * sum(tab["TRUE", ]) / sum(!is.na(sub_df[[col]]))
      cc_early <- 100 * tab["TRUE", "early_PT"]    / sum(tab[, "early_PT"])
      cc_not   <- 100 * tab["TRUE", "no_early_PT"] / sum(tab[, "no_early_PT"])
      cat(sprintf("\n%s,,%.2f%%,%.2f%%,%.2f%%,%.5f", lab, cc_all, cc_early, cc_not, p_value), file = con)
    } else {
      cat(sprintf("\n%s,,0.00%%,0.00%%,0.00%%,N/A", lab), file = con)
    }

  } else if (is.numeric(col_vals)) {

    cc_all   <- col_vals[!is.na(col_vals)]
    cc_early <- col_vals[block_df_final[[early_col]]  & !is.na(block_df_final[[early_col]])]
    cc_early <- cc_early[!is.na(cc_early)]
    cc_not   <- col_vals[!block_df_final[[early_col]] & !is.na(block_df_final[[early_col]])]
    cc_not   <- cc_not[!is.na(cc_not)]
    SMD <- calculate_smd(cc_early, cc_not)

    cat(sprintf("\n%s (Med & IQR),,%.2f  (%.2f - %.2f)",
                lab, median(cc_all), unname(quantile(cc_all, 0.25)), unname(quantile(cc_all, 0.75))),
        file = con)
    cat(sprintf(",%.2f  (%.2f - %.2f)",
                median(cc_early), unname(quantile(cc_early, 0.25)), unname(quantile(cc_early, 0.75))),
        file = con)
    cat(sprintf(",%.2f  (%.2f - %.2f)",
                median(cc_not), unname(quantile(cc_not, 0.25)), unname(quantile(cc_not, 0.75))),
        file = con)
    cat(sprintf(",%.5f", SMD), file = con)

  } else {

    cat(sprintf("\n%s,ERROR,,,,,", lab), file = con)

  }

  ## Missing data column
  mis_pct <- 100 * sum(is.na(col_vals)) / n_total
  cat(sprintf(",%.2f%%", mis_pct), file = con)
}

close(con)

## ------------------------------------------------------------------
## 8. CIF Graph
## ------------------------------------------------------------------

graphs_folder <- file.path(output_folder, "final", "graphs")
dir.create(graphs_folder, recursive = TRUE, showWarnings = FALSE)

yellow_list <- sort(block_df_final$yellow_time_eligibility_2h[!is.na(block_df_final$yellow_time_eligibility_2h)])
pt_list     <- sort(block_df_final$Time_first_PT[!is.na(block_df_final$Time_first_PT)])

cif_path <- file.path(graphs_folder, "CIF_Yellow_v_PT.png")
png(cif_path, width = 8, height = 6, units = "in", res = 150)
plot(yellow_list, seq_len(length(yellow_list)) - 1, type = "s", col = "orange",
     xlim = c(0, 72), ylim = c(0, max(length(yellow_list), length(pt_list)) - 1),
     xlab = "Hours from IMV Initiation", ylab = "Encounters",
     main = "Physioliogic Readiness versus PT Consult Order CIF")
lines(pt_list, seq_len(length(pt_list)) - 1, type = "s", col = "blue")
legend("bottomright", legend = c("Physiologic readiness (yellow)", "PT consult order"),
       col = c("orange", "blue"), lty = 1, bty = "n")
dev.off()

## ------------------------------------------------------------------
## 9. Time Bin Summary Graphs
## ------------------------------------------------------------------
## The closed time-bin data set was written out of 3_calculations.ipynb

time_bin_path <- file.path(output_folder, "intermediate", "time_bin_3_end.parquet")
time_bin_df   <- read_parquet(time_bin_path, stringsAsFactors = FALSE)

pt_time_bin_df <- time_bin_df %>%
  group_by(bin_start) %>%
  summarise(pt_order = sum(pt_order, na.rm = TRUE), .groups = "drop") %>%
  arrange(bin_start)

time_bin_graph_path <- file.path(graphs_folder, "Time_bin_PT.png")
png(time_bin_graph_path, width = 8, height = 6, units = "in", res = 150)
barplot(height = pt_time_bin_df$pt_order, names.arg = pt_time_bin_df$bin_start,
        col = "blue", xlab = "Hours from IMV Initiation", ylab = "Encounters",
        main = "Early PT Consult Order Prevalence Over Time")
dev.off()

## ------------------------------------------------------------------
## 10. Merging for CCW
## ------------------------------------------------------------------
## Create a merged block_df_final and time_bin_df to be used for stats.

mask_cols <- (column_order$name == "encounter_block") |
  (column_order$covariate %in% 1) |
  (column_order$outcome %in% 1) |
  !is.na(column_order$other)
stats_cols <- column_order[mask_cols, ]

stats_df <- block_df_final[, stats_cols$name]
stats_df <- inner_join(stats_df, time_bin_df, by = "encounter_block")

message(sprintf("Stats data set contains %d encounter_blocks.",
                 n_distinct(block_df_final$encounter_block)))
message(sprintf("Stats data set contains %d encounter_blocks.",
                 n_distinct(stats_df$encounter_block)))

stats_out_path <- file.path(output_folder, "intermediate", "block_and_time_bins_for_stats.parquet")
write_parquet(stats_df, stats_out_path)