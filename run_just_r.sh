#!/usr/bin/env bash
# ════════════════════════════════════════════════════════════════════════════════
#  run_just_r.sh — Execute only the R steps of the pipeline
#
#  Steps:
#    4. R       04_table_one.R
#    5. R       05_ccw.R
#
#  Usage:  bash run_just_r.sh
# ════════════════════════════════════════════════════════════════════════════════
set -euo pipefail

# ── colours ──────────────────────────────────────────────────────────────────
GREEN="\033[32m"; RED="\033[31m"; CYAN="\033[36m"; YELLOW="\033[33m"
BOLD="\033[1m"; RESET="\033[0m"

# ── paths ────────────────────────────────────────────────────────────────────
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_DIR="${PROJECT_ROOT}"
LOG_FILE="${LOG_DIR}/pipeline_${TIMESTAMP}.log"
mkdir -p "$LOG_DIR"

# ── logging ──────────────────────────────────────────────────────────────────
log() { echo -e "$1" | tee -a "$LOG_FILE"; }

log "${CYAN}${BOLD} CLIF PT Consult-order-tte Pipeline${RESET}"
log "Started: $(date)"
log "Log: ${LOG_FILE}"
log ""

# ── environment (renv) ─────────────────────────────────────────────────────────
if ! command -v Rscript >/dev/null 2>&1; then
  log "${RED}Rscript not found. Steps 4 & 5 could not be run.${RESET}"
  exit 1
else
  #Load renv
  cd "$PROJECT_ROOT"
  Rscript -e 'renv::restore(prompt = FALSE)' || { log "${RED}renv::restore failed${RESET}"; exit 1; }
  log "${CYAN}RENV Status:${RESET}"
  Rscript -e 'renv::status()'
  log "========== STARTING STEP 4: Table One =========="
  if Rscript code/4_table_one.R; then
    log "Step 4: Table One ran successfully"
  else
    log "${RED}Step 4: Table One FAILED${RESET}"
    exit 1
  fi

  log "========== STARTING STEP 5: CCW =========="
  if Rscript code/5_ccw.R; then
    log "Step 5: CCW ran successfully"
  else
    log "${RED}Step 5: CCW FAILED${RESET}"
    exit 1
  fi
fi

# ── output files ──────────────────────────────────────────────────────────────
log ""
log "Output files in ${PROJECT_ROOT}/output/final/:"
if [ -d "${PROJECT_ROOT}/output/final" ]; then
  # List generated files with sizes
  find "${PROJECT_ROOT}/output/final" -type f -newer "${LOG_FILE}" -exec ls -lh {} \; 2>/dev/null | \
    awk '{printf "  %-8s %s\n", $5, $NF}' | tee -a "$LOG_FILE" || true
  # If nothing newer, just list everything
  FILE_COUNT=$(find "${PROJECT_ROOT}/output/final" -type f | wc -l | tr -d ' ')
  log "  Total files: ${FILE_COUNT}"
else
  log "  (directory not yet created)"
fi

log ""
log "Full log: ${LOG_FILE}"
log "Finished: $(date)"
