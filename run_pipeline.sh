#!/usr/bin/env bash
# ════════════════════════════════════════════════════════════════════════════════
#  run_pipeline.sh — Execute the full CLIF-pt-consult-order-tte pipeline
#
#  Steps:
#    1. Python  01_cohort.ipynb -> Python
#    2. Python  02_data_gathering.ipynb -> Python
#    3. Python  03_calculations.ipynb -> Python
#    4. R       04_ccw.R
#
#  Usage:  bash run_pipeline.sh
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

# ── environment (uv) ─────────────────────────────────────────────────────────
if ! command -v uv >/dev/null 2>&1; then
  log "${RED}uv not found. Install it: https://docs.astral.sh/uv/getting-started/installation/${RESET}"
  exit 1
fi

log "Syncing dependencies with uv..."
uv sync --project "${PROJECT_ROOT}" 2>&1 | tee -a "$LOG_FILE"
log "${GREEN}Environment ready${RESET}"
log ""

# ── pipeline (cwd = code/ so relative paths work) ───────────────────────────
cd "${PROJECT_ROOT}/code"

# Python steps
log "========== STARTING STEP 1: COHORT =========="
uv run python 1_cohort.py
log "Step 1: Cohort ran"
log "========== STARTING STEP 2: DATA GATHERING =========="
uv run python 2_data_gathering.py
log "Step 2: Data Gathering ran"
log "========== STARTING STEP 3: CALCULATIONS =========="
uv run python 3_calculations.py
log "Step 3: Calculations ran"

# R steps
if command -v Rscript >/dev/null 2>&1; then
  run_step "04 CCW"       Rscript --vanilla 4_ccw.R
  log "Step 4: CCW ran"
else
  log "${YELLOW}Rscript not found — skipping R analysis.${RESET}"
  log "${YELLOW}Run manually: cd code && Rscript 4_ccw.R${RESET}"
fi

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
