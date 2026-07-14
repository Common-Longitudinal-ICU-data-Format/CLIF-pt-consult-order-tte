#!/usr/bin/env bash
#
# rebuild_uv_lock.sh
#
# Removes the existing uv.lock file and regenerates it based on the
# third-party packages imported by the project's Python files:
#   - pthelperfunctions.py
#   - 1_cohort.py
#   - 2_data_gathering.py
#   - 3_calculations.py
#
# Usage:
#   ./rebuild_uv_lock.sh
#
# This script always operates on the directory ONE LEVEL UP from wherever
# this script file itself lives (not the directory you happen to run it
# from). That parent directory is treated as the project root and is
# where pyproject.toml / uv.lock / .venv will live.
#
# e.g. if this script is at   /path/to/project/scripts/rebuild_uv_lock.sh
#      the project root is    /path/to/project
#
# Requires `uv` to be installed:
#   https://docs.astral.sh/uv/getting-started/installation/

set -euo pipefail

# ---------------------------------------------------------------------------
# Resolve project root: one directory above this script's location
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." &> /dev/null && pwd)"

echo "==> Script directory : ${SCRIPT_DIR}"
echo "==> Project root      : ${PROJECT_ROOT}"

cd "${PROJECT_ROOT}"

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

# Third-party (PyPI) packages detected via `import` statements in the
# project's .py files. Standard-library modules (os, sys, json, datetime,
# shutil, warnings, logging, etc.) and the local `pthelperfunctions` module
# are intentionally excluded.
PACKAGES=(
    "pandas"
    "numpy"
    "pyarrow"
    "pytz"
    "matplotlib"
    "clifpy"
)

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------

if ! command -v uv &> /dev/null; then
    echo "Error: 'uv' is not installed or not on PATH." >&2
    echo "Install it first: https://docs.astral.sh/uv/getting-started/installation/" >&2
    exit 1
fi

echo "==> Using $(uv --version)"

# If there's no pyproject.toml yet, initialize one so `uv add` has
# somewhere to record the dependencies.
if [[ ! -f "pyproject.toml" ]]; then
    echo "==> No pyproject.toml found — initializing a new uv project."
    uv init --no-readme --no-pin-python .
fi

# ---------------------------------------------------------------------------
# Remove the existing lock file
# ---------------------------------------------------------------------------

if [[ -f "uv.lock" ]]; then
    echo "==> Removing existing uv.lock"
    rm -f uv.lock
else
    echo "==> No existing uv.lock found — skipping removal"
fi

# ---------------------------------------------------------------------------
# Ensure required packages are declared as dependencies
# ---------------------------------------------------------------------------

echo "==> Adding/ensuring required packages in pyproject.toml:"
printf '    - %s\n' "${PACKAGES[@]}"

# `uv add` is idempotent: if a package is already a dependency it just
# leaves it alone (or updates version constraints if you re-run it).
uv add "${PACKAGES[@]}"

# ---------------------------------------------------------------------------
# Generate a fresh lock file
# ---------------------------------------------------------------------------

echo "==> Generating new uv.lock"
uv lock

echo "==> Done. New uv.lock has been created."
echo "    Run 'uv sync' to install the environment from the new lock file."