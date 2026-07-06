#!/bin/bash
 
# =============================================================
#  convert_notebooks.sh — Export all .ipynb files to .py
#  Usage: bash convert_notebooks.sh
# =============================================================
 
#!/bin/bash

# =============================================================
#  convert_notebooks.sh — Export all .ipynb files to .py
#  Usage: bash convert_notebooks.sh
# =============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

NOTEBOOKS=(
  "1_cohort.ipynb"
  "2_data_gathering.ipynb"
  "3_calculations.ipynb"
)

echo "📂  Working directory: $SCRIPT_DIR"

for notebook in "${NOTEBOOKS[@]}"; do
  FULL_PATH="$SCRIPT_DIR/$notebook"

  if [ ! -f "$FULL_PATH" ]; then
    echo "⚠️   Not found, skipping: $notebook"
    continue
  fi

  echo "🔄  Converting: $notebook"
  jupyter nbconvert --to script "$FULL_PATH" --output-dir "$SCRIPT_DIR"

  if [ $? -eq 0 ]; then
    echo "✅  Done: ${notebook%.ipynb}.py"
  else
    echo "❌  Failed: $notebook"
  fi
done

echo ""
echo "🎉  Conversion complete."