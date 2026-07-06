## Configuration

1. Rename `config_template.json` to `config.json`.
2. Fill in your site-specific settings:
   - `site_name` — your site's identifier (used in output file names).
   - `clif_folder` — path to the directory holding your CLIF table files
     (`clif_vitals.parquet`, `clif_labs.parquet`, …).
   - `mimic` — path to MIMIC folder, only used for MIMIC-CLIF "site".
   - `timezone` — your data's timezone, e.g. `"US/Eastern"` (required by clifpy).
   - `filetype` — `"csv"` or `"parquet"`.