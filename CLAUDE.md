# ds-precip-intercomparison — notes for Claude

Global, land-only, monthly intercomparison of observed precipitation products.
README has the product table and pipeline; this file is the non-obvious stuff.

## Layout
- `src/products.py` — one builder per product → 0.5° monthly cube (mm/month). `src/grid.py`
  regridding (exact block means for nesting grids; ERA5 uses a (0.25, 0.5, 0.25) kernel
  because its 0.25° grid is point-centred). `src/aux.py` — land fraction, Köppen, LW factors,
  CHIRPS v3 diagnostics.
- `scripts/ingest.py` — runs one builder, fingerprint-checks it, uploads to dev blob
  `projects/ds-precip-intercomparison/processed/{grid05,aux05}/`. Runs on Databricks
  (`databricks.yml`, job `precip_ingest`, unscheduled, prod target + dev data plane).
- `scripts/pack.py` → `cache/packed/*.npz` (land cells only; the laptop is disk-constrained).
- `scripts/analyze.py` → `site/report/data/*.json` + `cache/fields.npz`;
  `scripts/figures.py` → `site/report/fig/*.png`; `scripts/explorer_data.py` →
  `site/explorer/data/*.bin` (+ `meta.json`, `cells.bin`; layout in the script docstring).
  Generated site output is NOT committed: `scripts/site_blob.py upload` parks it on dev blob
  (skips unchanged files by sha256), the Pages workflow downloads and verifies it.
- The explorer computes everything client-side (correlation map, unit aggregation) from the
  two loaded product files; adding a statistic is a JS change, adding an aggregation is an
  `explorer_data.py` change + re-upload (~500 MB; Pages limit is 1 GB).

## Gotchas
- **Earthdata lockout**: repeated failed logins lock the NASA account for 10 min. The
  Databricks IMERG credentials are shared with the daily `Run IMERG` raster job — `fetch.py`
  never retries 401/403 and `build_imerg_final` probes auth single-threaded first. Keep it so.
- PSL `downloads.psl.noaa.gov` 502s; use the THREDDS `fileServer` URLs. THREDDS stalls
  mid-file — `download()` resumes only on a 206 reply.
- GPCC Monitoring `solid_p` is a **percentage**, not mm. `s` = gauges in the 1° cell.
- ERA5 COGs on blob are **mm/day** (monthly-mean rate) → × days in month.
- IMERG Final is calibrated to GPCC (Full Data to 2020, Monitoring from 2021) — agreement with
  GPCC is partly circular, and 2021 is a calibration break.
- CHIRPS v3 stations are undercatch-corrected (wetter than v2 by design); gaps filled with ERA5.
- CRU TS is ODbL: publish derived statistics only, never the regridded field.
- Reference choice: GPCC Full Data, stratified by its gauge count. No ensemble-median
  reference — most products share GPCC/GHCN stations, so a median is anchored to GPCC anyway.
