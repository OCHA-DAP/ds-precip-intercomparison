# ds-precip-intercomparison

Global, land-only, **monthly** comparison of observed precipitation products —
relative biases, temporal agreement, and trends — with a focus on what matters for
anticipatory-action (AA) triggers: gauge-sparse regions, dry-tercile agreement, and
whether apparent trends are real or artefacts of changing station networks.

## Products

| key | product | type | native grid | span | notes |
|---|---|---|---|---|---|
| `chirps_v2` | CHIRPS v2.0 | satellite + stations | 0.05°, 50S–50N | 1981– | used by JRC ASAP within ±50° |
| `chirp_v2` | CHIRP v2.0 | satellite + climatology, **no stations** | 0.05° | 1981– | isolates the station effect in v2 |
| `chirps_v3` | CHIRPS v3.0 | satellite + stations | 0.05°, 60S–60N | 1981– | stations **undercatch-corrected** → wetter than v2 by design; gaps filled with ERA5 |
| `chirp_v3` | CHIRP v3.0 | satellite + climatology, no stations | 0.05° | 1981– | |
| `imerg_final` | IMERG V07 Final, monthly (GPM_3IMERGM) | satellite, gauge-calibrated | 0.1° | 1998– | calibrated to GPCC over land → agreement with GPCC is partly circular |
| `imerg_late` | IMERG V07 Late, daily → monthly | satellite (NRT) | 0.1° | 1998– | the team's operational archive (blob `raster/imerg/daily/late/v7`) |
| `era5` | ERA5 monthly total precipitation | reanalysis (no gauge assimilation) | 0.25° | 1981– | team blob `raster/era5/monthly` |
| `gpcc_full` | GPCC Full Data Monthly v2022 | gauges | 0.5° | 1891–2020 | with `numgauge` per cell |
| `gpcc_monitoring` | GPCC Monitoring Product v2022 | gauges (near-real-time) | 1.0° | 1982– | with `numgauge`, GPCC undercatch `corr_fac`, `solid_frac` |
| `cru` | CRU TS 4.10 | gauges | 0.5° | 1901–2025 | `stn` = stations within correlation distance; relaxes to climatology where none (ODbL — gridded data not redistributed) |
| `cpc` | CPC Global Unified Gauge-Based (daily → monthly) | gauges | 0.5° | 1979– | network change in 2006 (retrospective → real-time) |
| `precl` | NOAA PREC/L | gauges (GHCN + CAMS) | 0.5° | 1948– | |
| `udel` | University of Delaware v5.01 | gauges | 0.5° | 1900–2017 | |

**ASAP blend** (derived in the analysis): CHIRPS v2 within ±50°, ERA5 outside — the
rainfall input of the EU JRC ASAP early-warning system.

Auxiliary layers (`aux05/`): Natural Earth land fraction, Köppen–Geiger 1991–2020
(Beck et al. 2023), Legates–Willmott undercatch factors (CHC distribution), CRU TS
temperature (snow masking), CHIRPS v3 ERA5-fill fraction and station density.

## Pipeline

1. **Ingest (Databricks)** — `scripts/ingest.py <product>` downloads the native data,
   converts to mm/month and regrids to a common 0.5° grid: exact block means for grids
   that nest in 0.5° cells, an area-weighted (0.25, 0.5, 0.25) kernel for ERA5's
   point-centred 0.25° grid, replication for GPCC Monitoring's 1°. A fingerprint check
   (desert / rainforest / Sahel / Zambia / England annual totals and peak months) fails
   the task on a units, orientation or time-offset error. Output: one NetCDF per
   product on the dev blob, `projects/ds-precip-intercomparison/processed/grid05/`.
   Job: `databricks.yml` → `precip_ingest` (unscheduled).
2. **Analysis (local)** — reads the 0.5° cubes and writes summary tables and figures
   for the site.

```bash
uv sync
uv run python scripts/ingest.py era5 --start 2010-01 --end 2010-12 --no-upload   # local test
databricks bundle deploy -t prod -p default && databricks bundle run precip_ingest -t prod -p default
```

Credentials: blob SAS via `DSCI_AZ_BLOB_*` (ocha-stratus); NASA Earthdata via
`IMERG_USERNAME`/`IMERG_PASSWORD` (Databricks Job Compute policy) or `~/.netrc`.
Downloads never retry 401/403 — repeated bad logins lock the Earthdata account,
which the team's daily IMERG pipeline shares.
