"""Shared constants: the common grid, blob layout, analysis windows."""

import numpy as np

PROJECT_PREFIX = "ds-precip-intercomparison"

# Common analysis grid: 0.5 deg, cell edges on multiples of 0.5, centres on x.25/x.75.
RES = 0.5
LAT = np.round(np.arange(-90 + RES / 2, 90, RES), 2)  # ascending, S -> N (360)
LON = np.round(np.arange(-180 + RES / 2, 180, RES), 2)  # ascending, W -> E (720)
NLAT, NLON = LAT.size, LON.size

# Ingest window (first month of every product we keep). Products that start later
# (IMERG 1998) or end earlier (GPCC Full Data 2020, UDel 2017) are stored on their
# own native span; the analysis aligns them.
INGEST_START = "1981-01"

# Blob layout (dev account, `projects` container).
GRID_DIR = f"{PROJECT_PREFIX}/processed/grid05"
AUX_DIR = f"{PROJECT_PREFIX}/processed/aux05"


def grid_blob(product: str) -> str:
    return f"{GRID_DIR}/{product}.nc"


def aux_blob(name: str) -> str:
    return f"{AUX_DIR}/{name}.nc"
