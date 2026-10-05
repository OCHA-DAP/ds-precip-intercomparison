"""Data for the interactive explorer (site/explorer/): per-pixel yearly totals.

    python scripts/explorer_data.py      # needs cache/packed (scripts/pack.py)

For every product and aggregation — calendar year, rainy season, and the twelve
3-month windows — writes one binary file holding each land cell's series of
yearly totals, so the browser can correlate any two products per pixel and
aggregate to countries / admin-1 units itself:

    site/explorer/data/<product>_<agg>.bin
        uint8  q[ncell * nyear]   cell-major; 255 = missing
        uint16 lo[ncell]          mm, per-cell minimum over years
        uint16 hi[ncell]          mm, per-cell maximum over years
        value = lo + q / 254 * (hi - lo)     (error <= (hi - lo) / 508)
    site/explorer/data/cells.bin  + meta.json (layout, years, names)

Rainy season: the calendar months whose GPCC 2001-2020 climatology is at least
the monthly mean (1/12 of the annual total), summed within a hydrological year
that starts the month after the climatologically driest month — so a Nov-Mar
season is not split across two years. Each season is labelled by the calendar
year of its last rainy month (a Jun-Sep season by that year, a Nov-Mar season
by the year of the March).
3-month windows crossing the year end (NDJ, DJF) are labelled by the year of
their last month.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import metrics as M  # noqa: E402
from src.adata import LABEL, ROOT, countries, load, products_present  # noqa: E402
from src.fetch import download  # noqa: E402

OUT = ROOT / "site" / "explorer" / "data"
CACHE = ROOT / "cache"
MONTHS = "JFMAMJJASOND"
AGGS = ["annual", "rainy"] + ["".join(MONTHS[(m + k) % 12] for k in range(3)) for m in range(12)]


def admin1(cells: dict) -> pd.DataFrame:
    """Majority Natural Earth admin-1 unit per 0.5 deg land cell."""
    import geopandas as gpd
    import rasterio.features as rf
    from rasterio.transform import from_origin

    out = CACHE / "cell_admin1.parquet"
    if out.exists():
        return pd.read_parquet(out)
    zp = download(
        "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_1_states_provinces.zip",
        CACHE / "ne_10m_admin_1_states_provinces.zip",
    )
    gdf = gpd.read_file(zp)[["name", "admin", "adm0_a3", "geometry"]].reset_index(drop=True)
    idx = rf.rasterize(((g, i + 1) for i, g in enumerate(gdf.geometry)), out_shape=(3600, 7200),
                       transform=from_origin(-180, 90, 0.05, 0.05), dtype="int32")[::-1]
    blocks = idx.reshape(360, 10, 720, 10).transpose(0, 2, 1, 3).reshape(360, 720, 100)
    b = blocks[cells["iy"], cells["ix"]]
    code = np.zeros(b.shape[0], dtype=np.int32)
    for k in range(b.shape[0]):
        row = b[k][b[k] > 0]
        if row.size:
            code[k] = np.bincount(row).argmax()
    df = pd.DataFrame({
        "admin1": np.where(code > 0, gdf.name.reindex(code - 1).fillna("").values, ""),
        "admin1_country": np.where(code > 0, gdf.admin.reindex(code - 1).fillna("").values, ""),
        "admin1_iso3": np.where(code > 0, gdf.adm0_a3.reindex(code - 1).fillna("").values, ""),
    })
    df.to_parquet(out)
    return df


def by_year(x: np.ndarray, time: pd.DatetimeIndex) -> tuple[np.ndarray, np.ndarray]:
    """(T, N) monthly -> (Y, 12, N) on whole calendar years, NaN-padded."""
    y0, y1 = time[0].year, time[-1].year
    years = np.arange(y0, y1 + 1)
    out = np.full((len(years) * 12, x.shape[1]), np.nan, dtype=np.float32)
    idx = (time.year.values - y0) * 12 + time.month.values - 1
    out[idx] = x
    return out.reshape(len(years), 12, -1), years


def aggregate(X: np.ndarray, years: np.ndarray, start: np.ndarray, rainy: np.ndarray) -> dict:
    """All aggregations -> {agg: (Y, N)} labelled by `years`; NaN unless every
    month of the window (or of the rainy season) is present."""
    Y, _, N = X.shape
    nxt = np.concatenate([X[1:], np.full((1, 12, N), np.nan, np.float32)], 0)
    ext = np.concatenate([X, nxt], 1)  # (Y, 24, N): this year + next year
    res = {}
    res["annual"] = np.where(np.isfinite(X).all(1), np.nansum(X, 1), np.nan)
    for m in range(12):
        w = ext[:, m:m + 3]
        tot = np.where(np.isfinite(w).all(1), np.nansum(w, 1), np.nan)
        if m + 2 >= 12:  # crosses the year end: label by the year of the last month
            tot = np.concatenate([np.full((1, N), np.nan, np.float32), tot[:-1]], 0)
        res[AGGS[2 + m]] = tot
    # rainy season over the hydrological year starting at `start` (0-based month)
    rain = np.full((Y, N), np.nan, dtype=np.float32)
    for s in range(12):
        cols = np.nonzero(start == s)[0]
        if cols.size == 0:
            continue
        mask = np.concatenate([rainy, rainy], 0)[s:s + 12][:, cols]  # (12, n) months in hydro order
        seg = ext[:, s:s + 12][:, :, cols]  # hydro year starting in year y
        ok = np.where(mask[None], np.isfinite(seg), True).all(1)
        tot = np.where(ok, np.nansum(np.where(mask[None], seg, 0), 1), np.nan)
        # Label by the calendar year the season ENDS in: only seasons with rainy
        # months before `s` (i.e. in the next calendar year) shift to y+1. A season
        # lying inside one calendar year keeps label y even when s > 0.
        crosses = rainy[:s, cols].any(0) if s > 0 else np.zeros(cols.size, bool)
        shifted = np.concatenate([np.full((1, cols.size), np.nan, np.float32), tot[:-1]], 0)
        rain[:, cols] = np.where(crosses[None, :], shifted, tot)
    res["rainy"] = rain
    return res


def quantise(A: np.ndarray) -> tuple[bytes, int, int]:
    """(Y, N) -> bytes of the .bin layout; trims leading/trailing all-NaN years."""
    ok_years = np.nonzero(np.isfinite(A).any(1))[0]
    y_first, y_last = ok_years[0], ok_years[-1]
    A = A[y_first:y_last + 1]
    with np.errstate(all="ignore"):
        lo = np.floor(np.nanmin(A, 0))
        hi = np.ceil(np.nanmax(A, 0))
    lo = np.where(np.isfinite(lo), np.clip(lo, 0, 65535), 0).astype(np.uint16)
    hi = np.where(np.isfinite(hi), np.clip(hi, 0, 65535), 0).astype(np.uint16)
    span = np.maximum(hi.astype(np.float32) - lo, 1e-6)
    q = np.where(np.isfinite(A), np.clip(np.rint((A - lo) / span * 254), 0, 254), 255).astype(np.uint8)
    blob = q.T.copy().tobytes() + lo.tobytes() + hi.tobytes()  # cell-major
    return blob, int(y_first), int(y_last)


def main() -> None:
    d = load()
    P = [p for p in products_present(d) if p != "asap"]
    cells = d.cells
    N = cells["lat"].size
    OUT.mkdir(parents=True, exist_ok=True)

    # rainy-season definition from GPCC 2001-2020 (ERA5 where GPCC has no climatology)
    tb = d.window(2001, 2020)
    clim = M.monthly_clim(d.precip["gpcc_full"][tb], d.months[tb])
    alt = M.monthly_clim(d.precip["era5"][tb], d.months[tb])
    clim = np.where(np.isfinite(clim).all(0)[None], clim, alt)
    rainy = clim >= np.nansum(clim, 0)[None] / 12
    start = (np.nanargmin(np.where(np.isfinite(clim), clim, np.inf), 0) + 1) % 12

    ctry = countries(cells, CACHE)
    a1 = admin1(cells)
    c_names = sorted({(i, n) for i, n in zip(ctry.iso3, ctry.country) if i})
    c_index = {iso: k for k, (iso, _) in enumerate(c_names)}
    a1_keys = sorted({(i, n) for i, n in zip(a1.admin1_iso3, a1.admin1) if n})
    a1_index = {key: k for k, key in enumerate(a1_keys)}
    cbin = (
        cells["iy"].astype(np.int16).tobytes()
        + cells["ix"].astype(np.int16).tobytes()
        + np.array([c_index.get(i, 65535) for i in ctry.iso3], dtype=np.uint16).tobytes()
        + np.array([a1_index.get((i, n), 65535) if n else 65535 for i, n in zip(a1.admin1_iso3, a1.admin1)],
                   dtype=np.uint16).tobytes()
        + (cells["area"] * cells["landfrac"]).astype(np.float32).tobytes()
        + (rainy * (1 << np.arange(12))[:, None]).sum(0).astype(np.uint16).tobytes()
        + start.astype(np.uint8).tobytes()
    )
    (OUT / "cells.bin").write_bytes(cbin)

    files = {}
    for p in P:
        X, years = by_year(d.precip[p], d.time)
        aggs = aggregate(X, years, start, rainy)
        for a, A in aggs.items():
            blob, y0, y1 = quantise(A)
            (OUT / f"{p}_{a}.bin").write_bytes(blob)
            files[f"{p}_{a}"] = {"first_year": int(years[y0]), "nyear": int(y1 - y0 + 1)}
        print(f"[explorer] {p}: {len(aggs)} aggregations", flush=True)

    meta = {
        "ncell": int(N),
        "cells_layout": ["iy:int16", "ix:int16", "country:uint16", "admin1:uint16", "weight:float32",
                         "rainy_mask:uint16", "hydro_start:uint8"],
        "products": [{"key": p, "label": LABEL[p]} for p in P],
        "aggs": [{"key": "annual", "label": "Calendar year"},
                 {"key": "rainy", "label": "Rainy season (hydrological year)"}]
                + [{"key": a, "label": f"{a} (3 months)"} for a in AGGS[2:]],
        "files": files,
        "countries": [{"iso3": i, "name": n} for i, n in c_names],
        "admin1": [{"iso3": i, "name": n} for i, n in a1_keys],
        "masked": {"imerg_late": "Oct-Nov 2024 (NASA calibration incident) excluded"},
    }
    (OUT / "meta.json").write_text(json.dumps(meta, separators=(",", ":")))
    size = sum(f.stat().st_size for f in OUT.glob("*")) / 1e6
    print(f"[explorer] wrote {len(files)} series files + cells.bin + meta.json ({size:.0f} MB)")


if __name__ == "__main__":
    main()
