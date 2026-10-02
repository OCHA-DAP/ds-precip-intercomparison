"""Auxiliary layers on the common 0.5 deg grid.

Static: land fraction (Natural Earth), Koppen-Geiger 1991-2020 class (Beck et al.
2023), Legates-Willmott monthly undercatch factors (as distributed by CHC).
Monthly: CHIRPS v3 ERA5-fill fraction and CHIRPS v3 station density.
"""

import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import rasterio
import requests
import xarray as xr
from rasterio import features
from rasterio.transform import from_origin

from src.constants import LAT, LON
from src.fetch import download
from src.grid import block_mean, empty_cube, fine_to_05
from src.products import CHC, months, read_fine, to_dataset, trim

NE = "https://naciscdn.org/naturalearth/10m/physical"
KOPPEN_ZIP = "https://ndownloader.figshare.com/files/61012822"


def _flip05(a: np.ndarray) -> np.ndarray:
    """North-up global 0.5 deg raster -> ascending-lat common grid."""
    if a.shape != (360, 720):
        raise ValueError(f"expected (360, 720), got {a.shape}")
    return a[::-1].copy()


def build_static(workdir: Path) -> xr.Dataset:
    import geopandas as gpd

    # Land fraction: Natural Earth 10m land minus 10m lakes, rasterised at 0.02 deg
    # in latitude bands to keep memory down, then block-averaged.
    land = gpd.read_file(download(f"{NE}/ne_10m_land.zip", workdir / "ne_10m_land.zip"))
    lakes = gpd.read_file(download(f"{NE}/ne_10m_lakes.zip", workdir / "ne_10m_lakes.zip"))
    res, factor = 0.02, 25
    nx = int(360 / res)
    band_rows = 450  # 9 deg bands
    rows = []
    for top in np.arange(90, -90, -band_rows * res):
        tr = from_origin(-180, top, res, res)
        m = features.rasterize(
            ((g, 1) for g in land.geometry), out_shape=(band_rows, nx), transform=tr, dtype="uint8"
        )
        lk = features.rasterize(
            ((g, 1) for g in lakes.geometry), out_shape=(band_rows, nx), transform=tr, dtype="uint8"
        )
        m = (m.astype(bool) & ~lk.astype(bool)).astype(np.float32)
        mean, _ = block_mean(m, factor)
        rows.append(mean)
    landfrac = _flip05(np.vstack(rows))

    # Koppen-Geiger 1991-2020 at 0.5 deg
    zp = download(KOPPEN_ZIP, workdir / "koppen_geiger_tif.zip")
    with zipfile.ZipFile(zp) as z:
        name = next(n for n in z.namelist() if n.endswith("1991_2020/koppen_geiger_0p5.tif"))
        z.extract(name, workdir)
    with rasterio.open(workdir / name) as src:
        if not np.isclose(src.transform.a, 0.5) or src.transform.c != -180 or src.transform.f != 90:
            raise ValueError(f"unexpected Koppen grid {src.transform}")
        kg = _flip05(src.read(1).astype(np.float32))

    # Legates-Willmott monthly correction factors (CHC distribution, 0.5 deg, 90N-90S)
    lw = np.full((12, 360, 720), np.nan, dtype=np.float32)
    for m in range(1, 13):
        url = f"{CHC}/CHIRPS/v3.0/diagnostics/legates-willmott_corrections/legates_corrfactor_0.5Deg.90N90S.{m:02d}.tif"
        with rasterio.open(download(url, workdir / Path(url).name)) as src:
            if not (np.isclose(src.transform.a, 0.5) and src.transform.f == 90):
                raise ValueError(f"unexpected LW grid {src.transform}")
            a = src.read(1).astype(np.float32)
            if src.nodata is not None:
                a[a == src.nodata] = np.nan
            a[a <= 0] = np.nan
            lw[m - 1] = _flip05(a)

    ds = xr.Dataset(
        {
            "landfrac": (("lat", "lon"), landfrac, {"long_name": "Natural Earth 10m land (minus lakes) fraction"}),
            "koppen": (("lat", "lon"), kg, {"long_name": "Koppen-Geiger class 1991-2020 (Beck et al. 2023), 0 = water"}),
            "lw_factor": (("month", "lat", "lon"), lw, {"long_name": "Legates-Willmott gauge undercatch factor"}),
        },
        coords={"lat": LAT, "lon": LON, "month": np.arange(1, 13)},
    )
    return ds


def build_chirps_v3_diag(start: str, end: str, workdir: Path, workers: int = 8) -> xr.Dataset:
    """Monthly ERA5-fill fraction (pentads filled / 6) and station count per cell."""
    ts = months(start, end)
    fill = empty_cube(len(ts))
    stn = empty_cube(len(ts))
    base = f"{CHC}/CHIRPS/v3.0/diagnostics"

    def one(i_t):
        i, t = i_t
        out = []
        for url, kind in (
            (f"{base}/fillmaps/monthly/era5-fill.mask.{t.year}.{t.month:02d}.tif", "fill"),
            (f"{base}/global_monthly_station_density/tifs/p25/v3.stn_density.{t.year}.{t.month:02d}.tif", "stn"),
        ):
            try:
                p = download(url, workdir / Path(url).name)
            except requests.HTTPError as exc:
                if exc.response is not None and exc.response.status_code == 404:
                    out.append(None)
                    continue
                raise
            a, top, left, res = read_fine(p)
            p.unlink()
            if kind == "fill":
                mean, _ = fine_to_05(a / 6.0, top, left, res)
                out.append(mean)
            else:
                mean, frac = fine_to_05(np.nan_to_num(a), top, left, res)
                out.append(np.where(np.isfinite(mean), mean * 4, np.nan))  # 4 x 0.25 cells
        return i, out

    with ThreadPoolExecutor(workers) as ex:
        for i, (f, s) in ex.map(one, enumerate(ts)):
            if f is not None:
                fill[i] = f
            if s is not None:
                stn[i] = s
    ds = to_dataset(
        ts,
        {
            "precip": (fill, {"long_name": "fraction of pentads filled with ERA5 (CHIRPS v3)"}),
            "stn": (stn, {"long_name": "CHIRPS v3 stations in cell"}),
        },
        {"source": "CHC CHIRPS v3.0 diagnostics (fillmaps, global_monthly_station_density p25)"},
    )
    ds = trim(ds).rename({"precip": "era5_fill_frac"})
    return ds
