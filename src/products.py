"""One builder per product: native source -> monthly totals (mm/month) on the
common 0.5 deg grid.

Every builder returns an ``xr.Dataset`` with dims (time, lat, lon) — ``time``
is the first of each month — holding at least ``precip`` [mm/month], plus any
product-specific diagnostics (gauge counts, valid fraction, ...). Builders do
not apply a land mask; the analysis does, from ``aux05/landfrac``.
"""

import calendar
import re
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import requests
import xarray as xr

from src.constants import LAT, LON
from src.fetch import download, earthdata_session
from src.grid import (
    assert_on_grid,
    block_mean,
    empty_cube,
    era5_to_05,
    fine_to_05,
    normalise,
    replicate,
)

CHC = "https://data.chc.ucsb.edu/products"
PSL = "https://psl.noaa.gov/thredds/fileServer/Datasets"
DWD = "https://opendata.dwd.de/climate_environment/GPCC"
CRU = "https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.10/cruts.2604091129.v4.10"
GES = "https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGM.07"


def months(start: str, end: str) -> pd.DatetimeIndex:
    return pd.date_range(start, end, freq="MS")


def ndays(t: pd.Timestamp) -> int:
    return calendar.monthrange(t.year, t.month)[1]


def to_dataset(times, variables: dict, attrs: dict) -> xr.Dataset:
    data_vars = {}
    for name, (arr, vattrs) in variables.items():
        dims = ("time", "lat", "lon") if arr.ndim == 3 else ("lat", "lon")
        data_vars[name] = (dims, arr, vattrs)
    ds = xr.Dataset(data_vars, coords={"time": pd.DatetimeIndex(times), "lat": LAT, "lon": LON})
    ds.attrs.update(attrs)
    return ds


def trim(ds: xr.Dataset) -> xr.Dataset:
    """Drop all-missing months at the end (product not yet published)."""
    ok = ds["precip"].notnull().any(("lat", "lon")).values
    if not ok.any():
        raise RuntimeError("no data at all")
    last = np.where(ok)[0][-1]
    return ds.isel(time=slice(0, last + 1))


# ---------------------------------------------------------------- CHC (0.05 deg)

CHC_URLS = {
    "chirps_v2": CHC + "/CHIRPS-2.0/global_monthly/cogs/chirps-v2.0.{y}.{m:02d}.cog",
    "chirp_v2": CHC + "/CHIRP/monthly/CHIRP.{y}.{m:02d}.tif",
    "chirps_v3": CHC + "/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.{y}.{m:02d}.tif",
    "chirp_v3": CHC + "/CHIRP-v3.0/monthly/global/tifs/chirp-v3.0.{y}.{m:02d}.tif",
}


# shape -> (top, left, res) for CHC rasters that arrive without a geotransform
UNGEOREFERENCED = {(2400, 7200): (60.0, -180.0, 0.05), (2000, 7200): (50.0, -180.0, 0.05)}


def read_fine(path: Path, scale: float = 1.0) -> tuple[np.ndarray, float, float, float]:
    """North-up float32 array with nodata/negatives as NaN, plus top/left/res."""
    with rasterio.open(path) as src:
        a = src.read(1).astype(np.float32)
        nod = src.nodata
        t = src.transform
        if t.is_identity and a.shape in UNGEOREFERENCED:
            # CHC has published some files without georeferencing (e.g. CHIRP v3
            # from 2024-10). Same pixel layout as the georeferenced months
            # (verified: r = 0.96 with CHIRPS v3 2024-10 as-is, 0.20 flipped).
            top, left, res = UNGEOREFERENCED[a.shape]
            print(f"[read_fine] {path.name}: no geotransform, assuming top={top} left={left} res={res}", flush=True)
        else:
            res = round(t.a, 4)
            top, left = round(t.f, 3), round(t.c, 3)
            if not np.isclose(-t.e, t.a, rtol=1e-4):
                raise ValueError(f"{path}: non-square pixels {t}")
    if nod is not None:
        a[a == nod] = np.nan
    a[a < 0] = np.nan  # -9999 / -99 style fills
    return a * scale, top, left, res


def build_chc(product: str, start: str, end: str, workdir: Path, workers: int = 6) -> xr.Dataset:
    tpl = CHC_URLS[product]
    ts = months(start, end)
    cube = empty_cube(len(ts))
    frac_max = np.zeros(cube.shape[1:], dtype=np.float32)

    def one(i_t):
        i, t = i_t
        url = tpl.format(y=t.year, m=t.month)
        dest = workdir / Path(url).name
        try:
            download(url, dest)
        except requests.HTTPError as exc:
            if exc.response is not None and exc.response.status_code == 404:
                return i, None, None
            raise
        a, top, left, res = read_fine(dest)
        dest.unlink()
        mean, frac = fine_to_05(a, top, left, res)
        return i, mean, frac

    with ThreadPoolExecutor(workers) as ex:
        for i, mean, frac in ex.map(one, enumerate(ts)):
            if mean is not None:
                cube[i] = mean
                np.fmax(frac_max, frac, out=frac_max)
            if i % 60 == 0:
                print(f"[{product}] {ts[i]:%Y-%m} done", flush=True)
    ds = to_dataset(
        ts,
        {
            "precip": (cube, {"units": "mm/month", "long_name": "monthly precipitation total"}),
            "valid_frac": (frac_max, {"long_name": "max fraction of 0.05 deg pixels with data"}),
        },
        {"source_url_template": tpl},
    )
    return trim(ds)


# ---------------------------------------------------------------- IMERG

def build_imerg_final(start: str, end: str, workdir: Path, workers: int = 4) -> xr.Dataset:
    """GPM_3IMERGM V07 Final monthly (mm/hr monthly-mean rate -> mm/month)."""
    import h5py

    sess = earthdata_session()
    ts = months(max(pd.Timestamp(start), pd.Timestamp("1998-01-01")).strftime("%Y-%m"), end)
    names: dict[tuple[int, int], str] = {}
    for y in sorted({t.year for t in ts}):
        r = sess.get(f"{GES}/{y}/", timeout=120)
        if r.status_code == 404:
            continue
        r.raise_for_status()
        for fn in set(re.findall(r"3B-MO\.MS\.MRG\.3IMERG\.\d{8}-S000000-E235959\.\d{2}\.V07[A-Z]\.HDF5(?!\.)", r.text)):
            ym = re.search(r"3IMERG\.(\d{4})(\d{2})01", fn)
            key = (int(ym.group(1)), int(ym.group(2)))
            # keep the latest version letter if a month has several
            if key not in names or fn > names[key]:
                names[key] = fn
    if not names:
        raise RuntimeError("no IMERG monthly files listed (directory listing failed?)")
    # One authenticated request before going parallel: a credential problem must
    # fail once, not N threads x retries (that locks the Earthdata account).
    first_key = min(names)
    probe = download(f"{GES}/{first_key[0]}/{names[first_key]}", workdir / names[first_key], session=sess, retries=1)
    print(f"[imerg_final] auth OK ({probe.stat().st_size / 1e6:.0f} MB probe)", flush=True)
    cube = empty_cube(len(ts))
    grw = empty_cube(len(ts))
    versions = [""] * len(ts)

    def one(i_t):
        i, t = i_t
        fn = names.get((t.year, t.month))
        if fn is None:
            return i, None, None, ""
        dest = download(f"{GES}/{t.year}/{fn}", workdir / fn, session=sess)
        with h5py.File(dest, "r") as f:
            p = f["Grid/precipitation"][0].astype(np.float32)  # (lon 3600, lat 1800)
            g = f["Grid/gaugeRelativeWeighting"][0].astype(np.float32)
            lat = f["Grid/lat"][:]
            lon = f["Grid/lon"][:]
        dest.unlink()
        if not (lat[0] < lat[-1] and np.isclose(lat[0], -89.95) and np.isclose(lon[0], -179.95)):
            raise ValueError(f"unexpected IMERG grid lat0={lat[0]} lon0={lon[0]}")
        p, g = p.T, g.T  # -> (lat ascending, lon)
        p[p < 0] = np.nan
        g[g < 0] = np.nan
        mean, _ = block_mean(p, 5)
        gm, _ = block_mean(g, 5)
        return i, mean * 24 * ndays(t), gm, fn.split(".")[-2]

    with ThreadPoolExecutor(workers) as ex:
        for i, p, g, v in ex.map(one, enumerate(ts)):
            if p is not None:
                cube[i], grw[i], versions[i] = p, g, v
            if i % 60 == 0:
                print(f"[imerg_final] {ts[i]:%Y-%m} done", flush=True)
    ds = to_dataset(
        ts,
        {
            "precip": (cube, {"units": "mm/month"}),
            "gauge_weight": (grw, {"units": "%", "long_name": "gaugeRelativeWeighting (block mean)"}),
        },
        {"source": "NASA GES DISC GPM_3IMERGM V07 (Final Run monthly)"},
    )
    ds["version"] = ("time", np.array(versions, dtype="U8"))
    return trim(ds)


def build_imerg_late(start: str, end: str, workdir: Path, workers: int = 8) -> xr.Dataset:
    """Monthly totals from the team's daily IMERG V07 Late COGs (prod raster blob)."""
    import ocha_stratus as stratus

    blobs = set(stratus.list_container_blobs("imerg/daily/late/v7/processed/", stage="prod", container_name="raster"))
    ts = months(max(pd.Timestamp(start), pd.Timestamp("1998-01-01")).strftime("%Y-%m"), end)
    cube = empty_cube(len(ts))
    ndays_found = np.zeros(len(ts), dtype=np.int16)
    cc = stratus.get_container_client(container_name="raster", stage="prod")

    def read_day(name):
        da = stratus.open_blob_cog(name, stage="prod", container_name="raster", container_client=cc)
        a = da.squeeze().values.astype(np.float32)
        y = da.y.values
        if not (y[0] > y[-1] and np.isclose(y[0], 89.95) and np.isclose(float(da.x[0]), -179.95)):
            raise ValueError(f"unexpected IMERG Late grid in {name}")
        a[a < 0] = np.nan
        return a

    for i, t in enumerate(ts):
        names = [
            f"imerg/daily/late/v7/processed/imerg-daily-late-{t.year}-{t.month:02d}-{d:02d}.tif"
            for d in range(1, ndays(t) + 1)
        ]
        have = [n for n in names if n in blobs]
        ndays_found[i] = len(have)
        if len(have) < len(names):
            continue  # incomplete month -> NaN
        with ThreadPoolExecutor(workers) as ex:
            tot = None
            for a in ex.map(read_day, have):
                tot = a if tot is None else tot + a  # NaN in any day -> NaN pixel
        mean, _ = block_mean(tot[::-1], 5)  # north-up -> ascending lat
        cube[i] = mean
        if i % 24 == 0:
            print(f"[imerg_late] {t:%Y-%m} done", flush=True)
    ds = to_dataset(ts, {"precip": (cube, {"units": "mm/month"})},
                    {"source": "team blob raster/imerg/daily/late/v7/processed (daily sums)"})
    ds["ndays"] = ("time", ndays_found)
    return trim(ds)


# ---------------------------------------------------------------- ERA5

def build_era5(start: str, end: str, workdir: Path, workers: int = 8) -> xr.Dataset:
    import ocha_stratus as stratus

    blobs = set(stratus.list_container_blobs("era5/monthly/processed/", stage="prod", container_name="raster"))
    ts = months(start, end)
    cube = empty_cube(len(ts))
    cc = stratus.get_container_client(container_name="raster", stage="prod")

    def one(i_t):
        i, t = i_t
        name = f"era5/monthly/processed/precip_reanalysis_v{t:%Y-%m}-01.tif"
        if name not in blobs:
            return i, None
        da = stratus.open_blob_cog(name, stage="prod", container_name="raster", container_client=cc).squeeze()
        if da.attrs.get("units") != "mm/day":
            raise ValueError(f"{name}: units {da.attrs.get('units')!r}, expected mm/day")
        out = era5_to_05(da.values, da.y.values, da.x.values)
        return i, out * ndays(t)

    with ThreadPoolExecutor(workers) as ex:
        for i, a in ex.map(one, enumerate(ts)):
            if a is not None:
                cube[i] = a
    ds = to_dataset(ts, {"precip": (cube, {"units": "mm/month"})},
                    {"source": "team blob raster/era5/monthly/processed (CDS reanalysis-era5-single-levels-monthly-means tp)"})
    return trim(ds)


# ---------------------------------------------------------------- gauge-only

def _slice(ds: xr.Dataset, start: str, end: str) -> xr.Dataset:
    ds = ds.sel(time=slice(start, end))
    ds = ds.assign_coords(time=pd.DatetimeIndex(ds.time.values).to_period("M").to_timestamp())
    return ds


def build_gpcc_full(start: str, end: str, workdir: Path, **_) -> xr.Dataset:
    # Per variable, float32 throughout: a to_array()/sortby() of the 4-decade
    # concat would peak at ~5 GB on a shared 14 GB driver.
    parts = {"precip": [], "numgauge": []}
    for d in range(1981, 2021, 10):
        url = f"{DWD}/full_data_monthly_v2022/05/full_data_monthly_v2022_{d}_{d + 9}_05.nc.gz"
        dest = download(url, workdir / f"gpcc_full_{d}.nc", gunzip=True)
        with xr.open_dataset(dest) as src:
            for v in parts:
                parts[v].append(normalise(src[v].astype("float32").load()))
        dest.unlink()
    ds = xr.Dataset({v: xr.concat(parts[v], "time") for v in parts})
    assert_on_grid(ds, "gpcc_full")
    ds = _slice(ds, start, end)
    ds["precip"].attrs["units"] = "mm/month"
    ds.attrs["source"] = "DWD GPCC Full Data Monthly v2022 0.5 deg"
    return ds


def build_gpcc_monitoring(start: str, end: str, workdir: Path, workers: int = 6) -> xr.Dataset:
    ts = months(max(pd.Timestamp(start), pd.Timestamp("1982-01-01")).strftime("%Y-%m"), end)
    out = {k: empty_cube(len(ts)) for k in ("precip", "numgauge", "corr_fac", "solid_frac")}

    def one(i_t):
        i, t = i_t
        url = f"{DWD}/monitoring_v2022/{t.year}/monitoring_v2022_10_{t.year}_{t.month:02d}.nc.gz"
        try:
            dest = download(url, workdir / f"gpcc_mon_{t:%Y%m}.nc", gunzip=True)
        except requests.HTTPError as exc:
            if exc.response is not None and exc.response.status_code == 404:
                return i, None
            raise
        ds = normalise(xr.open_dataset(dest).isel(time=0)[["p", "s", "corr_fac", "solid_p"]].to_array("v").load())
        dest.unlink()
        if ds.lat.size != 180 or not np.isclose(float(ds.lat[0]), -89.5):
            raise ValueError("unexpected GPCC monitoring grid")
        res = {
            "precip": ds.sel(v="p").values,
            "numgauge": ds.sel(v="s").values,
            "corr_fac": ds.sel(v="corr_fac").values,
            "solid_frac": ds.sel(v="solid_p").values / 100.0,  # solid_p is a percentage
        }
        return i, {k: replicate(v.astype(np.float32), 2) for k, v in res.items()}

    with ThreadPoolExecutor(workers) as ex:
        for i, res in ex.map(one, enumerate(ts)):
            if res is not None:
                for k, v in res.items():
                    out[k][i] = v
    ds = to_dataset(
        ts,
        {
            "precip": (out["precip"], {"units": "mm/month"}),
            "numgauge": (out["numgauge"], {"long_name": "gauges in the parent 1 deg cell"}),
            "corr_fac": (out["corr_fac"], {"long_name": "GPCC systematic gauge-error (undercatch) correction factor"}),
            "solid_frac": (out["solid_frac"], {"long_name": "solid / total precipitation (GPCC)"}),
        },
        {"source": "DWD GPCC Monitoring Product v2022 1.0 deg, replicated to 0.5 deg"},
    )
    return trim(ds)


def build_cru(start: str, end: str, workdir: Path, **_) -> xr.Dataset:
    url = f"{CRU}/pre/cru_ts4.10.1901.2025.pre.dat.nc.gz"
    dest = download(url, workdir / "cru_pre.nc", gunzip=True)
    with xr.open_dataset(dest) as src:
        ds = xr.Dataset({
            "precip": normalise(src["pre"].sel(time=slice(start, end)).astype("float32").load()),
            "stn": normalise(src["stn"].sel(time=slice(start, end)).astype("float32").load()),
        })
    dest.unlink()
    assert_on_grid(ds, "cru")
    ds = _slice(ds, start, end)
    ds["precip"].attrs["units"] = "mm/month"
    ds["stn"].attrs["long_name"] = "stations within correlation-decay distance (CRU)"
    ds.attrs["source"] = "CRU TS 4.10 (UEA), ODbL - do not redistribute the gridded data"
    return ds


def build_cru_tmp(start: str, end: str, workdir: Path, **_) -> xr.Dataset:
    url = f"{CRU}/tmp/cru_ts4.10.1901.2025.tmp.dat.nc.gz"
    dest = download(url, workdir / "cru_tmp.nc", gunzip=True)
    with xr.open_dataset(dest) as src:
        ds = xr.Dataset({"tmp": normalise(src["tmp"].sel(time=slice(start, end)).astype("float32").load())})
    dest.unlink()
    assert_on_grid(ds, "cru_tmp")
    ds = _slice(ds, start, end)
    ds.attrs["source"] = "CRU TS 4.10 tmp (UEA), ODbL"
    return ds


def build_cpc(start: str, end: str, workdir: Path, **_) -> xr.Dataset:
    y0, y1 = pd.Timestamp(start).year, pd.Timestamp(end).year
    parts = []
    for y in range(y0, y1 + 1):
        try:
            dest = download(f"{PSL}/cpc_global_precip/precip.{y}.nc", workdir / f"cpc_{y}.nc")
        except requests.HTTPError as exc:
            if exc.response is not None and exc.response.status_code == 404:
                break
            raise
        da = xr.open_dataset(dest)["precip"].astype("float32").load()
        # monthly sum only where every day of the month is present
        cnt = da.notnull().resample(time="MS").sum()
        print(f"[cpc] {y}: max valid days per month {cnt.max(('lat', 'lon')).values.tolist()}", flush=True)
        tot = da.resample(time="MS").sum(min_count=1)
        dim = xr.DataArray([ndays(pd.Timestamp(t)) for t in tot.time.values], coords={"time": tot.time})
        parts.append(tot.where(cnt == dim))
        dest.unlink()
        print(f"[cpc] {y} done", flush=True)
    da = normalise(xr.concat(parts, "time"))
    assert_on_grid(da, "cpc")
    ds = _slice(da.to_dataset(name="precip"), start, end)
    ds["precip"].attrs["units"] = "mm/month"
    ds.attrs["source"] = "NOAA CPC Global Unified Gauge-Based Analysis of Daily Precipitation (PSL), monthly sums"
    return trim(ds.astype("float32"))


def build_precl(start: str, end: str, workdir: Path, **_) -> xr.Dataset:
    dest = download(f"{PSL}/precl/0.5deg/precip.mon.mean.0.5x0.5.nc", workdir / "precl.nc")
    da = xr.open_dataset(dest)["precip"].sel(time=slice(start, end)).load()
    units = da.attrs.get("units", "")
    if units not in ("mm/day", "mm/day "):
        raise ValueError(f"PREC/L units {units!r}")
    dim = xr.DataArray([ndays(pd.Timestamp(t)) for t in da.time.values], coords={"time": da.time})
    da = normalise(da * dim)
    assert_on_grid(da, "precl")
    ds = _slice(da.to_dataset(name="precip"), start, end)
    ds["precip"].attrs["units"] = "mm/month"
    ds.attrs["source"] = "NOAA PREC/L 0.5 deg (PSL)"
    return trim(ds.astype("float32"))


def build_udel(start: str, end: str, workdir: Path, **_) -> xr.Dataset:
    dest = download(f"{PSL}/udel.airt.precip/precip.mon.total.v501.nc", workdir / "udel.nc")
    da = xr.open_dataset(dest)["precip"].sel(time=slice(start, end)).load()
    units = da.attrs.get("units", "")
    if not units.lower().startswith("cm"):
        raise ValueError(f"UDel units {units!r}, expected cm")
    da = normalise(da * 10.0)
    assert_on_grid(da, "udel")
    ds = _slice(da.to_dataset(name="precip"), start, end)
    ds["precip"].attrs["units"] = "mm/month"
    ds.attrs["source"] = "University of Delaware v5.01 (PSL)"
    return ds.astype("float32")


BUILDERS = {
    "chirps_v2": lambda s, e, w, **k: build_chc("chirps_v2", s, e, w, **k),
    "chirp_v2": lambda s, e, w, **k: build_chc("chirp_v2", s, e, w, **k),
    "chirps_v3": lambda s, e, w, **k: build_chc("chirps_v3", s, e, w, **k),
    "chirp_v3": lambda s, e, w, **k: build_chc("chirp_v3", s, e, w, **k),
    "imerg_final": build_imerg_final,
    "imerg_late": build_imerg_late,
    "era5": build_era5,
    "gpcc_full": build_gpcc_full,
    "gpcc_monitoring": build_gpcc_monitoring,
    "cru": build_cru,
    "cru_tmp": build_cru_tmp,
    "cpc": build_cpc,
    "precl": build_precl,
    "udel": build_udel,
}


def default_workdir() -> Path:
    base = Path("/local_disk0") if Path("/local_disk0").is_dir() else Path(tempfile.gettempdir())
    return base / "precip_ingest"
