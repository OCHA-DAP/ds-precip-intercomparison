"""Analysis-side data: packed land cells aligned on one monthly axis, plus
the region definitions (country, continent, Koppen, gauge density)."""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PACKED = ROOT / "cache" / "packed"

# key, label, family, colour slot (family-level identity; see site legend)
PRODUCTS = [
    ("gpcc_full", "GPCC Full Data v2022", "gauge"),
    ("gpcc_monitoring", "GPCC Monitoring v2022", "gauge"),
    ("cru", "CRU TS 4.10", "gauge"),
    ("cpc", "CPC Unified", "gauge"),
    ("precl", "PREC/L", "gauge"),
    ("udel", "UDel v5.01", "gauge"),
    ("chirps_v2", "CHIRPS v2.0", "satellite+gauge"),
    ("chirp_v2", "CHIRP v2.0 (no stations)", "satellite"),
    ("chirps_v3", "CHIRPS v3.0", "satellite+gauge"),
    ("chirp_v3", "CHIRP v3.0 (no stations)", "satellite"),
    ("imerg_final", "IMERG V07 Final", "satellite+gauge"),
    ("imerg_late", "IMERG V07 Late", "satellite"),
    ("era5", "ERA5", "reanalysis"),
    ("asap", "ASAP blend (CHIRPS v2 ±50°, ERA5 beyond)", "composite"),
]
LABEL = {k: lab for k, lab, _ in PRODUCTS}
FAMILY = {k: fam for k, _, fam in PRODUCTS}
GAUGE = [k for k, _, f in PRODUCTS if f == "gauge"]

KOPPEN_MAJOR = {
    "A Tropical": range(1, 4),
    "B Arid": range(4, 8),
    "C Temperate": range(8, 17),
    "D Cold": range(17, 29),
    "E Polar": range(29, 31),
}


@dataclass
class Data:
    time: pd.DatetimeIndex
    cells: dict
    precip: dict = field(default_factory=dict)  # product -> (T, N) float32
    extra: dict = field(default_factory=dict)  # name -> (T, N) or (N,)

    @property
    def years(self) -> np.ndarray:
        return self.time.year.values

    @property
    def months(self) -> np.ndarray:
        return self.time.month.values

    def window(self, y0: int, y1: int) -> np.ndarray:
        return (self.years >= y0) & (self.years <= y1)

    def span(self, p: str) -> tuple[str, str]:
        ok = np.isfinite(self.precip[p]).any(1)
        idx = np.nonzero(ok)[0]
        return f"{self.time[idx[0]]:%Y-%m}", f"{self.time[idx[-1]]:%Y-%m}"


def _align(t_src: np.ndarray, arr: np.ndarray, time: pd.DatetimeIndex) -> np.ndarray:
    out = np.full((len(time),) + arr.shape[1:], np.nan, dtype=np.float32)
    pos = time.get_indexer(pd.DatetimeIndex(pd.to_datetime(t_src)))
    ok = pos >= 0
    out[pos[ok]] = arr[ok]
    return out


def load(end: str | None = None) -> Data:
    cells = dict(np.load(PACKED / "cells.npz"))
    files = {p.stem: p for p in PACKED.glob("*.npz") if p.stem != "cells"}
    ends = []
    raw = {}
    for name, path in files.items():
        z = dict(np.load(path, allow_pickle=False))
        raw[name] = z
        ends.append(pd.Timestamp(str(z["time"][-1])))
    time = pd.date_range("1981-01-01", end or max(ends), freq="MS")
    d = Data(time=time, cells=cells)
    extras = {
        "gpcc_full": {"numgauge": "gpcc_numgauge"},
        "gpcc_monitoring": {"numgauge": "gpccmon_numgauge", "corr_fac": "gpcc_corr_fac", "solid_frac": "gpcc_solid_frac"},
        "cru": {"stn": "cru_stn"},
        "imerg_final": {"gauge_weight": "imerg_gauge_weight"},
        "chirps_v3_diag": {"era5_fill_frac": "chirps3_fill", "stn": "chirps3_stn"},
        "cru_tmp": {"tmp": "cru_tmp"},
    }
    for name, z in raw.items():
        t = z["time"]
        if "precip" in z and name not in ("cru_tmp", "chirps_v3_diag"):
            d.precip[name] = _align(t, z["precip"], time)
        for var, alias in extras.get(name, {}).items():
            if var in z:
                d.extra[alias] = _align(t, z[var], time)
        if "valid_frac" in z:
            d.extra[f"{name}_valid_frac"] = z["valid_frac"]
    # ASAP blend: CHIRPS v2 within +-50 deg, ERA5 beyond.
    if "chirps_v2" in d.precip and "era5" in d.precip:
        inside = np.abs(cells["lat"]) < 50
        d.precip["asap"] = np.where(inside[None, :], d.precip["chirps_v2"], d.precip["era5"])
    # CHIRPS cells with < half their 0.05 deg pixels valid (coasts) are dropped.
    for p in ("chirps_v2", "chirp_v2", "chirps_v3", "chirp_v3"):
        vf = d.extra.get(f"{p}_valid_frac")
        if p in d.precip and vf is not None:
            d.precip[p][:, vf < 0.5] = np.nan
    return d


def products_present(d: Data) -> list[str]:
    return [k for k, _, _ in PRODUCTS if k in d.precip]


# ------------------------------------------------------------------ regions

def koppen_major(cells: dict) -> np.ndarray:
    out = np.full(cells["koppen"].shape, "", dtype=object)
    for name, rng in KOPPEN_MAJOR.items():
        out[np.isin(cells["koppen"], list(rng))] = name
    return out


def gauge_class(numgauge: np.ndarray, sel_t: np.ndarray) -> np.ndarray:
    """Per-cell class from the median GPCC gauge count over the window."""
    med = np.nanmedian(numgauge[sel_t], 0)
    return np.where(np.isnan(med), "", np.where(med < 0.5, "0 gauges", np.where(med < 2.5, "1–2 gauges", "≥3 gauges")))


def countries(cells: dict, cache: Path) -> pd.DataFrame:
    """Majority Natural Earth admin-0 unit per 0.5 deg land cell."""
    import geopandas as gpd
    import rasterio.features as rf
    from rasterio.transform import from_origin

    from src.fetch import download

    out = cache / "cell_country.parquet"
    if out.exists():
        return pd.read_parquet(out)
    zp = download(
        "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries.zip",
        cache / "ne_10m_admin_0_countries.zip",
    )
    gdf = gpd.read_file(zp)[["ADM0_A3", "NAME", "CONTINENT", "SUBREGION", "geometry"]].reset_index(drop=True)
    res = 0.05
    idx = rf.rasterize(
        ((g, i + 1) for i, g in enumerate(gdf.geometry)),
        out_shape=(3600, 7200),
        transform=from_origin(-180, 90, res, res),
        dtype="int32",
    )
    idx = idx[::-1]  # south-up to match cell iy (ascending latitude)
    blocks = idx.reshape(360, 10, 720, 10).transpose(0, 2, 1, 3).reshape(360, 720, 100)
    b = blocks[cells["iy"], cells["ix"]]  # (N, 100)
    code = np.zeros(b.shape[0], dtype=np.int32)
    for k in range(b.shape[0]):
        row = b[k][b[k] > 0]
        if row.size:
            code[k] = np.bincount(row).argmax()
    df = pd.DataFrame({"code": code})
    df["iso3"] = np.where(code > 0, gdf.ADM0_A3.reindex(code - 1).values, "")
    df["country"] = np.where(code > 0, gdf.NAME.reindex(code - 1).values, "")
    df["continent"] = np.where(code > 0, gdf.CONTINENT.reindex(code - 1).values, "")
    df["subregion"] = np.where(code > 0, gdf.SUBREGION.reindex(code - 1).values, "")
    df = df.drop(columns="code").fillna("")
    df.to_parquet(out)
    return df
