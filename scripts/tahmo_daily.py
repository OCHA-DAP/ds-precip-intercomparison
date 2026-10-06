"""TAHMO curated station rainfall (5-15 min) -> daily totals cache.

    python scripts/tahmo_daily.py   # data/gauges/tahmo/*.zip -> cache/tahmo/daily.parquet + stations.parquet

A day's total is kept only if at least DAY_MIN_FRAC of its expected records (from
the station's logging step) have a value; otherwise the day is missing. Days are
UTC calendar days (TAHMO's convention; the files don't state a time zone, which
does not matter at monthly scale). Data stay local: TAHMO's research licence
forbids sharing the data or derived products (see data/gauges/tahmo/MANIFEST.md).
"""

import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "gauges" / "tahmo"
OUT = ROOT / "cache" / "tahmo"
DAY_MIN_FRAC = 0.95
STEPS = [5.0, 10.0, 15.0, 20.0, 30.0, 60.0]
SPARSE_DAYS = [0]  # days with no regular logging interval (too few records), across all stations


def station_daily(fh) -> tuple[pd.Series, int]:
    df = pd.read_csv(fh, usecols=["time", "precip_mm", "final_quality_flag"])
    if df.empty:
        return pd.Series(dtype=float), 0
    if not set(df.final_quality_flag.dropna().unique()) <= {1.0}:
        raise ValueError(f"unexpected quality flags {set(df.final_quality_flag.dropna().unique())}")
    t = pd.to_datetime(df.time)
    v = pd.to_numeric(df.precip_mm, errors="coerce")
    if (v < 0).any():
        raise ValueError("negative precipitation")
    # expected records per day from that day's own logging interval (stations can switch
    # interval): the most common gap between consecutive records among the regular ones
    # (STEPS). A day with no regular gap has too few records to be complete -> missing.
    dt = (t.shift(-1) - t).dt.total_seconds() / 60
    day = t.dt.floor("D")
    frame = pd.DataFrame({"day": day, "v": v, "dt": dt.where(dt.isin(STEPS))})
    day_step = frame.groupby("day").dt.agg(lambda x: x.mode().iloc[0] if x.notna().any() else np.nan)
    if day_step.isna().all():
        raise ValueError("no regular logging interval on any day")
    sparse = day_step.isna()
    g = frame.groupby("day").v
    tot, n = g.sum(min_count=1), g.count()
    step = int(day_step.mode().iloc[0])
    SPARSE_DAYS[0] += int(sparse.sum())
    return tot.where((n >= DAY_MIN_FRAC * 1440 / day_step) & ~sparse), step


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    zips = sorted(SRC.glob("TAHMO_Historical_precipitation_dataset_*.zip"))
    if len(zips) != 26:
        raise FileNotFoundError(f"expected 26 TAHMO zips in {SRC}, found {len(zips)}")
    metas, parts = [], []
    for zp in zips:
        with zipfile.ZipFile(zp) as z:
            meta = pd.read_csv(z.open("metadata.csv"))
            meta["source_zip"] = zp.name
            for sid in meta.id:
                with z.open(f"{sid}.csv") as fh:
                    s, step = station_daily(fh)
                meta.loc[meta.id == sid, "step_min"] = step
                meta.loc[meta.id == sid, "complete_days"] = int(s.notna().sum())
                parts.append(pd.DataFrame({"station": sid, "date": s.index, "precip": s.values}))
            metas.append(meta)
            print(f"[tahmo] {zp.name}: {len(meta)} stations", flush=True)
    st = pd.concat(metas, ignore_index=True)
    if st.id.duplicated().any():
        raise ValueError(f"duplicate station ids: {st.id[st.id.duplicated()].tolist()}")
    pd.concat(parts, ignore_index=True).to_parquet(OUT / "daily.parquet")
    st.to_parquet(OUT / "stations.parquet")
    print(f"[tahmo] {len(st)} stations -> {OUT}; {SPARSE_DAYS[0]} station-days with too few records to have a regular interval set missing")


if __name__ == "__main__":
    main()
