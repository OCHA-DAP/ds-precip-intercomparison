"""Validate the products, and the multivariate model's skill estimates, against
independent dense rain-gauge networks (AMMA-CATCH, Niger / Benin / Mali).

    # data: AMMA-CATCH zips in data/gauges/amma-catch/ (git-ignored; see its MANIFEST.md
    # for licence CC-BY 4.0 and the mandatory citations)
    python scripts/validate_gauges.py --overlap   # step 0: are these gauges already in GPCC / CHIRPS?
    python scripts/validate_gauges.py             # -> site/multivariate/data/gauges.json

Step 0 decides whether the gauges are an independent witness: per AMMA cell, GPCC's
and CHIRPS v3's own gauge counts against the AMMA station count, and AMMA stations
matched (within MATCH_KM) to the CHIRPS v3 station lists. Matched stations are
dropped from the gauge areal mean.

Gauge "truth" per 0.5 deg cell-month: mean of the complete station-months (every
day present) of the unmatched stations, at least MIN_STATIONS of them. Daily values
are 24 h ending 06:00 UTC and are assigned to the previous day.
"""

import argparse
import io
import re
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import multivariate as mv  # noqa: E402
from src import metrics as M  # noqa: E402
from src.adata import load  # noqa: E402
from src.fetch import download  # noqa: E402

GAUGE_DIR = mv.ROOT / "data" / "gauges" / "amma-catch"
CHIRPS_LIST = "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/diagnostics/list_of_stations_used/monthly/global.stationsUsed.v3.{y}.{m:02d}.csv"
# daily datasets used; the irregular reading-day records (CE.RainD_Nc, CE.RainD_Nct) and
# the secondary flux/met stations (AE.*, AL.*, PA.*) are left out
DAILY = ["CL.Rain_O", "CL.Rain_Od", "CL.Rain_Cotonou", "CL.Rain_G", "CL.Rain_GT", "CL.RainD_Gh",
         "CL.Rain_N", "CE.Rain_Nc", "CE.Rain_Nc_9004", "CL.Rain_Nig"]
MATCH_KM = 2.0
MIN_STATIONS = 3
KNOWN_FLAGS = {"checked data", "missing value", "Cumulative rainfall estimated from control bucket"}


def read_station(text: str) -> tuple[dict, pd.Series]:
    meta, rows = {}, []
    for line in text.splitlines():
        if line.startswith("#"):
            m = re.match(r"# (station_name|station_latitude[^;]*|station_longitude[^;]*);([^;]*)", line)
            if m:
                meta[m.group(1).split(" ")[0]] = m.group(2).strip()
            continue
        parts = line.split(";")
        if len(parts) > 2 and parts[0][:2] in ("19", "20"):
            rows.append((parts[0], parts[1], parts[2].strip()))
    flags = {f for _, _, f in rows}
    unknown = flags - KNOWN_FLAGS
    if unknown:
        raise ValueError(f"{meta.get('station_name')}: unknown quality flags {unknown}")
    # only quality-checked daily values; "missing value" and the control-bucket gap totals
    # (negative, repeated over the gap's days) are treated as missing
    s = pd.Series([float(v) if f == "checked data" else np.nan for _, v, f in rows], index=pd.to_datetime([t for t, _, _ in rows]))
    if (s < 0).any():
        raise ValueError(f"{meta.get('station_name')}: negative value flagged as checked data")
    s.index = (s.index - pd.Timedelta(hours=30)).normalize()  # 24 h ending 06:00 UTC on D -> day D-1
    return meta, s


def load_stations() -> pd.DataFrame:
    """Daily series of every station in the daily datasets: columns = station key."""
    recs, series = [], {}
    for code in DAILY:
        zf = next(GAUGE_DIR.glob(f"AMMA-CATCH_{code}_*.zip"), None)
        if zf is None:
            raise FileNotFoundError(f"{code}: zip not found in {GAUGE_DIR}")
        with zipfile.ZipFile(zf) as z:
            for name in z.namelist():
                if not name.endswith(".csv"):
                    continue
                meta, s = read_station(io.TextIOWrapper(z.open(name), encoding="utf-8", errors="replace").read())
                key = f"{code}|{meta['station_name']}"
                series.setdefault(key, []).append(s)
                recs.append({"key": key, "dataset": code, "name": meta["station_name"],
                             "lat": float(meta["station_latitude"]), "lon": float(meta["station_longitude"])})
    st = pd.DataFrame(recs).drop_duplicates("key").set_index("key")
    daily = pd.DataFrame({k: pd.concat(v).groupby(level=0).first() for k, v in series.items()})
    return st, daily


def monthly(daily: pd.DataFrame) -> pd.DataFrame:
    """Station-month totals; NaN unless every day of the month is present."""
    g = daily.groupby(daily.index.to_period("M"))
    tot, n = g.sum(min_count=1), g.count()
    days = pd.Series(tot.index.days_in_month, index=tot.index)
    return tot.where(n.eq(days, axis=0))


def dedupe(st: pd.DataFrame, mon: pd.DataFrame) -> pd.DataFrame:
    """Same station in several datasets (within 500 m): keep the one with the most complete months."""
    st = st.assign(nmon=mon.notna().sum().reindex(st.index).values).sort_values("nmon", ascending=False)
    keep = []
    for k, r in st.iterrows():
        if all(np.hypot((r.lat - st.at[j, "lat"]) * 111, (r.lon - st.at[j, "lon"]) * 111 * np.cos(np.radians(r.lat))) > 0.5 for j in keep):
            keep.append(k)
    return st.loc[keep]


def chirps_matches(st: pd.DataFrame) -> pd.Series:
    """Distance (km) from each AMMA station to the nearest CHIRPS v3 station in any July list 2001-2020."""
    pts = []
    for y in range(2001, 2021):
        f = download(CHIRPS_LIST.format(y=y, m=7), mv.CACHE / "chirps3_stations" / f"{y}07.csv")
        # station names contain unquoted commas: read only the first two fields
        lines = Path(f).read_text(errors="replace").splitlines()
        if not lines[0].lower().startswith("latitude,longitude"):
            raise ValueError(f"CHIRPS station list {y}-07: unexpected header {lines[0]!r}")
        xy = np.array([[float(v) for v in ln.split(",", 2)[:2]] for ln in lines[1:] if ln.strip()])
        pts.append(xy[(xy[:, 0] > 5) & (xy[:, 0] < 20) & (xy[:, 1] > -6) & (xy[:, 1] < 6)])
    pts = np.unique(np.vstack(pts).round(4), axis=0)
    d = [np.min(np.hypot((r.lat - pts[:, 0]) * 111, (r.lon - pts[:, 1]) * 111 * np.cos(np.radians(r.lat)))) for _, r in st.iterrows()]
    return pd.Series(d, index=st.index)


def cell_of(lat, lon):
    return np.floor((np.asarray(lat) + 90) / 0.5).astype(int), np.floor((np.asarray(lon) + 180) / 0.5).astype(int)


def overlap_table(st, mon, d) -> pd.DataFrame:
    tb = d.window(*mv.WINDOW)
    rainy = tb & np.isin(d.months, [6, 7, 8, 9])
    key = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(d.cells["iy"], d.cells["ix"]))}
    st = st.assign(iy=cell_of(st.lat, st.lon)[0], ix=cell_of(st.lat, st.lon)[1])
    per = mon.loc["2001-01":"2020-12"]
    rows = []
    for (iy, ix), grp in st.groupby(["iy", "ix"]):
        c = key.get((iy, ix))
        n_amma = per[grp.index].loc[per.index.month.isin([6, 7, 8, 9])].notna().sum(1)
        n_free = per[grp.index[~grp.in_chirps]].loc[per.index.month.isin([6, 7, 8, 9])].notna().sum(1)
        rows.append({"cell": f"{-89.75 + 0.5 * iy:.2f}N {-179.75 + 0.5 * ix:.2f}E", "stations": len(grp),
                     "in_chirps_v3": int(grp.in_chirps.sum()),
                     "amma_jjas_median": float(n_amma.median()), "amma_not_chirps_jjas_median": float(n_free.median()),
                     "gpcc_full_jjas_median": float(np.nanmedian(d.extra["gpcc_numgauge"][rainy, c])) if c is not None else np.nan,
                     "gpcc_mon_jjas_median": float(np.nanmedian(d.extra["gpccmon_numgauge"][rainy, c])) if c is not None else np.nan,
                     "chirps3_stn_jjas_median": float(np.nanmedian(d.extra["chirps3_stn"][rainy, c])) if c is not None and "chirps3_stn" in d.extra else np.nan})
    return pd.DataFrame(rows).sort_values("stations", ascending=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--overlap", action="store_true", help="only print the step-0 overlap table")
    a = ap.parse_args()
    st, daily = load_stations()
    mon = monthly(daily)
    print(f"[gauges] {len(st)} station records, {daily.shape[0]} days, {daily.index.min():%Y}-{daily.index.max():%Y}")
    st = dedupe(st, mon)
    mon = mon[st.index]
    print(f"[gauges] {len(st)} stations after removing duplicates across datasets")
    st["chirps_km"] = chirps_matches(st)
    st["in_chirps"] = st.chirps_km <= MATCH_KM
    print(f"[gauges] matched to a CHIRPS v3 station within {MATCH_KM} km: {st.in_chirps.sum()} "
          f"({', '.join(st[st.in_chirps].name)})")
    d = load()
    tab = overlap_table(st, mon, d)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(tab.to_string(index=False))
    if a.overlap:
        return
    score(st, mon, d, tab)


# Core-product pairs that are independent under BOTH shared-error explanations
# (satellites share mistakes; GPCC and ERA5 share mistakes): used to estimate how
# reliable the gauge areal mean itself is.
CLEAN_PAIRS = [("gpcc_full", "chirp_v3"), ("gpcc_full", "imerg_late"), ("chirp_v3", "era5"), ("imerg_late", "era5")]
N_BOOT = 500


def pooled_stats(X: np.ndarray, years: np.ndarray) -> dict:
    out = {}
    for y in np.unique(years):
        x = X[years == y]
        out[int(y)] = (len(x), x.sum(0), x.T @ x)
    return out


def cov_from(stats: dict, yrs) -> np.ndarray:
    n = sum(stats[y][0] for y in yrs)
    s = sum(stats[y][1] for y in yrs)
    ss = sum(stats[y][2] for y in yrs)
    mu = s / n
    return (ss - n * np.outer(mu, mu)) / (n - 1)


def evaluate(C: np.ndarray) -> dict:
    """C over [gauge] + P. Gauge-based rho^2 for every product, the gauge's own
    reliability, and the model's estimates under both explanations (products only)."""
    ix = {p: i + 1 for i, p in enumerate(mv.P)}
    with np.errstate(invalid="ignore", divide="ignore"):
        rg = [C[0, ix[a]] * C[0, ix[b]] / (C[ix[a], ix[b]] * C[0, 0]) for a, b in CLEAN_PAIRS]
        rho2_g = float(np.median(rg))
        r2 = {p: C[0, ix[p]] ** 2 / (C[0, 0] * C[ix[p], ix[p]]) for p in mv.P}
        gauge_based = {p: r2[p] / rho2_g for p in mv.P}
    model = {k: {p: float(v[0]) for p, v in mv.skill(C[None, 1:, 1:], pair=sc["pair"])["rho2"].items()}
             for k, sc in mv.SCENARIOS.items()}
    return {"rho2_gauge": rho2_g, "rho2_gauge_by_pair": [float(x) for x in rg], "r2_vs_gauge": r2,
            "gauge_based": gauge_based, "model": model}


def score(st, mon, d, tab) -> None:
    tb = d.window(*mv.WINDOW)
    times = d.time[tb]
    mo, yr = d.months[tb], d.years[tb]
    key = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(d.cells["iy"], d.cells["ix"]))}
    free = st[~st.in_chirps].copy()
    free["iy"], free["ix"] = cell_of(free.lat, free.lon)
    per = mon.reindex(pd.PeriodIndex(times, freq="M"))
    rows, X, Y, CELLS = [], [], [], []
    for (iy, ix), grp in free.groupby(["iy", "ix"]):
        c = key.get((iy, ix))
        if c is None:
            continue
        g = per[grp.index]
        n = g.notna().sum(1).to_numpy()
        gauge = np.where(n >= MIN_STATIONS, g.mean(1).to_numpy(), np.nan)
        if np.isfinite(gauge).sum() < 30:
            continue
        prods = np.stack([d.precip[p][tb][:, c] for p in mv.P], 1)
        clim = M.monthly_clim(d.precip[mv.REF_CLIM][tb][:, [c]], mo)[:, 0]
        ok = np.isfinite(gauge) & np.isfinite(prods).all(1) & (clim[mo - 1] >= mv.RAIN_RELEVANT)
        if ok.sum() < 30:
            continue
        raw = np.column_stack([gauge, prods]).astype(np.float64)
        raw[~ok] = np.nan
        z = np.column_stack([M.std_anom(raw[:, [j]], mo)[:, 0] for j in range(raw.shape[1])])
        z = z[ok]
        X.append(z)
        Y.append(yr[ok])
        CELLS.append(np.full(ok.sum(), len(rows)))
        rows.append({"cell": f"{-89.75 + 0.5 * iy:.2f}N {-179.75 + 0.5 * ix:.2f}E", "lat": -89.75 + 0.5 * iy, "lon": -179.75 + 0.5 * ix,
                     "stations": int(len(grp)), "months": int(ok.sum()), "median_stations_per_month": float(np.median(n[ok])),
                     "site": "Benin (Ouémé)" if iy < 200 else ("Mali (Gourma)" if ix < 360 else "Niger (Niamey)")})
    if not rows:
        raise ValueError("no cell has enough independent gauge-months")
    X, Y, CELLS = np.vstack(X), np.concatenate(Y), np.concatenate(CELLS)
    print(f"[gauges] {len(rows)} cells, {len(X)} cell-months scored")
    stats = pooled_stats(X, Y)
    yrs = sorted(stats)
    est = evaluate(cov_from(stats, yrs))
    rng = np.random.default_rng(0)
    boots = [evaluate(cov_from(stats, list(rng.choice(yrs, len(yrs))))) for _ in range(N_BOOT)]

    def ci(path):
        v = np.array([path(b) for b in boots], dtype=float)
        return [float(np.nanpercentile(v, 5)), float(np.nanpercentile(v, 95))]

    out = {"cells": rows, "n_cell_months": int(len(X)), "years": [int(min(yrs)), int(max(yrs))],
           "rho2_gauge": est["rho2_gauge"], "rho2_gauge_ci": ci(lambda b: b["rho2_gauge"]),
           "rho2_gauge_by_pair": dict(zip([f"{a}+{b}" for a, b in CLEAN_PAIRS], est["rho2_gauge_by_pair"])),
           "products": {}}
    for p in mv.P:
        out["products"][p] = {
            "r2_vs_gauge": float(est["r2_vs_gauge"][p]), "r2_ci": ci(lambda b: b["r2_vs_gauge"][p]),
            "gauge_based": float(est["gauge_based"][p]), "gauge_based_ci": ci(lambda b: b["gauge_based"][p]),
            **{f"model_{k}": est["model"][k][p] for k in mv.SCENARIOS},
            **{f"model_{k}_ci": ci(lambda b, k=k: b["model"][k][p]) for k in mv.SCENARIOS}}
    # which explanation do the gauges support? mean absolute gap over products
    gap = {k: float(np.mean([abs(out["products"][p][f"model_{k}"] - out["products"][p]["gauge_based"]) for p in mv.P]))
           for k in mv.SCENARIOS}
    ratio = lambda e, k: e[k]["era5"] / e[k]["imerg_late"]  # noqa: E731
    out["mean_abs_gap"] = gap
    out["era5_over_late"] = {"gauge": est["r2_vs_gauge"]["era5"] / est["r2_vs_gauge"]["imerg_late"],
                             "gauge_ci": ci(lambda b: b["r2_vs_gauge"]["era5"] / b["r2_vs_gauge"]["imerg_late"]),
                             **{k: ratio(est["model"], k) for k in mv.SCENARIOS}}
    # per-cell spread (no bootstrap)
    per_cell = []
    for i, r in enumerate(rows):
        sel = CELLS == i
        e = evaluate(np.cov(X[sel], rowvar=False))
        per_cell.append({"cell": r["cell"], "gauge_based": e["gauge_based"], "rho2_gauge": e["rho2_gauge"]})
    out["per_cell"] = per_cell
    out["overlap"] = tab.to_dict(orient="records")
    out["stations_in_chirps"] = st[st.in_chirps].name.tolist()
    mv.dump("gauges", out)


if __name__ == "__main__":
    main()
