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
SITE_OF = {"CL.Rain_O": "Benin, upper Ouémé", "CL.Rain_Od": "Benin, Donga (in the Ouémé)", "CL.Rain_Cotonou": "Benin, Cotonou",
           "CL.Rain_G": "Mali, Gourma (Hombori)", "CL.Rain_GT": "Mali, Gourma transect", "CL.RainD_Gh": "Mali, Gourma (Hombori, reading gauges)",
           "CL.Rain_N": "Niger, Niamey (long-term)", "CE.Rain_Nc": "Niger, Niamey (2004 on)", "CE.Rain_Nc_9004": "Niger, Niamey (1990-2004)",
           "CL.Rain_Nig": "Niger, other towns"}
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
        # exact name: a glob "AMMA-CATCH_CE.Rain_Nc_*" would also match the CE.Rain_Nc_9004 zip
        pat = re.compile(rf"AMMA-CATCH_{re.escape(code)}_(daily|Precipitation_[A-Za-z_]+)\.zip")
        hits = [p for p in GAUGE_DIR.iterdir() if pat.fullmatch(p.name)]
        if len(hits) != 1:
            raise FileNotFoundError(f"{code}: expected one zip in {GAUGE_DIR}, found {[h.name for h in hits]}")
        zf = hits[0]
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


def dedupe(st: pd.DataFrame, mon: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Same station in several datasets (within 500 m), e.g. a series split between
    CE.Rain_Nc_9004 and CE.Rain_Nc at 2004: keep the record with the most complete
    months and fill its missing months from the others. Where both have a month they
    should agree; the kept record's value is used and disagreements over 1 mm are printed."""
    st = st.assign(nmon=mon.notna().sum().reindex(st.index).values).sort_values("nmon", ascending=False)
    mon = mon.copy()
    keep, filled, merged = [], 0, {}
    for k, r in st.iterrows():
        near = [j for j in keep
                if np.hypot((r.lat - st.at[j, "lat"]) * 111, (r.lon - st.at[j, "lon"]) * 111 * np.cos(np.radians(r.lat))) <= 0.5]
        if not near:
            keep.append(k)
            continue
        j = near[0]
        merged[k] = j
        both = mon[j].notna() & mon[k].notna()
        bad = both & ((mon[j] - mon[k]).abs() > 1.0)
        if bad.any():
            print(f"[gauges] {j} vs {k}: {int(bad.sum())} months differ by >1 mm; kept {j}'s values")
        gap = mon[j].isna() & mon[k].notna()
        mon.loc[gap, j] = mon.loc[gap, k]
        filled += int(gap.loc["2001-01":"2020-12"].sum())
    print(f"[gauges] duplicates merged: {len(st) - len(keep)} records; {filled} station-months (2001-2020) filled from them")
    return st.loc[keep], mon[keep], merged


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


def viz_export(st, st_all, merged, mon, rows, series, n_cm_ge3, X, Y, CELLS, MO, dry, valid, times, d, mon_all) -> None:
    """Everything the explainer section draws, straight from this run -> gauge_viz.json."""
    import geopandas as gpd
    from shapely.geometry import box

    used_cells = {(s_["iy"], s_["ix"]) for s_ in series}
    used_st = {k for s_ in series for k in s_["stations"]}
    per = mon.loc["2001-01":"2020-12"]
    complete = per.notna().sum()
    stations = []
    for k, r in st_all.iterrows():
        iy, ix = cell_of(r.lat, r.lon)
        if k in merged:
            fate, why = "merged", f"duplicate of {st_all.at[merged[k], 'name']} ({merged[k].split('|')[0]}), merged into it"
        elif bool(st.at[k, "in_chirps"]):
            fate, why = "chirps", f"used by CHIRPS v3 (station {st.at[k, 'chirps_km']:.1f} km away) - removed"
        elif complete.get(k, 0) == 0:
            fate, why = "nodata", "no complete month in 2001-2020"
        elif k in used_st:
            fate, why = "used", "used"
        else:
            fate, why = "fewcells", "its 0.5° cell never has 3 complete stations in enough months"
        stations.append({"key": k, "name": r["name"], "dataset": r["dataset"], "lat": round(float(r.lat), 4), "lon": round(float(r.lon), 4),
                         "fate": fate, "why": why, "complete_months": int(complete.get(k, 0)) if k in complete.index else 0,
                         "cell": [int(iy), int(ix)]})
    n_any = int(sum(1 for x in stations if x["fate"] in ("used", "fewcells")))
    funnel = [
        {"step": "Gauge records in the 10 daily datasets", "n": len(st_all), "unit": "records"},
        {"step": "After merging duplicates (same spot, within 500 m)", "n": len(st), "unit": "stations"},
        {"step": "After removing gauges CHIRPS v3 uses", "n": int((~st.in_chirps).sum()), "unit": "stations"},
        {"step": "With at least one complete month in 2001-2020", "n": n_any, "unit": "stations"},
        {"step": "...and in a 0.5° cell that reaches 3 complete stations: used", "n": int(sum(1 for x in stations if x["fate"] == "used")), "unit": "stations"},
        {"step": "Cell-months with 3 or more complete stations", "n": n_cm_ge3, "unit": "cell-months"},
        {"step": "...that are rain-season months with every product present", "n": int(len(X)), "unit": "cell-months"},
    ]
    # availability: complete (non-CHIRPS) station-months per scored cell per year, 1990-2025
    years = list(range(1990, 2026))
    avail = []
    for s_ in series:
        sub = mon[[k for k in s_["stations"]]]
        cnt = sub.notna().sum(1)
        avail.append([int(cnt[cnt.index.year == y].sum()) for y in years])
    # coverage: complete station-months per gauge dataset and year (before any selection),
    # and each product's span, for the timeline
    yrs_all = list(range(1989, 2026))
    datasets = []
    for code in DAILY:
        cols = [k for k in st_all.index if st_all.at[k, "dataset"] == code]
        sub = mon_all[cols]
        per_year = [int(sub[sub.index.year == y].notna().sum().sum()) for y in yrs_all]
        have = [y for y, v in zip(yrs_all, per_year) if v]
        n_used = sum(1 for x in stations if x["dataset"] == code and x["fate"] == "used")
        datasets.append({"code": code, "site": SITE_OF[code], "stations": len(cols), "used": n_used,
                         "first": have[0] if have else None, "last": have[-1] if have else None, "per_year": per_year})
    prod_span = {p: list(d.span(p)) for p in mv.P}
    # basemap: country outlines in the gauge region
    bbox = (-6.5, 4.0, 10.5, 18.5)
    zp = mv.CACHE / "ne_10m_admin_0_countries.zip"
    if not zp.exists():
        raise FileNotFoundError(f"{zp} missing — run analyze.py once (src.adata.countries downloads it)")
    gdf = gpd.read_file(zp)[["NAME", "geometry"]]
    gdf = gdf[gdf.intersects(box(*bbox))].copy()
    gdf["geometry"] = gdf.geometry.intersection(box(*bbox)).simplify(0.02)
    rings = []
    for name, geom in zip(gdf.NAME, gdf.geometry):
        polys = [geom] if geom.geom_type == "Polygon" else list(getattr(geom, "geoms", []))
        for pg in polys:
            if pg.geom_type == "Polygon":
                rings.append({"name": name, "xy": np.round(np.asarray(pg.exterior.coords), 3).tolist()})
    labels = {n: [round(c.x, 2), round(c.y, 2)] for n, c in zip(gdf.NAME, gdf.geometry.representative_point())}
    mv.dump("gauge_viz", {
        "stations": stations, "funnel": funnel, "coverage": {"years": yrs_all, "datasets": datasets, "products": prod_span,
                                                          "window": list(mv.WINDOW)}, "bbox": bbox, "countries": rings, "country_labels": labels,
        "cells": [{**r, "iy": s_["iy"], "ix": s_["ix"]} for r, s_ in zip(rows, series)],
        "availability": {"years": years, "counts": avail},
        "times": [f"{t:%Y-%m}" for t in times],
        "series": [{k_: v for k_, v in s_.items() if k_ != "stations"} for s_ in series],
        "pooled": {"z": np.round(X, 2).tolist(), "cell": CELLS.astype(int).tolist(), "month": MO.astype(int).tolist(),
                   "year": Y.astype(int).tolist(), "dry": dry.astype(int).tolist(), "valid": valid.astype(int).tolist()},
        "products": mv.P, "match_km": MATCH_KM, "min_stations": MIN_STATIONS,
    })


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--overlap", action="store_true", help="only print the step-0 overlap table")
    a = ap.parse_args()
    st, daily = load_stations()
    mon = monthly(daily)
    print(f"[gauges] {len(st)} station records, {daily.shape[0]} days, {daily.index.min():%Y}-{daily.index.max():%Y}")
    st_all, mon_all = st.copy(), mon.copy()
    st, mon, merged = dedupe(st, mon)
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
    score(st, mon, d, tab, st_all, merged, mon_all)


# Core-product pairs that are independent under BOTH shared-error explanations
# (satellites share mistakes; GPCC and ERA5 share mistakes): used to estimate how
# reliable the gauge areal mean itself is.
CLEAN_PAIRS = [("gpcc_full", "chirp_v3"), ("gpcc_full", "imerg_late"), ("chirp_v3", "era5"), ("imerg_late", "era5")]
N_BOOT = 500
MIN_GROUP = 40  # minimum cell-months for a disaggregated score


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


def score(st, mon, d, tab, st_all, merged, mon_all) -> None:
    tb = d.window(*mv.WINDOW)
    times = d.time[tb]
    mo, yr = d.months[tb], d.years[tb]
    key = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(d.cells["iy"], d.cells["ix"]))}
    free = st[~st.in_chirps].copy()
    free["iy"], free["ix"] = cell_of(free.lat, free.lon)
    per = mon.reindex(pd.PeriodIndex(times, freq="M"))
    rows, X, Y, CELLS, MO = [], [], [], [], []
    series, n_cm_ge3 = [], 0
    for (iy, ix), grp in free.groupby(["iy", "ix"]):
        c = key.get((iy, ix))
        if c is None:
            continue
        g = per[grp.index]
        n = g.notna().sum(1).to_numpy()
        gauge = np.where(n >= MIN_STATIONS, g.mean(1).to_numpy(), np.nan)
        if np.isfinite(gauge).sum() < 30:
            continue
        n_ge3 = int(np.isfinite(gauge).sum())
        prods = np.stack([d.precip[p][tb][:, c] for p in mv.P], 1)
        clim = M.monthly_clim(d.precip[mv.REF_CLIM][tb][:, [c]], mo)[:, 0]
        ok = np.isfinite(gauge) & np.isfinite(prods).all(1) & (clim[mo - 1] >= mv.RAIN_RELEVANT)
        if ok.sum() < 30:
            continue
        n_cm_ge3 += n_ge3
        series.append({"gauge": np.round(np.where(np.isfinite(gauge), gauge, -1), 1).tolist(), "n": n.astype(int).tolist(),
                       "ok": ok.astype(int).tolist(), "relevant": (clim[mo - 1] >= mv.RAIN_RELEVANT).astype(int).tolist(),
                       "products": {p: np.round(np.nan_to_num(prods[:, i], nan=-1), 1).tolist() for i, p in enumerate(mv.P)},
                       "stations": grp.index.tolist(), "iy": int(iy), "ix": int(ix)})
        raw = np.column_stack([gauge, prods]).astype(np.float64)
        raw[~ok] = np.nan
        z = np.column_stack([M.std_anom(raw[:, [j]], mo)[:, 0] for j in range(raw.shape[1])])
        z = z[ok]
        X.append(z)
        Y.append(yr[ok])
        MO.append(mo[ok])
        CELLS.append(np.full(ok.sum(), len(rows)))
        rows.append({"cell": f"{-89.75 + 0.5 * iy:.2f}N {-179.75 + 0.5 * ix:.2f}E", "lat": -89.75 + 0.5 * iy, "lon": -179.75 + 0.5 * ix,
                     "stations": int(len(grp)), "months": int(ok.sum()), "median_stations_per_month": float(np.median(n[ok])),
                     "site": "Benin (Ouémé)" if iy < 200 else ("Mali (Gourma)" if ix < 360 else "Niger (Niamey)")})
    if not rows:
        raise ValueError("no cell has enough independent gauge-months")
    X, Y, CELLS, MO = np.vstack(X), np.concatenate(Y), np.concatenate(CELLS), np.concatenate(MO)
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
    # rank-based scores against the gauges: Spearman correlation of the anomalies, and the
    # dry-month hit rate (gauge month in its driest third for that cell and calendar month ->
    # how often the product also puts it in its own driest third) with the false-alarm share
    dry = np.zeros_like(X, dtype=bool)
    for c in np.unique(CELLS):
        for m in np.unique(MO[CELLS == c]):
            sel = (CELLS == c) & (MO == m)
            if sel.sum() < 6:
                continue
            q = np.nanpercentile(X[sel], 100 / 3, axis=0)
            dry[sel] = X[sel] <= q
    valid = np.zeros(len(X), bool)
    for c in np.unique(CELLS):
        for m in np.unique(MO[CELLS == c]):
            sel = (CELLS == c) & (MO == m)
            valid[sel] = sel.sum() >= 6

    def rank_scores(idx):
        Xi, Di = X[idx], dry[idx]
        R = pd.DataFrame(Xi).rank().to_numpy()
        sp = {p: float(np.corrcoef(R[:, 0], R[:, i + 1])[0, 1]) for i, p in enumerate(mv.P)}
        gd = Di[:, 0]
        hit = {p: float(Di[gd, i + 1].mean()) for i, p in enumerate(mv.P)}
        fa = {p: float((~gd[Di[:, i + 1]]).mean()) for i, p in enumerate(mv.P)}
        return sp, hit, fa

    vidx = np.nonzero(valid)[0]
    sp, hit, fa = rank_scores(vidx)
    bs = []
    for _ in range(N_BOOT):
        yy = rng.choice(yrs, len(yrs))
        idx = np.concatenate([vidx[Y[vidx] == y] for y in yy])
        bs.append(rank_scores(idx))
    for i, (name, val) in enumerate((("spearman", sp), ("dry_hit", hit), ("dry_false_alarm", fa))):
        for p in mv.P:
            v = np.array([b[i][p] for b in bs])
            out["products"][p][name] = val[p]
            out["products"][p][f"{name}_ci"] = [float(np.percentile(v, 5)), float(np.percentile(v, 95))]
    out["n_dry_months"] = int(dry[valid, 0].sum())

    # disaggregated: r^2 and dry-month hit rate against the gauges by site x calendar month and
    # site x 5-year period, each with a year-block bootstrap. Groups under MIN_GROUP values are
    # reported with n only (no score).
    site_of = np.array([rows[c]["site"] for c in CELLS])
    periods = [(2001, 2005), (2006, 2010), (2011, 2015), (2016, 2020)]

    def group_scores(idx):
        if len(idx) < MIN_GROUP:
            return None
        Xi, Di, Yi = X[idx], dry[idx], Y[idx]
        C = np.corrcoef(Xi, rowvar=False)
        r2 = {p: float(C[0, i + 1] ** 2) for i, p in enumerate(mv.P)}
        gd = Di[:, 0] & valid[idx]
        hit = {p: (float(Di[gd, i + 1].mean()) if gd.sum() >= 10 else None) for i, p in enumerate(mv.P)}
        bs = []
        for _ in range(200):
            yy = rng.choice(np.unique(Yi), len(np.unique(Yi)))
            j = np.concatenate([np.nonzero(Yi == y)[0] for y in yy])
            Cb = np.corrcoef(Xi[j], rowvar=False)
            bs.append([Cb[0, i + 1] ** 2 for i in range(len(mv.P))])
        lo, hi = np.nanpercentile(np.array(bs), [5, 95], axis=0)
        return {"n": int(len(idx)), "n_dry": int(gd.sum()),
                "products": {p: {"r2": r2[p], "lo": float(lo[i]), "hi": float(hi[i]), "dry_hit": hit[p]} for i, p in enumerate(mv.P)}}

    by_month, by_period = {}, {}
    for site in ["All"] + sorted(set(site_of)):
        sel_s = np.ones(len(X), bool) if site == "All" else (site_of == site)
        by_month[site] = {}
        for m in range(1, 13):
            idx = np.nonzero(sel_s & (MO == m))[0]
            by_month[site][m] = group_scores(idx) if len(idx) else None
        by_period[site] = {f"{a}–{b}": group_scores(np.nonzero(sel_s & (Y >= a) & (Y <= b))[0]) for a, b in periods}
    out["by_month"], out["by_period"] = by_month, by_period
    out["min_group"] = MIN_GROUP

    # per-cell spread (no bootstrap)
    per_cell = []
    for i, r in enumerate(rows):
        sel = CELLS == i
        e = evaluate(np.cov(X[sel], rowvar=False))
        per_cell.append({"cell": r["cell"], "gauge_based": e["gauge_based"], "rho2_gauge": e["rho2_gauge"]})
    out["per_cell"] = per_cell
    viz_export(st, st_all, merged, mon, rows, series, n_cm_ge3, X, Y, CELLS, MO, dry, valid, times, d, mon_all)
    out["overlap"] = tab.to_dict(orient="records")
    out["stations_in_chirps"] = st[st.in_chirps].name.tolist()
    mv.dump("gauges", out)


if __name__ == "__main__":
    main()
