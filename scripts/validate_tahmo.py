"""Score the products against TAHMO station rainfall (2015-2025), and calibrate what a
single gauge per cell costs, using AMMA-CATCH's dense cells.

    python scripts/tahmo_daily.py      # once: zips -> cache/tahmo/daily.parquet
    python scripts/validate_tahmo.py   # -> cache/tahmo/scores.json (aggregate scores only)

Differences from validate_gauges.py (AMMA-CATCH), and why:
- Window 2015-2025: TAHMO stations start in 2015 or later. GPCC Full ends in 2020, so it
  is left out here and GPCC Monitoring takes its place as the gauge yardstick in the pairs
  used to estimate the gauges' own reliability. "Rain-season month" still uses the
  long-term 2001-2020 GPCC Full climatology.
- Cells with 1 or 2 stations are scored too, as their own tiers: TAHMO is mostly one
  station per 0.5 deg cell. A single gauge is a point, not the cell's average, so it
  lowers every product's r^2; the AMMA calibration below measures by how much.
- Independence: stations within MATCH_KM of a CHIRPS v3 station in any July list
  2015-2025 are removed; GPCC Monitoring's gauge count in each cell is reported.
- Output is aggregate only (by region and tier; station and cell counts, no coordinates or
  series): TAHMO's research licence forbids passing on the data or derived products.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import multivariate as mv  # noqa: E402
import validate_gauges as vg  # noqa: E402
from src import metrics as M  # noqa: E402
from src.adata import load  # noqa: E402
from src.fetch import download  # noqa: E402

CACHE = mv.ROOT / "cache" / "tahmo"
WINDOW = (2015, 2025)
PRODS = [p for p in mv.P if p != "gpcc_full"]
CLEAN_PAIRS = [("gpcc_monitoring", "chirp_v3"), ("gpcc_monitoring", "imerg_late"), ("chirp_v3", "era5"), ("imerg_late", "era5")]
REGION = {**{c: "East Africa" for c in ["KE", "UG", "TZ", "RW", "ET", "BI"]},
          **{c: "West Africa" for c in ["GH", "NG", "ML", "TG", "BJ", "BF", "NE", "CI", "SN"]},
          **{c: "Central Africa" for c in ["CM", "CD", "TD"]},
          **{c: "Southern Africa" for c in ["MW", "ZM", "MZ", "ZA", "MG", "ZW", "LS"]}}
MIN_MONTHS = 24  # per cell, scored months
N_BOOT = 300
DROPPED_NO_VARIATION = [0]


def chirps_points(years) -> np.ndarray:
    pts = []
    for y in years:
        f = download(vg.CHIRPS_LIST.format(y=y, m=7), mv.CACHE / "chirps3_stations" / f"{y}07.csv")
        lines = Path(f).read_text(errors="replace").splitlines()
        if not lines[0].lower().startswith("latitude,longitude"):
            raise ValueError(f"CHIRPS station list {y}-07: unexpected header {lines[0]!r}")
        xy = np.array([[float(v) for v in ln.split(",", 2)[:2]] for ln in lines[1:] if ln.strip()])
        pts.append(xy[(xy[:, 0] > -36) & (xy[:, 0] < 38) & (xy[:, 1] > -20) & (xy[:, 1] < 55)])
    return np.unique(np.vstack(pts).round(4), axis=0)


def nearest_km(lat, lon, pts) -> np.ndarray:
    return np.array([np.min(np.hypot((a - pts[:, 0]) * 111, (b - pts[:, 1]) * 111 * np.cos(np.radians(a)))) for a, b in zip(lat, lon)])


def three_month(x: np.ndarray) -> np.ndarray:
    """3-month totals ending in each month (all three months present), along axis 0."""
    s = pd.DataFrame(x).rolling(3, min_periods=3).sum().to_numpy()
    return s[:, 0] if x.ndim == 1 else s


def zrows(gauge: np.ndarray, prods: np.ndarray, mo: np.ndarray, rel: np.ndarray):
    """Standardised anomalies (per calendar month) of gauge + products on the months
    where all are present and the month is rain-relevant."""
    ok = np.isfinite(gauge) & np.isfinite(prods).all(1) & rel
    if ok.sum() < MIN_MONTHS:
        return None, ok
    raw = np.column_stack([gauge, prods]).astype(np.float64)
    raw[~ok] = np.nan
    z = np.column_stack([M.std_anom(raw[:, [j]], mo)[:, 0] for j in range(raw.shape[1])])
    # a calendar month where any source has no year-to-year variation (e.g. always 0 mm)
    # has no anomaly to compare: drop those months (counted in DROPPED_NO_VARIATION)
    ok2 = ok & np.isfinite(z).all(1)
    DROPPED_NO_VARIATION[0] += int((ok & ~ok2).sum())
    if ok2.sum() < MIN_MONTHS:
        return None, ok2
    return z[ok2], ok2


def scores(X: np.ndarray, Y: np.ndarray, MO: np.ndarray, prods: list, rng) -> dict:
    """r^2, Spearman and dry-month hit rate against the gauge (column 0), the gauge's own
    reliability from CLEAN_PAIRS, and year-block bootstrap ranges."""
    ix = {p: i + 1 for i, p in enumerate(prods)}
    pairs = [(a, b) for a, b in CLEAN_PAIRS if a in ix and b in ix]

    def one(idx):
        C = np.cov(X[idx], rowvar=False)
        r2 = {p: C[0, ix[p]] ** 2 / (C[0, 0] * C[ix[p], ix[p]]) for p in prods}
        rg = float(np.median([C[0, ix[a]] * C[0, ix[b]] / (C[ix[a], ix[b]] * C[0, 0]) for a, b in pairs])) if pairs else np.nan
        R = pd.DataFrame(X[idx]).rank().to_numpy()
        sp = {p: float(np.corrcoef(R[:, 0], R[:, ix[p]])[0, 1]) for p in prods}
        return r2, sp, rg

    allidx = np.arange(len(X))
    r2, sp, rg = one(allidx)
    yrs = np.unique(Y)
    boots = []
    for _ in range(N_BOOT):
        yy = rng.choice(yrs, len(yrs))
        boots.append(one(np.concatenate([np.nonzero(Y == y)[0] for y in yy])))
    # dry-month hit rate: lowest third per (row group of) calendar month, pooled
    dry = np.zeros_like(X, dtype=bool)
    for m in np.unique(MO):
        sel = MO == m
        if sel.sum() >= 6:
            dry[sel] = X[sel] <= np.nanpercentile(X[sel], 100 / 3, axis=0)
    gd = dry[:, 0]
    out = {"n": int(len(X)), "years": [int(yrs.min()), int(yrs.max())], "rho2_gauge": rg,
           "rho2_gauge_ci": [float(np.nanpercentile([b[2] for b in boots], q)) for q in (5, 95)], "products": {}}
    for p in prods:
        b = np.array([x[0][p] for x in boots])
        out["products"][p] = {"r2": float(r2[p]), "r2_ci": [float(np.percentile(b, 5)), float(np.percentile(b, 95))],
                              "spearman": sp[p], "dry_hit": float(dry[gd, ix[p]].mean()) if gd.sum() >= 10 else None,
                              "dry_false_alarm": float((~gd[dry[:, ix[p]]]).mean()) if dry[:, ix[p]].sum() >= 10 else None,
                              "gauge_based": float(r2[p] / rg) if rg and np.isfinite(rg) else None}
    return out


def amma_single_station_calibration(d, rng) -> dict:
    """In AMMA's dense cells: score each product against single stations vs the cell mean,
    over the same months (2001-2020)."""
    st, daily = vg.load_stations()
    mon = vg.monthly(daily)
    st, mon, _ = vg.dedupe(st, mon)
    tb = d.window(*mv.WINDOW)
    mo, yr = d.months[tb], d.years[tb]
    per = mon.reindex(pd.PeriodIndex(d.time[tb], freq="M"))
    key = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(d.cells["iy"], d.cells["ix"]))}
    st["iy"], st["ix"] = vg.cell_of(st.lat, st.lon)
    cell_rows, single_rows = [], []
    for (iy, ix), grp in st.groupby(["iy", "ix"]):
        c = key.get((iy, ix))
        if c is None or len(grp) < 3:
            continue
        g = per[grp.index]
        n = g.notna().sum(1).to_numpy()
        cell_mean = np.where(n >= vg.MIN_STATIONS, g.mean(1).to_numpy(), np.nan)
        prods = np.stack([d.precip[p][tb][:, c] for p in PRODS], 1)
        clim = M.monthly_clim(d.precip[mv.REF_CLIM][tb][:, [c]], mo)[:, 0]
        rel = clim[mo - 1] >= mv.RAIN_RELEVANT
        for k in grp.index:
            one = g[k].to_numpy()
            both = np.isfinite(one) & np.isfinite(cell_mean)  # same months for both
            zc, okc = zrows(np.where(both, cell_mean, np.nan), prods, mo, rel)
            zs, oks = zrows(np.where(both, one, np.nan), prods, mo, rel)
            if zc is None or zs is None:
                continue
            cell_rows.append((zc, yr[okc], mo[okc]))
            single_rows.append((zs, yr[oks], mo[oks]))
    if not cell_rows:
        raise ValueError("no AMMA station-cell pairs for the calibration")
    Xc, Yc, Mc = (np.concatenate(v) for v in zip(*cell_rows))
    Xs, Ys, Ms = (np.concatenate(v) for v in zip(*single_rows))
    sc_cell, sc_single = scores(Xc, Yc, Mc, PRODS, rng), scores(Xs, Ys, Ms, PRODS, rng)
    ratio = {p: sc_single["products"][p]["r2"] / sc_cell["products"][p]["r2"] for p in PRODS}
    rank_c = pd.Series({p: sc_cell["products"][p]["r2"] for p in PRODS}).rank()
    rank_s = pd.Series({p: sc_single["products"][p]["r2"] for p in PRODS}).rank()
    return {"pairs": len(cell_rows), "n": int(len(Xs)), "cell_mean": sc_cell, "single_station": sc_single,
            "r2_ratio_single_over_cell": ratio, "rank_spearman": float(rank_c.corr(rank_s, method="spearman"))}


def main() -> None:
    rng = np.random.default_rng(0)
    d = load()
    st = pd.read_parquet(CACHE / "stations.parquet")
    daily = pd.read_parquet(CACHE / "daily.parquet").pivot(index="date", columns="station", values="precip").sort_index()
    daily.index = pd.DatetimeIndex(daily.index)
    mon = vg.monthly(daily.asfreq("D"))
    st = st[st.country_code.isin(REGION)].copy()
    st["region"] = st.country_code.map(REGION)
    print(f"[tahmo] {len(st)} African stations in the mapped regions")
    pts = chirps_points(range(WINDOW[0], WINDOW[1] + 1))
    st["chirps_km"] = nearest_km(st.latitude.values, st.longitude.values, pts)
    st["in_chirps"] = st.chirps_km <= vg.MATCH_KM
    print(f"[tahmo] within {vg.MATCH_KM} km of a CHIRPS v3 station: {int(st.in_chirps.sum())} of {len(st)}")
    st["iy"], st["ix"] = vg.cell_of(st.latitude.values, st.longitude.values)

    tb = d.window(*WINDOW)
    mo, yr = d.months[tb], d.years[tb]
    tb_clim = d.window(*mv.WINDOW)
    per = mon.reindex(pd.PeriodIndex(d.time[tb], freq="M"))
    key = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(d.cells["iy"], d.cells["ix"]))}
    gpccmon = d.extra["gpccmon_numgauge"]
    free = st[~st.in_chirps]
    groups, groups3 = {}, {}
    cells_info = []
    for (iy, ix), grp in free.groupby(["iy", "ix"]):
        c = key.get((int(iy), int(ix)))
        if c is None:
            continue
        g = per[[s_ for s_ in grp.id if s_ in per.columns]]
        n = g.notna().sum(1).to_numpy()
        gauge = np.where(n >= 1, g.mean(1).to_numpy(), np.nan)
        prods = np.stack([d.precip[p][tb][:, c] for p in PRODS], 1)
        clim = M.monthly_clim(d.precip[mv.REF_CLIM][tb_clim][:, [c]], d.months[tb_clim])[:, 0]
        rel = clim[mo - 1] >= mv.RAIN_RELEVANT
        z, ok = zrows(gauge, prods, mo, rel)
        if z is None:
            continue
        # 3-month totals: windows ending in each month, compared with the same window in other
        # years; rain-relevant if the long-term 3-month total is at least 3x the monthly threshold
        clim3 = M.monthly_clim(three_month(d.precip[mv.REF_CLIM][tb_clim][:, [c]]), d.months[tb_clim])[:, 0]
        z3, ok3 = zrows(three_month(gauge), three_month(prods), mo, clim3[mo - 1] >= 3 * mv.RAIN_RELEVANT)
        tier = "1 station" if np.median(n[ok]) < 1.5 else ("2 stations" if np.median(n[ok]) < 2.5 else "3+ stations")
        region = grp.region.mode().iloc[0]
        cells_info.append({"region": region, "tier": tier, "months": int(ok.sum()), "stations": len(grp),
                           "gpccmon_median": float(np.nanmedian(gpccmon[tb][:, c])), "gpccmon_frac_months_gt0": float(np.nanmean(gpccmon[tb][:, c] > 0))})
        for k_ in (("All", "All"), (region, "All"), ("All", tier), (region, tier)):
            groups.setdefault(k_, []).append((z, yr[ok], mo[ok]))
            if z3 is not None:
                groups3.setdefault(k_, []).append((z3, yr[ok3], mo[ok3]))
    print(f"[tahmo] scored cells: {len(cells_info)}")
    res = {}
    for (reg, tier), parts in sorted(groups.items()):
        X, Y, MO = (np.concatenate(v) for v in zip(*parts))
        if len(X) < 60:
            continue
        res[f"{reg}|{tier}"] = {"region": reg, "tier": tier, "cells": len(parts), **scores(X, Y, MO, PRODS, rng)}
        print(f"[tahmo] {reg:16} {tier:12} cells {len(parts):4} months {len(X):6}")
    res3 = {}
    for (reg, tier), parts in sorted(groups3.items()):
        X, Y, MO = (np.concatenate(v) for v in zip(*parts))
        if len(X) < 60:
            continue
        res3[f"{reg}|{tier}"] = {"region": reg, "tier": tier, "cells": len(parts), **scores(X, Y, MO, PRODS, rng)}
    ci = pd.DataFrame(cells_info)
    print("[tahmo] calibrating the single-station tier on AMMA-CATCH ...")
    cal = amma_single_station_calibration(d, rng)
    out = {"window": WINDOW, "products": PRODS, "n_stations_african": int(len(st)), "n_in_chirps": int(st.in_chirps.sum()),
           "in_chirps_by_country": st.groupby("country_code").in_chirps.agg(["sum", "count"]).astype(int).to_dict(orient="index"),
           "cells_by_region_tier": ci.groupby(["region", "tier"]).size().rename("cells").reset_index().to_dict(orient="records"),
           "gpccmon_in_cells": {"median_of_medians": float(ci.gpccmon_median.median()), "share_cells_median_gt0": float((ci.gpccmon_median > 0).mean()),
                                "by_tier": ci.groupby("tier").gpccmon_median.agg(["median", "mean"]).to_dict(orient="index")},
           "groups": res, "groups_3month": res3, "amma_calibration": cal}
    out["months_dropped_no_variation"] = DROPPED_NO_VARIATION[0]
    print(f"[tahmo] months dropped because a source had no year-to-year variation that calendar month: {DROPPED_NO_VARIATION[0]}")
    CACHE.mkdir(parents=True, exist_ok=True)
    (CACHE / "scores.json").write_text(json.dumps(mv.clean(out), separators=(",", ":")))
    print(f"[tahmo] wrote {CACHE / 'scores.json'}")


if __name__ == "__main__":
    main()
