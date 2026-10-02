"""Compute every metric from the packed cubes; write site JSON + per-cell fields.

    python scripts/pack.py        # once: blob -> cache/packed
    python scripts/analyze.py     # -> site/report/data/*.json, cache/fields.npz

Windows: BIAS 2001-2020 (UDel to 2017), LONG trends 1983-2020, RECENT trends
2001-2025. Reference for "bias vs gauges" is GPCC Full Data v2022 (stratified by
its gauge count); every product is also compared with every other (pairwise),
so distance to the station-free axes (CHIRP, ERA5) is always visible too.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import metrics as M  # noqa: E402
from src.adata import (  # noqa: E402
    FAMILY,
    G_MOST,
    G_NEVER,
    GAUGE_CLASSES,
    gauged_fraction,
    GAUGE,
    LABEL,
    ROOT,
    countries,
    gauge_class,
    koppen_major,
    load,
    products_present,
)

OUT = ROOT / "site" / "report" / "data"
CACHE = ROOT / "cache"
BIAS = (2001, 2020)
LONG = (1983, 2020)
RECENT = (2001, 2025)
REF = "gpcc_full"
MIN_REF_ANNUAL = 50.0  # mm/yr: below this, relative bias is not reported
RAIN_RELEVANT = 10.0  # mm/month in the reference climatology

# Triple-collocation partners: two products whose errors are plausibly
# independent of the target's and of each other (see methodology).
TC_PARTNERS = {
    "gpcc_full": ("chirp_v3", "era5"),
    "gpcc_monitoring": ("chirp_v3", "era5"),
    "cru": ("chirp_v3", "era5"),
    "cpc": ("chirp_v3", "era5"),
    "precl": ("chirp_v3", "era5"),
    "udel": ("chirp_v3", "era5"),
    "chirp_v2": ("gpcc_full", "era5"),
    "chirp_v3": ("gpcc_full", "era5"),
    "era5": ("gpcc_full", "chirp_v3"),
    "imerg_late": ("gpcc_full", "era5"),
    "imerg_final": ("chirp_v3", "era5"),
    "chirps_v2": ("imerg_late", "era5"),
    "chirps_v3": ("imerg_late", "era5"),
    "asap": ("imerg_late", "era5"),  # = CHIRPS v2 in the common domain; GPCC shares its stations
}

DIFF_PAIRS_LONG = [
    ("chirps_v2", "chirp_v2"),
    ("chirps_v3", "chirp_v3"),
    ("chirps_v2", "gpcc_full"),
    ("chirps_v3", "gpcc_full"),
    ("chirp_v3", "gpcc_full"),
    ("era5", "gpcc_full"),
    ("era5", "chirp_v3"),
    ("cru", "gpcc_full"),
    ("cpc", "gpcc_full"),
    ("precl", "gpcc_full"),
]
DIFF_PAIRS_RECENT = [
    ("imerg_final", "gpcc_monitoring"),
    ("imerg_late", "gpcc_monitoring"),
    ("imerg_late", "imerg_final"),
    ("chirps_v2", "gpcc_monitoring"),
    ("chirps_v3", "gpcc_monitoring"),
    ("chirp_v3", "gpcc_monitoring"),
    ("era5", "gpcc_monitoring"),
    ("cru", "gpcc_monitoring"),
    ("cpc", "gpcc_monitoring"),
]


def clean(o):
    """JSON-safe: NaN/inf -> None, numpy -> python."""
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else round(float(o), 4)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def dump(name: str, obj) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.json").write_text(json.dumps(clean(obj), separators=(",", ":")))
    print(f"[analyze] wrote {name}.json")


def regional_series(x: np.ndarray, w: np.ndarray, sel: np.ndarray, min_cov: float = 0.95) -> np.ndarray:
    """Area-weighted mean over cells `sel` per row; NaN if coverage < min_cov
    of the region's weight (so a missing tile can't fake a trend)."""
    xs, ws = x[:, sel], w[sel]
    fin = np.isfinite(xs)
    cov = (fin * ws).sum(1) / ws.sum()
    with np.errstate(invalid="ignore"):
        m = np.nansum(np.where(fin, xs, 0) * ws, 1) / (fin * ws).sum(1)
    return np.where(cov >= min_cov, m, np.nan)


def main() -> None:
    d = load()
    P = products_present(d)
    cells = d.cells
    N = cells["lat"].size
    w = cells["area"].astype(np.float64) * cells["landfrac"]
    months, years = d.months, d.years
    print(f"[analyze] {N} cells, {len(d.time)} months, products: {P}")

    ctry = countries(cells, CACHE)
    kg = koppen_major(cells)
    tb = d.window(*BIAS)
    gclass = gauge_class(d.extra["gpcc_numgauge"], tb) if "gpcc_numgauge" in d.extra else np.full(N, "")
    gauged_frac = gauged_fraction(d.extra["gpcc_numgauge"], tb) if "gpcc_numgauge" in d.extra else np.full(N, np.nan)
    lat = cells["lat"]
    latband = np.select(
        [lat < -23.5, lat < 0, lat < 23.5, lat < 50], ["23.5–60°S", "0–23.5°S", "0–23.5°N", "23.5–50°N"], "50–90°N"
    )

    # ---------------------------------------------------------------- climatology
    clim = {}
    for p in P:
        y0, y1 = BIAS if p != "udel" else (2001, 2017)
        sel = d.window(y0, y1)
        clim[p] = M.annual_clim(d.precip[p][sel], months[sel])
    common = np.all([np.isfinite(clim[p]) for p in P], 0)
    print(f"[analyze] common domain: {common.sum()} cells")
    regions = {"All products' common domain": common}
    for c in ["Africa", "Asia", "Europe", "North America", "South America", "Oceania"]:
        regions[f"{c}"] = (ctry.continent.values == c) & common
    for k in ["A Tropical", "B Arid", "C Temperate", "D Cold"]:
        regions[f"Köppen {k}"] = (kg == k) & common
    for g in GAUGE_CLASSES:
        regions[f"GPCC {g}"] = (gclass == g) & common
    africa = ctry.continent.values == "Africa"
    regions["Sahel (Africa 10–18°N)"] = africa & (lat >= 10) & (lat < 18) & common
    regions["Sahara margin (Africa 18–22°N)"] = africa & (lat >= 18) & (lat < 22) & common
    # Regions on each product's own domain (for products/regions outside 50S-50N)
    own_regions = {"Global land (each product's own domain)": np.ones(N, bool)}
    for b in ["50–90°N", "23.5–60°S"]:
        own_regions[f"Lat {b}"] = latband == b

    region_meta = {r: {"cells": int(m.sum()), "area_frac": float(w[m].sum() / w.sum())}
                   for r, m in {**regions, **own_regions}.items()}

    # rain-dominated months: CRU temperature >= 2 C where known
    tmp = d.extra.get("cru_tmp")
    rain_ok = np.isfinite(tmp) & (tmp >= 2.0) if tmp is not None else np.ones((len(d.time), N), bool)

    ref = d.precip[REF]
    ref_clim_m = M.monthly_clim(ref[tb], months[tb])  # (12, N)
    ref_ok = clim[REF] >= MIN_REF_ANNUAL

    # ---------------------------------------------------------------- bias vs GPCC
    fields = {"gpcc_gauged_frac": gauged_frac, "common": common.astype(np.float32)}  # per-cell fields for maps
    bias_tab = {}
    for p in P:
        x = d.precip[p][tb]
        r = ref[tb]
        cell_ratio = M.ratio_of_sums(x, r)
        fields[f"bias_{p}"] = np.where(ref_ok, cell_ratio, np.nan)
        fields[f"clim_{p}"] = clim[p]
        row = {}
        for reg, sel in {**regions, **own_regions}.items():
            sel2 = sel & ref_ok
            if sel2.sum() < 20:
                continue
            xr_ = np.where(rain_ok[tb], x, np.nan)
            lw = cells["lw_factor"][months[tb] - 1]  # (T, N)
            gauge_like = p in GAUGE
            row[reg] = {
                "all": M.regional_ratio(x, r, w, sel2),
                "rain": M.regional_ratio(xr_, r, w, sel2),
                # sensitivity: gauge products (incl. the reference) lifted by
                # Legates-Willmott; satellite/reanalysis products unchanged
                "lw": M.regional_ratio(x * (lw if gauge_like else 1.0), r * lw, w, sel2),
            }
        bias_tab[p] = row

    # pairwise ratio matrix (rows = product, cols = reference product)
    pair_bias = {}
    for reg in ["All products' common domain", f"GPCC {G_NEVER}", f"GPCC {G_MOST}", "Africa"]:
        sel = regions[reg] & ref_ok
        pair_bias[reg] = {a: {b: M.regional_ratio(d.precip[a][tb], d.precip[b][tb], w, sel) for b in P} for a in P}

    # undercatch premise check: CHIRPS v3 / v2 ratio vs LW factor
    if {"chirps_v2", "chirps_v3"} <= set(P):
        v3v2 = M.ratio_of_sums(d.precip["chirps_v3"][tb], d.precip["chirps_v2"][tb])
        lw_ann = np.nansum(cells["lw_factor"] * ref_clim_m, 0) / np.nansum(ref_clim_m, 0)
        ok = np.isfinite(v3v2) & np.isfinite(lw_ann) & common
        fields["v3_over_v2"] = v3v2
        fields["lw_annual"] = lw_ann
        undercatch_check = {
            "corr": float(np.corrcoef(v3v2[ok], lw_ann[ok])[0, 1]),
            "median_v3_over_v2": float(np.median(v3v2[ok])),
            "median_lw": float(np.median(lw_ann[ok])),
            "binned": [
                {"lw_lo": lo, "lw_hi": hi, "median_v3_over_v2": float(np.median(v3v2[ok & (lw_ann >= lo) & (lw_ann < hi)]))
                 if (ok & (lw_ann >= lo) & (lw_ann < hi)).sum() > 20 else None, "n": int((ok & (lw_ann >= lo) & (lw_ann < hi)).sum())}
                for lo, hi in [(1.0, 1.05), (1.05, 1.1), (1.1, 1.15), (1.15, 1.25), (1.25, 1.5), (1.5, 3.0)]
            ],
        }
    else:
        undercatch_check = None

    dump("bias", {"ref": REF, "window": BIAS, "regions": region_meta, "table": bias_tab,
                  "pairwise": pair_bias, "undercatch_check": undercatch_check})

    # ---------------------------------------------------------------- temporal agreement
    relevant = ref_clim_m[months[tb] - 1] >= RAIN_RELEVANT  # (T, N)
    z = {p: np.where(relevant, M.std_anom(d.precip[p][tb], months[tb]), np.nan).astype(np.float32) for p in P}

    def anom_of(p):
        x = d.precip[p][tb]
        return np.where(relevant, x - M.monthly_clim(x, months[tb])[months[tb] - 1], np.nan)

    anom_ref = anom_of(REF)
    agree = {}
    sp_ref = {}
    for p in P:
        sp_ref[p] = M.spearman_cols(z[p], z[REF])
        fields[f"spearman_{p}"] = sp_ref[p]
        ap = anom_of(p)
        sd_p = np.nanstd(np.where(np.isfinite(anom_ref), ap, np.nan), 0)
        sd_r = np.nanstd(np.where(np.isfinite(ap), anom_ref, np.nan), 0)
        with np.errstate(invalid="ignore", divide="ignore"):
            gamma = sd_p / sd_r
        agree[p] = {
            reg: {
                "spearman_median": float(np.nanmedian(sp_ref[p][sel])),
                "spearman_p25": float(np.nanpercentile(sp_ref[p][sel], 25)),
                "variability_ratio_median": float(np.nanmedian(gamma[sel & ref_ok])),
            }
            for reg, sel in regions.items() if sel.sum() >= 20
        }
    pair_sp = {}
    for reg in ["All products' common domain", f"GPCC {G_NEVER}", f"GPCC {G_MOST}"]:
        sel = regions[reg]
        pair_sp[reg] = {}
        for a in P:
            pair_sp[reg][a] = {}
            for b in P:
                if a == b:
                    pair_sp[reg][a][b] = 1.0
                elif b in pair_sp[reg] and a in pair_sp[reg][b]:
                    pair_sp[reg][a][b] = pair_sp[reg][b][a]
                else:
                    pair_sp[reg][a][b] = float(np.nanmedian(M.spearman_cols(z[a][:, sel], z[b][:, sel])))
    dump("agreement", {"window": BIAS, "ref": REF, "table": agree, "pairwise_spearman": pair_sp})

    # ---------------------------------------------------------------- triple collocation
    tc = {}
    for p in P:
        a, b = TC_PARTNERS.get(p, (None, None))
        if a not in P or b not in P:
            continue
        raw, _, _ = M.etc_rho2(z[p], z[a], z[b], clip=False)
        r2 = np.where(np.isfinite(raw), np.clip(raw, 0.0, 1.0), np.nan)
        fields[f"tc_{p}"] = r2
        fin = np.isfinite(raw) & common
        tc[p] = {"partners": [a, b],
                 "clipped_high_frac": float((raw[fin] > 1).mean()), "clipped_low_frac": float((raw[fin] < 0).mean()),
                 "regions": {
            reg: {"rho2_median": float(np.nanmedian(r2[sel])), "n": int(np.isfinite(r2[sel]).sum())}
            for reg, sel in regions.items() if sel.sum() >= 20}}
    dump("tc", {"window": BIAS, "table": tc})

    # ---------------------------------------------------------------- terciles
    terc = {}
    pair_kappa = {}
    sums3 = {p: pd.DataFrame(d.precip[p]).rolling(3).sum().to_numpy()[tb].astype(np.float32) for p in P}
    relevant3 = pd.DataFrame(np.where(np.isfinite(ref), ref, np.nan)).rolling(3).sum().to_numpy()[tb]
    ref3_clim = M.monthly_clim(relevant3, months[tb])
    rel3 = ref3_clim[months[tb] - 1] >= 3 * RAIN_RELEVANT
    cats = {}
    crefs = {}
    for scale, src, rel in (("1-month", {p: d.precip[p][tb] for p in P}, relevant), ("3-month", sums3, rel3)):
        cref = np.where(rel, M.tercile_cat(src[REF], months[tb]), np.nan)
        crefs[scale] = cref
        cats[scale] = {p: np.where(rel, M.tercile_cat(src[p], months[tb]), np.nan) for p in P}
        common_clim = {p: np.where(rel, M.tercile_cat(src[p], months[tb], ref=src[REF]), np.nan) for p in P}
        for p in P:
            terc.setdefault(p, {})[scale] = {
                reg: {
                    "own": M.dry_skill(cats[scale][p][:, sel], cref[:, sel]),
                    "common": M.dry_skill(common_clim[p][:, sel], cref[:, sel]),
                }
                for reg, sel in regions.items() if sel.sum() >= 20
            }
            if scale == "3-month":
                fields[f"pss3_{p}"] = M.dry_skill_cells(cats[scale][p], cref)
    for reg in ["All products' common domain", f"GPCC {G_NEVER}", f"GPCC {G_MOST}"]:
        sel = regions[reg]
        pair_kappa[reg] = {a: {b: M.dry_skill(cats["3-month"][a][:, sel], cats["3-month"][b][:, sel])["pss"] for b in P} for a in P}
    dump("terciles", {"window": BIAS, "ref": REF, "table": terc, "pairwise_pss_3month": pair_kappa})

    # ---------------------------------------------------------------- trends
    trends = {}
    trend_regions = dict(regions)
    trend_regions.update(own_regions)
    fixed = None
    if "gpcc_numgauge" in d.extra:
        ng = d.extra["gpcc_numgauge"]
        yrs = np.arange(LONG[0], LONG[1] + 1)
        has = np.stack([(ng[years == y] >= 1).sum(0) >= 6 for y in yrs])
        fixed = has.all(0)
        never = np.stack([(ng[years == y] >= 1).sum(0) == 0 for y in yrs]).mean(0) >= 0.9
        trend_regions["Fixed-gauge cells (≥1 GPCC gauge every year)"] = fixed & common
        trend_regions["Gauge-free cells (no GPCC gauge in ≥90% of years)"] = never & common
    for wname, (y0, y1) in (("long", LONG), ("recent", RECENT)):
        yrs = np.arange(y0, y1 + 1)
        trends[wname] = {"years": [y0, y1], "products": {}}
        for p in P:
            ann = M.annual_totals(d.precip[p], years, yrs)
            if np.isfinite(ann).mean(0).max() < 0.9:
                continue
            cell = M.sen_mk(ann)
            sig = M.bh_fdr(cell["p"], q=0.1)
            fields[f"trend_{wname}_{p}"] = cell["pct_dec"]
            fields[f"trendsig_{wname}_{p}"] = sig.astype(np.float32)
            reg_out = {}
            for reg, sel in trend_regions.items():
                if sel.sum() < 20:
                    continue
                s = regional_series(ann, w, sel)
                if np.isfinite(s).mean() < 0.9:
                    continue
                r = M.sen_mk(s[:, None])
                reg_out[reg] = {"pct_dec": float(r["pct_dec"][0]), "p": float(r["p"][0]),
                                "series": s, "frac_cells_sig": float(sig[sel].mean()),
                                "frac_cells_pos": float((cell["pct_dec"][sel] > 0).mean())}
            trends[wname]["products"][p] = {"regions": reg_out, "frac_cells_sig": float(np.nanmean(sig[np.isfinite(cell["p"])]))}
        # sign agreement across products available for the window
        # ASAP is CHIRPS v2 inside +-50 deg: counting both would double-vote CHIRPS v2
        avail = [p for p in P if f"trend_{wname}_{p}" in fields and p != "asap"]
        stack = np.stack([fields[f"trend_{wname}_{p}"] for p in avail])
        fields[f"agree_pos_{wname}"] = (stack > 0).sum(0).astype(np.float32)
        fields[f"agree_n_{wname}"] = np.isfinite(stack).sum(0).astype(np.float32)
        trends[wname]["available"] = avail
    # difference-series trends (product A - product B), % of B's mean per decade
    diff = {}
    for wname, (y0, y1), pairs in (("long", LONG, DIFF_PAIRS_LONG), ("recent", RECENT, DIFF_PAIRS_RECENT)):
        yrs = np.arange(y0, y1 + 1)
        diff[wname] = {}
        for a, b in pairs:
            if a not in P or b not in P:
                continue
            A = M.annual_totals(d.precip[a], years, yrs)
            B = M.annual_totals(d.precip[b], years, yrs)
            out = {}
            for reg, sel in trend_regions.items():
                if sel.sum() < 20:
                    continue
                sa, sb = regional_series(A, w, sel), regional_series(B, w, sel)
                ok = np.isfinite(sa) & np.isfinite(sb)
                if ok.mean() < 0.9:
                    continue
                r = M.sen_mk((sa - sb)[:, None])
                out[reg] = {"pct_dec_of_ref": float(100 * 10 * r["slope"][0] / np.nanmean(sb)), "p": float(r["p"][0])}
            dcell = M.sen_mk(A - B)
            sig = M.bh_fdr(dcell["p"], q=0.1)
            out["_cells"] = {"frac_sig": float(sig[np.isfinite(dcell["p"])].mean()) if np.isfinite(dcell["p"]).any() else None}
            fields[f"difftrend_{wname}_{a}_minus_{b}"] = np.where(np.isfinite(clim.get(b, np.nan)), 100 * 10 * dcell["slope"] / clim[b], np.nan)
            diff[wname][f"{a}|{b}"] = out
    trends["difference"] = diff
    dump("trends", trends)

    # ---------------------------------------------------------------- station networks + drift
    net = {"years": list(range(1981, int(years.max()) + 1))}
    yrs_all = np.arange(1981, int(years.max()) + 1)

    def per_year(arr, fn):
        return [fn(arr[years == y]) if (years == y).any() else None for y in yrs_all]

    if "gpcc_numgauge" in d.extra:
        ng = d.extra["gpcc_numgauge"]
        net["gpcc_full_gauges"] = per_year(ng, lambda a: float(np.nanmean(np.nansum(a, 1))) if np.isfinite(a).any() else None)
        net["gpcc_full_frac_cells_gauged"] = per_year(ng, lambda a: float(np.nanmean((a >= 1).mean(1))) if np.isfinite(a).any() else None)
    if "gpccmon_numgauge" in d.extra:
        ng = d.extra["gpccmon_numgauge"]
        # replicated from 1 deg: each 1 deg count appears in up to 4 land cells
        net["gpcc_monitoring_frac_cells_gauged"] = per_year(ng, lambda a: float(np.nanmean((a >= 1).mean(1))) if np.isfinite(a).any() else None)
    if "cru_stn" in d.extra:
        net["cru_frac_cells_with_station"] = per_year(d.extra["cru_stn"], lambda a: float(np.nanmean((a >= 1).mean(1))) if np.isfinite(a).any() else None)
    if "chirps3_stn" in d.extra:
        net["chirps_v3_stations"] = per_year(d.extra["chirps3_stn"], lambda a: float(np.nanmean(np.nansum(a, 1))) if np.isfinite(a).any() else None)
        net["chirps_v3_fill_frac"] = per_year(d.extra["chirps3_fill"], lambda a: float(np.nanmean(a[:, common])) if np.isfinite(a).any() else None)
    if "imerg_gauge_weight" in d.extra:
        net["imerg_final_gauge_weight"] = per_year(d.extra["imerg_gauge_weight"], lambda a: float(np.nanmean(a[:, common])) if np.isfinite(a).any() else None)
    # IMERG Late / Final, monthly, common domain — from the UNMASKED Late
    # archive, so the Oct 2024 calibration incident stays visible
    if {"imerg_late", "imerg_final"} <= set(P):
        zl = np.load(ROOT / "cache" / "packed" / "imerg_late.npz")
        tl = [str(t) for t in zl["time"]]
        XL = zl["precip"]
        xf = d.precip["imerg_final"]
        ser_t, ser_r = [], []
        for k, t in enumerate(d.time):
            ym = f"{t:%Y-%m}"
            if ym not in tl:
                continue
            a, b = XL[tl.index(ym)], xf[k]
            ok = common & np.isfinite(a) & np.isfinite(b)
            if ok.sum() > 0.9 * common.sum():
                ser_t.append(ym)
                ser_r.append(float(np.sum(a[ok] * w[ok]) / np.sum(b[ok] * w[ok])))
        net["imerg_late_over_final_monthly"] = {"time": ser_t, "ratio": ser_r}
    # annual product/reference ratios on the common domain (drift view)
    drift = {}
    ann_cache = {}

    def ann(p):  # annual totals per product, computed once
        if p not in ann_cache:
            ann_cache[p] = M.annual_totals(d.precip[p], years, yrs_all)
        return ann_cache[p]

    for refname in ("gpcc_full", "gpcc_monitoring"):
        if refname not in P:
            continue
        R = ann(refname)
        drift[refname] = {}
        for reg in ["All products' common domain", "Africa", "South America", "Asia", f"GPCC {G_NEVER}", f"GPCC {G_MOST}"]:
            sel = regions[reg]
            sr = regional_series(R, w, sel)
            drift[refname][reg] = {}
            for p in P:
                sp = regional_series(ann(p), w, sel)
                with np.errstate(invalid="ignore", divide="ignore"):
                    drift[refname][reg][p] = sp / sr
    # step tests
    steps = {}

    def step(p, refn, before, after, reg="All products' common domain"):
        if p not in P or refn not in P:
            return None
        sel = regions[reg]
        out = []
        for y0, y1 in (before, after):
            tw = d.window(y0, y1)
            out.append(M.regional_ratio(d.precip[p][tw], d.precip[refn][tw], w, sel))
        return {"before": before, "after": after, "ratio_before": out[0], "ratio_after": out[1], "change_pct": 100 * (out[1] / out[0] - 1)}

    steps["CPC 2006 network change (vs GPCC Full)"] = step("cpc", "gpcc_full", (1996, 2005), (2007, 2016))
    steps["IMERG Final TRMM→GPM 2014 (vs GPCC Monitoring)"] = step("imerg_final", "gpcc_monitoring", (2004, 2013), (2015, 2024))
    steps["IMERG Late TRMM→GPM 2014 (vs GPCC Monitoring)"] = step("imerg_late", "gpcc_monitoring", (2004, 2013), (2015, 2024))
    steps["IMERG Late drift since 2022 (vs IMERG Final)"] = step("imerg_late", "imerg_final", (2001, 2021), (2022, 2025))
    steps["IMERG Final 2021 calibration switch (vs GPCC Monitoring)"] = step("imerg_final", "gpcc_monitoring", (2011, 2020), (2021, 2025))
    steps["ERA5 1983–2000 → 2001–2020 (vs GPCC Full)"] = step("era5", "gpcc_full", (1983, 2000), (2001, 2020))
    steps["CHIRPS v2 1983–2000 → 2001–2020 (vs GPCC Full)"] = step("chirps_v2", "gpcc_full", (1983, 2000), (2001, 2020))
    steps["CHIRPS v3 1983–2000 → 2001–2020 (vs GPCC Full)"] = step("chirps_v3", "gpcc_full", (1983, 2000), (2001, 2020))
    dump("networks", {"networks": net, "drift": drift, "steps": steps})

    # ---------------------------------------------------------------- countries
    cdf = pd.DataFrame({"iso3": ctry.iso3.values, "country": ctry.country.values, "continent": ctry.continent.values})
    ng_med = np.nanmean(d.extra["gpcc_numgauge"][tb], 0) if "gpcc_numgauge" in d.extra else np.full(N, np.nan)
    ctab = []
    (OUT / "country").mkdir(parents=True, exist_ok=True)
    ann_all = {p: ann(p) for p in P}
    cm_all = {p: M.monthly_clim(d.precip[p][tb], months[tb]) for p in P}
    yrs_all_l = [int(y) for y in yrs_all]
    for iso, g in cdf.groupby("iso3"):
        if not iso or len(g) < 3:
            continue
        sel = np.zeros(N, bool)
        sel[g.index.values] = True
        row = {"iso3": iso, "country": g.country.iloc[0], "continent": g.continent.iloc[0], "cells": int(sel.sum()),
               "gauges_per_cell": float(np.nanmean(ng_med[sel])), "frac_gauge_free": float(np.nanmean(gauged_frac[sel] == 0)),
               "products": {}}
        series = {"years": yrs_all_l, "annual": {}, "clim": {}}
        for p in P:
            ok = sel & np.isfinite(clim[p])
            if ok.sum() < max(3, 0.5 * sel.sum()):
                continue
            s = regional_series(ann_all[p], w, sel, min_cov=0.9)
            series["annual"][p] = s
            cm = cm_all[p]
            series["clim"][p] = [float(np.nansum(cm[m][ok] * w[ok]) / w[ok].sum()) for m in range(12)]
            prow = {
                "clim": float(np.nansum(clim[p][ok] * w[ok]) / w[ok].sum()),
                "bias": M.regional_ratio(d.precip[p][tb], ref[tb], w, sel & ref_ok) if (sel & ref_ok).sum() >= 3 else None,
                "spearman": float(np.nanmedian(sp_ref[p][sel])),
                "pss3": M.dry_skill(cats["3-month"][p][:, sel], crefs["3-month"][:, sel])["pss"],
            }
            for wname, (y0, y1) in (("long", LONG), ("recent", RECENT)):
                ii = [yrs_all_l.index(y) for y in range(y0, y1 + 1) if y in yrs_all_l]
                sw = s[ii]
                if np.isfinite(sw).mean() >= 0.9 and len(ii) == (y1 - y0 + 1):
                    r = M.sen_mk(sw[:, None])
                    prow[f"trend_{wname}"] = float(r["pct_dec"][0])
                    prow[f"trend_{wname}_p"] = float(r["p"][0])
            row["products"][p] = prow
        ctab.append(row)
        (OUT / "country" / f"{iso}.json").write_text(json.dumps(clean(series), separators=(",", ":")))
    dump("countries", ctab)

    meta = {
        "products": [{"key": p, "label": LABEL[p], "family": FAMILY[p], "span": d.span(p)} for p in P],
        "windows": {"bias": BIAS, "long": LONG, "recent": RECENT},
        "reference": REF,
        "cells": N,
        "common_cells": int(common.sum()),
    }
    dump("meta", meta)
    np.savez_compressed(CACHE / "fields.npz", iy=cells["iy"], ix=cells["ix"], **{k: np.asarray(v, dtype=np.float32) for k, v in fields.items()})
    print(f"[analyze] fields: {len(fields)}")


if __name__ == "__main__":
    main()
