"""Which product where: multivariate agreement, skill without a reference, regimes.

    python scripts/pack.py           # once: blob -> cache/packed
    python scripts/multivariate.py   # -> site/multivariate/data/*.json

Four layers, each answering a different question (2001-2020, common domain,
rain-relevant months, standardised anomalies; monthly and 3-month sums):

1. Dependency  - which products move together, and how much of that is shared
                 error rather than shared signal (residual correlation after
                 removing the latent signal estimated in layer 3).
2. Consensus   - per-cell PCA of the product correlation matrix: share of the
                 variance on the first component, for all products and for the
                 four input-independent "core" products only. Agreement, not skill.
3. Skill       - squared correlation of each product with the unknown true
                 signal (rho^2). A one-factor model on four core products with no
                 shared inputs (an overdetermined triple collocation), then every
                 other product scored against the core members it shares no input
                 with. Same assumptions as triple collocation; the residual
                 correlations of layer 1 test them.
4. Regimes     - k-means on the per-cell rho^2 profile -> regions where the
                 ranking of products is similar; best product per regime, overall
                 and among products available in near-real time.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.vq import kmeans2
from scipy.spatial.distance import cdist, squareform

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import metrics as M  # noqa: E402
from src.adata import (  # noqa: E402
    G_MOST,
    G_NEVER,
    G_SOME,
    LABEL,
    ROOT,
    countries,
    gauge_class,
    koppen_major,
    load,
)

OUT = ROOT / "site" / "multivariate" / "data"
CACHE = ROOT / "cache"
WINDOW = (2001, 2020)
RAIN_RELEVANT = 10.0  # mm/month in the GPCC climatology, as in analyze.py
REF_CLIM = "gpcc_full"  # only defines which months are rain-relevant
MIN_MONTHS = 60  # per cell, months where every product is present

# Inputs each product is built from. Two products may be used to score each
# other only if these sets are disjoint (independent errors). G = the gauge
# network behind GPCC/GHCN/CAMS (all gauge analyses share stations); C = the
# CHIRP satellite + CHPclim base (v2 and v3 share IR and climatology);
# I = the IMERG satellite constellation; E = ERA5.
INPUTS = {
    "gpcc_full": {"G"},
    "gpcc_monitoring": {"G"},
    "cru": {"G"},
    "cpc": {"G"},
    "precl": {"G"},
    "chirp_v2": {"C"},
    "chirp_v3": {"C"},
    "chirps_v2": {"C", "G"},
    "chirps_v3": {"C", "G"},  # ERA5 gap-fill is 0.25% of the 50S-50N domain; cells with >1% are masked (FILL_MAX)
    "imerg_late": {"I"},
    "imerg_final": {"I", "G"},  # calibrated to GPCC
    "era5": {"E"},
}
P = list(INPUTS)  # UDel (ends 2017) and ASAP (= CHIRPS v2 here) left out
CORE = ["gpcc_full", "chirp_v3", "imerg_late", "era5"]  # one per input family
FILL_MAX = 0.01  # CHIRPS v3 cells with more ERA5 gap-fill than this are not scored

# The four input families are not fully independent: the four core products'
# triple-collocation estimates disagree in a pattern that one shared-error pair
# explains, and two candidate pairs explain it equally well (they cannot be told
# apart from four products). Every result is computed under both; the first is
# primary because the inconsistency grows with gauge density, which a shared
# satellite error predicts and a shared gauge/reanalysis error does not (and it is
# what the report's TC_PARTNERS already assume).
SCENARIOS = {
    "sat": {"label": "CHIRP and IMERG share satellite (infrared) error", "pair": ("C", "I")},
    "smooth": {"label": "GPCC and ERA5 share smoothing error", "pair": ("G", "E")},
}
PRIMARY = "sat"


def dependent(a: str, b: str, pair: tuple | None) -> bool:
    """True if products a and b may have correlated errors: a shared input, or
    inputs that form the scenario's shared-error pair."""
    ia, ib = INPUTS[a], INPUTS[b]
    if ia & ib:
        return True
    return bool(pair) and any({x, y} == set(pair) for x in ia for y in ib)
# Approximate latency of the product's first release, for the near-real-time cut.
LATENCY = {
    "gpcc_full": "years (not operational)",
    "gpcc_monitoring": "~2 months",
    "cru": "~months, yearly release",
    "cpc": "~1 day",
    "precl": "~1 month",
    "chirp_v2": "~2 days",
    "chirp_v3": "~2 days",
    "chirps_v2": "~3 weeks (prelim ~2 days)",
    "chirps_v3": "~3 weeks (prelim ~2 days)",
    "imerg_late": "~14 hours",
    "imerg_final": "~3.5 months",
    "era5": "~5 days (ERA5T)",
}
# Sensitivity cores. "gauge swap" replaces GPCC with CRU (same family, sparser
# network use); "calibrated core" deliberately breaks the independence assumption
# (IMERG Final is calibrated to GPCC) to show how shared error inflates scores.
ALT_CORES = {
    "gauge swap (CRU for GPCC)": {"core": ["cru", "chirp_v3", "imerg_late", "era5"]},
    "naive: all triplets, IMERG Final in core": {"core": ["gpcc_full", "chirp_v3", "imerg_final", "era5"], "ignore": True},
}
NRT = ["chirps_v2", "chirps_v3", "imerg_late", "era5", "cpc", "chirp_v3"]


def clean(o):
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
    return o


def dump(name: str, obj) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.json").write_text(json.dumps(clean(obj), separators=(",", ":")))
    print(f"[mv] wrote {name}.json")


def cell_cov(X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """X (T, N, K) with a joint mask (a month is NaN for all K or none).
    Returns per-cell covariance (N, K, K) and month count (N,)."""
    ok = np.isfinite(X[..., 0])
    n = ok.sum(0)
    Xc = np.where(ok[..., None], X - np.nanmean(X, 0, keepdims=True), 0.0)
    C = np.einsum("tnk,tnl->nkl", Xc, Xc, optimize=True) / np.maximum(n - 1, 1)[:, None, None]
    C[n < MIN_MONTHS] = np.nan
    return C, n


def skill(C: np.ndarray, core: list[str] = CORE, pair: tuple | None = None,
          allow_unscored: bool = False, ignore_dependence: bool = False) -> dict:
    """rho^2 of every product with the latent signal, per cell (var T = 1).

    Core loadings: a_i^2 = C_ij C_ik / C_jk over the triplets (i, j, k) of core
    products that are pairwise independent under `pair` (each is one triple
    collocation estimate); the median when several are valid. Others:
    b_p = C_pk / a_k, averaged over the core products k independent of p;
    rho^2_p = b_p^2 / C_pp. A product with no valid triplet or partner raises,
    unless `allow_unscored` (sensitivity runs), where it is NaN and reported.
    `ignore_dependence` uses every triplet regardless (the deliberately broken
    demonstration)."""
    ix = {p: i for i, p in enumerate(P)}
    dep = (lambda a, b: False) if ignore_dependence else (lambda a, b: dependent(a, b, pair))
    a2, spread, used = {}, {}, {}
    nan = np.full(C.shape[0], np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        for i in core:
            others = [k for k in core if k != i]
            est, trip = [], []
            for j_ in range(len(others)):
                for k_ in range(j_ + 1, len(others)):
                    j, k = others[j_], others[k_]
                    if dep(i, j) or dep(i, k) or dep(j, k):
                        continue
                    est.append(C[:, ix[i], ix[j]] * C[:, ix[i], ix[k]] / C[:, ix[j], ix[k]])
                    trip.append([j, k])
            if not est:
                raise ValueError(f"core {i}: no independent triplet under pair {pair}")
            est = np.stack(est)
            a2[i] = np.median(est, 0)
            spread[i] = (est.max(0) - est.min(0)) / C[:, ix[i], ix[i]]
            used[i] = trip
        rho2_raw = {}
        for p in P:
            if p in core:
                rho2_raw[p] = a2[p] / C[:, ix[p], ix[p]]
                continue
            partners = [k for k in core if not dep(p, k)]
            used[p] = partners
            if not partners:
                if not allow_unscored:
                    raise ValueError(f"{p}: no independent core partner under pair {pair}")
                rho2_raw[p] = nan.copy()
                continue
            b = np.mean([C[:, ix[p], ix[k]] / np.sqrt(np.clip(a2[k], 1e-9, None)) for k in partners], 0)
            rho2_raw[p] = b**2 / C[:, ix[p], ix[p]]
    rho2 = {p: np.clip(v, 0.0, 1.0) for p, v in rho2_raw.items()}
    clipped = {p: float(np.nanmean(rho2_raw[p] > 1)) for p in P}
    return {"rho2": rho2, "spread": spread, "clipped_high": clipped, "used": used}


def corr_from_cov(C: np.ndarray) -> np.ndarray:
    sd = np.sqrt(np.einsum("nkk->nk", C))
    with np.errstate(invalid="ignore", divide="ignore"):
        return C / sd[:, :, None] / sd[:, None, :]


def stack(src: dict, rel: np.ndarray, months: np.ndarray) -> np.ndarray:
    """Standardised anomalies, jointly masked across products -> (T, N, K)."""
    z = np.stack([np.where(rel, M.std_anom(src[p], months), np.nan) for p in P], -1).astype(np.float32)
    z[~np.isfinite(z).all(-1)] = np.nan
    return z


def analyse(z: np.ndarray, label: str, fill_bad: np.ndarray) -> dict:
    C, n = cell_cov(z)
    R = corr_from_cov(C)
    scen = {}
    for key, sc in SCENARIOS.items():
        r = skill(C, pair=sc["pair"])
        r["rho2"]["chirps_v3"] = np.where(fill_bad, np.nan, r["rho2"]["chirps_v3"])
        scen[key] = r
    sk = scen[PRIMARY]
    ok = np.isfinite(R).all((1, 2))
    # consensus: PC1 share of the correlation matrix (all products; core only)
    ev = np.full(len(n), np.nan)
    ev_core = np.full(len(n), np.nan)
    ci = [P.index(k) for k in CORE]
    ev[ok] = np.linalg.eigvalsh(R[ok])[:, -1] / len(P)
    ev_core[ok] = np.linalg.eigvalsh(R[ok][:, ci][:, :, ci])[:, -1] / len(CORE)
    # residual correlation after the latent signal: r_ij - rho_i rho_j
    rho = np.sqrt(np.stack([np.nan_to_num(sk["rho2"][p], nan=np.nan) for p in P], -1))
    resid = R - rho[:, :, None] * rho[:, None, :]
    resid[:, np.arange(len(P)), np.arange(len(P))] = np.nan  # 1 - rho^2 on the diagonal: not a dependency
    print(f"[mv] {label}: {ok.sum()} cells with all {len(P)} products, median n = {np.median(n[ok]):.0f}")
    # sensitivity: the same model with a different core set
    alt = {name: skill(C, cfg["core"], pair=SCENARIOS[PRIMARY]["pair"], allow_unscored=True,
                       ignore_dependence=cfg.get("ignore", False))["rho2"] for name, cfg in ALT_CORES.items()}
    return {"R": R, "resid": resid, "alt": alt, "scen": scen, "ok": ok, "n": n, "pc1": ev, "pc1_core": ev_core, **sk}


def med_matrix(A: np.ndarray, sel: np.ndarray) -> list:
    return np.nanmedian(A[sel], 0)


def main() -> None:
    d = load()
    missing = [p for p in P if p not in d.precip]
    if missing:
        raise FileNotFoundError(f"packed products missing: {missing} — run scripts/pack.py")
    cells = d.cells
    N = cells["lat"].size
    tb = d.window(*WINDOW)
    mo = d.months[tb]
    ctry = countries(cells, CACHE)
    kg = koppen_major(cells)
    gclass = gauge_class(d.extra["gpcc_numgauge"], tb)

    ref_clim = M.monthly_clim(d.precip[REF_CLIM][tb], mo)
    rel1 = ref_clim[mo - 1] >= RAIN_RELEVANT
    s3 = {p: pd.DataFrame(d.precip[p]).rolling(3).sum().to_numpy()[tb].astype(np.float32) for p in P}
    ref3 = M.monthly_clim(s3[REF_CLIM], mo)
    rel3 = ref3[mo - 1] >= 3 * RAIN_RELEVANT

    if "chirps3_fill" not in d.extra:
        raise FileNotFoundError("chirps_v3_diag not packed — needed to mask CHIRPS v3 cells with ERA5 gap-fill")
    with np.errstate(all="ignore"):
        fill = np.nanmean(d.extra["chirps3_fill"][tb], 0)
    fill_bad = ~(fill <= FILL_MAX)  # NaN fill counts as not verifiable -> not scored
    res = {
        "monthly": analyse(stack({p: d.precip[p][tb] for p in P}, rel1, mo), "monthly", fill_bad),
        "3-month": analyse(stack(s3, rel3, mo), "3-month", fill_bad),
    }
    common = res["monthly"]["ok"] & res["3-month"]["ok"]
    print(f"[mv] analysed cells (both scales): {common.sum()}")

    # split-half stability of the monthly skill ranking
    halves = {}
    for h, (y0, y1) in {"2001-2010": (2001, 2010), "2011-2020": (2011, 2020)}.items():
        th = d.window(y0, y1)
        moh = d.months[th]
        relh = ref_clim[moh - 1] >= RAIN_RELEVANT
        zh = stack({p: d.precip[p][th] for p in P}, relh, moh)
        Ch, _ = cell_cov(zh)
        hr = skill(Ch, pair=SCENARIOS[PRIMARY]["pair"])["rho2"]
        hr["chirps_v3"] = np.where(fill_bad, np.nan, hr["chirps_v3"])
        halves[h] = hr

    # ------------------------------------------------------------ regions
    regions = {"All": common}
    for g in (G_NEVER, G_SOME, G_MOST):
        regions[g] = common & (gclass == g)
    for c in ["Africa", "Asia", "Europe", "North America", "South America", "Oceania"]:
        regions[c] = common & (ctry.continent.values == c)
    for k in ["A Tropical", "B Arid", "C Temperate", "D Cold"]:
        regions[f"Köppen {k}"] = common & (kg == k)

    # ------------------------------------------------------------ layer 1: dependency
    dep = {}
    for scale, r in res.items():
        dep[scale] = {}
        for reg in ["All", G_NEVER, G_SOME, G_MOST, "Africa"]:
            sel = regions[reg]
            Rm = med_matrix(r["R"], sel)
            dist = np.clip(1 - Rm, 0, None)
            np.fill_diagonal(dist, 0)
            Z = linkage(squareform((dist + dist.T) / 2, checks=False), method="average")
            dn = dendrogram(Z, no_plot=True, labels=P)
            dep[scale][reg] = {"corr": Rm, "resid": med_matrix(r["resid"], sel),
                               "linkage": Z, "order": dn["ivl"], "cells": int(sel.sum())}

    # ------------------------------------------------------------ layer 3: skill tables
    skill_tab = {
        sc: {scale: {
            reg: {p: {"median": float(np.nanmedian(r["scen"][sc]["rho2"][p][sel])),
                      "p25": float(np.nanpercentile(r["scen"][sc]["rho2"][p][sel], 25)),
                      "p75": float(np.nanpercentile(r["scen"][sc]["rho2"][p][sel], 75))} for p in P}
            for reg, sel in regions.items() if sel.sum() >= 50} for scale, r in res.items()}
        for sc in SCENARIOS}
    # tetrad check that picked the primary scenario: r_ge r_cl - r_gc r_el (shared
    # error in g-e or c-l) and its control r_gc r_el - r_gl r_ec (neither), by gauge class
    Rm_ = res["monthly"]["R"]
    rr = lambda a, b: Rm_[:, P.index(a), P.index(b)]  # noqa: E731
    g_, e_, c_, l_ = CORE[0], CORE[3], CORE[1], CORE[2]
    tet = rr(g_, e_) * rr(c_, l_) - rr(g_, c_) * rr(e_, l_)
    ctl = rr(g_, c_) * rr(e_, l_) - rr(g_, l_) * rr(e_, c_)
    tetrad = {reg: {"delta": float(np.nanmedian(tet[regions[reg]])), "control": float(np.nanmedian(ctl[regions[reg]])),
                    "cells": int(regions[reg].sum())}
              for reg in [G_NEVER, G_SOME, G_MOST]}
    tetrad_africa = {}
    for g in (G_NEVER, G_SOME, G_MOST):
        sel = regions[g] & regions["Africa"]
        tetrad_africa[g] = {"delta": float(np.nanmedian(tet[sel])), "cells": int(sel.sum())}
    sensitivity = {
        scale: {name: {reg: {p: float(np.nanmedian(alt[p][regions[reg]])) for p in P}
                       for reg in ["All", G_NEVER, G_MOST, "Africa"]}
                for name, alt in r["alt"].items()}
        for scale, r in res.items()}
    core_spread = {scale: {k: float(np.nanmedian(r["spread"][k][common])) for k in CORE} for scale, r in res.items()}
    clipped = {scale: r["clipped_high"] for scale, r in res.items()}

    # ------------------------------------------------------------ layer 4: regimes
    # cluster on each cell's relative profile (rho^2 minus the cell's mean over
    # products): which products do better or worse here, not how well measured the
    # cell is overall. The level is reported per regime alongside.
    S_abs = np.stack([res["monthly"]["rho2"][p] for p in P], -1)[common]  # (Nc, K)
    S = S_abs - np.nanmean(S_abs, 1, keepdims=True)
    S = np.nan_to_num(S, nan=0.0)  # unscored CHIRPS v3 cells (ERA5 fill): neutral, not "bad"
    rng = np.random.default_rng(0)
    best_k, scores = None, {}
    samp = rng.choice(len(S), min(3000, len(S)), replace=False)
    D = cdist(S[samp], S[samp])
    for k in range(3, 11):
        cent, lab = kmeans2(S, k, seed=0, minit="++", iter=50)
        # silhouette on a sample (scipy has none): mean (b - a) / max(a, b)
        ls = lab[samp]
        a = np.array([D[i, ls == ls[i]].sum() / max((ls == ls[i]).sum() - 1, 1) for i in range(len(samp))])
        b = np.array([min(D[i, ls == c].mean() for c in np.unique(ls) if c != ls[i]) for i in range(len(samp))])
        scores[k] = float(np.mean((b - a) / np.maximum(a, b)))
    # smallest k within 0.02 of the best silhouette: fewer regimes, easier to read
    top = max(scores.values())
    best_k = min(k for k, s in scores.items() if s >= top - 0.02)
    cent, lab = kmeans2(S, best_k, seed=0, minit="++", iter=100)
    # order regimes by mean skill, highest first
    order = np.argsort(-cent.mean(1))
    remap = np.empty_like(order)
    remap[order] = np.arange(best_k)
    lab = remap[lab]
    cent = cent[order]
    regime = np.full(N, -1, np.int16)
    regime[common] = lab

    w = (cells["area"] * cells["landfrac"]).astype(np.float64)
    reg_tab = []
    for g in range(best_k):
        sel = regime == g
        med = {p: float(np.nanmedian(res["monthly"]["rho2"][p][sel])) for p in P}
        med3 = {p: float(np.nanmedian(res["3-month"]["rho2"][p][sel])) for p in P}
        rank = sorted(P, key=lambda p: -med[p] if np.isfinite(med[p]) else 9)
        nrt_rank = [p for p in rank if p in NRT]
        # best NRT product under each shared-error scenario; robust if they agree
        by_scen = {}
        for sc in SCENARIOS:
            m = {q: float(np.nanmedian(res["monthly"]["scen"][sc]["rho2"][q][sel])) for q in NRT}
            o = sorted(NRT, key=lambda q: -m[q])
            by_scen[sc] = {"best": o[0], "margin": m[o[0]] - m[o[1]], "rho2": m}
        # split-half: does each half pick the same best NRT product?
        half_best = {h: max(NRT, key=lambda p: float(np.nanmedian(hv[p][sel]))) for h, hv in halves.items()}
        cont = pd.Series(ctry.continent.values[sel]).value_counts(normalize=True).head(4)
        reg_tab.append({
            "id": g, "cells": int(sel.sum()), "area_frac": float(w[sel].sum() / w[common].sum()),
            "rho2_monthly": med, "rho2_3month": med3, "rank": rank, "best_nrt": nrt_rank[0],
            "best_nrt_margin": med[nrt_rank[0]] - med[nrt_rank[1]], "half_best_nrt": half_best, "best_nrt_by_scenario": by_scen,
            "robust": len({v["best"] for v in by_scen.values()}) == 1,
            "gauge_mix": {c: float(np.mean(gclass[sel] == c)) for c in (G_NEVER, G_SOME, G_MOST)},
            "koppen_mix": {k: float(np.mean(kg[sel] == k)) for k in ["A Tropical", "B Arid", "C Temperate", "D Cold"]},
            "continents": cont.to_dict(),
            "mean_rho2": float(np.nanmean([med[p] for p in P])),
            "pc1_all": float(np.nanmedian(res["monthly"]["pc1"][sel])),
            "pc1_core": float(np.nanmedian(res["monthly"]["pc1_core"][sel])),
        })

    # per-pixel best NRT product and its stability across halves
    nrt_stack = np.stack([res["monthly"]["rho2"][p] for p in NRT], -1)
    best_px = np.full(N, -1, np.int8)
    best_px[common] = np.nanargmax(np.where(np.isfinite(nrt_stack[common]), nrt_stack[common], -1), 1)
    srt = np.sort(np.nan_to_num(nrt_stack, nan=-1), 1)
    margin = np.where(common, srt[:, -1] - srt[:, -2], np.nan)
    hb = [np.argmax(np.stack([np.nan_to_num(hv[p], nan=-1) for p in NRT], -1), 1) for hv in halves.values()]
    px_stable = float(np.mean(hb[0][common] == hb[1][common]))
    alt_sc = [k for k in SCENARIOS if k != PRIMARY][0]
    alt_stack = np.stack([np.nan_to_num(res["monthly"]["scen"][alt_sc]["rho2"][p], nan=-1) for p in NRT], -1)
    best_alt = np.argmax(alt_stack, 1)
    px_robust = np.where(common, best_alt == best_px, False)

    # ------------------------------------------------------------ per-cell fields
    def u8(x, lo=0.0, hi=1.0):
        q = np.round((np.clip(x, lo, hi) - lo) / (hi - lo) * 254)
        return np.where(np.isfinite(x) & common, q, 255).astype(np.uint8)

    fields = {
        "iy": cells["iy"][common].tolist(), "ix": cells["ix"][common].tolist(),
        "regime": regime[common].tolist(), "best_nrt": best_px[common].tolist(),
        "margin": u8(margin, 0, 0.5)[common].tolist(),
        "pc1": u8(res["monthly"]["pc1"])[common].tolist(), "pc1_core": u8(res["monthly"]["pc1_core"])[common].tolist(),
        "gauge": pd.Categorical(gclass[common], categories=[G_NEVER, G_SOME, G_MOST]).codes.tolist(),
        "best_alt": best_alt[common].astype(np.int8).tolist(), "robust": px_robust[common].astype(np.int8).tolist(),
        "rho2": {sc: {s: {p: u8(r["scen"][sc]["rho2"][p])[common].tolist() for p in P} for s, r in res.items()}
                 for sc in SCENARIOS},
    }
    dump("fields", fields)
    dump("summary", {
        "window": WINDOW, "products": P, "labels": {p: LABEL[p] for p in P}, "inputs": {p: sorted(v) for p, v in INPUTS.items()},
        "core": CORE, "nrt": NRT, "latency": LATENCY, "cells": int(common.sum()),
        "regions": {r: int(s.sum()) for r, s in regions.items()},
        "dependency": dep, "skill": skill_tab, "core_spread": core_spread, "sensitivity": sensitivity, "alt_cores": {k: v["core"] for k, v in ALT_CORES.items()}, "clipped_high": clipped,
        "scenarios": {k: {"label": v["label"], "pair": v["pair"],
                          "used": {p: u for p, u in res["monthly"]["scen"][k]["used"].items()}} for k, v in SCENARIOS.items()},
        "primary": PRIMARY, "tetrad": tetrad, "tetrad_africa": tetrad_africa,
        "pixel_best_robust_frac": float(px_robust[common].mean()), "fill_max": FILL_MAX,
        "chirps_v3_unscored_cells": int((fill_bad & common).sum()),
        "regimes": reg_tab, "k_silhouette": scores, "k": best_k, "pixel_best_nrt_stable_frac": px_stable,
        "consensus": {s: {reg: {"pc1_all": float(np.nanmedian(r["pc1"][sel])), "pc1_core": float(np.nanmedian(r["pc1_core"][sel]))}
                          for reg, sel in regions.items() if sel.sum() >= 50} for s, r in res.items()},
    })


if __name__ == "__main__":
    main()
