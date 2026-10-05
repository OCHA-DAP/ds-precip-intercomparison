"""Where, when and which products sit apart: PCA of the products, per-cell outliers,
disagreement over time.

    python scripts/multivariate_outliers.py   # -> site/multivariate/data/outliers.json

1. PCA of the product correlation matrix per region (pooled over every cell-month
   of standardised, rain-relevant anomalies). Loadings on PC1-3 place each product
   in a common space; products that move together sit together. This is the
   agreement picture (shared signal and shared mistakes alike), not the skill one.
2. Per cell: each product's distance from the consensus of the other eleven,
   1 - corr(product, mean of the others), and which product is farthest.
3. Per month: spread of the 12 products' standardised anomalies (standard
   deviation across products), averaged over a region -> when they disagree.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import multivariate as mv  # noqa: E402
from src import metrics as M  # noqa: E402
from src.adata import G_MOST, G_NEVER, G_SOME, countries, gauge_class, load  # noqa: E402

P = mv.P
K = len(P)
CONTINENTS = ["Africa", "Asia", "Europe", "North America", "South America", "Oceania"]


def pooled_corr(z: np.ndarray, sel: np.ndarray) -> np.ndarray:
    X = z[:, sel].reshape(-1, K)
    X = X[np.isfinite(X).all(1)].astype(np.float64)
    return np.corrcoef(X, rowvar=False)


def pca(R: np.ndarray) -> dict:
    w, V = np.linalg.eigh(R)
    o = np.argsort(w)[::-1]
    w, V = w[o], V[:, o]
    L = V[:, :3] * np.sqrt(np.clip(w[:3], 0, None))
    # orientation: PC1 positive for most products; PC2/PC3 positive for GPCC Full
    if L[:, 0].sum() < 0:
        L[:, 0] *= -1
    for j in (1, 2):
        if L[P.index("gpcc_full"), j] < 0:
            L[:, j] *= -1
    return {"loadings": {p: L[i].tolist() for i, p in enumerate(P)}, "var_frac": (w / w.sum())[:6].tolist()}


def main() -> None:
    d = load()
    cells = d.cells
    tb = d.window(*mv.WINDOW)
    mo = d.months[tb]
    times = [f"{t:%Y-%m}" for t in d.time[tb]]
    in_dom = np.abs(cells["lat"]) < 50
    cont = countries(cells, mv.CACHE).continent.values.astype(object)[in_dom]
    gclass = gauge_class(d.extra["gpcc_numgauge"], tb)[in_dom]
    base = {p: d.precip[p][tb][:, in_dom] for p in P}
    clim = M.monthly_clim(base[mv.REF_CLIM], mo)
    rel = clim[mo - 1] >= mv.RAIN_RELEVANT
    z = mv.stack(base, rel, mo)  # (T, N, K), jointly masked
    regions = {"All": np.ones(len(cont), bool)} | {c: cont == c for c in CONTINENTS} | \
        {g: gclass == g for g in (G_NEVER, G_SOME, G_MOST)}

    # 1. PCA per region
    pcs = {reg: pca(pooled_corr(z, sel)) | {"cells": int(sel.sum())} for reg, sel in regions.items()}
    print("[outliers] PCA done")

    # 2. per-cell distance of each product from the mean of the others
    ok = np.isfinite(z).all(-1)
    n = ok.sum(0)
    tot = np.nansum(z, -1, keepdims=True)
    others = (tot - z) / (K - 1)  # (T, N, K) leave-one-out consensus
    dist = np.full((z.shape[1], K), np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        for k in range(K):
            dist[:, k] = 1 - M.corr_cols(np.where(ok, z[..., k], np.nan), np.where(ok, others[..., k], np.nan), min_n=mv.MIN_MONTHS)
    good = np.isfinite(dist).all(1)
    far = np.where(good, np.nanargmax(np.nan_to_num(dist, nan=-1), 1), -1)
    srt = np.sort(np.nan_to_num(dist, nan=-1), 1)
    gap = np.where(good, srt[:, -1] - srt[:, -2], np.nan)  # how much farther than the runner-up
    share = {reg: {p: float(np.mean(far[sel & good] == i)) for i, p in enumerate(P)} | {"_cells": int((sel & good).sum())}
             for reg, sel in regions.items()}
    med_dist = {reg: {p: float(np.nanmedian(dist[sel & good, i])) for i, p in enumerate(P)} for reg, sel in regions.items()}
    print(f"[outliers] per-cell outliers: {good.sum()} cells")

    # 3. disagreement per month: SD across products, median over the region's cells
    with np.errstate(invalid="ignore"):
        spread = np.where(ok, np.nanstd(z, -1), np.nan)  # (T, N)
    ts = {reg: [float(v) if np.isfinite(v) else None for v in np.nanmedian(spread[:, sel], 1)] for reg, sel in regions.items()}
    disagree_map = np.nanmean(spread, 0)

    # per-cell fields aligned with fields.json (common cells of multivariate.py)
    fields = __import__("json").loads((mv.OUT / "fields.json").read_text())
    dom_idx = np.full(len(in_dom), -1)
    dom_idx[np.nonzero(in_dom)[0]] = np.arange(in_dom.sum())
    # fields.json lists cells by (iy, ix); map back to domain index
    key = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(cells["iy"], cells["ix"]))}
    idx = np.array([dom_idx[key[(a, b)]] for a, b in zip(fields["iy"], fields["ix"])])
    if (idx < 0).any():
        raise ValueError("fields.json cell outside the common domain — rerun multivariate.py")

    def u8(x, hi):
        q = np.round(np.clip(x, 0, hi) / hi * 254)
        return np.where(np.isfinite(x), q, 255).astype(np.uint8).tolist()

    mv.dump("outliers", {
        "products": P, "pca": pcs, "share_farthest": share, "median_distance": med_dist,
        "spread_ts": {"times": times, "series": ts},
        "cell": {"farthest": far[idx].astype(int).tolist(), "gap": u8(gap[idx], 0.5),
                 "far_dist": u8(np.nanmax(dist, 1)[idx], 1.0), "disagree": u8(disagree_map[idx], 1.0)},
    })


if __name__ == "__main__":
    main()
