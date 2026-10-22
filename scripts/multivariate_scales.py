"""Does the multivariate picture change with scale, time, season and region?

    python scripts/multivariate.py          # first (writes summary.json used here for context)
    python scripts/multivariate_scales.py   # -> site/multivariate/data/scales.json

Same model as multivariate.py (one-factor extended triple collocation, CHIRP v3 and
IMERG Late treated as sharing mistakes), applied four more ways:

1. Aggregation scale - space: 0.5 deg cells, 1 deg, 2 deg, Natural Earth admin-1;
   time: 1-, 3- and 6-month sums. Skill per unit, median by continent.
2. Over time - rolling 5-year windows, 2001-2005 .. 2016-2020.
3. Season - month relative to each cell's climatological wettest month.
4. Grouping by region - the product clustering tree per region, compared with the
   global tree (cophenetic correlation) and each product's closest partner.

Checks 2 and 3 have too few months per cell for per-cell estimates (about 40 and
20), so they pool all cells of a region: covariances of the standardised
anomalies over every (cell, month) pair, then the same model. Uncertainty is a
bootstrap over whole years (200 resamples) of those pooled covariances, which
keeps each year's spatial and seasonal structure intact. Pooled numbers are not
directly comparable with the per-cell medians of multivariate.py; compare them
with each other.

Anomalies are standardised once over 2001-2020 (per cell and calendar month), so
a product whose variability drifts shows up in the windows instead of being
re-normalised away.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.cluster.hierarchy import cophenet, dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import multivariate as mv  # noqa: E402
from src import metrics as M  # noqa: E402
from src.adata import G_MOST, G_NEVER, G_SOME, countries, gauge_class, koppen_major, load  # noqa: E402

PAIR = mv.SCENARIOS[mv.PRIMARY]["pair"]
P = mv.P
CONTINENTS = ["Africa", "Asia", "Europe", "North America", "South America", "Oceania"]
N_BOOT = 200
RNG = np.random.default_rng(0)


# ------------------------------------------------------------------ aggregation

def aggregator(cells: dict, keys: np.ndarray) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Sparse (G, N) area-weight matrix grouping cells by `keys` (-1 = drop)."""
    w = (cells["area"] * cells["landfrac"]).astype(np.float64)
    ok = keys >= 0
    uniq, g = np.unique(keys[ok], return_inverse=True)
    A = sparse.csr_matrix((w[ok], (g, np.nonzero(ok)[0])), shape=(len(uniq), len(keys)))
    return A, uniq


def aggregate(x: np.ndarray, A: sparse.csr_matrix, min_cov: float = 0.8) -> np.ndarray:
    """Area-weighted mean of (T, N) over groups -> (T, G); NaN where less than
    `min_cov` of the group's weight has data that month."""
    fin = np.isfinite(x)
    num = (A @ np.where(fin, x, 0.0).T).T
    cov = (A @ fin.T.astype(np.float64)).T
    tot = np.asarray(A.sum(1)).ravel()
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(cov >= min_cov * tot, num / cov, np.nan).astype(np.float32)


def admin1_keys(cells: dict) -> tuple[np.ndarray, pd.DataFrame]:
    """Majority Natural Earth admin-1 unit per 0.5 deg cell (same method as
    src.adata.countries)."""
    import geopandas as gpd
    import rasterio.features as rf
    from rasterio.transform import from_origin

    from src.fetch import download

    out = mv.CACHE / "cell_admin1.parquet"
    if not out.exists():
        zp = download("https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_1_states_provinces.zip",
                      mv.CACHE / "ne_10m_admin_1_states_provinces.zip")
        gdf = gpd.read_file(zp)[["adm1_code", "name", "admin", "geometry"]].reset_index(drop=True)
        idx = rf.rasterize(((g, i + 1) for i, g in enumerate(gdf.geometry)), out_shape=(3600, 7200),
                           transform=from_origin(-180, 90, 0.05, 0.05), dtype="int32")[::-1]
        blocks = idx.reshape(360, 10, 720, 10).transpose(0, 2, 1, 3).reshape(360, 720, 100)
        b = blocks[cells["iy"], cells["ix"]]
        code = np.zeros(b.shape[0], dtype=np.int32)
        for k in range(b.shape[0]):
            row = b[k][b[k] > 0]
            if row.size:
                code[k] = np.bincount(row).argmax()
        df = pd.DataFrame({"unit": code - 1})
        meta = gdf.drop(columns="geometry")
        df["adm1_code"] = np.where(code > 0, meta.adm1_code.reindex(code - 1).values, "")
        df["name"] = np.where(code > 0, meta.name.reindex(code - 1).values, "")
        df["admin"] = np.where(code > 0, meta.admin.reindex(code - 1).values, "")
        df.to_parquet(out)
    df = pd.read_parquet(out)
    return df.unit.values.astype(np.int64), df


def group_label(labels: np.ndarray, A: sparse.csr_matrix) -> np.ndarray:
    """Majority (by weight) label of each group."""
    A = A.tocoo()
    df = pd.DataFrame({"g": A.row, "lab": labels[A.col], "w": A.data})
    df = df[df.lab != ""]
    top = df.groupby(["g", "lab"]).w.sum().reset_index().sort_values("w").groupby("g").tail(1)
    out = np.full(A.shape[0], "", dtype=object)
    out[top.g.values] = top.lab.values
    return out


def zstack(src: dict, months: np.ndarray, ref_months: int) -> np.ndarray:
    """Rain-relevant standardised anomalies, jointly masked: (T, N, K)."""
    clim = M.monthly_clim(src[mv.REF_CLIM], months)
    rel = clim[months - 1] >= mv.RAIN_RELEVANT * ref_months
    return mv.stack(src, rel, months)


def rho2_table(z: np.ndarray, fill_bad: np.ndarray, groups: dict) -> dict:
    C, _ = mv.cell_cov(z)
    r = mv.skill(C, pair=PAIR)["rho2"]
    r["chirps_v3"] = np.where(fill_bad, np.nan, r["chirps_v3"])
    return {g: {p: float(np.nanmedian(r[p][sel])) for p in P} | {"_units": int((sel & np.isfinite(r["era5"])).sum())}
            for g, sel in groups.items() if sel.sum() >= 5}


# ------------------------------------------------------------------ pooled + bootstrap

def year_stats(z: np.ndarray, years: np.ndarray, sel_cells: np.ndarray, sel_t: np.ndarray | None = None) -> dict:
    """Per-year sufficient statistics (n, sum, cross-products) of the pooled
    (cell, month) anomaly vectors of a region."""
    out = {}
    for y in np.unique(years):
        t = years == y
        if sel_t is not None:
            t = t & sel_t
        X = z[t][:, sel_cells].reshape(-1, z.shape[-1])
        X = X[np.isfinite(X).all(1)].astype(np.float64)
        out[int(y)] = (len(X), X.sum(0), X.T @ X)
    return out


def pooled_rho2(stats: dict, yrs: list) -> dict:
    n = sum(stats[y][0] for y in yrs)
    if n < 50:
        return {p: np.nan for p in P}
    s = sum(stats[y][1] for y in yrs)
    ss = sum(stats[y][2] for y in yrs)
    mu = s / n
    C = (ss - n * np.outer(mu, mu)) / (n - 1)
    r = mv.skill(C[None], pair=PAIR)["rho2"]
    return {p: float(r[p][0]) for p in P}


def pooled_with_ci(stats: dict, yrs: list) -> dict:
    est = pooled_rho2(stats, yrs)
    boots = np.array([[v for v in pooled_rho2(stats, list(RNG.choice(yrs, len(yrs)))).values()] for _ in range(N_BOOT)])
    lo, hi = np.nanpercentile(boots, [5, 95], axis=0)
    return {p: {"est": est[p], "lo": float(lo[i]), "hi": float(hi[i])} for i, p in enumerate(P)}


# ------------------------------------------------------------------ main

def main() -> None:
    d = load()
    cells = d.cells
    tb = d.window(*mv.WINDOW)
    mo, yr = d.months[tb], d.years[tb]
    ctry = countries(cells, mv.CACHE)
    cont = ctry.continent.values.astype(object)
    gclass = gauge_class(d.extra["gpcc_numgauge"], tb)
    kg = koppen_major(cells)
    with np.errstate(all="ignore"):
        fill = np.nanmean(d.extra["chirps3_fill"][tb], 0)
    lat = cells["lat"]
    base = {p: d.precip[p][tb] for p in P}
    in_dom = np.abs(lat) < 50  # the common domain (CHIRPS v2 limit)

    # ---------------------------------------------- 1. aggregation scale
    scale_out = {}
    a1, a1_df = admin1_keys(cells)
    iy64, ix64 = cells["iy"].astype(np.int64), cells["ix"].astype(np.int64)  # int16 overflows in the keys below
    spatial = {
        "0.5°": np.where(in_dom, np.arange(len(lat)), -1),
        "1°": np.where(in_dom, (iy64 // 2) * 1000 + ix64 // 2, -1),
        "2°": np.where(in_dom, (iy64 // 4) * 1000 + ix64 // 4, -1),
        "admin-1": np.where(in_dom & (a1 >= 0), a1, -1),
    }
    for sname, keys in spatial.items():
        A, _ = aggregator(cells, keys)
        agg = {p: aggregate(base[p], A) for p in P}
        fill_g = aggregate(fill[None], A, min_cov=0.0)[0]
        fill_bad = ~(fill_g <= mv.FILL_MAX)
        gc = group_label(cont, A)
        groups = {"All": np.ones(A.shape[0], bool)} | {c: gc == c for c in CONTINENTS}
        scale_out[sname] = {}
        for months_n in (1, 3, 6):
            src = agg if months_n == 1 else {p: pd.DataFrame(v).rolling(months_n).sum().to_numpy().astype(np.float32)
                                             for p, v in agg.items()}
            z = zstack(src, mo, months_n)
            scale_out[sname][f"{months_n}-month"] = rho2_table(z, fill_bad, groups)
            print(f"[scales] {sname} {months_n}-month: {A.shape[0]} units, {scale_out[sname][f'{months_n}-month']['All']['_units']} scored")

    # ---------------------------------------------- 2 + 3. pooled: time windows, season
    # restrict to the common domain from here on (memory)
    base = {p: v[:, in_dom] for p, v in base.items()}
    cont, gclass, kg = cont[in_dom], gclass[in_dom], kg[in_dom]
    z1 = zstack(base, mo, 1)
    z3 = zstack({p: pd.DataFrame(d.precip[p]).rolling(3).sum().to_numpy()[tb][:, in_dom].astype(np.float32) for p in P}, mo, 3)
    alln = np.ones(len(cont), bool)
    regions = {"All": alln, "Africa": cont == "Africa", "Asia": cont == "Asia",
               "South America": cont == "South America", G_NEVER: gclass == G_NEVER, G_MOST: gclass == G_MOST}
    windows = [(y, y + 4) for y in range(mv.WINDOW[0], mv.WINDOW[1] - 3)]
    time_out, season_out = {}, {}
    clim = M.monthly_clim(base[mv.REF_CLIM], mo)  # (12, N), common domain
    peak = np.nanargmax(np.nan_to_num(clim, nan=-1), 0) + 1  # wettest calendar month per cell
    pos = ((mo[:, None] - peak[None, :] + 6) % 12) - 6  # (T, N): -6..5, 0 = peak month
    for scale, z in (("monthly", z1), ("3-month", z3)):
        time_out[scale], season_out[scale] = {}, {}
        for reg, sel in regions.items():
            st = year_stats(z, yr, sel)
            time_out[scale][reg] = [{"window": f"{a}–{b}", **pooled_with_ci(st, list(range(a, b + 1)))} for a, b in windows]
            if scale == "monthly":
                season_out[scale][reg] = []
                zr, posr = z[:, sel], pos[:, sel]
                for k in range(-4, 5):
                    zk = np.where((posr == k)[..., None], zr, np.nan)
                    stk = year_stats(zk, yr, np.ones(zr.shape[1], bool))
                    season_out[scale][reg].append({"pos": k, **pooled_with_ci(stk, sorted(stk))})
            print(f"[scales] pooled {scale} {reg}")

    # ---------------------------------------------- 4. grouping by region
    C, _ = mv.cell_cov(z1)
    R = mv.corr_from_cov(C)
    resid = R.copy()
    r = mv.skill(C, pair=PAIR)["rho2"]
    rho = np.sqrt(np.stack([r[p] for p in P], -1))
    resid = R - rho[:, :, None] * rho[:, None, :]
    resid[:, np.arange(len(P)), np.arange(len(P))] = np.nan

    def tree(sel):
        Rm = np.nanmedian(R[sel], 0)
        dist = np.clip(1 - (Rm + Rm.T) / 2, 0, None)
        np.fill_diagonal(dist, 0)
        cd = squareform(dist, checks=False)
        Z = linkage(cd, method="average")
        return Z, cd

    Zg, cdg = tree(np.isfinite(R).all((1, 2)))
    group_regions = {"All": alln} | {c: cont == c for c in CONTINENTS} | \
        {f"Köppen {k}": kg == k for k in ["A Tropical", "B Arid", "C Temperate", "D Cold"]} | \
        {g: gclass == g for g in (G_NEVER, G_SOME, G_MOST)}
    groups_out = {}
    for reg, sel in group_regions.items():
        sel = sel & np.isfinite(R).all((1, 2))
        if sel.sum() < 200:
            continue
        Z, cd = tree(sel)
        coph_vs_global = float(np.corrcoef(cophenet(Z), cophenet(Zg))[0, 1])
        fam = fcluster(Z, t=4, criterion="maxclust")
        Rres = np.nanmedian(resid[sel], 0)
        partner = {p: P[int(np.nanargmax(np.where(np.arange(len(P)) == i, -9, Rres[i])))] for i, p in enumerate(P)}
        partner_r = {p: float(np.nanmax(np.where(np.arange(len(P)) == i, -9, Rres[i]))) for i, p in enumerate(P)}
        groups_out[reg] = {"cells": int(sel.sum()), "coph_vs_global": coph_vs_global,
                           "families": {p: int(f) for p, f in zip(P, fam)},
                           "order": dendrogram(Z, no_plot=True, labels=P)["ivl"], "linkage": Z,
                           "closest_partner": partner, "closest_partner_resid": partner_r}

    mv.dump("scales", {"products": P, "nrt": mv.NRT, "scale": scale_out, "time": time_out, "season": season_out,
                       "groups": groups_out, "windows": [f"{a}–{b}" for a, b in windows], "n_boot": N_BOOT,
                       "admin1_units": int(len(np.unique(a1[a1 >= 0])))})


if __name__ == "__main__":
    main()
