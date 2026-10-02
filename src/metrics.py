"""Vectorised metrics over (time, cell) arrays. NaN = missing throughout.

Conventions: `x`, `y` are float arrays of shape (T, N) (months x cells) on the
same time axis; `months` is an int array (T,) of calendar months 1..12; `w` are
cell area weights (N,).
"""

import numpy as np
import pandas as pd
from scipy import stats


# ------------------------------------------------------------------ basics

def paired(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    m = np.isfinite(x) & np.isfinite(y)
    return np.where(m, x, np.nan), np.where(m, y, np.nan)


def monthly_clim(x: np.ndarray, months: np.ndarray, min_frac: float = 0.8) -> np.ndarray:
    """(12, N) calendar-month means; NaN where < min_frac of years are present."""
    out = np.full((12, x.shape[1]), np.nan, dtype=np.float64)
    for m in range(1, 13):
        xm = x[months == m]
        ok = np.isfinite(xm).mean(0) >= min_frac
        with np.errstate(invalid="ignore"):
            out[m - 1] = np.where(ok, np.nanmean(xm, 0), np.nan)
    return out


def annual_clim(x: np.ndarray, months: np.ndarray) -> np.ndarray:
    """Mean annual total (N,) = sum of the 12 monthly climatologies."""
    c = monthly_clim(x, months)
    return np.where(np.isfinite(c).all(0), c.sum(0), np.nan)


def ratio_of_sums(x: np.ndarray, y: np.ndarray, axis=0) -> np.ndarray:
    """Sum(x)/Sum(y) over paired months (per cell for axis=0)."""
    xp, yp = paired(x, y)
    sx, sy = np.nansum(xp, axis), np.nansum(yp, axis)
    n = np.isfinite(xp).sum(axis)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where((n > 0) & (sy > 0), sx / sy, np.nan)


def regional_ratio(x, y, w, sel) -> float:
    """Area-weighted ratio of sums over cells `sel` (bool (N,)), paired months."""
    xp, yp = paired(x[:, sel], y[:, sel])
    sx = np.nansum(np.nansum(xp, 0) * w[sel])
    sy = np.nansum(np.nansum(yp, 0) * w[sel])
    return float(sx / sy) if sy > 0 else np.nan


def std_anom(x: np.ndarray, months: np.ndarray) -> np.ndarray:
    """Standardised anomalies per calendar month (own climatology)."""
    z = np.full(x.shape, np.nan, dtype=np.float32)
    for m in range(1, 13):
        sel = months == m
        mu = np.nanmean(x[sel], 0)
        sd = np.nanstd(x[sel], 0)
        with np.errstate(invalid="ignore", divide="ignore"):
            z[sel] = np.where(sd > 0, (x[sel] - mu) / sd, np.nan)
    return z


# ------------------------------------------------------------------ correlation

def _rank_cols(a: np.ndarray) -> np.ndarray:
    return pd.DataFrame(a).rank(axis=0, method="average").to_numpy()


def corr_cols(a: np.ndarray, b: np.ndarray, min_n: int = 24) -> np.ndarray:
    """Pearson correlation per column over paired finite rows."""
    a, b = paired(a, b)
    n = np.isfinite(a).sum(0)
    am = a - np.nanmean(a, 0)
    bm = b - np.nanmean(b, 0)
    num = np.nansum(am * bm, 0)
    den = np.sqrt(np.nansum(am**2, 0) * np.nansum(bm**2, 0))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = num / den
    return np.where(n >= min_n, r, np.nan)


def spearman_cols(a: np.ndarray, b: np.ndarray, min_n: int = 24) -> np.ndarray:
    a, b = paired(a, b)
    return corr_cols(_rank_cols(a), _rank_cols(b), min_n)


# ------------------------------------------------------------------ triple collocation

def etc_rho2(x, y, z, min_n: int = 36, clip: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extended triple collocation (McColl et al. 2014): squared correlation of
    each of x, y, z with the unknown truth, per column. Assumes mutually
    independent, zero-mean errors and a linear relation to the truth."""
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = (np.where(m, a, np.nan) for a in (x, y, z))
    n = m.sum(0)

    def cov(a, b):
        return np.nanmean((a - np.nanmean(a, 0)) * (b - np.nanmean(b, 0)), 0)

    qxx, qyy, qzz = cov(x, x), cov(y, y), cov(z, z)
    qxy, qxz, qyz = cov(x, y), cov(x, z), cov(y, z)
    with np.errstate(invalid="ignore", divide="ignore"):
        rx = qxy * qxz / (qxx * qyz)
        ry = qxy * qyz / (qyy * qxz)
        rz = qxz * qyz / (qzz * qxy)
    bad = n < min_n
    if not clip:
        return tuple(np.where(bad, np.nan, r) for r in (rx, ry, rz))
    # sampling noise pushes estimates outside [0, 1] where the true value is near
    # an end: clip to the bound (dropping them would bias the medians)
    return tuple(np.where(bad, np.nan, np.clip(r, 0.0, 1.0)) for r in (rx, ry, rz))


# ------------------------------------------------------------------ terciles

def tercile_cat(x: np.ndarray, months: np.ndarray, ref: np.ndarray | None = None) -> np.ndarray:
    """Categorise each value into 0 = dry, 1 = normal, 2 = wet tercile of its
    calendar month. Thresholds come from `ref` (same shape) if given — the
    common-climatology variant — else from `x` itself. NaN stays NaN."""
    src = x if ref is None else ref
    cat = np.full(x.shape, np.nan, dtype=np.float32)
    for m in range(1, 13):
        sel = months == m
        lo = np.nanpercentile(src[sel], 100 / 3, axis=0)
        hi = np.nanpercentile(src[sel], 200 / 3, axis=0)
        v = x[sel]
        c = np.where(v <= lo, 0.0, np.where(v > hi, 2.0, 1.0))
        cat[sel] = np.where(np.isfinite(v), c, np.nan)
    return cat


def dry_skill(cat_p: np.ndarray, cat_r: np.ndarray) -> dict[str, float]:
    """Pooled dry-tercile detection of product p against reference r.
    Returns hit rate, false-alarm rate (POFD), Peirce skill score and Cohen's
    kappa for the 3-category table."""
    m = np.isfinite(cat_p) & np.isfinite(cat_r)
    p, r = cat_p[m], cat_r[m]
    if p.size < 30:
        return {"n": int(p.size), "hit": np.nan, "pofd": np.nan, "far": np.nan, "base": np.nan, "pss": np.nan, "kappa": np.nan}
    obs_dry, fc_dry = r == 0, p == 0
    hit = (fc_dry & obs_dry).sum() / max(obs_dry.sum(), 1)
    pofd = (fc_dry & ~obs_dry).sum() / max((~obs_dry).sum(), 1)
    # false-alarm RATIO: share of the product's dry calls the reference does not call dry
    far = (fc_dry & ~obs_dry).sum() / max(fc_dry.sum(), 1)
    po = (p == r).mean()
    pe = sum((p == k).mean() * (r == k).mean() for k in (0, 1, 2))
    kappa = (po - pe) / (1 - pe) if pe < 1 else np.nan
    return {"n": int(p.size), "hit": float(hit), "pofd": float(pofd), "far": float(far), "base": float(obs_dry.mean()),
            "pss": float(hit - pofd), "kappa": float(kappa)}


def dry_skill_cells(cat_p: np.ndarray, cat_r: np.ndarray, min_n: int = 30) -> np.ndarray:
    """Per-cell Peirce skill score for dry-tercile detection."""
    m = np.isfinite(cat_p) & np.isfinite(cat_r)
    obs_dry = (cat_r == 0) & m
    fc_dry = (cat_p == 0) & m
    n_obs = obs_dry.sum(0)
    n_not = (m & ~obs_dry).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        h = (fc_dry & obs_dry).sum(0) / n_obs
        f = (fc_dry & ~obs_dry & m).sum(0) / n_not
    return np.where(m.sum(0) >= min_n, h - f, np.nan)


# ------------------------------------------------------------------ trends

def annual_totals(x: np.ndarray, years: np.ndarray, yrs: np.ndarray) -> np.ndarray:
    """(len(yrs), N) annual totals; NaN unless all 12 months present."""
    out = np.full((len(yrs), x.shape[1]), np.nan)
    for i, y in enumerate(yrs):
        xy = x[years == y]
        if xy.shape[0] == 12:
            out[i] = np.where(np.isfinite(xy).all(0), xy.sum(0), np.nan)
    return out


def sen_mk(y: np.ndarray, min_frac: float = 0.9) -> dict[str, np.ndarray]:
    """Sen slope (units/yr) and autocorrelation-corrected Mann-Kendall p-value per
    column of y (n_years, N). Columns with < min_frac finite years -> NaN.
    Missing years are dropped pairwise (slopes) / treated as ties-free gaps."""
    n, N = y.shape
    t = np.arange(n, dtype=np.float64)
    ok = np.isfinite(y).mean(0) >= min_frac
    i, j = np.triu_indices(n, 1)
    dy = y[j] - y[i]  # (pairs, N)
    dt = (j - i)[:, None].astype(np.float64)
    slope = np.nanmedian(dy / dt, axis=0)
    s = np.nansum(np.sign(dy), 0)
    nn = np.isfinite(y).sum(0).astype(np.float64)
    var_s = nn * (nn - 1) * (2 * nn + 5) / 18.0
    # Autocorrelation correction (Yue & Wang 2004 style): AR(1) effective sample
    # size from the lag-1 autocorrelation of the Sen-detrended series, with the
    # small-sample bias correction r1* = (n r1 + 1) / (n - 4). A significance
    # filter on r1 (Hamed & Rao) silently drops the correction at n ~ 25-40,
    # which doubles the false-positive rate under modest persistence.
    resid = y - slope[None, :] * t[:, None]
    rm = resid - np.nanmean(resid, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        r1 = np.nansum(rm[:-1] * rm[1:], 0) / np.nansum(rm**2, 0)
    r1 = np.clip((nn * r1 + 1) / np.maximum(nn - 4, 1), -0.9, 0.9)
    k = np.arange(1, n)[:, None]
    eta = 1 + 2 * np.nansum((1 - k / nn[None, :]) * np.where(r1 > 0, r1, 0.0)[None, :] ** k, 0)
    var_c = var_s * eta
    z = np.where(s > 0, (s - 1) / np.sqrt(var_c), np.where(s < 0, (s + 1) / np.sqrt(var_c), 0.0))
    p = 2 * stats.norm.sf(np.abs(z))
    mean = np.nanmean(y, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        pct_dec = 100 * 10 * slope / mean
    return {
        "slope": np.where(ok, slope, np.nan),
        "pct_dec": np.where(ok & (mean > 0), pct_dec, np.nan),
        "p": np.where(ok, p, np.nan),
        "mean": np.where(ok, mean, np.nan),
    }


def bh_fdr(p: np.ndarray, q: float = 0.1) -> np.ndarray:
    """Benjamini-Hochberg: bool mask of discoveries at FDR level q (Wilks 2016
    recommends q = 2 x the intended global alpha, i.e. 0.1 for 0.05)."""
    out = np.zeros(p.shape, dtype=bool)
    fin = np.isfinite(p)
    ps = p[fin]
    if ps.size == 0:
        return out
    order = np.argsort(ps)
    thresh = q * np.arange(1, ps.size + 1) / ps.size
    passed = ps[order] <= thresh
    if passed.any():
        k = np.nonzero(passed)[0].max()
        cutoff = ps[order][k]
        out[fin] = ps <= cutoff
    return out


def wet_season_start(clim: np.ndarray) -> np.ndarray:
    """Start month (1..12) of the wettest 3-consecutive-month window per cell
    from a (12, N) climatology (wraps around the year)."""
    ext = np.concatenate([clim, clim[:2]], 0)
    s3 = ext[:-2] + ext[1:-1] + ext[2:]
    s3 = np.where(np.isfinite(s3), s3, -np.inf)
    return np.argmax(s3, 0) + 1


def season_totals(x: np.ndarray, years: np.ndarray, months: np.ndarray, start: np.ndarray, yrs: np.ndarray) -> np.ndarray:
    """(len(yrs), N) totals of the 3-month season starting at `start` (per
    cell); a season wrapping the year end is assigned to the year it ends in."""
    out = np.full((len(yrs), x.shape[1]), np.nan)
    tindex = {(int(y), int(m)): k for k, (y, m) in enumerate(zip(years, months))}
    for s in range(1, 13):
        cols = np.nonzero(start == s)[0]
        if cols.size == 0:
            continue
        for i, yend in enumerate(yrs):
            ms = [(s + d - 1) % 12 + 1 for d in range(3)]
            ys = [yend - (1 if (s + 2 > 12 and m >= s) else 0) for m in ms]
            ks = [tindex.get((y, m)) for y, m in zip(ys, ms)]
            if any(k is None for k in ks):
                continue
            v = x[ks][:, cols]
            out[i, cols] = np.where(np.isfinite(v).all(0), v.sum(0), np.nan)
    return out
