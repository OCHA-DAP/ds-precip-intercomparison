"""Render the map figures (PNG) from cache/fields.npz into site/report/fig/.

Colour: drier = red, wetter = blue, neutral grey midpoint (diverging); one-hue
blue ramp for magnitudes (sequential). Land cells without a value are light
grey so "no data" never reads as "zero".
"""

import json
import sys
from pathlib import Path

import logging

import matplotlib

matplotlib.use("Agg")
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.adata import LABEL, ROOT  # noqa: E402

CACHE = ROOT / "cache"
FIG = ROOT / "site" / "report" / "fig"
DATA = ROOT / "site" / "report" / "data"

INK = "#1f2324"
MUTED = "#5e6a6b"
NODATA = "#e4e6e6"
LAT0, LAT1 = -60, 80  # map extent

# Diverging: red (drier / decreasing) <-> grey <-> blue (wetter / increasing)
DIV = LinearSegmentedColormap.from_list(
    "div", ["#8f1d1d", "#d03b3b", "#ef9a8a", "#f0efec", "#86b6ef", "#2a78d6", "#104281"]
)
SEQ = LinearSegmentedColormap.from_list("seq", ["#eef4fc", "#b7d3f6", "#6da7ec", "#2a78d6", "#1c5cab", "#0d366b"])

plt.rcParams.update({
    "font.family": ["Roboto", "DejaVu Sans"],
    "font.size": 8,
    "axes.edgecolor": "#c9cfcf",
    "axes.labelcolor": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "text.color": INK,
})


def grid(fields, key):
    a = np.full((360, 720), np.nan, dtype=np.float32)
    a[fields["iy"], fields["ix"]] = fields[key]
    return a


def landmask(fields):
    a = np.zeros((360, 720), bool)
    a[fields["iy"], fields["ix"]] = True
    return a


def draw(ax, arr, land, cmap, norm, title, alpha_mask=None):
    r0, r1 = (LAT0 + 90) * 2, (LAT1 + 90) * 2
    ext = (-180, 180, LAT0, LAT1)
    bg = np.where(land, 1.0, np.nan)[r0:r1]
    ax.imshow(bg, origin="lower", extent=ext, cmap=ListedColormap([NODATA]), interpolation="nearest")
    sub = arr[r0:r1]
    if alpha_mask is not None:
        faint = np.where(alpha_mask[r0:r1], np.nan, sub)
        strong = np.where(alpha_mask[r0:r1], sub, np.nan)
        ax.imshow(faint, origin="lower", extent=ext, cmap=cmap, norm=norm, interpolation="nearest", alpha=0.35)
        im = ax.imshow(strong, origin="lower", extent=ext, cmap=cmap, norm=norm, interpolation="nearest")
    else:
        im = ax.imshow(sub, origin="lower", extent=ext, cmap=cmap, norm=norm, interpolation="nearest")
    ax.set_title(title, fontsize=8.5, loc="left", color=INK, pad=3)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_aspect("equal")
    return im


def small_multiples(fields, keys, titles, cmap, bounds, cbar_label, out, ncol=4, sig_keys=None, extend="both",
                    tick_labels=None, tick_mid=False):
    land = landmask(fields)
    n = len(keys)
    nrow = int(np.ceil(n / ncol))
    if ncol == 1:
        fig, axes = plt.subplots(1, 1, figsize=(8.0, 3.9), squeeze=False)
    else:
        fig, axes = plt.subplots(nrow, ncol, figsize=(2.75 * ncol, 1.38 * nrow + 0.9), squeeze=False)
    norm = BoundaryNorm(bounds, cmap.N, extend=extend)
    im = None
    for ax, k, t, sk in zip(axes.flat, keys, titles, sig_keys or [None] * n):
        if k not in fields:
            ax.set_visible(False)
            continue
        mask = None
        if sk and sk in fields:
            mask = grid(fields, sk) > 0.5
        im = draw(ax, grid(fields, k), land, cmap, norm, t, alpha_mask=mask)
    for ax in list(axes.flat)[n:]:
        ax.set_visible(False)
    bottom = 0.2 if ncol == 1 else 0.9 / (1.38 * nrow + 0.9)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.93 if ncol == 1 else 0.96, bottom=bottom, wspace=0.03, hspace=0.2)
    if im is not None:
        cax = fig.add_axes([0.25, bottom * 0.55, 0.5, bottom * 0.13])
        cb = fig.colorbar(im, cax=cax, orientation="horizontal", extend=extend)
        cb.set_label(cbar_label, fontsize=8, color=INK)
        cb.outline.set_visible(False)
        if tick_labels:
            ticks = [(a + b) / 2 for a, b in zip(bounds[:-1], bounds[1:])] if tick_mid else bounds
            cb.set_ticks(ticks)
            cb.set_ticklabels(tick_labels)
        cb.ax.tick_params(labelsize=7, length=2)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / out, dpi=170)
    plt.close(fig)
    print(f"[figures] {out}")


def main() -> None:
    f = dict(np.load(CACHE / "fields.npz"))
    meta = json.loads((DATA / "meta.json").read_text())
    P = [p["key"] for p in meta["products"]]
    T = {p: ("ASAP blend (CHIRPS v2 / ERA5)" if p == "asap" else LABEL[p]) for p in P}

    # gauge density
    # never-gauged gets its own warm colour: it is the class that matters most
    # and must not be confused with the grey "no data"
    gcmap = ListedColormap(["#f0a58a", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#184f95"])
    small_multiples(f, ["gpcc_gauged_frac"], ["Share of months 2001–2020 with at least one GPCC gauge in the 0.5° cell"],
                    gcmap, [0, 1e-6, 0.1, 0.25, 0.5, 0.75, 1.0], "share of months gauged", "gauges.png", ncol=1,
                    extend="neither", tick_labels=["never", "<10%", "10–25%", "25–50%", "50–75%", "≥75%"], tick_mid=True)

    small_multiples(f, [f"clim_{p}" for p in P], [T[p] for p in P], SEQ,
                    [0, 100, 250, 500, 750, 1000, 1500, 2000, 3000, 5000], "mean annual precipitation 2001–2020 (mm)",
                    "climatology.png", extend="max")

    ratio_b = [0.5, 0.67, 0.8, 0.9, 0.95, 1.05, 1.11, 1.25, 1.5, 2.0]
    ratio_t = ["½", "⅔", "0.8", "0.9", "0.95", "1.05", "1.11", "1.25", "1.5", "2"]
    nb = [p for p in P if p != "gpcc_full"]
    small_multiples(f, [f"bias_{p}" for p in nb], [T[p] for p in nb], DIV, ratio_b,
                    "ratio to GPCC Full Data, 2001–2020 (sum of product / sum of GPCC); red = drier, blue = wetter",
                    "bias_vs_gpcc.png", tick_labels=ratio_t)

    small_multiples(f, [f"spearman_{p}" for p in nb], [T[p] for p in nb], SEQ,
                    [0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                    "Spearman correlation of monthly standardised anomalies with GPCC, 2001–2020", "spearman_vs_gpcc.png",
                    extend="min")

    tc = [p for p in P if f"tc_{p}" in f]
    small_multiples(f, [f"tc_{p}" for p in tc], [T[p] for p in tc], SEQ,
                    [0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                    "triple collocation: squared correlation with the unknown truth (monthly anomalies, 2001–2020)",
                    "tc_rho2.png", extend="neither")

    small_multiples(f, [f"pss3_{p}" for p in nb], [T[p] for p in nb], SEQ,
                    [0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0],
                    "dry-tercile detection vs GPCC, 3-month totals (Peirce skill score)", "pss3_vs_gpcc.png", extend="min")

    trend_b = [-30, -20, -10, -5, -2, 2, 5, 10, 20, 30]
    for w, yrs in (("long", "1983–2020"), ("recent", "2001–2025")):
        keys = [p for p in P if f"trend_{w}_{p}" in f]
        small_multiples(f, [f"trend_{w}_{p}" for p in keys], [T[p] for p in keys], DIV, trend_b,
                        f"trend in annual total {yrs} (% of mean per decade); faded = not significant (FDR 10%)",
                        f"trend_{w}.png", sig_keys=[f"trendsig_{w}_{p}" for p in keys])
        if f"agree_pos_{w}" in f:
            n = f[f"agree_n_{w}"]
            with np.errstate(invalid="ignore", divide="ignore"):
                f[f"agree_frac_{w}"] = np.where(n >= 3, f[f"agree_pos_{w}"] / n, np.nan)
            small_multiples(f, [f"agree_frac_{w}"], [f"Share of products with a positive trend, {yrs}"], DIV,
                            [0, 0.1, 0.25, 0.4, 0.6, 0.75, 0.9, 1.0], "share of products with an increasing annual total",
                            f"trend_agreement_{w}.png", ncol=1, extend="neither",
                            tick_labels=["0", "0.1", "0.25", "0.4", "0.6", "0.75", "0.9", "1"])

    dk = [k for k in f if k.startswith("difftrend_")]
    if dk:
        titles = [k.replace("difftrend_long_", "1983–2020: ").replace("difftrend_recent_", "2001–2025: ")
                  .replace("_minus_", " − ") for k in dk]
        short = {"gpcc_full": "GPCC", "gpcc_monitoring": "GPCC Mon.", "cru": "CRU", "cpc": "CPC", "precl": "PREC/L",
                 "imerg_final": "IMERG Final", "imerg_late": "IMERG Late", "chirps_v2": "CHIRPS v2", "chirp_v2": "CHIRP v2",
                 "chirps_v3": "CHIRPS v3", "chirp_v3": "CHIRP v3", "era5": "ERA5"}
        for p in sorted(short, key=len, reverse=True):
            titles = [t.replace(p, short[p]) for t in titles]
        small_multiples(f, dk, titles, DIV, trend_b,
                        "trend of the difference A − B (% of B's mean per decade)", "difference_trends.png", ncol=3)

    if "v3_over_v2" in f:
        small_multiples(f, ["v3_over_v2", "lw_annual"],
                        ["CHIRPS v3 / CHIRPS v2, 2001–2020", "Legates–Willmott undercatch factor (annual, GPCC-weighted)"],
                        DIV, ratio_b, "ratio", "undercatch_check.png", ncol=2, tick_labels=ratio_t)


if __name__ == "__main__":
    main()
