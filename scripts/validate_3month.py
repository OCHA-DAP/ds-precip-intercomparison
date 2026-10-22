"""AMMA-CATCH gauges against the products on 3-month totals (the scale drought triggers use).

    python scripts/validate_gauges.py   # first: writes site/multivariate/data/gauge_viz.json
    python scripts/validate_3month.py   # -> site/multivariate/data/gauges_3month.json

Uses the same cells and gauge cell-means as the monthly check (gauge_viz.json series).
A 3-month total needs all three months present (gauge: three consecutive cell-months with
at least 3 complete stations). Anomalies are standardised per end-month (Jun-Aug against
other Jun-Aug), rain-relevant windows are those whose long-term GPCC Full total is at least
30 mm, and consecutive windows overlap, which the year-block bootstrap accounts for.
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import multivariate as mv  # noqa: E402
import validate_tahmo as vt  # noqa: E402
from src import metrics as M  # noqa: E402


def main() -> None:
    v = json.loads((mv.OUT / "gauge_viz.json").read_text())
    P = v["products"]
    times = np.array(v["times"])
    mo = np.array([int(t[5:7]) for t in times])
    yr = np.array([int(t[:4]) for t in times])
    rng = np.random.default_rng(0)
    by_site, by_month = {}, {}
    for c, s in zip(v["cells"], v["series"]):
        g = np.array(s["gauge"], dtype=float)
        g[g < 0] = np.nan
        pr = np.column_stack([np.where(np.array(s["products"][p]) < 0, np.nan, s["products"][p]) for p in P])
        g3, p3 = vt.three_month(g), vt.three_month(pr)
        clim3 = M.monthly_clim(p3[:, [P.index(mv.REF_CLIM)]], mo)[:, 0]
        z, ok = vt.zrows(g3, p3, mo, clim3[mo - 1] >= 3 * mv.RAIN_RELEVANT)
        if z is None:
            continue
        for key in ("All", c["site"]):
            by_site.setdefault(key, []).append((z, yr[ok], mo[ok]))
            for m in np.unique(mo[ok]):
                sel = mo[ok] == m
                by_month.setdefault(key, {}).setdefault(int(m), []).append((z[sel], yr[ok][sel], mo[ok][sel]))
    out = {"products": P, "groups": {}, "by_end_month": {}}
    for key, parts in by_site.items():
        X, Y, MO = (np.concatenate(a) for a in zip(*parts))
        out["groups"][key] = {"cells": len(parts), **vt.scores(X, Y, MO, P, rng)}
    for key, months in by_month.items():
        out["by_end_month"][key] = {}
        for m, parts in sorted(months.items()):
            X, Y, MO = (np.concatenate(a) for a in zip(*parts))
            out["by_end_month"][key][m] = vt.scores(X, Y, MO, P, rng) if len(X) >= 40 else None
    mv.dump("gauges_3month", out)


if __name__ == "__main__":
    main()
