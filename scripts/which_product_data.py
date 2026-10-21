"""Public data for the "Which rainfall product where" page (site/which-product/).

    python scripts/which_product_data.py   # -> site/which-product/data/*.json
    python scripts/site_blob.py upload     # parks it on blob with the report and explorer

Inputs: site/multivariate/data/ (multivariate*.py, validate_gauges.py, validate_3month.py)
and cache/tahmo/ (validate_tahmo.py). The page is public, so two things are left out:
- TAHMO: only the regional/tier scores and the selection-funnel counts go out. Station
  locations, countries and records in cache/tahmo/viz.json stay local (licence, ADR 0002).
- CRU TS (ODbL): its monthly values at the gauge cells are dropped from the time series
  (the chart never draws them); only derived statistics remain.
Per-cell arrays are packed as base64 little-endian integers, which the page decodes with
the typed array named in PACK; that takes fields.json from 7 MB to 2.7 MB.
"""

import base64
import json
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MV = ROOT / "site" / "multivariate" / "data"
TAHMO = ROOT / "cache" / "tahmo"
OUT = ROOT / "site" / "which-product" / "data"
# numpy dtype per packed array, matching the page's dec(..., Int16Array / Int8Array / Uint8Array)
PACK = {"iy": "<i2", "ix": "<i2", "regime": "i1", "best_nrt": "i1", "best_alt": "i1", "robust": "i1",
        "gauge": "i1", "farthest": "i1", "margin": "u1", "pc1": "u1", "pc1_core": "u1", "rho2": "u1",
        "far_dist": "u1", "gap": "u1", "disagree": "u1"}


def _read(p: Path):
    if not p.exists():
        raise SystemExit(f"missing input {p.relative_to(ROOT)}: run the script that writes it first (see docstring)")
    return json.loads(p.read_text())


def _pack(values: list, dtype: str) -> str:
    info = np.iinfo(dtype)
    a = np.asarray(values)
    if a.dtype.kind not in "iu" or a.min() < info.min or a.max() > info.max:
        raise ValueError(f"cannot pack values {a.dtype} [{a.min()}, {a.max()}] as {dtype}")
    return base64.b64encode(a.astype(dtype).tobytes()).decode()


def _pack_tree(x, dtype: str):
    return _pack(x, dtype) if isinstance(x, list) else {k: _pack_tree(v, dtype) for k, v in x.items()}


def build() -> dict:
    summary = _read(MV / "summary.json")
    for scale in summary["dependency"].values():
        for region in scale.values():
            region.pop("linkage")  # clustering order is already in "order"; the page doesn't use it

    fields = _read(MV / "fields.json")
    packed = {"n": len(fields["iy"])}
    packed.update({k: _pack_tree(v, PACK[k]) for k, v in fields.items()})

    outliers = _read(MV / "outliers.json")
    outliers["cell"] = {k: _pack(v, PACK[k]) for k, v in outliers["cell"].items()}

    gauges = _read(MV / "gauges.json")
    gauges.pop("per_cell")  # not drawn

    viz = _read(MV / "gauge_viz.json")
    for s in viz["series"]:
        s["products"].pop("cru")

    scores = _read(TAHMO / "scores.json")
    tahmo_scores = {k: scores[k] for k in ("products", "groups", "groups_3month")}
    tahmo_funnel = {"funnel": _read(TAHMO / "viz.json")["funnel"]}

    return {"summary": summary, "fields": packed, "outliers": outliers, "gauges": gauges, "gauge_viz": viz,
            "gauges_3month": _read(MV / "gauges_3month.json"), "tahmo_scores": tahmo_scores, "tahmo_funnel": tahmo_funnel}


def main() -> None:
    out = build()
    for name in ("tahmo_scores", "tahmo_funnel"):
        if re.search(r'"(lat|lon|latitude|longitude|station_id|stations|country)"\s*:|TA\d{5}', json.dumps(out[name])):
            raise SystemExit(f"{name} carries station-level TAHMO fields; refusing to write a public file")
    OUT.mkdir(parents=True, exist_ok=True)
    for name, d in out.items():
        (OUT / f"{name}.json").write_text(json.dumps(d, separators=(",", ":"), ensure_ascii=False))
    total = sum((OUT / f"{n}.json").stat().st_size for n in out)
    print(f"wrote {len(out)} files, {total / 1e6:.1f} MB, to {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
