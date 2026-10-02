"""Build one product (or aux layer) on the common 0.5 deg grid and upload it.

    python scripts/ingest.py chirps_v2                    # full span -> dev blob
    python scripts/ingest.py era5 --start 2010-01 --end 2010-12 --no-upload --out /tmp/x.nc

Products: see src.products.BUILDERS; aux layers: `static`, `chirps_v3_diag`.
Writes processed/grid05/<product>.nc (or processed/aux05/<name>.nc) on the DEV
blob (`projects` container). A fingerprint check (desert/rainforest/monsoon
cells, hemispheric seasonality) runs before upload and fails the task on a
units, orientation or time-offset error.
"""

import argparse
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.aux import build_chirps_v3_diag, build_static  # noqa: E402
from src.constants import INGEST_START, PROJECT_PREFIX  # noqa: E402
from src.products import BUILDERS, default_workdir  # noqa: E402

AUX = {"static", "chirps_v3_diag"}

# (name, lat, lon, annual-mean range mm/yr, peak months or None)
CHECKS = [
    ("Sahara (Algeria/Niger)", 23.25, 8.75, (0, 120), None),
    ("Amazon (Manaus)", -3.25, -60.25, (1500, 3800), None),
    ("Sahel (Niamey)", 13.25, 2.25, (250, 900), {7, 8, 9}),
    ("Zambia (Lusaka)", -15.25, 28.25, (500, 1300), {12, 1, 2}),
    ("SE England", 51.25, -0.25, (400, 1200), None),
]


def fingerprint(ds: xr.Dataset, product: str) -> list[str]:
    pr = ds["precip"]
    years = pd.DatetimeIndex(pr.time.values).year
    full = [y for y in np.unique(years) if (years == y).sum() == 12]
    window = [y for y in full if 2001 <= y <= 2020] or full
    if not window:
        return [f"{product}: no complete year to check"]
    sub = pr.sel(time=np.isin(years, window))
    problems = []
    for name, lat, lon, (lo, hi), peak in CHECKS:
        s = sub.sel(lat=lat, lon=lon)
        if s.isnull().mean() > 0.2:
            print(f"  [check] {name}: not covered, skipped")
            continue
        clim = s.groupby("time.month").mean()
        annual = float(clim.sum())
        pk = int(clim.month[int(clim.argmax())])
        ok = lo <= annual <= hi and (peak is None or pk in peak)
        print(f"  [check] {name}: {annual:.0f} mm/yr (expect {lo}-{hi}), peak month {pk}"
              f"{'' if peak is None else f' (expect {sorted(peak)})'} {'OK' if ok else 'FAIL'}")
        if not ok:
            problems.append(f"{name}: {annual:.0f} mm/yr, peak {pk}")
    return problems


def check_static(ds: xr.Dataset) -> list[str]:
    """An inverted or empty land mask would silently corrupt every comparison."""
    lf = ds["landfrac"]
    probs = []
    for name, lat, lon, lo, hi in (
        ("Antarctica", -80.25, 0.25, 0.95, 1.0),
        ("Sahara", 23.25, 8.75, 0.95, 1.0),
        ("central Pacific", 0.25, -150.25, 0.0, 0.05),
    ):
        v = float(lf.sel(lat=lat, lon=lon))
        print(f"  [check] landfrac {name}: {v:.2f}")
        if not lo <= v <= hi:
            probs.append(f"landfrac {name} = {v:.2f}")
    m = float(lf.mean())
    print(f"  [check] landfrac global unweighted mean {m:.3f}")
    if not 0.30 <= m <= 0.38:
        probs.append(f"landfrac mean {m:.3f}")
    agree = float(((ds["koppen"] > 0) == (lf >= 0.5)).mean())
    print(f"  [check] koppen>0 vs landfrac>=0.5 agreement {agree:.3f}")
    if agree < 0.95:
        probs.append(f"koppen/landfrac agreement {agree:.3f}")
    return probs


def write_nc(ds: xr.Dataset, path: Path) -> None:
    enc = {}
    for v, da in ds.data_vars.items():
        if da.dtype.kind == "f" and da.ndim >= 2:
            chunks = tuple(min(n, c) for n, c in zip(da.shape, (12, 90, 180)[-da.ndim:]))
            enc[v] = {"zlib": True, "complevel": 4, "chunksizes": chunks, "dtype": "float32"}
    path.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(path, encoding=enc, engine="netcdf4")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("product", choices=sorted(set(BUILDERS) | AUX))
    ap.add_argument("--start", default="")
    ap.add_argument("--end", default="", help="last month YYYY-MM (default: last complete month)")
    ap.add_argument("--subdir", default="", help="blob subfolder override, e.g. 'smoke' for test runs")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--no-upload", action="store_true")
    ap.add_argument("--no-check", action="store_true")
    ap.add_argument("--out", type=Path, help="local output path (default: workdir/<product>.nc)")
    args = ap.parse_args()

    if not args.end:
        args.end = str(pd.Timestamp.today().to_period("M") - 1)
    if not args.start:
        args.start = INGEST_START
    stage = os.environ.get("STAGE", "dev")
    if stage != "dev" and not args.no_upload:
        raise SystemExit("this project writes the DEV blob only (STAGE=dev)")
    workdir = default_workdir() / args.product
    workdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    print(f"[ingest] {args.product} {args.start}..{args.end} workdir={workdir}", flush=True)

    aux = args.product in AUX
    sub = args.subdir or ("aux05" if aux else "grid05")
    blob = f"{PROJECT_PREFIX}/processed/{sub}/{args.product}.nc"
    if args.product == "static":
        ds = build_static(workdir)
    elif args.product == "chirps_v3_diag":
        ds = build_chirps_v3_diag(args.start, args.end, workdir, workers=args.workers)
    else:
        ds = BUILDERS[args.product](args.start, args.end, workdir, workers=args.workers)
        ds.attrs["product"] = args.product
        print(f"[ingest] {args.product}: {ds.sizes['time']} months "
              f"{pd.Timestamp(ds.time.values[0]):%Y-%m}..{pd.Timestamp(ds.time.values[-1]):%Y-%m}", flush=True)

    ds.attrs["built"] = pd.Timestamp.now(tz="UTC").isoformat()
    out = args.out or (workdir / f"{args.product}.nc")
    write_nc(ds, out)  # write BEFORE checking: a failed check still leaves an artefact
    print(f"[ingest] wrote {out} ({out.stat().st_size / 1e6:.0f} MB) in {time.time() - t0:.0f}s", flush=True)

    problems = []
    if not args.no_check:
        if args.product == "static":
            problems = check_static(ds)
        elif args.product not in ("cru_tmp", "chirps_v3_diag"):
            problems = fingerprint(ds, args.product)
    if problems:
        blob = f"{PROJECT_PREFIX}/processed/failed/{args.product}.nc"
    if not args.no_upload:
        import ocha_stratus as stratus

        with open(out, "rb") as f:
            stratus.upload_blob_data(f, blob, stage="dev")
        print(f"[ingest] uploaded -> dev projects/{blob}", flush=True)
        if args.out is None:
            shutil.rmtree(workdir, ignore_errors=True)
    if problems:
        raise SystemExit(f"check failed for {args.product} (file kept at {blob}): {problems}")


if __name__ == "__main__":
    main()
