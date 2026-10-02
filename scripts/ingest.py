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
from src.constants import INGEST_START, aux_blob, grid_blob  # noqa: E402
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
    ap.add_argument("--start", default=INGEST_START)
    ap.add_argument("--end", default=(pd.Timestamp.today() - pd.offsets.MonthBegin(1)).strftime("%Y-%m"))
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--no-upload", action="store_true")
    ap.add_argument("--no-check", action="store_true")
    ap.add_argument("--out", type=Path, help="local output path (default: workdir/<product>.nc)")
    args = ap.parse_args()

    stage = os.environ.get("STAGE", "dev")
    if stage != "dev" and not args.no_upload:
        raise SystemExit("this project writes the DEV blob only (STAGE=dev)")
    workdir = default_workdir() / args.product
    workdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    print(f"[ingest] {args.product} {args.start}..{args.end} workdir={workdir}", flush=True)

    if args.product == "static":
        ds, blob = build_static(workdir), aux_blob("static")
    elif args.product == "chirps_v3_diag":
        ds = build_chirps_v3_diag(args.start, args.end, workdir, workers=args.workers)
        blob = aux_blob("chirps_v3_diag")
    else:
        ds = BUILDERS[args.product](args.start, args.end, workdir, workers=args.workers)
        blob = grid_blob(args.product)
        ds.attrs["product"] = args.product
        print(f"[ingest] {args.product}: {ds.sizes['time']} months "
              f"{pd.Timestamp(ds.time.values[0]):%Y-%m}..{pd.Timestamp(ds.time.values[-1]):%Y-%m}", flush=True)
        if args.product != "cru_tmp" and not args.no_check:
            problems = fingerprint(ds, args.product)
            if problems:
                raise SystemExit(f"fingerprint check failed for {args.product}: {problems}")

    ds.attrs["built"] = pd.Timestamp.utcnow().isoformat()
    out = args.out or (workdir / f"{args.product}.nc")
    write_nc(ds, out)
    size = out.stat().st_size / 1e6
    print(f"[ingest] wrote {out} ({size:.0f} MB) in {time.time() - t0:.0f}s", flush=True)

    if not args.no_upload:
        import ocha_stratus as stratus

        with open(out, "rb") as f:
            stratus.upload_blob_data(f, blob, stage="dev")
        print(f"[ingest] uploaded -> dev projects/{blob}", flush=True)
        if args.out is None:
            shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
