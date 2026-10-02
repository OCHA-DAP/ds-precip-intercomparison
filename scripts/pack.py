"""Pull the 0.5 deg cubes from the dev blob and pack land cells into a local cache.

    python scripts/pack.py            # all products present on blob
    python scripts/pack.py era5 cru   # selected

Writes cache/packed/<product>.npz with arrays of shape (time, cell) for every
(time, lat, lon) variable, plus cache/packed/cells.npz describing the cells
(land fraction >= 0.5, south of 60S excluded). The full cubes are deleted after
packing: the laptop has little disk, the analysis only needs land cells.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import ocha_stratus as stratus  # noqa: E402

from src.constants import GRID_DIR, LAT, LON, aux_blob  # noqa: E402

CACHE = Path(__file__).resolve().parents[1] / "cache"
PACKED = CACHE / "packed"


def fetch(blob: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        return dest
    cc = stratus.get_container_client(stage="dev")
    tmp = dest.with_suffix(".part")
    with open(tmp, "wb") as f:
        cc.get_blob_client(blob).download_blob(max_concurrency=4).readinto(f)
    tmp.rename(dest)
    return dest


def build_cells() -> dict:
    static = xr.open_dataset(fetch(aux_blob("static"), CACHE / "static.nc")).load()
    lat2, lon2 = np.meshgrid(LAT, LON, indexing="ij")
    mask = (static.landfrac.values >= 0.5) & (lat2 >= -60)
    iy, ix = np.nonzero(mask)
    cells = {
        "iy": iy.astype(np.int16),
        "ix": ix.astype(np.int16),
        "lat": LAT[iy].astype(np.float32),
        "lon": LON[ix].astype(np.float32),
        "area": np.cos(np.deg2rad(LAT[iy])).astype(np.float32),
        "landfrac": static.landfrac.values[iy, ix].astype(np.float32),
        "koppen": static.koppen.values[iy, ix].astype(np.int16),
        "lw_factor": static.lw_factor.values[:, iy, ix].astype(np.float32),  # (12, cell)
    }
    PACKED.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(PACKED / "cells.npz", **cells)
    print(f"[pack] {iy.size} land cells")
    return cells


def pack(product: str, cells: dict, aux: bool = False) -> None:
    blob = aux_blob(product) if aux else f"{GRID_DIR}/{product}.nc"
    path = fetch(blob, CACHE / "nc" / f"{product}.nc")
    ds = xr.open_dataset(path)
    out = {"time": pd.DatetimeIndex(ds.time.values).values.astype("datetime64[M]").astype(str)}
    for v, da in ds.data_vars.items():
        if da.dims == ("time", "lat", "lon"):
            out[v] = da.values[:, cells["iy"], cells["ix"]].astype(np.float32)
        elif da.dims == ("lat", "lon"):
            out[v] = da.values[cells["iy"], cells["ix"]].astype(np.float32)
        elif da.dims == ("time",):
            out[v] = da.values
    out["attrs"] = np.array(str(ds.attrs))
    ds.close()
    np.savez_compressed(PACKED / f"{product}.npz", **out)
    path.unlink()
    print(f"[pack] {product}: {len(out['time'])} months {out['time'][0]}..{out['time'][-1]}")


def main(argv: list[str]) -> None:
    cells = build_cells() if not (PACKED / "cells.npz").exists() else dict(np.load(PACKED / "cells.npz"))
    on_blob = {
        Path(b).stem
        for b in stratus.list_container_blobs(name_starts_with=GRID_DIR + "/", stage="dev")
        if b.endswith(".nc")
    }
    wanted = argv or sorted(on_blob)
    for p in wanted:
        if p not in on_blob:
            print(f"[pack] {p}: not on blob yet, skipped")
            continue
        if (PACKED / f"{p}.npz").exists():
            continue
        pack(p, cells)
    if not (PACKED / "chirps_v3_diag.npz").exists():
        try:
            pack("chirps_v3_diag", cells, aux=True)
        except Exception as exc:  # noqa: BLE001
            print(f"[pack] chirps_v3_diag: {exc}")


if __name__ == "__main__":
    main(sys.argv[1:])
