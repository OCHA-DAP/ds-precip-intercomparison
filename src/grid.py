"""Regridding onto the common 0.5 deg grid.

Every source grid used here either nests inside the 0.5 deg cells (CHIRPS 0.05,
IMERG 0.1, GPCC/CRU/CPC/PREC-L/UDel 0.5, GPCC Monitoring 1.0 - edges on multiples
of the resolution) or is ERA5's point-centred 0.25 grid. Nesting grids are
block-averaged exactly; ERA5 uses its own area-weighted kernel (`era5_to_05`).
"""

import numpy as np
import xarray as xr

from src.constants import LAT, LON, NLAT, NLON, RES


def normalise(da: xr.DataArray) -> xr.DataArray:
    """Rename to lat/lon, lon to [-180, 180), both ascending."""
    ren = {}
    for name in da.dims:
        low = name.lower()
        if low in ("latitude", "y") and "lat" not in da.dims:
            ren[name] = "lat"
        if low in ("longitude", "x") and "lon" not in da.dims:
            ren[name] = "lon"
    da = da.rename(ren)
    if float(da.lon.max()) > 180:
        da = da.assign_coords(lon=((da.lon + 180) % 360) - 180)
    return da.sortby("lat").sortby("lon")


def assert_on_grid(da: xr.DataArray, name: str) -> None:
    """Fail loudly unless `da` is exactly on the common 0.5 deg grid."""
    lat = np.round(da.lat.values.astype(float), 3)
    lon = np.round(da.lon.values.astype(float), 3)
    if lat.size != NLAT or lon.size != NLON or not (
        np.allclose(lat, LAT) and np.allclose(lon, LON)
    ):
        raise ValueError(
            f"{name}: not on the 0.5 deg grid "
            f"(lat {lat[:2]}..{lat[-1:]} n={lat.size}, lon {lon[:2]}..{lon[-1:]} n={lon.size})"
        )


def block_mean(
    arr: np.ndarray, factor: int, min_valid: float = 0.0
) -> tuple[np.ndarray, np.ndarray]:
    """NaN-aware block mean of a (ny, nx) array, ny and nx multiples of `factor`.

    Returns (mean, valid_fraction). Cells whose valid fraction is <= `min_valid`
    are NaN.
    """
    ny, nx = arr.shape
    if ny % factor or nx % factor:
        raise ValueError(f"shape {arr.shape} not divisible by {factor}")
    b = arr.reshape(ny // factor, factor, nx // factor, factor)
    valid = np.isfinite(b)
    n = valid.sum(axis=(1, 3))
    s = np.where(valid, b, 0).sum(axis=(1, 3), dtype=np.float64)
    frac = n / factor**2
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where((n > 0) & (frac > min_valid), s / np.maximum(n, 1), np.nan)
    return mean.astype(np.float32), frac.astype(np.float32)


def place_fine(arr: np.ndarray, top: float, left: float, res: float) -> np.ndarray:
    """Embed a north-up fine-grid array into a global north-up array at `res`.

    `top`/`left` are the outer edges of the source raster. The result covers
    90N..90S, 180W..180E with NaN outside the source footprint.
    """
    ny_g, nx_g = int(round(180 / res)), int(round(360 / res))
    r0 = int(round((90 - top) / res))
    c0 = int(round((left + 180) / res))
    if c0 != 0 or arr.shape[1] != nx_g:
        raise ValueError(f"expected full-longitude raster, got left={left} nx={arr.shape[1]}")
    out = np.full((ny_g, nx_g), np.nan, dtype=np.float32)
    out[r0 : r0 + arr.shape[0], :] = arr
    return out


def fine_to_05(
    arr: np.ndarray, top: float, left: float, res: float, min_valid: float = 0.0
) -> tuple[np.ndarray, np.ndarray]:
    """North-up fine raster (nodata already NaN) -> (mean, valid_frac) on the
    common grid, both south-up (ascending lat) to match LAT."""
    factor = int(round(RES / res))
    if not np.isclose(factor * res, RES):
        raise ValueError(f"resolution {res} does not nest in {RES}")
    full = place_fine(arr, top, left, res)
    mean, frac = block_mean(full, factor, min_valid)
    return mean[::-1], frac[::-1]


def era5_to_05(arr: np.ndarray, lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """ERA5 0.25 deg point-centred field -> 0.5 deg cell means.

    ERA5 grid points sit on multiples of 0.25 and each represents the cell
    +-0.125 around it. A 0.5 deg target cell [c-0.25, c+0.25] (c on x.25/x.75)
    therefore contains the point at c fully and the points at c+-0.25 by half:
    1-D weights (0.25, 0.5, 0.25), separable in lat and lon. cos(lat) variation
    inside one cell is ignored (< 0.2% effect away from the poles).
    """
    da = normalise(xr.DataArray(arr, coords={"lat": lat, "lon": lon}, dims=("lat", "lon")))
    a = da.values.astype(np.float64)  # lat ascending -90..90 (721), lon -180..179.75 (1440)
    lat_s = np.round(da.lat.values, 3)
    lon_s = np.round(da.lon.values, 3)
    if lat_s.size != 721 or lon_s.size != 1440 or lat_s[0] != -90 or lon_s[0] != -180:
        raise ValueError(f"unexpected ERA5 grid lat {lat_s[:2]} lon {lon_s[:2]}")
    w = np.array([0.25, 0.5, 0.25])
    # Longitude: target centre c = -179.75 + 0.5 j sits on source index 1 + 2 j;
    # neighbours 2j, 2j+1, 2j+2 (index 1440 wraps to 0).
    a_wrap = np.concatenate([a, a[:, :1]], axis=1)  # 1441 columns
    lon_out = (
        w[0] * a_wrap[:, 0:-1:2] + w[1] * a_wrap[:, 1::2] + w[2] * a_wrap[:, 2::2]
    )  # (721, 720)
    # Latitude: target centre -89.75 + 0.5 i sits on source index 1 + 2 i.
    out = w[0] * lon_out[0:-1:2] + w[1] * lon_out[1::2] + w[2] * lon_out[2::2]  # (360, 720)
    return out.astype(np.float32)


def replicate(arr: np.ndarray, factor: int) -> np.ndarray:
    """Coarse grid nesting the 0.5 deg grid (e.g. 1 deg) -> repeat into 0.5 deg."""
    return np.repeat(np.repeat(arr, factor, axis=0), factor, axis=1)


def empty_cube(ntime: int) -> np.ndarray:
    return np.full((ntime, NLAT, NLON), np.nan, dtype=np.float32)
