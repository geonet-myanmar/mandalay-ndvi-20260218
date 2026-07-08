"""Download cloud-free Sentinel-2 Level-2A imagery and compute NDVI.

This module queries the Microsoft Planetary Computer STAC API for Sentinel-2
L2A scenes that intersect a user-supplied bounding box, selects the best
available acquisition (lowest cloud cover while still fully covering the box),
mosaics the required spectral bands (red B04 and near-infrared B08), and writes
a georeferenced NDVI raster.

Why mosaicking is needed
------------------------
A single Sentinel-2 granule is delivered in one UTM zone. Bounding boxes that
straddle a UTM zone boundary (e.g. the 96°E boundary between zones 46N and 47N)
are therefore covered by *several* granules. This tool automatically gathers all
granules for a given acquisition date and mosaics them into a single raster so
the whole bounding box is covered.

Usage
-----
    python src/sentinel2_ndvi.py \
        --bbox 95.964203 21.725697 96.330872 22.104726 \
        --date 2026-02-18 \
        --outdir output

Omit ``--date`` to automatically pick the least-cloudy acquisition that still
fully covers the bounding box.
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import planetary_computer as pc
import rasterio
from pystac_client import Client
from rasterio.transform import from_bounds
from rasterio.warp import Resampling, reproject
import shapely.geometry as sg
from shapely.ops import unary_union
from warnings import filterwarnings

filterwarnings("ignore")

# Microsoft Planetary Computer STAC endpoint and collection name.
STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
COLLECTION = "sentinel-2-l2a"

# S2 L2A surface reflectance is encoded as uint16 with a 1/10000 scale factor.
REFLECTANCE_SCALE = 10000.0


def connect() -> Client:
    """Open a connection to the Planetary Computer STAC API."""
    return Client.open(STAC_URL)


def fetch_items(catalog: Client, bbox, start: str, end: str):
    """Return all S2 L2A items intersecting ``bbox`` within [start, end]."""
    search = catalog.search(
        collections=[COLLECTION],
        bbox=list(bbox),
        datetime=f"{start}T00:00:00Z/{end}T23:59:59Z",
        max_items=5000,
    )
    return list(search.get_items())


def select_best_date(items, bbox, prefer_recent: bool = True):
    """Group items by acquisition date and return those fully covering ``bbox``.

    Returns a list of ``(date, items, mean_cloud, max_cloud)`` tuples sorted by
    ascending mean cloud cover, with ties broken by date (most recent first when
    ``prefer_recent`` is True, otherwise earliest first).
    """
    box = sg.box(*bbox)
    by_date: dict[str, list] = defaultdict(list)
    for it in items:
        by_date[it.properties["datetime"][:10]].append(it)

    candidates = []
    for date, day_items in by_date.items():
        union = unary_union([sg.shape(it.geometry) for it in day_items])
        # small buffer tolerates sub-pixel edge gaps between granules
        if union.buffer(1e-4).contains(box):
            clouds = [it.properties["eo:cloud_cover"] for it in day_items]
            candidates.append((date, day_items, float(np.mean(clouds)), float(np.max(clouds))))

    if not candidates:
        raise RuntimeError("No acquisition date fully covers the bounding box.")

    candidates.sort(key=lambda c: (c[2], -int(c[0].replace("-", "")) if prefer_recent else int(c[0].replace("-", ""))))
    return candidates


def build_target_grid(bbox, res_m: float = 10.0):
    """Build a common EPSG:4326 target grid (~``res_m`` metres per pixel).

    EPSG:4326 is used because it spans UTM zone boundaries, allowing granules
    from different zones to be mosaicked into one raster.
    """
    minx, miny, maxx, maxy = bbox
    lat_center = (miny + maxy) / 2.0
    m_per_deg_x = 111320.0 * np.cos(np.radians(lat_center))
    m_per_deg_y = 110540.0

    res_x = res_m / m_per_deg_x
    res_y = res_m / m_per_deg_y

    minx = np.floor(minx / res_x) * res_x
    maxx = np.ceil(maxx / res_x) * res_x
    miny = np.floor(miny / res_y) * res_y
    maxy = np.ceil(maxy / res_y) * res_y

    width = int(round((maxx - minx) / res_x))
    height = int(round((maxy - miny) / res_y))
    transform = from_bounds(minx, miny, maxx, maxy, width, height)
    return transform, width, height, "EPSG:4326"


def mosaic_band(items, asset_key: str, transform, width: int, height: int, crs: str) -> np.ndarray:
    """Reproject and mosaic one band across all granules into the target grid."""
    acc = np.full((height, width), np.nan, dtype=np.float32)
    for it in items:
        href = pc.sign(it.assets[asset_key].href)
        with rasterio.open(href) as src:
            tmp = np.full((height, width), np.nan, dtype=np.float32)
            reproject(
                source=rasterio.band(src, 1),
                destination=tmp,
                src_crs=src.crs,
                src_transform=src.transform,
                dst_crs=crs,
                dst_transform=transform,
                dst_nodata=np.nan,
                resampling=Resampling.bilinear,
            )
        tmp = tmp / REFLECTANCE_SCALE
        acc = np.where(np.isnan(acc), tmp, acc)
    return acc


def compute_ndvi(red: np.ndarray, nir: np.ndarray) -> np.ndarray:
    """NDVI = (NIR - Red) / (NIR + Red), clipped to [-1, 1]."""
    ndvi = np.full_like(red, np.nan, dtype=np.float32)
    denom = nir + red
    valid = (~np.isnan(red)) & (~np.isnan(nir)) & (denom > 1e-6)
    ndvi[valid] = (nir[valid] - red[valid]) / denom[valid]
    return np.clip(ndvi, -1.0, 1.0)


def write_geotiff(path: str, ndvi: np.ndarray, transform, crs: str) -> None:
    """Write the NDVI array to a compressed, georeferenced GeoTIFF."""
    profile = {
        "driver": "GTiff",
        "dtype": "float32",
        "count": 1,
        "width": ndvi.shape[1],
        "height": ndvi.shape[0],
        "crs": crs,
        "transform": transform,
        "nodata": np.nan,
        "compress": "deflate",
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(ndvi, 1)
        dst.set_band_description(1, "NDVI")


def run(bbox, outdir: str, date: str | None = None, res_m: float = 10.0,
        start: str = "2019-01-01", end: str = "2026-12-31"):
    os.makedirs(outdir, exist_ok=True)
    catalog = connect()
    items = fetch_items(catalog, bbox, start, end)

    if date:
        chosen = [it for it in items if it.properties["datetime"][:10] == date]
        if not chosen:
            raise RuntimeError(f"No scenes found for requested date {date}.")
        day_items = chosen
        used_date = date
        mean_cloud = float(np.mean([it.properties["eo:cloud_cover"] for it in chosen]))
        max_cloud = float(np.max([it.properties["eo:cloud_cover"] for it in chosen]))
    else:
        candidates = select_best_date(items, bbox)
        used_date, day_items, mean_cloud, max_cloud = candidates[0]

    print(f"Selected acquisition date: {used_date}")
    print(f"Granules used: {len(day_items)}  |  mean cloud {mean_cloud:.4f}%  max cloud {max_cloud:.4f}%")

    transform, width, height, crs = build_target_grid(bbox, res_m)
    print(f"Target grid: {width} x {height} px at ~{res_m:.0f} m (EPSG:4326)")

    print("Mosaicking B04 (red)...")
    b04 = mosaic_band(day_items, "B04", transform, width, height, crs)
    print("Mosaicking B08 (NIR)...")
    b08 = mosaic_band(day_items, "B08", transform, width, height, crs)

    valid = (~np.isnan(b04)) & (~np.isnan(b08))
    print(f"Valid coverage: {100 * valid.sum() / valid.size:.1f}% of bounding box")

    ndvi = compute_ndvi(b04, b08)
    print(f"NDVI: min={np.nanmin(ndvi):.3f} mean={np.nanmean(ndvi):.3f} max={np.nanmax(ndvi):.3f}")

    tif_path = os.path.join(outdir, f"NDVI_{used_date}.tif")
    write_geotiff(tif_path, ndvi, transform, crs)
    print(f"Saved GeoTIFF: {tif_path}")
    return tif_path


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Compute NDVI from Planetary Computer Sentinel-2 L2A.")
    p.add_argument("--bbox", nargs=4, type=float, required=True,
                   metavar=("MINX", "MINY", "MAXX", "MAXY"),
                   help="Bounding box in lon/lat: min_lon min_lat max_lon max_lat")
    p.add_argument("--date", type=str, default=None,
                   help="Acquisition date YYYY-MM-DD (default: auto-select least cloudy)")
    p.add_argument("--outdir", type=str, default="output")
    p.add_argument("--res", type=float, default=10.0, help="Output pixel size in metres")
    p.add_argument("--start", type=str, default="2019-01-01", help="Search start date")
    p.add_argument("--end", type=str, default="2026-12-31", help="Search end date")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    run(args.bbox, args.outdir, args.date, args.res, args.start, args.end)


if __name__ == "__main__":
    main()
