"""Export a publication-ready NDVI JPEG from a GeoTIFF produced by
``src/sentinel2_ndvi.py``.

Two products are written:

* ``<date>_publication.jpg`` — annotated figure with colour bar, scale bar,
  north arrow, title and lat/lon grid (best for papers / reports).
* ``<date>_colormap.jpg``  — plain native-resolution colormapped raster with no
  annotations (best as a raw data product / quick look).

Usage
-----
    python src/export_ndvi_jpeg.py --tif output/NDVI_2026-02-18.tif
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import rasterio
from matplotlib import cm, colors
from matplotlib import pyplot as plt
from PIL import Image
from warnings import filterwarnings

filterwarnings("ignore")

# Display range for NDVI and the vegetation-oriented colour map.
CMAP = "RdYlGn"          # red (low) -> yellow -> green (high vegetation)
VMIN, VMAX = -0.1, 0.9


def load(tif_path: str):
    with rasterio.open(tif_path) as ds:
        ndvi = ds.read(1)
        transform = ds.transform
        height, width = ndvi.shape
    # GeoTIFF transform: (c, a, b, f, d, e). Top-left is (c, f).
    minx = transform.c
    maxx = transform.c + transform.a * width
    maxy = transform.f
    miny = transform.f + transform.e * height
    return ndvi, transform, (minx, miny, maxx, maxy), width, height


def metres_per_degree(lat_center: float):
    return 111320.0 * np.cos(np.radians(lat_center)), 110540.0


def export(tif_path: str, outdir: str, dpi: int = 110):
    os.makedirs(outdir, exist_ok=True)
    ndvi, transform, (minx, miny, maxx, maxy), width, height = load(tif_path)
    lat_center = (miny + maxy) / 2.0
    mdx, mdy = metres_per_degree(lat_center)

    norm = colors.Normalize(vmin=VMIN, vmax=VMAX)
    cmap = plt.get_cmap(CMAP)
    date = os.path.basename(tif_path).replace("NDVI_", "").replace(".tif", "")

    # ---------------- Annotated publication figure ----------------
    fig = plt.figure(figsize=(width / dpi, height / dpi), dpi=dpi)
    ax = fig.add_axes([0.02, 0.10, 0.80, 0.85])
    im = ax.imshow(
        ndvi, cmap=cmap, norm=norm, extent=[minx, maxx, miny, maxy],
        origin="upper", interpolation="nearest",
    )
    ax.set_title(
        f"Sentinel-2 L2A NDVI — {date}\n"
        f"Bounding box {minx:.3f}°E–{maxx:.3f}°E, {miny:.3f}°N–{maxy:.3f}°N\n"
        "Source: Microsoft Planetary Computer (ESA Copernicus)",
        fontsize=11, weight="bold",
    )
    ax.grid(color="white", alpha=0.35, linestyle=":", linewidth=0.5)
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    for s in ax.spines.values():
        s.set_edgecolor("black")
        s.set_linewidth(0.8)

    cax = fig.add_axes([0.84, 0.10, 0.025, 0.85])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("NDVI", fontsize=11, weight="bold")
    cb.ax.tick_params(labelsize=9)

    # 5 km scale bar drawn in geographic (longitude) coordinates.
    km = 5.0
    bar_deg = km * 1000.0 / mdx
    x0 = minx + (maxx - minx) * 0.04
    y0 = miny + (maxy - miny) * 0.04
    ax.plot([x0, x0 + bar_deg], [y0, y0], color="black", lw=2.2)
    ax.plot([x0, x0], [y0 - bar_deg * 0.15, y0 + bar_deg * 0.15], color="black", lw=2.2)
    ax.plot([x0 + bar_deg, x0 + bar_deg], [y0 - bar_deg * 0.15, y0 + bar_deg * 0.15], color="black", lw=2.2)
    ax.text(x0 + bar_deg / 2, y0 + bar_deg * 0.35, f"{km} km", color="black",
            fontsize=10, ha="center", weight="bold",
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=1))

    # North arrow (image is north-up).
    nx = minx + (maxx - minx) * 0.93
    ny = maxy - (maxy - miny) * 0.06
    ax.annotate("", xy=(nx, ny - bar_deg * 0.9), xytext=(nx, ny),
                arrowprops=dict(arrowstyle="-|>", color="black", lw=2.2))
    ax.text(nx, ny + bar_deg * 0.25, "N", color="black", fontsize=12, ha="center",
            weight="bold", bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=1))

    pub_path = os.path.join(outdir, f"NDVI_{date}_publication.jpg")
    fig.savefig(pub_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {pub_path}")

    # ---------------- Plain colormap raster (native resolution) ----------------
    rgba = cmap(norm(ndvi))
    rgb = (rgba[:, :, :3] * 255).astype(np.uint8)
    rgb[np.isnan(ndvi)] = (220, 220, 220)
    col_path = os.path.join(outdir, f"NDVI_{date}_colormap.jpg")
    Image.fromarray(rgb, "RGB").save(col_path, "JPEG", quality=95, subsampling=0)
    print(f"Saved: {col_path}  ({rgb.shape[1]} x {rgb.shape[0]} px)")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Export a publication-ready NDVI JPEG.")
    p.add_argument("--tif", type=str, required=True, help="Input NDVI GeoTIFF")
    p.add_argument("--outdir", type=str, default="output")
    p.add_argument("--dpi", type=int, default=110, help="Figure resolution")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    export(args.tif, args.outdir, args.dpi)


if __name__ == "__main__":
    main()
