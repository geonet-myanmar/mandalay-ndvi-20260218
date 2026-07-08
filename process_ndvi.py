import warnings
warnings.filterwarnings("ignore")
import os, numpy as np
import rasterio
import rasterio.warp
from rasterio.warp import reproject, Resampling
import planetary_computer as pc
from pystac_client import Client

BBOX = [95.964203, 21.725697, 96.330872, 22.104726]   # minx,miny,maxx,maxy (lon,lat)
DATE = "2026-02-18"
OUTDIR = "/workspaces/experiment-hub/output"
os.makedirs(OUTDIR, exist_ok=True)

# ---------- 1. Query items for the chosen date ----------
catalog = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1")
search = catalog.search(collections=["sentinel-2-l2a"], bbox=BBOX,
                        datetime=f"{DATE}T00:00:00Z/{DATE}T23:59:59Z", max_items=50)
items = list(search.get_items())
print(f"Found {len(items)} tiles for {DATE}")

# ---------- 2. Build a common target grid in EPSG:4326 at 10 m ----------
minx, miny, maxx, maxy = BBOX
lat_center = (miny + maxy) / 2.0
m_per_deg_x = 111320.0 * np.cos(np.radians(lat_center))   # meters per degree longitude
m_per_deg_y = 110540.0                                     # approx meters per degree latitude
target_res_m = 10.0
res_x = target_res_m / m_per_deg_x   # degrees per pixel (x)
res_y = target_res_m / m_per_deg_y   # degrees per pixel (y)

# align to pixel grid
minx_a = np.floor(minx / res_x) * res_x
maxx_a = np.ceil(maxx / res_x) * res_x
miny_a = np.floor(miny / res_y) * res_y
maxy_a = np.ceil(maxy / res_y) * res_y

width = int(round((maxx_a - minx_a) / res_x))
height = int(round((maxy_a - miny_a) / res_y))
transform = rasterio.transform.from_bounds(minx_a, miny_a, maxx_a, maxy_a, width, height)
crs = "EPSG:4326"
print(f"Target grid: {width} x {height} px, res~10m, EPSG:4326")

# ---------- 3. Download & mosaic B04 and B08 ----------
def mosaic_band(asset_key):
    acc = np.full((height, width), np.nan, dtype=np.float32)
    for it in items:
        href = pc.sign(it.assets[asset_key].href)
        with rasterio.open(href) as src:
            tmp = np.full((height, width), np.nan, dtype=np.float32)
            reproject(
                source=rasterio.band(src, 1),
                destination=tmp,
                src_crs=src.crs, src_transform=src.transform,
                dst_crs=crs, dst_transform=transform,
                dst_nodata=np.nan,
                resampling=Resampling.bilinear,
            )
        # S2 L2A reflectance = DN / 10000
        tmp = tmp / 10000.0
        acc = np.where(np.isnan(acc), tmp, acc)
    return acc

print("Mosaicking B04 (red)...")
b04 = mosaic_band("B04")
print("Mosaicking B08 (NIR)...")
b08 = mosaic_band("B08")

valid = (~np.isnan(b04)) & (~np.isnan(b08))
print(f"Valid pixels: {valid.sum()} / {valid.size} ({100*valid.sum()/valid.size:.1f}%)")

# ---------- 4. Compute NDVI ----------
ndvi = np.full_like(b04, np.nan, dtype=np.float32)
denom = (b08 + b04)
safe = valid & (denom > 1e-6)
ndvi[safe] = (b08[safe] - b04[safe]) / denom[safe]
ndvi = np.clip(ndvi, -1.0, 1.0)

# Save the NDVI raster as GeoTIFF (full quality, georeferenced)
gtiff = os.path.join(OUTDIR, "NDVI_2026-02-18.tif")
profile = {"driver":"GTiff","dtype":"float32","count":1,"width":width,"height":height,
           "crs":crs,"transform":transform,"nodata":np.nan,"compress":"deflate"}
with rasterio.open(gtiff, "w", **profile) as dst:
    dst.write(ndvi, 1)
    dst.set_band_description(1, "NDVI")
print(f"Saved GeoTIFF: {gtiff}")

# Save NDVI stats
print("NDVI stats (valid): min=%.3f mean=%.3f max=%.3f" % (
    np.nanmin(ndvi), np.nanmean(ndvi), np.nanmax(ndvi)))

np.save(os.path.join(OUTDIR, "ndvi.npy"), ndvi)
np.save(os.path.join(OUTDIR, "transform.npy"), np.array(transform.column_vectors()).reshape(9) if False else np.array(transform))
