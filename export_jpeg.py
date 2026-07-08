import warnings
warnings.filterwarnings("ignore")
import numpy as np
import rasterio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.patches import Rectangle
import os

OUTDIR = "/workspaces/experiment-hub/output"
SRC = os.path.join(OUTDIR, "NDVI_2026-02-18.tif")
DATE = "2026-02-18"

with rasterio.open(SRC) as ds:
    ndvi = ds.read(1)
    transform = ds.transform
    crs = ds.crs

minx, maxx = transform.c, transform.c + transform.a * ds.width
miny, maxy = transform.f + transform.e * ds.height, transform.f
# note: transform.f is top; bottom = f + e*height (e negative)
maxy, miny = transform.f, transform.f + transform.e * ds.height

height, width = ndvi.shape
lat_center = (miny + maxy) / 2.0
m_per_deg_x = 111320.0 * np.cos(np.exp(np.log(np.e) * np.radians(lat_center)))  # =cos
m_per_deg_x = 111320.0 * np.cos(np.radians(lat_center))

# pixel size in meters (near-isotropic at 10 m)
px_x_m = abs(transform.a) * m_per_deg_x
px_y_m = abs(transform.e) * 110540.0
print(f"pixel size: {px_x_m:.2f} m x {px_y_m:.2f} m")

cmap = plt.get_cmap("RdYlGn")
vmin, vmax = -0.1, 0.9
norm = colors.Normalize(vmin=vmin, vmax=vmax)

# ---------------- Publication figure ----------------
fig = plt.figure(figsize=(width/110.0, height/110.0), dpi=110)
ax = fig.add_axes([0.02, 0.10, 0.80, 0.85])
im = ax.imshow(ndvi, cmap=cmap, norm=norm,
               extent=[minx, maxx, miny, maxy], origin="upper",
               interpolation="nearest")
ax.set_title(f"Sentinel-2 L2A NDVI — {DATE}\nBounding box 95.964°E–96.331°E, 21.726°N–22.105°N\n"
             "Source: Microsoft Planetary Computer (European Space Agency Copernicus)",
             fontsize=11, weight="bold")
ax.grid(color="white", alpha=0.35, linestyle=":", linewidth=0.5)
ax.set_xlabel("Longitude (°E)")
ax.set_ylabel("Latitude (°N)")
for s in ax.spines.values():
    s.set_edgecolor("black"); s.set_linewidth(0.8)

# Colorbar
cax = fig.add_axes([0.84, 0.10, 0.025, 0.85])
cb = fig.colorbar(im, cax=cax)
cb.set_label("NDVI", fontsize=11, weight="bold")
cb.ax.tick_params(labelsize=9)

# Scale bar (5 km) in data (lon) coordinates
km = 5.0
bar_deg = km * 1000.0 / m_per_deg_x
x0 = minx + (maxx - minx) * 0.04
y0 = miny + (maxy - miny) * 0.04
ax.plot([x0, x0 + bar_deg], [y0, y0], color="black", lw=2.2)
ax.plot([x0, x0], [y0 - bar_deg*0.15, y0 + bar_deg*0.15], color="black", lw=2.2)
ax.plot([x0 + bar_deg, x0 + bar_deg], [y0 - bar_deg*0.15, y0 + bar_deg*0.15], color="black", lw=2.2)
ax.text(x0 + bar_deg/2, y0 + bar_deg*0.35, f"{km} km", color="black",
        fontsize=10, ha="center", weight="bold",
        bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=1))

# North arrow (image is north-up)
nx = minx + (maxx - minx) * 0.93
ny = maxy - (maxy - miny) * 0.06
ax.annotate("", xy=(nx, ny - bar_deg*0.9), xytext=(nx, ny),
            arrowprops=dict(arrowstyle="-|>", color="black", lw=2.2))
ax.text(nx, ny + bar_deg*0.25, "N", color="black", fontsize=12,
        ha="center", weight="bold", bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=1))

pub = os.path.join(OUTDIR, "NDVI_2026-02-18_publication.jpg")
fig.savefig(pub, dpi=110, bbox_inches="tight")
print("Saved:", pub)

# ---------------- Plain colormap raster (native resolution) ----------------
rgb = (cmap(norm(ndvi))[:, :, :3] * 255).astype(np.uint8)
# mask invalid (none, but safe)
mask = np.isnan(ndvi)
rgb[mask] = (220, 220, 220)
from PIL import Image
img = Image.fromarray(rgb, "RGB")
col = os.path.join(OUTDIR, "NDVI_2026-02-18_colormap.jpg")
img.save(col, "JPEG", quality=95, subsampling=0)
print("Saved:", col, img.size)
