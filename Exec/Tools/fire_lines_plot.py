#!/usr/bin/env python3
"""Map frames, movie, line-temperature and tower-temperature plots of one demo run.

    fire_lines_plot.py <case> <run_dir> <basemap.tif> <out_dir> --domain <case_domain.json> [--title T]
                       [--every N] [--crop x0,x1,y0,y1] [--trange lo,hi]

<case> is a short label used in file names. --crop is a zoom window in km; --trange the conductor colour range
in deg F (fit it to the run). <case_domain.json> is the Real_Terrain domain file (UTM zone, SW corner, size).

Reads plt_fire_*, conductors/conductor_temperature.dat, conductors/tower_*_member_temperature.dat and
lines_meta.json in run_dir. Writes into out_dir:
  <case>_map_final.png   burned area by arrival time over imagery, the lines coloured by their peak temperature
  <case>_map.mp4         the same through the run, the lines coloured by their current temperature
  <case>_line_T.png      hottest span of each line and hottest tower member against time
  <case>_spans.csv       per span: case, line, span, length_m, peak_T_C, peak_fire_flux_W_m2, closest_burn_m,
                         fire_dx_m (the fire cell size)
  <case>_towers.csv      per tower: case, line, tower, x, y, peak_steel_T_C, closest_burn_m
"""
import argparse, glob, json, os, re, shutil, subprocess, sys
import numpy as np
import yt
yt.set_log_level(50)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling
from pyproj import Transformer

ap = argparse.ArgumentParser()
ap.add_argument("case"); ap.add_argument("run"); ap.add_argument("basemap"); ap.add_argument("out")
ap.add_argument("--domain", required=True, help="Real_Terrain <case>_domain.json")
ap.add_argument("--title", default="", help="figure title (default: the case label)")
ap.add_argument("--every", type=int, default=4, help="use every N-th plt_fire for the movie (N >= 1)")
ap.add_argument("--crop", default="")
ap.add_argument("--trange", default="32,300", help="conductor colour range lo,hi in deg F (fit it to the run)")
a = ap.parse_args()
if a.every < 1:
    ap.error("--every must be at least 1")
try:
    TLO_F, TMAX_F = (float(v) for v in a.trange.split(","))
except ValueError:
    ap.error("--trange must be two numbers lo,hi in deg F")
if not TLO_F < TMAX_F:
    ap.error("--trange must have lo < hi")
try:
    crop = [float(c) for c in a.crop.split(",")] if a.crop else None
except ValueError:
    ap.error("--crop must be four numbers x0,x1,y0,y1 in km")
if crop and not (len(crop) == 4 and crop[0] < crop[1] and crop[2] < crop[3]):
    ap.error("--crop must be four numbers x0,x1,y0,y1 in km with x0 < x1 and y0 < y1")
os.makedirs(a.out, exist_ok=True)
TITLE = a.title or a.case
DOM = json.load(open(a.domain))
L = float(DOM["domain_m"])
meta = json.load(open(os.path.join(a.run, "lines_meta.json")))
T_IGN = float([l for l in open(os.path.join(a.run, 'ignitions_spinup.csv')) if not l.startswith('#')][0].split()[0])
IGN = DOM.get("ignition") or DOM["ignitions"][0]

# ---------------------------------------------------------------- basemap on the model square
zone = DOM["utm_zone"]; utm = f"EPSG:{(32600 if zone.upper().endswith('N') else 32700) + int(zone[:-1])}"
x0, y0 = Transformer.from_crs("EPSG:4326", utm, always_xy=True).transform(*DOM["sw_corner_lonlat"])
res = 10.0; nb = int(round(L / res)); rgb = np.zeros((3, nb, nb), np.uint8)
with rasterio.open(a.basemap) as src:
    credit = src.tags().get("attribution", "Imagery: U.S. Geological Survey")
    for b in range(3):
        reproject(source=src.read(b + 1), destination=rgb[b], src_transform=src.transform, src_crs=src.crs,
                  dst_transform=from_origin(x0, y0 + L, res, res), dst_crs=utm, resampling=Resampling.bilinear)
rgb = rgb.transpose(1, 2, 0)
credit += "; lines (c) OpenStreetMap contributors; fuel: LANDFIRE; wind: NOAA HRRR"

# ---------------------------------------------------------------- conductor and tower temperatures
ct = np.genfromtxt(os.path.join(a.run, "conductors", "conductor_temperature.dat"), dtype=None, encoding=None,
                   names=["t", "line", "span", "Tm", "Tx", "node", "Ta", "w", "wd", "E", "Ex", "strain", "Tr"])
towers = {}
for fn in glob.glob(os.path.join(a.run, "conductors", "tower_*_member_temperature.dat")):
    m = re.search(r"tower_(\w+?)_t(\d+)_member_temperature", fn)
    d = np.loadtxt(fn, comments="#")
    towers[(m.group(1), int(m.group(2)))] = (d[:, 0], d[:, 1:].max(axis=1))


def span_T(t, field="Tx"):
    """{(line, span): value} at the diagnostics time closest to t (and at or before)."""
    ts = np.unique(ct["t"]); tt = ts[max(0, np.searchsorted(ts, t, side="right") - 1)]
    sel = ct[ct["t"] == tt]
    return {(r["line"], int(r["span"])): r[field] for r in sel}


peak = {}
for r in ct:
    k = (r["line"], int(r["span"]))
    p = peak.get(k, (-1e9, 0.0))
    peak[k] = (max(p[0], r["Tx"]), max(p[1], r["Ex"]))

# ---------------------------------------------------------------- fire plotfiles
files = sorted(glob.glob(os.path.join(a.run, "plt_fire_?????")))
FIRE_MIN = max(5.0, round((float(yt.load(files[-1]).current_time) - T_IGN) / 60.0))
def miles_axes(ax):
    """Label a km-coordinate map in miles (US audience); the data stay in km."""
    MI = 1.609344
    for get, setl, sett, lab in ((ax.get_xlim, ax.set_xlim, ax.set_xticks, ax.set_xlabel),
                                 (ax.get_ylim, ax.set_ylim, ax.set_yticks, ax.set_ylabel)):
        lo, hi = get(); step = 1 if (hi - lo) / MI <= 8 else 2
        t = np.arange(np.ceil(lo / MI / step) * step, hi / MI + 1e-9, step)
        sett(t * MI, [f"{v:g}" for v in t]); setl(lo, hi)
    ax.set_xlabel("miles east"); ax.set_ylabel("miles north")


def F(c):
    """deg C -> deg F (plots are in F for a US audience; the run and the csv files stay in C)."""
    return np.asarray(c) * 9.0 / 5.0 + 32.0
WIRE_FLAG_F = F(100.0)        # the chosen 100 C wire flag (fire_lines_fragility.py), not a rating
norm = plt.Normalize(TLO_F, TMAX_F); cmapL = plt.get_cmap("plasma")
ext = [0, L / 1000, 0, L / 1000]


def segs():
    out = []
    for m in meta:
        p = np.array(m["pts"]) / 1000.0
        for k in range(len(p) - 1):
            out.append(((m["name"], k), [p[k], p[k + 1]]))
    return out


SEG = segs()


def draw(ds, Tmap, fn, label):
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    at = np.array(g["fire_arrival_time"])[:, :, 0]; phi = np.array(g["fire_phi"])[:, :, 0]
    t = float(ds.current_time); nx = at.shape[0]; dx = L / nx
    xc = (np.arange(nx) + 0.5) * dx / 1000.0
    fig, ax = plt.subplots(figsize=(8.0, 6.6))
    ax.imshow(rgb, origin="upper", extent=ext, interpolation="bilinear")
    burned = np.ma.masked_where(at < 0, (at - T_IGN) / 60.0)
    im = ax.imshow(burned.T, cmap="YlOrRd_r", origin="lower", extent=ext, vmin=0, vmax=FIRE_MIN, alpha=0.65)
    if (phi < 0).any():
        ax.contour(xc, xc, phi.T, levels=[0.0], colors="red", linewidths=1.4)
    lc = LineCollection([s for _, s in SEG], colors=[cmapL(norm(F(Tmap.get(k, 0)))) for k, _ in SEG], linewidths=3.2,
                        path_effects=None)
    ax.add_collection(lc)
    for m in meta:
        p = np.array(m["pts"]) / 1000.0
        ax.plot(p[1:-1, 0], p[1:-1, 1], "s", ms=3.2, mfc="w", mec="k", mew=0.5)
        ax.text(p[0, 0], p[0, 1], " " + m["name"], color="w", fontsize=9, fontweight="bold")
    ax.plot(IGN["x"] / 1000, IGN["y"] / 1000, "*", color="cyan", ms=13, mec="k")
    if crop:
        ax.set_xlim(crop[0], crop[1]); ax.set_ylim(crop[2], crop[3])
    n_b = int((at >= 0).sum()); ha = n_b * dx * dx / 1e4
    ax.set_title(f"{TITLE}   {label}   burned {ha * 2.47105:,.0f} acres", fontsize=10)
    miles_axes(ax)
    ax.text(0.99, 0.01, credit, transform=ax.transAxes, ha="right", va="bottom", fontsize=5.5, color="w",
            bbox=dict(facecolor="k", alpha=0.5, pad=1.2, linewidth=0))
    c1 = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.01); c1.solids.set_alpha(1); c1.set_label("fire arrival [min after ignition]")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmapL)
    c2 = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.035, pad=0.09); c2.set_label("conductor peak temperature [°F]")
    fig.tight_layout(); fig.savefig(fn, dpi=130); plt.close(fig)
    return at, dx


fr = os.path.join(a.out, f"frames_{a.case}"); os.makedirs(fr, exist_ok=True)
for old in glob.glob(fr + "/f_*.png"):
    os.remove(old)
k = 0
for pf in files[::a.every]:
    ds = yt.load(pf); t = float(ds.current_time)
    if t < T_IGN:
        continue
    draw(ds, span_T(t), f"{fr}/f_{k:04d}.png", f"{(t - T_IGN) / 60:4.0f} min after ignition"); k += 1
ds = yt.load(files[-1])
at, dx = draw(ds, {kk: v[0] for kk, v in peak.items()}, os.path.join(a.out, f"{a.case}_map_final.png"),
              "lines: peak temperature")

# ---------------------------------------------------------------- distances from the burned area
bi, bj = np.where(at >= 0)
bx, by = (bi + 0.5) * dx, (bj + 0.5) * dx


def dist_seg(p, q):
    if len(bx) == 0:
        return np.inf
    s = np.linspace(0, 1, 25)[:, None]
    pts = p[None, :] + s * (q - p)[None, :]
    d = np.hypot(pts[:, 0][:, None] - bx[None, :], pts[:, 1][:, None] - by[None, :])
    return float(d.min())


with open(os.path.join(a.out, f"{a.case}_spans.csv"), "w") as f:
    f.write("case,line,span,length_m,peak_T_C,peak_fire_flux_W_m2,closest_burn_m,fire_dx_m\n")
    for m in meta:
        p = np.array(m["pts"])
        for k in range(len(p) - 1):
            Tx, Ex = peak.get((m["name"], k), (np.nan, np.nan))
            f.write(f"{a.case},{m['name']},{k},{np.hypot(*(p[k + 1] - p[k])):.0f},{Tx:.2f},{Ex:.1f},{dist_seg(p[k], p[k + 1]):.0f},{dx:.1f}\n")
with open(os.path.join(a.out, f"{a.case}_towers.csv"), "w") as f:
    f.write("case,line,tower,x,y,peak_steel_T_C,closest_burn_m\n")
    for (ln, k), (tt, Tm) in sorted(towers.items()):
        m = next(m for m in meta if m["name"] == ln); p = np.array(m["pts"][k])
        d = float(np.hypot(bx - p[0], by - p[1]).min()) if len(bx) else np.inf
        f.write(f"{a.case},{ln},{k},{p[0]:.0f},{p[1]:.0f},{Tm.max():.2f},{d:.0f}\n")

# ---------------------------------------------------------------- movie (optional: needs ffmpeg)
if shutil.which("ffmpeg") and k > 0:
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "4", "-i", f"{fr}/f_%04d.png",
                    "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20",
                    os.path.join(a.out, f"{a.case}_map.mp4")], check=True)
else:
    print(f"no movie (ffmpeg not found, or no frame after ignition); frames in {fr}", file=sys.stderr)

# ---------------------------------------------------------------- temperatures against time
fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
for m in meta:
    sel = ct[ct["line"] == m["name"]]
    ts = np.unique(sel["t"])
    Tx = np.array([sel["Tx"][sel["t"] == t].max() for t in ts])
    axs[0].plot((ts - T_IGN) / 60, F(Tx), label=f"{m['name']}: hottest span")
    tw = [v for (ln, _), v in towers.items() if ln == m["name"]]
    if tw:
        t0 = tw[0][0]; Tt = np.max([np.interp(t0, tt, Tm) for tt, Tm in tw], axis=0)
        axs[1].plot((t0 - T_IGN) / 60, F(Tt), label=f"{m['name']}: hottest tower member")
axs[0].axhline(WIRE_FLAG_F, color="k", ls="--", lw=1.2)
axs[0].text(0.02, WIRE_FLAG_F, f" damage flag {WIRE_FLAG_F:.0f} °F (chosen, not a rating)", transform=axs[0].get_yaxis_transform(),
            va="bottom", fontsize=8)
axs[0].set_ylim(top=max(axs[0].get_ylim()[1], WIRE_FLAG_F * 1.08))
axs[1].text(0.02, 0.97, f"steel loses strength above {F(400.0):.0f} °F (EN 1993-1-2)", transform=axs[1].transAxes,
            va="top", fontsize=8)
axs[0].set_ylabel("conductor temperature [°F]"); axs[1].set_ylabel("tower steel temperature [°F]")
for x in axs:
    x.set_xlabel("minutes after ignition"); x.set_xlim(left=-5); x.grid(alpha=0.3); x.legend(fontsize=8, loc="center left")
fig.suptitle(TITLE, fontsize=10); fig.tight_layout()
fig.savefig(os.path.join(a.out, f"{a.case}_line_T.png"), dpi=150)
print("peak conductor T", max(v[0] for v in peak.values()), "peak steel T", max(v[1].max() for v in towers.values()) if towers else None)
