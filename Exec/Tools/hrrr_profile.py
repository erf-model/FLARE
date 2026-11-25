#!/usr/bin/env python3
"""HRRR analysis column at points of a FLARE real-terrain domain -> ERF sounding + inflow file.

    hrrr_profile.py --date 2021-12-30 --hour 18 --domain marshall_domain.json --point-xy 0 11000 \
        --out-prefix marshall_hrrr18 [--top-km 8]

Fetches only the needed GRIB messages of the HRRR pressure-level analysis (product prs, fxx 0) with Herbie
(HGT, TMP, SPFH, UGRD, VGRD on isobaric levels) plus surface HGT and PRES (product sfc). U, V are rotated from
grid- to earth-relative (alpha = sin(Latin1)(lon - LoV)). Writes
  <prefix>_sounding.txt  ERF input sounding: p_floor[hPa] theta_sfc[K] qv_sfc[g/kg]; then z theta qv u v,
                         z above the local ground (use with erf.input_sounding_{wind,theta}_above_ground)
  <prefix>_inflow.txt    # z u v T above ground (xlo.inflow_profile = "file")
  <prefix>_info.json     ground height, 10 m wind, geostrophic estimate (mean wind 2-4 km AGL), points
p_floor is the HRRR surface pressure brought hydrostatically down to the domain floor (floor_elevation_m_asl).
"""
import argparse, json, math, os, warnings
import numpy as np
import pandas as pd
from herbie import Herbie
from pyproj import Transformer

ap = argparse.ArgumentParser()
ap.add_argument("--date", required=True)
ap.add_argument("--hour", type=int, required=True)
ap.add_argument("--domain", required=True)
ap.add_argument("--point-xy", nargs=2, type=float, action="append", required=True, metavar=("X", "Y"),
                help="domain coordinates (m); several points are averaged column by column")
ap.add_argument("--out-prefix", required=True)
ap.add_argument("--top-km", type=float, default=12.0)
ap.add_argument("--save-dir", default="hrrr_cache")
a = ap.parse_args()

info = json.load(open(a.domain))
zone = info["utm_zone"]
utm = f"EPSG:{(32600 if zone.upper().endswith('N') else 32700) + int(zone[:-1])}"
fwd = Transformer.from_crs("EPSG:4326", utm, always_xy=True)
inv = Transformer.from_crs(utm, "EPSG:4326", always_xy=True)
x0, y0 = fwd.transform(*info["sw_corner_lonlat"])
floor = float(info["floor_elevation_m_asl"])
pts = [inv.transform(x0 + x, y0 + y) for x, y in a.point_xy]

when = pd.Timestamp(f"{a.date} {a.hour:02d}:00")
SP = r":(?:HGT|TMP|SPFH|UGRD|VGRD):\d+ mb:"
SS = r":(?:HGT|PRES):surface:|:(?:UGRD|VGRD):10 m above ground:|:TMP:2 m above ground:|:SPFH:2 m above ground:"


def load(product, search):
    H = Herbie(when, model="hrrr", product=product, fxx=0, save_dir=a.save_dir, priority=["aws"], verbose=False)
    H.download(search)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        d = H.xarray(search, remove_grib=False)
    return d if isinstance(d, list) else [d]


def nearest(ds, lon, lat):
    la, lo = ds.latitude.values, ds.longitude.values
    lo = np.where(lo > 180, lo - 360, lo)
    k = np.argmin((la - lat) ** 2 + ((lo - lon) * math.cos(math.radians(lat))) ** 2)
    return np.unravel_index(k, la.shape)


prs = load("prs", SP)
sfc = load("sfc", SS)


def gridinfo(ds):
    for v in ds.data_vars:
        at = ds[v].attrs
        if "GRIB_LoVInDegrees" in at:
            return at["GRIB_LoVInDegrees"], at["GRIB_Latin1InDegrees"]
    return -97.5, 38.5


cols = []
for lon, lat in pts:
    c = {}
    for ds in prs:
        if "isobaricInhPa" not in ds.dims and "isobaricInhPa" not in ds.coords:
            continue
        ij = nearest(ds, lon, lat)
        lov, lat1 = gridinfo(ds)
        p = np.atleast_1d(ds["isobaricInhPa"].values)
        for v in ds.data_vars:
            c[v] = np.atleast_1d(ds[v].values[(..., ) + ij])
        c["p"] = p
    alpha = math.sin(math.radians(lat1)) * math.radians(lon - (lov - 360 if lov > 180 else lov))
    ug, vg = c["u"], c["v"]
    c["u_earth"] = ug * math.cos(alpha) + vg * math.sin(alpha)
    c["v_earth"] = -ug * math.sin(alpha) + vg * math.cos(alpha)
    s = {}
    for ds in sfc:
        ij = nearest(ds, lon, lat)
        for v in ds.data_vars:
            lev = ds[v].attrs.get("GRIB_typeOfLevel", "")
            s[f"{v}_{lev}"] = float(ds[v].values[ij])
    c["sfc"] = s
    u10, v10 = s.get("u10_heightAboveGround"), s.get("v10_heightAboveGround")
    c["u10e"] = u10 * math.cos(alpha) + v10 * math.sin(alpha)
    c["v10e"] = -u10 * math.sin(alpha) + v10 * math.cos(alpha)
    cols.append(c)

# average the columns on a common above-ground height grid
zg = np.concatenate([np.arange(5.0, 200.0, 10.0), np.arange(200.0, 1000.0, 25.0), np.arange(1000.0, a.top_km * 1000 + 1, 100.0)])
acc = {k: [] for k in ("th", "qv", "u", "v", "T")}
grounds, psfc, t2, q2, u10l, v10l = [], [], [], [], [], []
for c in cols:
    zs = c["sfc"]["orog_surface"] if "orog_surface" in c["sfc"] else c["sfc"]["gh_surface"]
    ps = c["sfc"]["sp_surface"]
    z = c["gh"] - zs
    keep = z > 2.0
    o = np.argsort(z[keep])
    zz = z[keep][o]
    T = c["t"][keep][o]; p = c["p"][keep][o]; q = c["q"][keep][o]
    th = T * (1000.0 / p) ** 0.2854
    T2 = c["sfc"]["t2m_heightAboveGround"]; Q2 = c["sfc"]["sh2_heightAboveGround"]
    th2 = T2 * (1e5 / ps) ** 0.2854
    # near-surface anchor: 2 m theta/q and the 10 m wind, then the isobaric levels
    zz = np.concatenate([[2.0, 10.0], zz]); th = np.concatenate([[th2, th2], th])
    q = np.concatenate([[Q2, Q2], q]); T = np.concatenate([[T2, T2], T])
    u = np.concatenate([[c["u10e"] * 0.6, c["u10e"]], c["u_earth"][keep][o]])
    v = np.concatenate([[c["v10e"] * 0.6, c["v10e"]], c["v_earth"][keep][o]])
    for k, arr in (("th", th), ("qv", q / (1 - q)), ("u", u), ("v", v), ("T", T)):
        acc[k].append(np.interp(zg, zz, arr))
    grounds.append(zs); psfc.append(ps); t2.append(T2); q2.append(Q2); u10l.append(c["u10e"]); v10l.append(c["v10e"])
prof = {k: np.mean(v, axis=0) for k, v in acc.items()}
ground = float(np.mean(grounds)); ps = float(np.mean(psfc)); T2 = float(np.mean(t2))
# hydrostatic to the floor with the 2 m temperature
p_floor = ps * math.exp(9.81 * (ground - floor) / (287.04 * T2))
theta_floor = T2 * (1e5 / p_floor) ** 0.2854

with open(f"{a.out_prefix}_sounding.txt", "w") as f:
    f.write(f"{p_floor / 100:.2f} {prof['th'][0]:.3f} {prof['qv'][0] * 1000:.4f}\n")
    for i, z in enumerate(zg):
        f.write(f"{z:10.3f} {prof['th'][i]:10.4f} {prof['qv'][i] * 1000:8.4f} {prof['u'][i]:10.4f} {prof['v'][i]:10.4f}\n")
with open(f"{a.out_prefix}_inflow.txt", "w") as f:
    f.write(f"# HRRR analysis {a.date} {a.hour:02d} UTC, mean of {len(pts)} column(s) at {a.point_xy}\n")
    f.write("# Heights are above the local ground; T is potential temperature (K)\n# z u v T\n")
    for i, z in enumerate(zg):
        f.write(f"{z:10.3f} {prof['u'][i]:10.5f} {prof['v'][i]:10.5f} {prof['th'][i]:10.4f}\n")
sel = (zg >= 2000) & (zg <= 4000)
ug, vg = float(prof["u"][sel].mean()), float(prof["v"][sel].mean())
u10, v10 = float(np.mean(u10l)), float(np.mean(v10l))
out = dict(date=a.date, hour_utc=a.hour, points_xy=a.point_xy, points_lonlat=pts, ground_m_asl=ground,
           psfc_hPa=ps / 100, p_floor_hPa=p_floor / 100, T2_K=T2,
           u10=u10, v10=v10, speed10=math.hypot(u10, v10), dir10_from=(math.degrees(math.atan2(-u10, -v10)) % 360),
           geo_2_4km=[ug, vg], geo_speed=math.hypot(ug, vg), geo_dir_from=(math.degrees(math.atan2(-ug, -vg)) % 360),
           hpbl_note="see hrrr_point_winds.py for HPBL")
json.dump(out, open(f"{a.out_prefix}_info.json", "w"), indent=1)
print(json.dumps(out, indent=1))
