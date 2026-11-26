#!/usr/bin/env python3
"""OSM power=line ways -> domain coordinates; list and plot them over the fuel map.

    osm_power_lines.py <osm_power.json> <case_domain.json> <fuel.asc> <out_prefix>

<osm_power.json> is an Overpass API answer for the domain's latitude-longitude box (S, W, N, E), for example
    q='[out:json][timeout:60];(way["power"="line"](S,W,N,E);way["power"="minor_line"](S,W,N,E);
       node["power"="tower"](S,W,N,E););out body;>;out skel qt;'
    curl -s --data-urlencode "data=$q" https://overpass-api.de/api/interpreter -o osm_power.json
OpenStreetMap data are (c) OpenStreetMap contributors and available under the ODbL
(https://www.openstreetmap.org/copyright); credit them wherever the lines are shown.

Writes <out_prefix>_lines.json (id, voltage, name, operator, [[x, y], ...] in domain metres, only ways with
at least one node inside the square) and <out_prefix>_lines.png.
"""
import json, sys
import numpy as np
from pyproj import Transformer
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

if len(sys.argv) != 5:
    sys.exit(__doc__.split("\n\n")[1])
osm, dom, fuel, out = sys.argv[1:5]
d = json.load(open(osm)); info = json.load(open(dom))
zone = info["utm_zone"]; L = float(info["domain_m"])
utm = f"EPSG:{(32600 if zone.upper().endswith('N') else 32700) + int(zone[:-1])}"
fwd = Transformer.from_crs("EPSG:4326", utm, always_xy=True)
x0, y0 = fwd.transform(*info["sw_corner_lonlat"])
nodes = {e["id"]: e for e in d["elements"] if e["type"] == "node"}
lines = []
for e in d["elements"]:
    if e["type"] != "way":
        continue
    t = e.get("tags", {})
    if t.get("power") not in ("line", "minor_line"):
        continue
    pts = []
    for n in e["nodes"]:
        if n in nodes:
            x, y = fwd.transform(nodes[n]["lon"], nodes[n]["lat"])
            pts.append([x - x0, y - y0, 1 if nodes[n].get("tags", {}).get("power") in ("tower", "pole") else 0])
    if not any(0 <= p[0] <= L and 0 <= p[1] <= L for p in pts):
        continue
    lines.append(dict(id=e["id"], power=t.get("power"), voltage=t.get("voltage"), name=t.get("name"),
                      operator=t.get("operator"), structure=t.get("line"), cables=t.get("cables"), pts=pts))
json.dump(lines, open(out + "_lines.json", "w"))

hdr = {}
with open(fuel) as f:
    for _ in range(6):
        k, v = f.readline().split(); hdr[k] = float(v)
F = np.loadtxt(fuel, skiprows=6)
burn = ~np.isin(F, [0, 91, 92, 93, 98, 99, -9999])
fig, ax = plt.subplots(figsize=(9, 9))
ax.imshow(burn, extent=[0, L, 0, L], cmap="Greens", alpha=0.5, origin="upper")
for l in lines:
    p = np.array(l["pts"])
    lw = 2.0 if l["power"] == "line" else 0.6
    ax.plot(p[:, 0], p[:, 1], "-", lw=lw, color="k" if l["power"] == "line" else "0.5")
    if l["power"] == "line":
        k = len(p) // 2
        ax.text(p[k, 0], p[k, 1], f'{l["id"]}\n{l["voltage"]}', fontsize=5, color="r")
ig = info.get("ignitions") or [info.get("ignition")]
for g in ig:
    ax.plot(g["x"], g["y"], "r*", ms=14)
ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_aspect("equal")
fig.savefig(out + "_lines.png", dpi=150, bbox_inches="tight")
print(len(lines), "ways in the domain;", sum(l["power"] == "line" for l in lines), "transmission")
