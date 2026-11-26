#!/usr/bin/env python3
"""Head distance from the ignition point against time, measured the way the observed anchors are.

    fire_head_distance.py <figs_dir> <case>=<run_dir> [<case>=<run_dir> ...]

Reads the last plt_fire_* of each <run_dir>. fire_arrival_time is model time (s, -1 = not burned). The head
distance at time t is the largest distance from the ignition point (ignitions_spinup.csv) among cells that arrived by t.
Writes <figs>/<case>_head.csv (min after ignition, head_km) and adds head_km (5-min list) to <figs>/summary.json.
Run after fire_lines_fragility.py, which rewrites summary.json.
"""
import csv, glob, json, os, sys
import numpy as np
import yt
yt.set_log_level(50)

if len(sys.argv) < 3 or any("=" not in arg for arg in sys.argv[2:]):
    sys.exit(__doc__.split("\n\n")[1])
FIG = sys.argv[1]
RUNS = dict(arg.split("=", 1) for arg in sys.argv[2:])
sp = os.path.join(FIG, "summary.json")
summ = json.load(open(sp)) if os.path.exists(sp) else {}
for c, run in RUNS.items():
    t_ign, x0, y0 = map(float, [l for l in open(os.path.join(run, "ignitions_spinup.csv"))
                                if not l.startswith("#")][0].split()[:3])
    ds = yt.load(sorted(glob.glob(os.path.join(run, "plt_fire_?????")))[-1])
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    at = np.array(g["fire_arrival_time"])[:, :, 0]
    lo, hi = ds.domain_left_edge.d, ds.domain_right_edge.d
    dx = (hi[0] - lo[0]) / at.shape[0]
    i, j = np.where(at >= 0)
    d = np.hypot(lo[0] + (i + 0.5) * dx - x0, lo[1] + (j + 0.5) * dx - y0) / 1000.0
    ta = (at[i, j] - t_ign) / 60.0
    t_end = (float(ds.current_time) - t_ign) / 60.0
    rows = [(m, float(d[ta <= m].max()) if (ta <= m).any() else 0.0) for m in range(0, int(t_end) + 1, 5)]
    with open(os.path.join(FIG, f"{c}_head.csv"), "w") as f:
        w = csv.writer(f); w.writerow(["min_after_ignition", "head_km"]); w.writerows(rows)
    summ.setdefault(c, {})["head_km"] = [[m, round(h, 3)] for m, h in rows]
    print(c, "t_end %.1f min" % t_end, " ".join(f"{m}:{h:.2f}" for m, h in rows))
json.dump(summ, open(sp, "w"), indent=1)
