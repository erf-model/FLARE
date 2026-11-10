#!/usr/bin/env python3
"""Compare the two fire-under-a-line runs: the fire's heat into the air (feedback on) and not (off).

Reads, from each run directory, fire_stats.csv, L1_span<k>.dat, conductors/conductor_temperature.dat,
conductors/tower_temperature.dat and conductors/towers.dat, and writes

  summary.csv          one row per quantity: the on and off values and their difference
  compare_runs.png     time series of the fire, the conductors and the towers in both runs

    python3 compare_runs.py run_on run_off [--out .] [--t0 600]

Times are from the start of the simulation; the fire is lit at --t0.
"""
import argparse
import csv
import os

import numpy as np


def read_table(path, comment_header=False):
    """Whitespace table with a header line (optionally starting with '#'): dict of column -> array or list."""
    with open(path) as f:
        head = f.readline().lstrip("#").split()
        rows = [ln.split() for ln in f if ln.strip() and not ln.startswith("#")]
    cols = {}
    for j, name in enumerate(head):
        vals = [r[j] for r in rows]
        try:
            cols[name] = np.array(vals, dtype=float)
        except ValueError:
            cols[name] = np.array(vals)
    return cols


def read_csv(path):
    if not os.path.exists(path):
        raise SystemExit(f"{path}: not found")
    with open(path) as f:
        rd = csv.DictReader(f)
        rows = list(rd)
    if not rows:
        raise SystemExit(f"{path}: no rows (did the run reach the fire?)")
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def load_run(d):
    run = {"fire": read_csv(os.path.join(d, "fire_stats.csv"))}
    run["spans"] = {}
    k = 1
    while os.path.exists(os.path.join(d, f"L1_span{k}.dat")):
        run["spans"][k] = read_table(os.path.join(d, f"L1_span{k}.dat"))
        k += 1
    if not run["spans"]:
        raise SystemExit(f"{d}: no L1_span<k>.dat")
    run["ct"] = read_table(os.path.join(d, "conductors", "conductor_temperature.dat"))
    run["tt"] = read_table(os.path.join(d, "conductors", "tower_temperature.dat"))
    run["tw"] = read_table(os.path.join(d, "conductors", "towers.dat"))
    return run


def per_key(table, key_col, key, *cols):
    m = table[key_col] == key
    return [table[c][m] for c in cols]


def summary_rows(run, t0):
    out = {}
    f = run["fire"]
    out["burned area at the end (ha)"] = f["burned_area_ha"][-1]
    # the head rate's mean over the first 10 min, before the front nears the outflow face (where the
    # head's estimate jumps as the front leaves the domain)
    w = (f["time_s"] > t0) & (f["time_s"] <= t0 + 600.0)
    if not np.any(w):
        raise SystemExit(f"fire_stats.csv has no rows in the 10 min after ignition at {t0} s (--t0?)")
    out["mean head rate of spread, first 10 min (m/s)"] = f["head_ros_ms"][w].mean()
    ct = run["ct"]
    for k in sorted(run["spans"]):
        t, tmax, ef = per_key(ct, "span", k - 1, "time", "T_max_C", "E_fire_max_W_m2")
        sel = t >= t0
        before = tmax[t < t0][-1] if np.any(t < t0) else tmax[0]
        out[f"span {k}: conductor temperature before the fire (C)"] = before
        out[f"span {k}: hottest conductor temperature (C)"] = tmax[sel].max()
        out[f"span {k}: time of the hottest temperature (s)"] = t[sel][np.argmax(tmax[sel])]
        out[f"span {k}: largest fire irradiance on the conductor (W/m2)"] = ef[sel].max()
        s = run["spans"][k]
        sel_s = s["time"] >= t0
        sag0 = s["mid_sag"][s["time"] < t0][-1] if np.any(s["time"] < t0) else s["mid_sag"][0]
        out[f"span {k}: largest added mid-span sag (m)"] = s["mid_sag"][sel_s].max() - sag0
        out[f"span {k}: smallest clearance to the ground (m)"] = s["min_clearance"][sel_s].min()
    tt = run["tt"]
    tw = run["tw"]
    for tower in sorted(set(tt["tower"])):
        t, tmax, kemin, reb = per_key(tt, "tower", tower, "time", "T_max_C", "kE_min", "rebuilds")
        sel = t >= t0
        out[f"tower {tower}: hottest member steel (C)"] = tmax[sel].max()
        out[f"tower {tower}: smallest stiffness factor k_E"] = kemin[sel].min()
        out[f"tower {tower}: frame rebuilds"] = reb[-1]
        u = tw[f"{tower}_utilisation"][tw["time"] >= t0]
        out[f"tower {tower}: largest member utilisation"] = u.max()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_on")
    ap.add_argument("run_off")
    ap.add_argument("--out", default=".")
    ap.add_argument("--t0", type=float, default=600.0, help="ignition time (s)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    on, off = load_run(a.run_on), load_run(a.run_off)
    s_on, s_off = summary_rows(on, a.t0), summary_rows(off, a.t0)

    with open(os.path.join(a.out, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "feedback_on", "feedback_off", "on_minus_off"])
        for q in s_on:
            w.writerow([q, f"{s_on[q]:.6g}", f"{s_off[q]:.6g}", f"{s_on[q] - s_off[q]:.6g}"])
    for q in s_on:
        print(f"{q:60s} {s_on[q]:12.4g} {s_off[q]:12.4g}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
    style = {"on": dict(color="#b4471c", lw=1.6), "off": dict(color="#1f6f8b", lw=1.6, ls="--")}
    for name, run in (("on", on), ("off", off)):
        lab = f"feedback {name}"
        f = run["fire"]
        ax[0, 0].plot((f["time_s"] - a.t0) / 60, f["burned_area_ha"], label=lab, **style[name])
        t, tmax = per_key(run["ct"], "span", 1, "time", "T_max_C")
        ax[0, 1].plot((t - a.t0) / 60, tmax, label=lab, **style[name])
        s = run["spans"][2]
        ax[1, 0].plot((s["time"] - a.t0) / 60, s["min_clearance"], label=lab, **style[name])
        for i, tower in enumerate(sorted(set(run["tt"]["tower"]))):
            t, tmax = per_key(run["tt"], "tower", tower, "time", "T_max_C")
            ax[1, 1].plot((t - a.t0) / 60, tmax, label=f"{tower}, {lab}", alpha=1.0 if i == 0 else 0.6,
                          **style[name])
    ax[0, 0].set(ylabel="burned area (ha)", title="Fire")
    ax[0, 1].set(ylabel="hottest conductor temperature (°C)", title="Middle span: conductor")
    ax[1, 0].set(ylabel="smallest ground clearance (m)", title="Middle span: clearance")
    ax[1, 1].set(ylabel="hottest member steel (°C)", title="Towers: steel")
    for x in ax.flat:
        x.set_xlabel("time since ignition (min)")
        x.set_xlim(left=-1)
        x.grid(alpha=0.3)
        x.legend(fontsize=8)
    fig.savefig(os.path.join(a.out, "compare_runs.png"), dpi=150)


if __name__ == "__main__":
    main()
