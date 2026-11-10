#!/usr/bin/env python3
"""Snapshots of the two runs, each when its middle span is hottest: the fire in plan view and the plume through
the middle span, side by side for feedback on and off.

    python3 make_figure.py run_on run_off [--out fire_under_line.png]

Each run directory holds the plotfiles plt<step> and plt_fire_<step> (every 60 s), L1_nodes.dat and
conductors/conductor_temperature.dat. Needs yt.
"""
import argparse
import glob
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import FancyArrow
import yt

yt.set_log_level(50)
LX, LY, LZ = 4000.0, 2000.0, 600.0
LINE_X, TOWERS_Y, SPAN_MID_Y = 1200.0, (850.0, 1150.0), 1000.0
T_IGN = 600.0


def table(path):
    with open(path) as f:
        head = f.readline().lstrip("#").split()
        rows = [ln.split() for ln in f if ln.strip() and not ln.startswith("#")]
    return head, rows


def hottest_time(run):
    head, rows = table(os.path.join(run, "conductors", "conductor_temperature.dat"))
    c = {n: i for i, n in enumerate(head)}
    mid = [(float(r[c["T_max_C"]]), float(r[c["time"]])) for r in rows if int(r[c["span"]]) == 1]
    return max(mid)[1]


def nearest_plotfile(run, prefix, t):
    files = sorted(glob.glob(os.path.join(run, f"{prefix}[0-9]*")))
    if not files:
        raise SystemExit(f"{run}: no {prefix}<step> plotfiles")
    times = [float(yt.load(f).current_time) for f in files]
    i = int(np.argmin(np.abs(np.array(times) - t)))
    return files[i], times[i]


def grid(path, *fields):
    ds = yt.load(path)
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    return [np.array(g[("boxlib", f)]) for f in fields], ds.domain_dimensions


def nodes_at(run, t):
    by_t = {}
    for ln in open(os.path.join(run, "L1_nodes.dat")):
        f = ln.split()
        if f[0] == "time":
            continue
        by_t.setdefault(float(f[0]), []).append((int(f[1]), float(f[2]), float(f[3]), float(f[4])))
    tn = min(by_t, key=lambda s: abs(s - t))
    return sorted(by_t[tn])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_on")
    ap.add_argument("run_off")
    ap.add_argument("--out", default="fire_under_line.png")
    a = ap.parse_args()

    fig, axs = plt.subplots(2, 2, figsize=(15, 7.6), gridspec_kw=dict(height_ratios=[1.0, 1.0]),
                            constrained_layout=True)
    hr_max = 0.0
    snaps = []
    for run in (a.run_on, a.run_off):
        f_fire, t = nearest_plotfile(run, "plt_fire_", hottest_time(run))
        f_air, _ = nearest_plotfile(run, "plt", t)
        (hr, arr), (nfx, nfy, _) = grid(f_fire, "fire_heat_release", "fire_arrival_time")
        (th, u, w), (nx, ny, nz) = grid(f_air, "theta", "x_velocity", "z_velocity")
        snaps.append((run, t, hr.squeeze(), arr.squeeze(), th, u, w, (nfx, nfy), (nx, ny, nz)))
        hr_max = max(hr_max, hr.max())

    norm_T = Normalize(30.0, 80.0)
    cmap_T = plt.get_cmap("coolwarm")
    for col, (run, t, hr, arr, th, u, w, (nfx, nfy), (nx, ny, nz)) in enumerate(snaps):
        label = "fire heat into the air (feedback on)" if col == 0 else "no fire heat into the air (feedback off)"
        ct_head, ct_rows = table(os.path.join(run, "conductors", "conductor_temperature.dat"))
        c = {n: i for i, n in enumerate(ct_head)}
        span_T = {int(r[c["span"]]): float(r[c["T_max_C"]]) for r in ct_rows if abs(float(r[c["time"]]) - t) < 0.26}
        x = (np.arange(nfx) + 0.5) * LX / nfx
        y = (np.arange(nfy) + 0.5) * LY / nfy
        ax = axs[0, col]
        burned = np.where((arr >= 0) & (arr <= t), 1.0, np.nan)
        ax.pcolormesh(x, y, burned.T, cmap="Greys", vmin=0, vmax=3, shading="auto", rasterized=True)
        hrm = np.ma.masked_less_equal(hr, 0.0)
        pc = ax.pcolormesh(x, y, hrm.T / 1e3, cmap="inferno", vmin=0, vmax=hr_max / 1e3, shading="auto",
                           rasterized=True)
        arrm = np.ma.masked_where((arr < 0) | (arr > t), arr)
        lv = np.arange(T_IGN + 120.0, t + 1e-6, 120.0)
        if len(lv):
            cs = ax.contour(x, y, arrm.T, levels=lv, colors="w", linewidths=0.7)
            ax.clabel(cs, fmt=lambda v: f"{(v - T_IGN) / 60:.0f} min", fontsize=7)
        ax.plot([600, 600], [700, 1300], color="tab:green", lw=2.5, label="ignition line")
        nd = nodes_at(run, t)
        seg = 21  # L1_nodes.dat numbers span k's nodes 21k .. 21k + 20 (20 segments), then the strings'
        for k in range(3):
            pts = [p for p in nd if k * seg <= p[0] < (k + 1) * seg]
            ax.plot([p[1] for p in pts], [p[2] for p in pts], color=cmap_T(norm_T(span_T.get(k, 30.0))), lw=4,
                    solid_capstyle="round")
        for ty in TOWERS_Y:
            ax.plot(LINE_X, ty, marker="s", ms=9, mfc="0.85", mec="k", mew=1.2)
        ax.add_patch(FancyArrow(150, 1880, 250, 0, width=8, head_width=40, head_length=40, color="tab:blue"))
        ax.text(150, 1790, "inflow 6 m/s at 10 m", color="tab:blue", fontsize=9)
        ax.set(xlim=(0, LX), ylim=(0, LY), aspect="equal", xlabel="x (m)", ylabel="y (m)")
        ax.set_title(f"{label}\n{(t - T_IGN) / 60:.1f} min after ignition, middle span at its hottest "
                     f"({span_T.get(1, float('nan')):.0f} °C)", fontsize=10.5)
        if col == 0:
            ax.legend(loc="lower right", fontsize=8)
        if col == 1:
            cb = fig.colorbar(pc, ax=axs[0, :], fraction=0.025, pad=0.01)
            cb.set_label("fire heat release (kW/m$^2$)")
            sm = plt.cm.ScalarMappable(norm=norm_T, cmap=cmap_T)
            cb2 = fig.colorbar(sm, ax=axs[0, :], fraction=0.025, pad=0.01)
            cb2.set_label("conductor, hottest node per span (°C)")

        ax2 = axs[1, col]
        j = int(round(SPAN_MID_Y / (LY / ny) - 0.5))   # the cell centre nearest the middle of the middle span
        xa = (np.arange(nx) + 0.5) * LX / nx
        za = (np.arange(nz) + 0.5) * LZ / nz
        dth = th[:, j, :] - 300.0
        pc2 = ax2.pcolormesh(xa, za, dth.T, cmap="magma", vmin=0, vmax=10.0, shading="auto", rasterized=True)
        kz = za < 300
        ax2.quiver(xa[::3], za[kz][::2], u[::3, j, :][:, kz][:, ::2].T, w[::3, j, :][:, kz][:, ::2].T,
                   color="w", scale=200, width=0.0018, alpha=0.7)
        zmid = [p[3] for p in nd if 21 <= p[0] <= 41]
        ax2.plot([LINE_X], [min(zmid)], marker="o", ms=7, color="cyan", mec="k", label="middle span, lowest point")
        ax2.set(xlim=(400, 3000), ylim=(0, 300), xlabel="x (m)", ylabel="z (m)")
        ax2.set_title(f"Air through the middle of the middle span (y = {(j + 0.5) * LY / ny:.1f} m)", fontsize=10.5)
        ax2.legend(loc="upper left", fontsize=8)
        if col == 1:
            cb3 = fig.colorbar(pc2, ax=axs[1, :], fraction=0.025, pad=0.01)
            cb3.set_label(r"potential temperature above ambient, $\theta - 300$ K")

    fig.suptitle("Grass fire under a transmission line: with and without the fire's heat in the air",
                 fontsize=14, weight="bold")
    fig.savefig(a.out, dpi=150)
    print(a.out)


if __name__ == "__main__":
    main()
