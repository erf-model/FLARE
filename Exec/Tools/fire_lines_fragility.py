#!/usr/bin/env python3
"""Illustrative fire-exposure fragility from the demo runs, plus the numbers for the deck.

    fire_lines_fragility.py <figs_dir> <case>=<run_dir> [<case>=<run_dir> ...]

Inputs per case (from fire_lines_plot.py): <figs>/<case>_spans.csv, <case>_towers.csv; <run_dir>/fire_stats.csv.
IM: peak fire radiation reaching the item (kW/m^2). For a tower, the larger peak of the two spans it carries.
Damage flags: wire T > 100 C (chosen), tower steel T > 400 C (EN 1993-1-2 Table 3.1: k_y < 1 above 400 C).
Fit: lognormal CDF by maximum likelihood (Baker 2015), when both outcomes occur; otherwise the figure shows
peak temperature against exposure with the two flags. Writes fragility_fire.png and
summary.json (per case: burned_ha, head_speed, peak_cond_T, peak_steel_T, spans_burned, spans_total; pool: n, fits).
"""
import csv, json, os, sys
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

if len(sys.argv) < 3 or any("=" not in arg for arg in sys.argv[2:]):
    sys.exit(__doc__.split("\n\n")[1])
FIG = sys.argv[1]
RUNS = dict(arg.split("=", 1) for arg in sys.argv[2:])
cases = list(RUNS)
T_WIRE, T_STEEL = 100.0, 400.0          # deg C; the margin plot shows deg F
F = lambda c: c * 9.0 / 5.0 + 32.0
summ, items = {}, []
for c in cases:
    sp = list(csv.DictReader(open(os.path.join(FIG, f"{c}_spans.csv"))))
    tw = list(csv.DictReader(open(os.path.join(FIG, f"{c}_towers.csv"))))
    fs = np.genfromtxt(os.path.join(RUNS[c], "fire_stats.csv"), delimiter=",", names=True)
    t_ign = float([l for l in open(os.path.join(RUNS[c], "ignitions_spinup.csv")) if not l.startswith("#")][0].split()[0])
    after = fs[fs["time_s"] >= t_ign]
    # head speed: burned-area-weighted is noisy; report the median of head_ros from 15 min after ignition on
    late = after[after["time_s"] >= t_ign + 900]
    head = float(np.median(late["head_ros_ms"])) if len(late) else None
    flux = {}
    for r in sp:
        flux[(r["line"], int(r["span"]))] = float(r["peak_fire_flux_W_m2"]) / 1000.0
        items.append(("wire", c, flux[(r["line"], int(r["span"]))], float(r["peak_T_C"]) > T_WIRE, float(r["peak_T_C"])))
    for r in tw:
        k = int(r["tower"]); ln = r["line"]
        im = max(flux.get((ln, k - 1), 0.0), flux.get((ln, k), 0.0))
        items.append(("tower", c, im, float(r["peak_steel_T_C"]) > T_STEEL, float(r["peak_steel_T_C"])))
    # a span counts as reached when burned cells lie within one fire cell of it
    dx_burn = float(sp[0]["fire_dx_m"])
    summ[c] = dict(burned_ha=float(fs["burned_area_ha"][-1]), head_speed=head,
                   fire_min=round((float(fs["time_s"][-1]) - t_ign) / 60.0),
                   peak_cond_T=max(float(r["peak_T_C"]) for r in sp),
                   peak_steel_T=max((float(r["peak_steel_T_C"]) for r in tw), default=None),
                   spans_burned=sum(float(r["closest_burn_m"]) <= dx_burn for r in sp), spans_total=len(sp),
                   peak_flux_kW=max(flux.values()))


def fit(x, y):
    x, y = np.asarray(x), np.asarray(y, bool)
    if y.all() or (~y).all() or (x <= 0).any():
        return None
    def nll(p):
        th, b = np.exp(p)
        P = np.clip(norm.cdf(np.log(x / th) / b), 1e-12, 1 - 1e-12)
        return -np.sum(np.where(y, np.log(P), np.log(1 - P)))
    r = minimize(nll, [np.log(np.median(x)), np.log(0.5)], method="Nelder-Mead")
    th, b = np.exp(r.x)
    return dict(median_kW=float(th), beta=float(b))


fits = {k: fit([max(i[2], 1e-3) for i in items if i[0] == k], [i[3] for i in items if i[0] == k])
        for k in ("wire", "tower") if any(i[0] == k for i in items)}
if not any(fits.values()):
    # Nothing (or everything) passed a flag, so no curve can be fitted: show how close each item came instead.
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 3.95), sharex=True)
    for ax, kind, col, lim, lab in ((axs[0], "wire", "#c0392b", T_WIRE, "wire flag (chosen)"),
                                    (axs[1], "tower", "#2f6b55", T_STEEL, "steel flag (EN 1993-1-2)")):
        for c, mk in zip(cases, "os^"):
            sel = [i for i in items if i[0] == kind and i[1] == c]
            if sel:
                ax.semilogx([max(i[2], 0.1) for i in sel], [F(i[4]) for i in sel], mk, ms=4.5, color=col, alpha=0.7,
                            mfc=col if mk == "o" else "none", label=f"{c.capitalize()} ({len(sel)})")
        ax.axhline(F(lim), color="k", ls="--", lw=1.2)
        ax.text(0.12, F(lim) * 0.93, f"{lab}: {F(lim):.0f} °F", fontsize=7, va="top")
        ax.set_ylim(0, F(lim) * 1.08); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=7, loc="center left")
        ax.set_xlabel("peak fire radiation reaching it [kW/m²]")
        ax.set_title({"wire": "Wires: peak temperature", "tower": "Tower steel: peak temperature"}[kind], fontsize=9)
    axs[0].set_ylabel("°F"); axs[0].set_xlim(0.08, 400)
    fig.text(0.01, 0.005, "Items the fire did not reach are drawn at 0.1 kW/m².", fontsize=6.5, color="0.4")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fragility_fire.png"), dpi=170)
    summ["pool"] = dict(n=len(items), fits=fits, flagged={k: int(sum(i[3] for i in items if i[0] == k)) for k in ("wire", "tower")})
    json.dump(summ, open(os.path.join(FIG, "summary.json"), "w"), indent=1)
    print(json.dumps(summ, indent=1)); sys.exit(0)
fig, ax = plt.subplots(figsize=(6.6, 3.95))
fits = {}
for kind, col, flag in (("wire", "#c0392b", f"wire above {F(T_WIRE):.0f} °F (chosen flag)"),
                        ("tower", "#2f6b55", f"tower steel above {F(T_STEEL):.0f} °F (EN 1993-1-2)")):
    sel = [i for i in items if i[0] == kind]
    if not sel:
        continue
    x = np.array([max(i[2], 1e-3) for i in sel]); y = np.array([i[3] for i in sel])
    jit = (np.random.default_rng(1).random(len(y)) - 0.5) * 0.04
    ax.semilogx(x, y + jit, "o", ms=4, color=col, alpha=0.6, label=f"{kind}s: {y.sum()} of {len(y)} flagged")  # flag text below in F
    f = fit(x, y); fits[kind] = f
    if f:
        xx = np.logspace(np.log10(max(x.min(), 0.1) / 2), np.log10(x.max() * 2), 200)
        ax.semilogx(xx, norm.cdf(np.log(xx / f["median_kW"]) / f["beta"]), "-", color=col, lw=2,
                    label=f"fit: median {f['median_kW']:.1f} kW/m², dispersion {f['beta']:.2f}")
ax.set_xlabel("peak fire radiation reaching the span or tower [kW/m²]")
ax.set_ylabel("chance of passing the flag")
ax.set_ylim(-0.08, 1.08); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=7.5, loc="center right")
ax.set_title("Illustrative fire fragility, " + " and ".join(c.capitalize() for c in cases) + " runs pooled", fontsize=9.5)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fragility_fire.png"), dpi=170)
summ["pool"] = dict(n=len(items), fits=fits,
                    flagged={k: int(sum(i[3] for i in items if i[0] == k)) for k in ("wire", "tower")})
json.dump(summ, open(os.path.join(FIG, "summary.json"), "w"), indent=1)
print(json.dumps(summ, indent=1))
