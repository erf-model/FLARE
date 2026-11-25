#!/usr/bin/env python3
"""FireBehaveFuelMap: BEHAVE's rate on a two-fuel map, against closed form.

    python3 check_behave_fuel_map.py              # every plt_fire_<variant>_* present
    python3 check_behave_fuel_map.py isotropic    # one variant

The decks run erf.fire.ros_model = "behave" with erf.fire.rothermel_per_fuel
= true on fuel_map_sh5_gr2.asc: Scott-Burgan SH5 (code 145) for x < 500 m and
GR2 (code 102) east of it. Dead moistures 3/4/5 %, live 60 %, a 6.5 m/s
midflame wind along +x (WAF off), flat ground, Rothermel's wind cap on.

The reference rates are the multi-class Rothermel (1972) model with Albini's
(1976) weighting, as Andrews (2018) gives it: characteristic SAV from the
area shares of all size classes, dead and live, and the live extinction
moisture from the per-class fine-fuel loads exp(-138/sigma) and
exp(-500/sigma). GR2's herbaceous load is two thirds cured at 60 % live
moisture (transfer window 0.30-1.20). The numbers come from an independent
script, not from the FLARE code:

    SH5 (145)  1.7584449850336146 m/s
    GR2 (102)  0.9489605773204224 m/s  (the 0.9 I_R wind cap binds)

Two line ignitions, 12.5 m half-width, at x = 12.5 m (SH5) and x = 512.5 m
(GR2) start two straight fronts at x = 25 m and x = 525 m with normal +x, so
each front crosses its own fuel at that fuel's head rate. The checks:

  SH5 half   at step 20 (10 s), every cell 100 m < x < 500 m has fire_ros
             = SH5 to 1e-8. fire_ros is the isotropic head rate, written on
             both paths. Before the per-fuel BEHAVE table every cell took
             the uniform fuel (erf.fire.fuel_model_id = 102), so this half
             read GR2's 0.949; with the dead-only SAV SH5 read 0.650.
  GR2 half   at step 20, every cell x > 500 m has fire_ros = GR2 to 1e-8.
  map read   fire_fuel_model_code is 145 west and 102 east of x = 500 m.
  SH5 front  at the last step (120 s), fire_arrival_time in the SH5 cells
             50 m < x < 225 m is (x - 25)/R_SH5 to half a cell crossing
             (7.1 s). On the level-set path this is the only check of the
             rate the front moves with (fill_directional_ros), which
             fire_ros does not show.
  GR2 front  fire_arrival_time in the GR2 cells 525 m < x < 625 m is
             (x - 525)/R_GR2 to half a cell crossing (13.2 s).

A pure-Python reader (erf_plotfile.py) reads the fire plotfile; no yt.
"""
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                                  # the CTest copy of the reader
sys.path.append(os.path.join(HERE, "..", "..", "..", "Submodules", "ERF", "Exec", "Tools"))  # from the source folder
import erf_plotfile  # noqa: E402

ROS_SH5 = 1.7584449850336146
ROS_GR2 = 0.9489605773204224
X_SEAM = 500.0
X_AHEAD = 100.0     # east of the 10 s front: the front normal is +x
REL = 1.0e-8
ROS_STEP = 20       # the plotfile the rate checks read
FIRE_DX = 25.0      # fire cell size [m]
# (fuel, rate, front start x, cell-centre window) for the arrival-time check
FRONTS = (("SH5", ROS_SH5, 25.0, (50.0, 225.0)),
          ("GR2", ROS_GR2, 525.0, (525.0, 625.0)))
VARIANTS = ("isotropic", "directional")

results = []


def check(name, ok, detail=""):
    ok = bool(ok)
    results.append(ok)
    print(f"  {name:20s} {'PASS' if ok else 'FAIL'}  {detail}")
    return ok


def check_arrival(v, plotfile):
    print(f"{v}: {plotfile} (arrival times)")
    hdr, data = erf_plotfile.read_fields(plotfile, ["fire_arrival_time"])
    at = data["fire_arrival_time"]
    x0, dx = hdr["prob_lo"][0], hdr["dx"][0]
    for key, ref, xs, (xa, xb) in FRONTS:
        tol = 0.5 * FIRE_DX / ref
        worst, n, unburnt = 0.0, 0, 0
        for i in range(len(at)):
            x = x0 + (i + 0.5) * dx
            if not (xa < x < xb):
                continue
            t_ref = (x - xs) / ref
            for j in range(len(at[i])):
                t = at[i][j][0]
                n += 1
                if t < 0.0:
                    unburnt += 1
                    continue
                worst = max(worst, abs(t - t_ref))
        check(f"{key} front", n > 0 and unburnt == 0 and worst <= tol,
              f"{n} cells, {unburnt} not reached, max |T - (x - {xs:g})/R| "
              f"{worst:.2f} s vs {tol:.2f} s")


def run_variant(v, plotfile=None):
    pfs = sorted(glob.glob(f"plt_fire_{v}_?????"))
    if plotfile is None:
        if not pfs:
            check(f"{v} plotfile", False, "missing")
            return
        last = pfs[-1]
        plotfile = f"plt_fire_{v}_{ROS_STEP:05d}"
        if not os.path.isdir(plotfile):
            check(f"{v} step {ROS_STEP}", False, "missing")
            return
        check_arrival(v, last)
    print(f"{v}: {plotfile}")
    hdr, data = erf_plotfile.read_fields(plotfile, ["fire_ros", "fire_fuel_model_code"])
    ros, code = data["fire_ros"], data["fire_fuel_model_code"]
    x0, dx = hdr["prob_lo"][0], hdr["dx"][0]
    worst = {"SH5": 0.0, "GR2": 0.0}
    seen = {"SH5": [], "GR2": []}
    codes_ok, n = True, {"SH5": 0, "GR2": 0}
    for i in range(len(ros)):
        x = x0 + (i + 0.5) * dx
        west = x < X_SEAM
        for j in range(len(ros[i])):
            c = int(round(code[i][j][0]))
            if c != (145 if west else 102):
                codes_ok = False
            if west and x < X_AHEAD:
                continue
            key, ref = ("SH5", ROS_SH5) if west else ("GR2", ROS_GR2)
            r = ros[i][j][0]
            worst[key] = max(worst[key], abs(r - ref) / ref)
            seen[key].append(r)
            n[key] += 1
    check("map read", codes_ok, "fire_fuel_model_code 145 | 102 across x = 500 m")
    for key, ref in (("SH5", ROS_SH5), ("GR2", ROS_GR2)):
        lo = min(seen[key]) if seen[key] else float("nan")
        hi = max(seen[key]) if seen[key] else float("nan")
        check(f"{key} half", n[key] > 0 and worst[key] <= REL,
              f"{n[key]} cells, fire_ros {lo:.10f}..{hi:.10f} vs {ref:.10f}, "
              f"max rel err {worst[key]:.2e}")


def main(argv):
    if len(argv) > 1:
        run_variant(argv[1], argv[2] if len(argv) > 2 else None)
    else:
        found = [v for v in VARIANTS if glob.glob(f"plt_fire_{v}_?????")]
        if not found:
            check("plotfiles", False, "no plt_fire_<variant>_* found")
        for v in found:
            run_variant(v)
    ok = bool(results) and all(results)
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
