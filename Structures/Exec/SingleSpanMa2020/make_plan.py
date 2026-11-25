#!/usr/bin/env python3
"""Write the multiple-stripes plan of the single span of Ma, Bocchini and Christou 2020 for the
fragility ensemble driver (Structures/Tools/fragility_ensemble) and its fit
(Structures/Tools/fragility_fit).

    python3 make_plan.py --out stripes --winds 40 50 60 70 80 --seeds 1 2 3 4 5
    python3 ../../Tools/fragility_ensemble/fragility_ensemble.py run --plan stripes/plan.json --exe path/to/flare --ranks 2
    python3 ../../Tools/fragility_ensemble/fragility_ensemble.py collect --plan stripes/plan.json
    python3 ../../Tools/fragility_fit/fragility_fit.py --plan stripes/fit_plan.json

The base deck is make_case.py --mode flow at the first wind (written to <out>/base). Every stripe
keeps the base's Courant number: its erf.fixed_dt is the base's times first wind / stripe, so the
shared spin-up of --spinup-flows flow-throughs of the domain at 10 m is the same number of steps
in every stripe; each member then runs 60 s more and the --average seconds the statistics are
taken over. The conductor is checked against Drake's rated strength, 140.18 kN (Ma 2020 Table 3);
nothing lets go.

The fit's intensity is the stripe: the inflow's wind at 10 m, a mean over the whole run. Ma et
al.'s intensity measure is V_ms, the highest 1-min mean wind at 10 m, which is higher than a
long mean by a gust factor that is not applied here. The run also records the wind at 10 m
averaged over --average s (conductors.intensity_measure = wind); that record is the largest such
mean over windows ending at or after conductors.capacity_check_start, so windows ending before the
run's end reach back into the settling period, and the 10 m wind of the first cell is below the
inflow's. It is kept for reference and not used by the fit.

The inflow profile reaches 1.6 times the wind at 10 m at the domain's top (240 m); the driver
takes profiles up to 150 m/s, so the winds are limited to 93 m/s.
"""
import argparse
import json
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LX = 600.0      # make_case.py --mode flow: the domain's length along the wind (m)
LZ = 240.0      # make_case.py --mode flow: the domain's height (m)
Z0 = 0.05       # make_case.py --mode flow: the default roughness length (m)
PROFILE_MAX = 150.0   # the driver's largest wind in a scaled profile (m/s)


def fmt(x):
    return float("%.6g" % x)


def deck_value(path, key):
    with open(path) as f:
        for ln in f:
            s = ln.split("#", 1)[0]
            if "=" in s and s.split("=", 1)[0].strip() == key:
                return s.split("=", 1)[1].split()
    raise SystemExit("make_plan: %s does not set %s" % (path, key))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="the plan's directory")
    ap.add_argument("--winds", type=float, nargs="+", required=True, help="stripes: the inflow's wind at 10 m (m/s), each in (0, 93]")
    ap.add_argument("--seeds", type=int, nargs="+", required=True, help="conductors.gust_seed of the members, each >= 0")
    ap.add_argument("--span", type=float, default=300.0, help="span (m), in [50, 1000]")
    ap.add_argument("--EA", type=float, default=0.0, help="axial stiffness (N), >= 0; 0: make_case.py's default")
    ap.add_argument("--without-fft", action="store_true", help="erf.use_fft = false, for an ERF built without FFT")
    ap.add_argument("--spinup-flows", type=float, default=3.0, help="flow-throughs of the spin-up at 10 m, in [1, 20]")
    ap.add_argument("--average", type=float, default=600.0, help="averaging time (s) of the wind at 10 m and of the statistics, in (0, 600]")
    a = ap.parse_args(argv)
    top = math.log((LZ + Z0) / Z0) / math.log((10.0 + Z0) / Z0)
    vmax = math.floor(PROFILE_MAX / top)
    if not all(0.0 < u <= vmax for u in a.winds) or len(set(a.winds)) != len(a.winds):
        ap.error("--winds: each in (0, %g] m/s (the profile's top, %.3f times the wind at 10 m, stays within the "
                 "driver's %g m/s), no repeats" % (vmax, top, PROFILE_MAX))
    if not 50.0 <= a.span <= 1000.0:
        ap.error("--span must be in [50, 1000] m")
    if a.EA < 0.0:
        ap.error("--EA must be >= 0 N")
    if not all(s >= 0 for s in a.seeds) or len(set(a.seeds)) != len(a.seeds):
        ap.error("--seeds: each >= 0, no repeats")
    if not 1.0 <= a.spinup_flows <= 20.0:
        ap.error("--spinup-flows must be in [1, 20]")
    if not 0.0 < a.average <= 600.0:
        ap.error("--average must be in (0, 600] s")
    winds = sorted(a.winds)
    v0 = winds[0]
    out = os.path.abspath(a.out)
    base = os.path.join(out, "base")
    spin_time = a.spinup_flows * LX / v0
    cmd = [sys.executable, "-I", os.path.join(HERE, "make_case.py"), "--out", base, "--mode", "flow", "--span", str(a.span),
           "--wind", str(v0), "--spinup", str(spin_time), "--average", str(a.average)]
    if a.EA > 0.0:
        cmd += ["--EA", str(a.EA)]
    if a.without_fft:
        cmd += ["--without-fft"]
    subprocess.run(cmd, check=True)
    deck = os.path.join(base, "inputs")
    dt0 = float(deck_value(deck, "erf.fixed_dt")[0])
    spin_steps = int(math.ceil(spin_time / dt0))
    speeds = []
    for v in winds:
        dt = fmt(dt0 * v0 / v)
        max_step = spin_steps + int(math.ceil((60.0 + a.average) / dt))
        start = math.floor((max_step * dt - a.average) * 1.0e3) / 1.0e3   # the statistics span at least --average s
        speeds.append({"speed": v, "inputs": {
            "erf.fixed_dt": dt,
            "max_step": max_step,
            "stop_time": fmt(max_step * dt + 10.0 * dt),
            "conductors.diagnostics_int": max(1, int(round(1.0 / dt))),
            "conductors.stats_start": start,
            "conductors.capacity_check_start": start}})
    plan = {"base_deck": os.path.join("base", "inputs"), "output_dir": "ensemble", "speeds": speeds,
            "seeds": a.seeds, "reference_height": 10,
            "spinup": {"max_step": spin_steps}, "inputs": {"conductors.failure_break": False}}
    with open(os.path.join(out, "plan.json"), "w") as f:
        json.dump(plan, f, indent=2)
        f.write("\n")
    fit = {"demand_table": os.path.join("ensemble", "demand_table_all.csv"), "output_dir": "fit", "intensity": "stripe"}   # the inflow's wind at 10 m (docstring)
    with open(os.path.join(out, "fit_plan.json"), "w") as f:
        json.dump(fit, f, indent=2)
        f.write("\n")
    print("%s: %d stripes x %d seeds; spin-up %d steps (%.1f s at %g m/s)" % (out, len(winds), len(a.seeds), spin_steps, spin_time, v0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
