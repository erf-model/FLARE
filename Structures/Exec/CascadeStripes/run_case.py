#!/usr/bin/env python3
"""Two wind stripes through the fragility ensemble driver with the line failures on, and the damage
states DS1-DS4 from the fragility fit.

    python3 run_case.py --exe path/to/flare --out case_dir [--mpiexec mpiexec] [--ranks 1]
                        [--uplift 20000] [--load-start 6 | --no-load-gate] [--check]

The lines are Tests/test_files/Conductors_RandomGusts (two lines over hills, three lattice towers,
random gusts). The driver runs two stripes, 10 and 25 m/s at 30 m, two gust seeds each, and one shared
spin-up of 5 steps of the flow alone per stripe; every member restarts from its stripe's checkpoint
with the lines, the failures (conductors.failure_break, ground_contact, tower_failure = release) and
runs on to step 15. The towers' footings hold a SYNTHETIC uplift allowable (--uplift, default 20 kN)
and every other capacity is 1e9, so the footing uplift alone decides the outcome. The fit then counts
each damage state per stripe (DS4: 2 or more towers lost).

The allowable comes from a calibration run with every capacity at 1e9 (README.md): the largest
footing uplift was 0 at 10 m/s in both seeds and 33 to 66 kN at 25 m/s. 20 kN is a value chosen
between the two for the test, not a design value.

With --check the case ends with return code 1 unless:
  1. every member restarted from its stripe's spin-up checkpoint, two members per checkpoint; the
     conductor model reports that the checkpoint holds no conductor state (the lines were made at
     the restart, not carried in it), and the conductor logs start at the end of the spin-up
     (5 x erf.fixed_dt);
  2. at 10 m/s no member exceeds a capacity, releases a line or loses a tower;
  3. at 25 m/s every member's first failure is a footing uplift above the allowable, and it releases
     a line and loses 2 or more towers;
  4. the fit's points give 0 of 2 members in every state at 10 m/s and 2 of 2 at 25 m/s.

Each member is one realisation: a falling span amplifies roundoff, so the number of towers lost at
25 m/s can differ between builds and rank counts; the checks ask only for 2 or more.
Needs an executable built with ERF_ENABLE_STRUCTURES (the MoorDyn stub or MoorDyn-C). Standard
library only.
"""
import argparse
import csv
import json
import os
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
DECK_DIR = os.path.join(ROOT, "Tests", "test_files", "Conductors_RandomGusts")
DECK = "Conductors_RandomGusts.i"
DRIVER = os.path.join(ROOT, "Structures", "Tools", "fragility_ensemble", "fragility_ensemble.py")
FIT = os.path.join(ROOT, "Structures", "Tools", "fragility_fit", "fragility_fit.py")

SPEEDS = (10, 25)
SEEDS = (7, 8)
SPINUP_STEPS = 5
MAX_STEP = 15
CASCADE_TOWERS = 2
STATES = ("DS1", "DS2", "DS3", "DS4")


def driver_plan(uplift):
    return {
        "base_deck": os.path.join(DECK_DIR, DECK),
        "output_dir": "ensemble",
        "speeds": list(SPEEDS),
        "seeds": list(SEEDS),
        "reference_height": 30,
        "support_files": ["terrain_hills.txt"],
        "spinup": {"max_step": SPINUP_STEPS},
        "inputs": {
            "max_step": MAX_STEP,
            "stop_time": 13.6,
            "erf.plot_int_1": -1,
            "conductors.failure_break": True,
            "conductors.ground_contact": True,
            "conductors.tower_failure": "release",
            "conductors.lattice.allowable_compression": 1.0e9,
            "conductors.lattice.allowable_uplift": uplift,
            "conductors.L1.rated_breaking_strength": 1.0e9,
            "conductors.L2.rated_breaking_strength": 1.0e9,
        },
    }


def fit_plan():
    return {
        "demand_table": "ensemble/demand_table_all.csv",
        "output_dir": "fit",
        "resamples": 0,
        "damage_states": {"table": "ensemble/damage_table_all.csv", "states": list(STATES),
                          "cascade_towers": CASCADE_TOWERS},
    }


def write_json(path, obj):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)
        f.write("\n")


def rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def data_rows(path):
    """The rows of a whitespace log after its header line (lines starting with # are notes)."""
    with open(path, encoding="utf-8") as f:
        lines = [ln.split() for ln in f if ln.strip() and not ln.startswith("#")]
    return lines[1:]


def deck_value(path, key):
    with open(path, encoding="utf-8") as f:
        for ln in f:
            parts = ln.split("#", 1)[0].split("=", 1)
            if len(parts) == 2 and parts[0].strip() == key:
                return parts[1].strip()
    return None


def check(out, uplift):
    bad = []
    ens = os.path.join(out, "ensemble")
    manifest = json.load(open(os.path.join(ens, "manifest.json"), encoding="utf-8"))
    members = manifest["members"]

    # 1. restart from the spin-up, lines made at the restart
    per_spinup = {}
    for m in members:
        d = os.path.join(ens, m["directory"])
        want = "../../spinup/%s/chk%05d" % (m["spinup"], SPINUP_STEPS)
        per_spinup[m["spinup"]] = per_spinup.get(m["spinup"], 0) + 1
        log = open(os.path.join(d, "run.log"), encoding="utf-8", errors="replace").read()
        if "Restart from native checkpoint " + want not in log:
            bad.append("%s did not restart from %s" % (m["member"], want))
        # a new member directory always starts its logs at the restart, so the logs alone cannot show
        # where the lines came from; the conductor model says so in the run log
        fresh = "conductors: the checkpoint %s holds no conductor state; the lines start at this time" % want
        if fresh not in log:
            bad.append("%s: the run log does not say '%s' (the lines came out of the checkpoint)"
                       % (m["member"], fresh))
        dt = float(deck_value(os.path.join(d, DECK), "erf.fixed_dt"))
        start = SPINUP_STEPS * dt
        towers = data_rows(os.path.join(d, "conductors", "towers.dat"))
        if not towers or abs(float(towers[0][0]) - start) > 1e-6 * start:
            bad.append("%s: towers.dat starts at %s s, not at the end of the spin-up (%g s)"
                       % (m["member"], towers[0][0] if towers else "no row", start))
    if sorted(per_spinup.values()) != [len(SEEDS)] * len(SPEEDS):
        bad.append("members per spin-up %s, expected %d spin-ups of %d" % (per_spinup, len(SPEEDS), len(SEEDS)))

    # 2 and 3. outcomes per stripe
    for r in rows(os.path.join(ens, "damage_table_all.csv")):
        name, stripe = r["member"], float(r["stripe"])
        ratio, released, lost = float(r["max_peak_ratio"]), int(r["lines_released"]), int(r["towers_lost"])
        if stripe == SPEEDS[0]:
            if ratio >= 1.0 or released or lost:
                bad.append("%s (10 m/s): peak ratio %g, %d line(s) released, %d tower(s) lost; expected none"
                           % (name, ratio, released, lost))
            continue
        if released < 1 or lost < CASCADE_TOWERS:
            bad.append("%s (25 m/s): %d line(s) released, %d tower(s) lost; expected at least 1 and %d"
                       % (name, released, lost, CASCADE_TOWERS))
        fails = data_rows(os.path.join(ens, "members", name, "conductors", "failures.dat"))
        if not fails:
            bad.append("%s (25 m/s): failures.dat has no row" % name)
            continue
        first = fails[0]        # time kind component demand capacity ratio
        if (first[1] != "footing_uplift" or abs(float(first[4]) - uplift) > 1e-9 * uplift
                or float(first[3]) <= uplift):
            bad.append("%s (25 m/s): first failure %s; expected a footing_uplift above %g N"
                       % (name, " ".join(first), uplift))

    # 4. the fit's counts
    points = {(p["curve"], float(p["stripe"])): p for p in rows(os.path.join(out, "fit", "fragility_points.csv"))}
    for state in STATES:
        for stripe, want in ((SPEEDS[0], 0), (SPEEDS[1], len(SEEDS))):
            p = points.get((state, float(stripe)))
            if p is None:
                bad.append("fit: no point for %s at %g m/s" % (state, stripe))
            elif int(p["members"]) != len(SEEDS) or float(p["exceedances"]) != want:
                bad.append("fit: %s at %g m/s has %s of %s members; expected %d of %d"
                           % (state, stripe, p["exceedances"], p["members"], want, len(SEEDS)))
    return bad


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, help="the ERF (flare) executable built with ERF_ENABLE_STRUCTURES")
    ap.add_argument("--out", required=True, help="directory for the plans, the members and the fit")
    ap.add_argument("--mpiexec", default="mpiexec", help="the MPI launcher; '' runs the executable directly")
    ap.add_argument("--ranks", type=int, default=1, help="MPI ranks of each run (default 1)")
    ap.add_argument("--uplift", type=float, default=20000.0,
                    help="the towers' synthetic footing uplift allowable, N (default 20000)")
    ap.add_argument("--load-start", type=float, default=6.0,
                    help="a run starts when the 1-minute load is below this (default 6, the driver's)")
    ap.add_argument("--no-load-gate", action="store_true", help="run without looking at the machine's load")
    ap.add_argument("--check", action="store_true", help="end with return code 1 unless the checks hold")
    a = ap.parse_args()
    if not a.uplift > 0.0:
        sys.exit("run_case.py: --uplift must be > 0 N, got %g" % a.uplift)
    out = os.path.abspath(a.out)
    for sub in ("ensemble", "fit"):
        if os.path.exists(os.path.join(out, sub)):
            # the driver skips members that finished, so an earlier run would be reused at another allowable
            sys.exit("run_case.py: %s already exists; remove it or give another --out" % os.path.join(out, sub))
    os.makedirs(out, exist_ok=True)
    write_json(os.path.join(out, "plan.json"), driver_plan(a.uplift))
    write_json(os.path.join(out, "fit_plan.json"), fit_plan())
    cmd = [sys.executable, "-I", "-u", DRIVER, "run", "--plan", os.path.join(out, "plan.json"), "--exe", a.exe,
           "--mpiexec", a.mpiexec, "--ranks", str(a.ranks), "--load-start", str(a.load_start)]
    if a.no_load_gate:
        cmd.append("--no-load-gate")
    rc = subprocess.call(cmd)
    if rc != 0:
        print("run_case.py: the ensemble driver ended with return code %d" % rc)
        return 1
    rc = subprocess.call([sys.executable, "-I", "-u", FIT, "--plan", os.path.join(out, "fit_plan.json")])
    if rc != 0:
        print("run_case.py: the fragility fit ended with return code %d" % rc)
        return 1
    with open(os.path.join(out, "fit", "fragility_points.csv"), encoding="utf-8") as f:
        print(f.read(), end="")
    if not a.check:
        return 0
    bad = check(out, a.uplift)
    for b in bad:
        print("FAIL " + b)
    print("CascadeStripes: %s" % ("FAILED" if bad else "all checks hold"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
