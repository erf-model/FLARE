#!/usr/bin/env python3
"""Seven wind stripes over the PowerLines hills through the fragility ensemble driver, with the line
failures on, and the damage states DS1-DS4 from the fragility fit.

    python3 run_case.py --exe path/to/flare --out case_dir [--speeds 11 12 13 14 15 16 17] [--seeds 1 2 3]
                        [--mpiexec mpiexec] [--mpiexec-args=...] [--ranks 4] [--load-start 10 | --no-load-gate]
                        [--dry-run]

The deck (powerlines_stripes.i) is Structures/Exec/PowerLines in one stage: five lines on lattice
suspension towers between six transformers over four hills, k-equation RANS, random gusts. The
script copies the deck and the case's terrain, inflow profile and sounding into --out/base, writes the
driver's plan (plan.json) and the fit's plan (fit_plan.json) and runs both.

The driver scales the inflow to each stripe's speed at 10 m. Per stripe one spin-up runs the flow alone
for 1715 steps (600 s); each seed restarts from it with the lines and runs to step 2058 (120 s more),
with conductors.failure_break, ground_contact and tower_failure = release. The capacities checked from
610 s on are the footings' SYNTHETIC allowables of the deck (uplift 50 kN, compression 150 kN) and the
insulators' specified mechanical load, 111,206 N (25,000 lb, Hubbell Hi*Lite catalog CA08051E); every
conductor's rated breaking strength is 1e9 N, so the conductors are not checked. The fit counts each
damage state per stripe (DS4: 2 or more towers lost).

A spin-up takes about 25 min and a member about 5 min on 4 ranks (MoorDyn-C), about 3.5 h for the
default 21 members. Standard library only.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CASE_FILES = os.path.join(ROOT, "Structures", "Exec", "PowerLines")
SUPPORT = ("terrain_hills.txt", "inflow_profile", "input_sounding")
DECK = "powerlines_stripes.i"
DRIVER = os.path.join(ROOT, "Structures", "Tools", "fragility_ensemble", "fragility_ensemble.py")
FIT = os.path.join(ROOT, "Structures", "Tools", "fragility_fit", "fragility_fit.py")

LINES = ("L1", "L2", "L3", "L4", "L5")
INSULATOR_RATED_LOAD = 111206.0   # N, 25,000 lb specified mechanical load, Hubbell Hi*Lite catalog CA08051E
SPINUP_STEPS = 1715               # 600 s at erf.fixed_dt = 0.35
CASCADE_TOWERS = 2
STATES = ("DS1", "DS2", "DS3", "DS4")


def driver_plan(speeds, seeds):
    inputs = {
        "conductors.gust_type": "random",
        "conductors.failure_break": True,
        "conductors.ground_contact": True,
        "conductors.tower_failure": "release",
        "conductors.capacity_check_start": 610.0,
        "conductors.capacity_check_interval": 0.01,
    }
    for ln in LINES:
        inputs["conductors.%s.insulator_rated_load" % ln] = INSULATOR_RATED_LOAD
    for ln in LINES:
        inputs["conductors.%s.rated_breaking_strength" % ln] = 1.0e9
    return {
        "base_deck": os.path.join("base", DECK),
        "output_dir": "ensemble",
        "speeds": list(speeds),
        "seeds": list(seeds),
        "reference_height": 10,
        "support_files": ["terrain_hills.txt"],
        "spinup": {"max_step": SPINUP_STEPS},
        "inputs": inputs,
    }


def fit_plan():
    return {
        "demand_table": "ensemble/demand_table_all.csv",
        "output_dir": "fit",
        "damage_states": {"table": "ensemble/damage_table_all.csv", "states": list(STATES),
                          "cascade_towers": CASCADE_TOWERS},
    }


def write_json(path, obj):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)
        f.write("\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", help="the ERF (flare) executable built with ERF_ENABLE_STRUCTURES and MoorDyn-C")
    ap.add_argument("--out", required=True, help="directory for the deck, the plans, the members and the fit")
    ap.add_argument("--speeds", type=float, nargs="+", default=[11, 12, 13, 14, 15, 16, 17],
                    help="the stripes, m/s at 10 m (default 11 to 17)")
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3], help="gust seeds per stripe (default 1 2 3)")
    ap.add_argument("--mpiexec", default="mpiexec", help="the MPI launcher; '' runs the executable directly")
    ap.add_argument("--mpiexec-args", default="", help="further launcher arguments in one string, passed to the "
                                                       "driver: --mpiexec-args='--oversubscribe'")
    ap.add_argument("--ranks", type=int, default=4, help="MPI ranks of each run (default 4)")
    ap.add_argument("--load-start", type=float, default=10.0,
                    help="a run starts when the 1-minute load is below this (default 10)")
    ap.add_argument("--no-load-gate", action="store_true", help="run without looking at the machine's load")
    ap.add_argument("--dry-run", action="store_true", help="write the decks, plans and manifest, run nothing")
    a = ap.parse_args()
    if not a.dry_run and not a.exe:
        sys.exit("run_case.py: --exe is needed unless --dry-run")
    if any(not v > 0.0 for v in a.speeds):
        sys.exit("run_case.py: --speeds must all be > 0 m/s, got %s" % a.speeds)
    out = os.path.abspath(a.out)
    os.makedirs(os.path.join(out, "base"), exist_ok=True)
    shutil.copyfile(os.path.join(HERE, DECK), os.path.join(out, "base", DECK))
    for f in SUPPORT:
        shutil.copyfile(os.path.join(CASE_FILES, f), os.path.join(out, "base", f))
    write_json(os.path.join(out, "plan.json"), driver_plan(a.speeds, a.seeds))
    write_json(os.path.join(out, "fit_plan.json"), fit_plan())
    os.chdir(out)
    cmd = [sys.executable, "-I", "-u", DRIVER, "run", "--plan", "plan.json",
           "--mpiexec", a.mpiexec, "--mpiexec-args=" + a.mpiexec_args, "--ranks", str(a.ranks),
           "--load-start", str(a.load_start)]
    if a.exe:
        cmd += ["--exe", os.path.abspath(a.exe)]
    if a.no_load_gate:
        cmd.append("--no-load-gate")
    if a.dry_run:
        cmd.append("--dry-run")
    rc = subprocess.call(cmd)
    if rc != 0 or a.dry_run:
        if rc != 0:
            print("run_case.py: the ensemble driver ended with return code %d" % rc)
        return 1 if rc else 0
    rc = subprocess.call([sys.executable, "-I", "-u", FIT, "--plan", "fit_plan.json"])
    if rc != 0:
        print("run_case.py: the fragility fit ended with return code %d" % rc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
