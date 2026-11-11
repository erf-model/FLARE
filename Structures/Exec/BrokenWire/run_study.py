#!/usr/bin/env python3
"""Run the broken-wire studies: decks from make_case.py through ERF, compared with the state at
rest of Structures/Tools/broken_wire/broken_wire.py.

    python3 run_study.py --exe path/to/erf_exec --out study_dir [--launcher "mpiexec -n 2"]
                         [--study residual overshoot cascade] [--check PERCENT]

Needs an ERF built with ERF_ENABLE_STRUCTURES and the real MoorDyn-C (the stub has no dynamics). The
runs are in still air on a coarse mesh and take seconds each; they run one after another.

  residual   one span cut on four sections; the towers' pulls, averaged over the last 240 s of
             600 s, against the force balance (residual.csv)
  overshoot  span 1 cut at tower 1, the middle span 0.5 to 2 times the others, cross-arms 100 m up
             and dead end a 90 m lower so that span 1 carries the largest tension: the largest load
             after the release over the load at rest and over the load before (overshoot.csv)
  cascade    eight spans of 300 m, span 1 cut at tower 1, the legs' bearing capacity from 1.5 to
             24 times their load before the release, checked once per step and every 0.01 s: the
             number of towers that fail (cascade.csv)

With --check the residual study ends with return code 1 if a tower's pull differs from the force
balance by more than that many percent of tower 1's largest pull (the CTest BrokenWire_ThroughERF
uses 0.2).

Every run is one realisation: a fall amplifies roundoff, so the largest loads and a cascade's
extent at a capacity close to a largest load differ between builds and rank counts.
Standard library only.
"""
import argparse
import csv
import importlib.util
import json
import math
import os
import shlex
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "broken_wire", os.path.join(HERE, "..", "..", "Tools", "broken_wire", "broken_wire.py"))
bw = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bw)


def run_case(a, name, options):
    """Write the decks of one case and run them; returns its directory."""
    d = os.path.join(a.out, name)
    if os.path.isdir(d):
        shutil.rmtree(d)
    subprocess.run([sys.executable, os.path.join(HERE, "make_case.py"), "--out", d] + [str(o) for o in options], check=True)
    for deck in ("inputs", "inputs_after"):
        if not os.path.isfile(os.path.join(d, deck)):
            continue
        with open(os.path.join(d, deck + ".log"), "w") as log:
            rc = subprocess.run(a.launch + [a.exe, deck], cwd=d, stdout=log, stderr=subprocess.STDOUT).returncode
        if rc != 0:
            raise SystemExit("run_study: %s: %s ended with return code %d (see %s.log)" % (name, deck, rc, os.path.join(d, deck)))
        if deck == "inputs" and os.path.isfile(os.path.join(d, "inputs_after")) and len(rows(os.path.join(d, "conductors", "breaks.dat"))) != 1:
            raise SystemExit("run_study: %s: the first deck must end with one release (see conductors/breaks.dat)" % name)
    return d


def rows(path):
    """A space-separated table with one header line, as a list of dictionaries of strings."""
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        head = f.readline().split()
        return [dict(zip(head, line.split())) for line in f if line.strip()]


def series(d):
    """towers.dat: the header's names and the rows as floats."""
    out = []
    for r in rows(os.path.join(d, "conductors", "towers.dat")):
        out.append({k: float(v) for k, v in r.items()})
    return out


def mean(table, name, t_from):
    v = [r[name] for r in table if r["time"] >= t_from]
    return sum(v) / len(v)


def records(d):
    """capacities.csv by component."""
    with open(os.path.join(d, "conductors", "capacities.csv")) as f:
        return {r["component"]: r for r in csv.DictReader(f)}


def plan_of(d, released=True):
    with open(os.path.join(d, "plan.json")) as f:
        p = json.load(f)
    if not released:
        p.pop("released", None)
    return p


def write(a, name, head, table):
    path = os.path.join(a.out, name)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(head)
        w.writerows(table)
    widths = [max(len(str(x)) for x in col) for col in zip(head, *table)]
    print("\n" + name)
    for r in [head] + table:
        print("  ".join(str(x).ljust(n) for x, n in zip(r, widths)))


def residual(a):
    cases = [("equal_2.5m", ["--spans", 300, 300, 300, "--release", "span_b"]),
             ("equal_4m", ["--spans", 300, 300, 300, "--string", 4.0, "--release", "span_b"]),
             ("unequal", ["--spans", 300, 500, 100, "--release", "span_b", "--height", 100, "--lower_a", 90]),
             ("dead_end", ["--spans", 300, 300, 300, "--release", "span_a"])]
    out = []
    for name, opt in cases:
        d = run_case(a, "residual_" + name, opt + ["--time", 600, "--release_time", 10])
        rest = bw.solve_section(plan_of(d))
        whole = bw.solve_section(plan_of(d, False))
        t = series(d)
        for i, s in enumerate(rest["strings"]):
            for q, k in (("along the line", 0), ("vertical", 1)):
                model = mean(t, "L_t%d_line_F%s" % (i + 1, "xz"[k]), 360.0)
                ref = s["pull_on_tower"][k]
                out.append([name, "tower %d" % (i + 1), q, "%.1f" % model, "%.1f" % ref,
                            "%.2f" % (100.0 * (model - ref) / max(abs(c) for c in rest["strings"][0]["pull_on_tower"]))])
        h_before = whole["spans"][1]["horizontal_tension"]
        h_after = mean(t, "L_t1_line_Fx", 360.0) if name != "dead_end" else rest["spans"][1]["horizontal_tension"]
        out.append([name, "span 2", "tension kept (%)", "%.1f" % (100.0 * h_after / h_before),
                    "%.1f" % (100.0 * rest["spans"][1]["horizontal_tension"] / h_before), ""])
    write(a, "residual.csv", ["case", "where", "quantity", "model (N)", "force balance (N)", "difference (% of tower 1's largest pull)"], out)
    worst = [r for r in out if r[5] and abs(float(r[5])) > a.check] if a.check > 0.0 else []
    for r in worst:
        print("run_study: %s, %s, %s: the model's %s N differs from the force balance's %s N by %s %%, more than %g %%"
              % (r[0], r[1], r[2], r[3], r[4], r[5], a.check))
    return 1 if worst else 0


def overshoot(a):
    out = []
    for ratio in (0.5, 1.0, 1.5, 2.0):
        d = run_case(a, "overshoot_%g" % ratio, ["--spans", 300, 300 * ratio, 300, "--release", "span_b", "--height", 100, "--lower_a", 90,
                                                 "--time", 300, "--release_time", 10, "--check_interval", 0.01])
        rest = bw.solve_section(plan_of(d))
        whole = bw.solve_section(plan_of(d, False))
        rec = records(d)
        t = series(d)
        loads = [("string 1", float(rec["L_string1"]["peak_demand"]), rest["strings"][0]["tension_top"], whole["strings"][0]["tension_top"]),
                 ("string 2", float(rec["L_string2"]["peak_demand"]), rest["strings"][1]["tension_top"], whole["strings"][1]["tension_top"]),
                 ("dead end b", float(rec["L_end_b"]["peak_demand"]), math.hypot(*rest["pull_on_end_b"]), math.hypot(*whole["pull_on_end_b"])),
                 ("tower 1 leg", float(rec["L_t1_footing_compression"]["peak_demand"]), mean(t, "L_t1_max_compression", 180.0), t[0]["L_t1_max_compression"]),
                 ("tower 2 leg", float(rec["L_t2_footing_compression"]["peak_demand"]), mean(t, "L_t2_max_compression", 180.0), t[0]["L_t2_max_compression"])]
        for what, peak, at_rest, before in loads:
            out.append(["%g" % ratio, what, "%.0f" % before, "%.0f" % at_rest, "%.0f" % peak, "%.2f" % (peak / at_rest), "%.2f" % (peak / before)])
    write(a, "overshoot.csv", ["span 2 / span 1", "load", "before (N)", "at rest after (N)", "largest (N)", "largest / at rest", "largest / before"], out)


def cascade(a):
    spans = ["--spans"] + [300] * 8
    d = run_case(a, "cascade_whole", spans + ["--release", "none", "--time", 1, "--release_time", 0])
    before = series(d)[0]["L_t4_max_compression"]
    out = []
    for factor in (1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 16.0, 24.0):
        row = ["%g" % factor, "%.0f" % (factor * before)]
        for interval in (0.0, 0.01):
            d = run_case(a, "cascade_%g_%g" % (factor, interval), spans + ["--release", "span_b", "--time", 120, "--release_time", 10,
                                                                         "--tower_capacity", factor * before, "--check_interval", interval])
            failed = rows(os.path.join(d, "conductors", "tower_failures.dat"))
            row += [len(failed), failed[-1]["time"] if failed else ""]
        out.append(row)
    write(a, "cascade.csv", ["capacity / load before", "capacity (N per leg)", "towers failed of 7, one check per step", "last at (s)",
                             "towers failed of 7, checks every 0.01 s", "last at (s)"], out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, help="the ERF executable (real MoorDyn-C)")
    ap.add_argument("--out", required=True, help="directory for the runs and the tables")
    ap.add_argument("--launcher", default="mpiexec -n 2", help="what starts the executable, with its rank count; \"\": nothing")
    ap.add_argument("--check", type=float, default=0.0, help="residual: largest difference allowed (%% of tower 1's largest pull); 0: no check")
    ap.add_argument("--study", nargs="+", choices=("residual", "overshoot", "cascade"), default=["residual", "overshoot", "cascade"])
    a = ap.parse_args(argv)
    a.exe = os.path.abspath(a.exe)
    if not os.path.isfile(a.exe):
        ap.error("--exe: no such file: " + a.exe)
    if not (a.check >= 0.0 and math.isfinite(a.check)):
        ap.error("--check must be finite and >= 0")
    a.launch = shlex.split(a.launcher)
    os.makedirs(a.out, exist_ok=True)
    rc = 0
    for s in a.study:
        rc = max(rc, {"residual": residual, "overshoot": overshoot, "cascade": cascade}[s](a) or 0)
    return rc


if __name__ == "__main__":
    sys.exit(main())
