#!/usr/bin/env python3
"""Write the deck of a broken-wire run: a straight conductor section over flat ground in still air,
on insulator strings at standing towers, in which one component is given a capacity below its load
so that it lets go at a chosen time (conductors.failure_break, conductors.capacity_check_start).

    python3 make_case.py --out run_dir [options]

Writes, into the output directory:
  inputs        the ERF deck (flow at rest; the lines take still air through
                conductors.prescribed_velocity). With a release it ends one step after the
                release and writes a checkpoint there
  inputs_after  with a release: the deck that restarts from that checkpoint and runs to --time
                without the capacity that set the release, and with --tower_capacity and
                --check_interval if given
  plan.json     the same section as a plan of Structures/Tools/broken_wire/broken_wire.py, with
                what is released first, for the state at rest the run is compared with

The section runs along x. Its spans' unstretched lengths give every span the same slack ratio at
300 m scaled with the cube of the span (--slack), so that the spans pull alike along the line and
the strings hang plumb. The first release is one of:

  --release span_a   span 1 lets go of dead end a: the hardware's rated load is set between the
                     pulls on the two dead ends, dead end a standing --raise_a higher than b
  --release span_b   span 1 lets go of the string at tower 1: the conductor's rated strength is
                     set between span 1's largest tension, at its higher end at the tower, and the
                     other spans' (dead end a stands --lower_a lower than the strings' feet)
  --release none     nothing is released: a run for the towers' loads on the whole line

Why two decks: a line has one rated load per kind of component, and the load that passes along
the line just after a release is well above the load at rest. A capacity set a few percent under
the first component's load would therefore release its neighbours too. The second deck replaces
that capacity by rated loads of ten times the largest tension at rest for the other kinds of
component, which release nothing and make capacities.csv record their largest loads after the
release (the strings and, after span_b, the attachments; after span_a, the spans), and the legs'.
What happens after the first release is then set by --tower_capacity alone: with it
the towers' footings are checked against that bearing capacity (N per leg) and a tower over it
lets go of its lines (conductors.tower_failure = release); without it nothing else lets go and the
section comes to the state at rest of plan.json.
Everything is standard library; the plan's state at rest comes from broken_wire.py.
"""
import argparse
import importlib.util
import json
import math
import os
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "broken_wire", os.path.join(HERE, "..", "..", "Tools", "broken_wire", "broken_wire.py"))
bw = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bw)

# 795 kcmil ACSR (Drake), the conductor of the model's tests
DIAMETER, MASS, EA = 0.0281, 1.628, 3.0e7
AIR = 1.2
CELL = 50.0      # m, the flow's cells: the lines put nothing into the flow


def parse(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--spans", type=float, nargs="+", default=[300.0, 300.0, 300.0], help="span lengths (m), at least two")
    ap.add_argument("--height", type=float, default=60.0, help="height of the cross-arms (m)")
    ap.add_argument("--string", type=float, default=2.5, help="insulator string length (m)")
    ap.add_argument("--string_mass", type=float, default=60.0, help="mass of one string (kg)")
    ap.add_argument("--slack", type=float, default=1.5, help="slack of a 300 m span (m); other spans by the cube of their length")
    ap.add_argument("--segments", type=int, default=20, help="line segments per span")
    ap.add_argument("--release", choices=("span_a", "span_b", "none"), default="span_a")
    ap.add_argument("--raise_a", type=float, default=30.0, help="span_a: dead end a above dead end b (m)")
    ap.add_argument("--lower_a", type=float, default=30.0, help="span_b: dead end a below the strings' feet (m)")
    ap.add_argument("--release_time", type=float, default=5.0, help="the checks start then (s): the first release")
    ap.add_argument("--tower_capacity", type=float, default=0.0, help="after the first release: bearing capacity of a tower's leg (N); 0: not checked")
    ap.add_argument("--tower_weight", type=float, default=60000.0, help="weight of a tower (N)")
    ap.add_argument("--base_width", type=float, default=6.0, help="width of a tower's base (m)")
    ap.add_argument("--check_interval", type=float, default=0.0, help="after the first release: conductors.capacity_check_interval (s); 0: once per step")
    ap.add_argument("--drag_coefficient", type=float, default=1.0, help="the conductor's drag coefficient")
    ap.add_argument("--dt", type=float, default=0.25, help="time step (s)")
    ap.add_argument("--time", type=float, default=60.0, help="length of the run (s)")
    a = ap.parse_args(argv)
    if len(a.spans) < 2 or any(not (s > 0.0) for s in a.spans):
        ap.error("--spans: at least two spans, each > 0")
    if not (a.string > 0.0 and a.string_mass > 0.0):
        ap.error("--string and --string_mass must be > 0: the section hangs from strings")
    if not a.height > a.string + (a.lower_a if a.release == "span_b" else 0.0) + 1.0:
        ap.error("--height: the strings' feet and the dead ends must stand above the ground")
    if not (a.slack > 0.0 and a.segments >= 2 and a.dt > 0.0 and a.time > a.release_time >= 0.0):
        ap.error("--slack > 0, --segments >= 2, --dt > 0 and --time > --release_time >= 0 are needed")
    if a.tower_capacity < 0.0 or a.check_interval < 0.0 or a.drag_coefficient < 0.0:
        ap.error("--tower_capacity, --check_interval and --drag_coefficient must be >= 0")
    if a.release == "none" and (a.tower_capacity > 0.0 or a.check_interval > 0.0):
        ap.error("--tower_capacity and --check_interval act after a first release: they need --release span_a or span_b")
    if abs(a.release_time / a.dt - round(a.release_time / a.dt)) > 1.0e-9:
        ap.error("--release_time must be a whole number of steps --dt")
    return a


def section(a):
    """The plan of the section and its x positions."""
    foot = a.height - a.string
    za = zb = foot
    if a.release == "span_a":
        za = foot + a.raise_a
    elif a.release == "span_b":
        za = foot - a.lower_a
    x0 = 4.0 * CELL
    xs = [x0]
    for s in a.spans:
        xs.append(xs[-1] + s)
    lengths = []
    for k, s in enumerate(a.spans):
        chord = math.hypot(s, (foot - za) if k == 0 else 0.0)
        lengths.append(round(chord + a.slack * (s / 300.0) ** 3, 6))
    plan = {"end_a": [xs[0], za], "end_b": [xs[-1], zb], "towers": [[x, a.height] for x in xs[1:-1]],
            "lengths": lengths, "diameter": DIAMETER, "mass_per_length": MASS, "axial_stiffness": EA,
            "insulator_length": a.string, "insulator_mass": a.string_mass, "air_density": AIR,
            "segments": a.segments}
    return plan, xs


def between(low, high, what):
    """A capacity between two loads that are far enough apart to tell."""
    if not high > 1.01 * low:
        raise SystemExit("make_case: %s: the loads %.1f and %.1f N are less than 1 %% apart; "
                         "raise --raise_a or --lower_a" % (what, low, high))
    return 0.5 * (low + high)


def main(argv=None):
    a = parse(argv)
    plan, xs = section(a)
    whole = bw.solve_section(plan)
    line = []
    released = None
    if a.release == "span_a":
        pa = math.hypot(*whole["pull_on_end_a"])
        pb = math.hypot(*whole["pull_on_end_b"])
        strings = max(t["tension_top"] for t in whole["strings"])
        cap = between(max(pb, strings), pa, "the pulls on the dead ends")
        line.append("conductors.L.hardware_rated_load = %.1f   # between the pulls on dead end b (%.1f N) and a (%.1f N)" % (cap, pb, pa))
        released = {"span": 1, "end": "a"}
    elif a.release == "span_b":
        first = whole["spans"][0]
        if not first["tension_b"] > first["tension_a"]:
            raise SystemExit("make_case: span 1's tension is not highest at tower 1; raise --lower_a")
        others = max(max(s["tension_a"], s["tension_b"]) for s in whole["spans"][1:])
        cap = between(others, first["tension_b"], "the spans' largest tensions")
        line.append("conductors.L.rated_breaking_strength = %.1f   # between the other spans' largest tension (%.1f N) and span 1's at tower 1 (%.1f N)"
                    % (cap, others, first["tension_b"]))
        released = {"span": 1, "end": "b"}
    if released is not None:
        plan["released"] = released

    lx = xs[-1] + 4.0 * CELL
    nx = int(math.ceil(lx / CELL))
    lx = nx * CELL
    nz = int(math.ceil((a.height + max(a.raise_a, 0.0) + 40.0) / 25.0))
    steps = int(round(a.time / a.dt))
    # the first deck ends one step after the release and writes a checkpoint there
    first = int(round(a.release_time / a.dt)) + 1 if released is not None else steps
    record = 10.0 * max(max(s["tension_a"], s["tension_b"]) for s in whole["spans"])
    y = 4.0 * CELL
    pts = lambda p: "%.3f %.3f %.3f" % (p[0], y, p[1])

    def deck(after):
        d = []
        d.append("# Broken wire: %d spans (%s m) on %.1f m strings, cross-arms %.1f m up, flat ground, still air." % (
            len(a.spans), " ".join("%g" % s for s in a.spans), a.string, a.height))
        if after:
            d.append("# Written by Structures/Exec/BrokenWire/make_case.py; restarts the deck 'inputs' one step after its release")
            d.append("# (%s at %g s) with the capacities of the rest of the run in place of the one that set the release." % (a.release, a.release_time))
        else:
            d.append("# Written by Structures/Exec/BrokenWire/make_case.py; first release: %s at %g s." % (a.release, a.release_time))
        d.append("max_step = %d" % (steps if after else first))
        d.append("stop_time = %g" % (a.time + 1.0))
        d.append("erf.fixed_dt = %g" % a.dt)
        d.append("")
        d.append("geometry.prob_lo     = 0. 0. 0.")
        d.append("geometry.prob_hi     = %g %g %g" % (lx, 8.0 * CELL, nz * 25.0))
        d.append("amr.n_cell           = %d 8 %d" % (nx, nz))
        d.append("amr.max_level        = 0")
        d.append("geometry.is_periodic = 1 1 0")
        d.append("zlo.type = \"SlipWall\"")
        d.append("zhi.type = \"SlipWall\"")
        d.append("")
        d.append("# the air is at rest and the lines put nothing into it")
        d.append("erf.anelastic = 1")
        d.append("erf.vert_implicit = true")
        d.append("erf.anelastic_type = MidPoint")
        d.append("erf.use_fft   = true")
        d.append("erf.molec_diff_type = \"None\"")
        d.append("erf.les_type        = \"None\"")
        d.append("erf.init_type = \"uniform\"")
        d.append("erf.prob_name = \"ABL\"")
        d.append("prob.rho_0 = 1.0")
        d.append("prob.T_0   = 300.0")
        if not after:
            d.append("prob.U_0   = 0.0")
        d.append("erf.sum_interval = -1")
        d.append("erf.plot_int_1   = -1")
        d.append("erf.check_file   = chk")
        if after:
            d.append("erf.check_int    = -1")
            d.append("erf.restart      = chk%05d" % first)
        else:
            d.append("erf.check_int    = %d" % (first if released is not None else -1))
        d.append("")
        d.append("erf.structures.model = conductors")
        d.append("conductors.lines = L")
        d.append("conductors.L.end_a  = " + pts(plan["end_a"]))
        d.append("conductors.L.towers = " + "  ".join(pts(t) for t in plan["towers"]))
        d.append("conductors.L.end_b  = " + pts(plan["end_b"]))
        d.append("conductors.L.length = " + " ".join("%.6f" % l for l in plan["lengths"]))
        d.append("conductors.L.diameter         = %g" % DIAMETER)
        d.append("conductors.L.mass_per_length  = %g" % MASS)
        d.append("conductors.L.axial_stiffness  = %g" % EA)
        d.append("conductors.L.segments         = %d" % a.segments)
        d.append("conductors.L.drag_coefficient = %g" % a.drag_coefficient)
        d.append("conductors.L.insulator_length = %g" % a.string)
        d.append("conductors.L.insulator_mass   = %g" % a.string_mass)
        d.append("conductors.L.tower_type       = lattice")
        if not after:
            d.extend(line)
        else:
            # ten times the largest tension of the whole section: nothing of the line lets go again, and
            # capacities.csv records the largest load after the release of the components of the other kinds.
            # The kind that set the release is left out: its records come with the checkpoint, and a largest
            # ratio taken against the first deck's capacity would hide every later one against this one
            if a.release != "span_b":
                d.append("conductors.L.rated_breaking_strength = %.1f   # ten times the largest tension at rest: records the peaks, releases nothing" % record)
            d.append("conductors.L.insulator_rated_load    = %.1f   # ten times the largest tension at rest: records the peaks, releases nothing" % record)
            if a.release != "span_a":
                d.append("conductors.L.hardware_rated_load     = %.1f" % record)
        d.append("conductors.tower_types = lattice")
        d.append("conductors.lattice.base_width = %g" % a.base_width)
        d.append("conductors.lattice.top_width  = 1.5")
        d.append("conductors.lattice.solidity   = 0.2")
        d.append("conductors.lattice.arm_length = 12")
        d.append("conductors.lattice.arm_depth  = 1.2")
        d.append("conductors.lattice.weight     = %g" % a.tower_weight)
        if after and a.tower_capacity > 0.0:
            d.append("conductors.lattice.allowable_compression = %.1f" % a.tower_capacity)
            d.append("conductors.tower_failure = release")
        elif after:
            d.append("conductors.lattice.allowable_compression = 1.0e9   # records the legs' largest load, releases nothing")
        d.append("conductors.air_density         = %g" % AIR)
        d.append("conductors.prescribed_velocity = 0. 0. 0.")
        d.append("conductors.diagnostics_dir     = conductors")
        d.append("conductors.ground_contact      = true")
        if released is not None:
            d.append("conductors.failure_break       = true")
            d.append("conductors.capacity_check_start = %g" % a.release_time)
        if after and a.check_interval > 0.0:
            d.append("conductors.capacity_check_interval = %g" % a.check_interval)
        return "\n".join(d) + "\n"

    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "inputs"), "w") as f:
        f.write(deck(False))
    if released is not None:
        with open(os.path.join(a.out, "inputs_after"), "w") as f:
            f.write(deck(True))
    with open(os.path.join(a.out, "plan.json"), "w") as f:
        json.dump(plan, f, indent=1)
        f.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
