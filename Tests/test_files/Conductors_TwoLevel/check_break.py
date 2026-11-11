#!/usr/bin/env python3
"""Check a run of the Conductors_TwoLevel deck with conductors.failure_break = true, conductors.ground_contact = true
and a rated_breaking_strength the span's tension exceeds.

The deck is one level span, 301.5 m of conductor between dead ends 300 m apart and 30 m above flat ground, in the
crosswind v = 15 + 0.1 z m/s, the line on a refined level. The checks:

  1. breaks.dat holds one row, at --time, the end of the anchor-level step at which failures.dat logs the
     conductor's first exceedance: the span S1 lets go at one of its ends (--end pins which; a level span's two end
     tensions are a tie, which releases end a), and the tension logged there is over the capacity and is S1.dat's at
     that end at that time;
  2. from the next row of S1.dat on, the tension at the released end stays below --free-fraction of the capacity
     (a free end carries its last segment's inertia and drag only);
  3. the span reaches the ground after the release (ground_contact.csv's first touch is later than --time, and the
     last row's smallest clearance is within the conductor's diameter of the terrain);
  4. at the end of the run the tension at the end still held is below --fallen-fraction of the capacity: the span
     has fallen and no longer pulls between its supports. With --rest-rtol it must also be the weight of the
     conductor that hangs from it, w h, with h the support's 30 m less the conductor's radius and w its weight per
     metre less buoyancy (closed form), within that tolerance. That is for the stub, which places the rest state.
     MoorDyn-C's span still swings then, and its 15 m segments lump the hanging length: its rest tension is the
     weight of the nodes clear of the ground, up to one segment's weight below w h.

    python3 check_break.py --time 0.125 --capacity 5000 [--end a] [--rest-rtol 1e-6] [--fallen-fraction 0.5]
                           [--free-fraction 0.2] [--no-rest]
"""
import argparse
import csv
import math
import sys


def table(path):
    with open(path) as f:
        head = f.readline().split()
        return head, [ln.split() for ln in f if ln.strip()]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--time", type=float, required=True)
    ap.add_argument("--capacity", type=float, required=True)
    ap.add_argument("--end", choices=["a", "b"], default=None)
    ap.add_argument("--rest-rtol", type=float, default=None)
    ap.add_argument("--fallen-fraction", type=float, default=0.5)
    ap.add_argument("--free-fraction", type=float, default=0.2)
    ap.add_argument("--no-rest", action="store_true", help="the run ends before the span is at rest: skip checks 3 and 4")
    a = ap.parse_args()
    bad = []

    # the deck's conductor and air
    diameter, mass, rho, g, height = 0.0281, 1.628, 1.0, 9.81, 30.0
    w = (mass - rho * 0.25 * math.pi * diameter ** 2) * g
    rest = w * (height - 0.5 * diameter)

    bh, breaks = table("conductors/breaks.dat")
    if bh != ["time", "line", "released", "component", "tension"]:
        bad.append(f"breaks.dat header: {bh}")
    fh, fails = table("conductors/failures.dat")
    sh, span = table("S1.dat")
    col = {n: i for i, n in enumerate(sh)}
    end = None
    if len(breaks) != 1:
        bad.append(f"breaks.dat holds {len(breaks)} rows, one expected")
    else:
        t, line, released, component, tension = breaks[0]
        if float(t) != a.time:
            bad.append(f"the release is at {t} s, expected {a.time}")
        if line != "S1" or component != "S1" or released not in ("span1_end_a", "span1_end_b"):
            bad.append(f"the release row: {breaks[0]}")
        else:
            end = released[-1]
        if a.end and end and end != a.end:
            bad.append(f"released at end {end}, expected end {a.end}")
        if not float(tension) > a.capacity:
            bad.append(f"the tension at the release, {tension} N, is not over the capacity {a.capacity} N")
        conductor = [r for r in fails if r[fh.index("kind")] == "conductor"] if "kind" in fh else []
        if len(conductor) != 1 or float(conductor[0][0]) != float(t):
            bad.append(f"failures.dat's conductor rows {conductor} do not hold one exceedance at the release's time {t}")
        at = [r for r in span if float(r[col["time"]]) == float(t)]
        if end and (len(at) != 1 or abs(float(at[0][col["tension_" + end]]) - float(tension)) > 1e-9 * float(tension)):
            bad.append(f"the tension at the release, {tension} N, is not S1.dat's tension_{end} at {t} s")
    worst = 0.0
    if end:
        after = [r for r in span if float(r[col["time"]]) > a.time]
        if not after:
            bad.append("S1.dat holds no row after the release")
        worst = max((float(r[col["tension_" + end]]) for r in after), default=0.0)
        if worst > a.free_fraction * a.capacity:
            bad.append(f"the released end carries up to {worst} N after the release, more than "
                       f"{a.free_fraction} of the capacity")
    held_tension = None
    if end and not a.no_rest:
        with open("conductors/ground_contact.csv", newline="") as f:
            ground = list(csv.DictReader(f))
        touch = ground[0]["first_touch_time"] if ground else ""
        if touch == "" or not float(touch) > a.time:
            bad.append(f"the span first touched the ground at {touch!r} s, not after the release at {a.time} s")
        last = span[-1]
        if not float(last[col["min_clearance"]]) <= diameter:
            bad.append(f"at the end the span's smallest clearance is {last[col['min_clearance']]} m: it is not on the ground")
        held = "b" if end == "a" else "a"
        held_tension = float(last[col["tension_" + held]])
        if not held_tension < a.fallen_fraction * a.capacity:
            bad.append(f"at the end the held end's tension is {held_tension} N, not below {a.fallen_fraction} of the capacity")
        if a.rest_rtol is not None and abs(held_tension - rest) > a.rest_rtol * rest:
            bad.append(f"at the end the held end's tension is {held_tension} N, expected the hanging weight {rest:.4f} N "
                       f"within {a.rest_rtol}")

    print(f"release: {breaks[0] if breaks else None}; the released end's largest tension afterwards {worst:.4g} N; "
          f"the held end's tension at the end {held_tension} N, closed form w h = {rest:.4f} N")
    for b in bad:
        print("FAIL:", b)
    print("PASS" if not bad else "FAILED")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
