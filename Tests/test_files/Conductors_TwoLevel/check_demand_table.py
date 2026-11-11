#!/usr/bin/env python3
"""Check conductors/demand_table.csv and intensity.csv of the Conductors_TwoLevel deck run with
conductors.intensity_measure = wind and capacities on its span and dead ends.

The deck's wind is the crosswind v = 15 + 0.1 z m/s along +y (u = 0), so the wind at the height h above the
flat ground is 15 + 0.1 h, direction 90 degrees, whatever the averaging time. The checks:

  1. every intensity point (the span's middle and its two ends) reads that speed within --rtol and that
     direction within --dtol, at the given height above the flat ground (z = height), stamped with the run's
     last time; its run maximum reads the same once the run is as long as the averaging time, and is empty
     before (a maximum is a full window's average);
  2. the demand table's peak ratio, demand, capacity and time, and its first exceedance, are capacities.csv's
     for every component, and the conductor's peak demand is the largest tension in the span's step log S1.dat;
  3. each component's wind at its peak and at its first exceedance is the closed form too, the first
     exceedance's stamped fields are filled only for a component that exceeded, and the point's maximum is not
     below the wind at the peak.

    python3 check_demand_table.py [--height 10] [--averaging 3] [--rtol 2e-3] [--dtol 0.05]
"""
import argparse
import csv
import sys


def rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--height", type=float, default=10.0)
    ap.add_argument("--averaging", type=float, default=3.0)
    ap.add_argument("--rtol", type=float, default=2e-3)
    ap.add_argument("--dtol", type=float, default=0.05)
    a = ap.parse_args()
    want = 15.0 + 0.1 * a.height
    bad = []

    def near(label, value, target, tol):
        if value == "" or abs(float(value) - target) > tol:
            bad.append(f"{label}: {value!r}, expected {target} +- {tol}")

    # the run's last time and the span's largest tension, from its step log
    with open("S1.dat") as f:
        head = f.readline().split()
        log = [ln.split() for ln in f if ln.strip()]
    t_end = float(log[-1][head.index("time")])
    peak_tension = max(float(r[head.index("max_tension")]) for r in log)

    pts = rows("conductors/intensity.csv")
    if sorted(p["point"] for p in pts) != ["S1", "S1_end_a", "S1_end_b"]:
        bad.append(f"intensity.csv points: {[p['point'] for p in pts]}")
    for p in pts:
        n = p["point"]
        near(f"{n} height", p["height"], a.height, 1e-9)
        near(f"{n} z", p["z"], a.height, 1e-9)
        near(f"{n} speed", p["speed"], want, a.rtol * want)
        near(f"{n} direction", p["direction"], 90.0, a.dtol)
        near(f"{n} time", p["time"], t_end, 1e-9)
        if t_end >= a.averaging:
            near(f"{n} max_speed", p["max_speed"], want, a.rtol * want)
            near(f"{n} max_direction", p["max_direction"], 90.0, a.dtol)
            if not float(p["max_time"]) >= a.averaging:
                bad.append(f"{n}: maximum at {p['max_time']} s, before the window of {a.averaging} s was full")
        elif p["max_speed"] or p["max_time"] or p["max_direction"]:
            bad.append(f"{n}: a run maximum after {t_end} s, less than the averaging time {a.averaging} s")
        if p["k"] != "" or p["max_k"] != "":
            bad.append(f"{n}: k given without the k-equation RANS")
        if p["max_speed"] and p["speed"] and float(p["max_speed"]) < float(p["speed"]):
            bad.append(f"{n}: the run's maximum {p['max_speed']} is below the last value {p['speed']}")
    by_point = {p["point"]: p for p in pts}

    cap = {r["component"]: r for r in rows("conductors/capacities.csv")}
    table = rows("conductors/demand_table.csv")
    if [r["component"] for r in table] != list(cap):
        bad.append("demand_table.csv does not list capacities.csv's components in order")
    for r in table:
        n = r["component"]
        c = cap.get(n)
        if c is None:
            continue
        for k in ("kind", "peak_ratio", "peak_demand", "peak_capacity", "peak_time", "exceeded", "first_time",
                  "first_demand", "first_capacity"):
            if r[k] != c[k]:
                bad.append(f"{n}: {k} is {r[k]!r} in the demand table and {c[k]!r} in capacities.csv")
        if r["point"] != n:
            bad.append(f"{n}: its intensity point is {r['point']!r}")
            continue
        near(f"{n} peak_speed", r["peak_speed"], want, a.rtol * want)
        near(f"{n} peak_direction", r["peak_direction"], 90.0, a.dtol)
        p = by_point[r["point"]]
        for k in ("max_speed", "max_direction", "max_time"):
            if r[k] != p[k]:
                bad.append(f"{n}: {k} is {r[k]!r} in the demand table and {p[k]!r} at its point")
        if r["exceeded"] == "1":
            near(f"{n} first_speed", r["first_speed"], want, a.rtol * want)
            near(f"{n} first_direction", r["first_direction"], 90.0, a.dtol)
        elif r["first_speed"] or r["first_direction"] or r["first_time"]:
            bad.append(f"{n}: first-exceedance fields without an exceedance")
    conductor = [r for r in table if r["kind"] == "conductor"]
    if len(conductor) != 1 or float(conductor[0]["peak_demand"]) != peak_tension:
        bad.append(f"the conductor's peak demand {[r['peak_demand'] for r in conductor]} is not S1.dat's largest "
                   f"tension {peak_tension!r}")
    if not any(r["exceeded"] == "1" for r in table) or not any(r["kind"] == "hardware" for r in table):
        bad.append("the case must hold an exceedance and the dead ends' hardware")

    print(f"wind at {a.height} m: closed form {want} m/s at 90 deg; points "
          + ", ".join(f"{p['point']} {p['speed']} m/s" for p in pts)
          + f"; conductor peak demand {conductor[0]['peak_demand'] if conductor else '?'} N, "
            f"largest tension in S1.dat {peak_tension} N")
    for b in bad:
        print("FAIL:", b)
    print("PASS" if not bad else "FAILED")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
