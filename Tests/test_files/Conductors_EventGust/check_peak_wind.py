#!/usr/bin/env python3
"""Check conductors/demand_table.csv and intensity.csv of the Conductors_EventGust deck run with
conductors.intensity_measure = wind: two lines on lattice towers over hills in a k-equation RANS wind, with a
travelling gust that makes the spans' tensions peak before the run ends.

  1. every component with a point keeps the wind of its peak's step: peak_wind_time is peak_time, and at least
     --mid components peak before the run's end with a wind that differs from their point's last value (a table
     that kept the last wind, or the first, fails);
  2. the points are one per span, one per attachment and one per tower; each attachment's point is sampled at
     the height asked for above the terrain there, as conductors/ground.dat gives it (the hills make the ground
     differ by tens of metres between points); a footing's point is its tower's, a conductor's its span's, the
     hardware's its attachment's;
  3. with the k-equation RANS every point carries a positive k, now, at its maximum and at the peaks;
  4. a point's maximum is a full window's average: not stamped before the averaging time, and not below the
     last value.

    python3 check_peak_wind.py [--height 10] [--averaging 3] [--mid 3]
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
    ap.add_argument("--mid", type=int, default=3)
    a = ap.parse_args()
    bad = []
    pts = {p["point"]: p for p in rows("conductors/intensity.csv")}
    table = rows("conductors/demand_table.csv")
    t_end = max(float(p["time"]) for p in pts.values())

    for n, p in pts.items():
        if float(p["height"]) != a.height:
            bad.append(f"{n}: height {p['height']}")
        for k in ("k", "max_k"):
            if p[k] == "" or not float(p[k]) > 0.0:
                bad.append(f"{n}: {k} is {p[k]!r}, not a positive k")
        if p["max_speed"] == "" or float(p["max_speed"]) < float(p["speed"]):
            bad.append(f"{n}: maximum {p['max_speed']!r} below the last value {p['speed']}")
        elif float(p["max_time"]) < a.averaging:
            bad.append(f"{n}: maximum at {p['max_time']} s, before the window of {a.averaging} s was full")
        if float(p["time"]) != t_end:
            bad.append(f"{n}: stamped {p['time']}, the run ends at {t_end}")

    # each attachment's point at the height above the terrain under it (ground.dat: line point x y ground z)
    grounds = []
    with open("conductors/ground.dat") as f:
        next(f)
        for ln in f:
            w = ln.split()
            if len(w) < 6 or w[1] == "transformer":
                continue
            which = {"a": "end_a", "b": "end_b"}.get(w[1], "tower" + w[1][1:])
            n = f"{w[0]}_{which}"
            if n not in pts:
                bad.append(f"ground.dat's attachment {w[0]} {w[1]} has no intensity point {n}")
                continue
            grounds.append(float(w[4]))
            if abs(float(pts[n]["z"]) - (float(w[4]) + a.height)) > 1e-6:
                bad.append(f"{n}: sampled at z = {pts[n]['z']}, the ground there is {w[4]} and the height {a.height}")
    if not grounds or max(grounds) - min(grounds) < 10.0:
        bad.append("the case must have ground that differs by 10 m or more between attachments")

    # one point per span and attachment of each line, one per tower
    lines = sorted({r["component"].rsplit("_span", 1)[0] for r in table if r["kind"] == "conductor"})
    spans = [n for n in pts if "_span" in n]
    ends = [n for n in pts if n.endswith("_end_a") or n.endswith("_end_b")]
    attach = [n for n in pts if "_tower" in n and not n.startswith("tower_")]
    towers = [n for n in pts if n.startswith("tower_")]
    if len(ends) != 2 * len(lines) or len(attach) != len(spans) - len(lines) or not towers:
        bad.append(f"points: {len(spans)} span(s), {len(ends)} end(s), {len(attach)} tower attachment(s), "
                   f"{len(towers)} tower(s) for {len(lines)} line(s)")

    mid = 0
    kinds = set()
    for r in table:
        n, kind = r["component"], r["kind"]
        kinds.add(kind)
        if kind == "conductor" or kind == "hardware":
            want = n
        elif kind.startswith("footing"):
            want = "tower_" + n.rsplit("_footing", 1)[0]
        else:
            continue
        if r["point"] != want or want not in pts:
            bad.append(f"{n}: its intensity point is {r['point']!r}, expected {want!r}")
            continue
        if r["peak_time"] == "":
            bad.append(f"{n}: never recorded")
            continue
        if r["peak_wind_time"] != r["peak_time"]:
            bad.append(f"{n}: the wind kept is stamped {r['peak_wind_time']}, its peak is at {r['peak_time']}")
        if r["peak_k"] == "" or not float(r["peak_k"]) > 0.0:
            bad.append(f"{n}: peak_k is {r['peak_k']!r}")
        if float(r["peak_time"]) < t_end and r["peak_speed"] != pts[want]["speed"]:
            mid += 1
    if mid < a.mid:
        bad.append(f"only {mid} component(s) peak before the end with a wind unlike the last; the case needs {a.mid}")
    for kind in ("conductor", "hardware", "footing_uplift", "footing_compression"):
        if kind not in kinds:
            bad.append(f"the case holds no {kind}")

    print(f"{len(pts)} point(s), {len(table)} component(s), {mid} peaking before the end at {t_end} s with their own wind")
    for b in bad:
        print("FAIL:", b)
    print("PASS" if not bad else "FAILED")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
