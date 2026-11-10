#!/usr/bin/env python3
"""CTest check of the benchmark's switch: erf.fire.fire_atm_feedback = 0 must keep the fire's heat out of the air.

Reads conductors/conductor_temperature.dat of a short run with the fire lit under the line (feedback off)
and of the same run with feedback on, and checks:

  on   the air temperature at some span rises by more than --rise_on K;
  off  the air temperature every span sees moves by less than --frac_off of that rise, while the fire's
       radiation reaches the line (E_fire_max > 0), so the conductors are heated by radiation alone.

The logged air temperature (T_air_C) is the air at the span's hottest node. As the fire heats one end of a span
more than its middle, that node moves along the span and up to 10 m in height, so with feedback off the logged
air still moves by a few hundredths of a kelvin in the neutral layer: the off limit is a fraction of the on rise,
not zero.

    python3 check_feedback.py off_dir on_dir [--t_ign 5] [--rise_on 1] [--frac_off 0.01]
"""
import argparse
import os
import sys


def spans(d):
    path = os.path.join(d, "conductors", "conductor_temperature.dat")
    if not os.path.exists(path):
        sys.exit(f"{path}: not found (did the run fail? see run.log there)")
    with open(path) as f:
        head = f.readline().lstrip("#").split()
        rows = [ln.split() for ln in f if ln.strip() and not ln.startswith("#")]
    col = {n: i for i, n in enumerate(head)}
    out = {}
    for r in rows:
        k = int(r[col["span"]])
        out.setdefault(k, []).append([float(r[col[c]]) for c in ("time", "T_air_C", "E_fire_max_W_m2")])
    if not out:
        sys.exit(f"{path}: no rows")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("off_dir")
    ap.add_argument("on_dir")
    ap.add_argument("--t_ign", type=float, default=5.0)
    ap.add_argument("--rise_on", type=float, default=1.0)
    ap.add_argument("--frac_off", type=float, default=0.01)
    a = ap.parse_args()
    ok = True

    def rises(runs, which):
        out = {}
        for k, v in sorted(runs.items()):
            before = [r for r in v if r[0] < a.t_ign]
            if not before:
                sys.exit(f"{which}: span {k} has no rows before ignition at {a.t_ign} s")
            out[k] = ([r[1] - before[-1][1] for r in v], max(r[2] for r in v))
        return out

    on = rises(spans(a.on_dir), "on")
    rise = max(max(d) for d, _ in on.values())
    for k, (d, _) in on.items():
        print(f"on: span {k}: air temperature rises by {max(d):.3f} K")
    if rise <= a.rise_on:
        ok = False
        print(f"FAIL: with feedback on the air at the line warmed by only {rise:.3f} K (needs > {a.rise_on} K)")

    off = rises(spans(a.off_dir), "off")
    lim = a.frac_off * max(rise, 0.0)
    for k, (d, e) in off.items():
        dev = max(abs(x) for x in d)
        print(f"off: span {k}: air temperature moves by at most {dev:.4f} K (limit {lim:.4f} K), "
              f"largest fire irradiance {e:.0f} W/m2")
        if dev > lim:
            ok = False
            print(f"FAIL: with feedback off the air at span {k} moved by {dev:.4f} K")
    if max(e for _, e in off.values()) <= 0.0:
        ok = False
        print("FAIL: the fire's radiation never reached the line with feedback off")

    print("PASS" if ok else "FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
