#!/usr/bin/env python3
"""Compare two fire plotfiles field by field, bit for bit.

    python3 compare_restart.py straight_00020 plt_fire_00020

Fields: fire_wind_ref_u, fire_wind_ref_v, fire_extract_z, fire_phi. Exits 1
if any differs. Standard library only (plus erf_plotfile.py).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import erf_plotfile  # noqa: E402

FIELDS = ["fire_wind_ref_u", "fire_wind_ref_v", "fire_extract_z", "fire_phi"]


def main():
    a_pf, b_pf = sys.argv[1], sys.argv[2]
    _, a = erf_plotfile.read_fields(a_pf, FIELDS)
    _, b = erf_plotfile.read_fields(b_pf, FIELDS)
    ok = True
    for f in FIELDS:
        diff = max(abs(x - y) for ra, rb in zip(a[f], b[f])
                   for ca, cb in zip(ra, rb) for x, y in zip(ca, cb))
        good = diff == 0.0
        ok = ok and good
        print(f"  {f:16s} {'PASS' if good else 'FAIL'}  max |restart - straight| = {diff:.3e}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
