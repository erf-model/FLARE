#!/usr/bin/env python3
"""Check the FARSITE gridded winds written by erf.fire.farsite_winds.

    python3 check_farsite_winds.py

Checks every output directory of this case that exists in the working
directory (winds_prescribed/ from inputs_prescribed, winds_atmosphere/ from
inputs_atmosphere) and exits 1 if any check fails. Standard library only
(plus erf_plotfile.py, copied next to it by CTest).

winds_prescribed (u = 3, v = 4 m/s, ENGLISH):
  atm      the .atm is WINDS / ENGLISH and three lines 12 30 1100, 1101, 1102
           naming files that exist, steps 1, 61 and 121
  header   ncols 60, nrows 20, cellsize 10, the map corner of the deck
  values   every speed is 11.18 (5 m/s in mph, 2 decimals) and every
           direction 217 (180 + atan2(3, 4) = 216.87, 0 decimals)
  prj      each grid has a .prj identical to utm13n.prj

winds_atmosphere (sounding u = 5, v = 2 m/s, METRIC, coarsen 2):
  atm      WINDS / METRIC and three lines
  header   ncols 30, nrows 10, cellsize 20, corner 0 0
  plotfile the step-61 grids equal the 2 x 2 block means of fire_wind_ref_u/v
           in plt_fire_00061 (speed to 0.006 km/h, direction to 0.51 deg:
           the rounding of the written values), north row first
  wind     the mean direction is within 10 deg of the sounding's 248.2 deg
           and the mean speed between 70 % and 100 % of its 19.39 km/h
           (a mph/km/h mix-up gives 62 %)
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

results = []


def check(name, ok, detail):
    results.append(bool(ok))
    print(f"  {name:9s} {'PASS' if ok else 'FAIL'}  {detail}")


def read_asc(path):
    with open(path) as f:
        lines = f.read().split("\n")
    hdr = {}
    for ln in lines[:6]:
        k, v = ln.split("\t")
        hdr[k] = float(v)
    rows = [[float(x) for x in ln.split()] for ln in lines[6:] if ln.strip()]
    return hdr, rows, lines


def read_atm(path):
    with open(path) as f:
        return [ln for ln in f.read().split("\n") if ln]


def direction_from(u, v):
    if u == 0.0 and v == 0.0:
        return 0.0
    d = math.degrees(math.atan2(-u, -v))
    return d + 360.0 if d < 0.0 else d


def ang_diff(a, b):
    return abs((a - b + 180.0) % 360.0 - 180.0)


def check_atm(d, prefix, unit, n_expected=3):
    atm = read_atm(os.path.join(d, prefix + ".atm"))
    times = ["12 30 1100", "12 30 1101", "12 30 1102"]
    steps = [1, 61, 121]
    want = ["WINDS", unit] + [f"{t} {prefix}_{s:06d}_vel.asc {prefix}_{s:06d}_ang.asc"
                               for t, s in zip(times, steps)]
    check("atm", atm == want, f"{os.path.join(d, prefix + '.atm')}: {atm}")
    files = [os.path.join(d, w) for ln in atm[2:] for w in ln.split()[3:]]
    check("files", all(os.path.isfile(f) for f in files) and len(files) == 2 * n_expected,
          f"{len(files)} grid files listed, all present")


def check_prescribed(d):
    print(f"{d}:")
    check_atm(d, "prescribed", "ENGLISH")
    mph = round(5.0 * 3600.0 / 1609.344, 2)
    ang = round(180.0 + math.degrees(math.atan2(3.0, 4.0)))
    with open(os.path.join(HERE, "utm13n.prj")) as f:
        prj_in = f.read()
    for s in (1, 61, 121):
        for kind, want, tol in (("vel", mph, 1e-9), ("ang", ang, 1e-9)):
            path = os.path.join(d, f"prescribed_{s:06d}_{kind}.asc")
            hdr, rows, lines = read_asc(path)
            ok_h = (hdr["ncols"] == 60 and hdr["nrows"] == 20 and hdr["cellsize"] == 10.0
                    and hdr["xllcorner"] == 482000.0 and hdr["yllcorner"] == 4423000.0
                    and hdr["NODATA_value"] == -9999.0)
            check("header", ok_h, f"{os.path.basename(path)}: {hdr}")
            vals = [x for r in rows for x in r]
            ok_v = (len(rows) == 20 and all(len(r) == 60 for r in rows)
                    and all(abs(x - want) <= tol for x in vals))
            check("values", ok_v, f"{os.path.basename(path)}: {len(vals)} values, "
                  f"min {min(vals)} max {max(vals)}, expected {want}")
            with open(path[:-4] + ".prj") as f:
                check("prj", f.read() == prj_in, os.path.basename(path[:-4] + ".prj"))


def check_atmosphere(d):
    print(f"{d}:")
    check_atm(d, "atm", "METRIC")
    hv, vel, _ = read_asc(os.path.join(d, "atm_000061_vel.asc"))
    ha, ang, _ = read_asc(os.path.join(d, "atm_000061_ang.asc"))
    ok_h = (hv["ncols"] == 30 and hv["nrows"] == 10 and hv["cellsize"] == 20.0
            and hv["xllcorner"] == 0.0 and hv["yllcorner"] == 0.0 and ha == hv)
    check("header", ok_h, str(hv))

    import erf_plotfile
    pf = "plt_fire_00061"
    if not os.path.isdir(pf):
        check("plotfile", False, f"{pf} missing")
        return
    hdr, data = erf_plotfile.read_fields(pf, ["fire_wind_ref_u", "fire_wind_ref_v"])
    U, V = data["fire_wind_ref_u"], data["fire_wind_ref_v"]
    nx, ny = len(U), len(U[0])
    worst_s, worst_a = 0.0, 0.0
    su = sv = 0.0
    for J in range(ny // 2):
        r = ny // 2 - 1 - J                 # file row r holds J; north first
        for I in range(nx // 2):
            u = sum(U[2 * I + a][2 * J + b][0] for a in (0, 1) for b in (0, 1)) / 4.0
            v = sum(V[2 * I + a][2 * J + b][0] for a in (0, 1) for b in (0, 1)) / 4.0
            su += u
            sv += v
            worst_s = max(worst_s, abs(vel[r][I] - 3.6 * math.hypot(u, v)))
            worst_a = max(worst_a, ang_diff(ang[r][I], direction_from(u, v)))
    check("plotfile", worst_s <= 0.006 and worst_a <= 0.51,
          f"step 61 vs {pf} 2x2 means: max speed diff {worst_s:.4f} km/h, "
          f"max direction diff {worst_a:.3f} deg")
    n = (nx // 2) * (ny // 2)
    mu, mv = su / n, sv / n
    d_mean = direction_from(mu, mv)
    s_mean = 3.6 * math.hypot(mu, mv)
    d_snd = direction_from(5.0, 2.0)
    s_snd = 3.6 * math.hypot(5.0, 2.0)
    check("wind", ang_diff(d_mean, d_snd) <= 10.0 and 0.7 * s_snd <= s_mean <= s_snd,
          f"mean from {d_mean:.1f} deg at {s_mean:.2f} km/h; sounding {d_snd:.1f} deg, {s_snd:.2f} km/h")


found = False
if os.path.isdir("winds_prescribed"):
    found = True
    check_prescribed("winds_prescribed")
if os.path.isdir("winds_atmosphere"):
    found = True
    check_atmosphere("winds_atmosphere")
if not found:
    check("output", False, "neither winds_prescribed/ nor winds_atmosphere/ exists")

print(f"{sum(results)}/{len(results)} checks passed")
sys.exit(0 if all(results) else 1)
