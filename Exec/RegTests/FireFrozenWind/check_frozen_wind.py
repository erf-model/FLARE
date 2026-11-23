#!/usr/bin/env python3
"""Check the frozen GUST wind of the FireFrozenWind case.

    python3 check_frozen_wind.py [--interp bilinear|nearest] [--zref Z]
                                 [--below clamp|log] [--zsample Z] [--z0 Z0]

Reads the fire plotfiles plt_fire_00010 and plt_fire_00020 in the working
directory and compares fire_wind_ref_u/v and fire_extract_z with the
reference wind sampled here, independently of FLARE, from the GUST plotfile
gust_hill_plt. The checkpoint deck is checked against the same plotfile: the
checkpoint holds the same field on faces. Exits 1 if any check fails.
Standard library only (plus erf_plotfile.py, copied next to it by CTest).

The sampling, as documented for erf.fire.frozen_wind.file:
  - fire cell centre (x, y), 10 m cells; GUST columns 40 m wide
  - the four GUST columns around (x, y), bilinear weights, indices clamped
    to the GUST domain
  - each column's ground is the mean of its four k = 0 node heights, read
    here from GUST's Level_0/z_phys_nd (they follow the analytic hill of
    make_hill_terrain.py to 3e-4 m: GUST interpolates the terrain file as
    ERF does)
  - each column is sampled at ground + 6.1 m, linearly between the cell
    centres that bracket it (6.1 m is above the first centre, 2.5 m)
  - the options give the variants of run_variants.sh: the nearest column
    (the one holding the fire cell centre); a height below the first centre
    with the wind scaled by log(z / z0) / log(z1 / z0), z1 the first centre
    above the ground; and a wind sampled at zsample and brought to zref by
    log(zref / z0) / log(zsample / z0)

Checks:
  wind     |fire_wind_ref - sampled| <= 1e-9 m/s for u and v, every fire cell
  height   |fire_extract_z - (ground + 6.1)| <= 1e-9 m
  gust     the wind is GUST's, not the atmosphere's: the mean v is > 0
           (GUST's Ekman turning) where the sounding has v = -5 m/s
  frozen   steps 10 and 20 hold the same reference wind, bit for bit
"""
import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import erf_plotfile  # noqa: E402

H, S, XC, YC = 40.0, 150.0, 640.0, 640.0   # make_hill_terrain.py
TINY = 1.0e-30                              # z_tiny of the fire layer's log factors
TOL_WIND = 1.0e-9
TOL_Z = 1.0e-9

results = []


def check(name, ok, detail):
    results.append(bool(ok))
    print(f"  {name:8s} {'PASS' if ok else 'FAIL'}  {detail}")


def hill(x, y):
    return H * math.exp(-((x - XC) ** 2 + (y - YC) ** 2) / (2.0 * S * S))


def ground_nodes(plotfile):
    """The k = 0 node heights {(i, j): z} of <plotfile>/Level_0/z_phys_nd."""
    level = os.path.join(plotfile, "Level_0")
    with open(os.path.join(level, "z_phys_nd_H")) as fh:
        fabs = [ln.split()[1:3] for ln in fh if ln.startswith("FabOnDisk:")]
    z0 = {}
    for name, off in fabs:
        lo, hi, data = erf_plotfile._read_fab(os.path.join(level, name), int(off), 1)
        nx, ny = hi[0] - lo[0] + 1, hi[1] - lo[1] + 1
        if lo[2] > 0 or hi[2] < 0:
            continue
        k = -lo[2]
        for j in range(ny):
            for i in range(nx):
                z0[(lo[0] + i, lo[1] + j)] = data[i + nx * (j + ny * k)]
    return z0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--interp", choices=("bilinear", "nearest"), default="bilinear")
    ap.add_argument("--zref", type=float, default=6.1)      # erf.fire.wind_ref_ht in inputs_base
    ap.add_argument("--below", choices=("clamp", "log"), default="clamp")
    ap.add_argument("--zsample", type=float, default=0.0)
    ap.add_argument("--z0", type=float, default=0.1)
    opt = ap.parse_args()
    Z_REF = opt.zref
    z_at = opt.zsample if opt.zsample > 0.0 else Z_REF
    z0s = abs(opt.z0) + TINY
    scale = (math.log((abs(Z_REF) + TINY) / z0s) / (math.log((abs(opt.zsample) + TINY) / z0s) + TINY)
             if opt.zsample > 0.0 else 1.0)
    print(f"interp {opt.interp}, wind_ref_ht {Z_REF}, below {opt.below}, "
          f"wind_sample_ht {opt.zsample}, z0 {opt.z0}")

    ghdr, g = erf_plotfile.read_fields("gust_hill_plt", ["u", "v", "z"])
    U, V, Z = g["u"], g["v"], g["z"]
    gnx, gny, gnz = len(U), len(U[0]), len(U[0][0])
    gx0, gy0 = ghdr["prob_lo"][0], ghdr["prob_lo"][1]
    gdx = (ghdr["prob_hi"][0] - gx0) / gnx
    gdy = (ghdr["prob_hi"][1] - gy0) / gny

    z0 = ground_nodes("gust_hill_plt")
    dh = max(abs(z0[(i, j)] - hill(gx0 + i * gdx, gy0 + j * gdy))
             for i in range(gnx + 1) for j in range(gny + 1))
    check("terrain", dh <= 1.0e-3,
          f"GUST's k = 0 nodes follow the analytic hill to {dh:.2e} m (tol 1e-3)")

    def ground(ia, ja):
        return sum(z0[(ia + a, ja + b)] for a in (0, 1) for b in (0, 1)) / 4.0

    def column(ia, ja, zt, zg):
        zc = Z[ia][ja]
        if zt <= zc[0]:
            f = 1.0
            if opt.below == "log":
                z1 = max(zc[0] - zg, z0s * 1.000001)
                f = math.log(max(min(zt - zg, z1), z0s) / z0s) / math.log(z1 / z0s)
            return f * U[ia][ja][0], f * V[ia][ja][0]
        k = 0
        while k < gnz - 2 and zc[k + 1] <= zt:
            k += 1
        a = min(1.0, max(0.0, (zt - zc[k]) / (zc[k + 1] - zc[k])))
        return (U[ia][ja][k] + a * (U[ia][ja][k + 1] - U[ia][ja][k]),
                V[ia][ja][k] + a * (V[ia][ja][k + 1] - V[ia][ja][k]))

    def sample(x, y):
        if opt.interp == "nearest":
            ia = max(0, min(math.floor((x - gx0) / gdx), gnx - 1))
            ja = max(0, min(math.floor((y - gy0) / gdy), gny - 1))
            zt = ground(ia, ja) + z_at
            uc, vc = column(ia, ja, zt, ground(ia, ja))
            return scale * uc, scale * vc, zt
        gx = (x - gx0) / gdx - 0.5
        gy = (y - gy0) / gdy - 0.5
        i0, j0 = math.floor(gx), math.floor(gy)
        wx, wy = gx - i0, gy - j0
        u = v = z = 0.0
        for di in (0, 1):
            for dj in (0, 1):
                ia = max(0, min(i0 + di, gnx - 1))
                ja = max(0, min(j0 + dj, gny - 1))
                w = (wx if di else 1.0 - wx) * (wy if dj else 1.0 - wy)
                zt = ground(ia, ja) + z_at
                uc, vc = column(ia, ja, zt, ground(ia, ja))
                u += w * uc
                v += w * vc
                z += w * zt
        return scale * u, scale * v, z

    fields = ["fire_wind_ref_u", "fire_wind_ref_v", "fire_extract_z"]
    winds = {}
    for pf in ("plt_fire_00010", "plt_fire_00020"):
        if not os.path.isdir(pf):
            check("wind", False, f"{pf} missing")
            continue
        hdr, d = erf_plotfile.read_fields(pf, fields)
        fu, fv, fz = (d[f] for f in fields)
        nx, ny = len(fu), len(fu[0])
        x0, y0 = hdr["prob_lo"][0], hdr["prob_lo"][1]
        dx = (hdr["prob_hi"][0] - x0) / nx
        dy = (hdr["prob_hi"][1] - y0) / ny
        du = dz = sv = 0.0
        for i in range(nx):
            for j in range(ny):
                u, v, z = sample(x0 + (i + 0.5) * dx, y0 + (j + 0.5) * dy)
                du = max(du, abs(fu[i][j][0] - u), abs(fv[i][j][0] - v))
                dz = max(dz, abs(fz[i][j][0] - z))
                sv += fv[i][j][0]
        print(f"{pf}: {nx} x {ny} fire cells")
        check("wind", du <= TOL_WIND,
              f"max |fire_wind_ref - GUST sampled| = {du:.3e} m/s (tol {TOL_WIND})")
        check("height", dz <= TOL_Z,
              f"max |fire_extract_z - (ground + {z_at})| = {dz:.3e} m (tol {TOL_Z})")
        check("gust", sv / (nx * ny) > 0.0,
              f"mean v = {sv / (nx * ny):.3f} m/s (the sounding has -5)")
        winds[pf] = (fu, fv)

    if len(winds) == 2:
        (a_u, a_v), (b_u, b_v) = winds["plt_fire_00010"], winds["plt_fire_00020"]
        diff = max(abs(a[i][j][0] - b[i][j][0])
                   for a, b in ((a_u, b_u), (a_v, b_v))
                   for i in range(len(a)) for j in range(len(a[0])))
        check("frozen", diff == 0.0, f"max |wind(step 20) - wind(step 10)| = {diff:.3e} m/s")

    ok = bool(results) and all(results)
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
