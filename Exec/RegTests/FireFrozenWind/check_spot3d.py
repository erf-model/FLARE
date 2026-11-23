#!/usr/bin/env python3
"""Check the ember drift of the FireFrozenWind spotting decks.

    python3 check_spot3d.py --drift frozen|profile [--plotfile plt_fire_00020]
                            [--interp bilinear|nearest] [--below clamp|log]

Reads the fire plotfile (step 20 by default) in the working directory and
replays, independently of FLARE, the flight of the brand each source cell
launched at that step. The replay's landing distance is compared with the
source's fire_spot_max_dist. Exits 1 if any check fails. Standard library
only (plus erf_plotfile.py and check_frozen_wind.py, next to it).

The flight, as documented for erf.fire.spotting.drift_wind:
  - the brand starts at the source cell centre, fire_lofting_height above
    the fire ground there; it falls at v_t (terminal_velocity) in
    n_traj_steps forward-Euler steps of dt = (z_start - z_min) / v_t / N,
    z_min the lowest fire ground
  - the fire ground is the GUST column ground (mean of its four k = 0
    nodes in gust_hill_plt/Level_0/z_phys_nd) of the column holding the fire
    cell centre (the ERF and GUST meshes coincide in this case); between
    fire cell centres it is bilinear, clamped at the edges
  - each step: the brand leaves the run if its cell-centre index leaves the
    interior of the fire grid; else it moves by (u, v) dt and falls v_t dt,
    and it lands when it is at or below the ground under it
  - frozen: (u, v) is GUST's wind at h = height above the fire ground, each
    of the four GUST columns around the brand sampled at its own ground + h,
    linearly between its cell centres, the lowest centre's wind below it and
    the top centre's wind above it, bilinear over the columns
    (erf.fire.wind_interp = bilinear); with --interp nearest only the column
    holding the brand; with --below log (erf.fire.wind_below_first_cell =
    log) the lowest centre's wind is scaled by ln(h' / z0) / ln(h1 / z0)
    below it, h' and h1 the heights above that column's ground and z0 =
    0.1 m (wind_sample_z0)
  - profile: (u, v) = fire_wind_ref x ln(h / z0) / ln(z_ref / z0), z_ref =
    6.1 m and z0 = 0.1 m (inputs_base), fire_wind_ref bilinear
  - the distance is capped at 200 m, the Scott cap of fuel model 1

Checks:
  sources   the step launched brands (fire_spot_brand_count > 0)
  uncapped  every landing is short of the 200 m cap, so the drift is tested
  replay    |fire_spot_max_dist - replay| <= 1e-6 m for every source
  other     the other drift's replay differs by more than 1 m somewhere
            (the check can tell the two drifts apart)
  options   frozen with a non-default --interp or --below: the replay with the
            default bilinear and clamp differs by more than 1 m somewhere, so
            the run read the option (and, for --below log, some brand sampled
            below a lowest centre)
"""
import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import erf_plotfile  # noqa: E402
from check_frozen_wind import ground_nodes  # noqa: E402

V_T = 2.5           # inputs_spot3d: terminal_velocity
N_TRAJ = 20         # inputs_spot3d: n_traj_steps
Z_REF, Z0 = 6.1, 0.1   # inputs_base: wind_ref_ht; wind_sample_z0 (default)
TINY = 1.0e-30
CAP = 200.0         # Scott cap of Anderson model 1
TOL = 1.0e-6

results = []


def check(name, ok, detail):
    results.append(bool(ok))
    print(f"  {name:8s} {'PASS' if ok else 'FAIL'}  {detail}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drift", choices=("frozen", "profile"), required=True)
    ap.add_argument("--plotfile", default="plt_fire_00020")
    ap.add_argument("--interp", choices=("bilinear", "nearest"), default="bilinear")
    ap.add_argument("--below", choices=("clamp", "log"), default="clamp")
    opt = ap.parse_args()

    # GUST columns
    ghdr, g = erf_plotfile.read_fields("gust_hill_plt", ["u", "v", "z"])
    U, V, Z = g["u"], g["v"], g["z"]
    gnx, gny, gnz = len(U), len(U[0]), len(U[0][0])
    gx0, gy0 = ghdr["prob_lo"][0], ghdr["prob_lo"][1]
    gdx = (ghdr["prob_hi"][0] - gx0) / gnx
    gdy = (ghdr["prob_hi"][1] - gy0) / gny
    zn = ground_nodes("gust_hill_plt")
    G = [[sum(zn[(i + a, j + b)] for a in (0, 1) for b in (0, 1)) / 4.0
          for j in range(gny)] for i in range(gnx)]

    z0s = abs(Z0) + TINY
    n_below = [0]

    def column(ia, ja, zt, below):
        zc = Z[ia][ja]
        if zt <= zc[0]:
            n_below[0] += 1
            f = 1.0
            if below == "log":
                zg = G[ia][ja]
                z1 = max(zc[0] - zg, z0s * 1.000001)
                f = math.log(max(min(zt - zg, z1), z0s) / z0s) / math.log(z1 / z0s)
            return f * U[ia][ja][0], f * V[ia][ja][0]
        k = 0
        while k < gnz - 2 and zc[k + 1] <= zt:
            k += 1
        a = min(1.0, max(0.0, (zt - zc[k]) / (zc[k + 1] - zc[k])))
        return (U[ia][ja][k] + a * (U[ia][ja][k + 1] - U[ia][ja][k]),
                V[ia][ja][k] + a * (V[ia][ja][k + 1] - V[ia][ja][k]))

    def frozen_wind(x, y, h, interp, below):
        if interp == "nearest":
            ia = max(0, min(math.floor((x - gx0) / gdx), gnx - 1))
            ja = max(0, min(math.floor((y - gy0) / gdy), gny - 1))
            return column(ia, ja, G[ia][ja] + h, below)
        gx = (x - gx0) / gdx - 0.5
        gy = (y - gy0) / gdy - 0.5
        i0, j0 = math.floor(gx), math.floor(gy)
        wx, wy = gx - i0, gy - j0
        u = v = 0.0
        for di in (0, 1):
            for dj in (0, 1):
                ia = max(0, min(i0 + di, gnx - 1))
                ja = max(0, min(j0 + dj, gny - 1))
                w = (wx if di else 1.0 - wx) * (wy if dj else 1.0 - wy)
                if w <= 0.0:
                    continue
                uc, vc = column(ia, ja, G[ia][ja] + h, below)
                u += w * uc
                v += w * vc
        return u, v

    # fire fields
    if not os.path.isdir(opt.plotfile):
        check("sources", False, f"{opt.plotfile} missing")
        print("FAIL")
        return 1
    fields = ["fire_wind_ref_u", "fire_wind_ref_v", "fire_lofting_height",
              "fire_spot_brand_count", "fire_spot_max_dist"]
    hdr, d = erf_plotfile.read_fields(opt.plotfile, fields)
    WU, WV, HZ, CNT, DST = ([[c[0] for c in row] for row in d[f]] for f in fields)
    nx, ny = len(WU), len(WU[0])
    x0, y0 = hdr["prob_lo"][0], hdr["prob_lo"][1]
    dx = (hdr["prob_hi"][0] - x0) / nx
    dy = (hdr["prob_hi"][1] - y0) / ny

    def col_of(x, y):
        return (max(0, min(math.floor((x - gx0) / gdx), gnx - 1)),
                max(0, min(math.floor((y - gy0) / gdy), gny - 1)))

    T = [[G[col_of(x0 + (i + 0.5) * dx, y0 + (j + 0.5) * dy)[0]]
           [col_of(x0 + (i + 0.5) * dx, y0 + (j + 0.5) * dy)[1]]
          for j in range(ny)] for i in range(nx)]
    zmin = min(min(r) for r in T)

    def bilin(f, x, y):
        fx = (x - x0) / dx - 0.5
        fy = (y - y0) / dy - 0.5
        i0, j0 = math.floor(fx), math.floor(fy)
        wx, wy = fx - i0, fy - j0

        def at(i, j):
            return f[max(0, min(nx - 1, i))][max(0, min(ny - 1, j))]
        return ((1 - wx) * (1 - wy) * at(i0, j0) + wx * (1 - wy) * at(i0 + 1, j0)
                + (1 - wx) * wy * at(i0, j0 + 1) + wx * wy * at(i0 + 1, j0 + 1))

    def inside(x, y):
        i0 = math.floor((x - x0) / dx - 0.5)
        j0 = math.floor((y - y0) / dy - 0.5)
        return 0 <= i0 < nx - 1 and 0 <= j0 < ny - 1

    def profile_factor(h):
        z0s = max(Z0, 1.0e-6)
        zr = max(Z_REF, z0s * 1.000001)
        return math.log(max(h, z0s) / z0s) / math.log(zr / z0s)

    def fly(si, sj, H, drift, interp=opt.interp, below=opt.below):
        cx, cy = x0 + (si + 0.5) * dx, y0 + (sj + 0.5) * dy
        zb = T[si][sj] + H
        dt = (zb - zmin) / V_T / N_TRAJ
        dz = V_T * dt
        px, py = cx, cy
        for _ in range(N_TRAJ):
            if not inside(px, py):
                return None
            h = zb - bilin(T, px, py)
            if drift == "frozen":
                u, v = frozen_wind(px, py, h, interp, below)
            else:
                f = profile_factor(h)
                u, v = f * bilin(WU, px, py), f * bilin(WV, px, py)
            px += u * dt
            py += v * dt
            zb -= dz
            if zb <= bilin(T, px, py):
                break
        return min(math.hypot(px - cx, py - cy), CAP)

    other = "profile" if opt.drift == "frozen" else "frozen"
    src = [(i, j) for i in range(nx) for j in range(ny) if CNT[i][j] > 0.0]
    check("sources", len(src) > 0, f"{len(src)} source cells at {opt.plotfile}")
    if not src:
        print("FAIL")
        return 1
    err = 0.0
    gap = 0.0
    gap_opt = 0.0
    below_seen = 0
    dmax = 0.0
    lost = 0
    for i, j in src:
        D = DST[i][j]
        dmax = max(dmax, D)
        n_below[0] = 0
        r = fly(i, j, HZ[i][j], opt.drift)
        below_seen += n_below[0]
        if r is None:
            lost += 1
            err = math.inf
            continue
        err = max(err, abs(D - r))
        if opt.drift == "frozen" and (opt.interp, opt.below) != ("bilinear", "clamp"):
            rd = fly(i, j, HZ[i][j], "frozen", "bilinear", "clamp")
            if rd is not None:
                gap_opt = max(gap_opt, abs(rd - r))
        ro = fly(i, j, HZ[i][j], other)
        if ro is not None:
            gap = max(gap, abs(ro - r))
    check("uncapped", dmax < CAP, f"longest landing {dmax:.3f} m (cap {CAP} m)")
    check("replay", err <= TOL and lost == 0,
          f"max |fire_spot_max_dist - {opt.drift} replay| = {err:.3e} m (tol {TOL}); "
          f"{lost} replays left the domain")
    check("other", gap > 1.0, f"the {other} replay differs by up to {gap:.3f} m")
    if opt.drift == "frozen" and (opt.interp, opt.below) != ("bilinear", "clamp"):
        ok_opt = gap_opt > 1.0 and (opt.below != "log" or below_seen > 0)
        check("options", ok_opt,
              f"the bilinear/clamp replay differs by up to {gap_opt:.3f} m; "
              f"{below_seen} samples below a lowest centre")

    ok = all(results)
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
