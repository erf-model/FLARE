#!/usr/bin/env python3
"""Checks on the FireWindColumnZ0 regression test.

    python3 check_firewindcolumnz0.py write_map
    python3 check_firewindcolumnz0.py [parity] [sample] [ft20] [m10] [fmap] [bil] [restart]

run_wind_column_z0.sh writes, for each variant v, run_<v>.log and
plt_fire_<v>_00004. Every pair of variants shares one atmosphere (one-way
coupling) and differs only in erf.fire.wind_sample_z0_source, so the ratio of
the two fire fields at a cell is a closed form of the cell's roughness z0c:
0.03 m where the atmospheric column index is below 32 (x < 160 m), 0.5 m above
(roughness_halves.txt), against the scalar erf.fire.wind_sample_z0 = 0.1 m.
With F(z_to, z_from, z0) = ln(z_to / z0) / ln(z_from / z0):

write_map  writes roughness_halves.txt, the (x, y, z0) node triples that
           erf.most.roughness_file_name reads, i fastest, for the 64 x 32 grid,
           and fuel_stripes.asc, a fire-grid fuel map of Anderson 1 and 3 in
           stripes 16 fire cells wide along x.
parity     inputs_base (z0 = 0.1 m everywhere) with both sources: every field
           the check reads is identical, bit for bit.
sample     inputs_map, the wind sampled at 20 m and brought to 6.1 m:
           fire_wind_ref (surface_layer) / fire_wind_ref (scalar)
               = F(6.1, 20, z0c) / F(6.1, 20, 0.1)        to 1e-12,
           for both components (the sounding's wind is u = 8, v = 4 m/s).
ft20       inputs_map with erf.fire.wind_ref_ht = 10: the wind sampled at 20 m
           and brought to 10 m (as in sample), then to 20 ft before the WAF:
           (fire_wind_eff / fire_wind_ref) (surface_layer)
             / (fire_wind_eff / fire_wind_ref) (scalar)
             = F(6.096, 10, z0c) / F(6.096, 10, 0.1)    to 1e-12.
m10        inputs_map with erf.fire.ros_model = cheney_gould, whose rate reads
           the reference wind brought to 10 m, U10 = |fire_wind_ref| F(10, 6.1, z0).
           Above 5 km/h and below 12 % moisture its rate is
           (1.4 + 0.838 (U - 5)^0.844) phi_M phi_C with U = 3.6 U10 [km/h]
           (Cheney, Gould and Catchpole 1998), so
           fire_ros (surface_layer) / fire_ros (scalar)
             = R(U10 with z0c) / R(U10 with 0.1 m)       to 1e-10,
           each U10 taken from its own run's fire_wind_ref. With
           moisture_dynamic = false the rate reads
           erf.fire.cheney_gould.moisture (default 10 %).
fmap       ft20 on the fuel map: the 20 ft factor of the spatial-fuel branch.
bil        inputs_map, wind_interp = bilinear, no two-height sampling, the
           2 m reference wind below the 2.5 m first centre, wind_below_first_cell
           = log: fire_wind_ref (surface_layer) / fire_wind_ref (scalar)
             = F(2, 2.5, z0c) / F(2, 2.5, 0.1)    to 1e-12,
           over the cells whose four stencil columns share one z0 (not in
           atmospheric columns 0, 31, 32 or 63, next to the two seams).
restart    bulk_coeff fluxes (z0 recovered every step): the restart from step
           2 gives the straight run's step-4 fire fields, bit for bit.

Standard library only (erf_plotfile.py is copied beside it by CTest).
"""
import glob
import math
import sys

NX, NY = 64, 32            # atmospheric columns
DX = 5.0                   # m
GRID_RATIO = 4             # fire cells per atmospheric column
Z0_SCALAR = 0.1            # erf.fire.wind_sample_z0
Z0_WEST, Z0_EAST = 0.03, 0.5
Z_SAMPLE = 20.0            # erf.fire.wind_sample_ht
results = []


def check(name, ok, detail):
    results.append(bool(ok))
    print(f"  {name:30s} {'PASS' if ok else 'FAIL'}  {detail}")


def column_z0(i_atm):
    return Z0_WEST if i_atm < NX // 2 else Z0_EAST


def F(z_to, z_from, z0):
    return math.log(z_to / z0) / math.log(z_from / z0)


def write_map():
    with open("roughness_halves.txt", "w") as f:
        for j in range(NY + 1):
            for i in range(NX + 1):
                f.write(f"{i * DX:.1f} {j * DX:.1f} {column_z0(i)}\n")
    nfx, nfy = NX * GRID_RATIO, NY * GRID_RATIO
    with open("fuel_stripes.asc", "w") as f:
        f.write(f"ncols {nfx}\nnrows {nfy}\nxllcorner 0.0\nyllcorner 0.0\n"
                f"cellsize {DX / GRID_RATIO}\nNODATA_value -9999\n")
        row = " ".join("1" if (i // 16) % 2 == 0 else "3" for i in range(nfx))
        for _ in range(nfy):
            f.write(row + "\n")


def plotfile(v):
    pfs = sorted(glob.glob(f"plt_fire_{v}_?????"))
    return pfs[-1] if pfs else None


def fields(v, names):
    import erf_plotfile as pf
    p = plotfile(v)
    if p is None:
        return None
    _, data = pf.read_fields(p, names)
    return {n: [[col[0] for col in row] for row in data[n]] for n in names}


def cells(d):
    first = next(iter(d.values()))
    for i in range(len(first)):
        for j in range(len(first[0])):
            yield i, j


def pair(sc, sl, names):
    a, b = fields(sc, names), fields(sl, names)
    if a is None or b is None:
        check(f"plotfiles {sc}/{sl}", False, "missing")
        return None, None
    return a, b


def parity():
    names = ["fire_wind_ref_u", "fire_wind_ref_v", "fire_wind_eff_u", "fire_wind_eff_v", "fire_ros", "fire_phi"]
    a, b = pair("uni_sc", "uni_sl", names)
    if a is None:
        return
    ndiff = sum(1 for n in names for i, j in cells(a) if a[n][i][j] != b[n][i][j])
    check("parity uniform z0", ndiff == 0, f"{ndiff} values differ over {len(names)} fields")


def ratio_check(name, a, b, expect, num, den, rtol, skip=lambda i: False):
    """max over cells of |num(b)/den(b) / (num(a)/den(a)) - expect| / expect"""
    worst, n, west, east = 0.0, 0, None, None
    for i, j in cells(a):
        if skip(i) or abs(den(a, i, j)) < 0.5 or abs(den(b, i, j)) < 0.5 or abs(num(a, i, j)) < 0.5:
            continue
        r = (num(b, i, j) / den(b, i, j)) / (num(a, i, j) / den(a, i, j))
        e = expect(i)
        worst = max(worst, abs(r / e - 1.0))
        n += 1
        if i // GRID_RATIO < NX // 2:
            west = r
        else:
            east = r
    check(name, n > 1000 and worst <= rtol and west is not None and east is not None,
          f"{n} cells, max rel. error {worst:.2e} (tol {rtol:g}); west {west}, east {east}")


def sample():
    a, b = pair("map_sc", "map_sl", ["fire_wind_ref_u", "fire_wind_ref_v"])
    if a is None:
        return
    expect = lambda i: F(6.1, Z_SAMPLE, column_z0(i // GRID_RATIO)) / F(6.1, Z_SAMPLE, Z0_SCALAR)
    one = lambda d, i, j: 1.0
    ratio_check("sampled wind, u", a, b, expect, lambda d, i, j: d["fire_wind_ref_u"][i][j], one, 1e-12)
    ratio_check("sampled wind, v", a, b, expect, lambda d, i, j: d["fire_wind_ref_v"][i][j], one, 1e-12)


def ft20(v="map10"):
    names = ["fire_wind_ref_u", "fire_wind_eff_u"]
    a, b = pair(f"{v}_sc", f"{v}_sl", names)
    if a is None:
        return
    expect_ref = lambda i: F(10.0, Z_SAMPLE, column_z0(i // GRID_RATIO)) / F(10.0, Z_SAMPLE, Z0_SCALAR)
    one = lambda d, i, j: 1.0
    ratio_check(f"{v}: sampled wind to 10 m, u", a, b, expect_ref, lambda d, i, j: d["fire_wind_ref_u"][i][j], one, 1e-12)
    expect = lambda i: F(6.096, 10.0, column_z0(i // GRID_RATIO)) / F(6.096, 10.0, Z0_SCALAR)
    ratio_check(f"{v}: 20 ft conversion, eff/ref", a, b, expect,
                lambda d, i, j: d["fire_wind_eff_u"][i][j], lambda d, i, j: d["fire_wind_ref_u"][i][j], 1e-12)


def fmap():
    ft20("fmap")


def bil():
    a, b = pair("bil_sc", "bil_sl", ["fire_wind_ref_u", "fire_wind_ref_v"])
    if a is None:
        return
    expect = lambda i: F(2.0, 2.5, column_z0(i // GRID_RATIO)) / F(2.0, 2.5, Z0_SCALAR)
    seam = lambda i: (i // GRID_RATIO) in (0, NX // 2 - 1, NX // 2, NX - 1)
    one = lambda d, i, j: 1.0
    ratio_check("bilinear below first, u", a, b, expect, lambda d, i, j: d["fire_wind_ref_u"][i][j], one, 1e-12, seam)
    ratio_check("bilinear below first, v", a, b, expect, lambda d, i, j: d["fire_wind_ref_v"][i][j], one, 1e-12, seam)


def restart():
    names = ["fire_wind_ref_u", "fire_wind_ref_v", "fire_wind_eff_u", "fire_wind_eff_v", "fire_ros", "fire_phi"]
    a, b = pair("rs_a", "rs_b", names)
    if a is None:
        return
    ndiff = sum(1 for n in names for i, j in cells(a) if a[n][i][j] != b[n][i][j])
    check("restart, bulk_coeff z0", ndiff == 0, f"{ndiff} values differ over {len(names)} fields")


def cg98(u10):
    U = 3.6 * u10
    return 1.4 + 0.838 * (U - 5.0) ** 0.844


def m10():
    names = ["fire_wind_ref_u", "fire_wind_ref_v", "fire_ros"]
    a, b = pair("cg_sc", "cg_sl", names)
    if a is None:
        return
    worst, n, umin = 0.0, 0, 1e30
    for i, j in cells(a):
        z0c = column_z0(i // GRID_RATIO)
        u_sc = math.hypot(a["fire_wind_ref_u"][i][j], a["fire_wind_ref_v"][i][j]) * F(10.0, 6.1, Z0_SCALAR)
        u_sl = math.hypot(b["fire_wind_ref_u"][i][j], b["fire_wind_ref_v"][i][j]) * F(10.0, 6.1, z0c)
        umin = min(umin, 3.6 * u_sc, 3.6 * u_sl)
        if 3.6 * min(u_sc, u_sl) <= 5.0 or a["fire_ros"][i][j] <= 0.0:
            continue
        r = b["fire_ros"][i][j] / a["fire_ros"][i][j]
        worst = max(worst, abs(r / (cg98(u_sl) / cg98(u_sc)) - 1.0))
        n += 1
    check("10 m wind, CG98 rate", n > 1000 and worst <= 1e-10,
          f"{n} cells, max rel. error {worst:.2e} (tol 1e-10); lowest U10 {umin:.1f} km/h")


if __name__ == "__main__":
    which = sys.argv[1:] or ["parity", "sample", "ft20", "m10", "fmap", "bil", "restart"]
    if which == ["write_map"]:
        write_map()
        sys.exit(0)
    for w in which:
        {"parity": parity, "sample": sample, "ft20": ft20, "m10": m10,
         "fmap": fmap, "bil": bil, "restart": restart}[w]()
    sys.exit(0 if results and all(results) else 1)
