#!/usr/bin/env python3
"""Check a Conductors_FireHeating run against an independent calculation.

The script recomputes, from the deck's numbers alone (no code shared with ERF), at every node of the span:

1. the air temperature at the node, T = theta - g z / c_p for the still, isentropic sounding (theta 300 K,
   surface pressure 1000 hPa), checked against the value the run logged at the hottest node;
2. the fire's irradiance at every step, from the fire cells inside the prescribed heat disc (point sources
   of chi_r q'' dA at the source height, within the cutoff) and the node positions the run wrote at the
   step's start (S1_nodes.dat);
3. the node's temperature, from its steady temperature without the fire, stepped with backward Euler and
   the IEEE 738-2023 terms written out here (Joule heating with a linear resistance, forced and natural
   convection, radiation to surroundings at the node's first air temperature), in still air;
4. from these, the hottest node's temperature, the span's mean (trapezoid) and reference temperatures, the
   thermal strain, and the sag at the end from the bundled MoorDyn stub's elastic parabola.

Pure Python (no numpy). Exit status 0 when every check passes.
"""
import argparse
import math
import sys

SIGMA_17_8 = 17.8          # IEEE 738-2023 eq 21: 17.8 D eps [((T+273)/100)^4 - ((Ta+273)/100)^4]


def air(t_film, elev):
    """IEEE 738-2023 eqs 13, 15, 17: viscosity, density, conductivity of air at the film temperature."""
    mu = 1.458e-6 * (t_film + 273.0) ** 1.5 / (t_film + 383.4)
    rho = (1.293 - 1.525e-4 * elev + 6.379e-9 * elev ** 2) / (1.0 + 0.00367 * t_film)
    k = 2.424e-2 + 7.477e-5 * t_film - 4.407e-9 * t_film ** 2
    return mu, rho, k


def k_angle(phi):
    """IEEE 738-2023 eq 19."""
    return 1.194 - math.cos(phi) + 0.194 * math.cos(2 * phi) + 0.368 * math.sin(2 * phi)


def convection(D, V, phi, T, Ta, elev):
    """The larger in magnitude of forced (larger of eqs 9, 10) and natural (eq 7) convection, signed by T - Ta."""
    mu, rho, k = air(0.5 * (T + Ta), elev)
    dT = T - Ta
    Re = D * rho * V / mu
    f1 = k_angle(phi) * (1.01 + 1.35 * Re ** 0.52) * k * dT
    f2 = k_angle(phi) * 0.754 * Re ** 0.6 * k * dT
    forced = f1 if abs(f1) >= abs(f2) else f2
    natural = math.copysign(3.645 * math.sqrt(rho) * D ** 0.75 * abs(dT) ** 1.25, dT) if dT != 0 else 0.0
    return forced if abs(forced) >= abs(natural) else natural


def radiation(D, eps, T, Ts):
    return SIGMA_17_8 * D * eps * (((T + 273.0) / 100.0) ** 4 - ((Ts + 273.0) / 100.0) ** 4)


def gain(a, T, Ta, V, phi, E, Ts=None):
    """Net heat gain [W/m]; the conductor radiates to surroundings at Ts (Ta when not given)."""
    Ts = Ta if Ts is None else Ts
    R = a.R[0] + (a.R[1] - a.R[0]) * (T - a.RT[0]) / (a.RT[1] - a.RT[0])
    return (a.current ** 2 * R + a.alpha * a.D * (a.solar + E)
            - convection(a.D, V, phi, T, Ta, a.elev) - radiation(a.D, a.eps, T, Ts))


def root(g, lo, hi):
    """Bisection on an increasing g with g(lo) < 0 < g(hi), to 1e-12 K."""
    glo, ghi = g(lo), g(hi)
    assert glo < 0 < ghi, (lo, hi, glo, ghi)
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if g(m) < 0:
            lo = m
        else:
            hi = m
        if hi - lo < 1e-12:
            break
    return 0.5 * (lo + hi)


def steady(a, Ta, V, phi, E):
    return root(lambda T: -gain(a, T, Ta, V, phi, E), Ta - 1.0, Ta + 1000.0)


def implicit(a, Told, dt, Ta, V, phi, E, Ts=None):
    return root(lambda T: a.mcp * (T - Told) - dt * gain(a, T, Ta, V, phi, E, Ts), Told - 200.0, Told + 1000.0)


def air_temperature(z):
    """The still, isentropic sounding (theta 300 K, 1000 hPa at z = 0) at height z [C]: T = theta - g z / c_p."""
    return 300.0 - 9.81 * z / 1004.5 - 273.15


def fire_cells(a):
    """The fire cells whose centres lie in the prescribed disc, as (x, y, Q)."""
    n = int(round(a.extent / a.fire_dx))
    out = []
    for j in range(n):
        y = (j + 0.5) * a.fire_dx
        if abs(y - a.cy) > a.radius + a.fire_dx:
            continue
        for i in range(n):
            x = (i + 0.5) * a.fire_dx
            if (x - a.cx) ** 2 + (y - a.cy) ** 2 <= a.radius ** 2:
                out.append((x, y, a.flux * a.fire_dx ** 2))
    return out


def irradiance(a, cells, node, prev, nxt):
    """chi_r sum Q sin(psi) / (4 pi max(R^2, dA/pi)) at node, the axis from its neighbours."""
    t = [nxt[c] - prev[c] for c in range(3)]
    tl = math.sqrt(sum(v * v for v in t))
    t = [v / tl for v in t]
    rmin2 = a.fire_dx ** 2 / math.pi
    s = 0.0
    for (x, y, Q) in cells:
        r = [node[0] - x, node[1] - y, node[2] - a.source_height]
        r2 = sum(v * v for v in r)
        if r2 > a.cutoff ** 2:
            continue
        cx = t[1] * r[2] - t[2] * r[1]
        cy = t[2] * r[0] - t[0] * r[2]
        cz = t[0] * r[1] - t[1] * r[0]
        sin_psi = min(1.0, math.sqrt((cx * cx + cy * cy + cz * cz) / r2))
        s += Q * sin_psi / (4.0 * math.pi * max(r2, rmin2))
    return a.chi * s


def parabola_sag(a, L):
    """The stub's span: chord c + 8 d^2 / (3 c) = L (1 + H / EA), H = w c^2 / (8 d)."""
    w = (a.mass - a.rho_air * 0.25 * math.pi * a.D ** 2) * 9.81
    c = a.chord

    def f(d):
        H = w * c * c / (8.0 * d)
        return c + 8.0 * d * d / (3.0 * c) - L * (1.0 + H / a.EA)
    return root(f, 1.0e-3, 0.5 * c)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dir", default=".")
    p.add_argument("--dt", type=float, required=True, help="the lines' step [s]")
    p.add_argument("--fire-lags", type=int, required=True, help="1: the first step sees no fire yet")
    p.add_argument("--fire-dx", type=float, required=True, help="fire cell size [m]")
    p.add_argument("--extent", type=float, default=1500.0)
    p.add_argument("--flux", type=float, default=2.0e5)
    p.add_argument("--radius", type=float, default=40.0)
    p.add_argument("--cx", type=float, default=750.0)
    p.add_argument("--cy", type=float, default=500.0)
    p.add_argument("--chi", type=float, default=0.3)
    p.add_argument("--cutoff", type=float, default=300.0)
    p.add_argument("--source-height", type=float, default=2.0)
    p.add_argument("--D", type=float, default=0.0281)
    p.add_argument("--alpha", type=float, default=0.5)
    p.add_argument("--eps", type=float, default=0.5)
    p.add_argument("--mcp", type=float, default=1310.0)
    p.add_argument("--current", type=float, default=600.0)
    p.add_argument("--R", type=float, nargs=2, default=[7.283e-5, 8.688e-5])
    p.add_argument("--RT", type=float, nargs=2, default=[25.0, 75.0])
    p.add_argument("--solar", type=float, default=0.0)
    p.add_argument("--elev", type=float, default=0.0)
    p.add_argument("--expansion", type=float, default=1.93e-5)
    p.add_argument("--length", type=float, default=301.5)
    p.add_argument("--chord", type=float, default=300.0)
    p.add_argument("--mass", type=float, default=1.628)
    p.add_argument("--EA", type=float, default=3.0e7)
    p.add_argument("--rho-air", type=float, default=1.17)
    p.add_argument("--min-rise", type=float, default=20.0, help="the hottest node must heat by at least this [K]")
    a = p.parse_args()

    rows = [l.split() for l in open(a.dir + "/conductors/conductor_temperature.dat") if not l.startswith("#")]
    nodes = {}
    for l in open(a.dir + "/S1_nodes.dat"):
        f = l.split()
        if f[0] == "time":
            continue
        nodes.setdefault(round(float(f[0]), 9), {})[int(f[1])] = [float(v) for v in f[2:5]]
    cells = fire_cells(a)
    nn = max(nodes[round(float(rows[0][0]) - a.dt, 9)]) + 1      # the span's nodes (no strings)
    bad = []
    T = Ts = None
    err = {"Ta": 0.0, "E": 0.0, "Tmax": 0.0, "Tmean": 0.0, "strain": 0.0}
    for n, r in enumerate(rows):
        t = float(r[0])
        start = nodes[round(t - a.dt, 9)]
        Ta = [air_temperature(start[m][2]) for m in range(nn)]
        if T is None:
            # the first step: each node at its steady temperature without the fire, radiating to its air then
            T = [steady(a, Ta[m], 0.0, math.pi / 2, 0.0) for m in range(nn)]
            Ts = list(Ta)
            Tref = (0.5 * (T[0] + T[-1]) + sum(T[1:-1])) / (nn - 1)
            T0 = list(T)
        fire_on = not (n == 0 and a.fire_lags)
        E = [irradiance(a, cells, start[m], start[max(m - 1, 0)], start[min(m + 1, nn - 1)]) if fire_on else 0.0
             for m in range(nn)]
        T = [implicit(a, T[m], a.dt, Ta[m], 0.0, math.pi / 2, E[m], Ts[m]) for m in range(nn)]
        hot = max(range(nn), key=lambda m: T[m])
        Tmean = (0.5 * (T[0] + T[-1]) + sum(T[1:-1])) / (nn - 1)
        strain = a.expansion * (Tmean - Tref)
        if int(r[5]) != hot:
            bad.append("row %d: the hottest node is %s, the calculation gives %d" % (n, r[5], hot))
            break
        err["Ta"] = max(err["Ta"], abs(float(r[6]) - Ta[hot]))
        err["E"] = max(err["E"], abs(float(r[9]) - E[hot]) / max(1.0, E[hot]))
        err["Tmax"] = max(err["Tmax"], abs(float(r[4]) - T[hot]))
        err["Tmean"] = max(err["Tmean"], abs(float(r[3]) - Tmean))
        err["strain"] = max(err["strain"], abs(float(r[11]) - strain))
    tol = {"Ta": 1e-4, "E": 1e-6, "Tmax": 1e-3, "Tmean": 1e-3, "strain": 1e-7}
    for k in tol:
        if err[k] > tol[k]:
            bad.append("%s: largest difference %.3e (tolerance %.0e)" % (k, err[k], tol[k]))
    if abs(float(rows[-1][12]) - Tref) > 1e-3:
        bad.append("reference temperature %s C, the calculation gives %.6f C" % (rows[-1][12], Tref))
    if T[hot] - T0[hot] < a.min_rise:
        bad.append("the hottest node heated by %.2f K only (at least %.1f)" % (T[hot] - T0[hot], a.min_rise))
    # the sag at the last row, from the stub's parabola at the length of that row's strain
    s1 = [l.split() for l in open(a.dir + "/S1.dat") if not l.startswith("time")]
    last = s1[-1]
    sag_py = parabola_sag(a, a.length * (1.0 + strain))
    sag_run, sag0 = float(last[4]), float(s1[0][4])
    if abs(float(last[0]) - float(rows[-1][0])) > 1e-9:
        bad.append("S1.dat ends at %s, the temperature log at %s" % (last[0], rows[-1][0]))
    if abs(sag_run - sag_py) > 1e-4:
        bad.append("sag %.9f m, the parabola gives %.9f m" % (sag_run, sag_py))
    if sag_run - sag0 < 0.05:
        bad.append("the span sagged by %.4f m only" % (sag_run - sag0))
    print("hottest node %d: %.4f C -> %.4f C over %d steps; span mean %.4f C (reference %.4f C); irradiance %.6g W/m2 at "
          "the end; largest differences: air %.1e K, E %.1e, T_max %.1e K, T_mean %.1e K, strain %.1e; sag %.6f -> %.6f m "
          "(parabola %.6f m)" % (hot, T0[hot], T[hot], len(rows), Tmean, Tref, E[hot], err["Ta"], err["E"], err["Tmax"],
                                 err["Tmean"], err["strain"], sag0, sag_run, sag_py))
    for b in bad:
        print("FAIL:", b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
