#!/usr/bin/env python3
"""Mean tension of one conductor span in a steady uniform crosswind, for the single-span case of
Ma, Bocchini and Christou 2020 (Struct. Saf. 82:101890): hinged supports at the same level, the
unstretched length set by the sag-to-span ratio under gravity alone.

    python3 static_reference.py --span 300 --wind 30 40 50

Two independent answers, which must agree:

  catenary  the elastic catenary of Structures/Tools/broken_wire/broken_wire.py, loaded by the
            resultant of the weight w and a uniform drag f = 0.5 rho C_D D U^2 per metre. It lies in
            the plane through the chord tilted by the swing angle atan(f / w) from the vertical.
            This takes the whole wind as normal to the cable everywhere.
  chain     a lumped chain of N straight elastic links (the model MoorDyn integrates), relaxed to
            rest in 3-D, whose drag acts only on the wind's component normal to each link and on
            the link's stretched length; its weight acts on the unstretched length.

The two differ in two ways that partly cancel. The cable's slope in its tilted plane takes the
wind's component along the cable off the chain's drag (the drag goes with the square of the
normal wind); on its own this lowers the end tension by about 0.3 % at 40 m/s, 0.63 % at 60 m/s and
1.28 % at 100 m/s. The chain's drag on the stretched length puts back about the strain times the
square of the swing's sine, 0.33 % at 100 m/s. Net, the chain's end tension is 0.03 % below the
catenary's at 20 m/s, 0.5 % at 60 m/s and 1.0 % at 100 m/s (python3 static_reference.py --span 300
--wind 20 60 100).

The last column, T_moordyn, is the chain with MoorDyn's two conventions (MoorDyn-C 2.7.1,
source/Line.cpp): the drag on the unstretched length and the weight less the displaced air's.
These take another 0.03 % (still air) to 0.33 % (100 m/s) off the end tension, and are what the
conductor model's run should reproduce (run_check.py).

The catenary and rest_length use the standard library only; the chain needs numpy.
"""
import argparse
import importlib.util
import math
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")        # a small solve: one thread, as the runs are load-gated
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "broken_wire", os.path.join(HERE, "..", "..", "Tools", "broken_wire", "broken_wire.py"))
bw = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bw)

# Ma et al. 2020: Table 3 (Drake 26/7: weight, diameter as the deck's 1107 mil OD) and Section 3
# (sag-to-span ratio, drag coefficient, air density)
DRAKE_WEIGHT = 15.95          # N/m, Ma 2020 Table 3
DIAMETER = 0.0281             # m, from Ma 2020 Table 3's strands: 3 x 3.5 mm (steel core) + 4 x 4.4 mm (aluminium)
SAG_RATIO = 1.0 / 50.0        # Ma 2020 Section 3
CD = 1.0                      # Ma 2020 Section 3
RHO = 1.226                   # kg/m^3, Ma 2020 Section 3
G = 9.81                      # m/s^2, as the conductor model hands MoorDyn
EA_DEFAULT = 3.0e7            # N, the conductor model's tests; not a supplier figure


def rest_length(span, sag_ratio, w, EA, segments=0):
    """The unstretched length whose low point hangs sag_ratio * span below the supports under its
    weight alone (bisection on the length)."""
    target = sag_ratio * span

    def low(L0):
        H, V0 = bw.solve_cable((span, 0.0), L0, w, EA, segments)
        return bw.cable_low_point(H, V0, L0, w, EA, segments)

    lo, hi = span * (1.0 - 1.0e-3), span * (1.0 + 8.0 / 3.0 * (2.0 * sag_ratio) ** 2)
    while low(hi) < target:
        hi = span + 2.0 * (hi - span)
    # a long span stretches more than 0.1 % under its weight: widen the lower end until it hangs higher
    for _ in range(40):
        if low(lo) <= target:
            break
        lo = span - 2.0 * (span - lo)
    else:
        raise SystemExit("rest_length: no unstretched length above %g m hangs higher than a sag of %g m" % (lo, target))
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if low(mid) < target:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1.0e-10 * span:
            break
    return 0.5 * (lo + hi)


def catenary(span, L0, w, f, EA, segments=0):
    """End tension (N), sag in the loaded plane (m) and swing angle (deg) under weight w and a
    uniform normal drag f (N/m)."""
    q = math.hypot(w, f)
    H, V0 = bw.solve_cable((span, 0.0), L0, q, EA, segments)
    Ta, _ = bw.cable_end_tensions(H, V0, L0, q, segments)
    return Ta, bw.cable_low_point(H, V0, L0, q, EA, segments), math.degrees(math.atan2(f, w))


def chain(span, L0, w, U, EA, segments, cdrag=0.5 * RHO * CD * DIAMETER, tol=1.0e-9, max_iter=50,
          moordyn=False):
    """End tension (N), sag of the middle node below the chord (m) and its swing angle (deg) of
    the lumped chain at rest in the uniform wind U (m/s) along +x; the span runs along +y and
    gravity along -z. Newton's method on the inner nodes' positions with a finite-difference
    Jacobian.

    moordyn=True takes MoorDyn's conventions (MoorDyn-C 2.7.1 source/Line.cpp, Line::getStateDeriv):
    the drag acts on each link's unstretched length, and the weight is less the buoyancy of the
    displaced air, RHO g pi DIAMETER^2 / 4 per metre. Needs numpy (the rest of this file does not)."""
    import numpy as np
    n = segments
    l0 = L0 / n
    if moordyn:
        w = w - RHO * G * math.pi * DIAMETER ** 2 / 4.0
    W = np.array([U, 0.0, 0.0])
    # start from the catenary in the tilted plane, sampled at equal unstretched arc length
    f0 = cdrag * U * U
    q = math.hypot(w, f0)
    H, V0 = bw.solve_cable((span, 0.0), L0, q, EA, n)
    P = np.zeros((n + 1, 3))
    zl = 0.0
    yl = 0.0
    for i in range(n):
        V = V0 + q * l0 * i
        T = math.hypot(H, V)
        yl += H * l0 * (1.0 / T + 1.0 / EA)
        zl += V * l0 * (1.0 / T + 1.0 / EA)
        P[i + 1] = [0.0, yl, zl]          # zl is the depth below the chord along the load
    P[n, 1] = span
    P[n, 2] = 0.0
    c, s = w / q, f0 / q                  # the load points along (s, 0, -c)
    depth = -P[:, 2].copy()
    P[:, 0] = s * depth
    P[:, 2] = -c * depth

    def forces(P):
        F = np.zeros_like(P)
        d = P[1:] - P[:-1]
        l = np.linalg.norm(d, axis=1)
        t = d / l[:, None]
        T = np.maximum(EA * (l - l0) / l0, 0.0)
        F[:-1] += T[:, None] * t
        F[1:] -= T[:, None] * t
        un = W[None, :] - (t @ W)[:, None] * t
        drag = cdrag * np.linalg.norm(un, axis=1)[:, None] * un * (np.full_like(l, l0) if moordyn else l)[:, None]
        F[:-1] += 0.5 * drag
        F[1:] += 0.5 * drag
        F[:, 2] -= w * l0 * np.r_[0.5, np.ones(n - 1), 0.5]
        return F, T

    def residual(x):
        Q = P.copy()
        Q[1:-1] = x.reshape(-1, 3)
        return forces(Q)[0][1:-1].ravel()

    x = P[1:-1].ravel().copy()
    m = x.size
    for it in range(max_iter):
        r = residual(x)
        Q = P.copy()
        Q[1:-1] = x.reshape(-1, 3)
        scale = max(np.max(forces(Q)[1]), w * l0)      # the largest link tension (N)
        if np.max(np.abs(r)) < tol * scale:
            break
        # banded Jacobian by central differences, three node colours at a time
        J = np.zeros((m, m))
        h = 1.0e-6 * span
        for colour in range(3):
            for k in range(3):
                e = np.zeros(m)
                idx = [3 * j + k for j in range(colour, n - 1, 3)]
                e[idx] = h
                dr = (residual(x + e) - residual(x - e)) / (2.0 * h)
                for j in range(colour, n - 1, 3):
                    rows = slice(max(0, 3 * (j - 1)), min(m, 3 * (j + 2)))
                    J[rows, 3 * j + k] = dr[rows]
        dx = np.linalg.solve(J, -r)
        step = 1.0
        r0 = np.linalg.norm(r)
        while step > 1.0e-4 and np.linalg.norm(residual(x + step * dx)) > (1.0 - 1.0e-4 * step) * r0:
            step *= 0.5
        x = x + step * dx
    else:
        raise SystemExit("static_reference: the chain did not come to rest in %d Newton steps" % max_iter)
    P[1:-1] = x.reshape(-1, 3)
    _, T = forces(P)
    mid = P[n // 2] if n % 2 == 0 else 0.5 * (P[n // 2] + P[n // 2 + 1])
    return T[0], math.hypot(mid[0], mid[2]), math.degrees(math.atan2(mid[0], -mid[2]))


def log_wind(v10, z, z0):
    """The log-law wind at height z given the wind at 10 m, in the form the case's inflow uses."""
    return v10 * math.log((z + z0) / z0) / math.log((10.0 + z0) / z0)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--span", type=float, default=300.0, help="span (m)")
    ap.add_argument("--wind", type=float, nargs="+", default=[0.0, 20.0, 40.0, 60.0],
                    help="uniform wind normal to the span at the conductor (m/s)")
    ap.add_argument("--EA", type=float, default=EA_DEFAULT, help="axial stiffness (N)")
    ap.add_argument("--weight", type=float, default=DRAKE_WEIGHT, help="weight per unstretched metre (N/m)")
    ap.add_argument("--segments", type=int, default=80, help="links of the chain (Ma 2020 use at least 80 elements)")
    a = ap.parse_args(argv)
    if not (a.span > 0.0 and a.EA > 0.0 and a.weight > 0.0 and a.segments >= 2):
        ap.error("--span, --EA and --weight must be > 0 and --segments >= 2")
    if any(u < 0.0 for u in a.wind):
        ap.error("--wind: every speed must be >= 0")
    L0 = rest_length(a.span, SAG_RATIO, a.weight, a.EA, a.segments)
    print("span %g m, unstretched length %.6f m (sag %g m at rest), w %g N/m, EA %g N, %d links"
          % (a.span, L0, SAG_RATIO * a.span, a.weight, a.EA, a.segments))
    print("%8s %10s %12s %12s %8s %10s %10s %8s %12s" % ("U", "f(N/m)", "T_cat(N)", "T_chain(N)", "diff%",
                                                       "sag_cat", "sag_chain", "swing", "T_moordyn(N)"))
    for U in a.wind:
        f = 0.5 * RHO * CD * DIAMETER * U * U
        Tc, sc, sw = catenary(a.span, L0, a.weight, f, a.EA, a.segments)
        Tn, sn, swn = chain(a.span, L0, a.weight, U, a.EA, a.segments)
        Tm, _, _ = chain(a.span, L0, a.weight, U, a.EA, a.segments, moordyn=True)
        print("%8.2f %10.3f %12.1f %12.1f %8.3f %10.3f %10.3f %8.2f %12.1f"
              % (U, f, Tc, Tn, 100.0 * (Tn - Tc) / Tc, sc, sn, swn, Tm))
    return 0


if __name__ == "__main__":
    sys.exit(main())
