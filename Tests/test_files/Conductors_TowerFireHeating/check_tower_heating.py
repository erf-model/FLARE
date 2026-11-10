#!/usr/bin/env python3
"""Check a Conductors_TowerFireHeating run against an independent calculation.

The script recomputes, from the deck's numbers and the frame the run wrote (conductors/frame_<tower>.dat, the
SubDyn file, and frame_<tower>_members.dat, the members' angles), with no code shared with ERF:

1. each member's three Gauss points and axis in ERF's axes: the tower stands on its ground.dat point, its
   frame's x axis along the line and its y axis along the cross-arm (normal to the line, to its left);
2. the air temperature at each point, T = theta - g z / c_p for the still, isentropic sounding (theta 300 K),
   and the fire's irradiance there from the fire cells inside the prescribed heat disc (point sources of
   chi_r q'' dA at the source height within the cutoff, projected on the member's axis), both as Gauss means;
3. each member's temperature from the air temperature at the first step, stepped with backward Euler in the
   steel's enthalpy (EN 1993-1-2 section 3.4.1.2, integrated here numerically) against convection, radiation to the
   first step's air and the absorbed fire radiation, for an equal-leg angle b x t: area t (2b - t), exchange
   perimeter 4b, intercepted width ((2 + sqrt 2) b + (2 - sqrt 2) t) / pi (the hull of the outline);
4. when the frame is rebuilt: whenever a member's EN 1993-1-2 k_y or k_E (Table 3.1) has moved by more than the
   tolerance, relative to the larger of the two, from the frame's, which then takes the members' temperatures.

Every member's temperature in tower_<tower>_member_temperature.dat and every row of tower_temperature.dat
(mean, largest, its member, largest irradiance, least k_E, rebuilds, the frame's largest temperature) must
match. Pure Python (no numpy). Exit status 0 when every check passes.
"""
import argparse
import bisect
import math
import sys

SIGMA = 5.670374419e-8
RHO_STEEL = 7850.0


def c_steel(th):
    """EN 1993-1-2 section 3.4.1.2 [J/(kg K)], the 20 C value below 20 C."""
    th = max(th, 20.0)
    if th < 600.0:
        return 425.0 + 0.773 * th - 1.69e-3 * th ** 2 + 2.22e-6 * th ** 3
    if th < 735.0:
        return 666.0 + 13002.0 / (738.0 - th)
    if th < 900.0:
        return 545.0 + 17820.0 / (th - 731.0)
    return 650.0


class Enthalpy:
    """H(th) - H(20) [J/kg] by Simpson's rule on a 0.01 K grid from -50 to 1250 C, split at the formula's breaks."""

    def __init__(self, lo=-50.0, hi=1250.0, h=0.01):
        n = int(round((hi - lo) / h))
        self.x = [lo + i * h for i in range(n + 1)]
        self.H = [0.0] * (n + 1)
        for i in range(n):
            a, b = self.x[i], self.x[i + 1]
            self.H[i + 1] = self.H[i] + (b - a) / 6.0 * (c_steel(a) + 4.0 * c_steel(0.5 * (a + b)) + c_steel(b))
        i20 = int(round((20.0 - lo) / h))
        h20 = self.H[i20]
        self.H = [v - h20 for v in self.H]

    def __call__(self, th):
        i = min(max(bisect.bisect_right(self.x, th) - 1, 0), len(self.x) - 2)
        f = (th - self.x[i]) / (self.x[i + 1] - self.x[i])
        return self.H[i] + f * (self.H[i + 1] - self.H[i])


KY = [(20, 1), (100, 1), (200, 1), (300, 1), (400, 1), (500, .78), (600, .47), (700, .23), (800, .11), (900, .06),
      (1000, .04), (1100, .02), (1200, 0)]
KE = [(20, 1), (100, 1), (200, .9), (300, .8), (400, .7), (500, .6), (600, .31), (700, .13), (800, .09), (900, .0675),
      (1000, .045), (1100, .0225), (1200, 0)]


def table(tab, th):
    if th <= tab[0][0]:
        return tab[0][1]
    for (a, fa), (b, fb) in zip(tab, tab[1:]):
        if th <= b:
            return fa + (fb - fa) * (th - a) / (b - a)
    return 0.0


def read_frame(path):
    """Joints {id: (x, y, z)} and members [(id, joint_a, joint_b)] of a SubDyn file."""
    lines = open(path).read().splitlines()
    joints, members = {}, []
    for i, l in enumerate(lines):
        if "NJoints" in l:
            n = int(l.split()[0])
            for r in lines[i + 3:i + 3 + n]:
                f = r.split()
                joints[int(f[0])] = tuple(float(v) for v in f[1:4])
        if "NMembers" in l:
            n = int(l.split()[0])
            for r in lines[i + 3:i + 3 + n]:
                f = r.split()
                members.append((int(f[0]), int(f[1]), int(f[2])))
    return joints, members


def read_angles(path):
    out = {}
    for l in open(path):
        if l.startswith("#") or not l.strip():
            continue
        f = l.split()
        out[int(f[0])] = (float(f[3]), float(f[4]))
    return out


def fire_cells(a):
    n = int(round(a.extent / a.fire_dx))
    out = []
    for j in range(int(round(a.extent_y / a.fire_dx))):
        y = (j + 0.5) * a.fire_dx
        for i in range(n):
            x = (i + 0.5) * a.fire_dx
            if (x - a.cx) ** 2 + (y - a.cy) ** 2 <= a.radius ** 2:
                out.append((x, y, a.flux * a.fire_dx ** 2))
    return out


def irradiance(a, cells, p, t):
    rmin2 = a.fire_dx ** 2 / math.pi
    s = 0.0
    for (x, y, Q) in cells:
        r = [p[0] - x, p[1] - y, p[2] - a.source_height]
        r2 = sum(v * v for v in r)
        if r2 > a.cutoff ** 2:
            continue
        cx = t[1] * r[2] - t[2] * r[1]
        cy = t[2] * r[0] - t[0] * r[2]
        cz = t[0] * r[1] - t[1] * r[0]
        sin_psi = min(1.0, math.sqrt((cx * cx + cy * cy + cz * cz) / r2))
        s += Q * sin_psi / (4.0 * math.pi * max(r2, rmin2))
    return a.chi * s


def air_temperature(z):
    return 300.0 - 9.81 * z / 1004.5 - 273.15


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dir", default=".")
    p.add_argument("--tower", default="S1_t1")
    p.add_argument("--dt", type=float, required=True)
    p.add_argument("--fire-lags", type=int, required=True, help="1: the first step sees no fire yet")
    p.add_argument("--fire-dx", type=float, required=True)
    p.add_argument("--extent", type=float, default=1500.0)
    p.add_argument("--extent-y", type=float, default=1000.0)
    p.add_argument("--flux", type=float, default=2.0e5)
    p.add_argument("--radius", type=float, default=20.0)
    p.add_argument("--cx", type=float, default=755.0)
    p.add_argument("--cy", type=float, default=510.0)
    p.add_argument("--chi", type=float, default=0.3)
    p.add_argument("--cutoff", type=float, default=300.0)
    p.add_argument("--source-height", type=float, default=2.0)
    p.add_argument("--line", type=float, nargs=4, default=[600.0, 500.0, 900.0, 500.0], help="end_a x y, end_b x y")
    p.add_argument("--eps", type=float, default=0.6)
    p.add_argument("--hc", type=float, default=20.0)
    p.add_argument("--tol", type=float, default=0.02, help="the reduction tolerance")
    p.add_argument("--min-rebuilds", type=int, default=1)
    p.add_argument("--atol", type=float, default=1.0e-3, help="K")
    a = p.parse_args()

    joints, members = read_frame("%s/conductors/frame_%s.dat" % (a.dir, a.tower))
    angles = read_angles("%s/conductors/frame_%s_members.dat" % (a.dir, a.tower))
    base = None
    for l in open(a.dir + "/conductors/ground.dat"):
        f = l.split()
        if len(f) == 6 and f[0] + "_" + f[1] == a.tower:
            base = (float(f[2]), float(f[3]), float(f[4]))
    assert base is not None, "the tower is not in ground.dat"
    # the frame's x axis along the line, y along the cross-arm, to the line's left (z x along)
    dx, dy = a.line[2] - a.line[0], a.line[3] - a.line[1]
    n = math.hypot(dx, dy)
    along, across = (dx / n, dy / n), (-dy / n, dx / n)

    def erf(x):
        return (base[0] + x[0] * along[0] + x[1] * across[0], base[1] + x[0] * along[1] + x[1] * across[1], base[2] + x[2])

    g = 0.5 * math.sqrt(0.6)
    gs, gw = (0.5 - g, 0.5, 0.5 + g), (5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0)
    cells = fire_cells(a)
    sec, Ta, E = [], [], []
    for (mid, ja, jb) in members:
        b, t = angles[mid]
        sec.append((t * (2 * b - t), 4.0 * b, ((2.0 + math.sqrt(2.0)) * b + (2.0 - math.sqrt(2.0)) * t) / math.pi))
        pa, pb = erf(joints[ja]), erf(joints[jb])
        L = math.sqrt(sum((pb[c] - pa[c]) ** 2 for c in range(3)))
        ax = [(pb[c] - pa[c]) / L for c in range(3)]
        pts = [[pa[c] + s * (pb[c] - pa[c]) for c in range(3)] for s in gs]
        Ta.append(sum(w * air_temperature(q[2]) for w, q in zip(gw, pts)))
        E.append(sum(w * irradiance(a, cells, q, ax) for w, q in zip(gw, pts)))
    H = Enthalpy()

    def step(m, th0, e):
        A, P, W = sec[m]
        h0 = H(th0)

        def f(th):
            K = 273.15
            gain = P * (a.hc * (Ta[m] - th) + a.eps * SIGMA * ((Ta[m] + K) ** 4 - (th + K) ** 4)) + a.eps * W * e
            return RHO_STEEL * A * (H(th) - h0) - a.dt * gain
        lo, hi = th0 - 50.0, th0 + 400.0
        assert f(lo) < 0 < f(hi)
        while hi - lo > 1e-9:
            mid = 0.5 * (lo + hi)
            if f(mid) < 0:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    member_rows = [l.split() for l in open("%s/conductors/tower_%s_member_temperature.dat" % (a.dir, a.tower)) if not l.startswith("#")]
    summary = [l.split() for l in open(a.dir + "/conductors/tower_temperature.dat") if not l.startswith("#") and l.split()[1] == a.tower]
    by_time = {round(float(r[0]), 6): r for r in member_rows}
    theta = list(Ta)
    compared = 0
    frame = [20.0] * len(members)
    rebuilds = 0
    bad, worst = [], {"T": 0.0, "E": 0.0}
    nsteps = int(round((float(summary[-1][0]) - float(summary[0][0])) / a.dt)) + 1
    t0 = float(summary[0][0]) - a.dt
    srow = {round(float(r[0]), 6): r for r in summary}
    for k in range(nsteps):
        fire_on = not (k == 0 and a.fire_lags)
        theta = [step(m, theta[m], E[m] if fire_on else 0.0) for m in range(len(members))]
        def rel(p_, q_):
            return abs(p_ - q_) / max(p_, q_) if max(p_, q_) > 0 else 0.0
        change = max(max(rel(table(KY, x), table(KY, y)), rel(table(KE, x), table(KE, y))) for x, y in zip(theta, frame))
        if change > a.tol:
            frame = list(theta)
            rebuilds += 1
        t = round(t0 + (k + 1) * a.dt, 6)
        if t in by_time:
            compared += 1
            run = [float(v) for v in by_time[t][1:]]
            d = max(abs(x - y) for x, y in zip(run, theta))
            worst["T"] = max(worst["T"], d)
            if d > a.atol:
                bad.append("t = %g: a member differs by %.3g K" % (t, d))
                break
        if t in srow:
            r = srow[t]
            hot = max(range(len(members)), key=lambda m: theta[m])
            emax = max(E) if fire_on else 0.0
            kmin = min(table(KE, x) for x in frame)
            worst["E"] = max(worst["E"], abs(float(r[5]) - emax) / max(emax, 1.0))
            checks = [(float(r[2]), sum(theta) / len(theta), a.atol, "mean"), (float(r[3]), theta[hot], a.atol, "largest"),
                      (float(r[6]), kmin, 1e-5, "least k_E"), (float(r[8]), max(frame), a.atol, "frame's largest")]
            for got, want, tol, what in checks:
                if abs(got - want) > tol:
                    bad.append("t = %g: %s %g, the calculation gives %g" % (t, what, got, want))
            if int(r[4]) != members[hot][0]:
                bad.append("t = %g: the hottest member is %s, the calculation gives %d" % (t, r[4], members[hot][0]))
            if int(r[7]) != rebuilds:
                bad.append("t = %g: %s rebuilds, the calculation gives %d" % (t, r[7], rebuilds))
            if abs(float(r[5]) - emax) > 1e-6 * max(emax, 1.0):
                bad.append("t = %g: largest irradiance %s, the calculation gives %g" % (t, r[5], emax))
            if bad:
                break
    hot = max(range(len(members)), key=lambda m: theta[m])
    print("members %d, fire cells %d; the hottest member %d at %.3f C (air %.3f C, irradiance %.1f W/m^2); %d rebuilds"
          % (len(members), len(cells), members[hot][0], theta[hot], Ta[hot], E[hot], rebuilds))
    print("largest differences: member temperature %.3g K, largest irradiance (relative) %.3g" % (worst["T"], worst["E"]))
    if not bad and compared != len(member_rows):
        bad.append("%d of the %d rows of the member log were compared: their times do not match the steps" % (compared, len(member_rows)))
    if rebuilds < a.min_rebuilds:
        bad.append("the frame was rebuilt %d time(s), fewer than %d: the case does not test the rebuild" % (rebuilds, a.min_rebuilds))
    if theta[hot] - Ta[hot] < 50.0:
        bad.append("the hottest member heated by %.1f K only" % (theta[hot] - Ta[hot]))
    for b in bad:
        print("FAIL: " + b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
