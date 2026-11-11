#!/usr/bin/env python3
"""Rest state of a conductor section after one span has let go of a support.

    python3 broken_wire.py plan.json [--json]

A section is a row of spans between two dead ends. At the towers between them the conductor
hangs from insulator strings (or is clamped). The tool finds the section's state at rest in still
air, whole or with one span released at one end, from the balance of forces alone: every span
and every string is a heavy elastic cable, and every string's lower end is a point at which the
forces of the cables meeting there add up to zero. It shares no code with ERF or MoorDyn and
serves as the independent reference for the model's state at rest after a release
(conductors.failure_break). See "Broken-wire state at rest" in
Docs/sphinx_doc/theory/Conductors.rst for the equations.

Python standard library only. The section lies in one vertical plane: x along the line, z up.
"""
import argparse
import json
import math
import sys

GRAVITY = 9.81            # m/s^2, the value ERF's conductor model is given in its tests and decks
STRING_STIFFNESS = 1.0e7  # N, ERF's LineInputs::insulator_axial_stiffness
STRING_SEGMENTS = 2       # ERF's LineInputs::insulator_segments


class PlanError(Exception):
    """A plan the tool refuses; the message names the key."""


class OutsideTool(Exception):
    """A state at rest the tool does not describe: a conductor still held at both ends, or one
    that has lost its string, would lie on the ground."""


# ---------------------------------------------------------------------------------------------
# One cable
# ---------------------------------------------------------------------------------------------
def cable_reach(H, V0, L, w, EA, segments=0):
    """End b's position relative to end a, (dx, dz), of a cable of unstretched length L (m),
    weight w (N per unstretched metre) and axial stiffness EA (N), whose tension at end a is the
    vector (H, V0) (N), pointing from a into the cable.

    segments = 0: the continuous cable. The piece between a and the unstretched arc length s is
    held by its two end tensions and its weight w s, so the tension vector at s is
    (H, V0 + w s); the cable lies along it and is stretched by the factor 1 + T / EA. Integrating
    dx/ds = H (1/T + 1/EA) and dz/ds = (V0 + w s) (1/T + 1/EA) over 0 <= s <= L gives the
    expressions below.

    segments = n > 0: the lumped chain MoorDyn integrates. n straight elastic links of
    unstretched length l = L / n; every inner node carries the weight w l, each end node w l / 2.
    Link i (1..n) carries the tension vector (H, V0 + w l (i - 1)).
    """
    if segments > 0:
        l = L / segments
        dx = dz = 0.0
        for i in range(segments):
            V = V0 + w * l * i
            T = math.hypot(H, V)
            f = l * (1.0 / T + 1.0 / EA)
            dx += H * f
            dz += V * f
        return dx, dz
    V1 = V0 + w * L
    T0 = math.hypot(H, V0)
    T1 = math.hypot(H, V1)
    a = abs(H)
    if a > 1.0e-12 * max(T0, T1):
        dx = H * L / EA + (H / w) * (math.asinh(V1 / a) - math.asinh(V0 / a))
    else:
        dx = 0.0
    dz = (V0 * L + 0.5 * w * L * L) / EA + (T1 - T0) / w
    return dx, dz


def cable_end_pulls(H, V0, L, w, segments=0):
    """The cable's pull (Fx, Fz) on what holds end a, and on what holds end b."""
    if segments > 0:
        half = 0.5 * w * L / segments
        return (H, V0 - half), (-H, -(V0 + w * L - half))
    return (H, V0), (-H, -(V0 + w * L))


def cable_end_tensions(H, V0, L, w, segments=0):
    """The tension at end a and at end b (N): of the continuous cable there, or of the lumped
    chain's first and last link."""
    if segments > 0:
        return math.hypot(H, V0), math.hypot(H, V0 + w * L * (segments - 1) / segments)
    return math.hypot(H, V0), math.hypot(H, V0 + w * L)


def cable_low_point(H, V0, L, w, EA, segments=0):
    """The lowest point of the cable below end a (m, >= 0)."""
    if segments > 0:
        l = L / segments
        z = low = 0.0
        for i in range(segments):
            V = V0 + w * l * i
            T = math.hypot(H, V)
            z += V * l * (1.0 / T + 1.0 / EA)
            low = min(low, z)
        return -low
    if V0 >= 0.0:
        return 0.0
    s = min(L, -V0 / w)   # where the tension is horizontal, or end b
    T0 = math.hypot(H, V0)
    Ts = math.hypot(H, V0 + w * s)
    return -((V0 * s + 0.5 * w * s * s) / EA + (Ts - T0) / w)


# ---------------------------------------------------------------------------------------------
# Linear algebra and Newton's method (standard library)
# ---------------------------------------------------------------------------------------------
def solve_linear(A, b):
    """Gaussian elimination with partial pivoting; A is a list of rows and is overwritten."""
    n = len(b)
    x = list(b)
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(A[r][c]))
        if A[p][c] == 0.0:
            raise ArithmeticError("singular system")
        if p != c:
            A[c], A[p] = A[p], A[c]
            x[c], x[p] = x[p], x[c]
        for r in range(c + 1, n):
            m = A[r][c] / A[c][c]
            if m != 0.0:
                for k in range(c, n):
                    A[r][k] -= m * A[c][k]
                x[r] -= m * x[c]
    for c in range(n - 1, -1, -1):
        s = x[c]
        for k in range(c + 1, n):
            s -= A[c][k] * x[k]
        x[c] = s / A[c][c]
    return x


def newton(residual, x0, steps, tol=1.0e-11, max_iter=200):
    """Solve residual(x) = 0 by Newton's method with a central-difference Jacobian and step
    halving. residual returns scaled, dimensionless values; steps are the difference steps of
    the unknowns. Raises ArithmeticError when it does not converge."""
    x = list(x0)
    r = residual(x)
    norm = max(abs(v) for v in r)
    n = len(x)
    for _ in range(max_iter):
        if norm < tol:
            return x
        J = [[0.0] * n for _ in range(n)]
        for j in range(n):
            xp = list(x); xp[j] += steps[j]
            xm = list(x); xm[j] -= steps[j]
            rp = residual(xp)
            rm = residual(xm)
            for i in range(n):
                J[i][j] = (rp[i] - rm[i]) / (2.0 * steps[j])
        dx = solve_linear(J, [-v for v in r])
        lam = 1.0
        while True:
            xn = [x[i] + lam * dx[i] for i in range(n)]
            try:
                rn = residual(xn)
                nn = max(abs(v) for v in rn)
            except (ValueError, ZeroDivisionError, OverflowError, ArithmeticError):
                nn = math.inf
            if nn < norm or lam < 1.0e-10:
                break
            lam *= 0.5
        if not nn < norm:
            if norm < 1.0e3 * tol:
                return x      # at the roundoff of the residual
            raise ArithmeticError("no descent: residual %.3e" % norm)
        x, r, norm = xn, rn, nn
    raise ArithmeticError("not converged: residual %.3e" % norm)


def solve_cable(reach, L, w, EA, segments=0):
    """The tension vector (H, V0) at end a of a cable whose end b is at reach = (dx, dz) from a."""
    dx, dz = reach
    chord = math.hypot(dx, dz)
    h = max(abs(dx), 1.0e-9 * L)
    sign = 1.0 if dx >= 0.0 else -1.0
    if L > chord:
        sag = math.sqrt(3.0 * chord * (L - chord) / 8.0)
        H = w * h * h / (8.0 * max(sag, 1.0e-6 * L))
    else:
        H = EA * (chord - L) / L * h / chord + w * h * h / (8.0 * 0.01 * L)
    H = max(H, 1.0e-3 * w * L)
    V0 = -0.5 * w * L + H * dz / h
    scale = max(L, 1.0)

    def residual(x):
        rx, rz = cable_reach(sign * x[0], x[1], L, w, EA, segments)
        return [(rx - dx) / scale, (rz - dz) / scale]

    F = max(w * L, H)
    x = newton(residual, [H, V0], [1.0e-6 * F, 1.0e-6 * F])
    return sign * x[0], x[1]


# ---------------------------------------------------------------------------------------------
# The plan
# ---------------------------------------------------------------------------------------------
def _number(d, key, where, lo=None, hi=None, lo_open=False, default=None):
    if key not in d:
        if default is None:
            raise PlanError("%s%s is required" % (where, key))
        return default
    v = d[key]
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        raise PlanError("%s%s must be a finite number" % (where, key))
    if lo is not None and (v < lo or (lo_open and v == lo)):
        raise PlanError("%s%s must be %s %g" % (where, key, ">" if lo_open else ">=", lo))
    if hi is not None and v > hi:
        raise PlanError("%s%s must be <= %g" % (where, key, hi))
    return float(v)


def _point(d, key, where):
    v = d.get(key)
    if (not isinstance(v, list) or len(v) != 2
            or any(isinstance(c, bool) or not isinstance(c, (int, float)) or not math.isfinite(c) for c in v)):
        raise PlanError("%s%s must be [x, z], two finite numbers (m)" % (where, key))
    return float(v[0]), float(v[1])


PLAN_KEYS = {"end_a", "end_b", "towers", "lengths", "diameter", "mass_per_length", "axial_stiffness",
             "insulator_length", "insulator_mass", "insulator_diameter", "air_density", "gravity",
             "ground", "segments", "released"}


def read_plan(plan):
    """Check a plan (a dict) and return it with defaults filled in. Raises PlanError."""
    if not isinstance(plan, dict):
        raise PlanError("the plan must be a JSON object")
    unknown = sorted(set(plan) - PLAN_KEYS)
    if unknown:
        raise PlanError("unknown key %s" % unknown[0])
    p = {}
    p["end_a"] = _point(plan, "end_a", "")
    p["end_b"] = _point(plan, "end_b", "")
    towers = plan.get("towers", [])
    if not isinstance(towers, list):
        raise PlanError("towers must be a list of [x, z]")
    p["towers"] = [_point({"t": t}, "t", "towers: ") for t in towers]
    xs = [p["end_a"][0]] + [t[0] for t in p["towers"]] + [p["end_b"][0]]
    if any(b <= a for a, b in zip(xs, xs[1:])):
        raise PlanError("end_a, towers and end_b must follow one another with increasing x")
    lengths = plan.get("lengths")
    if not isinstance(lengths, list) or len(lengths) != len(xs) - 1:
        raise PlanError("lengths must hold one unstretched length (m) per span: %d" % (len(xs) - 1))
    p["lengths"] = [_number({"l": l}, "l", "lengths: ", 0.0, lo_open=True) for l in lengths]
    p["diameter"] = _number(plan, "diameter", "", 0.0, 1.0, lo_open=True)
    p["mass_per_length"] = _number(plan, "mass_per_length", "", 0.0, lo_open=True)
    p["axial_stiffness"] = _number(plan, "axial_stiffness", "", 0.0, lo_open=True)
    p["insulator_length"] = _number(plan, "insulator_length", "", 0.0, default=0.0)
    if p["insulator_length"] > 0.0:
        p["insulator_mass"] = _number(plan, "insulator_mass", "", 0.0, lo_open=True)
        p["insulator_diameter"] = _number(plan, "insulator_diameter", "", 0.0, lo_open=True, default=0.254)
    else:
        for key in ("insulator_mass", "insulator_diameter"):
            if key in plan:
                raise PlanError("%s is given but insulator_length is 0: the conductor is clamped" % key)
    p["air_density"] = _number(plan, "air_density", "", 0.0, default=1.2)
    p["gravity"] = _number(plan, "gravity", "", 0.0, lo_open=True, default=GRAVITY)
    p["ground"] = _number(plan, "ground", "", default=0.0)
    seg = plan.get("segments", 0)
    if isinstance(seg, bool) or not isinstance(seg, int) or seg < 0 or seg > 100000:
        raise PlanError("segments must be an integer in [0, 100000]; 0 is the continuous cable")
    p["segments"] = seg
    nspans = len(p["lengths"])
    rel = plan.get("released")
    if rel is not None:
        if isinstance(rel, dict) and set(rel) == {"string"}:
            j = rel["string"]
            if p["insulator_length"] <= 0.0:
                raise PlanError("released: string needs insulator_length > 0")
            if isinstance(j, bool) or not isinstance(j, int) or j < 1 or j > nspans - 1:
                raise PlanError("released: string must be an integer in [1, %d]" % (nspans - 1))
            rel = None
            p["released"] = ("string", j)
    if rel is None:
        p.setdefault("released", None)
    else:
        if not isinstance(rel, dict) or set(rel) != {"span", "end"}:
            raise PlanError("released must be {\"span\": k, \"end\": \"a\" or \"b\"} or {\"string\": j}")
        k = rel["span"]
        if isinstance(k, bool) or not isinstance(k, int) or k < 1 or k > nspans:
            raise PlanError("released: span must be an integer in [1, %d]" % nspans)
        if rel["end"] not in ("a", "b"):
            raise PlanError("released: end must be \"a\" or \"b\"")
        p["released"] = (k, rel["end"])
    w = (p["mass_per_length"] - p["air_density"] * 0.25 * math.pi * p["diameter"] ** 2) * p["gravity"]
    if not w > 0.0:
        raise PlanError("air_density: the conductor must be heavier than the air it displaces")
    p["weight"] = w
    if p["insulator_length"] > 0.0:
        wi = (p["insulator_mass"] / p["insulator_length"]
              - p["air_density"] * 0.25 * math.pi * p["insulator_diameter"] ** 2) * p["gravity"]
        if not wi > 0.0:
            raise PlanError("air_density: the string must be heavier than the air it displaces")
        p["insulator_weight"] = wi
    lowest = min([p["end_a"][1], p["end_b"][1]] + [t[1] - p["insulator_length"] for t in p["towers"]])
    if not lowest > p["ground"] + 0.5 * p["diameter"]:
        raise PlanError("ground: every attachment must be above the ground")
    return p


# ---------------------------------------------------------------------------------------------
# The section at rest
# ---------------------------------------------------------------------------------------------
def hanging_load(height, L, w, EA, segments=0):
    """A span of unstretched length L that hangs from one support, 'height' (m) above the level
    its axis rests at on the ground. Returns the pull it puts on the support (N), the
    unstretched length that hangs (m) and the tension at the held end (N). The ground has no
    friction, so the pull is vertical.

    Continuous: the hanging length l stretches under its own weight to
    l (1 + w l / (2 EA)) = height, unless the whole span is shorter; pull and tension are w l.
    Lumped: the n nodes below the held end that hang clear of the ground weigh w L / segments
    each (the free end's node half of it); the first link's tension is their weight (the docs'
    Eq. br-lumped) and the support also carries the held end's own node, half a link's weight.
    The node at the foot is taken as resting on the ground.
    """
    if segments > 0:
        l = L / segments
        depth = 0.0
        clear = 0
        for i in range(1, segments + 1):
            # link i, counted from the held end, stretched by about the weight below its middle
            depth += l * (1.0 + w * l * (i - 0.5) / EA)
            if depth < height:
                clear = i
            else:
                break
        if clear == segments:
            return w * L, L, w * L - 0.5 * w * l
        return w * l * clear + 0.5 * w * l, l * clear, w * l * clear
    if L * (1.0 + 0.5 * w * L / EA) <= height:
        return w * L, L, w * L
    # l^2 w / (2 EA) + l - height = 0
    a = 0.5 * w / EA
    l = 2.0 * height / (1.0 + math.sqrt(1.0 + 4.0 * a * height))
    return w * l, l, w * l


def solve_section(plan):
    """The section's state at rest. plan: see README.md. Returns a dict of results."""
    p = read_plan(plan)
    w, EA, seg = p["weight"], p["axial_stiffness"], p["segments"]
    nspans = len(p["lengths"])
    strings = p["insulator_length"] > 0.0
    Li = p["insulator_length"]
    wi = p.get("insulator_weight", 0.0)
    sseg = STRING_SEGMENTS if seg > 0 else 0
    rest_level = p["ground"] + 0.5 * p["diameter"]

    # attachment a = 0 .. nspans: where the conductor is held. With strings the inner ones are
    # the strings' lower ends, free to move; their first guess is plumb below the cross-arm
    tops = [p["end_a"]] + p["towers"] + [p["end_b"]]
    free = [strings and 0 < a < nspans for a in range(nspans + 1)]
    pos = []
    for a, (x, z) in enumerate(tops):
        if free[a]:
            # a string under half of each span beside it, stretched by that load
            load = 0.5 * w * (p["lengths"][a - 1] + p["lengths"][a]) + wi * Li
            pos.append((x, z - Li * (1.0 + load / STRING_STIFFNESS)))
        else:
            pos.append((x, z))

    released = p["released"]
    attached = [True] * nspans            # span k (0-based) is attached at both ends
    held_at = None                        # the attachment the released span still hangs from
    dropped = None                        # the tower whose string has left its cross-arm
    if released is not None and released[0] == "string":
        dropped = released[1]
    elif released is not None:
        k, end = released
        attached[k - 1] = False
        held_at = k if end == "a" else k - 1
    unknown = [a for a in range(nspans + 1) if free[a]]
    index = {a: 2 * i for i, a in enumerate(unknown)}
    scale_f = w * max(p["lengths"])

    def state(x, lam=1.0, kept=None):
        """Positions, the cables' end vectors and the residual forces on the free points.
        lam < 1 is a step of the continuation from the whole section (lam = 0) to the released
        one (lam = 1): the element that lets go keeps the fraction 1 - lam of the pull 'kept' it
        had in the whole section, and what replaces it acts with the fraction lam."""
        q = list(pos)
        for a in unknown:
            q[a] = (x[index[a]], x[index[a] + 1])
        spans = [None] * nspans
        force = {a: [0.0, 0.0] for a in unknown}
        pulls = [[0.0, 0.0] for _ in range(nspans + 1)]     # the conductor's pull on attachment a
        for k in range(nspans):
            if not attached[k]:
                if lam < 1.0:
                    for a in (k, k + 1):
                        pulls[a][0] += (1.0 - lam) * kept[a][0]
                        pulls[a][1] += (1.0 - lam) * kept[a][1]
                continue
            reach = (q[k + 1][0] - q[k][0], q[k + 1][1] - q[k][1])
            H, V0 = solve_cable(reach, p["lengths"][k], w, EA, seg)
            spans[k] = (H, V0)
            fa, fb = cable_end_pulls(H, V0, p["lengths"][k], w, seg)
            pulls[k][0] += fa[0]; pulls[k][1] += fa[1]
            pulls[k + 1][0] += fb[0]; pulls[k + 1][1] += fb[1]
        hang = None
        if held_at is not None:
            k = released[0] - 1
            load, length, tension = hanging_load(q[held_at][1] - rest_level, p["lengths"][k], w, EA, seg)
            pulls[held_at][1] -= lam * load
            hang = (load, length, tension)
        strs = {}
        for a in unknown:
            if a == dropped:
                # the string hangs from the conductor: its weight, and no support
                strs[a] = None
                force[a][0] = pulls[a][0]
                force[a][1] = pulls[a][1] - lam * wi * Li
                if lam < 1.0:
                    force[a][0] += (1.0 - lam) * kept[a][0]
                    force[a][1] += (1.0 - lam) * kept[a][1]
                continue
            top = tops[a]
            reach = (q[a][0] - top[0], q[a][1] - top[1])
            # a string is close to plumb: solve for its end vector from the vertical guess
            Hs, Vs = solve_string(reach, Li, wi, sseg)
            strs[a] = (Hs, Vs)
            _, fb = cable_end_pulls(Hs, Vs, Li, wi, sseg)
            force[a][0] = pulls[a][0] + fb[0]
            force[a][1] = pulls[a][1] + fb[1]
        return q, spans, strs, pulls, hang, force

    def solve(x0, lam, kept, kappa=0.0):
        def residual(x):
            force = state(x, lam, kept)[5]
            out = []
            for a in unknown:
                # kappa: a spring (N/m) that holds the foot under its cross-arm, for solve_whole
                out += [(force[a][0] - kappa * (x[index[a]] - pos[a][0])) / scale_f, force[a][1] / scale_f]
            return out
        return newton(residual, x0, [1.0e-6] * len(x0), tol=1.0e-10)

    def solve_whole(x0):
        """The whole section from plumb strings. Where the spans pull unequally the strings lean
        far and Newton's method from plumb may find no descent: the feet are then held under
        their cross-arms by springs along the line, which are slackened in steps to nothing."""
        try:
            return solve(x0, 1.0, None)
        except ArithmeticError:
            pass
        kappa = 100.0 * scale_f / max(Li, 1.0)
        x = solve(x0, 1.0, None, kappa)
        factor = 0.25
        floor = 1.0e-9 * kappa
        while kappa > 0.0:
            trial = kappa * factor if kappa * factor > floor else 0.0
            try:
                x = solve(x, 1.0, None, trial)
                kappa = trial
            except ArithmeticError:
                factor = math.sqrt(factor)
                if factor > 0.999:
                    raise
        return x

    x = []
    for a in unknown:
        x += [pos[a][0], pos[a][1]]
    if unknown and released is None:
        x = solve_whole(x)
    elif unknown:
        # the whole section first, then the release in steps: each step starts from the last
        # one's answer; a step that finds no answer is halved
        rel_attached, rel_dropped = list(attached), dropped
        attached = [True] * nspans
        dropped = None
        x = solve_whole(x)
        _, spans0, strs0, _, _, _ = state(x)
        kept = {}
        if rel_dropped is not None:
            Hs, Vs = strs0[rel_dropped]
            kept[rel_dropped] = cable_end_pulls(Hs, Vs, Li, wi, sseg)[1]
        else:
            k = released[0] - 1
            H, V0 = spans0[k]
            kept[k], kept[k + 1] = cable_end_pulls(H, V0, p["lengths"][k], w, seg)
        attached, dropped = rel_attached, rel_dropped
        lam, step = 0.0, 0.25
        while lam < 1.0:
            trial = min(1.0, lam + step)
            try:
                x = solve(x, trial, kept)
                lam = trial
                step = min(0.25, 2.0 * step)
            except ArithmeticError:
                step *= 0.5
                if step < 1.0e-6:
                    raise
    q, spans, strs, pulls, hang, _ = state(x)

    out = {"weight": w, "segments": seg, "spans": [], "strings": [], "attachments": [list(c) for c in q]}
    for k in range(nspans):
        L = p["lengths"][k]
        if spans[k] is None:
            load, length, tension = hang
            out["spans"].append({"span": k + 1, "attached": released[1] == "a" and "b" or "a",
                                 "hanging_length": length, "tension_held_end": tension,
                                 "pull_on_support": load, "on_ground": L - length})
            continue
        H, V0 = spans[k]
        if q[k][1] - cable_low_point(H, V0, L, w, EA, seg) < rest_level:
            raise OutsideTool("span %d, held at both ends, would lie on the ground (its lowest point %.3f m "
                              "below the level the conductor rests at)"
                              % (k + 1, rest_level - q[k][1] + cable_low_point(H, V0, L, w, EA, seg)))
        Ta, Tb = cable_end_tensions(H, V0, L, w, seg)
        out["spans"].append({"span": k + 1, "attached": "both", "horizontal_tension": H,
                             "tension_a": Ta, "tension_b": Tb,
                             "low_point_below_a": cable_low_point(H, V0, L, w, EA, seg),
                             "lowest_height": q[k][1] - cable_low_point(H, V0, L, w, EA, seg) - p["ground"]})
    for a in range(1, nspans):
        if a == dropped:
            out["strings"].append({"tower": a, "released": True, "angle_deg": None, "tension_top": 0.0,
                                   "pull_on_tower": [0.0, 0.0], "foot": [q[a][0], q[a][1]]})
        elif free[a]:
            Hs, Vs = strs[a]
            fa, _ = cable_end_pulls(Hs, Vs, Li, wi, sseg)
            dx, dz = q[a][0] - tops[a][0], q[a][1] - tops[a][1]
            out["strings"].append({"tower": a, "angle_deg": math.degrees(math.atan2(dx, -dz)),
                                   "tension_top": math.hypot(Hs, Vs),
                                   "pull_on_tower": [fa[0], fa[1]], "foot": [q[a][0], q[a][1]]})
        else:
            out["strings"].append({"tower": a, "angle_deg": 0.0, "tension_top": math.hypot(*pulls[a]),
                                   "pull_on_tower": list(pulls[a]), "foot": [q[a][0], q[a][1]]})
    out["pull_on_end_a"] = list(pulls[0])
    out["pull_on_end_b"] = list(pulls[nspans])
    return out


def solve_string(reach, L, w, segments):
    """A string's end vector at its top for a lower end at 'reach' from the top. The string is
    short, stiff and close to straight, so the first guess is the straight link."""
    dx, dz = reach
    dist = math.hypot(dx, dz)
    EA = STRING_STIFFNESS
    T = max(EA * (dist - L) / L, 0.0) + 0.5 * w * L + 1.0
    H, V0 = T * dx / dist, T * dz / dist - 0.5 * w * L
    scale = L

    def residual(x):
        rx, rz = cable_reach(x[0], x[1], L, w, EA, segments)
        return [(rx - dx) / scale, (rz - dz) / scale]

    F = max(T, w * L)
    x = newton(residual, [H, V0], [1.0e-7 * F, 1.0e-7 * F], tol=1.0e-12)
    return x[0], x[1]


# ---------------------------------------------------------------------------------------------
# Command
# ---------------------------------------------------------------------------------------------
def report(out):
    lines = ["conductor weight %.6f N/m; %s" % (out["weight"], "continuous cables" if out["segments"] == 0
                                                 else "lumped chains, %d links per span" % out["segments"])]
    lines.append("span  attached  H (N)        T_a (N)      T_b (N)      lowest point above ground (m)")
    for s in out["spans"]:
        if s["attached"] == "both":
            lines.append("%-5d both      %-12.3f %-12.3f %-12.3f %.4f" % (
                s["span"], s["horizontal_tension"], s["tension_a"], s["tension_b"], s["lowest_height"]))
        else:
            lines.append("%-5d end %s     hangs %.3f m, %.3f m on the ground; tension at the held end %.3f N" % (
                s["span"], s["attached"], s["hanging_length"], s["on_ground"], s["tension_held_end"]))
    lines.append("tower angle (deg)  T_top (N)    pull on the tower: along the line, vertical (N)")
    for t in out["strings"]:
        if t.get("released"):
            lines.append("%-5d released; the conductor there is at x = %.3f m, z = %.3f m" % (
                t["tower"], t["foot"][0], t["foot"][1]))
            continue
        lines.append("%-5d %-12.4f %-12.3f %-12.3f %.3f" % (
            t["tower"], t["angle_deg"], t["tension_top"], t["pull_on_tower"][0], t["pull_on_tower"][1]))
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Rest state of a conductor section, whole or after a release.")
    ap.add_argument("plan", help="JSON plan (see README.md)")
    ap.add_argument("--json", action="store_true", help="print the results as JSON")
    args = ap.parse_args(argv)
    try:
        with open(args.plan) as f:
            plan = json.load(f)
        out = solve_section(plan)
    except (OSError, ValueError) as e:
        print("broken_wire: %s" % e, file=sys.stderr)
        return 2
    except PlanError as e:
        print("broken_wire: %s" % e, file=sys.stderr)
        return 2
    except ArithmeticError as e:
        print("broken_wire: no state at rest found (%s)" % e, file=sys.stderr)
        return 3
    except OutsideTool as e:
        print("broken_wire: %s: the tool describes conductors that hang clear of the ground" % e, file=sys.stderr)
        return 4
    print(json.dumps(out, indent=1) if args.json else report(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
