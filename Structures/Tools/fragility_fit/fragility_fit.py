#!/usr/bin/env python3
"""Fragility curves from an ensemble of conductor runs (standard library only).

Reads the table an ensemble gathers (demand_table_all.csv of
Structures/Tools/fragility_ensemble: one row per member and component with the wind recorded
near the component and its peak demand-to-capacity ratio) and fits, for every component and for
every group of components, the lognormal fragility function

    P(exceedance | V = x) = Phi(ln(x / theta) / beta)

by maximum likelihood on the members' outcomes (Baker 2015, Eqs. 1 and 8 to 11: multiple stripes
analysis, doi:10.1193/021113EQS025M), with the spread of the two estimates found by simulating
outcomes from the fitted function and estimating again (the same paper's assessment procedure).

With damage_states it also reads the driver's damage_table_all.csv (one row per member) and fits one
curve per damage state of the line network against the nominal stripe: DS1 a peak ratio above 1,
DS2 a line released, DS3 a tower lost, DS4 at least cascade_towers towers lost (N has no default).

    python3 fragility_fit.py --plan plan.json

The plan and the outputs are described in README.md beside this file and in
Docs/sphinx_doc/theory/Conductors.rst ("Fragility curves from an ensemble").
"""
import argparse
import csv
import hashlib
import json
import math
import os
import random
import re
import sys

TOOL = "fragility_fit"
PLAN_KEYS = ("demand_table", "output_dir", "intensity", "by_direction", "capacity", "groups",
             "resamples", "seed", "confidence", "damage_states")
INTENSITIES = ("max_speed", "peak_speed", "stripe")
KINDS = ("conductor", "insulator", "hardware", "footing_uplift", "footing_compression", "member")
NEEDED = ("member", "stripe", "nominal_direction", "seed", "kind", "component", "peak_ratio", "peak_speed",
          "max_speed")
CURVES = "fragility_curves.csv"
POINTS = "fragility_points.csv"
CURVE_COLUMNS = ("curve", "kind", "direction", "intensity", "members", "exceedances", "status", "theta", "beta",
                 "theta_low", "theta_high", "beta_low", "beta_high", "theta_sd", "beta_sd", "resamples_used",
                 "bracket_low", "bracket_high", "log_likelihood")
POINT_COLUMNS = ("curve", "direction", "stripe", "members", "mean_intensity", "exceedances", "fraction", "fitted")
# the damage states of a line network, each read from the member's damage summary
# (damage_table_all.csv of the ensemble driver); a member is in the highest state it reaches
STATES = ("DS1", "DS2", "DS3", "DS4")
STATE_LOG = {"DS1": None, "DS2": "breaks_log", "DS3": "tower_log", "DS4": "tower_log"}
DAMAGE_NEEDED = ("member", "stripe", "nominal_direction", "seed", "max_peak_ratio", "breaks_log", "lines_released",
                 "tower_log", "towers_lost")
MAX_CASCADE = 100000
MAX_RESAMPLES = 100000
FEW = 1.0e-6           # fewer expected exceedances (or survivals) than this over all members count as none
FLAT_RISE = 0.05       # a fitted curve that rises by less than this over the data's winds is not reported


class Stop(Exception):
    """An input the tool refuses; the message names the option."""


def stop(message):
    raise Stop(message)


def fmt(x):
    if x is None:
        return ""
    s = "%.10g" % x
    return "0" if float(s) == 0.0 else s


def is_number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


# ---------------------------------------------------------------------------------------------
# The fragility function and its fit
# ---------------------------------------------------------------------------------------------

def Phi(t):
    """The standard normal distribution function."""
    return 0.5 * math.erfc(-t / math.sqrt(2.0))


def log_Phi(t):
    """ln Phi(t), accurate in both tails: 1 - Phi is taken from erfc of the other sign, and below
    t = -37, where erfc underflows, from the asymptotic series of the normal tail."""
    if t >= 0.0:
        return math.log1p(-0.5 * math.erfc(t / math.sqrt(2.0)))
    if t > -37.0:
        return math.log(0.5 * math.erfc(-t / math.sqrt(2.0)))
    t2 = t * t
    return -0.5 * t2 - math.log(-t) - 0.5 * math.log(2.0 * math.pi) + math.log1p(-1.0 / t2 + 3.0 / (t2 * t2))


def hazard(t):
    """phi(t) / Phi(t), the derivative of ln Phi."""
    return math.exp(-0.5 * t * t - 0.5 * math.log(2.0 * math.pi) - log_Phi(t))


def fragility(x, theta, beta):
    """Baker (2015) Eq. 1: the probability of exceedance at the intensity x."""
    return Phi(math.log(x / theta) / beta)


def grouped(data):
    """[(ln x, n, z, x)] with the observations of one intensity added up: n observations of which z
    (a sum of outcomes, each 0, 1 or a probability in between) exceeded; ascending in x."""
    by_x = {}
    for x, y in data:
        n, z = by_x.get(x, (0.0, 0.0))
        by_x[x] = (n + 1.0, z + y)
    return [(math.log(x), n, z, x) for x, (n, z) in sorted(by_x.items())]


def log_likelihood(groups, a, b):
    """Baker (2015) Eq. 11 without its constant binomial term, for Phi(a + b ln x)."""
    total = 0.0
    for lx, n, z, _ in groups:
        t = a + b * lx
        if z > 0.0:
            total += z * log_Phi(t)
        if n - z > 0.0:
            total += (n - z) * log_Phi(-t)
    return total


def identify(groups):
    """Why a set of observations has no maximum-likelihood fit, or None when it may have one.
    Returns (status, bracket_low, bracket_high): the largest intensity at which no observation
    exceeds and the smallest at which all do (None where the data have no such intensity)."""
    total_n = sum(g[1] for g in groups)
    total_z = sum(g[2] for g in groups)
    none = [x for _, n, z, x in groups if z < FEW]
    every = [x for _, n, z, x in groups if n - z < FEW]
    if total_z < FEW:
        return "no_exceedance", groups[-1][3], None
    if total_n - total_z < FEW:
        return "all_exceed", None, groups[0][3]
    if len(groups) < 2:
        return "one_intensity", None, None
    # the outcomes separate when no intensity with a non-exceedance lies above one with an exceedance:
    # the likelihood then grows as beta goes to zero and has no maximum
    lowest_exceeding = min(x for _, n, z, x in groups if z > 0.0)
    highest_surviving = max(x for _, n, z, x in groups if z < n)
    if highest_surviving <= lowest_exceeding:
        low = max([x for x in none if x <= lowest_exceeding] or [None], key=lambda v: -1.0 if v is None else v)
        high = min([x for x in every if x >= highest_surviving] or [None], key=lambda v: math.inf if v is None else v)
        return "separated", low, high
    return None, None, None


def fit(data):
    """The maximum-likelihood lognormal fragility function of data = [(intensity, outcome)].
    Returns a dict with status "fitted" and theta, beta, log_likelihood, or another status
    (no_exceedance, all_exceed, one_intensity, separated, not_increasing, flat, not_converged) with
    the bracket of the median where the data give one.

    The fit is Baker (2015) Eq. 11 written as the probit regression Phi(a + b ln x) the paper names
    as its equivalent (theta = exp(-a / b), beta = 1 / b), whose log-likelihood is concave in
    (a, b). It is maximised by Fisher scoring with step halving, in the centred and scaled
    variable u = (ln x - mean) / spread so that close or large intensities do not spoil it."""
    groups = grouped(data)
    status, low, high = identify(groups)
    out = {"status": status, "theta": None, "beta": None, "log_likelihood": None,
           "bracket_low": low, "bracket_high": high}
    if status is not None:
        return out
    n = sum(g[1] for g in groups)
    mean = sum(g[0] * g[1] for g in groups) / n
    spread = math.sqrt(sum(g[1] * (g[0] - mean) ** 2 for g in groups) / n)
    if spread <= 0.0:
        out["status"] = "one_intensity"
        return out
    scaled = [((lx - mean) / spread, m, z, x) for lx, m, z, x in groups]
    alpha, gamma = 0.0, 1.0
    ll = log_likelihood(scaled, alpha, gamma)
    settled = False
    for _ in range(500):
        g0 = g1 = h00 = h01 = h11 = 0.0
        for u, m, z, _ in scaled:
            t = alpha + gamma * u
            up, down = hazard(t), hazard(-t)      # phi / Phi and phi / (1 - Phi)
            score = z * up - (m - z) * down
            w = m * up * down                     # the Fisher information of one observation
            g0 += score
            g1 += score * u
            h00 += w
            h01 += w * u
            h11 += w * u * u
        if abs(g0) + abs(g1) <= 1.0e-9 * n:
            settled = True
            break
        det = h00 * h11 - h01 * h01
        if not det > 1.0e-14 * h00 * h11:
            break
        da = (h11 * g0 - h01 * g1) / det
        db = (h00 * g1 - h01 * g0) / det
        step, new = 1.0, None
        while step > 1.0e-12:
            new = log_likelihood(scaled, alpha + step * da, gamma + step * db)
            if new >= ll:
                break
            step *= 0.5
        else:
            break
        alpha, gamma, ll = alpha + step * da, gamma + step * db, new
        if abs(gamma) > 1.0e8:
            break
    if not settled or not all(math.isfinite(v) for v in (alpha, gamma, ll)):
        out["status"] = "not_converged"
        return out
    if gamma <= 1.0e-8:        # the probit does not rise over the spread of the data
        out["status"] = "not_increasing"
        return out
    b = gamma / spread
    a = alpha - gamma * mean / spread
    rise = Phi(a + b * groups[-1][0]) - Phi(a + b * groups[0][0])
    if rise < FLAT_RISE:
        out["status"] = "flat"
        return out
    try:
        theta = math.exp(-a / b)
    except OverflowError:
        out["status"] = "flat"
        return out
    out.update({"status": "fitted", "theta": theta, "beta": 1.0 / b, "log_likelihood": ll})
    return out


def reestimate(data, theta, beta, resamples, confidence, rng):
    """The spread of the estimates, by Baker's (2015) assessment procedure applied to the fitted
    function: outcomes are simulated at the data's intensities from the function with theta and
    beta, the function is estimated from them again, and this is repeated. Simulated sets without
    a fit (no exceedance, all exceeding, separated) are left out and counted. Returns the standard
    deviations and the central confidence interval of theta and beta over the sets that fitted."""
    xs = [x for x, _ in data]
    ps = [fragility(x, theta, beta) for x in xs]
    thetas, betas = [], []
    for _ in range(resamples):
        f = fit([(x, 1.0 if rng.random() < p else 0.0) for x, p in zip(xs, ps)])
        if f["status"] == "fitted":
            thetas.append(f["theta"])
            betas.append(f["beta"])
    out = {"resamples_used": len(thetas)}
    for name, values in (("theta", thetas), ("beta", betas)):
        if len(values) < 2:
            out.update({name + "_sd": None, name + "_low": None, name + "_high": None})
            continue
        mean = sum(values) / len(values)
        out[name + "_sd"] = math.sqrt(sum((v - mean) ** 2 for v in values) / (len(values) - 1))
        values.sort()
        out[name + "_low"] = quantile(values, 0.5 * (1.0 - confidence))
        out[name + "_high"] = quantile(values, 0.5 * (1.0 + confidence))
    return out


def quantile(ordered, q):
    """The q quantile of ascending values, by linear interpolation between order statistics."""
    pos = q * (len(ordered) - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(ordered) - 1)
    return ordered[lo] + (pos - lo) * (ordered[hi] - ordered[lo])


# ---------------------------------------------------------------------------------------------
# Capacity: nominal, or lognormal about the nominal value
# ---------------------------------------------------------------------------------------------

def exceedance_probability(ratio, capacity):
    """The probability that the peak demand exceeds the capacity, given the peak demand-to-nominal
    capacity ratio. Without a capacity distribution it is the model's own rule, demand > capacity.
    With one, the capacity is lognormal with mean mean_factor x nominal and coefficient of
    variation cov, so ln(capacity / nominal) is normal with standard deviation
    zeta = sqrt(ln(1 + cov^2)) and mean ln(mean_factor) - zeta^2 / 2."""
    if capacity is None:
        return 1.0 if ratio > 1.0 else 0.0
    if ratio <= 0.0:
        return 0.0
    zeta = math.sqrt(math.log1p(capacity["cov"] ** 2))
    return Phi((math.log(ratio) - math.log(capacity["mean_factor"]) + 0.5 * zeta * zeta) / zeta)


# ---------------------------------------------------------------------------------------------
# The plan and the table
# ---------------------------------------------------------------------------------------------

class Plan:
    def __init__(self, path):
        self.path = os.path.abspath(path)
        self.directory = os.path.dirname(self.path)
        try:
            with open(self.path, encoding="utf-8") as f:
                p = json.load(f)
        except OSError as e:
            stop("--plan: cannot read %s (%s)" % (path, e.strerror))
        except ValueError as e:
            stop("--plan: %s is not JSON (%s)" % (path, e))
        if not isinstance(p, dict):
            stop("--plan: %s must hold one JSON object" % path)
        unknown = sorted(k for k in p if k not in PLAN_KEYS)
        if unknown:
            stop("--plan: unknown key(s) %s; the plan's keys are %s" % (", ".join(unknown), ", ".join(PLAN_KEYS)))

        self.demand_table = p.get("demand_table")
        if not isinstance(self.demand_table, str) or not self.demand_table.strip():
            stop("demand_table must name the gathered demand table (a path relative to the plan)")
        self.output_dir = p.get("output_dir", "fragility")
        if not isinstance(self.output_dir, str) or not self.output_dir.strip():
            stop("output_dir must name a directory (a path relative to the plan)")
        self.intensity = p.get("intensity", "max_speed")
        if self.intensity not in INTENSITIES:
            stop("intensity must be one of %s, got %s" % (", ".join(INTENSITIES), json.dumps(self.intensity)))
        self.by_direction = p.get("by_direction", False)
        if not isinstance(self.by_direction, bool):
            stop("by_direction must be true or false, got %s" % json.dumps(self.by_direction))

        self.capacity = {}
        cap = p.get("capacity", {})
        if not isinstance(cap, dict):
            stop("capacity must be an object {kind: {\"cov\": .., \"mean_factor\": ..}}")
        for kind, c in cap.items():
            if kind not in KINDS:
                stop("capacity: '%s' is not a kind of component; the kinds are %s" % (kind, ", ".join(KINDS)))
            if not isinstance(c, dict) or sorted(k for k in c if k not in ("cov", "mean_factor")) or "cov" not in c:
                stop("capacity: %s needs {\"cov\": .., \"mean_factor\": ..} (mean_factor is 1 when left out)" % kind)
            cov, mean = c["cov"], c.get("mean_factor", 1.0)
            if not is_number(cov) or not 0.001 <= cov <= 1.0:
                stop("capacity: %s.cov must be in [0.001, 1], got %s" % (kind, json.dumps(cov)))
            if not is_number(mean) or not 0.1 <= mean <= 10.0:
                stop("capacity: %s.mean_factor must be in [0.1, 10], got %s" % (kind, json.dumps(mean)))
            self.capacity[kind] = {"cov": float(cov), "mean_factor": float(mean)}

        self.groups = []
        groups = p.get("groups", {})
        if not isinstance(groups, dict):
            stop("groups must be an object {name: {\"components\": regex, \"kinds\": [..]}}")
        for name, g in groups.items():
            if not name.strip() or "," in name:
                stop("groups: '%s' is not a usable group name" % name)
            if not isinstance(g, dict) or sorted(k for k in g if k not in ("components", "kinds")) or not g:
                stop("groups: %s needs \"components\" (a regular expression), \"kinds\" (a list) or both" % name)
            pattern = g.get("components", ".*")
            try:
                rx = re.compile(pattern) if isinstance(pattern, str) else None
            except re.error:
                rx = None
            if rx is None:
                stop("groups: %s.components must be a regular expression, got %s" % (name, json.dumps(pattern)))
            kinds = g.get("kinds", list(KINDS))
            if not isinstance(kinds, list) or not kinds or any(k not in KINDS for k in kinds):
                stop("groups: %s.kinds must be a list of %s" % (name, ", ".join(KINDS)))
            self.groups.append((name, rx, tuple(kinds)))

        self.resamples = p.get("resamples", 1000)
        if not isinstance(self.resamples, int) or isinstance(self.resamples, bool) \
                or not 0 <= self.resamples <= MAX_RESAMPLES:
            stop("resamples must be an integer in [0, %d], got %s" % (MAX_RESAMPLES, json.dumps(self.resamples)))
        self.seed = p.get("seed", 1)
        if not isinstance(self.seed, int) or isinstance(self.seed, bool) or not 0 <= self.seed <= 2147483647:
            stop("seed must be an integer in [0, 2147483647], got %s" % json.dumps(self.seed))
        self.confidence = p.get("confidence", 0.9)
        if not is_number(self.confidence) or not 0.5 <= self.confidence < 1.0:
            stop("confidence must be in [0.5, 1), got %s" % json.dumps(self.confidence))

        self.damage = None
        if "damage_states" in p:
            ds = p["damage_states"]
            keys = ("table", "states", "cascade_towers")
            if not isinstance(ds, dict) or sorted(k for k in ds if k not in keys) or "table" not in ds \
                    or "states" not in ds:
                stop("damage_states must be {\"table\": .., \"states\": [..], \"cascade_towers\": N} (cascade_towers "
                     "with DS4 only)")
            if not isinstance(ds["table"], str) or not ds["table"].strip():
                stop("damage_states: table must name the gathered damage table (a path relative to the plan)")
            states = ds["states"]
            if not isinstance(states, list) or not states or any(x not in STATES for x in states) \
                    or len(set(states)) != len(states):
                stop("damage_states: states must list some of %s, each once, got %s"
                     % (", ".join(STATES), json.dumps(states)))
            n = ds.get("cascade_towers")
            if "DS4" in states and n is None:
                stop("damage_states: DS4 (a cascade over N towers) needs cascade_towers = N; it has no default, since "
                     "N is the user's choice")
            if n is not None and "DS4" not in states:
                stop("damage_states: cascade_towers is used by DS4 only, which states does not list")
            if n is not None and (not isinstance(n, int) or isinstance(n, bool) or not 2 <= n <= MAX_CASCADE):
                stop("damage_states: cascade_towers must be an integer in [2, %d] (one tower lost is DS3), got %s"
                     % (MAX_CASCADE, json.dumps(n)))
            self.damage = {"table": ds["table"], "states": [x for x in STATES if x in states], "cascade_towers": n}

    def path_of(self, name):
        return os.path.normpath(os.path.join(self.directory, name))


def read_table(path):
    """(rows, notes): the rows of the gathered demand table as dicts with numbers where the fit needs
    them, and what was left out of it."""
    try:
        with open(path, newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            header = reader.fieldnames or []
            rows = list(reader)
    except OSError as e:
        stop("demand_table: cannot read %s (%s)" % (path, e.strerror))
    except (ValueError, csv.Error) as e:
        stop("demand_table: %s is not a readable CSV table (%s)" % (path, e))
    if not rows:
        stop("demand_table: %s holds no rows" % path)
    lacking = [c for c in NEEDED if c not in header]
    if lacking:
        stop("demand_table: %s lacks the column(s) %s" % (path, ", ".join(lacking)))
    out, notes, unrecorded, disagree = [], [], 0, 0
    for n, r in enumerate(rows, start=2):
        if None in r or any(v is None for v in r.values()):
            stop("demand_table: line %d of %s does not have the header's number of fields" % (n, path))
        row = {k: v.strip() for k, v in r.items()}
        if not row["peak_ratio"]:        # the model writes an empty ratio for a component it never recorded
            unrecorded += 1
            continue
        for col in ("stripe", "nominal_direction", "peak_ratio"):
            try:
                row[col] = float(row[col])
            except ValueError:
                stop("demand_table: line %d of %s has no number in %s" % (n, path, col))
            if not math.isfinite(row[col]):
                stop("demand_table: line %d of %s is not finite in %s" % (n, path, col))
        if row["peak_ratio"] < 0.0 or row["stripe"] <= 0.0:
            stop("demand_table: line %d of %s has a negative peak_ratio or a stripe that is not positive" % (n, path))
        for col in ("peak_speed", "max_speed"):      # empty when the run recorded none
            if not row[col]:
                row[col] = None
                continue
            try:
                row[col] = float(row[col])
            except ValueError:
                stop("demand_table: line %d of %s has no number in %s" % (n, path, col))
            if not (math.isfinite(row[col]) and row[col] > 0.0):
                stop("demand_table: line %d of %s has a %s that is not positive and finite" % (n, path, col))
        if row["kind"] not in KINDS:
            stop("demand_table: line %d of %s has the unknown kind '%s'" % (n, path, row["kind"]))
        if row.get("exceeded") in ("0", "1") and (row["exceeded"] == "1") != (row["peak_ratio"] > 1.0):
            disagree += 1
        out.append(row)
    seen, members = set(), {}
    for r in out:
        key = (r["member"], r["component"])
        if key in seen:
            stop("demand_table: member %s lists the component %s twice" % key)
        seen.add(key)
        same = members.setdefault(r["member"], (r["stripe"], r["nominal_direction"], r["seed"]))
        if same != (r["stripe"], r["nominal_direction"], r["seed"]):
            stop("demand_table: the rows of member %s do not agree in stripe, nominal_direction and seed" % r["member"])
    if not out:
        stop("demand_table: %s holds no row with a peak_ratio" % path)
    if unrecorded:
        notes.append("NOTE %d row(s) have no peak_ratio (a component the run never recorded) and are left out" % unrecorded)
    if disagree:
        notes.append("NOTE %d row(s) have an 'exceeded' flag that differs from peak_ratio > 1 (the ratio is printed "
                     "to 10 digits); the fit uses the ratio" % disagree)
    return out, notes


# ---------------------------------------------------------------------------------------------
# The curves
# ---------------------------------------------------------------------------------------------

def intensity_of(row, which):
    return row["stripe"] if which == "stripe" else row[which]


def component_curves(rows, plan):
    """One set of observations per component: [(curve name, kind, direction, observations, members
    left out)] with observations = [(member, stripe, intensity, probability of exceedance)]."""
    by = {}
    for r in rows:
        direction = r["nominal_direction"] if plan.by_direction else None
        by.setdefault((r["component"], r["kind"], direction), []).append(r)
    curves = []
    for (name, kind, direction), rs in by.items():
        obs, skipped = [], []
        for r in rs:
            x = intensity_of(r, plan.intensity)
            if x is None:
                skipped.append(r["member"])
                continue
            obs.append((r["member"], r["stripe"], x, exceedance_probability(r["peak_ratio"], plan.capacity.get(kind))))
        curves.append((name, kind, direction, obs, skipped))
    return curves


def group_curves(rows, plan):
    """One set of observations per group: a member exceeds when any component of the group does.
    With capacity distributions the components' capacities are taken as independent given the run's
    demands, so the member's probability is 1 - prod(1 - p_c) (the form of the series rule of
    Panteli et al. 2017, Eq. 4). The group's intensity is the member's stripe, since its components
    stand at different points. Every member must hold every component of the group."""
    components = set(r["component"] for r in rows)
    for name, _, _ in plan.groups:
        if name in components:
            stop("groups: the group name '%s' is also a component's name" % name)
    curves = []
    for name, rx, kinds in plan.groups:
        chosen = [r for r in rows if r["kind"] in kinds and rx.search(r["component"])]
        if not chosen:
            stop("groups: %s matches no component of the table" % name)
        wanted = set(r["component"] for r in chosen)
        held = {}
        for r in rows:
            held.setdefault(r["member"], set())
        for r in chosen:
            held[r["member"]].add(r["component"])
        for member in sorted(held):
            if held[member] != wanted:
                stop("groups: member %s lacks %d of the %d components of the group %s (first: %s); a group needs "
                     "every component in every member" % (member, len(wanted - held[member]), len(wanted), name,
                                                          sorted(wanted - held[member])[0]))
        by = {}
        for r in chosen:
            direction = r["nominal_direction"] if plan.by_direction else None
            by.setdefault(direction, {}).setdefault(r["member"], []).append(r)
        for direction, members in by.items():
            obs = []
            for member, rs in members.items():
                survive = 1.0
                for r in rs:
                    survive *= 1.0 - exceedance_probability(r["peak_ratio"], plan.capacity.get(r["kind"]))
                obs.append((member, rs[0]["stripe"], rs[0]["stripe"], 1.0 - survive))
            curves.append((name, "group", direction, obs, []))
    return curves


def read_damage_table(path):
    """The members' damage summaries as dicts, with numbers where the states need them."""
    try:
        with open(path, newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            header = reader.fieldnames or []
            rows = list(reader)
    except OSError as e:
        stop("damage_states: cannot read %s (%s)" % (path, e.strerror))
    except (ValueError, csv.Error) as e:
        stop("damage_states: %s is not a readable CSV table (%s)" % (path, e))
    if not rows:
        stop("damage_states: %s holds no rows" % path)
    lacking = [c for c in DAMAGE_NEEDED if c not in header]
    if lacking:
        stop("damage_states: %s lacks the column(s) %s" % (path, ", ".join(lacking)))
    out, seen = [], set()
    for n, r in enumerate(rows, start=2):
        if None in r or any(v is None for v in r.values()):
            stop("damage_states: line %d of %s does not have the header's number of fields" % (n, path))
        row = {k: v.strip() for k, v in r.items()}
        if row["member"] in seen:
            stop("damage_states: member %s is listed twice in %s" % (row["member"], path))
        seen.add(row["member"])
        for col in ("stripe", "nominal_direction"):
            try:
                row[col] = float(row[col])
            except ValueError:
                stop("damage_states: line %d of %s has no number in %s" % (n, path, col))
            if not math.isfinite(row[col]):
                stop("damage_states: line %d of %s is not finite in %s" % (n, path, col))
        if row["stripe"] <= 0.0:
            stop("damage_states: line %d of %s has a stripe that is not positive" % (n, path))
        for log, col in (("breaks_log", "lines_released"), ("tower_log", "towers_lost")):
            if row[log] not in ("yes", "no"):
                stop("damage_states: line %d of %s has %s = '%s'; it is yes or no" % (n, path, log, row[log]))
            if row[log] == "no":
                row[col] = None
                continue
            try:
                row[col] = int(row[col])
            except ValueError:
                stop("damage_states: line %d of %s has no whole number in %s" % (n, path, col))
            if row[col] < 0:
                stop("damage_states: line %d of %s has a negative %s" % (n, path, col))
        if row["max_peak_ratio"]:
            try:
                row["max_peak_ratio"] = float(row["max_peak_ratio"])
            except ValueError:
                stop("damage_states: line %d of %s has no number in max_peak_ratio" % (n, path))
            if not (math.isfinite(row["max_peak_ratio"]) and row["max_peak_ratio"] >= 0.0):
                stop("damage_states: line %d of %s has a max_peak_ratio that is negative or not finite" % (n, path))
        else:
            row["max_peak_ratio"] = None
        out.append(row)
    return out


def damage_curves(rows, plan, notes):
    """One set of observations per damage state: a member is in the highest of the listed states it reaches,
    and counts as exceeding every listed state up to that one. DS1 is a peak ratio above 1 (the rule of the
    component curves without a capacity distribution), DS2 a line released (breaks.dat), DS3 a tower lost and
    DS4 at least cascade_towers towers lost (tower_failures.dat). Every member must hold the log a listed
    state is read from; the intensity is the member's stripe."""
    states, n = plan.damage["states"], plan.damage["cascade_towers"]
    for state in states:
        log = STATE_LOG[state]
        lacking = [r["member"] for r in rows if (r["max_peak_ratio"] is None if log is None else r[log] != "yes")]
        if lacking:
            why = "no peak ratio" if log is None else "no %s (its option was off)" % (
                "breaks.dat" if log == "breaks_log" else "tower_failures.dat")
            stop("damage_states: %s needs every member's record, but %d member(s) have %s: %s%s"
                 % (state, len(lacking), why, ", ".join(lacking[:5]), ", ..." if len(lacking) > 5 else ""))
    reached = {"DS1": lambda r: r["max_peak_ratio"] > 1.0, "DS2": lambda r: r["lines_released"] >= 1,
               "DS3": lambda r: r["towers_lost"] >= 1, "DS4": lambda r: r["towers_lost"] >= n}
    by, odd = {}, []
    for r in rows:
        flags = [reached[s](r) for s in states]
        level = max([i + 1 for i, f in enumerate(flags) if f] or [0])
        if not all(flags[:level]):
            odd.append(r["member"])
        direction = r["nominal_direction"] if plan.by_direction else None
        by.setdefault(direction, []).append((r["member"], r["stripe"], level))
    if odd:
        notes.append("NOTE %d member(s) reach a damage state without every lower one listed (%s%s); each counts in "
                     "the highest it reaches and every one below" % (len(odd), ", ".join(odd[:5]),
                                                                     ", ..." if len(odd) > 5 else ""))
    curves = []
    for direction, members in by.items():
        for i, state in enumerate(states):
            obs = [(m, x, x, 1.0 if level > i else 0.0) for m, x, level in members]
            curves.append((state, "damage_state", direction, obs, []))
    return curves


def curve_seed(seed, name, kind, direction):
    """A seed of the curve's own, so that its spread does not depend on the other curves or their order."""
    text = "%d|%s|%s|%s" % (seed, name, kind, direction)
    return int(hashlib.sha256(text.encode("utf-8")).hexdigest()[:16], 16)


def run(plan, log=print):
    rows, notes = read_table(plan.path_of(plan.demand_table))
    curves = component_curves(rows, plan) + group_curves(rows, plan)
    if plan.damage is not None:
        names = set(r["component"] for r in rows) | set(g[0] for g in plan.groups)
        for state in plan.damage["states"]:
            if state in names:
                stop("damage_states: %s is also the name of a component or group" % state)
        curves += damage_curves(read_damage_table(plan.path_of(plan.damage["table"])), plan, notes)
        if plan.capacity:
            notes.append("NOTE the damage states are the runs' own outcomes, at the nominal capacities; capacity "
                         "applies to the component and group curves only")
        if plan.damage["cascade_towers"] is not None:
            notes.append("NOTE DS4 is %d or more towers lost (cascade_towers, the user's choice)"
                         % plan.damage["cascade_towers"])
    curve_rows, point_rows = [], []
    for name, kind, direction, obs, skipped in curves:
        d = "all" if direction is None else fmt(direction)
        used = "stripe" if kind in ("group", "damage_state") else plan.intensity
        if skipped:
            notes.append("NOTE %s (direction %s): %d member(s) recorded no %s and are left out: %s%s"
                         % (name, d, len(skipped), plan.intensity, ", ".join(skipped[:5]),
                            ", ..." if len(skipped) > 5 else ""))
        if not obs:
            curve_rows.append([name, kind, d, used, "0", "0", "no_intensity"] + [""] * 12)
            continue
        data = [(x, y) for _, _, x, y in obs]
        f = fit(data)
        spread = {}
        if f["status"] == "fitted" and plan.resamples:
            if all(y in (0.0, 1.0) for _, y in data):
                spread = reestimate(data, f["theta"], f["beta"], plan.resamples, plan.confidence,
                                    random.Random(curve_seed(plan.seed, name, kind, d)))
                if spread["resamples_used"] < 0.9 * plan.resamples:
                    notes.append("NOTE %s (direction %s): %d of %d simulated sets had a fit; its spread describes "
                                 "those sets only" % (name, d, spread["resamples_used"], plan.resamples))
            else:
                notes.append("NOTE %s (direction %s): no spread is given for outcomes that are probabilities "
                             "(a capacity distribution)" % (name, d))
        curve_rows.append([name, kind, d, used, str(len(obs)), fmt(sum(y for _, y in data)), f["status"],
                           fmt(f["theta"]), fmt(f["beta"]),
                           fmt(spread.get("theta_low")), fmt(spread.get("theta_high")),
                           fmt(spread.get("beta_low")), fmt(spread.get("beta_high")),
                           fmt(spread.get("theta_sd")), fmt(spread.get("beta_sd")),
                           "" if not spread else str(spread["resamples_used"]),
                           fmt(f["bracket_low"]), fmt(f["bracket_high"]), fmt(f["log_likelihood"])])
        stripes = {}
        for _, stripe, x, y in obs:
            st = stripes.setdefault(stripe, [0, 0.0, 0.0])
            st[0] += 1
            st[1] += x
            st[2] += y
        for stripe in sorted(stripes):
            n, sx, sy = stripes[stripe]
            fitted = fragility(sx / n, f["theta"], f["beta"]) if f["status"] == "fitted" else None
            point_rows.append([name, d, fmt(stripe), str(n), fmt(sx / n), fmt(sy), fmt(sy / n), fmt(fitted)])
    out = plan.path_of(plan.output_dir)
    try:
        os.makedirs(out, exist_ok=True)
        files = ((CURVES, CURVE_COLUMNS, curve_rows), (POINTS, POINT_COLUMNS, point_rows))
        for file_name, header, body in files:
            with open(os.path.join(out, file_name + ".tmp"), "w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh, lineterminator="\n")
                w.writerow(header)
                w.writerows(body)
        for file_name, _, _ in files:        # both written before either replaces an earlier run's
            os.replace(os.path.join(out, file_name + ".tmp"), os.path.join(out, file_name))
    except OSError as e:
        stop("output_dir: cannot write to %s (%s)" % (out, e.strerror))
    counts = {}
    for r in curve_rows:
        counts[r[6]] = counts.get(r[6], 0) + 1
    log("%d curves in %s: %s" % (len(curve_rows), os.path.join(out, CURVES),
                                 ", ".join("%d %s" % (counts[k], k) for k in sorted(counts))))
    for note in notes:
        log(note)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog=TOOL + ".py", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plan", required=True, help="the plan (JSON)")
    a = ap.parse_args(sys.argv[1:] if argv is None else argv)
    try:
        return run(Plan(a.plan))
    except Stop as e:
        print("%s: %s" % (TOOL, e), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
