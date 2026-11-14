#!/usr/bin/env python3
"""Tests of fragility_fit.py (standard library only; no ERF run).

    python3 test_fragility_fit.py            every test
    python3 test_fragility_fit.py Published  one class

Published: numbers printed in Baker (2015), "Efficient analytical fragility function fitting using
dynamic structural analysis", Earthquake Spectra 31(1), 579-599, doi:10.1193/021113EQS025M.
ByHand: fits whose answer follows from the normal distribution's values Phi(0) = 1/2 and
Phi(1) = 0.8413447461, and a brute-force search of the likelihood that shares no code with the
tool's optimiser. Table: the stored tables ensemble.csv and two_stripes.csv through the command.
Refused: every plan key out of range.
"""
import contextlib
import csv
import importlib.util
import io
import json
import math
import os
import random
import shutil
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("fragility_fit", os.path.join(HERE, "..", "fragility_fit.py"))
ff = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ff)

PHI_1 = 0.8413447461        # Phi(1)
LN_1P5 = 0.4054651081       # ln(30 / 20)


def simulated_spread(levels, per_level, sets, seed):
    """Baker's assessment procedure with the assumed function theta = 1, beta = 0.4: the mean and
    standard deviation of the estimates over simulated sets of outcomes."""
    rng = random.Random(seed)
    thetas, betas = [], []
    for _ in range(sets):
        data = [(x, 1.0 if rng.random() < ff.fragility(x, 1.0, 0.4) else 0.0) for x in levels for _ in range(per_level)]
        f = ff.fit(data)
        if f["status"] == "fitted":
            thetas.append(f["theta"])
            betas.append(f["beta"])

    def stats(v):
        m = sum(v) / len(v)
        return m, math.sqrt(sum((a - m) ** 2 for a in v) / (len(v) - 1))
    return stats(thetas), stats(betas), len(thetas)


class Published(unittest.TestCase):
    def test_the_standard_deviation_of_the_median_with_three_stripes(self):
        # p. 589: 40 ground motions at each of IM = {0.6, 1, 1.5} "resulted in a standard deviation of theta of
        # 0.056"; with 20 "the standard deviation of theta increases to 0.078". 20,000 simulated sets give 0.054
        # and 0.079 with this estimator; the 4,000 used here have a sampling error near 0.001.
        (mean, sd), _, used = simulated_spread((0.6, 1.0, 1.5), 40, 4000, seed=2015)
        self.assertGreater(used, 3980)      # a few sets rise too little to be reported
        self.assertAlmostEqual(sd, 0.056, delta=0.003)
        self.assertAlmostEqual(mean, 1.0, delta=0.01)       # "essentially 1": unbiased
        (_, sd), _, _ = simulated_spread((0.6, 1.0, 1.5), 20, 4000, seed=2015)
        self.assertAlmostEqual(sd, 0.078, delta=0.004)

    def test_the_coefficients_of_variation_of_table_3(self):
        # Table 3, multiple stripes with 45 ground motions at each of IM = [0.4, 0.8, 1.2]: COV of theta 0.06, COV
        # of beta 0.20, printed to two digits from 1,000 simulated sets (a sampling error near 0.005 for beta);
        # 4,000 sets here give 0.057 and 0.19
        (mt, st), (mb, sb), _ = simulated_spread((0.4, 0.8, 1.2), 45, 4000, seed=3)
        self.assertAlmostEqual(st / mt, 0.06, delta=0.006)
        self.assertAlmostEqual(sb / mb, 0.20, delta=0.02)
        # The column with 30 ground motions per stripe prints 0.08 and 0.26. It is not reproduced: 20,000 sets
        # give 0.070 and 0.23 with this estimator and with an independent one, well outside the sampling error
        # (0.001), so it is recorded here and not asserted as agreement.
        (mt, st), (mb, sb), _ = simulated_spread((0.4, 0.8, 1.2), 30, 4000, seed=3)
        self.assertAlmostEqual(st / mt, 0.070, delta=0.005)
        self.assertAlmostEqual(sb / mb, 0.23, delta=0.02)
        # the two-stripe column (IM = [0.5, 1.2], COV of beta 0.40) is not compared: about one set in seven has no
        # collapse at 0.5, and such a set has no maximum (the paper: "if one produces no collapses there is no
        # unique solution"). This tool names those sets "separated" and leaves them out; the paper does not say
        # how it treated them.
        _, _, used = simulated_spread((0.5, 1.2), 45, 2000, seed=3)
        self.assertTrue(0.80 * 2000 < used < 0.90 * 2000, used)      # (1 - Phi(ln 0.5 / 0.4))^45 = 0.148 left out

    def test_equation_1(self):
        # the median has an exceedance probability of one half; one beta above it, Phi(1)
        self.assertAlmostEqual(ff.fragility(30.0, 30.0, 0.4), 0.5, places=12)
        self.assertAlmostEqual(ff.fragility(30.0 * math.exp(0.4), 30.0, 0.4), PHI_1, places=10)
        self.assertAlmostEqual(ff.fragility(30.0 * math.exp(-0.4), 30.0, 0.4), 1.0 - PHI_1, places=10)


class ByHand(unittest.TestCase):
    def test_two_stripes_are_fitted_exactly(self):
        # one half exceed at 20 m/s and Phi(1) at 30 m/s: the curve through both has theta = 20 and
        # beta = ln(30 / 20) / 1
        f = ff.fit([(20.0, 1.0), (20.0, 0.0), (30.0, PHI_1)])
        self.assertEqual(f["status"], "fitted")
        self.assertAlmostEqual(f["theta"], 20.0, places=6)
        self.assertAlmostEqual(f["beta"], LN_1P5, places=7)
        # counts: 5 of 10 at 20 m/s; at 40 m/s the fraction Phi(1): beta = ln 2
        f = ff.fit([(20.0, 1.0)] * 5 + [(20.0, 0.0)] * 5 + [(40.0, PHI_1)] * 10)
        self.assertAlmostEqual(f["theta"], 20.0, places=6)
        self.assertAlmostEqual(f["beta"], math.log(2.0), places=7)

    def test_the_optimum_is_the_brute_force_optimum(self):
        # Eq. 11 searched on a grid of theta and beta, with its own likelihood sum
        data = [(19.5, 0.0), (20.5, 1.0), (29.5, 0.0), (30.5, 1.0), (39.5, 1.0), (40.5, 1.0),
                (25.0, 0.0), (25.0, 0.0), (35.0, 1.0), (35.0, 0.0), (45.0, 1.0)]

        def loglik(theta, beta):
            total = 0.0
            for x, y in data:
                p = min(max(0.5 * math.erfc(-math.log(x / theta) / beta / math.sqrt(2.0)), 1.0e-300), 1.0 - 1.0e-16)
                total += math.log(p) if y else math.log(1.0 - p)
            return total
        best = max((loglik(20.0 + 0.02 * i, 0.05 + 0.002 * j), 20.0 + 0.02 * i, 0.05 + 0.002 * j)
                   for i in range(1000) for j in range(600))
        f = ff.fit(data)
        self.assertEqual(f["status"], "fitted")
        self.assertAlmostEqual(f["theta"], best[1], delta=0.02)
        self.assertAlmostEqual(f["beta"], best[2], delta=0.002)
        self.assertGreaterEqual(f["log_likelihood"], best[0] - 1.0e-9)
        self.assertAlmostEqual(f["log_likelihood"], loglik(f["theta"], f["beta"]), places=9)

    def test_a_shallow_curve_is_found_too(self):
        # fractions that hardly rise with the intensity (26 of 40, 22 of 40, 3 of 3, 3 of 5): the first full
        # step of the optimiser overshoots here, so this fit needs its shortened steps. The optimum is again
        # that of a grid search.
        counts = ((60.49, 40, 26), (60.63, 40, 22), (64.64, 3, 3), (75.68, 5, 3))
        data = [(x, 1.0 if i < z else 0.0) for x, n, z in counts for i in range(n)]

        def loglik(theta, beta):
            total = 0.0
            for x, n, z in counts:
                p = 0.5 * math.erfc(-math.log(x / theta) / beta / math.sqrt(2.0))
                total += z * math.log(p) + (n - z) * math.log(1.0 - p)
            return total
        best = max((loglik(25.0 + 0.05 * i, 0.8 + 0.005 * j), 25.0 + 0.05 * i, 0.8 + 0.005 * j)
                   for i in range(600) for j in range(400))
        f = ff.fit(data)
        self.assertEqual(f["status"], "fitted")
        self.assertAlmostEqual(f["theta"], best[1], delta=0.06)
        self.assertAlmostEqual(f["beta"], best[2], delta=0.006)
        self.assertGreaterEqual(f["log_likelihood"], best[0] - 1.0e-9)

    def test_outcomes_without_a_fit_are_named_with_their_bracket(self):
        cases = (
            ([(20.0, 0.0), (30.0, 0.0), (40.0, 0.0)], ("no_exceedance", 40.0, None)),
            ([(20.0, 1.0), (30.0, 1.0)], ("all_exceed", None, 20.0)),
            ([(30.0, 1.0), (30.0, 0.0)], ("one_intensity", None, None)),
            # none exceed up to 29.5 m/s and all from 30.5 m/s: the median lies between
            ([(19.5, 0.0), (29.5, 0.0), (30.5, 1.0), (40.5, 1.0)], ("separated", 29.5, 30.5)),
            # the same with both outcomes at one intensity, 30 m/s: the likelihood still grows as beta goes to zero;
            # the bracket is the largest wind without an exceedance and the smallest at which all exceed
            ([(20.0, 0.0), (30.0, 0.0), (30.0, 1.0), (40.0, 1.0)], ("separated", 20.0, 40.0)),
            # none at 0.5 and 30 of 45 at 1.2 (the paper's two-stripe case): no wind at which all exceed
            ([(0.5, 0.0)] * 45 + [(1.2, 1.0)] * 30 + [(1.2, 0.0)] * 15, ("separated", 0.5, None)),
            # probabilities: none at 20 and 30 m/s, 0.3 at 40 m/s
            ([(20.0, 0.0), (30.0, 0.0), (40.0, 0.3)], ("separated", 30.0, None)),
            # the same fraction at every intensity: the curve does not rise
            ([(20.0, 1.0), (20.0, 0.0), (30.0, 1.0), (30.0, 0.0)], ("not_increasing", None, None)),
            ([(20.0, 0.3), (30.0, 0.3)], ("not_increasing", None, None)),
            # 1 of 20 at the middle intensity only: a best curve exists but rises by less than 0.05 over the data
            ([(20.0, 0.0)] * 20 + [(30.0, 1.0)] + [(30.0, 0.0)] * 19 + [(40.0, 0.0)] * 20, ("flat", None, None)),
            # outcomes too small to count as an exceedance (a capacity far above the demand)
            ([(20.0, 1.0e-30), (30.0, 2.0e-30), (40.0, 4.0e-9)], ("no_exceedance", 40.0, None)),
            # more exceed at the lower intensity: no increasing curve fits
            ([(20.0, 1.0), (20.0, 1.0), (20.0, 0.0), (40.0, 1.0), (40.0, 0.0), (40.0, 0.0)], ("not_increasing", None, None)))
        for data, (status, low, high) in cases:
            with self.subTest(status=status, data=data):
                f = ff.fit(data)
                self.assertEqual((f["status"], f["bracket_low"], f["bracket_high"]), (status, low, high))
                self.assertIsNone(f["theta"])
                self.assertIsNone(f["beta"])

    def test_a_lognormal_capacity_gives_the_probability_of_exceedance(self):
        # cov = sqrt(exp(0.04) - 1) makes zeta = 0.2, and mean_factor = exp(0.02) puts the median capacity at the
        # nominal value: a peak ratio of 1 exceeds with probability 1/2, a ratio of exp(0.2) with Phi(1)
        cap = {"cov": math.sqrt(math.exp(0.04) - 1.0), "mean_factor": math.exp(0.02)}
        self.assertAlmostEqual(ff.exceedance_probability(1.0, cap), 0.5, places=12)
        self.assertAlmostEqual(ff.exceedance_probability(math.exp(0.2), cap), PHI_1, places=10)
        self.assertAlmostEqual(ff.exceedance_probability(math.exp(-0.2), cap), 1.0 - PHI_1, places=10)
        self.assertEqual(ff.exceedance_probability(0.0, cap), 0.0)
        # a mean capacity 10 % above nominal moves the ratio with probability 1/2 to 1.1 exp(-0.02)
        self.assertAlmostEqual(ff.exceedance_probability(1.1 * math.exp(-0.02), {"cov": cap["cov"], "mean_factor": 1.1}),
                               0.5, places=12)
        # the nominal capacity: the model's own rule, demand above capacity
        self.assertEqual([ff.exceedance_probability(r, None) for r in (0.999, 1.0, 1.001)], [0.0, 0.0, 1.0])

    def test_the_spread_is_that_of_the_assessment_procedure(self):
        # 40 observations at each of 0.6, 1 and 1.5 with theta = 1 and beta = 0.4: the standard deviation the paper
        # prints, 0.056, from the tool's own routine; the interval is the central 90 % of the estimates
        data = [(x, 0.0) for x in (0.6, 1.0, 1.5) for _ in range(40)]
        s = ff.reestimate(data, 1.0, 0.4, 3000, 0.9, random.Random(7))
        self.assertGreater(s["resamples_used"], 2990)
        self.assertAlmostEqual(s["theta_sd"], 0.056, delta=0.003)
        self.assertLess(s["theta_low"], 1.0 - 1.2 * s["theta_sd"])
        self.assertGreater(s["theta_high"], 1.0 + 1.2 * s["theta_sd"])
        self.assertLess(s["beta_low"], 0.4)
        self.assertGreater(s["beta_high"], 0.4)
        wide = ff.reestimate(data, 1.0, 0.4, 3000, 0.99, random.Random(7))
        self.assertLess(wide["theta_low"], s["theta_low"])
        self.assertGreater(wide["theta_high"], s["theta_high"])
        self.assertEqual(ff.quantile([1.0, 2.0, 3.0, 4.0, 5.0], 0.05), 1.2)
        self.assertEqual(ff.quantile([1.0, 2.0, 3.0, 4.0, 5.0], 0.5), 3.0)

    def test_separated_sets_count_as_beta_zero(self):
        # 5 stripes of 3 with 1, 2, 3, 3, 3 exceeding (the PowerLines DS1 counts, 13 to 17 m/s): most simulated
        # sets separate. Their likelihood's supremum is at beta -> 0, so they count as beta = 0 in beta's interval
        # and the lower bound is 0; left out, the bound would sit at the fitted beta itself
        data = [(x, 1.0 if k < n else 0.0) for x, n in zip((13.0, 14.0, 15.0, 16.0, 17.0), (1, 2, 3, 3, 3))
                for k in range(3)]
        f = ff.fit(data)
        self.assertEqual(f["status"], "fitted")
        s = ff.reestimate(data, f["theta"], f["beta"], 1000, 0.9, random.Random(3))
        self.assertGreater(s["resamples_separated"], 300)
        self.assertLessEqual(s["resamples_used"] + s["resamples_separated"], 1000)
        self.assertEqual(s["beta_low"], 0.0)
        self.assertGreater(s["beta_high"], f["beta"])
        self.assertGreater(s["theta_high"], s["theta_low"])
        # no separated set: the interval is over the fitted sets alone, as before
        wide = [(x, 0.0) for x in (0.6, 1.0, 1.5) for _ in range(40)]
        t = ff.reestimate(wide, 1.0, 0.4, 200, 0.9, random.Random(7))
        self.assertEqual(t["resamples_separated"], 0)
        self.assertGreater(t["beta_low"], 0.2)

    def test_crossing_damage_state_curves_are_named(self):
        # DS3 (14.028, 0.0913) and DS4 (14.480, 0.0501) over stripes 13 to 17 cross where
        # ln x = (b3 ln t4 - b4 ln t3) / (b3 - b4) = (0.0913 ln 14.480 - 0.0501 ln 14.028) / 0.0412 = 2.71134, x = 15.049
        notes = ff.crossing_notes({"all": [("DS3", 14.028, 0.0913, 13.0, 17.0), ("DS4", 14.480, 0.0501, 13.0, 17.0)]})
        self.assertEqual(len(notes), 1)
        self.assertIn("DS3 and DS4 (direction all): the curves cross at 15.049", notes[0])
        # the PowerLines fits cross at 11.25 m/s, below the stripes: no note
        self.assertEqual(ff.crossing_notes({"all": [("DS3", 13.43696145, 0.05893949262, 13.0, 17.0),
                                                    ("DS4", 13.74233833, 0.06638969011, 13.0, 17.0)]}), [])
        # a higher state's median below a lower one's is always named
        notes = ff.crossing_notes({"0": [("DS2", 20.0, 0.1, 10.0, 30.0), ("DS3", 19.0, 0.1, 10.0, 30.0)]})
        self.assertEqual(len(notes), 1)
        self.assertIn("DS3 (direction 0): its median 19 lies below DS2's median 20", notes[0])

    def test_tails_and_close_intensities(self):
        # ln Phi far in the tail, against the series -t^2/2 - ln(-t) - ln(2 pi)/2 + ln(1 - 1/t^2 + 3/t^4) worked
        # here at t = -40: -800 - 3.688879454 - 0.918938533 - 0.000624024 = -804.608442011
        self.assertAlmostEqual(ff.log_Phi(-40.0), -804.608442011, places=6)
        self.assertAlmostEqual(ff.log_Phi(0.0), math.log(0.5), places=14)
        self.assertAlmostEqual(ff.log_Phi(-1.0), math.log(1.0 - PHI_1), places=9)
        self.assertAlmostEqual(ff.log_Phi(10.0), -7.619853024e-24, delta=1e-30)     # -(1 - Phi(10))
        # winds a hundred-thousandth apart about 1000 m/s: the fit is that of the same pattern about 1 with the
        # spacing scaled, theta at the middle and a beta 1000 times smaller in ln x
        pattern = ((0.0, 0.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0))
        near = ff.fit([(1000.0 + 1.0e-5 * k, y) for k, y in pattern])
        wide = ff.fit([(math.exp(k), y) for k, y in pattern])
        self.assertEqual((near["status"], wide["status"]), ("fitted", "fitted"))
        self.assertAlmostEqual(math.log(wide["theta"]), 1.5, places=6)              # the pattern's centre
        self.assertAlmostEqual(near["theta"], 1000.000015, places=6)
        self.assertAlmostEqual(near["beta"] / wide["beta"], 1.0e-8, delta=1.0e-11)


class Case(unittest.TestCase):
    maxDiff = None

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="fragility_fit_")
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def plan(self, table="ensemble.csv", **keys):
        shutil.copy(os.path.join(HERE, table), self.tmp)
        p = {"demand_table": table, "resamples": 0}
        p.update(keys)
        p = {k: v for k, v in p.items() if v is not None}
        path = os.path.join(self.tmp, "plan.json")
        with open(path, "w") as f:
            json.dump(p, f)
        return path

    def call(self, path):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = ff.main(["--plan", path])
        return code, out.getvalue(), err.getvalue()

    def table(self, name, directory="fragility"):
        with open(os.path.join(self.tmp, directory, name), newline="") as f:
            return list(csv.DictReader(f))


class Table(Case):
    def test_recorded_winds_separate_the_conductor_and_fit_the_hardware(self):
        # ensemble.csv: the recorded wind is 19.5, 20.5, 29.5, 30.5, 39.5, 40.5 m/s. The conductor S1 exceeds its
        # nominal capacity (peak ratio above 1) from 30.5 m/s on and never below: separated between 29.5 and
        # 30.5. The hardware exceeds at 20.5, 30.5, 39.5 and 40.5 and not at 19.5 and 29.5: a fit exists.
        code, out, err = self.call(self.plan())
        self.assertEqual((code, err), (0, ""))
        self.assertIn("2 curves", out)
        self.assertIn("1 fitted, 1 separated", out)
        rows = {r["curve"]: r for r in self.table("fragility_curves.csv")}
        self.assertEqual(sorted(rows), ["S1", "S1_end_a"])
        s1 = rows["S1"]
        self.assertEqual((s1["kind"], s1["direction"], s1["members"], s1["exceedances"], s1["status"]),
                         ("conductor", "all", "6", "3", "separated"))
        self.assertEqual((s1["bracket_low"], s1["bracket_high"], s1["theta"], s1["beta"]), ("29.5", "30.5", "", ""))
        self.assertEqual(s1["intensity"], "max_speed")
        hw = rows["S1_end_a"]
        self.assertEqual((hw["kind"], hw["members"], hw["exceedances"], hw["status"]), ("hardware", "6", "4", "fitted"))
        data = [(19.5, 0), (20.5, 1), (29.5, 0), (30.5, 1), (39.5, 1), (40.5, 1)]

        def loglik(theta, beta):
            total = 0.0
            for x, y in data:
                p = min(max(0.5 * math.erfc(-math.log(x / theta) / beta / math.sqrt(2.0)), 1.0e-300), 1.0 - 1.0e-16)
                total += math.log(p) if y else math.log(1.0 - p)
            return total
        best = max((loglik(15.0 + 0.02 * i, 0.1 + 0.005 * j), 15.0 + 0.02 * i, 0.1 + 0.005 * j)
                   for i in range(1000) for j in range(500))
        self.assertAlmostEqual(float(hw["theta"]), best[1], delta=0.03)
        self.assertAlmostEqual(float(hw["beta"]), best[2], delta=0.006)
        self.assertEqual(hw["theta_sd"], "")        # resamples = 0: no spread

        points = [(p["curve"], p["stripe"], p["members"], p["mean_intensity"], p["exceedances"], p["fraction"])
                  for p in self.table("fragility_points.csv")]
        self.assertEqual(points, [
            ("S1", "20", "2", "20", "0", "0"), ("S1", "30", "2", "30", "1", "0.5"), ("S1", "40", "2", "40", "2", "1"),
            ("S1_end_a", "20", "2", "20", "1", "0.5"), ("S1_end_a", "30", "2", "30", "1", "0.5"),
            ("S1_end_a", "40", "2", "40", "2", "1")])
        fitted = [p["fitted"] for p in self.table("fragility_points.csv")]
        self.assertEqual(fitted[:3], ["", "", ""])
        self.assertAlmostEqual(float(fitted[4]), ff.fragility(30.0, float(hw["theta"]), float(hw["beta"])), places=9)

    def test_the_stripe_as_intensity_puts_both_outcomes_at_30(self):
        # with the nominal stripe the conductor has 0 of 2, 1 of 2 and 2 of 2 at 20, 30 and 40 m/s: still
        # separated; none exceed at 20 m/s and all at 40 m/s
        self.call(self.plan(intensity="stripe"))
        s1 = {r["curve"]: r for r in self.table("fragility_curves.csv")}["S1"]
        self.assertEqual((s1["status"], s1["bracket_low"], s1["bracket_high"], s1["intensity"]),
                         ("separated", "20", "40", "stripe"))
        # a group is fitted against the stripe whatever the components use
        self.call(self.plan(groups={"both": {"kinds": ["conductor", "hardware"]}}))
        used = {r["curve"]: r["intensity"] for r in self.table("fragility_curves.csv")}
        self.assertEqual(used, {"S1": "max_speed", "S1_end_a": "max_speed", "both": "stripe"})
        # the wind at the peak demand is the same column of numbers in this table
        self.call(self.plan(intensity="peak_speed"))
        s1 = {r["curve"]: r for r in self.table("fragility_curves.csv")}["S1"]
        self.assertEqual((s1["bracket_low"], s1["bracket_high"], s1["intensity"]), ("29.5", "30.5", "peak_speed"))

    def test_a_capacity_distribution_and_groups(self):
        # two_stripes.csv: peak ratios 1 at 20 m/s and exp(0.2) at 30 m/s for the conductor, 1 at both for the
        # hardware. With zeta = 0.2 about the nominal capacity the conductor exceeds with probability 1/2 and
        # Phi(1): theta = 20 and beta = ln 1.5. The hardware keeps its nominal capacity and never exceeds, so the
        # group of both is the conductor's curve.
        cap = {"cov": math.sqrt(math.exp(0.04) - 1.0), "mean_factor": math.exp(0.02)}
        plan = self.plan("two_stripes.csv", intensity="stripe", capacity={"conductor": cap}, resamples=50,
                         groups={"line_S1": {"components": "^S1"}, "ends": {"kinds": ["hardware"]}})
        code, out, err = self.call(plan)
        self.assertEqual((code, err), (0, ""))
        rows = {r["curve"]: r for r in self.table("fragility_curves.csv")}
        self.assertEqual(sorted(rows), ["S1", "S1_end_a", "ends", "line_S1"])
        for name in ("S1", "line_S1"):
            self.assertEqual(rows[name]["status"], "fitted")
            self.assertAlmostEqual(float(rows[name]["theta"]), 20.0, places=5)
            self.assertAlmostEqual(float(rows[name]["beta"]), LN_1P5, places=6)
            self.assertAlmostEqual(float(rows[name]["exceedances"]), 2 * 0.5 + 2 * PHI_1, places=8)
            self.assertEqual((rows[name]["theta_sd"], rows[name]["resamples_used"]), ("", ""))
        self.assertEqual((rows["line_S1"]["kind"], rows["line_S1"]["members"], rows["line_S1"]["intensity"]),
                         ("group", "4", "stripe"))
        self.assertIn("NOTE S1 (direction all): no spread is given for outcomes that are probabilities", out)
        self.assertEqual((rows["S1_end_a"]["status"], rows["ends"]["status"]), ("no_exceedance", "no_exceedance"))
        self.assertEqual(rows["ends"]["bracket_low"], "30")

        # with the same distribution on the hardware, each end exceeds with probability 1/2 at both stripes, and
        # the line with 1 - (1/2)(1/2) = 0.75 at 20 m/s and 1 - (1 - Phi(1))(1/2) = 0.9206723731 at 30 m/s:
        # theta and beta follow from the two normal quantiles, found here by bisection
        plan = self.plan("two_stripes.csv", intensity="stripe", capacity={"conductor": cap, "hardware": cap},
                         groups={"line_S1": {"components": "^S1"}})
        self.assertEqual(self.call(plan)[0], 0)
        line = {r["curve"]: r for r in self.table("fragility_curves.csv")}["line_S1"]

        def quantile(p):
            lo, hi = -8.0, 8.0
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if 0.5 * math.erfc(-mid / math.sqrt(2.0)) < p else (lo, mid)
            return 0.5 * (lo + hi)
        t20, t30 = quantile(0.75), quantile(0.9206723731)
        self.assertAlmostEqual(t20, 0.6744897502, places=9)      # the upper quartile of the normal distribution
        beta = LN_1P5 / (t30 - t20)
        self.assertAlmostEqual(float(line["beta"]), beta, places=6)
        self.assertAlmostEqual(float(line["theta"]), 20.0 * math.exp(-t20 * beta), places=5)
        self.assertAlmostEqual(float(line["exceedances"]), 2 * 0.75 + 2 * 0.9206723731, places=8)

    def test_by_direction_and_the_spread_columns(self):
        # the ensemble's rows again as a second direction with every hardware ratio below 1: one curve per
        # direction, and the spread of the fitted one is repeatable with the seed
        with open(os.path.join(HERE, "ensemble.csv")) as f:
            lines = f.read().splitlines()
        extra = []
        for ln in lines[1:]:
            c = ln.split(",")
            c[0], c[2] = c[0].replace("_d0_", "_d90_"), "90"
            if c[5] == "S1_end_a":
                c[6] = "0.5"
            extra.append(",".join(c))
        path = self.plan(by_direction=True, resamples=200, seed=11, confidence=0.8)
        with open(os.path.join(self.tmp, "ensemble.csv"), "w") as f:
            f.write("\n".join(lines + extra) + "\n")
        self.assertEqual(self.call(path)[0], 0)
        rows = {(r["curve"], r["direction"]): r for r in self.table("fragility_curves.csv")}
        self.assertEqual(sorted(rows), [("S1", "0"), ("S1", "90"), ("S1_end_a", "0"), ("S1_end_a", "90")])
        self.assertEqual(rows[("S1_end_a", "90")]["status"], "no_exceedance")
        hw = rows[("S1_end_a", "0")]
        self.assertEqual(hw["status"], "fitted")
        self.assertTrue(0 < int(hw["resamples_used"]) <= 200)
        self.assertLess(float(hw["theta_low"]), float(hw["theta_high"]))
        self.assertGreater(float(hw["theta_sd"]), 0.0)
        first = dict(hw)
        self.assertEqual(self.call(path)[0], 0)
        again = {(r["curve"], r["direction"]): r for r in self.table("fragility_curves.csv")}[("S1_end_a", "0")]
        self.assertEqual(again, first)

        def hardware(rows_first=None, **keys):
            p = self.plan(by_direction=True, resamples=200, **keys)
            body = lines[1:] + extra
            if rows_first:
                body = sorted(body, key=lambda ln: rows_first not in ln)
            with open(os.path.join(self.tmp, "ensemble.csv"), "w") as f:
                f.write("\n".join(lines[:1] + body) + "\n")
            self.assertEqual(self.call(p)[0], 0)
            return {(r["curve"], r["direction"]): r for r in self.table("fragility_curves.csv")}[("S1_end_a", "0")]
        # the same curve whatever the order of the table's rows; another seed gives another spread; a higher
        # confidence a wider interval from the same simulated sets
        self.assertEqual(hardware(rows_first="S1_end_a", seed=11, confidence=0.8), first)
        self.assertNotEqual(hardware(seed=12, confidence=0.8)["theta_sd"], first["theta_sd"])
        wide = hardware(seed=11, confidence=0.98)
        self.assertEqual(wide["theta_sd"], first["theta_sd"])
        self.assertLess(float(wide["theta_low"]), float(first["theta_low"]))
        self.assertGreater(float(wide["theta_high"]), float(first["theta_high"]))
        # pooled over the directions the hardware has 12 members
        path = self.plan(by_direction=False)
        with open(os.path.join(self.tmp, "ensemble.csv"), "w") as f:
            f.write("\n".join(lines + extra) + "\n")
        self.call(path)
        self.assertEqual({r["curve"]: r["members"] for r in self.table("fragility_curves.csv")},
                         {"S1": "12", "S1_end_a": "12"})

    def test_members_without_a_recorded_wind_are_left_out_and_named(self):
        path = self.plan()
        with open(os.path.join(self.tmp, "ensemble.csv")) as f:
            lines = f.read().splitlines()
        c = lines[1].split(",")
        c[15] = ""                      # max_speed of the first member's conductor: a run shorter than the window
        lines[1] = ",".join(c)
        with open(os.path.join(self.tmp, "ensemble.csv"), "w") as f:
            f.write("\n".join(lines) + "\n")
        code, out, _ = self.call(path)
        self.assertEqual(code, 0)
        self.assertIn("NOTE S1 (direction all): 1 member(s) recorded no max_speed and are left out: m001_s20_d0_seed1", out)
        self.assertEqual({r["curve"]: r["members"] for r in self.table("fragility_curves.csv")},
                         {"S1": "5", "S1_end_a": "6"})
        # the 20 m/s stripe of S1 keeps the member that recorded 20.5 m/s: the mean is the recorded wind's
        first = self.table("fragility_points.csv")[0]
        self.assertEqual((first["curve"], first["stripe"], first["members"], first["mean_intensity"]),
                         ("S1", "20", "1", "20.5"))


# Twelve members, four seeds at each of 20, 30 and 40 m/s. Each member's line S1 has the peak ratio
# below, its end S1_end_a 0.5; lines released and towers lost are counted from its logs.
#   20 m/s: ratios 0.8 1.2 1.0 1.1   released 0 0 0 0   towers 0 0 0 0
#   30 m/s: ratios 1.3 1.5 1.4 0.95  released 1 2 0 0   towers 0 1 0 0
#   40 m/s: ratios 1.6 1.8 1.7 1.9   released 1 3 2 4   towers 1 3 0 2
# A ratio of exactly 1 is not an exceedance (the model's rule is demand > capacity). Exceedances per
# stripe, with cascade_towers = 2:
#   DS1 (ratio > 1)      2 3 4      DS2 (a line released)  0 2 4
#   DS3 (a tower lost)   0 1 3      DS4 (2 or more towers) 0 0 2
RATIOS = [0.8, 1.2, 1.0, 1.1, 1.3, 1.5, 1.4, 0.95, 1.6, 1.8, 1.7, 1.9]
RELEASED = [0, 0, 0, 0, 1, 2, 0, 0, 1, 3, 2, 4]
TOWERS = [0, 0, 0, 0, 0, 1, 0, 0, 1, 3, 0, 2]
EXPECTED = {"DS1": (2, 3, 4), "DS2": (0, 2, 4), "DS3": (0, 1, 3), "DS4": (0, 0, 2)}


def member_name(i):
    return "m%03d_s%d_d0_seed%d" % (i + 1, 20 + 10 * (i // 4), i % 4 + 1)


def damage_text(breaks="yes", towers="yes", released=RELEASED, ratios=RATIOS):
    rows = ["member,stripe,nominal_direction,seed,max_peak_ratio,components_exceeded,failures_log,failures,"
            "breaks_log,lines_released,releases,tower_log,towers_lost"]
    for i, r in enumerate(ratios):
        rows.append("%s,%d,0,%d,%s,%d,yes,%d,%s,%s,%s,%s,%s" % (
            member_name(i), 20 + 10 * (i // 4), i % 4 + 1, r, 1 if r > 1 else 0, 1 if r > 1 else 0, breaks,
            released[i] if breaks == "yes" else "", released[i] if breaks == "yes" else "", towers,
            TOWERS[i] if towers == "yes" else ""))
    return "\n".join(rows) + "\n"


def demand_text(ratios=RATIOS):
    rows = ["member,stripe,nominal_direction,seed,kind,component,peak_ratio,peak_speed,max_speed"]
    for i, r in enumerate(ratios):
        x = 20 + 10 * (i // 4)
        rows.append("%s,%d,0,%d,conductor,S1,%s,%d,%d" % (member_name(i), x, i % 4 + 1, r, x, x))
        rows.append("%s,%d,0,%d,hardware,S1_end_a,0.5,%d,%d" % (member_name(i), x, i % 4 + 1, x, x))
    return "\n".join(rows) + "\n"


class DamageStates(Case):
    def damage_plan(self, damage=None, demand=None, **keys):
        for name, text in (("demand.csv", demand or demand_text()), ("damage.csv", damage or damage_text())):
            with open(os.path.join(self.tmp, name), "w") as f:
                f.write(text)
        p = {"demand_table": "demand.csv", "resamples": 0,
             "damage_states": {"table": "damage.csv", "states": ["DS1", "DS2", "DS3", "DS4"], "cascade_towers": 2}}
        p.update(keys)
        p = {k: v for k, v in p.items() if v is not None}
        path = os.path.join(self.tmp, "plan.json")
        with open(path, "w") as f:
            json.dump(p, f)
        return path

    def test_each_state_counts_the_members_that_reach_it_or_a_higher_one(self):
        code, out, err = self.call(self.damage_plan(groups={"network": {"components": ".*"}}))
        self.assertEqual((code, err), (0, ""), out)
        curves = {r["curve"]: r for r in self.table("fragility_curves.csv")}
        points = self.table("fragility_points.csv")
        for state, counts in EXPECTED.items():
            got = [(p["stripe"], p["members"], p["exceedances"]) for p in points if p["curve"] == state]
            self.assertEqual(got, [("20", "4", str(counts[0])), ("30", "4", str(counts[1])),
                                   ("40", "4", str(counts[2]))], state)
            c = curves[state]
            self.assertEqual((c["kind"], c["intensity"], c["members"], c["direction"]),
                             ("damage_state", "stripe", "12", "all"))
            # the curve is the fit of the outcomes above, at the stripes
            data = [(x, 1.0 if k < n else 0.0) for x, n in zip((20.0, 30.0, 40.0), counts) for k in range(4)]
            f = ff.fit(data)
            self.assertEqual(c["status"], f["status"], state)
            if f["status"] == "fitted":
                self.assertAlmostEqual(float(c["theta"]), f["theta"], places=6)
                self.assertAlmostEqual(float(c["beta"]), f["beta"], places=6)
        # DS1 is the group of every component of the demand table (no capacity distribution): the same curve
        for col in ("status", "exceedances", "theta", "beta", "log_likelihood"):
            self.assertEqual(curves["DS1"][col], curves["network"][col], col)
        self.assertIn("NOTE DS4 is 2 or more towers lost (cascade_towers, the user's choice)", out)
        self.assertNotIn("reach a damage state without", out)

    def test_members_the_driver_did_not_collect_are_named(self):
        status = os.path.join(self.tmp, "members_status.csv")
        header = "member,stripe,nominal_direction,seed,status,returncode,demand_rows,intensity_rows,detail\n"
        rows = "".join("%s,20,0,1,complete,0,2,1,\n" % member_name(i) for i in range(12))
        with open(status, "w") as f:
            f.write(header + rows + "m013_s40_d0_seed5,40,0,5,failed,1,,,amrex::Abort line S1 left the domain\n")
        code, out, err = self.call(self.damage_plan())
        self.assertEqual((code, err), (0, ""), out)
        self.assertIn("NOTE damage_states: 1 of 13 member(s) of the ensemble are not in the tables "
                      "(m013_s40_d0_seed5 (failed))", out)
        with open(status, "w") as f:
            f.write(header + rows)
        code, out, _ = self.call(self.damage_plan())
        self.assertEqual(code, 0)
        self.assertNotIn("not in the tables", out)
        self.assertNotIn("members_status.csv", out)
        os.remove(status)
        code, out, _ = self.call(self.damage_plan())
        self.assertIn("NOTE damage_states: no members_status.csv beside", out)

    def test_a_cascade_threshold_moves_members_out_of_ds4(self):
        # with cascade_towers = 3 only the 40 m/s member with three towers lost is in DS4
        code, _, _ = self.call(self.damage_plan(damage_states={"table": "damage.csv", "states": ["DS3", "DS4"],
                                                               "cascade_towers": 3}))
        self.assertEqual(code, 0)
        points = self.table("fragility_points.csv")
        self.assertEqual([p["exceedances"] for p in points if p["curve"] == "DS4"], ["0", "0", "1"])
        self.assertEqual([p["exceedances"] for p in points if p["curve"] == "DS3"], ["0", "1", "3"])
        self.assertFalse(any(p["curve"] in ("DS1", "DS2") for p in points))

    def test_a_member_that_skips_a_state_is_named_and_counts_below(self):
        # m008 releases a line with its peak ratio at 0.95: it is in DS2, so it counts in DS1 too
        released = list(RELEASED)
        released[7] = 1
        code, out, _ = self.call(self.damage_plan(damage=damage_text(released=released)))
        self.assertEqual(code, 0)
        self.assertIn("NOTE 1 member(s) reach a damage state without every lower one (m008_s30_d0_seed4)", out)
        points = self.table("fragility_points.csv")
        self.assertEqual([p["exceedances"] for p in points if p["curve"] == "DS1"], ["2", "4", "4"])
        self.assertEqual([p["exceedances"] for p in points if p["curve"] == "DS2"], ["0", "3", "4"])

    def test_a_state_does_not_depend_on_the_other_states_listed(self):
        # m008 (ratio 0.95, a line released) counts in DS1 whether or not DS2 is listed: the states are
        # cumulative over all four, read from every log the runs wrote
        released = list(RELEASED)
        released[7] = 1
        for states in (["DS1"], ["DS1", "DS3"], ["DS1", "DS2"]):
            code, out, _ = self.call(self.damage_plan(damage=damage_text(released=released),
                                                      damage_states={"table": "damage.csv", "states": states}))
            self.assertEqual(code, 0, states)
            self.assertIn("(m008_s30_d0_seed4)", out)
            points = self.table("fragility_points.csv")
            self.assertEqual([p["exceedances"] for p in points if p["curve"] == "DS1"], ["2", "4", "4"], states)

    def test_only_ds1_needs_no_failure_log(self):
        # a capacity distribution changes the component curves only: DS1 is the runs' own outcome
        code, out, err = self.call(self.damage_plan(damage=damage_text(breaks="no", towers="no"),
                                                    damage_states={"table": "damage.csv", "states": ["DS1"]},
                                                    capacity={"conductor": {"cov": 0.1}}))
        self.assertEqual((code, err), (0, ""))
        self.assertIn("NOTE the damage states are the runs' own outcomes, at the nominal capacities; capacity applies "
                      "to the component and group curves only", out)
        self.assertEqual([p["exceedances"] for p in self.table("fragility_points.csv") if p["curve"] == "DS1"],
                         ["2", "3", "4"])

    def test_tables_from_other_members_are_refused(self):
        lines = damage_text().splitlines(True)
        cases = (("1 member(s) of the other table are not in the damage table: m012_s40_d0_seed4",
                  "".join(lines[:-1])),
                 ("1 member(s) of the other table are not in the demand table: m099",
                  "".join(lines) + lines[1].replace("m001_s20_d0_seed1", "m099")),
                 ("member(s) m002_s20_d0_seed2 differ between the two tables in stripe, nominal_direction or seed",
                  "".join(lines).replace("m002_s20_d0_seed2,20,0,2,", "m002_s20_d0_seed2,20,0,7,")),
                 ("member(s) m005_s30_d0_seed1 differ between the two tables",
                  "".join(lines).replace("m005_s30_d0_seed1,30,0,", "m005_s30_d0_seed1,35,0,")))
        for expect, text in cases:
            code, _, err = self.call(self.damage_plan(damage=text))
            self.assertEqual(code, 2, expect)
            self.assertIn("damage_states: " if "member(s) m0" in expect else "damage_states: the two tables are "
                          "from different members: ", err)
            self.assertIn(expect, err)

    def test_refused(self):
        def refuse(expect, damage=None, **keys):
            code, _, err = self.call(self.damage_plan(damage=damage, **keys))
            self.assertEqual(code, 2, expect)
            self.assertIn(expect, err)

        def ds(**d):
            base = {"table": "damage.csv", "states": ["DS1", "DS2", "DS3", "DS4"], "cascade_towers": 2}
            base.update(d)
            return {k: v for k, v in base.items() if v is not None}

        refuse("damage_states: DS4 (a cascade over N towers) needs cascade_towers = N; it has no default, since N is "
               "the user's choice", damage_states=ds(cascade_towers=None))
        refuse("damage_states: cascade_towers is used by DS4 only, which states does not list",
               damage_states=ds(states=["DS1", "DS3"]))
        for n in (1, 0, 2.5, "3", True, 100001):
            refuse("damage_states: cascade_towers must be an integer in [2, 100000] (one tower lost is DS3), got",
                   damage_states=ds(cascade_towers=n))
        for states in ([], ["DS5"], ["DS1", "DS1"], "DS1"):
            refuse("damage_states: states must list some of DS1, DS2, DS3, DS4, each once",
                   damage_states=ds(states=states))
        refuse("damage_states must be {", damage_states=ds(level=1))
        refuse("damage_states must be {", damage_states=ds(table=None))
        refuse("damage_states must be {", damage_states=["DS1"])
        refuse("damage_states: table must name the gathered damage table", damage_states=ds(table=" "))
        refuse("damage_states: DS2 needs every member's record, but 12 member(s) have no breaks.dat (its option was "
               "off): m001_s20_d0_seed1, m002_s20_d0_seed2, m003_s20_d0_seed3, m004_s20_d0_seed4, m005_s30_d0_seed1, "
               "...", damage=damage_text(breaks="no"))
        refuse("damage_states: DS3 needs every member's record, but 12 member(s) have no tower_failures.dat",
               damage=damage_text(towers="no"), damage_states=ds(states=["DS3"], cascade_towers=None))
        refuse("damage_states: DS1 is also the name of a component or group", groups={"DS1": {"kinds": ["conductor"]}})
        refuse("damage_states: ", damage=damage_text().replace("towers_lost", "towers"))
        refuse("damage_states: line 2 of", damage=damage_text().replace(",yes,0\n", ",maybe,0\n", 1))
        refuse("has a negative towers_lost", damage=damage_text().replace(",yes,0\n", ",yes,-1\n", 1))
        refuse("is listed twice", damage=damage_text() + damage_text().splitlines()[1] + "\n")


class Refused(Case):
    def test_plan_keys(self):
        for expect, keys in (
                ("demand_table must name the gathered demand table", {"demand_table": ""}),
                ("demand_table: cannot read", {"demand_table": "missing.csv"}),
                ("output_dir must name a directory", {"output_dir": ""}),
                ("intensity must be one of max_speed, peak_speed, stripe, got \"gust\"", {"intensity": "gust"}),
                ("by_direction must be true or false, got 1", {"by_direction": 1}),
                ("capacity: 'pole' is not a kind of component", {"capacity": {"pole": {"cov": 0.1}}}),
                ("capacity: member needs {\"cov\"", {"capacity": {"member": {"mean_factor": 1.0}}}),
                ("capacity: member needs {\"cov\"", {"capacity": {"member": {"cov": 0.1, "median": 1.0}}}),
                ("capacity: member.cov must be in [0.001, 1], got 0", {"capacity": {"member": {"cov": 0}}}),
                ("capacity: member.cov must be in [0.001, 1], got 1e-09", {"capacity": {"member": {"cov": 1e-9}}}),
                ("capacity: member.cov must be in [0.001, 1], got 1.5", {"capacity": {"member": {"cov": 1.5}}}),
                ("capacity: member.mean_factor must be in [0.1, 10], got 0.05",
                 {"capacity": {"member": {"cov": 0.1, "mean_factor": 0.05}}}),
                ("capacity: member.mean_factor must be in [0.1, 10], got 11",
                 {"capacity": {"member": {"cov": 0.1, "mean_factor": 11}}}),
                ("groups: line needs \"components\"", {"groups": {"line": {}}}),
                ("groups: line.components must be a regular expression, got \"(\"", {"groups": {"line": {"components": "("}}}),
                ("groups: line.kinds must be a list of", {"groups": {"line": {"kinds": ["pole"]}}}),
                ("groups: line matches no component of the table", {"groups": {"line": {"components": "^L9"}}}),
                ("groups: the group name 'S1' is also a component's name", {"groups": {"S1": {"components": "^S1"}}}),
                ("resamples must be an integer in [0, 100000], got -1", {"resamples": -1}),
                ("resamples must be an integer in [0, 100000], got 100001", {"resamples": 100001}),
                ("resamples must be an integer in [0, 100000], got 2.5", {"resamples": 2.5}),
                ("seed must be an integer in [0, 2147483647], got -3", {"seed": -3}),
                ("confidence must be in [0.5, 1), got 1", {"confidence": 1}),
                ("confidence must be in [0.5, 1), got 0.4", {"confidence": 0.4}),
                ("--plan: unknown key(s) speeds;", {"speeds": [10]})):
            with self.subTest(expect=expect):
                code, _, err = self.call(self.plan(**keys))
                self.assertEqual(code, 2, err)
                self.assertIn(expect, err)
                self.assertFalse(os.path.exists(os.path.join(self.tmp, "fragility")), "a refused plan writes nothing")

    def test_the_table(self):
        path = self.plan()
        with open(os.path.join(self.tmp, "ensemble.csv")) as f:
            good = f.read()
        lines = good.splitlines()
        for expect, text in (
                ("holds no rows", lines[0] + "\n"),
                ("lacks the column(s) peak_ratio", good.replace(",peak_ratio,", ",ratio,")),
                ("line 2 of", good.replace(",S1,0.8,", ",S1,high,")),
                ("has no number in peak_ratio", good.replace(",S1,0.8,", ",S1,high,")),
                ("is not finite in peak_ratio", good.replace(",S1,0.8,", ",S1,inf,")),
                ("has a negative peak_ratio", good.replace(",S1,0.8,", ",S1,-0.8,")),
                ("has the unknown kind 'pole'", good.replace("conductor,S1,0.8", "pole,S1,0.8")),
                ("line 2 of", lines[0] + "\n" + ",".join(lines[1].split(",")[:7]) + "\n"),
                ("does not have the header's number of fields",
                 lines[0] + "\n" + ",".join(lines[1].split(",")[:7]) + "\n"),
                ("has a max_speed that is not positive and finite", good.replace(",19.5,0,,8,0,", ",-3,0,,8,0,", 1)),
                ("the rows of member m001_s20_d0_seed1 do not agree in stripe",
                 good.replace("m001_s20_d0_seed1,20,0,1,hardware", "m001_s20_d0_seed1,40,0,1,hardware")),
                # a field longer than the csv module reads (131072 characters) on every Python version
                ("is not a readable CSV table", good.replace("conductor,S1,0.8", "conductor," + "S" * 200000 + ",0.8")),
                ("holds no row with a peak_ratio", lines[0] + "\n" + lines[1].replace(",S1,0.8,", ",S1,,") + "\n"),
                ("member m001_s20_d0_seed1 lists the component S1 twice", good + lines[1] + "\n")):
            with self.subTest(expect=expect):
                with open(os.path.join(self.tmp, "ensemble.csv"), "w") as f:
                    f.write(text)
                code, _, err = self.call(path)
                self.assertEqual(code, 2, err)
                self.assertIn("demand_table:", err)
                self.assertIn(expect, err)


    def test_rows_left_out_and_groups_that_do_not_hold_every_component(self):
        path = self.plan(groups={"line": {"components": "^S1"}})
        with open(os.path.join(self.tmp, "ensemble.csv")) as f:
            good = f.read()

        def with_table(text):
            with open(os.path.join(self.tmp, "ensemble.csv"), "w", encoding="utf-8") as f:
                f.write(text)
            return self.call(path)
        # a byte-order mark, as a spreadsheet writes it, is read
        self.assertEqual(with_table("\ufeff" + good)[0], 0)
        # a component the run never recorded has an empty ratio: the row is left out and counted, and the group,
        # which needs every component in every member, is refused by name
        code, _, err = with_table(good.replace(",S1_end_a,0.7,", ",S1_end_a,,"))
        self.assertEqual(code, 2)
        self.assertIn("groups: member m001_s20_d0_seed1 lacks 1 of the 2 components of the group line (first: S1_end_a)",
                      err)
        path = self.plan()
        code, out, _ = with_table(good.replace(",S1_end_a,0.7,", ",S1_end_a,,"))
        self.assertEqual(code, 0)
        self.assertIn("NOTE 1 row(s) have no peak_ratio (a component the run never recorded) and are left out", out)
        self.assertEqual({r["curve"]: r["members"] for r in self.table("fragility_curves.csv")},
                         {"S1": "6", "S1_end_a": "5"})
        # an 'exceeded' flag the ratio does not bear out is named; the ratio decides
        code, out, _ = with_table(good.replace(",S1,0.9,", ",S1,1,").replace("20.5,0,,8,0,,,,,,,", "20.5,0,,8,1,,,,,,,", 1))
        self.assertEqual(code, 0)
        self.assertIn("NOTE 1 row(s) have an 'exceeded' flag that differs from peak_ratio > 1", out)
        self.assertEqual({r["curve"]: r["exceedances"] for r in self.table("fragility_curves.csv")}["S1"], "3")
        # an earlier run's output is replaced whole
        self.assertEqual(sorted(os.listdir(os.path.join(self.tmp, "fragility"))),
                         ["fragility_curves.csv", "fragility_points.csv"])


if __name__ == "__main__":
    unittest.main(verbosity=1)
