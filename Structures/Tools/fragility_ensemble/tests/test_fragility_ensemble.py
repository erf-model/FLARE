#!/usr/bin/env python3
"""Tests of fragility_ensemble.py (standard library only; no ERF run).

    python3 test_fragility_ensemble.py            every test
    python3 test_fragility_ensemble.py DryRun     one class

Every expected deck line, wind-file row, manifest entry and collected row below is worked out by
hand from the stored base decks (base_inflow, base_periodic, base_prescribed) and the stored member
tables (collect); nothing is regenerated with the driver.

base_inflow: the inflow profile has u = 6, 10, 14 m/s at z = 0, 20, 100 m, so 8 m/s along +x at the
reference height of 10 m. A 10 m/s stripe multiplies the wind by 10/8 = 1.25 and a 20 m/s stripe by
2.5; xlo.KE = 2.0 becomes 2.0 * 1.25^2 = 3.125 and 2.0 * 2.5^2 = 12.5. Towards -135 degrees the
2.5-fold wind has u = v = -2.5 cos(45 deg) u0 = -1.767766953 u0: -10.60660172, -17.67766953 and
-24.74873734 m/s. The profile's top row carries w = 0.4 m/s, which becomes 0.5 and 1 m/s. The
mesh is 50 m by 50 m and the step 0.5 s (0.25 s for the 20 m/s stripe), so the largest
(|u|/dx + |v|/dy) dt is 17.5/50 x 0.5 = 0.175 along an axis and 2 x 12.37436867/50 x 0.5 =
0.2474873734 towards -135 degrees, for both stripes.
"""
import contextlib
import csv
import importlib.util
import io
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import time
import unittest

sys.dont_write_bytecode = True     # the driver is loaded from the source tree
HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("fragility_ensemble", os.path.join(HERE, "..", "fragility_ensemble.py"))
fe = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fe)


def call(argv, gate=None):
    """(return code, stdout, stderr) of the driver's main()."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = fe.main(argv, gate=gate)
    return code, out.getvalue(), err.getvalue()


def read(path):
    with open(path) as f:
        return f.read()


class Case(unittest.TestCase):
    maxDiff = None

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="fragility_ensemble_")
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def stage(self, name):
        self.staged = getattr(self, "staged", 0) + 1
        d = os.path.join(self.tmp, "%d" % self.staged, name)
        shutil.copytree(os.path.join(HERE, name), d)
        return d

    def edit_plan(self, d, **changes):
        path = os.path.join(d, "plan.json")
        plan = json.loads(read(path))
        for k, v in changes.items():
            if v is None:
                plan.pop(k, None)
            else:
                plan[k] = v
        with open(path, "w") as f:
            json.dump(plan, f)
        return path

    def changed_lines(self, base_text, member_text):
        """(lines only in the base deck, lines only in the member's), each in its own order."""
        base, member = base_text.splitlines(), member_text.splitlines()
        return [ln for ln in base if ln not in member], [ln for ln in member if ln not in base]


# the base deck with lines.inputs read in place, as every member's deck carries it
INFLOW_FLAT = """\
# A deck for the fragility ensemble's tests: it is rewritten, never run. One span in a wind
# that enters through x-low as the profile of inflow_profile (8 m/s along +x at 10 m).
max_step = 10
erf.fixed_dt = 0.5          # s
geometry.prob_lo     = 0.    0.    0.
geometry.prob_hi     = 1500. 1000. 300.
amr.n_cell           = 30    20    12
geometry.is_periodic = 0 0 0
xlo.type           = "Inflow"
xlo.dirichlet_file = "inflow_profile"   # z u v w theta
xlo.KE             = 2.0
xhi.type           = "Outflow"
ylo.type           = "SlipWall"
yhi.type           = "SlipWall"
zlo.type           = "SlipWall"
zhi.type           = "SlipWall"
erf.init_type           = "input_sounding"
erf.input_sounding_file = "input_sounding"
erf.terrain_file_name   = "terrain.txt"
conductors.diagnostics_dir = conductors
# fragility_ensemble: the lines of lines.inputs
# the span and its gusts
erf.structures.model = conductors
conductors.lines     = S1
conductors.gust_type = random
conductors.gust_seed = 7
"""

INFLOW_MEMBERS = [
    "m001_s10_d0_seed3", "m002_s10_d0_seed8", "m003_s10_d90_seed3", "m004_s10_d90_seed8",
    "m005_s10_dm135_seed3", "m006_s10_dm135_seed8", "m007_s20_d0_seed3", "m008_s20_d0_seed8",
    "m009_s20_d90_seed3", "m010_s20_d90_seed8", "m011_s20_dm135_seed3", "m012_s20_dm135_seed8"]


class DryRun(Case):
    def run_inflow(self):
        d = self.stage("base_inflow")
        code, out, err = call(["run", "--plan", os.path.join(d, "plan.json"), "--dry-run"])
        self.assertEqual((code, err), (0, ""))
        return d, out

    def test_the_member_set_is_stripes_by_directions_by_seeds(self):
        d, out = self.run_inflow()
        self.assertEqual(sorted(os.listdir(os.path.join(d, "out", "members"))), INFLOW_MEMBERS)
        self.assertEqual(sorted(os.listdir(os.path.join(d, "out"))), ["manifest.json", "members"])
        self.assertIn("12 members", out)
        self.assertIn("8 m/s towards 0 degrees at 10 m in inflow_profile", out)
        for m in INFLOW_MEMBERS:   # the deck, the two wind files and the support file; nothing was run
            self.assertEqual(sorted(os.listdir(os.path.join(d, "out", "members", m))),
                             ["inflow_profile", "input_sounding", "span.i", "terrain.txt"])
            self.assertEqual(read(os.path.join(d, "out", "members", m, "terrain.txt")),
                             read(os.path.join(d, "terrain.txt")))

    def test_a_stripe_along_the_base_wind_changes_only_the_named_inputs(self):
        d, _ = self.run_inflow()
        m = os.path.join(d, "out", "members", "m001_s10_d0_seed3")
        want = INFLOW_FLAT.replace("xlo.KE             = 2.0", "xlo.KE = 3.125") \
                          .replace("conductors.gust_seed = 7", "conductors.gust_seed = 3") + (
            "\n# set by fragility_ensemble\n"
            "conductors.intensity_measure = wind\n"
            "conductors.S1.rated_breaking_strength = 50000\n")
        self.assertEqual(read(os.path.join(m, "span.i")), want)
        self.assertEqual(read(os.path.join(m, "input_sounding")),
                         "1000.0 300.0 0.0\n0.0 300.0 0.0 7.5 0\n20.0 300.0 0.0 12.5 0\n100.0 300.0 0.0 17.5 0\n")
        self.assertEqual(read(os.path.join(m, "inflow_profile")),
                         "0.0 7.5 0 0 300.0\n20.0 12.5 0 0 300.0\n100.0 17.5 0 0.5 300.0\n")
        # the other seed differs in the seed's line only
        other = read(os.path.join(d, "out", "members", "m002_s10_d0_seed8", "span.i"))
        self.assertEqual(other, want.replace("conductors.gust_seed = 3", "conductors.gust_seed = 8"))

    def test_a_wind_along_y_moves_the_inflow_to_the_y_faces(self):
        d, _ = self.run_inflow()
        m = os.path.join(d, "out", "members", "m003_s10_d90_seed3")
        gone, new = self.changed_lines(INFLOW_FLAT, read(os.path.join(m, "span.i")))
        self.assertEqual(gone, [
            'xlo.type           = "Inflow"',
            'xlo.dirichlet_file = "inflow_profile"   # z u v w theta',
            'xlo.KE             = 2.0',
            'xhi.type           = "Outflow"',
            'ylo.type           = "SlipWall"',
            'yhi.type           = "SlipWall"',
            'conductors.gust_seed = 7'])
        self.assertEqual(new, [
            'xlo.type = "SlipWall"',
            'xhi.type = "SlipWall"',
            'ylo.type = "Inflow"',
            'yhi.type = "Outflow"',
            'conductors.gust_seed = 3',
            '',
            '# set by fragility_ensemble',
            'ylo.dirichlet_file = "inflow_profile"',
            'ylo.KE = 3.125',
            'conductors.intensity_measure = wind',
            'conductors.S1.rated_breaking_strength = 50000'])
        self.assertEqual(read(os.path.join(m, "input_sounding")),
                         "1000.0 300.0 0.0\n0.0 300.0 0.0 0 7.5\n20.0 300.0 0.0 0 12.5\n100.0 300.0 0.0 0 17.5\n")
        self.assertEqual(read(os.path.join(m, "inflow_profile")),
                         "0.0 0 7.5 0 300.0\n20.0 0 12.5 0 300.0\n100.0 0 17.5 0.5 300.0\n")

    def test_an_oblique_wind_enters_through_two_faces_and_a_stripe_carries_its_own_inputs(self):
        d, _ = self.run_inflow()
        m = os.path.join(d, "out", "members", "m011_s20_dm135_seed3")
        gone, new = self.changed_lines(INFLOW_FLAT, read(os.path.join(m, "span.i")))
        self.assertEqual(gone, [
            'max_step = 10',
            'erf.fixed_dt = 0.5          # s',
            'xlo.type           = "Inflow"',
            'xlo.dirichlet_file = "inflow_profile"   # z u v w theta',
            'xlo.KE             = 2.0',
            'xhi.type           = "Outflow"',
            'ylo.type           = "SlipWall"',
            'yhi.type           = "SlipWall"',
            'conductors.gust_seed = 7'])
        self.assertEqual(new, [
            'max_step = 20',
            'erf.fixed_dt = 0.25',
            'xlo.type = "Outflow"',
            'xhi.type = "Inflow"',
            'ylo.type = "Outflow"',
            'yhi.type = "Inflow"',
            'conductors.gust_seed = 3',
            '',
            '# set by fragility_ensemble',
            'xhi.dirichlet_file = "inflow_profile"',
            'xhi.KE = 12.5',
            'yhi.dirichlet_file = "inflow_profile"',
            'yhi.KE = 12.5',
            'conductors.intensity_measure = wind',
            'conductors.S1.rated_breaking_strength = 50000'])
        self.assertEqual(read(os.path.join(m, "input_sounding")),
                         "1000.0 300.0 0.0\n"
                         "0.0 300.0 0.0 -10.60660172 -10.60660172\n"
                         "20.0 300.0 0.0 -17.67766953 -17.67766953\n"
                         "100.0 300.0 0.0 -24.74873734 -24.74873734\n")
        self.assertEqual(read(os.path.join(m, "inflow_profile")),
                         "0.0 -10.60660172 -10.60660172 0 300.0\n"
                         "20.0 -17.67766953 -17.67766953 0 300.0\n"
                         "100.0 -24.74873734 -24.74873734 1 300.0\n")

    def test_the_manifest_lists_every_member_and_its_settings(self):
        d, _ = self.run_inflow()
        man = json.loads(read(os.path.join(d, "out", "manifest.json")))
        members = man.pop("members")
        sources = man.pop("sources")    # a SHA-256 of every file the members were written from
        self.assertEqual(sorted(sources), ["inflow_profile", "input_sounding", "span.i", "terrain.txt"])
        self.assertTrue(all(len(v) == 64 and set(v) <= set("0123456789abcdef") for v in sources.values()))
        self.assertEqual(len(set(sources.values())), 4)
        self.assertEqual(man, {
            "tool": "fragility_ensemble", "format": 1, "base_deck": "span.i", "deck": "span.i",
            "diagnostics_dir": "conductors", "support_files": ["terrain.txt"],
            "wind": {"mode": "sounding", "base_speed": 8.0, "base_direction": 0.0, "reference_height": 10,
                     "reference_file": "inflow_profile", "files": ["inflow_profile", "input_sounding"]}})
        self.assertEqual([m["member"] for m in members], INFLOW_MEMBERS)
        self.assertEqual([m["directory"] for m in members], ["members/" + m for m in INFLOW_MEMBERS])
        self.assertEqual([m["stripe"] for m in members], [10.0] * 6 + [20.0] * 6)
        self.assertEqual([m["nominal_direction"] for m in members], [0.0, 0.0, 90.0, 90.0, -135.0, -135.0] * 2)
        self.assertEqual([m["seed"] for m in members], [3, 8] * 6)
        self.assertEqual([m["scale"] for m in members], [1.25] * 6 + [2.5] * 6)
        self.assertEqual([m["rotation"] for m in members], [0.0, 0.0, 90.0, 90.0, -135.0, -135.0] * 2)
        self.assertEqual([m["max_wind"] for m in members], [17.5] * 6 + [35.0] * 6)
        self.assertEqual([m["courant"] for m in members], [0.175] * 4 + [0.2474873734] * 2 + [0.175] * 4 + [0.2474873734] * 2)
        self.assertEqual(members[0]["inputs"], {
            "xlo.KE": "3.125", "conductors.intensity_measure": "wind", "conductors.gust_seed": "3",
            "conductors.S1.rated_breaking_strength": "50000"})
        self.assertEqual(members[3]["inputs"], {
            "xlo.type": '"SlipWall"', "xlo.dirichlet_file": None, "xlo.KE": None, "xhi.type": '"SlipWall"',
            "ylo.type": '"Inflow"', "ylo.dirichlet_file": '"inflow_profile"', "ylo.KE": "3.125",
            "yhi.type": '"Outflow"', "conductors.intensity_measure": "wind", "conductors.gust_seed": "8",
            "conductors.S1.rated_breaking_strength": "50000"})
        self.assertEqual(members[10]["inputs"], {
            "xlo.type": '"Outflow"', "xlo.dirichlet_file": None, "xlo.KE": None,
            "xhi.type": '"Inflow"', "xhi.dirichlet_file": '"inflow_profile"', "xhi.KE": "12.5",
            "ylo.type": '"Outflow"', "yhi.type": '"Inflow"', "yhi.dirichlet_file": '"inflow_profile"',
            "yhi.KE": "12.5", "conductors.intensity_measure": "wind", "conductors.gust_seed": "3",
            "conductors.S1.rated_breaking_strength": "50000", "erf.fixed_dt": "0.25", "max_step": "20"})

    def test_a_periodic_deck_turns_its_sounding_only(self):
        # v = 15 and 45 m/s at 0 and 300 m: 16 m/s along +y at 10 m; 32 m/s towards +x doubles and turns it
        d = self.stage("base_periodic")
        code, _, err = call(["run", "--plan", os.path.join(d, "plan.json"), "--dry-run"])
        self.assertEqual((code, err), (0, ""))
        m = os.path.join(d, "out", "members", "m001_s32_d0")
        self.assertEqual(read(os.path.join(m, "span.i")), read(os.path.join(d, "span.i")) + (
            "\n# set by fragility_ensemble\nconductors.intensity_measure = wind\n"))
        self.assertEqual(read(os.path.join(m, "input_sounding")),
                         "1000.0 300.0 0.0\n0.0 300.0 0.0 30 0\n300.0 300.0 0.0 90 0\n")
        man = json.loads(read(os.path.join(d, "out", "manifest.json")))
        self.assertEqual(sorted(man["sources"]), ["input_sounding", "span.i"])
        self.assertEqual(man["wind"], {"mode": "sounding", "base_speed": 16.0, "base_direction": 90.0,
                                       "reference_height": 10.0, "reference_file": "input_sounding",
                                       "files": ["input_sounding"]})
        self.assertEqual(man["members"], [{
            "member": "m001_s32_d0", "stripe": 32.0, "nominal_direction": 0.0, "seed": None,
            "directory": "members/m001_s32_d0", "scale": 2.0, "rotation": -90.0,
            "max_wind": 90.0, "courant": None,      # 2 x 45 m/s at 300 m; the deck gives no mesh or step
            "inputs": {"conductors.intensity_measure": "wind"}}])

    def test_without_directions_the_base_direction_is_kept_and_the_base_stripe_keeps_the_files(self):
        d = self.stage("base_periodic")
        plan = self.edit_plan(d, directions=None, speeds=[16, 24])
        self.assertEqual(call(["run", "--plan", plan, "--dry-run"])[0], 0)
        base = os.path.join(d, "out", "members", "m001_s16_d90")
        self.assertEqual(read(os.path.join(base, "input_sounding")), read(os.path.join(d, "input_sounding")))
        self.assertEqual(read(os.path.join(d, "out", "members", "m002_s24_d90", "input_sounding")),
                         "1000.0 300.0 0.0\n0.0 300.0 0.0 0 22.5\n300.0 300.0 0.0 0 67.5\n")

    def test_a_prescribed_wind_is_scaled_and_turned_in_the_deck(self):
        # 15 m/s along +y, 1.5 m/s upwards: 30 m/s towards 180 degrees is (-30, 0), towards 45 degrees
        # 30 cos 45 = 21.21320344 in u and v; the vertical component doubles with the speed
        d = self.stage("base_prescribed")
        code, _, err = call(["run", "--plan", os.path.join(d, "plan.json"), "--dry-run"])
        self.assertEqual((code, err), (0, ""))
        base = read(os.path.join(d, "span.i"))
        for name, wind in (("m001_s30_d180", "-30 0 3"), ("m002_s30_d45", "21.21320344 21.21320344 3")):
            m = os.path.join(d, "out", "members", name)
            self.assertEqual(os.listdir(m), ["span.i"])
            self.assertEqual(read(os.path.join(m, "span.i")),
                             base.replace("conductors.prescribed_velocity = 0. 15. 1.5",
                                          "conductors.prescribed_velocity = " + wind)
                             + "\n# set by fragility_ensemble\nconductors.intensity_measure = wind\n")

    def test_a_second_dry_run_keeps_the_directories_and_another_plan_is_refused(self):
        d = self.stage("base_inflow")
        plan = os.path.join(d, "plan.json")
        self.assertIn("(12 directories written now)", call(["run", "--plan", plan, "--dry-run"])[1])
        marker = os.path.join(d, "out", "members", "m001_s10_d0_seed3", "kept")
        open(marker, "w").close()
        self.assertIn("(0 directories written now)", call(["run", "--plan", plan, "--dry-run"])[1])
        self.assertTrue(os.path.exists(marker))
        self.edit_plan(d, seeds=[3, 9])
        code, _, err = call(["run", "--plan", plan, "--dry-run"])
        self.assertEqual(code, 2)
        self.assertIn("output_dir:", err)
        self.assertIn("was written from another plan", err)

    def test_a_changed_source_file_or_a_lost_manifest_is_refused(self):
        d = self.stage("base_inflow")
        plan = os.path.join(d, "plan.json")
        self.assertEqual(call(["run", "--plan", plan, "--dry-run"])[0], 0)
        manifest = read(os.path.join(d, "out", "manifest.json"))
        # the base deck, an included file, a wind file's theta and a support file: none changes a member's listed
        # inputs, and each would leave the member directories behind the files they were written from
        for name, old, new in (("span.i", "max_step = 10", "max_step = 11"),
                               ("lines.inputs", "conductors.lines     = S1", "conductors.lines     = S2"),
                               ("input_sounding", "100.0 300.0", "100.0 301.0"),
                               ("inflow_profile", "0.4 300.0", "0.4 301.0"),
                               ("terrain.txt", "0. 0. 0. 0.", "0. 0. 0. 1.")):
            with self.subTest(name=name):
                text = read(os.path.join(d, name))
                self.assertIn(old, text)
                with open(os.path.join(d, name), "w") as f:
                    f.write(text.replace(old, new))
                code, _, err = call(["run", "--plan", plan, "--dry-run"])
                self.assertEqual(code, 2)
                self.assertIn("was written from another plan, base deck, wind file or support file", err)
                with open(os.path.join(d, name), "w") as f:
                    f.write(text)
                self.assertEqual(call(["run", "--plan", plan, "--dry-run"])[0], 0)
        self.assertEqual(read(os.path.join(d, "out", "manifest.json")), manifest)
        os.remove(os.path.join(d, "out", "manifest.json"))
        code, _, err = call(["run", "--plan", plan, "--dry-run"])
        self.assertEqual(code, 2)
        self.assertIn("holds member directories but no manifest.json", err)

    def test_the_notes_name_a_stripe_too_fast_for_the_step_and_an_unused_inflow_k(self):
        # 100 m/s at 10 m is 12.5 times the base wind: 175 m/s at the top is refused; 60 m/s gives 105 m/s there
        # and 105/50 x 0.5 = 1.05 along x
        d = self.stage("base_inflow")
        plan = self.edit_plan(d, speeds=[10, 60], directions=[0], seeds=[3])
        code, out, _ = call(["run", "--plan", plan, "--dry-run"])
        self.assertEqual(code, 0)
        self.assertIn("NOTE m002_s60_d0_seed3: (|u|/dx + |v|/dy) erf.fixed_dt reaches 1.05 on level 0 "
                      "(largest wind 105 m/s)", out)
        self.assertNotIn("NOTE m001", out)
        self.assertIn("NOTE the inflow face's KE is scaled with the stripe, but ERF uses it only together with "
                      "<face>.density", out)


class Refused(Case):
    """Every option's range check stops the driver with a message that names the option."""

    def refuse(self, fixture, expect, args=(), **plan_changes):
        d = self.stage(fixture)
        plan = self.edit_plan(d, **plan_changes)
        code, _, err = call(["run", "--plan", plan, "--dry-run"] + list(args))
        self.assertEqual(code, 2, err)
        self.assertIn(expect, err)
        self.assertFalse(os.path.exists(os.path.join(d, "out")), "a refused plan must write nothing")

    def test_plan_options(self):
        for expect, changes in (
                ("speeds must be a list of at least one wind speed", {"speeds": []}),
                ("speeds must be in (0, 150] m/s, got 0", {"speeds": [0]}),
                ("speeds must be in (0, 150] m/s, got 151", {"speeds": [10, 151]}),
                ("speeds must be in (0, 150] m/s, got \"fast\"", {"speeds": ["fast"]}),
                ("speeds: 10 m/s is listed twice", {"speeds": [10, 10.0]}),
                ("speeds: unknown key(s) dt in a stripe", {"speeds": [{"speed": 10, "dt": 1}]}),
                ("directions must be in (-180, 180] degrees", {"directions": [-180]}),
                ("directions must be in (-180, 180] degrees", {"directions": [270]}),
                ("directions: 90 degrees is listed twice", {"directions": [90, 90]}),
                ("directions must be a list of at least one direction", {"directions": []}),
                ("seeds must be integers in [0, 2147483647], got -1", {"seeds": [-1]}),
                ("seeds must be integers in [0, 2147483647], got 1.5", {"seeds": [1.5]}),
                ("seeds: a seed is listed twice", {"seeds": [4, 4]}),
                ("reference_height must be in (0, 500] m, got 0", {"reference_height": 0}),
                ("reference_height must be in (0, 500] m, got 501", {"reference_height": 501}),
                ("reference_height must be in [0, 100] m, the heights of", {"reference_height": 120}),
                ("side_type must be one of SlipWall, NoSlipWall, Symmetry, got \"Open\"", {"side_type": "Open"}),
                ("inputs must be an object of deck inputs", {"inputs": [1]}),
                ("inputs: 'a b' is not a deck key", {"inputs": {"a b": 1}}),
                ("inputs: the value of erf.x must be", {"inputs": {"erf.x": {}}}),
                ("inputs: conductors.intensity_measure cannot be given", {"inputs": {"conductors.intensity_measure": "none"}}),
                ("inputs: conductors.prescribed_velocity cannot be given", {"inputs": {"conductors.prescribed_velocity": [1, 0, 0]}}),
                ("inputs: erf.init_type cannot be given", {"inputs": {"erf.init_type": "uniform"}}),
                ("inputs: erf.input_sounding_file cannot be given", {"inputs": {"erf.input_sounding_file": "other"}}),
                ("inputs: geometry.is_periodic cannot be given", {"inputs": {"geometry.is_periodic": [1, 1, 0]}}),
                ("inputs: amr.restart cannot be given", {"inputs": {"amr.restart": "chk00010"}}),
                ("inputs: conductors.diagnostics_dir cannot be given", {"inputs": {"conductors.diagnostics_dir": "d"}}),
                ("speeds: a stripe's inputs: erf.restart cannot be given",
                 {"speeds": [{"speed": 10, "inputs": {"erf.restart": "chk00010"}}]}),
                ("inputs: yhi.type cannot be given: the driver sets the lateral faces", {"inputs": {"yhi.type": "Open"}}),
                ("inputs: the value of erf.x must be", {"inputs": {"erf.x": "a # b"}}),
                ("inputs: the value of erf.x must be", {"inputs": {"erf.x": "a = b"}}),
                ("seeds must be integers in [0, 2147483647], got 2147483648", {"seeds": [2147483648]}),
                ("speeds: the stripe 100 m/s at 10 m makes the profile's largest wind 175 m/s", {"speeds": [100]}),
                ("support_files: the base deck's erf.terrain_file_name names the file terrain.txt",
                 {"support_files": []}),
                ("speeds x directions x seeds gives 108000 members; at most 100000",
                 {"speeds": [0.01 * i for i in range(1, 1501)], "directions": [2.5 * i for i in range(-71, 1)],
                  "seeds": None}),
                ("inputs: conductors.gust_seed is set by the driver", {"inputs": {"conductors.gust_seed": 2}}),
                ("inputs: xlo.KE cannot be given: the driver sets the lateral faces", {"inputs": {"xlo.KE": 1.0}}),
                ("speeds: a stripe's inputs: xlo.KE cannot be given: the driver sets the lateral faces",
                 {"speeds": [{"speed": 10, "inputs": {"xlo.KE": 1.0}}]}),
                ("support_files: missing.txt is not a file", {"support_files": ["missing.txt"]}),
                ("support_files: '../x' must lie in the base deck's directory", {"support_files": ["../x"]}),
                ("support_files: a file is listed twice", {"support_files": ["terrain.txt", "terrain.txt"]}),
                ("base_deck: cannot read", {"base_deck": "missing.i"}),
                ("output_dir must name a directory", {"output_dir": ""}),
                ("--plan: unknown key(s) speed;", {"speed": [10]}),
                ("seeds need conductors.gust_type = random", {"inputs": {"conductors.gust_type": "factor"}})):
            with self.subTest(expect=expect):
                self.refuse("base_inflow", expect, **changes)

    def test_run_options(self):
        for expect, args in (
                ("--ranks must be in [1, --max-ranks = 4], got 5", ["--ranks", "5"]),
                ("--ranks must be in [1, --max-ranks = 4], got 0", ["--ranks", "0"]),
                ("--ranks must be in [1, --max-ranks = 2], got 3", ["--ranks", "3", "--max-ranks", "2"]),
                ("--max-ranks must be in [1, 1048576], got 0", ["--max-ranks", "0"]),
                ("--max-ranks must be in [1, 1048576], got 1048577", ["--max-ranks", "1048577"]),
                ("--exe: nowhere/erf_exec is not an executable file", ["--exe", "nowhere/erf_exec"]),
                ("--load-start must be in (0, 100000], got 0.0", ["--load-start", "0"]),
                ("--load-pause must be in (0, 100000], got -1.0", ["--load-pause", "-1"]),
                ("--load-resume must be in (0, 100000], got nan", ["--load-resume", "nan"]),
                ("--load-resume must be below --load-pause = 12, got 12", ["--load-resume", "12"]),
                ("--load-start must not be above --load-pause = 12, got 13", ["--load-start", "13"]),
                ("--load-poll must be in (0, 3600] s, got 0.0", ["--load-poll", "0"]),
                ("--load-poll must be in (0, 3600] s, got 4000.0", ["--load-poll", "4000"])):
            with self.subTest(expect=expect):
                self.refuse("base_inflow", expect, args=args)

    def test_the_deck_and_the_wind(self):
        # a direction along x on a deck whose only lateral faces are the y inflow and outflow has no side face to copy
        d = self.stage("base_inflow")
        deck = os.path.join(d, "span.i")
        text = read(deck)

        def with_deck(new_text, expect, **changes):
            with open(deck, "w") as f:
                f.write(new_text)
            code, _, err = call(["run", "--plan", self.edit_plan(d, **changes), "--dry-run"])
            self.assertEqual(code, 2, err)
            self.assertIn(expect, err)
            self.assertFalse(os.path.exists(os.path.join(d, "out")))

        with_deck(text.replace('"SlipWall"\nyhi', '"NoSlipWall"\nyhi'),
                  "base_deck: the side faces ylo and yhi carry different inputs")
        with_deck(text.replace("geometry.is_periodic = 0 0 0", "geometry.is_periodic = 0 1 0")
                      .replace('ylo.type           = "SlipWall"\nyhi.type           = "SlipWall"\n', ""),
                  "side_type: the direction 90 degrees blows along the xlo and xhi faces")
        with_deck(text.replace('xlo.type           = "Inflow"', 'xlo.type = "SlipWall"')
                      .replace('xhi.type           = "Outflow"', 'xhi.type = "SlipWall"')
                      .replace('xlo.dirichlet_file = "inflow_profile"   # z u v w theta\nxlo.KE             = 2.0\n', ""),
                  "directions: the direction 0 degrees blows through the xlo and xhi faces, which the base deck closes")
        with_deck(text.replace('erf.init_type           = "input_sounding"', 'erf.init_type = "uniform"'),
                  "base_deck: the driver sets the wind through erf.init_type = input_sounding or "
                  "conductors.prescribed_velocity")
        with_deck(text.replace('xlo.dirichlet_file = "inflow_profile"   # z u v w theta', "xlo.velocity = 8. 0. 0."),
                  "base_deck: xlo.velocity is not supported")
        with_deck(text + "erf.input_bndry_planes = 1\n", "base_deck: erf.input_bndry_planes = 1 is not supported")
        with_deck(text + "amr.restart = chk01000\n", "base_deck: amr.restart is not supported")
        with_deck(text + "erf.restart = chk01000\n", "base_deck: erf.restart is not supported")
        with_deck(text + "erf.abl_pressure_grad = -0.001 0. 0.\n", "base_deck: erf.abl_pressure_grad is not supported")
        with_deck(text + "erf.abl_geo_wind = 10. 0. 0.\n", "base_deck: erf.abl_geo_wind is not supported")
        with_deck(text.replace('"SlipWall"\nyhi', '"Open"\nyhi'), "base_deck: ylo.type = Open is not supported")
        with_deck(text.replace("conductors.diagnostics_dir = conductors", "conductors.diagnostics_dir = ../shared"),
                  "base_deck: conductors.diagnostics_dir = ../shared must lie inside the run directory")
        with_deck(text + "FILE = span.i\n", "includes itself through FILE")
        with_deck(text.replace("conductors.diagnostics_dir = conductors", "a = 1 b = 2"),
                  "the driver needs one 'key = value' per line")
        with_deck(text.replace("FILE = lines.inputs", "conductors.lines = S1"),
                  "does not set erf.structures.model = conductors")

        def with_wind(name, new_text, expect, **changes):
            old = read(os.path.join(d, name))
            with open(os.path.join(d, name), "w") as f:
                f.write(new_text)
            with_deck(text, expect, **changes)
            with open(os.path.join(d, name), "w") as f:
                f.write(old)

        with_wind("inflow_profile", "0.0 6.0 0.0 0.0\n20.0 10.0 0.0 0.0\n20.0 14.0 0.0 0.0\n", "must increase")
        with_wind("inflow_profile", "0.0 6.0 0.0\n20.0 10.0 0.0\n", "has fewer than 4 columns")
        with_wind("input_sounding", "1000.0 300.0 0.0\n0.0 300.0 0.0 6.0\n", "(5 columns expected)")
        with_wind("inflow_profile", "0.0 6.0 0.0 0.0\n20.0 nan 0.0 0.0\n", "is not finite")
        # a log-law file starts with no wind at the ground: the wind 5 m up is not half the wind at 10 m
        with_wind("inflow_profile", "0.0 0.0 0.0 0.0\n10.0 8.0 0.0 0.0\n100.0 14.0 0.0 0.0\n",
                  "reference_height = 5 m lies between the row at 0 m, where the wind is zero, and the row at 10 m",
                  reference_height=5)
        with_deck(text, "seeds need", reference_height=10, inputs={"conductors.gust_type": "factor"})
        self.edit_plan(d, inputs={"conductors.S1.rated_breaking_strength": 50000})
        with_deck(text + "conductors.prescribed_velocity = 0. 0. 1.\n", "base_deck: the base wind is zero",
                  reference_height=None, seeds=None)
        with_deck(text + "conductors.prescribed_velocity = 0. 5. 0.\n",
                  "reference_height is not used with conductors.prescribed_velocity", reference_height=10)

    def test_a_run_needs_an_executable(self):
        d = self.stage("base_prescribed")
        plan = os.path.join(d, "plan.json")
        for expect, args in (
                ("--exe is needed to run the members", []),
                ("--exe: nowhere/erf_exec is not an executable file", ["--exe", "nowhere/erf_exec"]),
                ("--exe: %s is not an executable file" % plan, ["--exe", plan]),
                ("--ranks must be 1 without a launcher", ["--exe", sys.executable, "--mpiexec", ""]),
                ("--mpiexec-args needs a launcher",
                 ["--exe", sys.executable, "--mpiexec", "", "--ranks", "1", "--mpiexec-args=-x"]),
                ("--mpiexec: 'no_such_launcher' is not on the PATH",
                 ["--exe", sys.executable, "--mpiexec", "no_such_launcher"])):
            with self.subTest(expect=expect):
                code, _, err = call(["run", "--plan", plan] + args)
                self.assertEqual(code, 2, err)
                self.assertIn(expect, err)
                self.assertFalse(os.path.exists(os.path.join(d, "out")))


# A stand-in for erf_exec: it reads the member's deck, notes that it ran, and writes the two
# tables as a run would. Seed 8 fails until the file "repaired" stands beside the plan; seed 5
# finishes without a demand table. With FAKE_WAIT set it starts a child process, writes both
# process ids to pids.txt, and the two wait until the file "go" stands beside the plan: a run of
# several processes that lasts as long as the test needs.
FAKE = """\
import os, sys, time
deck = dict((k.strip(), v.strip()) for k, _, v in (ln.partition("=") for ln in open(sys.argv[1]) if "=" in ln))
with open("runs.txt", "a") as f:
    f.write("run\\n")
seed = deck.get("conductors.gust_seed", "")
go = os.path.join("..", "..", "..", "go")
wait = "import os, time\\nfor i in range(1200):\\n    if os.path.exists(%r): break\\n    time.sleep(0.05)\\n" % go
if os.environ.get("FAKE_WAIT"):
    import subprocess
    # in its own session and process group, as Open MPI starts a rank
    child = subprocess.Popen([sys.executable, "-c", wait], start_new_session=True)
    with open("pids.txt", "w") as f:
        f.write("%d %d\\n" % (os.getpid(), child.pid))
    exec(wait)
    child.wait()
if seed == "8" and not os.path.exists(os.path.join("..", "..", "..", "repaired")):
    print("amrex::Abort::0::the stand-in fails seed 8 !!!")
    sys.exit(3)
os.makedirs("conductors", exist_ok=True)
if seed != "5":
    with open("conductors/demand_table.csv", "w") as f:
        f.write("kind,component,peak_ratio\\nconductor,S1,%s\\n" % deck["conductors.prescribed_velocity"].split()[0])
with open("conductors/intensity.csv", "w") as f:
    f.write("point,speed\\nS1,%s\\n" % deck["conductors.prescribed_velocity"].split()[0])
"""


class Run(Case):
    def stage_run(self, seeds):
        d = self.stage("base_prescribed")
        with open(os.path.join(d, "span.i"), "a") as f:
            f.write("conductors.gust_type = random\n")
        plan = self.edit_plan(d, seeds=seeds, directions=[0], speeds=[20, 40])
        with open(os.path.join(d, "fake_erf.py"), "w") as f:
            f.write(FAKE)
        exe = os.path.join(d, "fake_erf")
        with open(exe, "w") as f:
            # not exec'ed: the shell is the launcher, the stand-in its child, as a rank is mpiexec's
            f.write("#!/bin/sh\n'%s' '%s' \"$@\"\nexit $?\n" % (sys.executable, os.path.join(d, "fake_erf.py")))
        os.chmod(exe, os.stat(exe).st_mode | stat.S_IXUSR)
        return d, ["run", "--plan", plan, "--exe", exe, "--mpiexec", "", "--ranks", "1"]

    def runs(self, d, member):
        path = os.path.join(d, "out", "members", member, "runs.txt")
        return len(read(path).splitlines()) if os.path.exists(path) else 0

    def test_finished_members_are_skipped_and_failed_ones_run_again(self):
        d, argv = self.stage_run([3, 8, 5])
        names = ["m001_s20_d0_seed3", "m002_s20_d0_seed8", "m003_s20_d0_seed5",
                 "m004_s40_d0_seed3", "m005_s40_d0_seed8", "m006_s40_d0_seed5"]
        code, out, err = call(argv + ["--no-load-gate", "--load-poll", "0.05"])
        self.assertEqual((code, err), (1, ""))     # two members failed, two lack a table
        self.assertEqual([self.runs(d, n) for n in names], [1] * 6)
        status = [json.loads(read(os.path.join(d, "out", "members", n, "member_status.json"))) for n in names]
        self.assertEqual([s["status"] for s in status], ["done", "failed", "done"] * 2)
        self.assertEqual([s["returncode"] for s in status], [0, 3, 0] * 2)
        self.assertIn("NOT COLLECTED m002_s20_d0_seed8: failed (amrex::Abort::0::the stand-in fails seed 8 !!!)", out)
        self.assertIn("NOT COLLECTED m003_s20_d0_seed5: no_demand_table (conductors/demand_table.csv: missing)", out)
        self.assertIn("collected 2 of 6 members", out)
        self.assertEqual(read(os.path.join(d, "out", "demand_table_all.csv")),
                         "member,stripe,nominal_direction,seed,kind,component,peak_ratio\n"
                         "m001_s20_d0_seed3,20,0,3,conductor,S1,20\n"
                         "m004_s40_d0_seed3,40,0,3,conductor,S1,40\n")

        # again, with seed 8 repaired: only the failed members run; the finished ones are not touched
        open(os.path.join(d, "repaired"), "w").close()
        code, out, _ = call(argv + ["--no-load-gate", "--load-poll", "0.05"])
        self.assertEqual(code, 1)                  # the two members without a demand table remain
        self.assertEqual([self.runs(d, n) for n in names], [1, 2, 1, 1, 2, 1])
        self.assertIn("m001_s20_d0_seed3: finished earlier, skipped", out)
        self.assertIn("collected 4 of 6 members", out)
        self.assertEqual(read(os.path.join(d, "out", "intensity_all.csv")),
                         "member,stripe,nominal_direction,seed,point,speed\n"
                         "m001_s20_d0_seed3,20,0,3,S1,20\n"
                         "m002_s20_d0_seed8,20,0,8,S1,20\n"
                         "m004_s40_d0_seed3,40,0,3,S1,40\n"
                         "m005_s40_d0_seed8,40,0,8,S1,40\n")

    def test_a_rerun_does_not_take_an_earlier_attempts_tables(self):
        # seed 5 finishes without a demand table; one left in its directory by an earlier attempt must not be
        # gathered as this run's
        d, argv = self.stage_run([5])
        self.assertEqual(call(argv + ["--dry-run"])[0], 0)
        for name in ("m001_s20_d0_seed5", "m002_s40_d0_seed5"):
            stale = os.path.join(d, "out", "members", name, "conductors")
            os.makedirs(stale)
            with open(os.path.join(stale, "demand_table.csv"), "w") as f:
                f.write("kind,component,peak_ratio\nconductor,S1,99\n")
        code, out, _ = call(argv + ["--no-load-gate", "--load-poll", "0.05"])
        self.assertEqual(code, 1)
        self.assertIn("collected 0 of 2 members", out)
        self.assertEqual(read(os.path.join(d, "out", "demand_table_all.csv")),
                         "member,stripe,nominal_direction,seed\n")

    def test_a_second_driver_on_the_same_directory_is_refused(self):
        d, argv = self.stage_run([3])
        self.assertEqual(call(argv + ["--dry-run"])[0], 0)
        lock = os.path.join(d, "out", "run.lock")
        with open(lock, "w") as f:
            f.write("1\n")
        code, _, err = call(argv + ["--no-load-gate"])
        self.assertEqual(code, 2)
        self.assertIn("another driver is running these members", err)
        self.assertEqual(self.runs(d, "m001_s20_d0_seed3"), 0)
        os.remove(lock)
        self.assertEqual(call(argv + ["--no-load-gate", "--load-poll", "0.05"])[0], 0)
        self.assertFalse(os.path.exists(lock))

    def states(self, pids):
        out = []
        for pid in pids:
            ps = subprocess.run(["ps", "-o", "stat=", "-p", str(pid)], stdout=subprocess.PIPE, universal_newlines=True)
            out.append(ps.stdout.strip()[:1])      # "" for a process that is gone; "Z" once it ended and waits
        return out

    def wait_for(self, what, seconds=20.0):
        t0 = time.time()
        while not what():
            self.assertLess(time.time() - t0, seconds, "the stand-in run did not get there")
            time.sleep(0.02)

    def test_the_gate_waits_for_a_quiet_machine_and_pauses_a_run_while_it_is_busy(self):
        if fe.descendants(os.getpid()) is None:
            self.skipTest("ps is not available")
        d, argv = self.stage_run([3])
        os.environ["FAKE_WAIT"] = "1"
        self.addCleanup(os.environ.pop, "FAKE_WAIT", None)
        pids_file = os.path.join(d, "out", "members", "m001_s20_d0_seed3", "pids.txt")
        # the load the gate sees: busy twice before the first run starts, then quiet; above the pause level
        # at the first look during the run, still above the resume level at the second, below it at the third
        loads = [9.0, 6.0, 5.9, 12.1, 7.0, 6.9]
        seen, states, log = [], [], []

        def load():
            value = loads[len(seen)] if len(seen) < len(loads) else 1.0
            seen.append(value)
            if len(seen) == 4:      # the run's two processes must be there before the pause
                self.wait_for(lambda: os.path.exists(pids_file) and len(read(pids_file).split()) == 2)
            if len(seen) == 5:      # the look after the pause: every process of the run is stopped
                states.extend(self.states(read(pids_file).split()))
            if len(seen) == 6:      # the run may end once it goes on
                open(os.path.join(d, "go"), "w").close()
            return value

        gate = fe.Gate(6.0, 12.0, 7.0, 0.1, load=load, log=log.append)
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = fe.main(argv, gate=gate)
        self.assertEqual(code, 0, out.getvalue())
        self.assertEqual([ln for ln in log if ln.startswith("m001_s20_d0_seed3")], [
            "m001_s20_d0_seed3: waiting for the 1-minute load to fall below 6",
            "m001_s20_d0_seed3: paused, load 12.1 above 12",
            "m001_s20_d0_seed3: resumed, load 6.9 below 7"])
        self.assertEqual(seen[:6], loads)
        self.assertEqual(states, ["T", "T"])
        self.assertEqual(self.runs(d, "m001_s20_d0_seed3"), 1)
        self.assertEqual(self.runs(d, "m002_s40_d0_seed3"), 1)

    def test_a_load_exactly_at_the_pause_level_does_not_pause(self):
        d, argv = self.stage_run([3])
        os.environ["FAKE_WAIT"] = "1"
        self.addCleanup(os.environ.pop, "FAKE_WAIT", None)
        seen, log = [], []

        def load():
            seen.append(1)
            if len(seen) == 3:
                open(os.path.join(d, "go"), "w").close()
            return 12.0 if len(seen) in (2, 3) else 1.0      # at the level during the first run only

        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(fe.main(argv, gate=fe.Gate(6.0, 12.0, 7.0, 0.1, load=load, log=log.append)), 0)
        self.assertEqual([ln for ln in log if "paused" in ln], [])

    def start_alone(self, d, argv):
        """The driver in its own process, running the first member's stand-in, which waits."""
        driver = os.path.join(HERE, "..", "fragility_ensemble.py")
        env = dict(os.environ, FAKE_WAIT="1")
        proc = subprocess.Popen([sys.executable, "-I", driver] + argv + ["--no-load-gate", "--load-poll", "0.1"],
                                env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        pids_file = os.path.join(d, "out", "members", "m001_s20_d0_seed3", "pids.txt")
        try:
            self.wait_for(lambda: os.path.exists(pids_file) and len(read(pids_file).split()) == 2)
        except BaseException:
            proc.kill()
            raise
        return proc, [int(p) for p in read(pids_file).split()]

    def test_a_driver_told_to_end_takes_its_run_with_it(self):
        # SIGTERM to the driver: the run's processes end too, the lock goes, and the member has no status, so
        # it runs again next time
        if fe.descendants(os.getpid()) is None:
            self.skipTest("ps is not available")
        d, argv = self.stage_run([3])
        proc, pids = self.start_alone(d, argv)
        proc.terminate()
        self.assertEqual(proc.wait(timeout=30), 128 + 15)
        self.wait_for(lambda: all(s in ("", "Z") for s in self.states(pids)))
        member = os.path.join(d, "out", "members", "m001_s20_d0_seed3")
        self.assertFalse(os.path.exists(os.path.join(member, "member_status.json")))
        self.assertFalse(os.path.exists(os.path.join(d, "out", "run.lock")))
        self.assertEqual(self.runs(d, "m002_s40_d0_seed3"), 0)
        open(os.path.join(d, "go"), "w").close()
        self.assertEqual(call(argv + ["--no-load-gate", "--load-poll", "0.05"])[0], 0)
        self.assertEqual(self.runs(d, "m001_s20_d0_seed3"), 2)

    def test_ending_a_paused_run_leaves_no_stopped_process(self):
        if fe.descendants(os.getpid()) is None:
            self.skipTest("ps is not available")
        d, argv = self.stage_run([3])
        self.assertEqual(call(argv + ["--dry-run"])[0], 0)
        member = os.path.join(d, "out", "members", "m001_s20_d0_seed3")
        proc = subprocess.Popen([argv[argv.index("--exe") + 1], "span.i"], cwd=member, start_new_session=True,
                                env=dict(os.environ, FAKE_WAIT="1"))
        try:
            pids_file = os.path.join(member, "pids.txt")
            self.wait_for(lambda: os.path.exists(pids_file) and len(read(pids_file).split()) == 2)
            pids = [int(p) for p in read(pids_file).split()]
            self.assertEqual(sorted(fe.descendants(proc.pid)), sorted(pids))
            fe.signal_all(pids + [proc.pid], fe.signal.SIGSTOP)      # the launcher too
            self.wait_for(lambda: self.states(pids + [proc.pid]) == ["T", "T", "T"])
        finally:
            t0 = time.time()
            fe.end_run(proc)
        # the launcher was continued and ended on SIGTERM, without the 10 s wait for a process that cannot answer
        self.assertLess(time.time() - t0, 8.0)
        self.wait_for(lambda: all(s in ("", "Z") for s in self.states(pids)))
        self.assertIsNotNone(proc.returncode)


# A stand-in for a sounding deck. A deck without erf.structures.model is a spin-up: it writes the
# checkpoint chk<max_step> (a Header file). A member must restart from a checkpoint that is there; it
# writes its demand table (the 10 m/s stripe's hardware exactly at its capacity, which is not an
# exceedance: the model's rule is demand > capacity) and, for the 20 m/s stripe, a release log with two releases of one line
# and a tower-failure log with one tower.
FAKE_SPIN = """\
import os, sys
deck = dict((k.strip(), v.strip()) for k, _, v in (ln.partition("=") for ln in open(sys.argv[1]) if "=" in ln))
with open("runs.txt", "a") as f:
    f.write("run\\n")
fast = deck["xlo.KE"] == "12.5"
if "erf.structures.model" not in deck:
    if fast and not os.path.exists(os.path.join("..", "..", "..", "repaired")):
        print("amrex::Abort::0::the stand-in fails the 20 m/s spin-up !!!")
        sys.exit(5)
    if any(k.startswith("conductors.") for k in deck) or "amr.restart" in deck:
        sys.exit(6)
    chk = "chk%05d" % int(deck["max_step"])
    os.makedirs(chk, exist_ok=True)
    open(os.path.join(chk, "Header"), "w").close()
    sys.exit(0)
if not os.path.isfile(os.path.join(deck["amr.restart"], "Header")):
    print("amrex::Abort::0::no checkpoint !!!")
    sys.exit(4)
os.makedirs("conductors", exist_ok=True)
with open("conductors/demand_table.csv", "w") as f:
    f.write("kind,component,peak_ratio\\nconductor,S1,%s\\nhardware,S1_end_a,%s\\n" % (
        "1.4" if fast else "0.6", "1.1" if fast else "1"))
with open("conductors/intensity.csv", "w") as f:
    f.write("point,speed\\nS1,%s\\n" % ("20" if fast else "10"))
with open("conductors/failures.dat", "w") as f:
    f.write("time kind component demand capacity ratio\\n")
    if fast:
        f.write("2.5 conductor S1 70000 50000 1.4\\n")
if fast:
    with open("conductors/breaks.dat", "w") as f:
        f.write("time line released component tension\\n2.5 S1 span0_end1 S1 70000\\n2.5 S1 span1_end0 T1 0\\n")
    with open("conductors/tower_failures.dat", "w") as f:
        f.write("time tower component releases\\n2.5 T1 T1_leg 2\\n")
"""


class Spinup(Case):
    def stage_spinup(self, **plan):
        d = self.stage("base_inflow")
        path = self.edit_plan(d, spinup={"max_step": 4}, **plan)
        return d, path

    def test_members_differing_only_in_the_seed_share_a_spin_up_of_the_flow_alone(self):
        d, plan = self.stage_spinup()
        code, out, err = call(["run", "--plan", plan, "--dry-run"])
        self.assertEqual((code, err), (0, ""), out)
        manifest = json.loads(read(os.path.join(d, "out", "manifest.json")))
        names = ["u001_s10_d0", "u002_s10_d90", "u003_s10_dm135", "u004_s20_d0", "u005_s20_d90", "u006_s20_dm135"]
        self.assertEqual([g["spinup"] for g in manifest["spinups"]], names)
        for i, g in enumerate(manifest["spinups"]):
            self.assertEqual(g["directory"], "spinup/" + names[i])
            self.assertEqual((g["max_step"], g["checkpoint"]), (4, "chk00004"))
            # members m001..m012 run seeds 3 and 8 of each stripe and direction in turn
            self.assertEqual(g["members"], [manifest["members"][2 * i]["member"], manifest["members"][2 * i + 1]["member"]])
        for i, m in enumerate(manifest["members"]):
            g = names[i // 2]
            self.assertEqual(m["spinup"], g)
            self.assertEqual(m["inputs"]["amr.restart"], "../../spinup/%s/chk00004" % g)
            member = read(os.path.join(d, "out", m["directory"], "span.i"))
            self.assertIn("amr.restart = ../../spinup/%s/chk00004\n" % g, member)
            spin = read(os.path.join(d, "out", "spinup", g, "span.i"))
            # the spin-up is the member's deck without the conductor keys and the restart, run to step 4
            # with a checkpoint there
            only_member, only_spin = self.changed_lines(member, spin)
            self.assertEqual(sorted(only_spin), ["erf.check_int = 4", "max_step = 4"])
            self.assertEqual(sorted(only_member), sorted(
                ["amr.restart = ../../spinup/%s/chk00004" % g, "conductors.diagnostics_dir = conductors",
                 "erf.structures.model = conductors", "conductors.lines     = S1",
                 "conductors.gust_type = random", "conductors.gust_seed = %d" % m["seed"],
                 "conductors.S1.rated_breaking_strength = 50000", "conductors.intensity_measure = wind",
                 "max_step = %d" % (20 if m["stripe"] == 20 else 10)]))
            for name in ("inflow_profile", "input_sounding", "terrain.txt"):
                self.assertEqual(read(os.path.join(d, "out", "spinup", g, name)),
                                 read(os.path.join(d, "out", m["directory"], name)))

    def test_without_spinup_the_manifest_has_no_spin_ups(self):
        d = self.stage("base_inflow")
        self.assertEqual(call(["run", "--plan", os.path.join(d, "plan.json"), "--dry-run"])[0], 0)
        manifest = json.loads(read(os.path.join(d, "out", "manifest.json")))
        self.assertNotIn("spinups", manifest)
        self.assertFalse(any("spinup" in m or "amr.restart" in m["inputs"] for m in manifest["members"]))
        self.assertFalse(os.path.exists(os.path.join(d, "out", "spinup")))

    def test_spin_up_options_out_of_range_are_refused(self):
        for value, expect in (([4], "spinup must be {\"max_step\": n}"),
                              ({"max_step": 4, "steps": 4}, "spinup must be {\"max_step\": n}"),
                              ({"max_step": 0}, "spinup: max_step must be an integer in [1, 100000000], got 0"),
                              ({"max_step": 2.5}, "spinup: max_step must be an integer in [1, 100000000], got 2.5"),
                              ({"max_step": True}, "spinup: max_step must be an integer"),
                              ({"max_step": 10}, "spinup: max_step = 10 must be below m001_s10_d0_seed3's max_step = 10")):
            d = self.stage("base_inflow")
            code, _, err = call(["run", "--plan", self.edit_plan(d, spinup=value), "--dry-run"])
            self.assertEqual(code, 2, value)
            self.assertIn(expect, err)
        d = self.stage("base_inflow")
        code, _, err = call(["run", "--plan", self.edit_plan(d, spinup={"max_step": 4},
                                                             inputs={"erf.check_int": 2}), "--dry-run"])
        self.assertEqual(code, 2)
        self.assertIn("spinup: m001_s10_d0_seed3's inputs set erf.check_int; the driver sets the spin-up's checkpoint "
                      "itself", err)
        d = self.stage("base_prescribed")
        code, _, err = call(["run", "--plan", self.edit_plan(d, spinup={"max_step": 4}), "--dry-run"])
        self.assertEqual(code, 2)
        self.assertIn("spinup needs a flow to spin up: conductors.prescribed_velocity runs no flow", err)

    def stage_run(self):
        d, plan = self.stage_spinup(directions=[0])
        with open(os.path.join(d, "fake_erf.py"), "w") as f:
            f.write(FAKE_SPIN)
        exe = os.path.join(d, "fake_erf")
        with open(exe, "w") as f:
            f.write("#!/bin/sh\n'%s' '%s' \"$@\"\nexit $?\n" % (sys.executable, os.path.join(d, "fake_erf.py")))
        os.chmod(exe, os.stat(exe).st_mode | stat.S_IXUSR)
        return d, ["run", "--plan", plan, "--exe", exe, "--mpiexec", "", "--ranks", "1", "--no-load-gate",
                   "--load-poll", "0.05"]

    def runs(self, d, rel):
        path = os.path.join(d, "out", rel, "runs.txt")
        return len(read(path).splitlines()) if os.path.exists(path) else 0

    def test_each_spin_up_runs_once_and_a_failed_one_holds_back_its_members(self):
        d, argv = self.stage_run()
        members = ["m001_s10_d0_seed3", "m002_s10_d0_seed8", "m003_s20_d0_seed3", "m004_s20_d0_seed8"]
        code, out, err = call(argv)
        self.assertEqual((code, err), (1, ""), out)
        self.assertEqual([self.runs(d, "spinup/" + g) for g in ("u001_s10_d0", "u002_s20_d0")], [1, 1])
        self.assertEqual([self.runs(d, "members/" + m) for m in members], [1, 1, 0, 0])
        self.assertIn("u002_s20_d0: spin-up failed (return code 5)", out)
        self.assertIn("m004_s20_d0_seed8: not run, its spin-up u002_s20_d0 failed", out)
        self.assertIn("NOT COLLECTED m003_s20_d0_seed3: failed (not run: its spin-up spinup/u002_s20_d0 failed "
                      "(return code 5))", out)

        # repaired: the finished spin-up and members are not run again; the failed spin-up runs once more,
        # then its two members. m001 is made to run again with a release log an earlier attempt left: the
        # rerun writes none (failure_break off), so the old one must not be gathered as its
        open(os.path.join(d, "repaired"), "w").close()
        m001 = os.path.join(d, "out", "members", "m001_s10_d0_seed3")
        os.remove(os.path.join(m001, "member_status.json"))
        with open(os.path.join(m001, "conductors", "breaks.dat"), "w") as f:
            f.write("time line released component tension\n1 S1 span0_end0 S1 1\n")
        code, out, _ = call(argv)
        self.assertEqual(code, 0, out)
        self.assertEqual([self.runs(d, "spinup/" + g) for g in ("u001_s10_d0", "u002_s20_d0")], [1, 2])
        self.assertEqual([self.runs(d, "members/" + m) for m in members], [2, 1, 1, 1])
        self.assertIn("collected 4 of 4 members", out)

        # the failure logs, gathered with the key columns, and each member's damage summary
        self.assertEqual(read(os.path.join(d, "out", "breaks_all.csv")),
                         "member,stripe,nominal_direction,seed,time,line,released,component,tension\n"
                         "m003_s20_d0_seed3,20,0,3,2.5,S1,span0_end1,S1,70000\n"
                         "m003_s20_d0_seed3,20,0,3,2.5,S1,span1_end0,T1,0\n"
                         "m004_s20_d0_seed8,20,0,8,2.5,S1,span0_end1,S1,70000\n"
                         "m004_s20_d0_seed8,20,0,8,2.5,S1,span1_end0,T1,0\n")
        self.assertEqual(read(os.path.join(d, "out", "tower_failures_all.csv")),
                         "member,stripe,nominal_direction,seed,time,tower,component,releases\n"
                         "m003_s20_d0_seed3,20,0,3,2.5,T1,T1_leg,2\n"
                         "m004_s20_d0_seed8,20,0,8,2.5,T1,T1_leg,2\n")
        self.assertEqual(read(os.path.join(d, "out", "damage_table_all.csv")),
                         "member,stripe,nominal_direction,seed,max_peak_ratio,components_exceeded,failures_log,failures,"
                         "breaks_log,lines_released,releases,tower_log,towers_lost\n"
                         "m001_s10_d0_seed3,10,0,3,1,0,yes,0,no,,,no,\n"
                         "m002_s10_d0_seed8,10,0,8,1,0,yes,0,no,,,no,\n"
                         "m003_s20_d0_seed3,20,0,3,1.4,2,yes,1,yes,1,2,yes,1\n"
                         "m004_s20_d0_seed8,20,0,8,1.4,2,yes,1,yes,1,2,yes,1\n")

        # a log with another header is named, and its member is not collected
        path = os.path.join(d, "out", "members", "m003_s20_d0_seed3", "conductors", "breaks.dat")
        with open(path, "w") as f:
            f.write("time line component tension\n")
        code, out, _ = call(["collect", "--plan", argv[2]])
        self.assertEqual(code, 1)
        self.assertIn("NOT COLLECTED m003_s20_d0_seed3: bad_breaks_log (conductors/breaks.dat: its header is not "
                      "'time line released component tension')", out)

    def test_a_finished_spin_up_whose_checkpoint_is_gone_runs_again(self):
        d, argv = self.stage_run()
        open(os.path.join(d, "repaired"), "w").close()
        self.assertEqual(call(argv)[0], 0)
        spin = os.path.join(d, "out", "spinup", "u001_s10_d0")
        shutil.rmtree(os.path.join(spin, "chk00004"))
        os.remove(os.path.join(d, "out", "members", "m001_s10_d0_seed3", "member_status.json"))
        code, out, _ = call(argv)
        self.assertEqual(code, 0, out)
        self.assertEqual(self.runs(d, "spinup/u001_s10_d0"), 2)
        self.assertEqual(self.runs(d, "spinup/u002_s20_d0"), 1)

    def test_a_spin_up_that_writes_no_checkpoint_has_failed(self):
        d, argv = self.stage_run()
        open(os.path.join(d, "repaired"), "w").close()
        # the stand-in writes chk<max_step>; a deck whose check file is another name leaves no chk00004
        text = read(os.path.join(d, "fake_erf.py")).replace('chk = "chk%05d"', 'chk = "other%05d"')
        with open(os.path.join(d, "fake_erf.py"), "w") as f:
            f.write(text)
        code, out, _ = call(argv)
        self.assertEqual(code, 1)
        self.assertIn("u001_s10_d0: spin-up failed (return code 0)", out)
        status = json.loads(read(os.path.join(d, "out", "spinup", "u001_s10_d0", "member_status.json")))
        self.assertEqual(status["detail"], "the run ended without writing chk00004 (does its stop_time end it before "
                                           "max_step?)")
        self.assertEqual(self.runs(d, "members/m001_s10_d0_seed3"), 0)


DEMAND_HEADER = ("kind,component,peak_ratio,peak_demand,peak_capacity,peak_time,point,peak_speed,peak_direction,"
                 "peak_k,peak_wind_time,max_speed,max_direction,max_k,max_time,exceeded,first_time,first_demand,"
                 "first_capacity,first_speed,first_direction,first_k,first_wind_time")


class Collect(Case):
    def stage(self, name):
        # the failed member's log (a *.log file is not kept in the repository, so it is written here)
        d = super().stage(name)
        with open(os.path.join(d, "members", "m003_s30_d0_seed1", "run.log"), "w") as f:
            f.write("Step 3 starts\namrex::Abort::0::conductors: line S1 left the domain !!!\nSIGABRT\n")
        return d

    def test_complete_members_are_gathered_and_the_others_reported(self):
        d = self.stage("collect")
        code, out, err = call(["collect", "--plan", os.path.join(d, "plan.json")])
        self.assertEqual((code, err), (1, ""))       # four of the six members are not complete
        self.assertEqual(read(os.path.join(d, "demand_table_all.csv")),
                         "member,stripe,nominal_direction,seed," + DEMAND_HEADER + "\n"
                         "m001_s20_d0_seed1,20,0,1,conductor,S1,0.41,20500,50000,7.5,S1,20.4,0.5,,7.5,20.6,0.4,,8,0,,,,,,,\n"
                         "m001_s20_d0_seed1,20,0,1,hardware,S1_end_a,0.52,20800,40000,7.5,S1_end_a,20.1,0.2,,7.5,20.3,0.1,,8,0,,,,,,,\n"
                         "m005_s40_d0_seed1,40,0,1,conductor,S1,1.64,82000,50000,6,S1,40.8,0.5,,6,41.2,0.4,,8,1,2.5,50100,50000,31.5,0.6,,2.5\n"
                         "m005_s40_d0_seed1,40,0,1,hardware,S1_end_a,2.08,83200,40000,6,S1_end_a,40.2,0.2,,6,40.6,0.1,,8,1,2,40050,40000,28.4,0.3,,2\n")
        self.assertEqual(read(os.path.join(d, "intensity_all.csv")),
                         "member,stripe,nominal_direction,seed,point,x,y,z,height,speed,direction,k,time,max_speed,"
                         "max_direction,max_k,max_time\n"
                         "m001_s20_d0_seed1,20,0,1,S1,450,500,10,10,20.2,0.5,,10,20.6,0.4,,8\n"
                         "m005_s40_d0_seed1,40,0,1,S1,450,500,10,10,40.4,0.5,,10,41.2,0.4,,8\n")
        with open(os.path.join(d, "members_status.csv"), newline="") as f:
            status = list(csv.reader(f))
        self.assertEqual(status, [
            ["member", "stripe", "nominal_direction", "seed", "status", "returncode", "demand_rows", "intensity_rows",
             "detail"],
            ["m001_s20_d0_seed1", "20", "0", "1", "complete", "0", "2", "1", ""],
            ["m002_s20_d0_seed2", "20", "0", "2", "no_demand_table", "0", "", "",
             "conductors/demand_table.csv: missing"],
            ["m003_s30_d0_seed1", "30", "0", "1", "failed", "134", "", "",
             "amrex::Abort::0::conductors: line S1 left the domain !!!"],
            ["m004_s30_d0_seed2", "30", "0", "2", "not_run", "", "", "", ""],
            ["m005_s40_d0_seed1", "40", "0", "1", "complete", "0", "2", "1", ""],
            ["m006_s40_d0_seed2", "40", "0", "2", "bad_demand_table", "0", "", "",
             "conductors/demand_table.csv: its header differs from m001_s20_d0_seed1's"]])
        for line in ("collected 2 of 6 members",
                     "NOT COLLECTED m002_s20_d0_seed2: no_demand_table",
                     "NOT COLLECTED m003_s30_d0_seed1: failed",
                     "NOT COLLECTED m004_s30_d0_seed2: not_run",
                     "NOT COLLECTED m006_s40_d0_seed2: bad_demand_table"):
            self.assertIn(line, out)

    def test_a_whole_ensemble_returns_zero_and_bad_tables_are_named(self):
        d = self.stage("collect")
        members = os.path.join(d, "members")
        table = os.path.join("conductors", "demand_table.csv")
        good = read(os.path.join(members, "m001_s20_d0_seed1", table))
        done = read(os.path.join(members, "m001_s20_d0_seed1", "member_status.json"))
        point = read(os.path.join(members, "m001_s20_d0_seed1", "conductors", "intensity.csv"))
        for name in ("m002_s20_d0_seed2", "m003_s30_d0_seed1", "m004_s30_d0_seed2", "m006_s40_d0_seed2"):
            os.makedirs(os.path.join(members, name, "conductors"), exist_ok=True)
            for rel, text in ((table, good), ("member_status.json", done), ("conductors/intensity.csv", point)):
                with open(os.path.join(members, name, rel), "w") as f:
                    f.write(text)
        code, out, _ = call(["collect", "--plan", os.path.join(d, "plan.json")])
        self.assertEqual(code, 0, out)
        self.assertIn("collected 6 of 6 members", out)
        self.assertEqual(len(read(os.path.join(d, "demand_table_all.csv")).splitlines()), 1 + 6 * 2)

        def status_with(name, rel, text):
            with open(os.path.join(members, name, rel), "w") as f:
                f.write(text)
            self.assertEqual(call(["collect", "--plan", os.path.join(d, "plan.json")])[0], 1)
            with open(os.path.join(d, "members_status.csv"), newline="") as f:
                return {r[0]: (r[4], r[8]) for r in csv.reader(f)}[name]

        self.assertEqual(status_with("m002_s20_d0_seed2", table, good.splitlines()[0] + "\n"),
                         ("bad_demand_table", "conductors/demand_table.csv: no rows"))
        self.assertEqual(status_with("m003_s30_d0_seed1", table, good + "conductor,S2,0.4\n"),
                         ("bad_demand_table", "conductors/demand_table.csv: a row's length differs from the header's"))
        self.assertEqual(status_with("m004_s30_d0_seed2", table, good.replace("kind,component", "kind,part")),
                         ("bad_demand_table", "conductors/demand_table.csv: no column component"))
        self.assertEqual(status_with("m006_s40_d0_seed2", "conductors/intensity.csv", ""),
                         ("bad_intensity_table", "conductors/intensity.csv: empty"))
        self.assertEqual(status_with("m005_s40_d0_seed1", "conductors/intensity.csv",
                                     point.replace("point,x", "point,stripe")),
                         ("bad_intensity_table", "conductors/intensity.csv: has the key column stripe"))

    def test_a_directory_without_a_manifest_is_refused(self):
        d = self.stage("collect")
        os.remove(os.path.join(d, "manifest.json"))
        code, _, err = call(["collect", "--plan", os.path.join(d, "plan.json")])
        self.assertEqual(code, 2)
        self.assertIn("output_dir:", err)
        self.assertIn("holds no manifest.json written by this driver", err)


if __name__ == "__main__":
    unittest.main(verbosity=1)
