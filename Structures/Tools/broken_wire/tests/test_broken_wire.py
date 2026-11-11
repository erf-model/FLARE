#!/usr/bin/env python3
"""Tests of broken_wire.py (standard library only; no ERF run).

    python3 test_broken_wire.py            every test
    python3 test_broken_wire.py Cable      one class

Cable: the cable's reach against a numerical integration of its shape that shares no formula with
the tool, the lumped chain against the continuous cable, and the level span against the elastic
catenary of the model's documentation (Docs/sphinx_doc/theory/Conductors.rst, "Verification").
Balance: solved sections against checks that do not use the tool's balance at the strings' feet:
the supports together carry the weight of what hangs from them and no net pull along the line; the
potential energy, summed by numerical integration, does not change to first order when a string's
foot is moved; a string points along the force it carries. Limits: states known without the tool.
Robust: a section whose strings lean far, and a state the tool refuses. Reference: the numbers the
model's unit test BrokenWire.* is pinned to. Refused: every plan key out of range.
"""
import contextlib
import importlib.util
import io
import json
import math
import os
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("broken_wire", os.path.join(HERE, "..", "broken_wire.py"))
bw = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bw)

# 795 kcmil ACSR (Drake) as in the model's tests: diameter, mass, axial stiffness; air of 1.2 kg/m^3
DRAKE = {"diameter": 0.0281, "mass_per_length": 1.628, "axial_stiffness": 3.0e7, "air_density": 1.2}
W = (1.628 - 1.2 * 0.25 * math.pi * 0.0281 ** 2) * 9.81     # N/m


def section(**more):
    """Spans of 300, 500 and 100 m on 2.5 m strings of 60 kg, every attachment 60 m up, the spans'
    slack in proportion to the cube of their length (the section of the model's test
    Conductors.AFailedTowerPassesItsLoadToTheNext)."""
    plan = dict(DRAKE, end_a=[100.0, 60.0], end_b=[1000.0, 60.0], towers=[[400.0, 60.0], [900.0, 60.0]],
                lengths=[301.5, 506.944, 100.056], insulator_length=2.5, insulator_mass=60.0)
    plan.update(more)
    return plan


def simpson(f, a, b, n=20000):
    """Composite Simpson's rule on n (even) intervals."""
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4.0 if i % 2 else 2.0) * f(a + i * h)
    return s * h / 3.0


def reach_by_integration(H, V0, L, w, EA):
    """End b from end a by integrating the stretched tangent along the unstretched length."""
    def T(s):
        return math.hypot(H, V0 + w * s)
    dx = simpson(lambda s: H * (1.0 / T(s) + 1.0 / EA), 0.0, L)
    dz = simpson(lambda s: (V0 + w * s) * (1.0 / T(s) + 1.0 / EA), 0.0, L)
    return dx, dz


def cable_energy(a, b, L, w, EA):
    """Potential energy (J) of a continuous cable between the points a and b: the strain energy
    T^2 / (2 EA) and the weight times the height, integrated along the unstretched length."""
    H, V0 = bw.solve_cable((b[0] - a[0], b[1] - a[1]), L, w, EA)

    def T(s):
        return math.hypot(H, V0 + w * s)

    def height(s):
        # z(s) - z_a, by the same closed integral the tool uses for dz; the check of that
        # integral is Cable.test_reach_against_integration
        return (V0 * s + 0.5 * w * s * s) / EA + (T(s) - T(0.0)) / w

    strain = simpson(lambda s: T(s) ** 2 / (2.0 * EA), 0.0, L, 4000)
    weight = simpson(lambda s: w * (a[1] + height(s)), 0.0, L, 4000)
    return strain + weight


class Cable(unittest.TestCase):
    def test_reach_against_integration(self):
        cases = [(13251.0, -2406.0, 301.5, W, 3.0e7),       # a level span
                 (9000.0, -6000.0, 420.0, W, 3.0e7),        # end b well below the low point's side
                 (5000.0, 300.0, 200.0, W, 1.0e6),          # rising from end a, soft
                 (-7000.0, -900.0, 150.0, 2.0 * W, 3.0e7),  # end b to the left
                 (40.0, -580.0, 2.5, 234.8, 1.0e7)]         # a string
        for H, V0, L, w, EA in cases:
            dx, dz = bw.cable_reach(H, V0, L, w, EA)
            ix, iz = reach_by_integration(H, V0, L, w, EA)
            self.assertAlmostEqual(dx, ix, delta=1.0e-9 * L)
            self.assertAlmostEqual(dz, iz, delta=1.0e-9 * L)

    def test_lumped_chain_tends_to_the_cable(self):
        H, V0, L, EA = 13251.0, -2406.0, 301.5, 3.0e7
        dx, dz = bw.cable_reach(H, V0, L, W, EA)
        last = None
        for n in (20, 200, 2000):
            # the same pull on the support at end a: the chain's first link starts half a
            # link's weight further along
            lx, lz = bw.cable_reach(H, V0 + 0.5 * W * L / n, L, W, EA, n)
            err = math.hypot(lx - dx, lz - dz)
            if last is not None:
                self.assertLess(err, last / 50.0)      # second order in the link length
            last = err
        self.assertLess(last, 1.0e-5)
        # the supports of a lumped chain carry its whole weight
        fa, fb = bw.cable_end_pulls(H, V0, L, W, 20)
        self.assertAlmostEqual(fa[1] + fb[1], -W * L, delta=1.0e-9 * W * L)
        self.assertAlmostEqual(fa[0] + fb[0], 0.0, delta=1.0e-12)
        # a level span is symmetric: each support carries half of the chain's weight
        out = bw.solve_section(dict(DRAKE, end_a=[0.0, 30.0], end_b=[300.0, 30.0], lengths=[301.5], segments=20))
        self.assertAlmostEqual(out["pull_on_end_a"][1], -0.5 * W * 301.5, delta=1.0e-6)
        self.assertAlmostEqual(out["pull_on_end_b"][1], -0.5 * W * 301.5, delta=1.0e-6)
        self.assertAlmostEqual(out["pull_on_end_a"][0], -out["pull_on_end_b"][0], delta=1.0e-9)

    def test_level_span_is_the_elastic_catenary_of_the_docs(self):
        # c = H L / EA + (2 H / w) asinh(w L / 2 H); sag and end tension as printed there
        out = bw.solve_section(dict(DRAKE, end_a=[0.0, 30.0], end_b=[300.0, 30.0], lengths=[301.5]))
        s = out["spans"][0]
        H = s["horizontal_tension"]
        L, EA = 301.5, 3.0e7
        self.assertAlmostEqual(H * L / EA + 2.0 * H / W * math.asinh(W * L / (2.0 * H)), 300.0, delta=1.0e-8)
        self.assertAlmostEqual(H, 13251.040, delta=0.001)
        self.assertAlmostEqual(s["tension_a"], 13467.784, delta=0.001)
        self.assertAlmostEqual(s["tension_b"], s["tension_a"], delta=1.0e-6)
        self.assertAlmostEqual(s["low_point_below_a"], 13.583619, delta=1.0e-6)
        self.assertAlmostEqual(s["low_point_below_a"],
                               H / W * (math.sqrt(1.0 + (W * L / (2.0 * H)) ** 2) - 1.0) + W * L * L / (8.0 * EA),
                               delta=1.0e-9)

    def test_inclined_span_closes(self):
        H, V0 = bw.solve_cable((250.0, 40.0), 256.0, W, 3.0e7)
        ix, iz = reach_by_integration(H, V0, 256.0, W, 3.0e7)
        self.assertAlmostEqual(ix, 250.0, delta=1.0e-7)
        self.assertAlmostEqual(iz, 40.0, delta=1.0e-7)
        # the higher end carries more: by the weight of the height between them, to first order
        Ta, Tb = bw.cable_end_tensions(H, V0, 256.0, W)
        self.assertAlmostEqual(Tb - Ta, W * 40.0, delta=0.01 * W * 40.0)


RELEASES = [None, {"span": 1, "end": "b"}, {"span": 1, "end": "a"}, {"span": 2, "end": "a"},
            {"span": 2, "end": "b"}, {"span": 3, "end": "a"}, {"string": 1}, {"string": 2}]


class Balance(unittest.TestCase):
    def test_supports_carry_the_weight_and_no_net_pull(self):
        wi_total = (60.0 - 1.2 * 0.25 * math.pi * 0.254 ** 2 * 2.5) * 9.81     # one string (N)
        for seg in (0, 20):
            for rel in RELEASES:
                plan = section(segments=seg)
                if rel is not None:
                    plan["released"] = rel
                out = bw.solve_section(plan)
                fx = out["pull_on_end_a"][0] + out["pull_on_end_b"][0]
                fz = out["pull_on_end_a"][1] + out["pull_on_end_b"][1]
                for t in out["strings"]:
                    fx += t["pull_on_tower"][0]
                    fz += t["pull_on_tower"][1]
                weight = 2.0 * wi_total
                for k, (s, L) in enumerate(zip(out["spans"], plan["lengths"])):
                    if s["attached"] == "both":
                        weight += W * L
                        continue
                    # counted here, not taken from the tool: what hangs between the attachment the
                    # span is still on and the level the conductor rests at
                    held = out["attachments"][k if s["attached"] == "a" else k + 1]
                    height = held[1] - 0.5 * 0.0281
                    if seg == 0:
                        hangs = min(L, height * (1.0 - 0.5 * W * height / 3.0e7))
                        weight += W * hangs
                    else:
                        l = L / seg
                        nodes = sum(1 for i in range(1, seg) if i * l * (1.0 + W * i * l / 6.0e7) < height)
                        weight += W * l * (nodes + 0.5)
                self.assertAlmostEqual(fx, 0.0, delta=1.0e-5, msg="%s %s" % (seg, rel))
                self.assertAlmostEqual(-fz, weight, delta=1.0e-6 * weight, msg="%s %s" % (seg, rel))

    def test_energy_is_stationary_at_the_strings_feet(self):
        # span 1 released at tower 1 (it hangs from dead end a and takes no further part): the
        # energy of string 1, span 2, string 2 and span 3 as a function of the two feet
        plan = section(released={"span": 1, "end": "b"})
        out = bw.solve_section(plan)
        feet = [list(t["foot"]) for t in out["strings"]]
        tops = plan["towers"]
        wi = (60.0 / 2.5 - 1.2 * 0.25 * math.pi * 0.254 ** 2) * 9.81

        def energy(f):
            e = cable_energy(tops[0], f[0], 2.5, wi, bw.STRING_STIFFNESS)
            e += cable_energy(tops[1], f[1], 2.5, wi, bw.STRING_STIFFNESS)
            e += cable_energy(f[0], f[1], plan["lengths"][1], W, 3.0e7)
            e += cable_energy(f[1], plan["end_b"], plan["lengths"][2], W, 3.0e7)
            return e

        e0 = energy(feet)
        # along the strings their stiffness (4e6 N/m) dominates, so the foot is moved across them
        for i in range(2):
            ang = math.radians(out["strings"][i]["angle_deg"])
            across = (math.cos(ang), math.sin(ang))
            def moved(d):
                return energy([[f[0] + d * across[0] * (j == i), f[1] + d * across[1] * (j == i)]
                               for j, f in enumerate(feet)])

            # a straight move across a string stretches it at second order, which leaves an error
            # of the central difference in proportion to d^2: two steps remove it
            d = 1.0e-3
            ep, em = moved(d), moved(-d)
            slope = (4.0 * (ep - em) / (2.0 * d) - (moved(2.0 * d) - moved(-2.0 * d)) / (4.0 * d)) / 3.0
            # N: the net force across the string at the foot
            self.assertLess(abs(slope), 0.05, "foot %d: %g N against a tension of 12 kN" % (i + 1, slope))
            self.assertGreater(ep + em - 2.0 * e0, 0.0, "a minimum of the energy")

    def test_a_string_points_along_its_load(self):
        # moments about the string's top: the pull along the line at its foot times its height
        # equals the vertical load at its foot plus half its own weight, times its offset
        wi_total = (60.0 - 1.2 * 0.25 * math.pi * 0.254 ** 2 * 2.5) * 9.81
        for rel in (None, {"span": 1, "end": "b"}, {"span": 3, "end": "a"}):
            plan = section()
            if rel is not None:
                plan["released"] = rel
            out = bw.solve_section(plan)
            for t in out["strings"]:
                fx, fz = t["pull_on_tower"]
                foot_down = -fz - wi_total            # the vertical load at the foot
                expect = math.degrees(math.atan2(fx, foot_down + 0.5 * wi_total))
                # a straight string; its own sag changes the angle by a few thousandths of a degree
                self.assertAlmostEqual(t["angle_deg"], expect, delta=0.01, msg=str(rel))


class Limits(unittest.TestCase):
    def test_whole_section_strung_alike_hangs_its_strings_plumb(self):
        # three equal level spans: by symmetry the strings are plumb and each carries one span
        plan = dict(DRAKE, end_a=[0.0, 57.5], end_b=[900.0, 57.5], towers=[[300.0, 60.0], [600.0, 60.0]],
                    lengths=[301.5, 301.5, 301.5], insulator_length=2.5, insulator_mass=60.0)
        out = bw.solve_section(plan)
        wi_total = (60.0 - 1.2 * 0.25 * math.pi * 0.254 ** 2 * 2.5) * 9.81
        for t in out["strings"]:
            self.assertAlmostEqual(t["angle_deg"], 0.0, delta=0.01)
            self.assertAlmostEqual(-t["pull_on_tower"][1], W * 301.5 + wi_total, delta=1.0e-3 * W * 301.5)
        self.assertAlmostEqual(out["spans"][1]["horizontal_tension"], 13251.0, delta=15.0)

    def test_a_released_string_leaves_one_span_with_a_weight_on_it(self):
        # two equal spans, the string between them released: one cable of twice the length with
        # the string's weight at its middle; each half ends there with half that weight
        plan = dict(DRAKE, end_a=[0.0, 80.0], end_b=[600.0, 80.0], towers=[[300.0, 82.5]],
                    lengths=[301.5, 301.5], insulator_length=2.5, insulator_mass=60.0, released={"string": 1})
        out = bw.solve_section(plan)
        wi_total = (60.0 - 1.2 * 0.25 * math.pi * 0.254 ** 2 * 2.5) * 9.81
        self.assertAlmostEqual(out["attachments"][1][0], 300.0, delta=1.0e-6)
        H = out["spans"][0]["horizontal_tension"]
        self.assertAlmostEqual(H, out["spans"][1]["horizontal_tension"], delta=1.0e-6 * H)
        self.assertAlmostEqual(-out["pull_on_end_a"][1], W * 301.5 + 0.5 * wi_total, delta=1.0e-6 * W * 301.5)
        # without the string's weight: the level span of 603 m between the dead ends
        light = dict(plan, insulator_mass=0.2)
        one = bw.solve_section(dict(DRAKE, end_a=[0.0, 80.0], end_b=[600.0, 80.0], lengths=[603.0]))
        two = bw.solve_section(light)
        self.assertAlmostEqual(two["spans"][0]["horizontal_tension"], one["spans"][0]["horizontal_tension"],
                               delta=2.0e-3 * one["spans"][0]["horizontal_tension"])

    def test_a_released_span_hangs_by_its_weight(self):
        # the model's Eq. br-rest and br-lumped: 31 m of conductor from a support 30 m up
        h = 30.0 - 0.5 * 0.0281
        load, length, tension = bw.hanging_load(h, 31.0, W, 3.0e7)
        self.assertEqual(load, tension)
        self.assertLess(load, W * h)                                # the hanging length stretches by 0.24 mm
        self.assertAlmostEqual(load, W * h, delta=0.005)            # 478.68 N
        self.assertAlmostEqual(load, 478.68, delta=0.01)
        load, length, tension = bw.hanging_load(h, 31.0, W, 3.0e7, 31)
        self.assertAlmostEqual(tension, W * 29.0, delta=1.0e-9)     # 462.94 N: the 29 nodes clear of the ground
        self.assertAlmostEqual(load, W * 29.5, delta=1.0e-9)        # the support also carries the held end's node
        # shorter than the support's height: its whole weight (340.8 N in the docs)
        load, length, tension = bw.hanging_load(h, 21.35, W, 3.0e7)
        self.assertAlmostEqual(load, W * 21.35, delta=1.0e-9)
        self.assertAlmostEqual(load, 340.8, delta=0.05)
        self.assertEqual(length, 21.35)
        # in links, hanging clear: the support carries all of it, the first link all but the held end's node
        load, length, tension = bw.hanging_load(h, 21.0, W, 3.0e7, 3)
        self.assertAlmostEqual(load, W * 21.0, delta=1.0e-9)
        self.assertAlmostEqual(tension, W * 17.5, delta=1.0e-9)

    def test_a_released_span_swings_its_string_into_the_span_that_stays(self):
        # span 1 released at tower 1: the string there carries span 2 alone and leans towards
        # it; the tension in span 2 falls, and tower 1 now takes it along the line
        whole = bw.solve_section(section())
        cut = bw.solve_section(section(released={"span": 1, "end": "b"}))
        t = cut["strings"][0]
        self.assertGreater(t["angle_deg"], 45.0)
        self.assertLess(cut["spans"][1]["horizontal_tension"], whole["spans"][1]["horizontal_tension"])
        self.assertAlmostEqual(t["pull_on_tower"][0], cut["spans"][1]["horizontal_tension"], delta=1.0e-6)
        self.assertLess(abs(whole["strings"][0]["pull_on_tower"][0]), 0.01 * whole["spans"][1]["horizontal_tension"])


class Robust(unittest.TestCase):
    def test_spans_strung_unalike_lean_their_strings_far(self):
        # short strings between spans that pull unequally (4.2, 11.8 and 8.7 kN): the strings lean
        # 52 and -24 degrees, far from the plumb start, which Newton's method alone does not reach
        plan = dict(DRAKE, end_a=[0.0, 48.103], end_b=[903.387, 47.603], towers=[[150.467, 48.103], [639.34, 48.103]],
                    lengths=[152.71, 494.886, 265.893], axial_stiffness=2591081.68, insulator_length=0.5,
                    insulator_mass=150.0)
        out = bw.solve_section(plan)
        self.assertAlmostEqual(out["strings"][0]["angle_deg"], 52.2, delta=0.1)
        self.assertAlmostEqual(out["strings"][1]["angle_deg"], -23.96, delta=0.1)
        # checked without the tool's balance: no net pull along the line on the supports, and the
        # difference of the tensions either side of a string is what its lean carries
        fx = out["pull_on_end_a"][0] + out["pull_on_end_b"][0] + sum(t["pull_on_tower"][0] for t in out["strings"])
        self.assertAlmostEqual(fx, 0.0, delta=1.0e-5)
        H = [sp["horizontal_tension"] for sp in out["spans"]]
        self.assertAlmostEqual(out["strings"][0]["pull_on_tower"][0], H[1] - H[0], delta=1.0e-5)
        self.assertAlmostEqual(out["strings"][1]["pull_on_tower"][0], H[2] - H[1], delta=1.0e-5)
        # each span closes on its two attachments, by integration of its shape
        att = out["attachments"]
        for k, L in enumerate(plan["lengths"]):
            reach = (att[k + 1][0] - att[k][0], att[k + 1][1] - att[k][1])
            Hk, V0 = bw.solve_cable(reach, L, W, plan["axial_stiffness"])
            self.assertAlmostEqual(Hk, H[k], delta=1.0e-6 * H[k])
            ix, iz = reach_by_integration(Hk, V0, L, W, plan["axial_stiffness"])
            self.assertAlmostEqual(ix, reach[0], delta=1.0e-6)
            self.assertAlmostEqual(iz, reach[1], delta=1.0e-6)

    def test_a_conductor_on_the_ground_is_outside_the_tool(self):
        # 300, 400 and 200 m on 30 m towers: with the string at tower 1 released the 700 m span
        # would lie on the ground while held at both ends
        plan = dict(DRAKE, end_a=[100.0, 30.0], end_b=[1000.0, 30.0], towers=[[400.0, 30.0], [800.0, 30.0]],
                    lengths=[301.5, 402.0, 201.0], insulator_length=2.5, insulator_mass=60.0)
        bw.solve_section(plan)                      # whole: clear of the ground
        with self.assertRaises(bw.OutsideTool) as cm:
            bw.solve_section(dict(plan, released={"string": 1}))
        self.assertIn("on the ground", str(cm.exception))
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "plan.json")
            with open(path, "w") as f:
                json.dump(dict(plan, released={"string": 1}), f)
            with contextlib.redirect_stderr(io.StringIO()) as err:
                self.assertEqual(bw.main([path]), 4)
            self.assertIn("clear of the ground", err.getvalue())


class Reference(unittest.TestCase):
    """The numbers the model's unit test BrokenWire.* (Structures/Tests/Unit/ERF_GTestBrokenWire.cpp)
    compares MoorDyn-C at rest with. A change here must be made there too."""

    def test_printed(self):
        buf = io.StringIO()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "plan.json")
            with open(path, "w") as f:
                json.dump(section(released={"span": 1, "end": "b"}), f)
            with contextlib.redirect_stdout(buf):
                self.assertEqual(bw.main([path]), 0)
            self.assertIn("hangs", buf.getvalue())
            with contextlib.redirect_stdout(io.StringIO()) as js:
                self.assertEqual(bw.main([path, "--json"]), 0)
            self.assertEqual(json.loads(js.getvalue())["spans"][0]["attached"], "a")


class Refused(unittest.TestCase):
    def refuse(self, text, **change):
        plan = section()
        for k, v in change.items():
            if v is None:
                plan.pop(k, None)
            else:
                plan[k] = v
        with self.assertRaises(bw.PlanError) as cm:
            bw.solve_section(plan)
        self.assertIn(text, str(cm.exception))

    def test_every_key(self):
        self.refuse("end_a", end_a=[1.0])
        self.refuse("end_b", end_b="x")
        self.refuse("towers", towers=[[400.0]])
        self.refuse("increasing x", towers=[[900.0, 60.0], [400.0, 60.0]])
        self.refuse("lengths", lengths=[301.5, 506.944])
        self.refuse("lengths", lengths=[301.5, 0.0, 100.0])
        self.refuse("diameter", diameter=0.0)
        self.refuse("diameter", diameter=None)
        self.refuse("mass_per_length", mass_per_length=-1.0)
        self.refuse("axial_stiffness", axial_stiffness=float("nan"))
        self.refuse("insulator_length", insulator_length=-2.5)
        self.refuse("insulator_mass", insulator_mass=0.0)
        self.refuse("insulator_mass", insulator_length=0.0)
        self.refuse("insulator_diameter", insulator_diameter=0.0)
        self.refuse("air_density", air_density=-1.0)
        self.refuse("air_density", air_density=3000.0)          # lighter than the air displaced
        self.refuse("gravity", gravity=0.0)
        self.refuse("ground", ground=57.6)                      # above the strings' feet
        self.refuse("segments", segments=-1)
        self.refuse("segments", segments=2.5)
        self.refuse("segments", segments=100001)
        self.refuse("released", released={"span": 4, "end": "a"})
        self.refuse("released", released={"span": 1, "end": "c"})
        self.refuse("released", released={"span": 1})
        self.refuse("released", released={"string": 3})
        self.refuse("unknown key", span_lengths=[1.0])

    def test_string_release_needs_strings(self):
        plan = dict(DRAKE, end_a=[0.0, 60.0], end_b=[600.0, 60.0], towers=[[300.0, 60.0]],
                    lengths=[301.5, 301.5], released={"string": 1})
        with self.assertRaises(bw.PlanError) as cm:
            bw.solve_section(plan)
        self.assertIn("insulator_length", str(cm.exception))

    def test_command_exit_codes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "plan.json")
            with open(path, "w") as f:
                json.dump(section(diameter=-1.0), f)
            with contextlib.redirect_stderr(io.StringIO()) as err:
                self.assertEqual(bw.main([path]), 2)
            self.assertIn("diameter", err.getvalue())
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(bw.main([os.path.join(tmp, "missing.json")]), 2)


if __name__ == "__main__":
    names = [a for a in sys.argv[1:] if not a.startswith("-")]
    if names:
        suite = unittest.TestSuite()
        for name in names:
            suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(globals()[name]))
        sys.exit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
    unittest.main(verbosity=2)
