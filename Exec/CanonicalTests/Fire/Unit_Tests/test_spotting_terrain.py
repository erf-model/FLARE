#!/usr/bin/env python3
"""
Terrain-aware firebrand trajectory unit tests.

A firebrand is lofted a height H_z above the ground it starts on and falls at
its terminal velocity while the wind carries it. The ground it lands on is not
the height it left: a brand drifting downslope falls further and lands further
out, and one drifting upslope meets the rising ground sooner. These tests cover
the descent, its reduction to the flat-earth trajectory when there is no
terrain, and a step that depends on the brand's own loft, not on the lowest
ground in the domain.

Reference implementation extracted from:
  - ERF_AlbiniSpotting.H: compute_albini_spotting() trajectory loop,
                          albini_lofting_height()

Run: python3 test_spotting_terrain.py
"""

import sys
import math


def lofting_height(I_B_kW, I_B_min=100.0, U10=5.0, A=83.9, B=-1.040):
    """
    Albini (1983, INT-309) lofting height, H_z = 0.173 (I_B A U10^B)^(1/2),
    U10 held to 5-30 m/s, zero below threshold (default A, B: Anderson 4).
    """
    if I_B_kW < I_B_min:
        return 0.0
    U = min(max(U10, 5.0), 30.0)
    return 0.173 * math.sqrt(I_B_kW * A * U ** B)


def fly(H_z, z_src, ground_of_x, u_wind, w_terminal, n_steps=20,
        z_ground_min=0.0, x0=0.0):
    """
    Integrate a firebrand from its source to the ground.

    Mirrors the trajectory loop: each step falls H_z / n_steps, and over
    lower ground the brand steps on, up to its fall to the lowest ground in
    the domain (n_fall steps), so it always lands. In the step where its
    height above the ground goes from a0 > 0 to a1 < 0 it lands at the
    fraction a0 / (a0 - a1) of the step (until Nov 2026 the step was the fall
    to the lowest ground over n_steps, and the brand landed at the step end).

    Returns (landing_x, flight_time, steps_taken).
    """
    z_brand = z_src + H_z
    dz = H_z / n_steps
    dt = dz / w_terminal
    n_fall = max(n_steps, math.ceil((z_brand - z_ground_min) / dz))

    x = x0
    for ns in range(n_fall):
        x_start = x
        a0 = z_brand - ground_of_x(x)
        x += u_wind * dt
        z_brand -= dz
        a1 = z_brand - ground_of_x(x)
        if a1 <= 0.0:
            frac = a0 / (a0 - a1) if a1 < 0.0 < a0 else 1.0
            x = x_start + (x - x_start) * frac
            return x, (ns + frac) * dt, ns + 1
    return x, n_fall * dt, n_fall


def check(failures, name, condition, detail=""):
    status = "✓" if condition else "✗"
    print(f"{status} {name}")
    if not condition:
        if detail:
            print(f"    {detail}")
        failures.append(name)


def test_flat_matches_the_old_flight_time():
    """With no terrain the trajectory is the flat-earth one, unchanged."""
    failures = []
    H_z, u, w = 60.0, 8.0, 1.5
    flat = lambda x: 0.0

    x_land, t_flight, steps = fly(H_z, 0.0, flat, u, w, n_steps=20)

    # The superseded model flew for exactly H_z / w_terminal and drifted u * t
    t_expected = H_z / w
    check(failures, "Flat terrain reproduces the flat-earth flight time",
          abs(t_flight - t_expected) < 1.0e-12,
          f"got {t_flight}, expected {t_expected}")
    check(failures, "Flat terrain reproduces the flat-earth landing point",
          abs(x_land - u * t_expected) < 1.0e-9,
          f"got {x_land}, expected {u * t_expected}")
    check(failures, "Flat terrain uses every trajectory step", steps == 20,
          f"took {steps} steps")
    return failures


def test_downslope_flies_further():
    """Drifting onto lower ground extends the fall and the reach."""
    failures = []
    H_z, u, w = 60.0, 8.0, 1.5

    # Ground descending at 10% in the direction of drift
    downslope = lambda x: -0.10 * x
    flat = lambda x: 0.0

    x_flat, t_flat, _ = fly(H_z, 0.0, flat, u, w, n_steps=400, z_ground_min=0.0)
    x_down, t_down, _ = fly(H_z, 0.0, downslope, u, w, n_steps=400,
                            z_ground_min=-200.0)

    check(failures, "Downslope drift lands further than flat",
          x_down > x_flat + 1.0,
          f"downslope {x_down:.1f} m vs flat {x_flat:.1f} m")
    check(failures, "Downslope drift stays aloft longer",
          t_down > t_flat,
          f"{t_down:.1f} s vs {t_flat:.1f} s")
    return failures


def test_upslope_lands_sooner():
    """Rising ground cuts the fall short."""
    failures = []
    H_z, u, w = 60.0, 8.0, 1.5

    upslope = lambda x: 0.10 * x
    flat = lambda x: 0.0

    x_flat, t_flat, _ = fly(H_z, 0.0, flat, u, w, n_steps=400)
    x_up, t_up, _ = fly(H_z, 0.0, upslope, u, w, n_steps=400)

    check(failures, "Upslope drift lands closer than flat",
          x_up < x_flat - 1.0, f"upslope {x_up:.1f} m vs flat {x_flat:.1f} m")
    check(failures, "Upslope drift spends less time aloft",
          t_up < t_flat, f"{t_up:.1f} s vs {t_flat:.1f} s")
    return failures


def test_ridge_to_valley_reach():
    """A brand leaving a ridge crest reaches well beyond the flat estimate."""
    failures = []
    # 60 m loft from a 150 m crest, falling into ground 150 m lower
    H_z, u, w = 60.0, 12.0, 1.5
    z_src = 150.0
    valley = lambda x: max(0.0, 150.0 - 0.30 * x)

    x_land, t_land, _ = fly(H_z, z_src, valley, u, w, n_steps=800,
                            z_ground_min=0.0)
    x_flat, t_flat, _ = fly(H_z, 0.0, lambda x: 0.0, u, w, n_steps=800)

    check(failures, "Ridge-to-valley brand outruns the flat-earth estimate",
          x_land > 1.5 * x_flat,
          f"terrain {x_land:.0f} m vs flat {x_flat:.0f} m")
    return failures


def test_step_does_not_depend_on_the_lowest_ground():
    """On a plateau the flight is the flat one, however low the lowest ground."""
    failures = []
    H_z, u, w = 60.0, 8.0, 1.5
    n = 20
    plateau = lambda x: 500.0

    # The step is H_z / n, so the brand falls to the plateau in n steps; until
    # Nov 2026 it was (500 + H_z) / n here and the brand took 3 steps of 28 m
    x_land, t_land, steps = fly(H_z, 500.0, plateau, u, w, n_steps=n,
                                z_ground_min=0.0)
    check(failures, "Plateau flight takes the flat step count", steps == n,
          f"took {steps} of {n}")
    check(failures, "Plateau flight lands at the flat-earth point",
          abs(x_land - u * H_z / w) < 1.0e-9,
          f"got {x_land}, expected {u * H_z / w}")

    # A brand from the plateau over a cliff down to the lowest ground steps on
    # past n until it lands there
    cliff = lambda x: 500.0 if x < 10.0 else 0.0
    x_cliff, _, steps_cliff = fly(H_z, 500.0, cliff, u, w, n_steps=n,
                                  z_ground_min=0.0)
    check(failures, "Brand over a cliff steps on and lands on the low ground",
          steps_cliff > n and abs(x_cliff - u * (500.0 + H_z) / w) < 1.0e-6,
          f"took {steps_cliff} steps, landed at {x_cliff}")
    return failures


def test_lofting_height_threshold():
    """No loft below the minimum intensity; square-root growth above it."""
    failures = []
    check(failures, "No loft below the intensity threshold",
          lofting_height(50.0) == 0.0)
    h1 = lofting_height(100.0)
    h8 = lofting_height(800.0)
    h4 = lofting_height(400.0)
    check(failures, "Loft grows as the square root of intensity",
          abs(h4 / h1 - 2.0) < 1.0e-12, f"ratio {h4 / h1}")
    check(failures, "Loft falls as the wind rises (B < 0)",
          lofting_height(800.0, U10=20.0) < h8)
    return failures


def main():
    print("=" * 70)
    print("Terrain-Aware Firebrand Trajectory Unit Tests")
    print("=" * 70)

    all_failures = []
    all_failures.extend(test_flat_matches_the_old_flight_time())
    all_failures.extend(test_downslope_flies_further())
    all_failures.extend(test_upslope_lands_sooner())
    all_failures.extend(test_ridge_to_valley_reach())
    all_failures.extend(test_step_does_not_depend_on_the_lowest_ground())
    all_failures.extend(test_lofting_height_threshold())

    print("\n" + "=" * 70)
    if not all_failures:
        print("FINAL RESULT: All tests PASSED")
        return 0
    print(f"FINAL RESULT: {len(all_failures)} test(s) FAILED")
    for failure in all_failures:
        print(f"  ✗ {failure}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
