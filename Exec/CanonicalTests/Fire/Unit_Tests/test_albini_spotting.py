#!/usr/bin/env python3
"""
Unit tests for the Albini (1983) firebrand spotting model.

Checks the lofting height against the two worked examples of INT-309, the
trajectory distance scaling and the intensity threshold.

Reference:
  Albini, F.A. (1983). Potential spotting distance from wind-driven
    surface fires. USDA Forest Service Research Paper INT-309.
    doi:10.2737/INT-RP-309. Eq. 2: H = 0.173 E^(1/2); eqs. 7-8:
    E = I A U^B with A, B from Table 4 and U the 10 m wind.
    Example 1: short grass, 2000 kW/m, 5 m/s -> E = 30 660 kJ/m, H = 30 m.
    Example 2: mature chaparral, 50 000 kW/m, 20 m/s -> E = 186 000 kJ/m, H = 75 m.

Run: python3 test_albini_spotting.py
"""

import math
import sys

# INT-309 Table 4 (A [s], B), Anderson models; 8 is not in the table
TABLE4 = {1: (162.0, -1.456), 2: (166.0, -1.282), 3: (129.0, -1.238), 9: (154.0, -1.151),
          4: (83.9, -1.040), 5: (94.4, -0.881), 6: (74.1, -0.906), 7: (75.1, -0.884),
          10: (66.1, -0.810), 11: (63.6, -0.776), 12: (69.6, -0.818), 13: (71.8, -0.835)}


def albini_lofting_height(I_B_kW, U10=5.0, model=4, I_B_min=100.0, AB=None):
    """H = 0.173 (I A U^B)^(1/2), U held to the 5-30 m/s fit range (INT-309 eqs. 2, 7, 8)"""
    if I_B_kW < I_B_min:
        return 0.0
    A, B = AB if AB is not None else TABLE4[model]
    U = min(max(U10, 5.0), 30.0)
    return 0.173 * math.sqrt(I_B_kW * A * U ** B)


def albini_flight_time(H_z_m, v_t_ms=0.5):
    """t_f = H_z / v_t"""
    return H_z_m / v_t_ms if v_t_ms > 0 else 0.0


def albini_spot_distance_uniform_wind(I_B_kW, u_ms, v_t_ms=0.5):
    """
    Analytical spotting distance for uniform horizontal wind.
    d = u * t_f = u * H_z / v_t
    """
    H_z = albini_lofting_height(I_B_kW)
    t_f = albini_flight_time(H_z, v_t_ms)
    return u_ms * t_f


def test_lofting_height_threshold():
    """H_z = 0 below I_B_min threshold."""
    failures = []
    if albini_lofting_height(50.0, I_B_min=100.0) != 0.0:
        failures.append("I_B=50 kW/m should give H_z=0 (below threshold)")
    if albini_lofting_height(99.9, I_B_min=100.0) != 0.0:
        failures.append("I_B=99.9 kW/m should give H_z=0 (below threshold)")
    if albini_lofting_height(100.0, I_B_min=100.0) <= 0.0:
        failures.append("I_B=100 kW/m should give H_z > 0 (at threshold)")
    return failures


def test_lofting_height_formula():
    """
    INT-309's worked examples. Example 1 uses B = -1.465 for short grass
    (Table 4 prints -1.456) and gives 30 m; example 2 gives 75 m.
    """
    failures = []
    cases = [("example 1", albini_lofting_height(2000.0, 5.0, AB=(162.0, -1.465)), 30.0),
             ("example 2", albini_lofting_height(50000.0, 20.0, model=4), 75.0)]
    for name, H_z, H_ref in cases:
        if abs(H_z - H_ref) >= 0.5:
            failures.append(f"{name}: H_z={H_z:.2f} m, INT-309 gives {H_ref} m")
    return failures


def test_lofting_height_monotone():
    """H_z must increase monotonically with I_B."""
    I_values = [100.0, 200.0, 500.0, 1000.0, 5000.0]
    H_values = [albini_lofting_height(I) for I in I_values]
    failures = []
    for k in range(len(H_values) - 1):
        if H_values[k] >= H_values[k + 1]:
            failures.append(
                f"H_z not monotone: H({I_values[k]})={H_values[k]:.2f} >= "
                f"H({I_values[k+1]})={H_values[k+1]:.2f}"
            )
    return failures


def test_spot_distance_uniform_wind():
    """
    Spotting distance in uniform wind: d = u * H_z / v_t.
    For I_B=500 kW/m in FM4 at 5 m/s, v_t=0.5 m/s:
      H_z = 0.173 (500 83.9 5^-1.04)^(1/2) = 15.3 m, t_f = 30.7 s, d = 153 m
    """
    failures = []
    d = albini_spot_distance_uniform_wind(500.0, u_ms=5.0, v_t_ms=0.5)
    if abs(d - 153.4) >= 1.0:
        failures.append(f"Spot distance {d:.1f} m, expected 153.4 m")
    return failures


def test_spot_distance_scales_with_wind():
    """Spotting distance scales linearly with wind speed."""
    failures = []
    d1 = albini_spot_distance_uniform_wind(500.0, u_ms=2.0)
    d2 = albini_spot_distance_uniform_wind(500.0, u_ms=4.0)
    ratio = d2 / d1
    if abs(ratio - 2.0) >= 0.01:
        failures.append(
            f"Distance ratio {ratio:.4f} should be 2.0 (linear in wind speed)"
        )
    return failures


def test_spot_distance_zero_below_threshold():
    """No spotting when I_B below threshold."""
    failures = []
    d = albini_spot_distance_uniform_wind(50.0, u_ms=10.0, v_t_ms=0.5)
    if d != 0.0:
        failures.append(f"Spot distance below threshold should be 0, got {d}")
    return failures


def test_terminal_velocity_effect():
    """Higher terminal velocity reduces flight time and spot distance."""
    failures = []
    d_slow = albini_spot_distance_uniform_wind(500.0, u_ms=5.0, v_t_ms=0.3)
    d_fast = albini_spot_distance_uniform_wind(500.0, u_ms=5.0, v_t_ms=1.0)
    if d_slow <= d_fast:
        failures.append(
            f"Slow fall ({d_slow:.1f} m) should give farther spotting than "
            f"fast fall ({d_fast:.1f} m)"
        )
    return failures


def main():
    print("=" * 70)
    print("Albini (1983) Spotting Model Unit Tests")
    print("Reference: Albini INT-309, USDA Forest Service, 1983")
    print("=" * 70)

    tests = [
        ("test_lofting_height_threshold",    test_lofting_height_threshold),
        ("test_lofting_height_formula",      test_lofting_height_formula),
        ("test_lofting_height_monotone",     test_lofting_height_monotone),
        ("test_spot_distance_uniform_wind",  test_spot_distance_uniform_wind),
        ("test_spot_distance_scales_with_wind", test_spot_distance_scales_with_wind),
        ("test_spot_distance_zero_below_threshold", test_spot_distance_zero_below_threshold),
        ("test_terminal_velocity_effect",    test_terminal_velocity_effect),
    ]

    failed = []
    for name, test_fn in tests:
        try:
            failures = test_fn()
            if failures:
                for msg in failures:
                    print(f"✗ {name} FAILED: {msg}")
                failed.append((name, failures[0]))
            else:
                print(f"✓ {name} PASSED")
        except Exception as e:
            print(f"✗ {name} ERROR: {e}")
            failed.append((name, str(e)))

    print("=" * 70)
    if not failed:
        print(f"All {len(tests)} tests PASSED")
        return 0
    else:
        print(f"{len(failed)} test(s) FAILED:")
        for name, error in failed:
            print(f"  - {name}: {error}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
