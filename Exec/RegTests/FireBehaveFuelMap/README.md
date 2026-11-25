# FireBehaveFuelMap

## Purpose

Checks that `erf.fire.ros_model = "behave"` with
`erf.fire.rothermel_per_fuel = true` spreads each cell of a fuel map with that
cell's own fuel model. Before this test, BEHAVE used the uniform
`erf.fire.fuel_model_id` on every cell, whatever the map said.

The test also pins the multi-class Rothermel form of BEHAVE. The fuel-bed
SAV is the characteristic SAV over the dead and live categories, and W' and
M'_f are summed per size class. With the earlier dead-only SAV, Scott-Burgan
SH5 spread at 0.65 m/s instead of 1.76 m/s.

The map has Scott-Burgan SH5 (code 145) west of x = 500 m and GR2 (code 102)
east of it. The uniform fuel is GR2, so a run that ignores the map gives GR2's
rate on the SH5 half. The atmosphere is still, and the fire sees a prescribed
6.5 m/s wind along +x with the wind adjustment factor off. The ground is
flat. The moistures are fixed: 3, 4 and 5 % dead, 60 % live.

Two straight ignition lines, one at the western wall (SH5) and one at
x = 512.5 m (GR2), make the level-set function the distance to two straight
fronts that start at x = 25 m and x = 525 m. East of each line the front
normal is +x in every cell, so each front crosses its own fuel at that
fuel's head rate on both paths. `fire_ros` holds the isotropic head rate on
both paths; the rate the level-set front moves with on the directional path
shows only in `fire_arrival_time`, so the test checks both.

## Files

| File | Role |
|------|------|
| `inputs_base` | Shared settings: domain, fuel map, moistures, prescribed wind |
| `inputs_isotropic` | `erf.fire.directional_ros = false` (`fill_behave_ros`) |
| `inputs_directional` | `erf.fire.directional_ros = true` (`fill_directional_ros`, projection shape) |
| `fuel_map_sh5_gr2.asc` | 40 x 20 fire-grid map, 25 m cells: 145 for x < 500 m, 102 east of it |
| `ignition_line.csv` | North-south ignition polyline at x = 12.5 m (SH5) |
| `ignition_line_gr2.csv` | North-south ignition polyline at x = 512.5 m (GR2) |
| `check_behave_fuel_map.py` | Reads `fire_ros` and `fire_fuel_model_code` at step 20 and `fire_arrival_time` at the last step, and compares them with the reference values |
| `reference_rates.py` | Independent closed-form multi-class Rothermel rates for SH5 and GR2 (prints the values the check pins) |

## Data and acknowledgements

No external data. The fuel models are Scott and Burgan (2005), *Standard Fire
Behavior Fuel Models: A Comprehensive Set for Use with Rothermel's Surface
Fire Spread Model*, USDA Forest Service RMRS-GTR-153, already tabulated in
`FireEngine/ERF_FuelModels.H`.

The reference rates come from `reference_rates.py`, a short evaluation of
Rothermel (1972, INT-115) with Albini's (1976, INT-30) multi-class
weighting, as Andrews (2018, RMRS-GTR-371) sets it out. It was written apart
from the FLARE code:

| Fuel | Code | Rate at 6.5 m/s (m/s) | Note |
|------|------|-----------------------|------|
| SH5 | 145 | 1.7584449850336146 | |
| GR2 | 102 | 0.9489605773204224 | The 0.9 I_R wind cap binds |

GR2's live herbaceous load is two thirds cured at 60 % live moisture
(transfer window 0.30-1.20).

## How to run

```
mpirun -np 2 /path/to/flare inputs_isotropic
python3 check_behave_fuel_map.py isotropic
mpirun -np 2 /path/to/flare inputs_directional
python3 check_behave_fuel_map.py directional
```

Each run takes 240 steps of 0.5 s (120 s) and writes a fire plotfile every
20 steps. The check reads `plt_fire_<variant>_00020` for the rates and the
last plotfile for the arrival times. `erf_plotfile.py` (from
`Submodules/ERF/Exec/Tools`) is found from the source folder; CTest copies it beside the script.
`python3 reference_rates.py` prints the reference rates.

## Expected results

| Check | Isotropic | Directional |
|-------|-----------|-------------|
| `fire_fuel_model_code` is 145 west and 102 east of x = 500 m | PASS | PASS |
| SH5 half (320 cells, x from 100 to 500 m): `fire_ros` = 1.7584449850 | rel. err 6e-16 | rel. err 6e-16 |
| GR2 half (400 cells): `fire_ros` = 0.9489605773 | rel. err 1e-15 | rel. err 1e-15 |
| SH5 front (140 cells, x from 50 to 225 m): arrival (x - 25)/R_SH5, tolerance 7.1 s | max err 3.2 s | max err 3.0 s |
| GR2 front (80 cells, x from 525 to 625 m): arrival (x - 525)/R_GR2, tolerance 13.2 s | max err 2.7 s | max err 1.8 s |

The rate tolerance is 1e-8; the arrival tolerance is half the time the
front takes to cross a 25 m fire cell. On the code before this test, both
variants fail the rate checks:

| Code | SH5 half | GR2 half |
|------|----------|----------|
| No per-fuel BEHAVE table | 0.9490 (the uniform GR2) | 0.9490 |
| Per-fuel table, dead-only characteristic SAV | 0.6264 | 0.9602 |

With the per-fuel table on the isotropic path but not on the level-set path
(`st.bs_tbl = nullptr` in `make_directional_state`), the directional
variant still passes the rate checks but fails the SH5 front: 80 of the 140
SH5 cells are not reached in 120 s, as the front moves at GR2's rate.

## CTest

```
ctest -R FireBehaveFuelMap
```

`FireBehaveFuelMap_isotropic` and `FireBehaveFuelMap_directional` run on two
ranks, one box each, through `add_test_fire_check`. The unit tests
`BehaveFuelMap.*` (`Tests/Unit/Fire/ERF_GTestBehaveFuelMap.cpp`) pin the same
reference values for the state, the slot table and the three kernels.
