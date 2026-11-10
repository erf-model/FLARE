# Conductors_TowerFireHeating

## Purpose

Checks the towers' steel temperature under a fire (`conductors.thermal.towers = true`, theory in
`Docs/sphinx_doc/theory/Conductors.rst`, "Tower steel temperature under a fire"). One line of two
150 m spans hangs from a suspension tower at (750, 500) between two fixed ends, 30 m above flat
ground, in still air at 300 K. The tower bends as a frame ERF generates from its type: 8 panels of
equal-leg angles below the cross-arm (legs 150 x 12 mm, the rest 100 x 8 mm), crossed bracing, a
cross-arm truss each side; 103 nodes, 310 members, 11.1 t. FLARE's fire engine releases a
prescribed 2e5 W/m^2 over a 20 m disc under the tower from t = 0 (`erf.fire.prescribed_heat.*`;
nothing else burns). The disc's centre stands 5 m along the line and 10 m along the cross-arm from
the tower's axis, so no two mirrored members heat alike and a member put on the wrong side of the
tower shows. The coupling is passive, so the air keeps the sounding's temperature, and each
member's temperature follows its EN 1993-1-2 heat balance alone. Each burning fire cell radiates
0.3 of its heat release from 2 m above the ground. The members near the fire pass 100 C, where the
steel's stiffness starts to fall, and the tower's frame is rebuilt at their temperatures. The deck
sets the steel's emissivity (0.6), the convection coefficient (20 W/(m^2 K)) and the rebuild
tolerance (0.02) away from their defaults, so the check shows that each key is read.

The conductors' own thermal model (`conductors.thermal_model = ieee738`) runs as in
Conductors_FireHeating; the towers' heating needs it.

## Files

- `Conductors_TowerFireHeating.i`: the deck, 1500 x 1000 x 300 m on 30 x 20 x 12 cells, 150 steps of
  2 s, fire cells of 10 m (grid ratio 5).
- `input_sounding`: still air, potential temperature 300 K, surface pressure 1000 hPa.
- `check_tower_heating.py`: the check. Pure Python, no numpy, and no code shared with ERF. From the
  deck's numbers and the frame the run wrote (`conductors/frame_S1_t1.dat`, a SubDyn file, and
  `conductors/frame_S1_t1_members.dat`, the members' angles) it recomputes:
  1. each member's three Gauss points and axis in ERF's axes (the tower on its `ground.dat` point, the
     frame's x axis along the line, its y axis along the cross-arm);
  2. the air temperature there, T = theta - g z / c_p, and the fire's irradiance from the fire cells
     whose centres lie in the disc, as Gauss means;
  3. every member's temperature, stepped with backward Euler in the steel's enthalpy, which the script
     integrates numerically from EN 1993-1-2's specific heat (section 3.4.1.2);
  4. the frame's rebuilds, by the rule of the theory section with EN 1993-1-2's Table 3.1.

  Every member's temperature in `tower_S1_t1_member_temperature.dat` must match to 1e-3 K at every
  step, and every row of `tower_temperature.dat` (mean, largest, its member, largest irradiance,
  least k_E, rebuilds, the frame's largest temperature) must match. The run must rebuild the frame
  at least once and heat the hottest member by at least 50 K, so the checks cannot pass on a cold
  tower.

## Data and acknowledgements

No external data. The steel's heating method and reduction factors are those of EN 1993-1-2 as the
open-access JRC report of Vassart et al. (2014), EUR 26698 EN, doi:10.2788/85432, reproduces them (with
thanks to the JRC and the authors); the tower's dimensions are the generic ones of
Conductors_GeneratedTowers.

## How to run

Build FLARE with `-DERF_ENABLE_STRUCTURES=ON -DERF_MOORDYN_USE_STUB=ON -DERF_ENABLE_FIRE=ON`, FFT and
MPI, then

    mpiexec -n 2 erf_exec Conductors_TowerFireHeating.i
    python3 check_tower_heating.py --dt 2 --fire-lags 1 --fire-dx 10

The two-level variants add the options listed in `Structures/Tests/ERFStructuresTests.cmake` (the
same as Conductors_FireHeating's) and pass `--dt 1 --fire-lags 0 --fire-dx 10` (the tower on level 1)
or `--dt 2 --fire-lags 1 --fire-dx 5` (the fire on level 1).

## Expected results

| Run | Fire / tower | Hottest member at 300 s | Members' mean | Largest member irradiance | Rebuilds | Least k_E | Largest difference |
|---|---|---|---|---|---|---|---|
| `Conductors_TowerFireHeating` | level 0 / level 0 | 76: 26.8 to 130.73 C | 57.35 C | 46907 W/m^2 | 1 (at 260 s) | 0.980 | 1.0e-4 K |
| `…_TowerOnFinerLevel` | level 0 / level 1 | 76: 26.8 to 131.30 C | 57.52 C | 46907 W/m^2 | 1 | 0.980 | 2.3e-5 K |
| `…_FireOnFinerLevel` | level 1 (5 m cells) / level 0 | 76: 26.8 to 164.63 C | 61.68 C | 59574 W/m^2 | 3 | 0.940 | 9.2e-5 K |

Member 76 is a 100 x 8 mm diagonal of the bottom panel. Its section factor is 260 m^-1, and its time
constant is about 450 s near ambient and 350 s at 250 C, so it is still heating at 300 s. The irradiance agrees to 1e-10.
The largest difference in the temperatures comes from the air temperature: ERF's equation of state,
extrapolated below the first cell centre, against the formula, up to 1e-4 K; it does not grow over
the run.

The fire on the refined level heats the members more. Its 5 m cells place 52 sources over the disc
instead of 12, nearer the bottom members, and floor the squared distance at dA/pi = 8 m^2 instead of
32 m^2. The point-source model's near field depends on the fire's cell size; the docs list it among
the limitations.

## CTest

`Conductors_TowerFireHeating`, `Conductors_TowerFireHeating_TowerOnFinerLevel`,
`Conductors_TowerFireHeating_FireOnFinerLevel` (deck and check), `Conductors_TowerFireHeating_Restart`
(a checkpoint at 270 s, after the first rebuild: the temperature logs and `towers.dat` must equal the
uninterrupted run's), `Conductors_TowerFireHeating_RestartUnheatedAbort` (a restart without the
heating from that checkpoint stops, naming `conductors.thermal.towers`), and `Conductors_Thermal_TowerSteelTemperatureAbort`,
`Conductors_Thermal_TowerRangeAbort` and `Conductors_Thermal_TowerMeltsAbort` (the deck with a type's
`steel_temperature`, a zero emissivity, and a fire 500 times stronger, which takes a member past
1200 C), all with the stub, in `Structures/Tests/ERFStructuresTests.cmake`.
