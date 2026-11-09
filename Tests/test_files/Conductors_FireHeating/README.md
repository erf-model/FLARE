# Conductors_FireHeating

## Purpose

Checks the conductors' temperature under a fire (`conductors.thermal_model = ieee738`, theory in
`Docs/sphinx_doc/theory/Conductors.rst`, "Conductor temperature under a fire"). A 300 m span of a
Drake-like ACSR conductor (28.1 mm, 1.628 kg/m, 1.5 m slack, 600 A) hangs between two points 30 m
above flat ground in still air at 300 K. FLARE's fire engine releases a prescribed 2e5 W/m^2 over a
40 m disc under the middle of the span from t = 0 (`erf.fire.prescribed_heat.*`; nothing else
burns). The coupling is passive, so the air stays at rest and at the sounding's temperature, and
the conductor's temperature follows its heat balance alone. Each burning fire cell radiates 0.3 of
its heat release from 2 m above the ground. The conductor starts at its steady temperature without
the fire, heats under it, lengthens and sags.

The same deck runs in three layouts of the fire and the lines on one or two levels. They also
differ in which fire step the lines see first through ERF's fire view (`StructureStep::fire`,
erf-model/ERF#5018).

## Files

- `Conductors_FireHeating.i`: the deck, 1500 x 1000 x 300 m on 30 x 20 x 12 cells, 200 steps of
  0.5 s, fire cells of 10 m (grid ratio 5).
- `input_sounding`: still air, potential temperature 300 K, surface pressure 1000 hPa.
- `check_fire_heating.py`: the check. Pure Python, no numpy, and no code shared with ERF. At every span node
  it recomputes, from the deck's numbers:
  1. the air temperature, T = theta - g z / c_p for the still, isentropic sounding, checked against the
     value the run logs at the hottest node (this catches a slip in the equation of state);
  2. the fire's irradiance at every step, from the fire cells whose centres lie in the disc (within the
     300 m cutoff) and the node positions the run wrote at the step's start (`S1_nodes.dat`);
  3. the temperature, from its steady value without the fire, stepped with backward Euler and the
     IEEE 738-2012 terms written out in the script, radiating to surroundings at the node's first air
     temperature;
  4. from these, the hottest node and its temperature, the span's mean and reference temperatures, the
     thermal strain of every row, and the sag at the end from the bundled MoorDyn stub's elastic parabola.

## Data and acknowledgements

No external data. The conductor's properties are typical of a 795 kcmil 26/7 ACSR conductor; the
resistances are test values, not a manufacturer's data.

## How to run

Build FLARE with `-DERF_ENABLE_STRUCTURES=ON -DERF_MOORDYN_USE_STUB=ON -DERF_ENABLE_FIRE=ON`, FFT and
MPI, then

    mpiexec -n 2 erf_exec Conductors_FireHeating.i
    python3 check_fire_heating.py --dt 0.5 --fire-lags 1 --fire-dx 10

The two-level variants add the options listed in `Structures/Tests/ERFStructuresTests.cmake` and
pass `--dt 0.25 --fire-lags 0 --fire-dx 10` (lines on level 1) or `--dt 0.5 --fire-lags 1 --fire-dx 5`
(fire on level 1).

## Expected results

| Run | Fire / lines | First step sees the fire | Hottest node | Span mean (reference) at 100 s | Irradiance at 100 s | Largest differences: air, E, T_max, T_mean | Sag |
|---|---|---|---|---|---|---|---|
| `Conductors_FireHeating` | level 0 / level 0 | no | 58.11 to 85.83 C | 65.68 C (58.06 C) | 27199 W/m^2 | 2.0e-5 K, 3.9e-10, 2.1e-5 K, 3.1e-5 K | 13.5547 to 13.7306 m |
| `Conductors_FireHeating_LinesOnFinerLevel` | level 0 / level 1 | yes | 58.11 to 85.97 C | 65.72 C (58.06 C) | 27201 W/m^2 | 1.2e-5 K, 4.1e-10, 5.8e-6 K, 5.4e-6 K | 13.5547 to 13.7314 m |
| `Conductors_FireHeating_FireOnFinerLevel` | level 1 (5 m cells) / level 0 | no | 58.11 to 85.82 C | 65.66 C (58.06 C) | 27183 W/m^2 | 2.9e-5 K, 4.0e-10, 3.0e-5 K, 4.0e-5 K | 13.5547 to 13.7301 m |

The sag matches the stub's parabola to the printed digits in every run, and the strain to 1e-11.

Tolerances: air temperature 1e-4 K (ERF's discrete hydrostatic balance differs from the formula by 2e-5 K),
irradiance 1e-6 (the node positions carry ten digits), temperatures 1e-3 K, strain 1e-7, sag 1e-4 m. The run must heat the hottest node by at least 20 K and sag the
span by at least 0.05 m, so that the checks cannot pass on a cold line.

The steady temperature of 58.1 C is the 600 A Joule heating against natural convection and
radiation in still air at 26.7 C (the air at 16 m). Under the fire the hottest node gains
alpha D E = 0.5 x 0.0281 m x 27 kW/m^2 = 380 W/m, about 0.29 K/s at the start.

## CTest

`Conductors_FireHeating`, `Conductors_FireHeating_LinesOnFinerLevel`,
`Conductors_FireHeating_FireOnFinerLevel` (deck and check), `Conductors_FireHeating_Restart` (a
checkpoint at 10 s, with the fire on: the temperature log and `S1.dat` must equal the uninterrupted
run's) and `Conductors_Thermal_FireNeedsFractionAbort` (the deck without `radiative_fraction`), all
with the stub, in `Structures/Tests/ERFStructuresTests.cmake`.
