# Fire under a transmission line

## Purpose

A benchmark of FLARE's fire coupled to a transmission line: a grass fire runs under three spans of conductor
on two lattice towers, and the conductors and the tower steel heat up from the fire's radiation and from the
air the fire warms. The same case runs twice, once with the fire's heat put into the air
(`erf.fire.fire_atm_feedback = 1`) and once without it (`= 0`), so the comparison shows what the two-way
coupling changes for the line.

It is a demonstration on a simple configuration, not a validation: no measured data of a fire under a line are
used. Each model it joins is checked on its own elsewhere (see the docs page
`Docs/sphinx_doc/FireUnderLine.rst`):

- the fire: FLARE's level-set spread with Rothermel's rate on Anderson fuel model 3 (tall grass), with the wind
  brought to mid-flame height by Andrews' wind adjustment factor, and lagged coupling to the atmosphere
  (`Docs/sphinx_doc/theory/fire_coupling.rst`);
- the flow: ERF's k-equation RANS closure (`Docs/sphinx_doc/theory/RANS.rst`) over flat ground, anelastic with
  MidPoint time stepping and the FFT pressure solve;
- the conductors: MoorDyn-C lines on insulator strings, their temperature from the heat balance of IEEE Std
  738-2023 with the fire's radiation and the local air, and the thermal sag that follows ("Conductor temperature
  under a fire" in `Docs/sphinx_doc/theory/Conductors.rst`; unit test
  `LineThermal.TheTermsMatchIEEE7382023sSampleCalculations`);
- the towers: lattice frames whose members heat as unprotected steel (EN 1993-1-2 lumped heating), the frame
  rebuilt as the steel softens ("Tower steel temperature under a fire", same page; unit test
  `TowerThermal.TheJRCWorkedExampleReachesItsCriticalTemperature`).

## Files

- `make_inflow.py` writes `inflow_profile` and `input_sounding`: a neutral log law, 6 m/s at 10 m over
  z0 = 0.05 m (u* = 0.464 m/s), under a capping inversion at 400 m, with a row at every cell centre so that ERF's
  linear interpolation gives each cell the log law's value. It prints the inflow turbulent kinetic energy
  `xlo.KE`, the k-equation closure's neutral wall value u*²/(c_μ⁰)². `python3 make_inflow.py` gives the
  committed files.
- `flow.inputs` holds everything the stages share: the flow, the fire, the line.
- `ignition_line.csv`: the 600 m ignition line at x = 600 m, across the wind, centred on the middle span.
- `inputs_spinup`: stage 1, the flow and the line without the fire, 590 s, to `chk01180`.
- `inputs_feedback_on`, `inputs_feedback_off`: stage 2, each restarting from `chk01180`; the fire is lit at
  600 s and followed to 2100 s (25 min).
- `compare_runs.py`: reads the two runs' logs and writes `summary.csv` and `compare_runs.png`.
- `make_figure.py`: snapshots of the two runs when their middle span is hottest (plan view of the fire, the
  air through the middle span); needs yt.
- `check_feedback.py`: the CTest's check (below).

The domain is 4000 m along the wind, 2000 m across it (periodic) and 600 m high, on 25 m × 25 m × 10 m cells
(768 000 cells), with 5 m fire cells and a fixed step of 0.5 s. The 8 × 4 boxes of 20 × 20 columns each divide by
the fire grid ratio; up to 32 ranks get whole boxes. The line runs north-south at x = 1200 m: dead ends at
y = 550 m and 1450 m, two 30 m lattice towers (310 members each) at y = 850 m and 1150 m, three 300 m spans of
Drake ACSR at 600 A strung to 20 kN.

Where the line's values come from:

| value | source |
|---|---|
| diameter 28.1 mm, resistance at 25 °C and 75 °C, heat capacity 1310 J/(m K) | IEEE Std 738-2023, 4.6 (Drake) |
| absorptivity and emissivity 0.5 | chosen (a weathered conductor; the standard's sample uses 0.8) |
| mass 1.628 kg/m, axial stiffness 3.0e7 N, thermal expansion 1.93e-5 /K | the repository's conductor test values (`Structures/Tests`) |
| radiative fraction 0.3 | chosen share of the fire's heat release that radiates |
| lattice geometry and member sections | the repository's generated lattice tower (`Docs/sphinx_doc/theory/TowerFrame.rst`) |

## Data and acknowledgements

No external data: the inflow is the synthetic log law of `make_inflow.py`, and the fuel is a standard fuel model
built into FLARE. Sources of the methods, each one we have access to:

- Anderson, H. E. (1982). *Aids to determining fuel models for estimating fire behavior.* USDA Forest Service
  General Technical Report INT-122. https://www.fs.usda.gov/rm/pubs_int/int_gtr122.pdf (public domain).
- Rothermel, R. C. (1972). *A mathematical model for predicting fire spread in wildland fuels.* USDA Forest
  Service Research Paper INT-115. https://www.fs.usda.gov/rm/pubs_int/int_rp115.pdf
- Andrews, P. L. (2012). *Modeling wind adjustment factor and midflame wind speed for Rothermel's surface fire
  spread model.* USDA Forest Service General Technical Report RMRS-GTR-266.
  https://www.fs.usda.gov/rm/pubs/rmrs_gtr266.pdf
- IEEE Std 738-2023, *IEEE Standard for Calculating the Current-Temperature Relationship of Bare Overhead
  Conductors*. doi:10.1109/IEEESTD.2023.10382442. Used with institutional access; cited, not reproduced.
- Vassart, O., Zhao, B., Cajot, L. G., Robert, F., Meyer, U. and Frangi, A. (2014). *Eurocodes: Background &
  Applications. Structural Fire Design. Worked examples.* JRC Science and Policy Report EUR 26698 EN.
  doi:10.2788/85432. Open access, reproduction authorised with acknowledgement; we thank its authors and the
  European Commission's Joint Research Centre.
- Axell and Liungman (2001), *Environmental Fluid Mechanics* 1, 71-106: the k-equation closure and its wall
  value, as `Docs/sphinx_doc/theory/RANS.rst` cites it.

We thank the USDA Forest Service for keeping the fire reports freely available.

## How to run

Build with `-DERF_ENABLE_STRUCTURES=ON` against MoorDyn-C (`Structures/setup_moordyn.sh`), `-DERF_ENABLE_FIRE=ON`,
`-DERF_ENABLE_FFT=ON` and MPI. Then, in an empty run directory, with `CASE` this folder and `EXE` the build's
`Exec/erf_exec`:

    mkdir spin on off
    for d in spin on off; do cp $CASE/flow.inputs $CASE/inflow_profile $CASE/input_sounding $CASE/ignition_line.csv $d/; done
    cp $CASE/inputs_spinup spin/; cp $CASE/inputs_feedback_on on/; cp $CASE/inputs_feedback_off off/
    (cd spin && mpiexec -n 4 $EXE inputs_spinup)
    for d in on off; do
      cp -r spin/chk01180 spin/conductors spin/L1_* $d/
      (cd $d && mpiexec -n 4 $EXE inputs_feedback_$d)
    done
    python3 $CASE/compare_runs.py on off
    python3 $CASE/make_figure.py on off

Copying the spin-up's `conductors/` and `L1_*` logs keeps those logs continuous from t = 0: a restart trims
them to the checkpoint's time and appends. `fire_stats.csv` starts at the restart, since there is no fire before
it. With 4 ranks on an Apple-silicon Mac the spin-up takes about 8 minutes and each stage-2 run about 20
minutes.

## Expected results

Built against MoorDyn-C 2.7.1, 4 ranks (`compare_runs.py on off`, rounded):

| quantity | feedback on | feedback off |
|---|---|---|
| mean head rate of spread, first 10 min after ignition | **1.58 m/s** | 0.65 m/s |
| burned area 25 min after ignition | 256 ha | 122 ha |
| middle span: conductor before the fire | 31.5 °C | 31.5 °C |
| middle span: hottest conductor node | **77.5 °C**, 7.1 min after ignition | 51.1 °C, 17.3 min after ignition |
| middle span: largest fire irradiance on the conductor | 48.8 kW/m² | 24.7 kW/m² |
| middle span: largest added sag | 0.93 m | 0.56 m |
| middle span: smallest ground clearance (18.50 m before) | 17.57 m | 17.95 m |
| outer spans: hottest conductor node | 70.0 °C | 51.8 °C |
| towers: hottest member steel | **116 °C** | 98 °C |
| towers: smallest stiffness factor k_E; frame rebuilds | 0.990; 2 | 1.000; 0 |
| towers: largest member utilisation | 0.072 | 0.060 |

With feedback on, the fire's heat changes the wind it spreads in: the head runs 2.5 times faster, and the line is
reached about 6 minutes after ignition instead of about 16. The plume carries hot air downwind of the front, so the
air at the conductor reaches 60 °C to 69 °C from about 3.7 min, before the front arrives; the conductor is at about
56 °C when the fire's radiation peaks (41 kW/m² at its hottest node at 6.3 min) and reaches 77.5 °C. With feedback
off, the air at the line stays at 26.7 °C and only the radiation heats the conductor, from a slower, less intense
front. The towers' steel follows the same pattern; with feedback on the hottest members pass the 1 % change in the
stiffness factor that rebuilds the frame twice. All members stay far below their capacity (utilisation under 0.08).

After the 590 s spin-up the wind at the line is 4.7 m/s at 5 m and 5.9 m/s at 15 m, about 10 % below the inflow's
log law: the RANS flow adjusts over the 1200 m from the inflow face. The fire stays inside the domain in both runs
(x from 465 m to 2985 m, y from 360 m to 1640 m at 25 min). FLARE's start-up warning that the fire "can reach the
west edge" uses the head's rate for the backing fire, which spreads far more slowly; it reaches x = 465 m.

These numbers come from MoorDyn-C. The CTest uses the MoorDyn stub, which hangs each span as a parabola instead of
solving the line, so its sag differs.

## CTest

`FireUnderLine_Feedback` (with the MoorDyn stub and the fire on) runs this directory's decks for 80 s with a
450 m line fire lit at 5 s, 50 m upwind of the line, once with feedback off and once with it on.
`check_feedback.py` requires the air at the line to warm by more than 1 K with feedback on, and the air at every
span to move by less than 1 % of that rise with feedback off while the fire's radiation reaches the line. The
logged air is the air at each span's hottest node; that node moves along the span (up to 10 m in height) as the
fire heats one end more than the middle, so with feedback off the logged air still moves by a few hundredths of a
kelvin.

With the decks here the run gives a rise of 34 K with feedback on and at most 0.08 K with feedback off (limit
0.34 K). A mutant that injects the fire's heat whatever the multiplier (the `fire_atm_feedback > 0` gate and the
multiplier dropped in `Source/FireCoupling/ERF_FireCoupler.cpp`) fails it: the air with feedback off then warms by
31 K to 34 K.
