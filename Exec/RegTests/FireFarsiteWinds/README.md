# FireFarsiteWinds

## Purpose

Checks `erf.fire.farsite_winds`, which writes the fire layer's reference wind as
gridded winds for FARSITE and FlamMap. The layout is the one WindNinja writes:
`_vel.asc` and `_ang.asc` ESRI ASCII grids, plus a `.atm` file that lists them by
date and time. The format is described in the Fire output docs
(`Docs/sphinx_doc/theory/fire_output.rst`, section "FARSITE and FlamMap gridded
winds").

The domain is 600 m x 200 m with 20 m cells, split into three boxes. The fire
grid is 60 x 20 cells of 10 m. Both regression tests run on two ranks, so each
rank owns a different number of boxes and the writer's gather onto the I/O rank
is tested. Grids are written every 60 s, at the steps that start at t = 0, 60 and
120 s (steps 1, 61 and 121).

## Files

| File | Role |
|------|------|
| `inputs_base` | What the decks share: domain, sounding, fire settings |
| `inputs_prescribed` | Prescribed u = 3, v = 4 m/s; ENGLISH; map corner and `.prj` |
| `inputs_atmosphere` | Wind interpolated from the ERF columns at 10 m; METRIC; `coarsen = 2`; fire plotfile at step 61 |
| `inputs_no_start` | Abort deck: no `start_datetime` |
| `inputs_key_only` | Abort deck: a `farsite_winds` key without `interval` |
| `input_sounding` | Uniform u = 5, v = 2 m/s, theta = 300 K |
| `utm13n.prj` | A UTM zone 13N projection, copied beside each grid |
| `check_farsite_winds.py` | The check script (standard library and `erf_plotfile.py`) |

## Data and acknowledgements

No external data is used. The output layout follows WindNinja (USDA Forest
Service, Rocky Mountain Research Station, Missoula Fire Sciences Laboratory;
https://github.com/firelab/windninja). We thank the WindNinja authors for
publishing the source. No WindNinja code is included here: the RMRS code is in
the public domain, and outside contributions are BSD-licensed.

## How to run

```bash
mpirun -np 2 erf_exec inputs_prescribed
mpirun -np 2 erf_exec inputs_atmosphere
python3 check_farsite_winds.py      # needs erf_plotfile.py from Exec/Tools beside it
```

To convert the output for viewing:

```bash
python3 ../../Tools/farsite_winds_convert.py winds_prescribed/prescribed.atm --to vtk shp kmz pdf --stride 4
```

## Expected results

| Check | Prescribed (ENGLISH) | Atmosphere (METRIC) |
|-------|----------------------|---------------------|
| `.atm` | `WINDS`, `ENGLISH`, lines `12 30 1100`, `1101`, `1102` | `WINDS`, `METRIC`, same times |
| Grid | 60 x 20 cells of 10 m, corner 482000, 4423000 | 30 x 10 cells of 20 m, corner 0, 0 |
| Speed | every cell 11.18 mph (5 m/s) | 2 x 2 means of `fire_wind_ref` to 0.006 km/h |
| Direction | every cell 217 (180 + atan2(3, 4) = 216.87) | the same means to 0.51 deg |
| Other | each grid has a `.prj` equal to `utm13n.prj` | mean from 248.2 deg (the sounding's, within 10 deg), 17.2 km/h (70-100 % of the sounding's) |

The atmospheric mean speed (17.2 km/h) is below the sounding's 19.4 km/h
because the surface layer slows the wind at 10 m.

## CTest

| Test | What it checks |
|------|----------------|
| `FireFarsiteWinds_prescribed` | Closed-form values, `.atm`, headers, `.prj` (2 ranks) |
| `FireFarsiteWinds_atmosphere` | Grids against the fire plotfile, coarsening (2 ranks) |
| `FireFarsiteWinds_height_abort` | `units = metric` with a 6.096 m reference wind stops the run |
| `FireFarsiteWinds_no_start_abort` | Missing `start_datetime` stops the run |
| `FireFarsiteWinds_interval_abort` | `interval = 30` stops the run (the `.atm` keeps minutes) |
| `FireFarsiteWinds_key_without_interval_abort` | A `farsite_winds` key without `interval` stops the run |
| `FireFarsiteWinds_coarsen_abort` | `coarsen = 7` does not divide 60 x 20 and stops the run |
| `FarsiteWindsConvert_SelfTest` | The converter round-trips a synthetic grid |

The unit tests `FarsiteWinds.*` (in `Tests/Unit/Fire/ERF_GTestFarsiteWindOutput.cpp`)
check the direction convention, the unit factors, the time stamps, north-first
row order on a field that varies in both directions, the cadence and the
restart append.
