# FireWindColumnZ0

## Purpose

Checks `erf.fire.wind_sample_z0_source`, which sets where the roughness of the
fire layer's neutral log profiles comes from:

- `scalar` (default): `erf.fire.wind_sample_z0` everywhere, as before.
- `surface_layer`: the roughness of each atmospheric column in ERF's surface
  layer (`erf.most.z0`, or the map of `erf.most.roughness_file_name`).

The roughness is used in four places. The sampler brings the wind sampled at
`erf.fire.wind_sample_ht` (or the wind of the first cell centre) down to the
reference height. The 10 m wind of FBP, Cheney-Gould and Cruz is computed from
the reference wind. The wind is brought to 20 ft before the wind adjustment
factor (WAF). The ember drift lifts the reference wind to the brand's height.
The test checks the first three; the unit test
`EmberDrift.PerCellRoughnessReplacesTheScalar` checks the drift.

The test also checks restart parity with `bulk_coeff` fluxes, whose z0 the
surface layer recovers every step.

The domain is 320 m x 160 m x 160 m with 5 m cells, in four boxes on two ranks.
The fire grid has four fire cells per atmospheric column. The wind is sampled
at 20 m with `erf.fire.wind_interp = "nearest"`, so each fire cell takes its
own column's wind and roughness. The coupling is one-way, so both roughness
sources see the same atmosphere. The ratio of the two runs' fire fields at a
cell is then a closed form of that cell's roughness.

## Files

| File | Role |
|------|------|
| `inputs_base` | Grass fire on flat ground; uniform `erf.most.z0 = 0.1`, the scalar's value |
| `inputs_map` | `inputs_base` with `erf.most.roughness_file_name = roughness_halves.txt` |
| `input_sounding` | Uniform u = 8, v = 4 m/s, theta = 300 K |
| `run_wind_column_z0.sh` | Writes the maps, runs the fourteen runs, runs the checks |
| `check_firewindcolumnz0.py` | Writes the maps (`write_map`) and checks the plotfiles (standard library and `erf_plotfile.py`) |

`roughness_halves.txt` and `fuel_stripes.asc` are written by the script, not
committed. The first holds the (x, y, z0) node triples of the 64 x 32 grid:
z0 = 0.03 m for columns i < 32 (x < 160 m) and 0.5 m for i >= 32. The second
is a fire-grid fuel map of Anderson 1 and 3 in stripes 16 fire cells wide.

## Data and acknowledgements

No external data. The Cheney, Gould and Catchpole (1998) grassland rate used
in the 10 m check is the formula already in `FireEngine/ERF_CheneyGouldModel.H`.

## How to run

```
MPIRUN="mpirun -np 2" ./run_wind_column_z0.sh /path/to/flare
SKIP_RUN=1 ./run_wind_column_z0.sh x     # only rerun the checks
```

Each run takes 4 steps and writes the step-4 fire plotfile. Each variant is
run with `scalar` (`_sc`) and with `surface_layer` (`_sl`):

| Variant | Deck | Extra inputs |
|---------|------|--------------|
| `uni` | `inputs_base` | none |
| `map` | `inputs_map` | none |
| `map10` | `inputs_map` | `erf.fire.wind_ref_ht = 10` |
| `cg` | `inputs_map` | `erf.fire.ros_model = cheney_gould` |
| `fmap` | `inputs_map` | as `map10`, with `erf.fire.fuel_map.file = fuel_stripes.asc` |
| `bil` | `inputs_map` | `wind_interp = bilinear`, `wind_sample_ht = 0`, `wind_below_first_cell = log`, `wind_ref_ht = 2` |

The restart pair runs `inputs_base` with `surface_layer` and `bulk_coeff`
fluxes: `rs_a` straight to step 4 with a checkpoint at step 2, and `rs_b` from
that checkpoint to step 4.

## Expected results

With F(z_to, z_from, z0) = ln(z_to / z0) / ln(z_from / z0), z0c the cell's
column roughness and 0.1 m the scalar:

| Check | Quantity, surface_layer over scalar | Expected | Measured max. rel. error |
|-------|-------------------------------------|----------|--------------------------|
| parity | the six fields, `uni` | identical | 0 values differ |
| sample | `fire_wind_ref` u and v, `map` | F(6.1, 20, z0c) / F(6.1, 20, 0.1) | 2.2e-16 |
| ft20 | `fire_wind_ref` u, `map10` | F(10, 20, z0c) / F(10, 20, 0.1) | 2.2e-16 |
| ft20 | `fire_wind_eff / fire_wind_ref`, `map10` | F(6.096, 10, z0c) / F(6.096, 10, 0.1) | 4.4e-16 |
| m10 | `fire_ros`, `cg` | R(U10 with z0c) / R(U10 with 0.1) | 4.4e-16 |
| fmap | `fire_wind_eff / fire_wind_ref`, `fmap` | F(6.096, 10, z0c) / F(6.096, 10, 0.1) | 4.4e-16 |
| bil | `fire_wind_ref` u and v, `bil` | F(2, 2.5, z0c) / F(2, 2.5, 0.1) | 4.4e-16 |
| restart | the six fields, `rs_b` over `rs_a` | identical | 0 values differ |

In bil the 2 m reference wind is below the 2.5 m first cell centre, so it is
the centre's wind times F(2, 2.5, z0). The check skips the fire cells in
atmospheric columns 0, 31, 32 and 63, whose bilinear stencil reaches a column
of the other roughness (the x = 160 m step and the periodic seam).

For m10, U10 = |`fire_wind_ref`| F(10, 6.1, z0), taken from each run's own
reference wind. Above 5 km/h and below 12 % moisture the CG98 rate is
R = (1.4 + 0.838 (3.6 U10 - 5)^0.844) phi_M phi_C, and phi_M phi_C cancels in
the ratio. With `erf.fire.moisture_dynamic = false` the rate reads
`erf.fire.cheney_gould.moisture` (default 10 %), not `moisture_1hr`. The lowest
U10 in the test is 26 km/h.

The ratios are about 1.05 over the western half and 0.87 over the eastern half
for the sampled wind. The tolerances are 1e-12 (1e-10 for the rate).

Each check fails against a one-line change to the code (measured 2026-11-23):

| Change | Check that fails |
|--------|------------------|
| sampler given no column roughness | sample (error 0.14) |
| 20 ft conversion with the scalar | ft20 eff/ref (error 0.069) |
| 10 m wind with the scalar | m10 (error 0.060) |
| column roughness scaled by 1.001 | parity (196608 values differ) |
| 20 ft conversion with the scalar on a fuel map | fmap eff/ref (error 0.069) |
| roughness copied only when the sea roughness varies | restart (196608 values differ) |

## CTest

| Test | What it runs |
|------|--------------|
| `FireWindColumnZ0` | `run_wind_column_z0.sh` on two ranks (MPI builds, not Windows) |
| `FireWindColumnZ0_value_abort` | `wind_sample_z0_source = column`: not one of the two values |
| `FireWindColumnZ0_frozen_abort` | `surface_layer` with `erf.fire.frozen_wind.file` |
| `FireWindColumnZ0_prescribed_abort` | `surface_layer` with `erf.fire.prescribed_wind = true` |
| `FireWindColumnZ0_z0_range_abort` | `surface_layer` with `erf.most.z0 = 2.4`, above 6.096 m / e |
| `FireWindColumnZ0_no_effect_abort` | `surface_layer` where no log profile is applied (no sampling height, a 10 m reference wind, no WAF) |
