# FireFrozenWind

## Purpose

Checks `erf.fire.frozen_wind.file`. The fire layer reads the wind once, at
start-up, from a GUST-CFD plotfile or checkpoint and holds it for the whole run.
The fire's own state (front, fuel, moisture) still changes in time. ERF still
advances the atmosphere, but the fire does not read its wind.

It also checks `erf.fire.spotting.drift_wind = frozen`: ember brands drift on
u and v of the frozen 3-D GUST field at their height above the ground, instead
of on the reference wind scaled by a log profile (`profile`, the default).
With `drift_wind = frozen_uvw` they also rise or fall with GUST's vertical wind
w.

The case is a grass fire next to a 40 m Gaussian hill. GUST ran steady RANS over
the hill once (`inputs_gust_hill`), and its plotfile and checkpoint are kept here
(`gust_hill_plt`, `gust_hill_chk`). They can serve as a reference for other
frozen-wind setups.

- Domain: 1280 m x 1280 m x 400 m, periodic in x and y.
- GUST mesh: 32 x 32 x 16 cells. Cells are 40 m wide; the first layer is 5 m,
  stretched to the top.
- ERF mesh: 32 x 32 x 16 cells, on the same terrain file.
- Fire grid: 128 x 128 cells of 10 m.
- GUST wind: geostrophic (10, 0) m/s at latitude 90. Near the ground it is
  about 2.5 to 4.6 m/s, turned towards +y.
- ERF sounding: 5 m/s towards -y. This differs from GUST's wind on purpose, so
  a fire that reads the atmosphere instead of the frozen field shows it.

## Files

| File | Role |
|------|------|
| `make_hill_terrain.py` | Writes `terrain_gaussian_hill.txt`: z = 40 exp(-r^2 / (2 x 150^2)) m, centred at (640, 640) m, nodes every 20 m |
| `terrain_gaussian_hill.txt` | The terrain file, read by both GUST (`terrain.type = erf_grid`) and ERF (`erf.terrain_file_name`) |
| `inputs_gust_hill` | The GUST deck that made `gust_hill_plt` and `gust_hill_chk` |
| `run_gust.sh` | Remakes the terrain file and the two GUST directories (needs the `gust` executable) |
| `gust_hill_plt` | GUST plotfile (u, v, w, p, theta, tke, z, detJ and `Level_0/z_phys_nd`), 1.2 MB |
| `gust_hill_chk` | GUST checkpoint (face velocities, p, theta, tke, eddy viscosity, `z_phys_nd`), 3.6 MB |
| `inputs_base` | What the decks share: domain, terrain, sounding, fire settings |
| `inputs_plt` | Frozen wind from `gust_hill_plt` |
| `inputs_chk` | Frozen wind from `gust_hill_chk`; must give the same reference wind as `inputs_plt` |
| `inputs_nearest_log` | Variant: nearest column, 2 m above the ground (below every first cell centre), log profile below it |
| `inputs_sample_ht` | Variant: checkpoint wind sampled at 20 m and brought to 6.1 m with a log profile (z0 = 0.1 m, GUST's) |
| `inputs_spot3d` | Ember spotting with `drift_wind = frozen`: every front cell launches a brand every step (fixed seed) |
| `inputs_spot_profile` | The same run with `drift_wind = profile` (the default drift) |
| `inputs_spot_nearest_log` | `inputs_spot3d` with `wind_interp = nearest` and `wind_below_first_cell = log`: the brands drift in the nearest GUST column, with the log profile below its lowest cell centre |
| `inputs_spot_uvw` | `inputs_spot3d` with `drift_wind = frozen_uvw`: the brands also move with GUST's w |
| `inputs_empty_abort` | An empty `erf.fire.frozen_wind.file` (abort test) |
| `input_sounding` | Uniform 5 m/s towards -y, theta 300 K |
| `check_frozen_wind.py` | Samples `gust_hill_plt` independently of FLARE and compares the fire plotfiles |
| `run_variants.sh` | Runs the two variant decks and checks each with `check_frozen_wind.py` and its options |
| `check_spot3d.py` | Replays the flight of every brand launched at step 20 and compares its landing distance with `fire_spot_max_dist` |
| `run_spot3d.sh` | Runs the four spotting decks and checks each with `check_spot3d.py` |
| `run_restart.sh`, `compare_restart.py` | Restart from step 10 and compare with a straight run |

## Data and acknowledgements

`gust_hill_plt` and `gust_hill_chk` were made for this case with GUST-CFD
(erf-model/gust-abl, commit 547e376, BSD 2-Clause licence) from
`inputs_gust_hill` and `terrain_gaussian_hill.txt`. Thanks to the GUST-CFD
developers. Both directories are redistributed here in full. No external
dataset is used: the terrain is analytic.

## How to run

```bash
mpirun -np 2 /path/to/flare inputs_plt max_step=20
python3 check_frozen_wind.py
mpirun -np 2 /path/to/flare inputs_spot3d max_step=20
python3 check_spot3d.py --drift frozen
mpirun -np 2 /path/to/flare inputs_spot_profile max_step=20
python3 check_spot3d.py --drift profile
mpirun -np 2 /path/to/flare inputs_spot_nearest_log max_step=20
python3 check_spot3d.py --drift frozen --interp nearest --below log
mpirun -np 2 /path/to/flare inputs_spot_uvw max_step=20
python3 check_spot3d.py --drift frozen_uvw
```

`run_spot3d.sh /path/to/flare` runs the four spotting pairs.

`check_frozen_wind.py` needs `erf_plotfile.py` (from `Submodules/ERF/Exec/Tools`)
next to it. CTest copies it there.

To remake the GUST files (about 1.5 minutes on 4 ranks):

```bash
MPIRUN="mpirun -np 4" ./run_gust.sh /path/to/gust
```

## Expected results

The check samples the GUST plotfile in Python, the way the inputs docs describe:

- the four GUST columns around each fire cell centre, with bilinear weights;
- each column's ground is the mean of its four k = 0 node heights;
- each column is sampled 6.1 m above its own ground, linearly between the two
  cell centres that bracket that height.

The variants follow the same rules with `--interp nearest` (the column that
holds the fire cell centre), `--below log` (below the first centre the wind is
scaled by log(z/z0)/log(z1/z0), z1 the first centre above the ground) and
`--zsample 20` (the wind at 20 m times log(6.1/z0)/log(20/z0)).

| Check | Result on the Mac (Debug, 2 ranks) |
|-------|-----------------------------------|
| `fire_wind_ref` vs the Python sampling, plotfile deck | 1.8e-15 m/s (tol 1e-9) |
| `fire_wind_ref` vs the Python sampling, checkpoint deck | 1.8e-15 m/s (tol 1e-9) |
| `fire_extract_z` vs ground + 6.1 m | 1.4e-14 m (tol 1e-9) |
| Mean v of the reference wind | +1.91 m/s (the sounding has -5) |
| Wind at step 20 vs step 10 | 0 (frozen) |
| Restart from step 10 vs straight run, step 20 | 0 for the wind, height and level set |
| Nearest column, 2 m, log below the first centre | 4.0e-15 m/s, 7.1e-15 m |
| Checkpoint sampled at 20 m, brought to 6.1 m | 1.8e-15 m/s, 1.4e-14 m |

The deck turns off `use_waf` and `use_terrain_wind`. With them on, they would
still act on the frozen wind, and `use_terrain_wind` would count the hill a
second time on top of GUST's flow.

The start-up line also prints the largest difference between the ground of the
GUST column under each fire cell and ERF's column ground there (0 m here: both
read the same terrain file on the same mesh). This is information only: the
fire keeps ERF's terrain for slope, and the reference wind is sampled above
GUST's ground, so `fire_extract_z` is GUST's ground + height.

The check was shown to fail on four mutants:

- sampling each column at 6.1 m above z = 0 instead of above its ground
  (height off by 38.6 m);
- running `inputs_base` without the frozen file (wind off by 8.3 m/s, mean
  v = -4.96 m/s, and the wind changes between steps 10 and 20);
- no log factor below the first centre and no 20 m to 6.1 m factor (wind off
  by 0.20 m/s and 1.27 m/s in the two variants);
- the nearest column shifted by half a GUST cell (wind off by 0.46 m/s,
  height by 6.2 m).

### Ember spotting

`check_spot3d.py` replays, in Python, the flight of the brand each source cell
launched at step 20:

- the brand starts `fire_lofting_height` above the ground at the source cell
  centre and falls at the terminal velocity (2.5 m/s) in 20 forward-Euler
  steps, sized from the fall down to the lowest ground of the domain;
- `frozen`: each step it moves with u, v of the GUST field at its height h
  above the ground. Each of the four GUST columns around it is sampled at its
  own ground + h, linearly between cell centres, and the columns are weighted
  bilinearly (`inputs_spot_nearest_log`: the column holding the brand, and
  below its lowest centre that centre's wind times
  ln(h'/z0)/ln(h1/z0), h' and h1 the heights above the column's ground);
- `frozen_uvw` (`inputs_spot_uvw`): as `frozen`, and it also moves with
  GUST's w, sampled as u and v, so it falls (2.5 m/s - w) per second. It keeps
  the step size, steps on until it lands (up to 200 steps) and is dropped if
  still aloft then;
- `profile`: it moves with `fire_wind_ref` times ln(h/z0)/ln(6.1/z0), with z0 =
  0.1 m;
- it lands at or below the ground under it, and the distance is capped at
  200 m (fuel model 1). The terminal velocity of 2.5 m/s keeps every landing
  short of that cap, so the drift sets every distance.

| Check | `inputs_spot3d` (frozen) | `inputs_spot_profile` (profile) | `inputs_spot_nearest_log` (frozen) | `inputs_spot_uvw` (frozen_uvw) |
|-------|--------------------------|---------------------------------|------------------------------------|--------------------------------|
| Source cells at step 20 | 726 | 752 | 719 | 703 |
| Longest landing | 151.1 m | 132.5 m | 144.6 m | 125.8 m |
| `fire_spot_max_dist` vs the replay | 2.6e-13 m (tol 1e-6) | 2.1e-13 m (tol 1e-6) | 2.1e-13 m (tol 1e-6) | 2.0e-13 m (tol 1e-6) |
| The other drift's replay differs by up to | 25.3 m | 25.3 m | 23.2 m | 25.8 m (frozen) |
| The bilinear, clamped replay differs by up to | | | 4.5 m (671 samples below a lowest centre) | |
| Flights with more / fewer steps than the frozen replay | | | | 71 / 252 (longest 20 steps, none aloft) |

The brands fly to about 47 m, below the highest GUST cell centres, so no
brand is counted above them here. The unit test
`FrozenDrift.AboveTheTopTheBrandTakesTheTopWind` covers that case: it checks
the landing and that every brand is counted. The unit tests also pin the shear
closed form and the gather of the GUST columns under the fire region
(`FrozenDrift.GatherCopiesTheColumnsUnderTheFireRegion`, on 8 x 8 column
boxes).

Each check fails when it is run on the other deck's output (landing off by
25.3 m), so it detects a drift that ignores `drift_wind`. The
`inputs_spot_nearest_log` replay fails by 3.3 m when FLARE drifts bilinearly
there, and by 4.5 m when it drops the log profile, so it detects a drift that
ignores `wind_interp` or `wind_below_first_cell`.

In `inputs_spot_uvw` GUST's w lies between -0.54 and +0.63 m/s, a quarter of
the 2.5 m/s terminal velocity at most, so w shortens or lengthens flights by
a few steps and no brand needs more than N = 20 steps or stays aloft. A drift
that ignores w fails the replay by 25.8 m. The cap of 10 N steps and the drop
of brands still aloft are pinned by the unit tests
(`FrozenDrift.VerticalWindChangesTheFallClosedForm`,
`FrozenDrift.BrandsStillAloftAreDroppedAndCounted`), and
`FrozenDrift.VerticalWindOfZeroLandsAsFrozen` checks that w = 0 lands every
brand where `frozen` does, bit for bit, over 400 sloping grounds. The
`inputs_spot_uvw` run with `erf.fire.frozen_wind.file = gust_hill_chk` (w
averaged from the checkpoint's z faces) passes the same replay, which reads
the plotfile's w, to 2.0e-13 m.

## CTest

| Test | What it checks |
|------|----------------|
| `FireFrozenWind_plt` | Reference wind and sampling height from the plotfile against the independent sampling; frozen in time (2 ranks) |
| `FireFrozenWind_chk` | The same from the checkpoint (2 ranks) |
| `FireFrozenWind_variants` | Nearest column with the log profile below the first centre, and the 20 m resampling from the checkpoint, against the independent sampling (2 ranks) |
| `FireFrozenWind_spot3d` | Ember landings with `drift_wind = frozen`, `frozen_uvw` and `profile`, and frozen with the nearest column and the log profile, against the replayed flights (2 ranks) |
| `FireFrozenWind_restart` | The restart re-reads the file; step 20 equals the straight run bit for bit (2 ranks) |
| `FireFrozenWind_prescribed_abort` | `prescribed_wind = true` with a frozen file stops the run |
| `FireFrozenWind_no_header_abort` | A path that is not a GUST directory stops the run |
| `FireFrozenWind_coverage_abort` | A fire region wider than GUST's domain (1300 m vs 1280 m) stops the run |
| `FireFrozenWind_above_top_abort` | `wind_ref_ht = 340` m stops the run: above the highest cell centre of the hill columns (330.6 m above their ground), below that of the flat columns (365.9 m) |
| `FireFrozenWind_empty_abort` | An empty file name stops the run (leave the key unset to turn it off) |
| `FireFrozenWind_open_columns_abort` | `structures.wind_open_columns` with a frozen file stops the run |
| `FireFrozenWind_spot_no_file_abort` | `spotting.drift_wind = frozen` without a frozen file stops the run |
| `FireFrozenWind_spot_uvw_no_file_abort` | `spotting.drift_wind = frozen_uvw` without a frozen file stops the run |
| `FireFrozenWind_spot_drift_value_abort` | A `spotting.drift_wind` value other than `profile`, `frozen` or `frozen_uvw` stops the run |

`FireFrozenWind_plt`, `_chk`, `_variants` and `_spot3d` run in double-precision builds
only: they compare to 1e-9, and `gust_hill_chk` holds 8-byte reals, which a
single-precision build refuses (use the plotfile there).
