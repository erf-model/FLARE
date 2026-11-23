# FireFrozenWind

## Purpose

Checks `erf.fire.frozen_wind.file`. The fire layer reads the wind once, at
start-up, from a GUST-CFD plotfile or checkpoint and holds it for the whole run.
The fire's own state (front, fuel, moisture) still changes in time. ERF still
advances the atmosphere, but the fire does not read its wind.

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
| `inputs_empty_abort` | An empty `erf.fire.frozen_wind.file` (abort test) |
| `input_sounding` | Uniform 5 m/s towards -y, theta 300 K |
| `check_frozen_wind.py` | Samples `gust_hill_plt` independently of FLARE and compares the fire plotfiles |
| `run_variants.sh` | Runs the two variant decks and checks each with `check_frozen_wind.py` and its options |
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
```

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

The check was shown to fail on two mutants:

- sampling each column at 6.1 m above z = 0 instead of above its ground
  (height off by 38.6 m);
- running `inputs_base` without the frozen file (wind off by 8.3 m/s, mean
  v = -4.96 m/s, and the wind changes between steps 10 and 20);
- no log factor below the first centre and no 20 m to 6.1 m factor (wind off
  by 0.20 m/s and 1.27 m/s in the two variants);
- the nearest column shifted by half a GUST cell (wind off by 0.46 m/s,
  height by 6.2 m).

## CTest

| Test | What it checks |
|------|----------------|
| `FireFrozenWind_plt` | Reference wind and sampling height from the plotfile against the independent sampling; frozen in time (2 ranks) |
| `FireFrozenWind_chk` | The same from the checkpoint (2 ranks) |
| `FireFrozenWind_variants` | Nearest column with the log profile below the first centre, and the 20 m resampling from the checkpoint, against the independent sampling (2 ranks) |
| `FireFrozenWind_restart` | The restart re-reads the file; step 20 equals the straight run bit for bit (2 ranks) |
| `FireFrozenWind_prescribed_abort` | `prescribed_wind = true` with a frozen file stops the run |
| `FireFrozenWind_no_header_abort` | A path that is not a GUST directory stops the run |
| `FireFrozenWind_coverage_abort` | A fire region wider than GUST's domain (1300 m vs 1280 m) stops the run |
| `FireFrozenWind_above_top_abort` | `wind_ref_ht = 340` m stops the run: above the highest cell centre of the hill columns (330.6 m above their ground), below that of the flat columns (365.9 m) |
| `FireFrozenWind_empty_abort` | An empty file name stops the run (leave the key unset to turn it off) |
| `FireFrozenWind_open_columns_abort` | `structures.wind_open_columns` with a frozen file stops the run |

`FireFrozenWind_plt`, `_chk` and `_variants` run in double-precision builds
only: they compare to 1e-9, and `gust_hill_chk` holds 8-byte reals, which a
single-precision build refuses (use the plotfile there).
