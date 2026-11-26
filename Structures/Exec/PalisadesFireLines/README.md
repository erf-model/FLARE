# Palisades Fire and two transmission lines

## Purpose

A demonstration of the full FLARE chain on a real fire: the Palisades Fire of 7 January 2025 in the Santa
Monica Mountains, Los Angeles, run from public data with two real transmission lines in the domain. It is
the companion of `Structures/Exec/MarshallFireLines` on steep chaparral instead of grass:

- **Weather:** one HRRR analysis hour (19 UTC), the column at the ignition point, as the sounding and the
  inflow profile (`Exec/Tools/hrrr_profile.py`). The wind is from the north-northeast.
- **Fuel:** the LANDFIRE LF2024 FBFM40 map on the fire grid (`Exec/Tools/make_landfire_fuel_map.py`), with
  the multi-class BEHAVE spread model on each fuel cell (`erf.fire.ros_model = "behave"`).
- **Fuel moisture:** from the Los Angeles County Fire Department Daily Fire Danger Analysis of 7 January 2025
  (see the deck comments and the table below).
- **Lines:** two OpenStreetMap `power=line` ways, L1 (voltage not tagged in OpenStreetMap) and the 230 kV
  Tarzana - Olympic line with the electrode line of the Pacific DC Intertie (L2), with towers picked along them
  over the terrain (`Exec/Tools/osm_power_lines.py`, `Exec/Tools/fire_lines_towers.py`).
- **Exposure:** the conductor temperature of every span (IEEE Std 738-2023 heat balance with the fire's
  radiation and the local air) and the steel temperature of L1's lattice towers (EN 1993-1-2 lumped heating).

It is a demonstration, not a validation. The head speed is compared with a satellite estimate; the burned area
and the perimeter are not compared. The atmosphere is dry: the HRRR moisture is
set to zero because the run has no moisture model. In this run the fire spreads downwind toward the coast,
away from both lines, so the lines stay near the air temperature; the case shows the chain working, not a
line exposure.

The setup:

| Item | Value |
|---|---|
| Domain | 24 km square, 80 x 80 x 102 cells, 10 m first cell stretched by 1.03 to a 6.5 km top |
| Terrain | `terrain_palisades.txt` and `terrain_palisades_fire.txt`, made by `Exec/CanonicalTests/Fire/Real_Terrain/Palisades_Fire/gen_palisades.py` |
| Atmosphere | anelastic, FFT pressure solve, YSU boundary layer with 2D Smagorinsky, Rayleigh damping in the top 1.5 km |
| Inflow | north and east faces, HRRR 2025-01-07 19 UTC (pressure levels from 20 m above the ground up, the 10 m wind below them, a log law with z0 = 0.1 m to 2 m); outflow on the west and south faces |
| Fire | level set on a grid 10 times finer (30 m cells) on the fire terrain, BEHAVE multi-class spread per fuel cell, Andrews WAF, lagged two-way coupling, spotting |
| Ignition | one fire of 100 m radius at (9000, 8000) m, lit at 600 s after a 10 min wind spin-up |
| Run length | 6900 s: 10 min spin-up, then 105 min of fire |
| L1 | 18 spans of Drake ACSR, 92-511 m long, lattice towers 30-33 m high with steel temperature |
| L2 | 17 spans of Drake ACSR, 208-541 m long, rigid point towers 30-48 m high (wire temperature only) |

Fuel moisture:

| Class | Value | Source |
|---|---|---|
| 1 h dead | 5 % | NFDRS equilibrium moisture (Cohen and Deeming 1985) at the Malibu Hills RAWS forecast of 22 % RH and 66 F (LACoFD, 7 Jan 2025): 4.8 %, rounded |
| 10 h dead | 8 % | LACoFD, 7 Jan 2025, Malibu Hills RAWS 10 h dead fuel moisture |
| 100 h dead | 10 % | chosen (not in the LACoFD table), above the 10 h value |
| Live | 62 % | LACoFD, 7 Jan 2025, Santa Monica Mountains live fuel moisture |

## Files

| File | What it is |
|---|---|
| `inputs` | the deck; the conductor block (from the `TRANSMISSION LINES` banner on) is the output of `fire_lines_towers.py` |
| `lines_config.json` | the two OpenStreetMap ways, the box each is cut to, the shortest span, and which line has frame towers |
| `osm_lines.json` | OpenStreetMap power lines in domain metres (output of `osm_power_lines.py`; ODbL, see below) |
| `lines_meta.json` | the picked towers: points, heights, span lengths and still-air clearances |
| `inflow_hrrr.txt` | inflow profile `z u v T` above the ground (output of `hrrr_profile.py`) |
| `sounding_hrrr.txt` | ERF input sounding from the same column, with the moisture column set to zero |
| `ignitions_spinup.csv` | the ignition schedule (`erf.fire.ignition.schedule_file`) |
| `case.py` | `prepare` builds a run directory; `post` makes the figures and the summary |

Shared tools in `Exec/Tools/`: `hrrr_profile.py`, `osm_power_lines.py`, `fire_lines_towers.py`,
`make_landfire_fuel_map.py`, `fire_lines_plot.py`, `fire_lines_fragility.py`, `fire_head_distance.py`.

## Data and acknowledgements

We thank the providers of the open data this case is built from:

- **NOAA High-Resolution Rapid Refresh (HRRR)** analysis, 2025-01-07 19 UTC, product `prs` and `sfc`, fetched
  from the NOAA Open Data Dissemination program on AWS (https://registry.opendata.aws/noaa-hrrr-pds/) with
  Herbie (B. K. Blaylock, https://github.com/blaylockbk/Herbie). NOAA data are open. Redistributed here: only
  the derived `inflow_hrrr.txt` and `sounding_hrrr.txt`. The GRIB files are not.
- **LANDFIRE** LF2024 FBFM40, U.S. Department of the Interior, Geological Survey, and U.S. Department of
  Agriculture, https://landfire.gov, fetched from the LANDFIRE Product Service (https://lfps.usgs.gov). Public
  domain. Not redistributed: `case.py prepare` downloads the window (about 2 MB) and maps it.
- **OpenStreetMap**, (c) OpenStreetMap contributors, https://www.openstreetmap.org/copyright. The power lines
  were read with the Overpass API (overpass-api.de; Overpass 0.7.62) for the box
  34.0063 to 34.2259 N, 118.6226 to 118.3663 W (query in `osm_power_lines.py`). `osm_lines.json` and
  `lines_meta.json` are derived databases and are available under the Open Database License (ODbL) 1.0; please
  credit "(c) OpenStreetMap contributors" wherever the lines are shown.
- **Terrain:** NASA Shuttle Radar Topography Mission 1 arc-second elevation (SRTM1), public domain, fetched with
  the `elevation` Python package and made into the two terrain files by `gen_palisades.py`; source and steps in
  `Exec/CanonicalTests/Fire/Real_Terrain/Palisades_Fire/README.md`. Not redistributed.
- **Imagery** under the maps (`case.py post` needs it): U.S. Geological Survey, The National Map, public
  domain, fetched by `Exec/Tools/make_map_movie.py --fetch usgs`. Not redistributed.
- **Fuel moisture:** Los Angeles County Fire Department, Daily Fire Danger Analysis, forecast for 7 January
  2025, https://fire.lacounty.gov/wp-content/uploads/2025/01/FDA-01-07-2025.pdf (Malibu Hills RAWS: 66 F,
  22 % RH, 8 % 10 h dead fuel moisture; Santa Monica Mountains live fuel moisture 62 %). Cited only; nothing is
  redistributed. We thank LACoFD for publishing it.
- **Observations:** F. Li et al., "Fire Spread, Intensity, and Emissions Observations by Multiple Satellites:
  The Southern California Wildfires of January 2025", *AGU Advances* (2026),
  https://doi.org/10.1029/2025AV002064, for the satellite spread speed (1.0 to 3.7 km/h, 0.6 to 2.3 mph) and the
  observed relative humidity below 20 %. Cited only; nothing is redistributed.

Sources of the methods: Cohen, J. D. and Deeming, J. E. (1985), *The National Fire-Danger Rating System: basic
equations*, USDA Forest Service General Technical Report PSW-82, for the equilibrium moisture; IEEE Std 738-2023
(doi:10.1109/IEEESTD.2023.10382442, cited, not reproduced) for the Drake conductor and its heat balance; EN
1993-1-2 for the steel; Rothermel (1972) and Andrews (2012) for the spread, as listed in
`Structures/Exec/FireUnderLine/README.md`; Scott and Burgan (2005), *Standard Fire Behavior Fuel Models: A
Comprehensive Set for Use with Rothermel's Surface Fire Spread Model*, USDA Forest Service RMRS-GTR-153, for the
FBFM40 fuel models; Andrews (2018), RMRS-GTR-371, for the BEHAVE multi-class spread (see
`Exec/RegTests/FireBehaveFuelMap/README.md`).

## How to run

Needs Python 3 with numpy, rasterio, pyproj, matplotlib, yt and scipy, ffmpeg for the movie (skipped
without it), and FLARE built with `-DERF_ENABLE_STRUCTURES=ON` against MoorDyn-C
(`Structures/setup_moordyn.sh`; the conductor model), `-DERF_ENABLE_FFT=ON` (the deck sets `erf.use_fft = 1`)
and MPI. Herbie and pandas are needed only to
remake the HRRR files; `gen_palisades.py` lists its own needs.

```bash
python3 case.py prepare                        # runs gen_palisades.py in run/ if the terrain is missing
cd run && mpiexec -n 4 <FLARE build>/flare inputs
python3 ../../../../Exec/Tools/make_map_movie.py --fetch usgs   # once: basemap_palisades_z14.tif (USGS imagery) and a map movie
cd .. && python3 case.py post --basemap run/basemap_palisades_z14.tif
```

`post` needs a basemap (any RGB GeoTIFF over the domain); `make_map_movie.py --fetch usgs` downloads the USGS
imagery tiles once and writes it.

`prepare` also checks that the deck's conductor block and `lines_meta.json` are the ones `fire_lines_towers.py`
builds from `lines_config.json` over the terrain; it stops if either differs. A line's last span may be shorter
than `min_span_m`, and when no node meets `min_span_m` and the clearance the tool takes the next node and prints a
warning (15 times on this steep terrain; the towers are then raised where a span still clears by less than 8 m). The 550 m and 8 m limits
hold for every span, and the tool stops if one cannot meet them. The towers
stand on the run's mesh surface (the terrain file at the bottom grid nodes, which the tool rebuilds from the deck),
and the 8 m is measured over the higher of that surface and the terrain file; see the note under Expected results.
The terrain files are not committed: `prepare` makes them with `gen_palisades.py` from the SRTM1 tiles. The tower
clearances depend on them, so on another machine (another tile copy, or other rasterio or scipy versions) the
check can stop with "towers differ" even though nothing is wrong. Then remake the conductor block and
`lines_meta.json` over the new terrain as described below and look at the clearances it prints before using them.
`post` writes into `run/figs/`: the burned-area map with the lines coloured by their peak temperature, a movie
(needs ffmpeg), the line and tower temperatures against time, per-span and per-tower CSV files,
`summary.json`, the fragility figure and the head distance.

To remake the inputs from the sources (not needed to run the case), in `run/` after `prepare` (the `T` path
assumes the default `--run`):

```bash
cd run
T=../../../../Exec/Tools
python3 $T/hrrr_profile.py --date 2025-01-07 --hour 19 --domain palisades_domain.json \
    --point-xy 9000 8000 --out-prefix pal19_ign --z-min 20 --z0 0.1
cp pal19_ign_inflow.txt inflow_hrrr.txt
awk 'NF { $3 = "0.0"; print }' pal19_ign_sounding.txt > sounding_hrrr.txt    # dry run: qv = 0
# osm_power_palisades.json: the Overpass query in osm_power_lines.py, box under Data and acknowledgements
python3 $T/osm_power_lines.py osm_power_palisades.json palisades_domain.json fuel_palisades_fbfm40.asc osm
python3 $T/fire_lines_towers.py ../lines_config.json --osm osm_lines.json --terrain terrain_palisades.txt \
    --fire-terrain terrain_palisades_fire.txt --deck inputs \
    --meta lines_meta_new.json > conductor_block.txt
```

To use the remade files, put `conductor_block.txt` in place of the conductor block of `../inputs` (from the
`TRANSMISSION LINES` banner to the end) and copy `lines_meta_new.json` to `../lines_meta.json`,
`inflow_hrrr.txt` and `sounding_hrrr.txt` to `..`; `prepare` then checks them.

The HRRR answer can change slightly if NOAA reprocesses the archive; the committed files are the ones the
results below were made with.

## Expected results

From the run made with this deck (4 ranks, 2 h 29 min wall time). That run wrote the fire plotfiles every 30 s
(`erf.fire.plot_per = 30`); the deck now writes them every 60 s to halve the output. The setting only sets the
output cadence, so the numbers below do not depend on it. The deck also spells three inputs by their current names
(`erf.compressibility_type`, `erf.timestepping_type`, `erf.most.pot_surf_temp`) where the run used the
deprecated ones; ERF maps them to the same settings.

| Quantity | Run | Observed |
|---|---|---|
| Head distance from the ignition at 60 min | 2.51 km (1.6 miles), so about 2.5 km/h (1.6 mph) | satellite spread 1.0 to 3.7 km/h (0.6 to 2.3 mph), Li et al. 2026 |
| Burned area at 105 min | 390 ha (965 acres) | not compared |
| Spans the fire reached | 0 of 35 | not compared |
| Peak conductor temperature | 21.9 C (71 F) | – |
| Peak tower steel temperature (L1) | 15.0 C (59 F) | – |

The mean head speed over the first hour, from the head distance, is inside the satellite range. It is not
steady: the head moves 2.2 km in the first 30 min (4.4 km/h, above the range), then stalls from 40 to 55 min. The fire runs south-west, away from both
lines, so no span is reached and the line and tower temperatures stay near the air temperature. The head
distance is the farthest burned cell from the ignition point in any direction, compared with the distance to
one observed place; it does not check the direction of spread. It is measured from the centre of the 100 m
radius ignition. `head_ros_median_ms` in `summary.json` is the median front speed from `fire_stats.csv`, not
the head distance over time; it is not compared.

Clearance: the towers stand on the run's terrain-following mesh surface (`terrain_palisades.txt` at nodes 300 m
apart), where the run also puts the fire's radiating sources. Every span clears the higher of that surface and
`terrain_palisades.txt` by at least 8 m in still air; towers are raised up to 33 m on L1 and 48 m on L2 for it.
The fire spreads on the sharper `terrain_palisades_fire.txt`, which lies up to 52 m above the smoothed
atmosphere terrain. Over it the lowest clearance is -17.7 m on L1 and -14.6 m on L2: parts of 12 of 18 L1 spans
and 10 of 17 L2 spans lie below the fire-terrain ground (`clearance_fire_terrain` in `lines_meta.json`). This
does not enter the heat transfer, which places towers, wire and fire sources on the mesh surface. The fire does
not reach any span in this run, and the case is not set up for a line exposure study.

## CTest

None. The full run takes about two hours on 4 ranks, longer than the regression suite allows. The case was
checked as follows instead:

- `case.py prepare` rebuilds the conductor block and the tower file from `lines_config.json`; it stops unless
  both match the committed `inputs` and `lines_meta.json` exactly.
- `Exec/Tools/test_fire_lines_towers.py` (CTest `FireLinesTowers_Tool`, registered when Python has numpy)
  checks `fire_lines_towers.py`: the terrain read on a plane over a non-square grid, the mesh surface, the
  clearance over flat ground against the closed form $30 - 2.5 - W L^2/(8T)$, the picks along a straight line,
  the tower raising over a higher surface, and the stops on a repeated node and on fewer than 3 towers. It fails
  on a transposed read, a flipped sag, a missing insulator offset and a missing guard.
- The deck runs 4 steps from a fresh `run/` with `amrex.parmparse.verbose = 1` and reports no unused inputs.
- `case.py post` on the finished run gives per-span, per-tower, head-distance and summary files identical to
  those behind the numbers above.
