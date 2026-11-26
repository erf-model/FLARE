# Marshall Fire and two transmission lines

## Purpose

A demonstration of the full FLARE chain on a real fire: the Marshall Fire of 30 December 2021 in Boulder
County, Colorado, run from public data with two real transmission lines in the domain. It shows what a user
can set up from open sources in about an hour of work and about an hour of computing on 4 ranks:

- **Weather:** one HRRR analysis hour (18 UTC), the mean of three columns on the upwind (west) face, as the
  sounding and the inflow profile (`Exec/Tools/hrrr_profile.py`).
- **Fuel:** the LANDFIRE 2.0.0 (2016) FBFM40 map on the fire grid (`Exec/Tools/make_landfire_fuel_map.py`).
- **Lines:** two OpenStreetMap `power=line` ways, the 230 kV Leggett to Plains End / Valmont to Simms line
  (L1) and the 115 kV Valmont to Ralston / Eldorado line (L2), with towers picked along them over the terrain
  (`Exec/Tools/osm_power_lines.py`, `Exec/Tools/fire_lines_towers.py`).
- **Exposure:** the conductor temperature of every span (IEEE Std 738-2023 heat balance with the fire's
  radiation and the local air) and the steel temperature of L1's lattice towers (EN 1993-1-2 lumped heating).

It is a demonstration, not a validation. One number is compared with an observation: how far the fire's head
is from the ignition point 75 min after ignition. The dead and live fuel moistures are **chosen** values
(6/7/8 % dead, 30 % live, as in `Exec/CanonicalTests/Fire/idealized_fire/inputs_marshall_ideal_ysunew`), not
taken from a station record. The atmosphere is dry: the HRRR moisture is set to zero because the run has no
moisture model.

The setup:

| Item | Value |
|---|---|
| Domain | 25.6 km square, 80 x 80 x 50 cells, 10 m first cell stretched by 1.1 to an 11.6 km top |
| Terrain | `Exec/CanonicalTests/Fire/Real_Terrain/Marshall_Fire/terrain_marshall.txt` |
| Atmosphere | anelastic, FFT pressure solve, YSU boundary layer with 2D Smagorinsky, Rayleigh damping in the top 3 km |
| Inflow | west and north faces, HRRR 2021-12-30 18 UTC (pressure levels from 100 m above the ground up, the 10 m wind below them, a log law with z0 = 0.1 m to 2 m); outflow on the east and south faces |
| Fire | level set on a grid 8 times finer (40 m cells), Rothermel per fuel cell, Andrews WAF, lagged two-way coupling, spotting |
| Ignition | one 100 m circle at (5000, 11000) m, lit at 600 s after a 10 min wind spin-up |
| Run length | 5400 s: 10 min spin-up, then 80 min of fire |
| L1 | 6 spans of Drake ACSR, 253-511 m long, lattice towers 30 m high with steel temperature |
| L2 | 14 spans of Drake ACSR, 210-512 m long, rigid point towers (wire temperature only) to keep the run short |

## Files

| File | What it is |
|---|---|
| `inputs` | the deck; the conductor block (from the `TRANSMISSION LINES` banner on) is the output of `fire_lines_towers.py` |
| `lines_config.json` | the two OpenStreetMap ways, the box each is cut to, the shortest span, and which line has frame towers |
| `osm_lines.json` | OpenStreetMap power lines in domain metres (output of `osm_power_lines.py`; ODbL, see below) |
| `lines_meta.json` | the picked towers: points, heights, span lengths and still-air clearances |
| `inflow_hrrr.txt` | inflow profile `z u v T` above the ground (output of `hrrr_profile.py`, mean of 3 columns) |
| `sounding_hrrr.txt` | ERF input sounding from the same columns, with the moisture column set to zero |
| `ignitions_spinup.csv` | the ignition schedule (`erf.fire.ignition.schedule_file`) |
| `case.py` | `prepare` builds a run directory; `post` makes the figures and the summary |

Shared tools in `Exec/Tools/`: `hrrr_profile.py`, `osm_power_lines.py`, `fire_lines_towers.py`,
`make_landfire_fuel_map.py`, `fire_lines_plot.py`, `fire_lines_fragility.py`, `fire_head_distance.py`.

## Data and acknowledgements

We thank the providers of the open data this case is built from:

- **NOAA High-Resolution Rapid Refresh (HRRR)** analysis, 2021-12-30 18 UTC, product `prs` and `sfc`, fetched
  from the NOAA Open Data Dissemination program on AWS (https://registry.opendata.aws/noaa-hrrr-pds/) with
  Herbie (B. K. Blaylock, https://github.com/blaylockbk/Herbie). NOAA data are open. Redistributed here: only
  the derived `inflow_hrrr.txt` and `sounding_hrrr.txt` (8 KB each). The GRIB files are not.
- **LANDFIRE** 2.0.0 (2016) FBFM40, U.S. Department of the Interior, Geological Survey, and U.S. Department of
  Agriculture, https://landfire.gov, fetched from the LANDFIRE Product Service (https://lfps.usgs.gov). Public
  domain. Not redistributed: `case.py prepare` downloads the window (about 2 MB) and maps it.
- **OpenStreetMap**, (c) OpenStreetMap contributors, https://www.openstreetmap.org/copyright. The power lines
  were read with the Overpass API (overpass-api.de; Overpass 0.7.62) for the box
  39.8538 to 40.0848 N, 105.2901 to 104.9908 W (query in `osm_power_lines.py`). `osm_lines.json` and
  `lines_meta.json` are derived databases and are available under the Open Database License (ODbL) 1.0; please
  credit "(c) OpenStreetMap contributors" wherever the lines are shown.
- **Terrain:** the committed Real_Terrain file; its source is given in
  `Exec/CanonicalTests/Fire/Real_Terrain/Marshall_Fire/README.md`.
- **Imagery** under the maps (optional, `case.py post`): U.S. Geological Survey, The National Map, public
  domain, fetched by `Exec/Tools/make_map_movie.py --fetch usgs`. Not redistributed.
- **Observation:** The Colorado Sun, Marshall Fire timeline, 6 January 2022,
https://coloradosun.com/2022/01/06/marshall-fire-boulder-county-timeline/: the fire reached the Superior
Costco, about 5.5 km (3.4 miles) from the origin area, about 75 min after it started.

Sources of the methods: IEEE Std 738-2023 (doi:10.1109/IEEESTD.2023.10382442, cited, not reproduced) for the
Drake conductor and its heat balance; EN 1993-1-2 for the steel; Rothermel (1972) and Andrews (2012) for the
spread, as listed in `Structures/Exec/FireUnderLine/README.md`.

## How to run

Needs Python 3 with numpy, rasterio, pyproj, matplotlib, yt and scipy, ffmpeg for the movie (skipped
without it), and a FLARE build with MoorDyn-C (the conductor model). Herbie and pandas are needed only to
remake the HRRR files.

```bash
python3 case.py prepare                        # run/ with the deck, terrain, small files and the LANDFIRE fuel map
cd run && mpiexec -n 4 <FLARE build>/flare inputs
cd .. && python3 case.py post --basemap <RGB GeoTIFF over the domain>
```

`prepare` also checks that the deck's conductor block and `lines_meta.json` are the ones `fire_lines_towers.py`
builds from `lines_config.json` over the terrain; it stops if either differs. A line's last span may be shorter
than `min_span_m`; the 550 m and 8 m limits hold for every span, and the tool stops if one cannot meet them.
`post` writes into `run/figs/`: the burned-area map with the lines coloured by their peak temperature, a movie
(needs ffmpeg), the line and tower temperatures against time, per-span and per-tower CSV files,
`summary.json`, the fragility figure and the head distance.

To remake the inputs from the sources (not needed to run the case), in `run/` after `prepare`:

```bash
cd run
T=../../../../Exec/Tools
python3 $T/hrrr_profile.py --date 2021-12-30 --hour 18 --domain marshall_domain.json \
    --point-xy 0 11000 --point-xy 0 8000 --point-xy 0 14000 --out-prefix marshall_hrrr18 \
    --z-min 100 --z0 0.1
cp marshall_hrrr18_inflow.txt inflow_hrrr.txt
awk 'NF { $3 = "0.0"; print }' marshall_hrrr18_sounding.txt > sounding_hrrr.txt    # dry run: qv = 0
# osm_power_marshall.json: the Overpass query in osm_power_lines.py, box under Data and acknowledgements
python3 $T/osm_power_lines.py osm_power_marshall.json marshall_domain.json fuel_marshall_fbfm40.asc osm
python3 $T/fire_lines_towers.py ../lines_config.json --osm osm_lines.json --terrain terrain_marshall.txt \
    > conductor_block.txt
```

The HRRR answer can change slightly if NOAA reprocesses the archive; the committed files are the ones the
results below were made with.

## Expected results

From the run made with this deck (FLARE `af41c74b0`, 4 ranks, 47 min wall time):

| Quantity | Run | Observed |
|---|---|---|
| Head distance from the ignition at 75 min | 4.75 km (3.0 miles) | about 5.5 km (3.4 miles), Colorado Sun |
| Burned area at 80 min | 452 ha (1,117 acres) | not compared |
| Spans the fire reached | 6 of 20 | not compared |
| Peak conductor temperature | 48.7 C (120 F) | – |
| Peak tower steel temperature (L1) | 29.3 C (85 F) | – |

The head reaches 86 % of the observed distance at 75 min. The burned area then ends at L1 and the edge of the
town, which LANDFIRE maps as non-burnable. The head distance is the farthest burned cell from the ignition
point in any direction, compared with the distance to one observed place; it does not check the direction of
spread. No span passes the chosen 100 C wire flag and no tower the 400 C steel flag
(`fire_lines_fragility.py`), so no fragility curve is fitted; the figure shows the margins instead.

## CTest

None. The full run takes about an hour on 4 ranks, longer than the regression suite allows. The case was
checked as follows instead:

- `case.py prepare` rebuilds the conductor block and the tower file from `lines_config.json`; it stops unless
  both match the committed `inputs` and `lines_meta.json` exactly.
- `Exec/Tools/test_fire_lines_towers.py` checks the terrain read of `fire_lines_towers.py` on a plane over a
  non-square grid; it fails on a transposed read.
- The deck runs 4 steps from a fresh `run/` with `amrex.parmparse.verbose = 1` and reports no unused inputs.
- `case.py post` on the finished run gives per-span, per-tower, head-distance and summary files identical to
  those behind the numbers above.
