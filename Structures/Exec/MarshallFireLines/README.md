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
| Inflow | west and north faces, HRRR 2021-12-30 18 UTC (pressure levels from 20 m above the ground up, the 10 m wind below them, a log law with z0 = 0.1 m to 2 m); outflow on the east and south faces |
| Fire | level set on a grid 8 times finer (40 m cells), Rothermel per fuel cell, Andrews WAF, lagged two-way coupling, spotting |
| Ignition | one 100 m circle at (5000, 11000) m, lit at 600 s after a 10 min wind spin-up |
| Run length | 5400 s: 10 min spin-up, then 80 min of fire |
| L1 | 6 spans of Drake ACSR, 253-511 m long, lattice towers 30 m high with steel temperature |
| L2 | 14 spans of Drake ACSR, 240-510 m long, rigid point towers (wire temperature only) to keep the run short |

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
  the derived `inflow_hrrr.txt` and `sounding_hrrr.txt` (about 7 KB and 6 KB). The GRIB files are not.
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
- **Reference** for the wind bias: R. G. Fovell, M. J. Brewer and R. J. Garmong, "The December 2021 Marshall
  Fire: Predictability and Gust Forecasts from Operational Models", Atmosphere 13, 765 (2022),
  https://doi.org/10.3390/atmos13050765, open access under CC BY 4.0. Cited only; nothing is redistributed.
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
than `min_span_m`, and when no node meets `min_span_m` and the clearance the tool takes the next node and prints a
warning. The 550 m and 8 m limits hold for every span, and the tool stops if one cannot meet them. The towers
stand on the run's mesh surface (the terrain file at the bottom grid nodes, which the tool rebuilds from the deck),
and the 8 m is measured over the higher of that surface and the terrain file; see the note under Expected results.
`post` writes into `run/figs/`: the burned-area map with the lines coloured by their peak temperature, a movie
(needs ffmpeg), the line and tower temperatures against time, per-span and per-tower CSV files,
`summary.json`, the fragility figure and the head distance.

To remake the inputs from the sources (not needed to run the case), in `run/` after `prepare`:

```bash
cd run
T=../../../../Exec/Tools
python3 $T/hrrr_profile.py --date 2021-12-30 --hour 18 --domain marshall_domain.json \
    --point-xy 0 11000 --point-xy 0 8000 --point-xy 0 14000 --out-prefix marshall_hrrr18 \
    --z-min 20 --z0 0.1
cp marshall_hrrr18_inflow.txt inflow_hrrr.txt
awk 'NF { $3 = "0.0"; print }' marshall_hrrr18_sounding.txt > sounding_hrrr.txt    # dry run: qv = 0
# osm_power_marshall.json: the Overpass query in osm_power_lines.py, box under Data and acknowledgements
python3 $T/osm_power_lines.py osm_power_marshall.json marshall_domain.json fuel_marshall_fbfm40.asc osm
python3 $T/fire_lines_towers.py ../lines_config.json --osm osm_lines.json --terrain terrain_marshall.txt --deck inputs \
    --meta lines_meta_new.json > conductor_block.txt
```

To use the remade files, put `conductor_block.txt` in place of the conductor block of `../inputs` (from the
`TRANSMISSION LINES` banner to the end) and copy `lines_meta_new.json` to `../lines_meta.json`,
`inflow_hrrr.txt` and `sounding_hrrr.txt` to `..`; `prepare` then checks them.

The HRRR answer can change slightly if NOAA reprocesses the archive; the committed files are the ones the
results below were made with.

## Expected results

From the run made with this deck (4 ranks, 50 min wall time):

| Quantity | Run | Observed |
|---|---|---|
| Head distance from the ignition at 75 min | 4.79 km (3.0 miles) | about 5.5 km (3.4 miles), Colorado Sun |
| Burned area at 80 min | 460 ha (1,137 acres) | not compared |
| Spans the fire reached | 6 of 20 | not compared |
| Peak conductor temperature | 48.6 C (119 F) | – |
| Peak tower steel temperature (L1) | 29.2 C (85 F) | – |

The head reaches 87 % of the observed distance at 75 min. The burned area then ends just past L1, near L2 and
the edge of the town, which LANDFIRE maps as non-burnable. The head distance is the farthest burned cell from the ignition
point in any direction, compared with the distance to one observed place; it does not check the direction of
spread. No span passes the chosen 100 C wire flag and no tower the 400 C steel flag
(`fire_lines_fragility.py`), so no fragility curve is fitted; the figure shows the margins instead. The head
distance is measured from the centre of the 100 m ignition circle.

Why the head falls short. The run is not tuned to close the gap; the likely reasons are:

- **The HRRR wind is a 3 km mean.** The fire burned in a downslope windstorm, and the strongest winds stayed
  within about 25 km of the foothills. Fovell et al. (2022) note that operational models tend to underpredict
  the winds at the windiest sites, that the HRRR missed the fastest gusts of this event, and that downslope
  winds depend strongly on grid spacing. The inflow here is one HRRR hour, so it may carry the same low bias.
- **The wind is steady.** The inflow is a mean profile held for the whole run. Gusts up to 51 m/s (115 mph) were
  reported about 11 km south of the ignition (Fovell et al. 2022); they are not in it, and the spread rate
  grows faster than linearly with wind.
- **The fuel moistures are chosen**, not measured (see Purpose).

How much wind would matter, as a back-of-the-envelope estimate. The spread rate is Rothermel's
$R = R_0 (1 + \phi_w + \phi_s)$ with $\phi_w = C\,U^B (\beta/\beta_{op})^{-E}$ and $B = 0.02526\,\sigma^{0.54}$
(Rothermel 1972, eqs. 47-50; $U$ the midflame wind, $\sigma$ the surface-area-to-volume ratio in 1/ft). The
wind limit is off in this deck (`erf.fire.use_wind_limit = false`). At the head the run's midflame wind is
about 6 m/s over GR2 grass (median of the burning cells 20, 40 and 60 min after ignition), so
$\phi_w \approx 100$ against a median $\phi_s$ of 0.2 to 0.4, and multiplying $U$ by $k$ multiplies $R$ by $k^B$
to within about 1 %. For the grass fuels here ($\sigma$ = 1800 to 2000 /ft), $B$ = 1.45 to 1.53:

| Wind multiplied by | Spread rate multiplied by |
|---|---|
| 1.10 | 1.15 to 1.16 |
| 1.20 | 1.30 to 1.32 |
| 1.36 | 1.56 to 1.60 |

If the head distance grew with the spread rate, a wind about 10 % higher, $(5.5/4.79)^{1/1.5} = 1.10$, would close
the gap, and a wind 36 % higher would carry the head to 7.5 to 7.7 km, past the observed 5.5 km. This ignores
spotting, the town edge and the slope, so it gives the size of the effect only. It shows the result is
sensitive to the wind: a small low bias in the mean wind over the fire's path is enough to explain the gap.

Clearance: the towers stand on the run's terrain-following mesh surface (`terrain_marshall.txt` at nodes 320 m
apart), where the run also puts the fire's radiating sources. Every span clears the higher of that surface and
`terrain_marshall.txt`, which the fire also burns on, by at least 8 m in still air: 10.2 m for L1 and 9.5 m for L2
at the lowest. All towers are 30 m.

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
