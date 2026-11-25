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

It is a demonstration, not a validation. The head speed is compared with a satellite estimate and the burned
area with the first aerial map; the perimeter is not compared. The atmosphere is dry: the HRRR moisture is
set to zero because the run has no moisture model. In this run the fire spreads downwind toward the coast,
away from both lines, so the lines stay near the air temperature; the case shows the chain working, not a
line exposure.

The setup:

| Item | Value |
|---|---|
| Domain | 24 km square, 80 x 80 x 102 cells, 10 m first cell stretched by 1.03 to a 6.5 km top |
| Terrain | `terrain_palisades.txt` and `terrain_palisades_fire.txt`, made by `Exec/CanonicalTests/Fire/Real_Terrain/Palisades_Fire/gen_palisades.py` |
| Atmosphere | anelastic, FFT pressure solve, YSU boundary layer with 2D Smagorinsky, Rayleigh damping in the top 1.5 km |
| Inflow | north and east faces, HRRR 2025-01-07 19 UTC; outflow on the west and south faces |
| Fire | level set on a grid 10 times finer (30 m cells) on the fire terrain, BEHAVE multi-class spread per fuel cell, Andrews WAF, lagged two-way coupling, spotting |
| Ignition | one 100 m circle at (9000, 8000) m, lit at 600 s after a 10 min wind spin-up |
| Run length | 6900 s: 10 min spin-up, then 105 min of fire |
| L1 | 13 spans of Drake ACSR, 142-541 m long, lattice towers 30 m high with steel temperature |
| L2 | 17 spans of Drake ACSR, 208-541 m long, rigid point towers (wire temperature only); six towers raised to 37-52 m so that every span clears the terrain by 8 m |

Fuel moisture:

| Class | Value | Source |
|---|---|---|
| 1 h dead | 4 % | NFDRS equilibrium moisture (Cohen and Deeming 1985) at the Malibu Hills RAWS forecast of 22 % RH and 66 F (LACoFD, 7 Jan 2025): 4.8 % |
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
- **Terrain:** SRTM1 elevation, made into the two terrain files by `gen_palisades.py`; source and steps in
  `Exec/CanonicalTests/Fire/Real_Terrain/Palisades_Fire/README.md`. Not redistributed.
- **Imagery** under the maps (optional, `case.py post`): U.S. Geological Survey, The National Map, public
  domain, fetched by `Exec/Tools/make_map_movie.py --fetch usgs`. Not redistributed.
- **Fuel moisture:** Los Angeles County Fire Department, Daily Fire Danger Analysis, 7 January 2025. Cited, not
  redistributed. We thank LACoFD for publishing it.
- **Observations:** Li et al. (2026), *AGU Advances*, doi:10.1029/2025AV002064, for the satellite spread speed
  (1.0 to 3.7 km/h, 0.6 to 2.3 mph) and the observed relative humidity below 20 %; the Los Angeles Fire
  Department after-action report on the Palisades Fire, for the FIRIS aerial map of 771 acres at about
  14:11 PST, 3 h 40 min after ignition.

Sources of the methods: Cohen, J. D. and Deeming, J. E. (1985), *The National Fire-Danger Rating System:
basic equations*, USDA Forest Service General Technical Report PSW-82, for the equilibrium moisture; IEEE Std
738-2023 (doi:10.1109/IEEESTD.2023.10382442, cited, not reproduced) for the Drake conductor and its heat
balance; EN 1993-1-2 for the steel; Rothermel (1972) and Andrews (2012) for the spread, as listed in
`Structures/Exec/FireUnderLine/README.md`.

## How to run

Needs Python 3 with numpy, rasterio, pyproj, matplotlib, yt and scipy (Herbie only to remake the HRRR files;
`gen_palisades.py` lists its own needs), and a FLARE build with MoorDyn-C (the conductor model).

```bash
python3 case.py prepare                        # runs gen_palisades.py if the terrain is missing, then run/
cd run && mpiexec -n 4 <FLARE build>/flare inputs
cd .. && python3 case.py post --basemap <RGB GeoTIFF over the domain>
```

`prepare` also checks that the deck's conductor block is the one `fire_lines_towers.py` builds from
`lines_config.json`; it stops if they differ. `post` writes into `run/figs/`: the burned-area map with the
lines coloured by their peak temperature, a movie, the line and tower temperatures against time, per-span and
per-tower CSV files, `summary.json`, the fragility figure and the head distance.

To remake the inputs from the sources (not needed to run the case), in `run/` after `prepare`:

```bash
cd run
T=../../../../Exec/Tools
python3 $T/hrrr_profile.py --date 2025-01-07 --hour 19 --domain palisades_domain.json \
    --point-xy 9000 8000 --out-prefix pal19_ign
cp pal19_ign_inflow.txt inflow_hrrr.txt
awk 'NF { $3 = "0.0"; print }' pal19_ign_sounding.txt > sounding_hrrr.txt    # dry run: qv = 0
# osm_power_palisades.json: the Overpass query in osm_power_lines.py, box under Data and acknowledgements
python3 $T/osm_power_lines.py osm_power_palisades.json palisades_domain.json fuel_palisades_fbfm40.asc osm
python3 $T/fire_lines_towers.py ../lines_config.json --osm osm_lines.json --terrain terrain_palisades.txt \
    > conductor_block.txt
```

The HRRR answer can change slightly if NOAA reprocesses the archive; the committed files are the ones the
results below were made with.

## Expected results

From the run made with this deck (FLARE `af41c74b0`, 4 ranks, 1 h 56 min wall time):

| Quantity | Run | Observed |
|---|---|---|
| Head distance from the ignition at 60 min | 2.56 km (1.6 miles), so about 2.6 km/h (1.6 mph) | satellite spread 1.0 to 3.7 km/h (0.6 to 2.3 mph), Li et al. 2026 |
| Burned area at 105 min | 403 ha (997 acres) | FIRIS map 771 acres at 3 h 40 min, LAFD |
| Spans the fire reached | 0 of 30 | not compared |
| Peak conductor temperature | 23.7 C (75 F) | – |
| Peak tower steel temperature (L1) | 15.2 C (59 F) | – |

The head speed in the first hour is inside the satellite range. The burned area is ahead of the first aerial
map: the run passes 771 acres well before 3 h 40 min. Possible reasons (one HRRR hour for the whole run, a dry
atmosphere, the fuel map's resolution) are not tested here. The fire runs south-west, away from both lines,
so no span is reached and the line and tower temperatures stay near the air temperature.

## CTest

None. The full run takes about two hours on 4 ranks, longer than the regression suite allows. The case was
checked as follows instead:

- `case.py prepare` rebuilds the conductor block and the tower file from `lines_config.json`; both match the
  committed `inputs` and `lines_meta.json` exactly.
- The deck runs 4 steps from a fresh `run/` with `amrex.parmparse.verbose = 1` and reports no unused inputs.
- `case.py post` on the finished run gives per-span, per-tower, head-distance and summary files identical to
  those behind the numbers above.
