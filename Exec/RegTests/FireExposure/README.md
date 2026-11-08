# FireExposure

What each building experiences as the fire passes: the per-structure
exposure diagnostics (`erf.fire.exposure.*`) on the FireHybridObstacles
scenario. Balbi 2020 on the reference wind, the three boxes and the street
non-burnable, the level-set wall extrapolation on, one exposure row per box
every 25 s. It is a demonstration deck, not a validation.

```
./run_exposure.sh /path/to/erf_exec
```

`MPIRUN="mpirun -np 4" ./run_exposure.sh ...` runs each variant in
parallel; `SKIP_RUN=1` only re-tabulates existing CSVs.

## The scenario

320 x 160 x 160 m at 5 m, fire grid at 1.25 m, neutral 8 m/s sounding with
a MOST surface layer, FM1 short grass at 6% moisture, passive coupling,
240 s. Three 20 m wide, 10 m tall boxes stand across the wind at
x = 180-200 m, at y = 30, 80 and 115 m (ids 1-3 in scan order); the fire is
ignited as a 15 m disc on the middle box's centreline, 35 m upwind of its
face. The wall band is the ring of burnable cells one fire cell (1.25 m)
wide around each footprint.

## The variants

- `noib`: flat ground, the boxes exist only in the fire grid.
- `ib`: the same heightmap drives immersed-forcing buildings, with the
  open-column wind weights on.
- `noib_spotting`: flat ground with Albini spotting on a fixed seed, so
  brands land on the footprints and the embers column is exercised.
- `noib_spotting_front`: the same with `erf.fire.spotting.launch_from =
  front`, so brands come from the fireline rather than the whole burned
  area.

## What the columns mean

For each box the script prints the last exposure row: the fraction of its
wall band burned, the first and last arrival of the front in the band and
their difference (how long the front spent passing the box), the largest
peak fireline intensity in the band, the mean and largest accumulated heat
load there in MJ/m², and the number of embers that landed on the footprint.
Arrival times are -1 for a box the front never reached.

## What to expect

- The middle box (id 2), on the ignition centreline, is reached first and
  carries the largest heat load; its band burns completely as the front
  wraps around it. The outer boxes are reached later by the flanks and
  their residence times are longer because the front passes them obliquely.
- With immersed-forcing buildings everything arrives later (the wind the
  fire reads is slowed) and the intensities are lower.
- Embers appear only in the spotting rows; the counts depend on the seed.
  Launching from the front gives far fewer brands at the same probability,
  because the fireline is about a thousand cells against tens of thousands
  behind it.

## Reference table

```
variant              id      x      y  burned t_first  t_last  resid  peak_kWm   HL_mean    HL_max  embers  landed
------------------- --- ------ ------ ------- ------- ------- ------ --------- --------- --------- ------- -------
noib                  1    190     30    0.61     107     220    113    2805.2     1.865     3.089       0       0
noib                  2    190     80    1.00      28     132    104    3368.5     3.089     3.089       0       0
noib                  3    190    115    1.00      40     195    156    3256.2     3.089     3.089       0       0
ib                    1    190     30    0.74     106     223    117    1435.0     2.220     3.089       0       0
ib                    2    190     80    1.00      40     165    125    1542.8     3.089     3.089       0       0
ib                    3    190    115    0.98      50     220    170    1461.1     2.997     3.089       0       0
noib_spotting         1    190     30    0.70      75     222    147    2984.0     2.111     3.089       1      34
noib_spotting         2    190     80    1.00      13      91     78    3543.3     3.089     3.089       3      34
noib_spotting         3    190    115    1.00      23     186    163    3425.9     3.089     3.089       1      34
noib_spotting_front   1    190     30    0.70      76     222    147    2983.3     2.111     3.089       2      10
noib_spotting_front   2    190     80    1.00      13     126    113    3543.3     3.089     3.089       1      10
noib_spotting_front   3    190    115    1.00      40     198    158    3256.2     3.088     3.089       0      10
```

Four ranks, 240 s, last row per box (written at 224.9 s), regenerated on
2026-10-09 on the validated code (the previous table, of 2026-09-11, is
reproduced by the code before it). Things to read out of it.

**The street moved on 2026-09-11.** Until hgopalan/ERF#411 the fuel map
reader put the first row of `fuel_map_street.asc` on the south edge, so the
non-burnable street ran at y = 30-35 m, across the first box, instead of at
125-130 m. That street cut the first box's band: 45 % of it burned, and the
front had passed it by 164 s. With the street in its place 62 % has burned
when the run ends and the front is still passing.

**The middle box is the one the head fire hits.** It is reached at 28 s, a
second after the arrival-time probe on its upwind face, its whole wall band
burns as the front wraps around it, and the front leaves its downwind side
at 128 s, two seconds after the downwind probe reports. The outer boxes are
reached by the flanks, later and obliquely, so their residence times are
longer; the flank is still passing the first box when the run ends.

**The heat load approaches the fuel's energy.** The largest heat load on
every box is 3.09 MJ/m², 98 % of the FM1 fuel load (0.166 kg/m²) times its
heat content (18.6 MJ/kg): a band cell that burned has released all but the
load still burning out at 225 s. The 3.16 MJ/m² the table carried until
2026-10 was the whole load, counted by a flux that ran ahead of the burnout.
The mean is lower where part of the band never burned.

**Immersed-forcing buildings halve the intensity at the walls** (1.4-1.5
against 2.8-3.4 MW/m) because the wind the Balbi model reads is slowed at
the walls, and every arrival is 10 s or so later.

**Spotting now lands brands inside the domain.** The lofting height is
Albini (1983, INT-309) eq. 2 with the Table 4 fit for FM1 and the reference
wind brought to 10 m (2026-11). At 3-3.5 MW/m this gives brands that land
22-31 m downwind of their source, where the earlier 12.2 I^(1/3) height
(about 150 m) sent them 175-191 m and mostly out of the domain.
`noib_spotting` lands 34 brands in 240 s, 5 of them on footprints;
`noib_spotting_front` lands 10, 3 on footprints. The spot fires bring the
arrivals forward: box 2 is reached at 13 s against 28 s without spotting,
and box 1 at 75 s against 107 s. The spot cells have their fuel capped at
5 % of the initial load. The table above was re-measured on 2026-11-08
with 2 ranks; the `noib` and `ib` rows, which have no spotting, were
not re-run.
