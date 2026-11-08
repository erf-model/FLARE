# FireCustomFuel

A fire on flat ground whose fuel model is described in the input deck
(`erf.fire.custom_fuel.*`) instead of coming from the Anderson 13 or the Scott
and Burgan 40 compiled tables. Fuel codes 1000-1015 take the fuel table slots
above the published sets, so every rate-of-spread model, the heat flux, the
fuel load and the burnout read them exactly as they read a published model.

Properties are given in SI and converted once to the units
`FuelModelParams` carries internally.

## Decks

| Deck | What it is |
| --- | --- |
| `inputs_base` | 320 x 160 x 160 m, 5 m cells, 1.25 m fire cells, one-way coupling, Anderson short grass, Rothermel |
| `inputs_anderson1` | the baseline: short grass from the compiled table |
| `inputs_custom_grass` | the same fuel written out in SI as code 1000; must reproduce the baseline |
| `inputs_custom_heavy` | a heavy bed of coarse woody material, which no published set can express |
| `inputs_custom_map` | a raster mixing Scott-Burgan GR2, a deck-defined block and the non-burnable codes, each cell starting with its own model's load |
| `inputs_bad_code` | a code outside 1000-1015 |
| `inputs_bad_missing` | a block missing a required property |
| `inputs_bad_depth` | a bed depth under the 0.01 m floor, where the Balbi models silently return zero spread |
| `inputs_bad_heat` | a heat content left in BTU/lb instead of J/kg |
| `inputs_bad_burnout` | `burnout_model = sfire` on a deck-defined model with no burn time |
| `inputs_bad_undeclared` | a raster holding a custom code the deck never defines |
| `inputs_bad_uniform_code` | `fuel_model_id` in the custom range with no block defining it |
| `inputs_custom_map_altid` | the mixed map with Anderson 13 as the uniform `fuel_model_id`, which a per-fuel run over a `load_from_map` raster must ignore |
| `inputs_undeclared_nonburnable` | the undeclared-code raster with that code declared non-burnable, which must run |

`fuel_map_mixed.asc` and `fuel_map_undeclared.asc` come from
`python3 make_fuel_map.py`.

## Running

```
MPIRUN="mpirun -np 2" FCOMPARE=/path/to/amrex_fcompare ./run_custom_fuel.sh /path/to/erf_exec
```

`SKIP_RUN=1` reruns only the checks. `FCOMPARE` adds the box-parity comparison.

## What the checks establish

`check_custom_fuel.py --selftest` exercises the checker's own pass/fail logic
first, with values just inside and just outside each band.

* **identity** — `custom_grass` reproduces `anderson1` over all 1920
  `[FIRE DEBUG]` numbers with a worst relative difference of 0. The fire
  plotfiles agree to 8.5e-14 absolute and 5.4e-16 relative, which is the
  round-off of the eight-digit SI values the deck carries, not bitwise
  equality: `0.16600262` is the written form of `0.034 lb/ft2 x 4.88243`.
* **summary** — the start-up line reports the deck's own SI back, at slot 54,
  which is above every published slot.
* **fuel** — the map run starts with 193430.6 kg on the grid, the sum of each
  cell's model load; 96% of it comes from the deck-defined block, so a run that
  ignored the block could not pass.
* **crossed** — the ignition circle straddles the block edge, so one front
  burns in both fuels. At 60 s it has burned 322 cells of the deck-defined code
  at 0.06003 m/s and 246 of GR2 at 0.2753 m/s: the deck's properties reach the
  per-cell kernel, not just the uniform path. The front is asymmetric because
  the sounding's wind carries the head fire east into the block and leaves a
  backing fire in GR2, so the burned counts are a wind effect, not a fuel one;
  the rates are the fuel one. Needs `yt`.

  The deck-defined rate is 0.05287 m/s on the binary before the wind adjustment
  factor was taken per cell: the block's 2 m bed gives 0.560 against the 1 ft
  uniform model's 0.362, so its midflame wind rose by 55% and its rate by 3.4%.
  GR2's own bed is 1 ft, so its rate is unchanged.
* **slower** — the coarse bed spreads at 0.06003 m/s against grass's 0.4625 m/s.
  The deck's properties drive the spread; they are not a published model in
  disguise.
* **box parity** — the map run on the default decomposition and on one box
  agree bitwise, which covers the per-cell slot lookup and the cross-rank
  reduction in the raster check.
* **uniform id is not read** — `custom_map_altid` differs from `custom_map`
  only in `erf.fire.fuel_model_id`, which a run with
  `erf.fire.rothermel_per_fuel` over a `load_from_map` raster must ignore
  entirely. Every number agrees over 1920 lines, with **no exemptions**, and
  `fcompare` finds no difference in any fire plotfile field.

  The deck uses Anderson 13, the furthest the Anderson set gets from the deck's
  Anderson 1: 3 ft of bed against 1 ft and 13.03 kg/m2 of load against
  0.166 kg/m2. Each of the three paths that used to read the uniform model
  fails this check on its own:

  | Path | How it failed | Where |
  | --- | --- | --- |
  | Rothermel coefficients | the level-set front followed the uniform rate: 608 against 704 cells at 30 s | fixed in #444 |
  | wind adjustment factor | built from the uniform bed depth: 0.362 at 1 ft against 0.459 at 3 ft, 27% of the midflame wind, so the two runs burned 558 against 566 cells at 60 s | fixed here |
  | Byram intensity and flame length | built from the uniform initial load: `I_B_max` differed by 37076 kW/m, the whole grid clamped to exactly 0 with Anderson 1 | fixed here |

  The Byram row is why the deck could not use Anderson 13 before: #444 had to
  fall back to Anderson 10, which shares Anderson 1's 1 ft depth, to isolate the
  coefficients from the wind adjustment factor, and still exempt `I_B_max` and
  `L_max` by name.
* **aborts** — each of the six bad decks stops the run with a message naming
  the input. `bad_undeclared` was also run on one and two ranks: its reduction
  aborts on both rather than hanging on one.

* **model sweep** — ten 10 s runs pair grass against the coarse deck fuel for
  every rate-of-spread model. Each has to give a different rate, since each
  reads its fuel through the same `FuelModelParams`:

  | `ros_model` | grass [m/s] | deck fuel [m/s] |
  | --- | --- | --- |
  | `rothermel` | 0.6746785108 | 0.06343788554 |
  | `balbi` | 0.4534711647 | 0.005781985511 |
  | `behave` | 0.6746785108 | 0.08237936874 |
  | `macarthur` | 1.793275906 | 6 |
  | `cheney_gould` | 1.343976782 | 1.343976782 |

  MacArthur rises rather than falls: it is a grassland correlation driven by
  the load and the bed depth, so a deeper, heavier bed speeds it up where
  Rothermel's reaction term and Balbi's radiation, which read the
  surface-area-to-volume ratio, slow down. Cheney-Gould 1998 reads no fuel bed
  property at all: its rate comes from the 10 m wind, the curing, the moisture
  and `erf.fire.cheney_gould.pasture`. The check therefore requires it to give
  the **same** rate for both fuels. Before this was changed the script asked it
  to differ, and that check failed on the base binary too (1.2112 m/s for both
  fuels); its rate is now higher because the reference wind is brought from
  6.1 m to 10 m.

## Measured numbers

| variant | active fire cells at 60 s | max ROS [m/s] | fuel at 0 s [kg] | fuel at 60 s [kg] |
| --- | --- | --- | --- | --- |
| `anderson1` | 926 | 0.462504577 | 8497.1 | 8268.3 |
| `custom_grass` | 926 | 0.462504577 | 8497.1 | 8268.3 |
| `custom_heavy` | 552 | 0.06002690628 | 1126391.0 | 1122149.7 |
| `custom_map` | 568 | 0.2752908476 | 193425.1 | 191018.7 |
| `custom_map_altid` | 568 | 0.2752908476 | 193425.1 | 191018.7 |

These numbers were re-measured with the unit-audit changes. The grass rates
were already higher than the earlier table on the base binary (0.6747 m/s for
grass in the 10 s Rothermel sweep on both binaries), so that part of the change
predates the audit. The deck-defined rows moved because the deck's 1 h
surface-area-to-volume ratio is now 328.084 1/m (100 1/ft, the lowest value
the range check accepts) instead of 200 1/m.

`custom_map` is re-measured twice over. With the level-set path reading the
per-cell coefficients (#444) it burned 558 cells at 60 s where it burned 774
before, since the front no longer spreads at the uniform grass rate inside the
coarse block. With the wind adjustment factor taken per cell it burned 564 (568 with the current deck): the
block's 2 m bed sees more of the wind than the 1 ft uniform model allowed it.

`custom_map`'s max ROS is GR2's rate, the fastest fuel on its grid; the
deck-defined block spreads at 0.06003 m/s in the same run.

`custom_map_altid` is listed because its agreeing with `custom_map` to the last
digit is the point of the deck; it differed in the two Byram diagnostics and in
the cell count until the two leaks above were closed.
