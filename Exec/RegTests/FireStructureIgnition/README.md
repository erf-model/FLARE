# FireStructureIgnition

Structure ignition and house-to-house spread (`erf.fire.structures.ignition.*`)
on three houses: a grass fire is lit against the wall of house A, house B stands
20 m downwind of it and house C 90 m further. Forty steps of 0.5 s, one rank.
It is the old-against-new regtest of the option and the smoke test CTest runs;
the thresholds and the burn curve are scaled to the 20 s run and are not the
documented defaults (the `ignition` variant of
`Exec/CanonicalTests/Fire/WUI_Subdivision` runs those).

```
./run_structure_ignition.sh /path/to/erf_exec          # four runs, then the checks
SKIP_RUN=1 ./run_structure_ignition.sh x               # checks only
python3 gen_structure_ignition.py                      # rasters (already committed)
```

## The scenario

320 x 160 x 120 m, 10 m atmosphere cells, 5 m fire cells; a neutral 10 m/s
sounding entering at the west face and leaving at the east face, periodic in
y; MOST surface layer, Smagorinsky; FM1 short grass at 6% moisture; passive
atmosphere. Three 20 m square, 8 m tall houses on the centreline y = 70-90 m:
A at x = 100-120 m, B at 140-160 m, C at 250-270 m (ids 1-3 in scan order).
The heightmap is nodal on the 10 m atmosphere grid and sampled onto the fire
cells by nearest node, so each footprint spans 2.5 m beyond its edges; A and
B are two node spacings apart so that they stay separate footprints (adjacent
roof nodes are one 4-connected footprint). The ignition disc (centre x = 78 m,
radius 20 m) ends on A's wall band, so the band burns from the first step.
Seeded spotting launches from the fire front and from burning houses.

| deck | what is on |
|---|---|
| `inputs_off` | the old path: exposure diagnostics, structure ignition off |
| `inputs_on` | ignition on: heat load 0.12 MJ/m² (default 6), 5 embers (default 50), intensity criterion off; curve 200 kW/m² peak in 2 s, 1.5 MJ/m² load, half of it in the decay (test values; the defaults are placeholders), radiant fraction 0.3 within 60 m |
| `inputs_radiation` | `inputs_on` without spotting and with the ember criterion off: B can only ignite from A's radiation |

The script also restarts `inputs_on` from its step-20 checkpoint into its own
outputs, and runs the checker on `inputs_off`, where it must fail.

## The checks (`check_structure_ignition.py`)

- **A ignites from the heat at its wall**, cause `heat`, at 0.5 s: the grass at
  the wall meets the 0.12 MJ/m² threshold on the first step.
- **A second house ignites later.** In `on`, B at 2 s and C at 12 s, both by
  heat. In `radiation`, B at 8 s by heat while its wall band has never burned
  (`t_first_s = -1`, `wall_burned_frac = 0`): the ignition came from A's
  radiation alone, whose largest incident flux on B's band was 22.0 kW/m².
  C, 90 m further, stays unignited.
- **A burns**: its release in the CSV reaches the 200 kW/m² peak, it burns out
  at 13.5 s (the 1.5 MJ/m² curve lasts 12.6 s) and its release is zero after.
  No house ignites by embers in these runs any more; the ember criterion's
  precedence is covered by the unit tests.
- **The plotfile carries the state**: A's footprint holds state 2 and its
  ignition time; cells off the footprints hold 0 and -1.
- **The radiation kernel**: the incident flux in the step-10 plotfile equals,
  to 1e-9, the point-source sum `chi_r q dA / (2 pi max(r^2, (dx/2)^2))`
  recomputed in the checker from the footprints and the releases the CSV
  reports for that step.
- **Restart**: the run restarted from step 20 reproduces the last CSV row of
  every structure exactly and the structure fields of the step-40 plotfile
  exactly; the other fields to round-off (the inflow atmosphere of this deck
  restarts to about 1e-11 relative with or without structure ignition).
- **The old path fails the check**: the exposure CSV of `inputs_off` has no
  `state` column, so the first assertion fails and the script requires that.
  The old path also prints no `[FIRE STRUCTURE]` line.

## Reference results

RelWithDebInfo build, one rank, 2026-11-10, with the test curve above:

| variant | A ignites | B ignites | C ignites | A burns out | largest incident flux on B [kW/m²] |
|---|---|---|---|---|---|
| off | never | never | never | - | - |
| on | 0.5 s, heat | 2.0 s, heat | 12.0 s, heat | 13.5 s | - |
| radiation | 0.5 s, heat | 8.0 s, heat | never | 13.5 s | 22.0 |

The table of 2026-09-16 (A at 1.0 s, B at 3.5 s and C at 12 s by embers) was
already out of date before the curve's values changed: the same decks with
the earlier curve (250 kW/m², 2 MJ/m², 30 % in the decay) on the 2026-11-10
code give A at 0.5 s, B at 2 s and C at 8 s, all by heat, A burning out at
12.5 s, and B at 6.5 s in `radiation`. The fire code changed in between
(among others, the lofting height of the spotting model in November 2026).

## The CTests

`FireStructureIgnition_on` runs `inputs_on` for 40 steps on one rank (the
smoke test); `FireStructureIgnition` (MPI builds, not Windows) runs the
script; `FireStructureIgnition_no_exposure_abort` and
`FireStructureIgnition_curve_abort` pass when the run stops at start-up
without `erf.fire.exposure.enable` and with a growth time whose t-squared
phase alone would release more than (1 - decay_fraction) of the fuel load;
`FireStructureIgnition_decay_abort` when `decay_fraction` is 1.
