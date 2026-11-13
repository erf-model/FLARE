# Cascade stripes

## Purpose

Two wind stripes through the fragility tools with the line failures on. The case checks the path
from ERF runs to damage-state counts end to end:

1. **Shared spin-up.** The ensemble driver (`Structures/Tools/fragility_ensemble`, plan key
   `spinup`) runs the flow alone once per stripe, and every gust seed of that stripe restarts from
   its checkpoint with the lines added. Do the members restart from it, and do the lines start at
   the restart rather than coming out of the checkpoint?
2. **Failures.** With a footing allowable between the stripes' demands, does the low stripe hold
   and the high stripe lose towers (`conductors.tower_failure = release`) and release spans?
3. **Damage states.** Does the fragility fit (`Structures/Tools/fragility_fit`, plan key
   `damage_states`) count DS1 (first exceedance), DS2 (a line released), DS3 (a tower lost) and DS4
   (2 or more towers lost) per stripe as the runs give them?

The damage states are defined in "Damage states" in `Structures/Tools/fragility_fit/README.md`.

## Files

- `run_case.py` writes the driver's plan (`plan.json`) and the fit's plan (`fit_plan.json`) in
  `--out`, runs the driver (members in `--out/ensemble`), runs the fit (`--out/fit`) and, with
  `--check`, checks the four results listed in its docstring. `python3 run_case.py --help` lists the
  options.

The lines, the terrain and the flow are those of `Tests/test_files/Conductors_RandomGusts`: two
lines over hills, three lattice towers, random gusts (`conductors.gust_type = random`). The driver
scales the deck's inflow to 10 and 25 m/s at 30 m. Each member runs 5 steps of spin-up and 10 steps
with the lines (0.9 s steps).

## Capacities

Every capacity is 1e9 N except the towers' footing uplift, which is a **synthetic** allowable of
20 kN (`--uplift`). It is a value chosen for the test between the two stripes' demands, not a
design value or a supplier's.

It comes from a calibration run of the same plan with every capacity at 1e9 N, which only records
the peak demands (MoorDyn-C, one rank):

| Component | 10 m/s (seeds 7, 8) | 25 m/s (seeds 7, 8) |
|---|---|---|
| Span tension, every span | 20.2 to 20.5 kN | 20.3 to 22.0 kN |
| Footing compression | 24 to 33 kN | 63 to 103 kN |
| Footing uplift | 0 | 33 to 66 kN |

The span tension is mostly the conductors' own weight and moves by about 2 kN between the stripes,
so a conductor rating cannot separate them in this deck. DS2 here comes from the spans that failed
towers release.

## Data and acknowledgements

No data set is used or redistributed. The deck is the repository's `Conductors_RandomGusts` test.

## How to run

Build with `-DERF_ENABLE_STRUCTURES=ON` and MPI (the MoorDyn stub or MoorDyn-C). Then:

```bash
python3 run_case.py --exe path/to/flare --out cascade --check
```

The driver waits for the machine's 1-minute load to fall below 6 before each run (`--load-start`,
or `--no-load-gate` to run at once). The six runs (two spin-ups, four members) take 1 to 2 s each
on one rank; the CTest takes about 8 s. `--out` must not already hold `ensemble` or `fit`: the driver skips
members that finished, so an earlier run would be reused at another allowable.

## Expected results

Numbers from MoorDyn-C on one rank per run (November 13, 2026). Every member restarted from its
stripe's `chk00005` at 4.5 s, the conductor model reported that the checkpoint holds no conductor
state, and the conductor logs start at 4.5 s.

| Member | Stripe | Largest peak ratio | Towers lost | Lines released | First failure |
|---|---|---|---|---|---|
| seed 7 | 10 m/s | 3.3e-5 | 0 | 0 | none |
| seed 8 | 10 m/s | 3.0e-5 | 0 | 0 | none |
| seed 7 | 25 m/s | 3.30 | 3 | 2 | L1_t1 footing uplift, 65.9 kN at 5.4 s |
| seed 8 | 25 m/s | 2.91 | 3 | 2 | L1_t1 footing uplift, 42.3 kN at 5.4 s |

At 25 m/s seed 7 loses its towers at 5.4, 6.3 and 8.1 s; seed 8 loses all three at 5.4 s. In
the calibration run (no failures) the footing uplift at 25 m/s stays above 20 kN at most steps on
every tower, with peaks between 5.4 and 13.5 s, so the failures do not hang on the first step
after the restart alone. The towers in this deck are rigid.

The MoorDyn stub build gives the same outcomes: no failure at 10 m/s, 3 towers lost and 2 lines
released at 25 m/s in both seeds (largest peak ratios 3.41 and 2.91).

The fit gives DS1 to DS4 at 0 of 2 members at 10 m/s and 2 of 2 at 25 m/s. With outcomes
separated between the stripes the curves have no maximum-likelihood estimate (status
`separated`); the case tests the counts, not a fitted median.

Each member is one realisation: a falling span amplifies roundoff, so the times and the number of
towers lost at 25 m/s can differ between builds and rank counts. The check asks for 2 or more.

**The check fails when it should.** Two deliberate changes, each run through the whole case:

| Change | What the check reports |
|---|---|
| `--uplift 1e9` (no capacity reached) | at 25 m/s no line released, no tower lost, no row in `failures.dat`; the fit's DS1 to DS4 at 0 of 2 |
| The driver keeps the conductor keys in the spin-up deck (`SPINUP_DROPPED = ()`) | every member's run log lacks "holds no conductor state"; four spin-ups of one member each instead of two of two |

The second change is the one the spin-up design guards against. The spin-up then runs the lines
too, its towers fail during the spin-up (from 0.9 s at 25 m/s), and the members restart with that
state and its logs. Their damage table still shows 3 towers lost at 25 m/s and none at 10 m/s,
the same as the correct run. Outcomes alone cannot tell the two apart, so the check reads the
conductor model's message at the restart.

## CTest

`CascadeStripes_DamageStates` (labels `regression`, `conductors`, `failure`, `fragility`) runs the
case with `--no-load-gate --check` on one rank per run and passes when the case prints
`CascadeStripes: all checks hold`.
