# Damage-state curves for the power lines over hills

## Purpose

Fragility curves of the damage states DS1 to DS4 for the line network of `Structures/Exec/PowerLines`,
from an ensemble of 7 wind stripes x 3 gust seeds run through ERF with the line failures on. The case
shows the whole chain on a full-size deck: the fragility ensemble driver's shared spin-up per stripe,
the members restarting with the lines, the failure logs, the damage table, and the fragility fit's
damage-state curves. The footing allowables are **synthetic**, so the curves are a demonstration of
the method, not the fragility of a real line.

## Files

- `powerlines_stripes.i` is `Structures/Exec/PowerLines` in one stage (`flow.inputs`,
  `network.inputs` and `inputs_lines` in one deck, the swaying towers kept, no restart).
- `run_case.py` copies the deck and the PowerLines case's `terrain_hills.txt`, `inflow_profile` and
  `input_sounding` into `--out/base`, writes the driver's plan (`plan.json`) and the fit's plan
  (`fit_plan.json`), and runs both. `python3 run_case.py --help` lists the options.

The plan, as `run_case.py` writes it:

| Setting | Value |
|---|---|
| Stripes | 11 to 17 m/s at 10 m, one per m/s, direction 0 degrees |
| Seeds | 1, 2, 3 (`conductors.gust_type = random`) |
| Spin-up | 1715 steps (600 s) of the flow alone per stripe |
| Member | restarts from the spin-up with the lines, runs to step 2058 (120 s more, 0.35 s steps) |
| Failures | `conductors.failure_break`, `ground_contact`, `tower_failure = release` |
| Capacity checks | from 610 s, every 0.01 s |
| DS4 | 2 or more towers lost (`cascade_towers = 2`) |

### Capacities

| Component | Capacity | Source |
|---|---|---|
| Tower footing uplift | 50 kN | **synthetic**, the deck's value from `make_case.py` |
| Tower footing compression | 150 kN | **synthetic**, the deck's value from `make_case.py` |
| Insulator string | 111,206 N (25,000 lb specified mechanical load) | Hubbell Hi\*Lite catalog CA08051E |
| Conductor | 1e9 N (not checked) | none |

No conductor rating is used: the conductors' peak tension is 20 to 23 kN at every stripe, mostly
their own weight, so a rating cannot separate the stripes in this deck. DS2 comes from the spans
that failed towers release.

Before the ensemble, a calibration with every capacity at 1e9 N (seed 1, 13 and 17 m/s) recorded
the largest footing uplift: 37.3 kN at 13 m/s and 76.1 kN at 17 m/s. 50 kN lies between them.
Insulators stayed below 12.7 kN. The transformers' allowables are in the deck but the
transformers are not in the demand table.

## Data and acknowledgements

No external data set is used or redistributed. The terrain, inflow profile and sounding are the
synthetic files of `Structures/Exec/PowerLines`, written by its `make_case.py`. The insulator rating
is taken from Hubbell Incorporated's Hi\*Lite catalog CA08051E, with thanks to Hubbell for making it
public; only the rated load is used.

## How to run

Build with `-DERF_ENABLE_STRUCTURES=ON` against MoorDyn-C (`Structures/setup_moordyn.sh`),
`-DERF_ENABLE_FFT=ON` and MPI. Then:

```bash
python3 run_case.py --exe path/to/flare --out powerlines_stripes --ranks 4
```

A spin-up takes about 25 min and a member about 5 min on 4 ranks, about 3.5 h for the 7 spin-ups
and 21 members. The driver waits for the 1-minute load to fall below 10 before each run
(`--load-start`, or `--no-load-gate`). Giving the same command again skips the runs that finished.
`--dry-run` writes the decks, the plans and the manifest and runs nothing.

## Expected results

Numbers from MoorDyn-C on 4 ranks per run (November 13 to 14, 2026). All 21 members finished
(return code 0). Every member restarted from its stripe's `chk01715`, and every run log says
the checkpoint holds no conductor state and the lines start at the restart.

Towers lost per member (seeds 1, 2, 3) and the largest peak ratio:

| Stripe (m/s) | 11 | 12 | 13 | 14 | 15 | 16 | 17 |
|---|---|---|---|---|---|---|---|
| Towers lost | 0, 0, 0 | 0, 0, 0 | 0, 2, 0 | 0, 3, 1 | 3, 8, 6 | 8, 8, 8 | 8, 8, 8 |
| Largest peak ratio | 0.44 to 0.65 | 0.58 to 0.85 | 0.75 to 1.60 | 0.92 to 1.62 | 2.35 to 2.73 | 2.44 to 2.56 | 2.43 to 2.75 |

Members reaching each state (of 3):

| Stripe (m/s) | 11 | 12 | 13 | 14 | 15 | 16 | 17 |
|---|---|---|---|---|---|---|---|
| DS1, DS2, DS3 | 0 | 0 | 1 | 2 | 3 | 3 | 3 |
| DS4 (2 or more towers) | 0 | 0 | 1 | 1 | 3 | 3 | 3 |

The fitted curves (1000 parametric-bootstrap sets each):

| State | Median theta (m/s) | Dispersion beta | theta, 95 % | beta, 95 % | Sets fitted / separated |
|---|---|---|---|---|---|
| **DS1 to DS3** | **13.48** | **0.054** | 12.8 to 14.1 | 0 to 0.084 | 400 / 600 |
| **DS4** | **13.77** | **0.062** | 13.2 to 14.5 | 0 to 0.096 | 504 / 496 |

A brute-force grid over (theta, beta) gives the same maxima (13.478, 0.054 and 13.772, 0.0625).

What the numbers show:

- **DS1 to DS3 are the same curve.** Only the footings have a capacity that the wind reaches, and a
  footing over its allowable fails its tower, which releases its spans. The first exceedance is
  therefore always a tower lost and a line released.
- **The losses saturate at 8 of the 10 towers.** Two towers are never lost, at any stripe.
- **The spread is wide for its few members.** Half or more of the simulated sets separate (no
  maximum-likelihood estimate), so beta's lower bound is 0. Three seeds per stripe resolve the
  median to about +/- 0.6 m/s and beta only loosely.
- **The tool notes that DS3 and DS4 cross at 11.78 m/s**, inside the stripes. Both curves are at
  0.6 % there; below it the DS4 curve, with its larger beta, lies above DS3. The stripes at 11 and
  12 m/s have no damage, so the crossing is in the fitted tails and not in the data.

Each member is one realisation: a falling span amplifies roundoff, so the times and the number of
towers lost can differ between builds and rank counts.

## CTest

This case runs for hours on 4 ranks and is not a CTest. `CascadeStripes_DamageStates`
(`Structures/Exec/CascadeStripes`) runs the same chain on a small deck in about 8 s.
