# Broken wire

## Purpose

A straight conductor section in still air in which one span lets go of a support, run through ERF
with the conductor model on MoorDyn-C. The case answers three questions about the line failure
path (`conductors.failure_break`, `conductors.tower_failure`, `conductors.capacity_check_interval`,
`conductors.capacity_check_start`):

1. **Residual tension.** After a span lets go, the insulator strings lean towards the spans that
   stay and the tension there falls. Does the model come to the state at rest that the balance of
   forces gives? The reference is `Structures/Tools/broken_wire/broken_wire.py`, which shares no
   code with ERF or MoorDyn.
2. **Overshoot.** How much larger than the load at rest is the largest load on the strings, the
   dead end and the towers' legs just after the release?
3. **Cascade.** On a line of equal spans, how many towers fail as their legs' bearing capacity is
   varied, with the capacities checked once per step and every 0.01 s?

The equations of the state at rest are in "State at rest after a release" in
`Docs/sphinx_doc/theory/Conductors.rst`. No design standard's broken-wire load is reproduced; see
"Data and acknowledgements".

## Files

- `make_case.py` writes the decks of one section: `inputs` runs to one step after the release and
  writes a checkpoint, `inputs_after` restarts from it and runs on, and `plan.json` is the same
  section for `broken_wire.py`. `python3 make_case.py --help` lists the options (spans, cross-arm
  height, string length, which end of span 1 lets go, the towers' leg capacity, the check interval).
- `run_study.py` writes and runs the cases of the three studies and prints their tables
  (`residual.csv`, `overshoot.csv`, `cascade.csv`).

Why two decks: a line has one rated load per kind of component, and the load that passes along the
line just after a release is two to four times the load at rest. A capacity set a few percent under
the first component's load therefore releases its neighbours within a second. The second deck
replaces that capacity, so that what happens after the first release is set by the towers' leg
capacity alone, or by nothing.

The flow is at rest on a coarse mesh (50 m cells) and the lines put nothing into it
(`conductors.prescribed_velocity = 0 0 0`); the flow solver only carries the time step.

## Data and acknowledgements

No data set is used or redistributed. The conductor is the 795 kcmil ACSR of the model's other
tests, with the same diameter, mass and axial stiffness as there. The strings' mass (60 kg), the
towers' weight (60 kN) and base width (6 m) are round numbers chosen for the case, not a supplier's.

The residual static load of a broken wire is tabulated in design guides (ASCE Manual of Practice
74; EPRI's work on cascades). We had no access to them, so none of their numbers or load factors
is used or compared here. The comparison is with the balance of forces alone.

## How to run

Build with `-DERF_ENABLE_STRUCTURES=ON` against MoorDyn-C (`Structures/setup_moordyn.sh`),
`-DERF_ENABLE_FFT=ON` and MPI. One case by hand:

```bash
python3 make_case.py --out run --release span_b --time 600 --release_time 10
```

```bash
cd run && mpiexec -n 2 erf_exec inputs && mpiexec -n 2 erf_exec inputs_after
```

```bash
python3 ../../Tools/broken_wire/broken_wire.py run/plan.json
```

All three studies (25 cases, about four minutes on two ranks):

```bash
python3 run_study.py --exe path/to/erf_exec --out study
```

## Expected results

Numbers from MoorDyn-C 2.7.1 on two ranks, 20 links per span, coupling step 0.25 s. Each run is
one realisation: a falling conductor amplifies roundoff, so the largest loads and the extent of a
cascade near a threshold differ between builds and rank counts. The loads at rest do not.

**1. Residual tension.** Span 1 let go at tower 1 (or at its dead end), 600 s, the towers' pulls
averaged over the last 240 s (the section keeps swaying by tens of newtons for minutes).

| Section | Span 2 keeps (model) | (force balance) | Tower 1 along the line: model (N) | force balance (N) |
|---|---|---|---|---|
| 3 x 300 m, 2.5 m strings | 70.7 % | 70.6 % | 9344.9 | 9342.8 |
| 3 x 300 m, 4 m strings | 62.5 % | 62.4 % | 8260.0 | 8258.0 |
| 300, 500, 100 m, 2.5 m strings | 87.0 % | 87.0 % | 11748.3 | 11746.9 |
| 3 x 300 m, cut at the dead end | 71.2 % | 71.2 % | 9413.5 | 9411.6 |

Every tower pull, along the line and vertical, is within 0.02 % of tower 1's largest pull.

**2. Overshoot.** Span 1 let go at tower 1, spans of 300 m either side of a middle span of 150 to
600 m, cross-arms 100 m up, capacities sampled every 0.01 s. Largest load after the release over
the load at rest after it:

| Middle span / outer spans | String 1 | String 2 | Dead end b | Tower 1 leg | Tower 2 leg |
|---|---|---|---|---|---|
| 0.5 | 4.2 | 3.9 | 3.4 | 3.4 | 3.8 |
| 1 | 4.4 | 3.5 | 3.6 | 4.0 | 2.5 |
| 1.5 | 3.9 | 3.6 | 4.2 | 3.6 | 2.6 |
| 2 | 4.1 | 3.3 | 4.4 | 3.9 | 2.4 |

Over the load before the release the string next to the cut sees 3.2 to 4.7 times, the dead end
2.5 to 4.2 times. The towers' legs see far more, 13 to 26 times at tower 1 and 2.7 to 10 times at
tower 2, because the pull along the line acts 100 m above a base 6 m wide; before the release the
legs carry the weight alone. `overshoot.csv` has every load in newtons.

**3. Cascade.** Eight spans of 300 m (seven towers), cross-arms 60 m up, span 1 let go at tower 1
at 10 s. The legs' bearing capacity is a multiple of their load before the release (16.3 kN).

| Capacity / load before | Towers failed of 7, one check per step | Checks every 0.01 s |
|---|---|---|
| 1.5 | 7 | 7 |
| 2 | 7 | 7 |
| 3 | 7 | 7 |
| 5 | 7 | 7 |
| 8 | 2 | 3 |
| 12 | 0 | 1 |
| 16 | 0 | 0 |
| 24 | 0 | 0 |

Up to five times the load before the release every tower fails, the last one 2 to 9 s after the release. At 8 and
12 times the checks inside the step find a failure that one check per step misses. From 16 times
on the line stands.

Limits of these numbers: the towers are rigid and stay standing after they let go of their lines,
the conductor has no structural damping beyond MoorDyn's, and a tower that has failed releases its
strings at the next coupling step. A real tower that yields takes up some of the load first.

## CTest

`BrokenWire_ThroughERF` (labels `regression`, `conductors`, `failure`; MoorDyn-C builds only) runs
study 1 and fails if a tower's pull differs from the force balance by more than 0.2 %:

```bash
python3 run_study.py --exe path/to/erf_exec --out study --study residual --check 0.2
```

The unit tests `BrokenWire.*` make the same comparison on the conductor line alone, to 0.03 %, and
`BrokenWire_Cable` to `BrokenWire_Refused` test the force balance itself.
