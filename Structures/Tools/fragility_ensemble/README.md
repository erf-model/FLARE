# Fragility ensemble driver

## Purpose

A fragility curve for a line component needs many pairs of wind intensity and peak demand. One
run with `conductors.intensity_measure = wind` gives one pair per component
(`conductors/demand_table.csv`, `conductors/intensity.csv`). `fragility_ensemble.py` writes, runs
and gathers a set of such runs from one base deck and a plan: wind-speed stripes x wind directions
x gust seeds. It is Python with the standard library only and changes no model. Fitting curves to
the gathered tables is a later step and is not done here. See "An ensemble of runs for fragility
curves" in `Docs/sphinx_doc/theory/Conductors.rst`.

## Files

- `fragility_ensemble.py`: the driver (subcommands `run` and `collect`).
- `examples/two_level/plan.json`: a plan on the smallest conductor deck,
  `Tests/test_files/Conductors_TwoLevel` (one span, periodic box, sheared crosswind).
- `tests/test_fragility_ensemble.py`: the tests, with their stored base decks (`tests/base_*`)
  and stored member tables (`tests/collect`).
- `FragilityEnsembleTests.cmake`: the CTest registrations, included by `Tests/CTestList.cmake`.

## Data and acknowledgements

No external data. The example's base deck and the tests' decks and tables are synthetic.

## The plan

A JSON object. Paths are relative to the plan file.

| Key | Default | Range | Effect |
| --- | --- | --- | --- |
| `base_deck` | none (required) | a readable deck that sets `erf.structures.model = conductors` | the deck every member's deck is written from; its `FILE` includes are written into the member's deck |
| `output_dir` | `ensemble` | a directory name | where `manifest.json`, `members/` and the gathered tables go |
| `speeds` | none (required) | each in (0, 150] m/s, no repeats; a number or `{"speed": .., "inputs": {..}}`; the scaled profile must stay within 150 m/s at every height | the stripes: the wind speed at the reference height; a stripe's `inputs` are deck inputs for its members only |
| `directions` | the base deck's direction | each in (-180, 180] degrees, no repeats | the direction the wind blows towards, from +x towards +y (the convention of `demand_table.csv`) |
| `seeds` | none (one member per stripe and direction) | integers in [0, 2147483647], no repeats; needs `conductors.gust_type = random` | `conductors.gust_seed` of each member |
| `reference_height` | 10 | (0, 500] m, within the wind profile's heights and not between a first row without wind and the next row; not with `conductors.prescribed_velocity` | the height at which the base wind's speed and direction are read (linear interpolation between the file's rows, as ERF reads it); independent of `conductors.intensity_height` |
| `side_type` | the type of the base deck's side faces | `SlipWall`, `NoSlipWall` or `Symmetry`; not with `conductors.prescribed_velocity` | the type of a pair of lateral faces the member's wind blows along |
| `inputs` | `{}` | deck keys with a number, `true`/`false`, a string (without `#` or `=`) or a list; not an input the members' wind is derived from (see below) | deck inputs set for every member (a capacity, `conductors.gust_type`, `conductors.intensity_averaging_time`) |
| `support_files` | `[]` | files in the base deck's directory or below it, no repeats; must hold every file the deck names that stands beside it | files copied into every member's directory (a terrain file, a member design file) |
| `spinup` | none (each member runs from the initial state) | `{"max_step": n}`, n an integer in [1, 100000000] below every member's `max_step`; not with `conductors.prescribed_velocity`; the members' inputs may not set `erf.check_int` or `erf.check_per` | the flow of each stripe and direction runs n steps once, without the structure model, and its members restart from that checkpoint with the lines (see below) |

The plan gives at most 100 000 members. These inputs cannot be given in `inputs` or in a stripe's
`inputs`: `conductors.intensity_measure`, `conductors.prescribed_velocity`,
`conductors.diagnostics_dir`, `erf.init_type`, `erf.input_sounding_file`, `geometry.is_periodic`,
`geometry.prob_lo`, `geometry.prob_hi`, `amr.n_cell`, `amr.restart`, `erf.restart`,
`erf.input_bndry_planes`, `erf.abl_driver_type`, `erf.sounding_type`, with a sounding any input
of the `xlo`, `xhi`, `ylo` and `yhi` faces, and with `seeds` `conductors.gust_seed`.

How a member's wind is set:

- **Sounding** (`erf.init_type = input_sounding`): the sounding and, when a lateral face is an
  `Inflow` with a `dirichlet_file`, that profile are scaled by stripe / base speed and turned by
  direction - base direction at every height; the base speed and direction are read at
  `reference_height` from the inflow profile when there is one, otherwise from the sounding. In a
  direction that is not periodic the face the wind enters becomes the `Inflow`, the face it leaves
  the `Outflow`, and a pair of faces the wind blows along become side faces; an inflow face's `KE`
  is multiplied by (stripe / base speed)^2 (ERF uses `KE` only together with `<face>.density`;
  the driver prints a note when the density is missing).
- **Prescribed wind** (`conductors.prescribed_velocity`): that velocity is scaled and turned in the
  deck; the flow and its boundaries are not changed.

Every member also gets `conductors.intensity_measure = wind`. The base deck (or the plan's
`inputs`) must give at least one capacity, which the intensity record needs; the driver does not
check this, and without one every member stops in ERF and is reported as failed. The time step,
the run's length and the averaging time are not changed: choose the base deck's step for the
fastest member, or give fast stripes their own `erf.fixed_dt` and `max_step`. The manifest holds
each member's largest wind (`max_wind`) and its largest (|u|/dx + |v|/dy) dt on level 0
(`courant`), and the driver prints a note for a member above 1. A base deck that sets or drives
its wind another way (`<face>.velocity`, `<face>.inflow_profile`, boundary planes, a restart,
another `erf.init_type`, `erf.abl_pressure_grad`, `erf.abl_geo_wind`), has a lateral face that is
not `Inflow`, `Outflow`, `SlipWall`, `NoSlipWall`, `Symmetry` or periodic, or writes its
diagnostics outside the run directory is refused with the input's name. Other inputs that hold a
wind speed (`conductors.asce74_wind`, the event gust's inputs) are not scaled. The driver reads
one `key = value` per line; a `FILE` include is looked for beside the file that names it.

### Shared spin-up

The seed sets only `conductors.gust_seed`, which the flow does not see. Members of one stripe and
direction therefore share one flow. With `spinup`, the driver writes one spin-up run per distinct
flow under `spinup/<name>` (`u001_s10_d0`: the flow's number, stripe and direction). Its deck is
the member's deck without any `conductors.` or `erf.structures.` input, with `max_step = n` and
`erf.check_int = n`. Members are grouped by a SHA-256 of that deck and of their wind files, so a
stripe's own `inputs` that change the flow give it its own spin-up. Each member's deck gets
`amr.restart = ../../spinup/<name>/<check_file>NNNNN` (`erf.check_file`, default `chk`, and n
written with five digits) and runs from step n to its own `max_step`.

The spin-up has no structure model because a restart keeps the conductor model's peak demands,
first exceedances, intensity records and failure logs from the checkpoint. A checkpoint without
them makes the lines afresh at the restart time, so every member starts its lines, records and
logs on the spun-up flow.

`run` runs each spin-up once, before its members, and skips one whose status is `done` and whose
checkpoint `Header` is there. A spin-up that ends with return code 0 but without the checkpoint
(its `stop_time` ended it before step n) is marked failed. When a spin-up fails, its members are
marked failed without a run, with the detail `not run: its spin-up spinup/<name> failed (return
code c)`; giving `run` again tries the spin-up again. The manifest lists the spin-ups
(`spinups`) and each member's (`spinup`).

## How to run

    python3 fragility_ensemble.py run --plan plan.json --dry-run
    python3 fragility_ensemble.py run --plan plan.json --exe <build>/Exec/erf_exec --ranks 2
    python3 fragility_ensemble.py collect --plan plan.json

`run` writes the member directories and `manifest.json`, runs the members one after another and
gathers their tables. `--dry-run` stops after the directories and the manifest. Giving `run`
again skips the members that finished (return code 0) and runs the others from their start. The
manifest holds a SHA-256 of the base deck, the wind files and the support files: a directory that
was written from another plan, or from files that have changed since, is refused, and so is a
directory another driver is running (`run.lock`). When the driver is interrupted (Ctrl-C, SIGTERM,
SIGHUP) it ends the run's processes with it; that member has no status and runs again next time.
`collect` gathers the tables again without running.

| Option | Default | Range | Effect |
| --- | --- | --- | --- |
| `--plan` | none (required) | a JSON plan | the plan |
| `--dry-run` | off | | write the decks and the manifest, run nothing |
| `--exe` | none | an executable file; needed unless `--dry-run` (checked when given) | the ERF executable, built with `-DERF_ENABLE_STRUCTURES=ON` |
| `--mpiexec` | `mpiexec` | a launcher on the PATH, or `''` for none (then `--ranks 1`) | each member runs as `<mpiexec> -n <ranks> <mpiexec-args> <exe> <deck>` |
| `--mpiexec-args` | none | one string, written `--mpiexec-args='--oversubscribe'` | further launcher arguments |
| `--ranks` | 4 | [1, `--max-ranks`] | MPI ranks of each run |
| `--max-ranks` | 4 | [1, 1048576] | the most ranks `--ranks` may ask for |
| `--load-start` | 6 | (0, 100000], not above `--load-pause` | a run starts when the 1-minute load is below this |
| `--load-pause` | 12 | (0, 100000] | a run's processes are stopped (SIGSTOP) when the load passes this |
| `--load-resume` | 7 | (0, 100000], below `--load-pause` | a stopped run continues (SIGCONT) when the load falls below this |
| `--load-poll` | 20 | (0, 3600] s | time between looks at the load |
| `--no-load-gate` | off | | run without looking at the load |

An option out of its range stops the driver with a message that names it (return code 2).

The load gate reads the 1-minute load of the machine the driver runs on. It finds the run's
processes as the launcher's descendants with `ps` and stops all of them, also those that start
during a pause. Ranks on other nodes, or started by a launcher that does not keep them as its
descendants (a scheduler's `srun`), are not paused: use `--no-load-gate` there. When `ps` is not
available the driver says so and runs on without pausing.

## Outputs

In `output_dir`:

- `manifest.json`: the base deck, the base wind (speed, direction, reference height and file),
  a SHA-256 of every source file, and for every member its name, stripe, direction, seed,
  directory, scale factor, rotation, largest wind, Courant number and the inputs in which its deck
  differs from the base deck (`null` marks an input removed).
- `members/<member>/`: the member's deck, wind files and support files; after its run `run.log`,
  `member_status.json` (status, return code, ranks, wall time) and the run's own output.
- `demand_table_all.csv`, `intensity_all.csv`: the rows of every complete member's
  `demand_table.csv` and `intensity.csv` behind the key columns `member`, `stripe`,
  `nominal_direction` and `seed` (the component or point is the table's own column).
- `members_status.csv`: one row per member: `complete`, `not_run`, `failed`, `no_demand_table`,
  `no_intensity_table`, `bad_demand_table` or `bad_intensity_table`, with the return code, the
  rows gathered and, for a failed run, the abort message of its log. A member goes into both
  gathered tables or into neither. `run` and `collect` return 1 when a member is not complete.
- `breaks_all.csv`, `tower_failures_all.csv`, `failures_all.csv`: the rows of every complete
  member's `breaks.dat` (`time line released component tension`), `tower_failures.dat`
  (`time tower component releases`) and `failures.dat`
  (`time kind component demand capacity ratio`) behind the same key columns. A run writes these
  logs only with the options that make the events (`conductors.failure_break`,
  `conductors.tower_failure = release`); a log with another header or a row with the wrong number
  of fields gives the member the status `bad_breaks_log`, `bad_tower_failures_log` or
  `bad_failures_log`, and it is not gathered.
- `damage_table_all.csv`: one row per complete member behind the key columns:

  | Column | Value |
  | --- | --- |
  | `max_peak_ratio` | the largest `peak_ratio` of its demand table |
  | `components_exceeded` | the number of components with `peak_ratio` > 1 |
  | `failures_log`, `breaks_log`, `tower_log` | `yes` when the run wrote that log, `no` when it did not |
  | `failures` | rows of `failures.dat` |
  | `lines_released` | distinct lines in `breaks.dat` |
  | `releases` | rows of `breaks.dat` |
  | `towers_lost` | distinct towers in `tower_failures.dat` |

  A count is empty when its log was not written. A tower that fails lets go of its spans, and
  those releases are also written to `breaks.dat`, so a member with a tower lost also has a line
  released. `fragility_fit.py` reads this table for damage-state curves.

## Expected results

`examples/two_level/plan.json` runs `Tests/test_files/Conductors_TwoLevel` (one 300 m span along
x, periodic box, crosswind v = 15 + 0.1 z m/s, so 16 m/s towards 90 degrees at 10 m) at two
stripes and two directions, with the span's and dead ends' capacities of the
`Conductors_Intensity_TwoLevel` test. Copy the plan to a work directory, set `base_deck` to the
deck's path, and

    python3 fragility_ensemble.py run --plan plan.json --exe <build>/Exec/erf_exec --ranks 2

On 2 ranks with the bundled MoorDyn stub (2026-11-10, one run by hand; each member takes about a
second) all four members were complete, and `demand_table_all.csv` held 12 rows:

| Member | Wind recorded under the span's middle (m/s, degrees) | Span's peak tension (N) | Peak ratio | Exceeded |
| --- | --- | --- | --- | --- |
| `m001_s16_d90` | 16.0008, 90.0 | 13892.6 | 1.0002 | yes |
| `m002_s16_d0` | 16.0002, 0.0 | 13464.8 | 0.9694 | no |
| `m003_s24_d90` | 24.0009, 90.0 | 15527.5 | 1.1179 | yes |
| `m004_s24_d0` | 24.0004, 0.0 | 13464.8 | 0.9694 | no |

The recorded wind is the run's largest 1-second average (`max_speed`, `max_direction`). It equals
the stripe within 0.01 % and the direction within 0.001 degrees at the three intensity points, so
the scaled and turned sounding reaches the lines as planned. A wind along the span adds no
tension. Giving the command again skipped all four members.

A second trial checked an oblique wind through two inflow faces and the gust seeds on
`Tests/test_files/Conductors_RandomGusts` (two lines on lattice towers over hills, k-equation
RANS, inflow through x-low): 18 m/s at 30 m, directions 0 and 30 degrees, seeds 7 and 8, 10 steps.
All four members ran. At the 15 intensity points the wind blew towards -6.9 to 3.1 degrees in the
0 degree members and towards 28.9 to 34.9 degrees in the 30 degree members (the hills turn it).
Members that differ only in the seed recorded the same wind and different peak tensions (for the
first span of L1 at 30 degrees, 23 992 N and 35 285 N). This is a run of 9 s from the initial
state, long enough to show that the rotated inflow runs and reaches the lines; it is not a
spun-up flow.


## CTest

`FragilityEnsemble_DryRun`, `FragilityEnsemble_Refused`, `FragilityEnsemble_Run`,
`FragilityEnsemble_Collect` and `FragilityEnsemble_Spinup` (labels `unit`, `conductors`, `fragility`) run
`tests/test_fragility_ensemble.py`; they need Python 3 and run no ERF executable. The expected
decks, wind files, manifest and gathered tables are written out by hand in the test file.

    python3 tests/test_fragility_ensemble.py
