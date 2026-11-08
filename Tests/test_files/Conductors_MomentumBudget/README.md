# Conductors_MomentumBudget

## Purpose

Checks that the drag of a conductor line reaches the flow whole on one level and on two, through
ERF's structure coupling (`erf.structures.model = conductors`). A periodic box of uniform air
(density 1, 300 K, no gravity, Coriolis force or diffusion, slip walls) blows at 10 m/s across one
line. Nothing but the line's drag changes the air's y-momentum, so its change over the run must
equal the impulse the model reports putting into the air. With a finer level over the drag of
lines on a coarser level, averaging down would remove the drag (the review point on
erf-model/ERF#5009); ERF's structure coupling stops such a run, which the abort tests check.

## Files

- `Conductors_MomentumBudget.i`: the deck, 600 x 400 x 150 m on 24 x 16 x 6 cells, ten steps of
  0.5 s. The line is 6 m thick (a test value) so that its impulse is about a thousandth of the air's
  momentum, which `fvolumesum`'s digits resolve.

## Data and acknowledgements

No external data.

## How to run

Build FLARE with `-DERF_ENABLE_STRUCTURES=ON` (the stub or MoorDyn-C) and MPI, then

    mpiexec -n 2 erf_exec Conductors_MomentumBudget.i
    amrex_fvolumesum -v y_velocity plt00000
    amrex_fvolumesum -v y_velocity plt00010

and sum `source_y` times the step over the rows of `conductors/total_load.dat`. The two-level
variants add the options listed in `Structures/Tests/ERFStructuresTests.cmake`.

## Expected results

The change of the volume integral of `y_velocity` (density 1) equals the sum of `source_y` times
the step:

| Run | Levels | Stub | MoorDyn-C |
|---|---|---|---|
| `Conductors_MomentumBudget` | one | 5e-7 | 3.7e-7 |
| `Conductors_MomentumBudget_FinestLevel` | the line on level 1, a slab over the whole box up to 150 m | 1.3e-5 | 1.0e-3 |
| `Conductors_MomentumBudget_FinestLevelPatch` | the line on level 1, a patch over the whole width in y | 1.9e-6 | 5.3e-4 |
| `Conductors_MomentumBudget_FinerLevelBelow` | the line on level 0, level 1 below the drag's reach | 1.1e-6 | 9.8e-5 |

The stub's line holds its quasi-static shape; MoorDyn-C's swings and bobs, and the velocities it
sends through the coarse/fine faces carry momentum fluxes the levels do not synchronise, so its
two-level budgets close less tightly (tolerances 3e-3, 2e-3 and 5e-4).

With ERF's check switched off and level 1 over the whole box, a line on level 0 lost all of its
drag. ERF does not synchronise momentum fluxes at coarse/fine faces: with level 1 a patch whose
faces normal to y were 75 to 150 m from the line, the y-momentum crossing them (mostly the
pressure term) differed between the levels and the budget missed by 23 % to 111 %; with faces
only normal to x it closed to 2e-6. The finer levels here have no faces normal to y. Every run
also checks the first step's source against the drag of the line at rest, 1/2 rho Cd D L U^2 =
9e4 N, since a wrong drag would move the applied and the reported source together: to 1e-6 with
the stub, whose line holds its quasi-static shape, and between half of it and all of it with
MoorDyn-C, whose line starts to move downwind within the step (6.18e4 N measured).

## CTest

`Conductors_MomentumBudget`, `Conductors_MomentumBudget_FinestLevel`, `Conductors_MomentumBudget_FinestLevelPatch`,
`Conductors_MomentumBudget_FinerLevelBelow`, `Conductors_Abort_FinerLevelCoversLines`,
`Conductors_Abort_FinerLevelCoversDragTail` and `Conductors_OneWayUnderFinerLevel`
(`Structures/Tests/ERFStructuresTests.cmake`, label `conductors`).
