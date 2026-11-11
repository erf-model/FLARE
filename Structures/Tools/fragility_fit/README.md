# Fragility fit

## Purpose

`fragility_fit.py` fits wind fragility curves to the table an ensemble of conductor runs gathers
(`demand_table_all.csv` of `Structures/Tools/fragility_ensemble`). For every component, and for
every group of components named in the plan, it estimates the median wind speed and the
dispersion of the lognormal fragility function

    P(exceedance | V = x) = Phi(ln(x / theta) / beta)

by maximum likelihood on the members' outcomes: the multiple stripes analysis of Baker (2015).
It is Python with the standard library only, reads one CSV file and runs no model. See
"Fragility curves from an ensemble" in `Docs/sphinx_doc/theory/Conductors.rst` for the
equations.

## Files

- `fragility_fit.py`: the tool.
- `tests/test_fragility_fit.py`: the tests, with two stored demand tables (`tests/ensemble.csv`,
  `tests/two_stripes.csv`).
- `FragilityFitTests.cmake`: the CTest registrations, included by `Tests/CTestList.cmake`.

## Data and acknowledgements

No data set is redistributed; the tests' tables are synthetic. The method and the numbers the
tests are pinned to come from these publications, whose authors are thanked:

- Baker, J. W. (2015). Efficient analytical fragility function fitting using dynamic structural
  analysis. *Earthquake Spectra* 31(1), 579–599. doi:10.1193/021113EQS025M. The fit (its Eqs. 1
  and 8 to 11), the procedure for the spread of the estimates, and the printed standard
  deviations the tests reproduce. Read through the institution's library access; the author
  also posts the paper and a code supplement (https://purl.stanford.edu/sw589ts9300), which was
  not used.
- Panteli, M., Pickering, C., Wilkinson, S., Dawson, R., and Mancarella, P. (2017). Power system
  resilience to extreme weather: fragility modeling, probabilistic impact assessment, and
  adaptation measures. *IEEE Transactions on Power Systems* 32(5), 3747–3757.
  doi:10.1109/TPWRS.2016.2641463 (Creative Commons Attribution 3.0). The series rule for a group
  of components (its Eq. 4).
- Liu, C., and Yan, Z. (2022). Fragility analysis of wind-induced collapse of a transmission tower
  considering corrosion. *Buildings* 12, 1500. doi:10.3390/buildings12101500 (open access, CC BY
  4.0). Cited for an example of a capacity distribution: steel yield strength, lognormal,
  coefficient of variation 0.07, the standard strength being the 0.05 quantile (its Table 2).
- Mazaheri, P., and Alipour, A. (2025). Fragility functions for weather events to enhance
  resiliency of transmission lines: from experiments to statistics. *Electrical Transmission and
  Substation Structures 2025*, 443–451, ASCE. Cited for the distributions it uses for conductor
  and insulator strength (lognormal, coefficient of variation 0.1; the ratio of mean to nominal is
  not stated). Read from a copy the user supplied.
- Zhu, C., Yang, Q., Huang, G., Zhang, X., and Wang, D. (2024). Fragility analysis and wind
  directionality-based failure probability evaluation of transmission tower under strong winds.
  *Journal of Wind Engineering and Industrial Aerodynamics* 246, 105668.
  doi:10.1016/j.jweia.2024.105668. Cited for an example of a capacity distribution (steel yield
  strength, lognormal, coefficient of variation 0.07, its Table 3) and for the treatment of
  directions this tool does not do. Read through the institution's library access.

## The plan

A JSON object. Paths are relative to the plan file.

| Key | Default | Range | Effect |
| --- | --- | --- | --- |
| `demand_table` | none (required) | a readable `demand_table_all.csv` | the gathered table: one row per member and component |
| `output_dir` | `fragility` | a directory name | where the two output tables go |
| `intensity` | `max_speed` | `max_speed`, `peak_speed` or `stripe` | the wind a member contributes for a component: the largest averaged wind its point recorded, the wind at the component's peak demand, or the member's nominal stripe |
| `by_direction` | `false` | `true` or `false` | one curve per nominal direction instead of one over all directions |
| `capacity` | `{}` (nominal capacities) | `{kind: {"cov": [0.001, 1], "mean_factor": [0.1, 10]}}`; kinds `conductor`, `insulator`, `hardware`, `footing_uplift`, `footing_compression`, `member`; `mean_factor` is 1 when left out | the kind's capacity is lognormal with mean `mean_factor` x nominal and coefficient of variation `cov`; a member's outcome is then the probability that the capacity lies below its peak demand |
| `groups` | `{}` | `{name: {"components": regular expression, "kinds": [..]}}`, at least one of the two; a name that is no component's and holds no comma; the group must match a component | a curve for the group: a member exceeds when any of the group's components does; the expression is searched in the name (`^L1` for names that start with L1), every member must hold every component of the group, and the curve is fitted against the nominal stripe whatever `intensity` is |
| `resamples` | 1000 | integer in [0, 100000] | simulated sets for the spread of the estimates; 0 leaves the spread out; the spread is given only for curves whose outcomes are all 0 or 1 (no capacity distribution) |
| `seed` | 1 | integer in [0, 2147483647] | makes the simulated sets repeatable; each curve has its own sequence from the seed and its name |
| `confidence` | 0.9 | [0.5, 1) | probability of the central interval reported for theta and beta |

There are no default capacity statistics: `cov` and `mean_factor` are the user's, with their
source. None of the publications above gives measured statistics of conductor, insulator,
hardware or footing capacity against a rated value. A key out of its range stops the tool with a message that names it (return code 2), and
so does a table that lacks a column, has a short row, holds a value that is not a number, a wind
that is not positive, an unknown kind, one member's component twice, or a member whose rows
disagree in stripe, direction or seed. A row without a peak ratio (a component the run never
recorded) is left out and counted in a note; so is a member without the chosen wind, by name.

## How to run

    python3 fragility_fit.py --plan plan.json

with, for example,

    {
     "demand_table": "ensemble/demand_table_all.csv",
     "capacity": {"member": {"cov": 0.07}},
     "groups": {"line_L1": {"components": "^L1"}, "towers": {"kinds": ["member", "footing_uplift", "footing_compression"]}}
    }

## Outputs

In `output_dir`:

- `fragility_curves.csv`: one row per component and per group (per direction with
  `by_direction`): `curve`, `kind` (`group` for a group), `direction`, `intensity` (the wind the
  curve was fitted against), `members`, `exceedances`
  (the sum of the outcomes), `status`, `theta` (m/s), `beta`, the interval `theta_low`,
  `theta_high`, `beta_low`, `beta_high`, the standard deviations `theta_sd`, `beta_sd`,
  `resamples_used`, `bracket_low`, `bracket_high` and `log_likelihood`.
- `fragility_points.csv`: one row per curve and stripe, for a plot: `members`, `mean_intensity`,
  `exceedances`, `fraction` and the `fitted` probability at the mean intensity.

`status` is `fitted`, or the reason there is no fit:

| Status | Meaning |
| --- | --- |
| `no_exceedance` | the outcomes add up to less than 1e-6 exceedances; `bracket_low` is the largest wind of the data |
| `all_exceed` | all but less than 1e-6 members exceed; `bracket_high` is the smallest wind of the data |
| `separated` | no wind with a member that does not exceed lies above a wind with one that does, so the likelihood has no maximum; `bracket_low` is the largest wind at which no member exceeds, `bracket_high` the smallest at which all do |
| `one_intensity` | all members recorded one wind speed |
| `not_increasing` | the best curve does not rise with the wind |
| `flat` | the best curve rises by less than 0.05 over the winds of the data |
| `not_converged` | the search did not settle |
| `no_intensity` | no member recorded the chosen wind |

The brackets say where the data change; they are not bounds on the median.

`separated` is the expected result for a nominal capacity in a set of RANS runs without gust
seeds: every member of a stripe has the same flow, so the demand crosses the capacity once. A
dispersion needs gust seeds, a capacity distribution or resolved turbulence.

## Expected results

One ensemble was run with the ensemble driver and fitted (2026-11-10, by hand): nine stripes
from 13 to 20 m/s on `Tests/test_files/Conductors_TwoLevel` (one 300 m span across a sheared
crosswind, MoorDyn stub, 2 ranks), with the rated breaking strength 13890.45 N and the hardware
rated load 13892 N of the `Conductors_Intensity_TwoLevel` test. The span's peak tension rises
from 13 653 N at 13 m/s to 14 491 N at 20 m/s, a peak ratio of 0.983 to 1.043: the weight of the
conductor carries most of the tension.

With nominal capacities (`{"demand_table": "ensemble/demand_table_all.csv"}`) every curve is
`separated`: no member exceeds up to 15.50 m/s and all do from 16.00 m/s on, for the span and for
both dead ends. This is the RANS result described above.

With a lognormal capacity of the conductor, `"capacity": {"conductor": {"cov": 0.1}}` (the
coefficient of variation Mazaheri and Alipour 2025 use; a mean equal to the nominal value is an
assumption made for this example), the span's curve is fitted:

| Curve | Status | theta (m/s) | beta |
| --- | --- | --- | --- |
| `S1` (conductor, recorded wind) | fitted | 15.00 | 0.748 |
| `S1_end_a`, `S1_end_b` (hardware, nominal) | separated, 15.50 to 16.00 m/s | | |

The members' outcomes are 0.451 at 13 m/s, 0.506 at 15.5 m/s and 0.682 at 20 m/s, and the fitted
curve gives 0.424, 0.518 and 0.650 there. The dispersion is large because the tension, and with
it the probability, changes little over these wind speeds; no spread is reported for a curve
fitted to probabilities.


## CTest

`FragilityFit_Published`, `FragilityFit_ByHand`, `FragilityFit_Table` and `FragilityFit_Refused`
(labels `unit`, `conductors`, `fragility`) run `tests/test_fragility_fit.py`; they need Python 3
and run no ERF executable.

    python3 tests/test_fragility_fit.py
