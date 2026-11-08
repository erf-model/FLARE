# Wind_Adjustment_Factor

## Purpose
Runs the wind adjustment factor that reduces the reference (20 ft) wind to
midflame height before it enters the spread calculation. A second, linear
"behaviorplus" formula (0.36 + 0.004 H) was removed in November 2026: it has no
published source (BehavePlus uses the Albini and Baughman 1979 form, Andrews 2012).

## Physics / Model Features Exercised
- Wind adjustment factor applied to the effective fire wind
- Andrews (2012) unsheltered formulation

## Expected Results
The midflame wind is below the reference wind (WAF 0.362 for the 1 ft bed of
fuel model 1), and the spread rate below an unadjusted run. An unrecognised
`waf_formula` aborts at start-up.

## Key Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `erf.fire.use_waf` | `true` | Enables the wind adjustment factor. |
| `erf.fire.waf_formula` | `"andrews"` | The formulation (the only one). |
| `erf.fire.fuel_model_id` | `1` | Fuel bed depth from this model drives the WAF. |

## References
- Andrews 2012, Modeling wind adjustment factor and midflame wind speed for Rothermel's surface fire spread model.
