# Flame_Temperature

## Purpose
Covers the McAlpine flame temperature, which no other deck uses, and the flame
tilt diagnostic. The `byram_radiant` default is used throughout the other fire
cases. A third method, `nelson_emissivity`, was removed in November 2026: it put
the fireline intensity [W/m] into the Stefan-Boltzmann law as a flux [W/m²]; its
deck became `inputs_fire_flame_tilt`.

## Physics / Model Features Exercised
- Flame temperature diagnostic from fireline intensity
- McAlpine heat-balance formulation (`inputs_fire_flame_temp_mcalpine`)
- Flame tilt diagnostic with the default temperature (`inputs_fire_flame_tilt`)

## Expected Results
Flame temperature should be finite and above ambient wherever the fire is
burning, and equal to ambient in unburned cells. Peak temperature should track
fireline intensity. The flame tilt is finite, in [0, 90) degrees, and zero in
unburned cells.

## Key Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `erf.fire.flame_temp_method` | `"mcalpine_heat"` / `"byram_radiant"` | Selects the formulation. |
| `erf.fire.compute_flame_tilt` | `true` (tilt case) | Also computes the flame tilt angle. |

## References
- McAlpine & Wakimoto 1991, The acceleration of fire from point source to equilibrium spread.
