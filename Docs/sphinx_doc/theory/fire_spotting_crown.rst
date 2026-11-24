.. role:: cpp(code)
   :language: c++

.. _sec:FireSpottingCrown:

Ember Spotting and Crown Fire
=============================

Both extensions are off by default and add fields, plotfile variables and
checkpoint state only when enabled.

Ember spotting
--------------

:cpp:`erf.fire.spotting.enable` turns on a stochastic firebrand model after
Albini (1983). Every :cpp:`erf.fire.spotting.spotting_interval` fire
steps (one fire step per atmospheric step; default every step), with the
same random draws on every rank whatever the decomposition:

1. The launch set is identified. With :cpp:`erf.fire.spotting.launch_from =
   "burned"` (default) it is every burned cell (:math:`\phi < 0`) with a
   positive rate of spread, so the number of brands grows with the burned
   area and the launch probability has to be retuned for each domain. With
   ``"front"`` it is only the burned cells that have unburned fuel within
   :cpp:`erf.fire.spotting.front_band` fire cells (default 2), the fireline
   Albini's model describes, so the count follows the fireline length.
2. The Byram fireline intensity :math:`I_B` [kW/m] of each front cell is
   formed from the rate of spread, the fuel load and the heat content
   (:math:`h w_0 R` on every burning cell until its load is 1 % of the
   initial, about :math:`4.6\,\tau` behind the front; the launch band of
   ``"burned"`` and the :cpp:`P_base` calibration moved with this in
   October 2026).
3. Cells below :cpp:`erf.fire.spotting.I_B_min` (default 100 kW/m) launch
   nothing. Otherwise the lofting height is

   .. math::

      H_z = 0.173\, E^{1/2}, \qquad E = I_B\, A\, U_{10}^{B}

   (Albini 1983, INT-309, equations 2, 7 and 8), with :math:`H_z` in m,
   :math:`E` the strength of the line thermals in kJ/m, :math:`I_B` the
   fireline intensity in kW/m and :math:`U_{10}` the mean 10 m wind in m/s:
   the reference wind brought to 10 m on the same log profile as the drift
   in step 5. :math:`A` [s] and :math:`B` come from INT-309's Table 4 for the
   launch cell's Anderson model:

   ========================  =======  =======
   Anderson model            A [s]    B
   ========================  =======  =======
   1 short grass             162      -1.456
   2 grassy understory       166      -1.282
   3 tall grass              129      -1.238
   9 hardwood litter         154      -1.151
   4 mature chaparral        83.9     -1.040
   5 young chaparral         94.4     -0.881
   6 dormant brush           74.1     -0.906
   7 southern rough          75.1     -0.884
   10 overgrown slash        66.1     -0.810
   11 light conifer slash    63.6     -0.776
   12 medium conifer slash   69.6     -0.818
   13 heavy conifer slash    71.8     -0.835
   ========================  =======  =======

   Model 8 (closed timber litter) is not in the table. A launch cell whose
   code has no entry (Anderson 8, Scott-Burgan codes, custom fuels) and every
   burning structure use :cpp:`erf.fire.spotting.thermal_fuel_model`. When
   that is 0, the run stops at start-up, naming the input, if
   :cpp:`spotting.fuel_system = "40"`, if structure ignition is on, if
   :cpp:`fuel_model_id` has no entry, or if the fuel map has a code without
   one on a cell with fuel.
   The fits were made over 10 m winds of 5 to 30 m/s (INT-309 Tables 1-3);
   since :math:`B < 0` they grow without bound as the wind falls, so
   :math:`U_{10}` is held to that range. INT-309's example 1 (short grass,
   2000 kW/m, 5 m/s) gives :math:`E` = 30 660 kJ/m and :math:`H_z` = 30 m
   (with :math:`B = -1.465`; Table 4 prints -1.456, which gives 30.5 m), and
   example 2 (mature chaparral, 50 000 kW/m, 20 m/s) gives 186 000 kJ/m and
   75 m.

   Until November 2026 the height was :math:`12.2\, I_B^{1/3}`, cited as
   INT-309's equation 1. That form is not in the paper (its equation 1 is the
   intensity at which a plume separates from the surface,
   :math:`I = 17.5\,U^3`); it has no wind dependence and gave 154 m at
   example 1, five times the paper's height.
4. A brand is launched with probability :cpp:`erf.fire.spotting.P_base` per
   front cell per application (default 0.01).
5. The brand falls at :cpp:`erf.fire.spotting.terminal_velocity` (default
   0.5 m/s) while drifting with the wind at its height: the reference wind
   at the height it was sampled (:cpp:`erf.fire.wind_ref_ht`, or 6.096 m for
   every cell with :cpp:`use_per_fuel_wind_ht`) scaled by the neutral log profile
   :math:`\ln(z/z_0) / \ln(z_{ref}/z_0)` with :math:`z_0` =
   :cpp:`erf.fire.wind_sample_z0`, so a brand lofted to 200 m rides a wind
   about twice the 6.1 m wind (until October 2026 it drifted on the
   WAF-reduced midflame wind, about a quarter of that). The fall is integrated by
   forward Euler in sub-steps of height :math:`\Delta z = H_z / N` and time
   :math:`\Delta t = \Delta z / v_t`, with :math:`N` =
   :cpp:`erf.fire.spotting.n_traj_steps` (default 20): :math:`N` sub-steps
   bring the brand from :math:`H_z` to the ground of its source. The terrain
   surface is interpolated under the brand, so drift downslope lands farther
   away; over ground lower than the source the brand steps on, up to

   .. math::

      N_{fall} = \max\left(N, \left\lceil \frac{z_s + H_z - z_{min}}{\Delta z} \right\rceil\right)

   sub-steps, with :math:`z_s` the source's ground and :math:`z_{min}` the
   lowest ground in the domain, so every brand lands. In the sub-step where
   the brand's height above the ground goes from :math:`a_0 > 0` to
   :math:`a_1 \le 0`, it lands at the fraction :math:`a_0 / (a_0 - a_1)` of
   that sub-step (the height above the ground taken linear over it). Until
   November 2026 :math:`\Delta z` was :math:`(z_s + H_z - z_{min}) / N`, the
   fall to the lowest ground in the domain, and the brand landed at the end
   of its last sub-step: a 30 m loft from ground 500 m above that lowest
   ground (5 m/s at 6.1 m, z0 = 0.1 m) took sub-steps of 26.5 m and landed
   74 % too far (310 % at 2000 m).
   With :math:`N = 20` the remaining forward Euler error of the distance is
   about 2.5 % for that 30 m loft on flat ground and about halves each time
   :math:`N` doubles. A brand with :math:`N_{fall} > 10^6` (a tiny loft over
   high relief, or a very large :math:`N`) stops the run with a message
   naming :cpp:`n_traj_steps` and :cpp:`I_B_min`.

   A brand stays in flight up to the faces of the fire grid. In the outer
   half of an edge cell the wind, ground and roughness are those of the edge
   cell; a brand that crosses a non-periodic face is lost, in flight or in
   the step in which it meets the ground. Across a periodic face the brand
   flies on from the other side: the fields are read at its position wrapped
   into the domain, the landing distance (and the cap of the source's fuel)
   is that of the unwrapped flight, the landing cell is found from the
   wrapped position, and the spot disc is measured to the nearest image of
   the landing, so a disc next to the seam ignites cells on both sides of
   it. Until November 2026 a brand past the outermost cell centre was
   dropped half a cell before the face, every face, periodic or not, was
   treated as an exit, and a brand that crossed a face in its last step
   landed in the edge cell.

   With :cpp:`erf.fire.spotting.drift_wind = frozen` (and
   :cpp:`erf.fire.frozen_wind.file` set) the brand drifts instead on
   :math:`u, v` of the frozen 3-D GUST-CFD field at its height
   :math:`h = z_b - z_g(x, y)` above the ground under it. Each of the four
   GUST columns around the brand is sampled at its own ground :math:`+ h`,
   linearly between its cell centres, and the four are weighted bilinearly
   (or the column holding the brand is used, with
   :cpp:`erf.fire.wind_interp = nearest`), the same sampling as the
   reference wind; at GUST's lateral edges the edge column is used (no
   periodic wrap: the brand's position is wrapped as above, but the frozen
   stencil does not blend the two edge columns across a periodic seam, so
   in the outer half of an edge cell it takes the edge column's wind). Below the first cell centre the lowest centre's wind is
   used, or the log profile of roughness :cpp:`erf.fire.wind_sample_z0`
   with :cpp:`erf.fire.wind_below_first_cell = log`. The reference height
   :cpp:`erf.fire.wind_sample_ht` and the per-fuel reference heights do not
   apply to the drift. Above the highest cell centre the brand takes that
   centre's wind; such brands are counted and reported at each spotting
   step, and the run goes on. The lofting height still uses the reference
   wind brought to 10 m (step 3). With a frozen wind that reference wind is
   GUST's own wind at each fire cell, so the terrain speed-up is already in
   it; only the step from the reference height to 10 m follows the
   log profile rather than GUST's profile. In the ``FireFrozenWind`` hill
   case GUST's 10 m wind differs from the log-scaled one by -2 % to +8 %;
   with the 5 m/s floor of the fit the lofting height would change by less
   than 2 %, in 54 of the 1024 GUST columns. With a 10 m reference height
   there is no log step. With
   :cpp:`frozen` the fall stays at the terminal velocity: the
   vertical wind :math:`w` is not used. Over terrain :math:`w` can be a
   sizeable fraction of :math:`v_t` (up to 0.6 m/s in the ``FireFrozenWind``
   hill case), so the flight time on a windward slope can differ from the
   one computed here; :cpp:`frozen_uvw`, below, includes :math:`w`.

   With :cpp:`erf.fire.spotting.drift_wind = frozen_uvw` the brand also
   moves with GUST's vertical wind :math:`w`, sampled as :math:`u, v` (below
   the first cell centre it is clamped or scaled with them, so the wind
   keeps its direction there), and its height changes by
   :math:`(w - v_t)\,\Delta t` per sub-step. :math:`\Delta t` stays the step of
   the :math:`v_t` fall, so in an updraft the brand needs more sub-steps: it
   steps on until it lands, up to :math:`10 N_{fall}` (:math:`10 N` on flat
   ground), and a brand still aloft then
   is dropped, counted and reported. In uniform :math:`u = U`, :math:`w = W
   < v_t` on flat ground it lands after

   .. math::

      n = \left\lceil \frac{N v_t}{v_t - W} \right\rceil

   sub-steps, :math:`\Delta t = H / (v_t N)`. It lands inside the last of
   them, where its height above the ground, linear over the sub-step, is zero
   (as with ``frozen`` and ``profile``); in uniform :math:`u` and :math:`w`
   this is the exact :math:`D = U H / (v_t - W)`, which the unit tests
   ``FrozenDrift.*`` pin (until November 2026 it landed at the end of that
   sub-step, up to one sub-step further). With :math:`W \ge 0.9\,v_t` it
   needs more than :math:`10 N` sub-steps and is dropped.

   The share of :math:`w` in the height change is summed apart from the
   :math:`v_t` fall. When :math:`w` has moved the brand, it lands within
   :math:`10^{-3}\,v_t \Delta t` above the ground, so the round-off of the
   sub-steps cannot add one; this can land it one sub-step before the
   :math:`n` above when :math:`N v_t / (v_t - W)` lies within
   :math:`10^{-3}\,v_t/(v_t - W)` above an integer. A brand that :math:`w`
   has not lifted overall lands after :math:`N_{fall}` sub-steps at the
   latest, as with ``frozen``. With :math:`w = 0` it therefore lands every brand where
   ``frozen`` does, bit for bit, in double and single precision.

   With the default :cpp:`erf.fire.wind_below_first_cell = clamp`, :math:`w`
   below the first cell centre keeps that centre's value down to the
   ground. On a steep windward slope where that value is close to or above
   :math:`v_t` a brand low over the ground can hover there and be dropped
   instead of landing; with :cpp:`log`, :math:`w` falls to zero at the
   ground with :math:`u` and :math:`v`.

   In a wind :math:`u = a\,h` on flat ground, with every sub-step between
   the first and the highest cell centres, the scheme lands a brand lofted
   to :math:`H` at

   .. math::

      D = \frac{a H^2}{2 v_t}\,\frac{N + 1}{N},

   :math:`v_t` the terminal velocity and :math:`N` the number of sub-steps;
   the unit tests ``FrozenDrift.*`` pin this value. The exact distance is
   :math:`a H^2 / (2 v_t)`: each sub-step takes the wind at its start
   height, so the scheme overestimates by :math:`1/N` (5 % at the default
   :math:`N = 20`). The same holds for the profile drift. The sub-step does
   not depend on the ground elsewhere in the domain (until November 2026 a
   brand on raised ground took fewer, larger sub-steps, sized for a fall to
   the lowest ground, and the error grew to 9 % at :math:`N = 20` in the
   ``FireFrozenWind`` hill case).
6. The landing distance is capped at the maximum spotting distance of the
   source cell's fuel (its code on a fuel map, the domain fuel otherwise; a
   burning structure's footprint cell takes the map's code under it, the
   tables' 500 m for a non-burnable code) from Scott (2006) and Albini, for
   the fuel system selected by
   :cpp:`erf.fire.spotting.fuel_system`, ``"13"`` (Anderson) or ``"40"``
   (Scott and Burgan).
7. The brand ignites with probability :cpp:`erf.fire.spotting.P_catch`
   (default 1) and only if the landing cell still holds at least
   :cpp:`erf.fire.spotting.reentry_fuel_thresh` (default 0.05) of its initial
   fuel load, so brands landing in the burned area or on consumed fuel are
   discarded.
8. A disc of radius :cpp:`erf.fire.spotting.spot_radius` (default 10 m) is
   stamped negative into :math:`\phi` at the landing point, in the
   convention of the propagation path: normalised by the radius on the
   FARSITE path, as a signed distance in metres on the level-set path. The
   landing removes no fuel: the spot burns its cells' load like any other
   burning cell (until October 2026 it stripped 95 % of the load from every
   cell of the disc, losing its heat). A brand that lands on a non-burnable
   cell is dropped and, with the exposure diagnostics on, counted on that
   cell.

.. note::

   Until the exposure work the spotting step clamped :math:`\phi` to
   :math:`[-1, 1]` on every call regardless of the path. On the level-set
   path, where :math:`\phi` is a distance in metres, that squashed the
   distance function every spotting step and changed the spread even when
   no brand landed. Level-set runs with spotting on therefore differ from
   earlier builds; FARSITE runs do not.

:cpp:`erf.fire.spotting.random_seed` fixes the draws for reproducibility;
zero seeds from the clock. The four-component diagnostic field
``fire_albini_data`` records, per cell, the lofting height at source cells,
the number of brands launched, the maximum landing distance reached from
that cell, and a flag marking cells that received a spot ignition; the
statistics CSV adds the number of spot fires per step and the largest
landing distance. The spotting sources and landings are on the fire grid
and use the fire-grid wind, so brands do not see the three-dimensional
plume; the terrain reader supplies the ground under the trajectory.

Crown fire
----------

:cpp:`erf.fire.crown.enable` adds the transition from surface fire to crown
fire and an active crown rate of spread. The canopy is described by uniform
parameters: base height :cpp:`erf.fire.crown.canopy_base_ht` (m), bulk
density :cpp:`erf.fire.crown.canopy_bulk_den` (kg/m³), depth
:cpp:`erf.fire.crown.canopy_depth` (m), foliar moisture
:cpp:`erf.fire.crown.foliar_moisture` (fraction) and heat content
:cpp:`erf.fire.crown.h_crown_BTU_lb`. Cruz's regression was fitted to North
American conifer stands and is not meant for grass or shrub fuel models.

**Initiation.** A surface fire crowns where its Byram intensity exceeds the
Van Wagner (1977) critical intensity,

.. math::

   I_{0} = \bigl(0.010\, C_{BH}\, (460 + 25.9\, m)\bigr)^{3/2} \quad [\mathrm{kW/m}],

with :math:`C_{BH}` the canopy base height [m] and :math:`m` the foliar
moisture content in percent of dry mass (:cpp:`erf.fire.crown.foliar_moisture`
times 100). At the defaults (5 m, 100 %) :math:`I_0` = 1883 kW/m. Until
September 2026 the 3/2 power was missing and the moisture entered as the excess
over a "critical" value, so the threshold was 16x lower and dry foliage could
not crown. Once a cell has crowned it stays crowned.

**Active crown rate of spread.** :cpp:`erf.fire.crown.ros_model` selects

- ``"cruz"`` (default), Cruz, Alexander and Wakimoto (2005),

  .. math::

     R_{crown} = \frac{11.02}{60}\, U_{10}^{0.90}\, C_{BD}^{0.19}\, e^{-0.17\, M_{10}} \quad [\mathrm{m/s}],

  with :math:`U_{10}` the 10 m wind in km/h (the reference wind at
  :cpp:`erf.fire.wind_ref_ht` brought to 10 m with the neutral log profile of
  :cpp:`erf.fire.wind_sample_z0`, an open wind, not the WAF-reduced midflame
  wind; or fixed by :cpp:`erf.fire.crown.wind_10m_kmh` when positive; until
  November 2026 the 6.1 m wind was used as the 10 m wind, 11 % low) and
  :math:`M_{10}` the 10-hour dead moisture in percent;
- ``"rothermel1991"``, :math:`R_{crown} = 3.34\, R_{surface}`;
- ``"van_wagner_proxy"``, :math:`R_{crown} = (3 / C_{BD})\, f(M_{fol}) / 60`.

With :cpp:`erf.fire.crown.use_passive_blend` the transition is continuous
through the crowning fraction :math:`CF = \min\bigl((I_B/I_{crit})^{2/3}, 1\bigr)`,
:math:`R = (1 - CF)\, R_{surface} + CF\, R_{crown}`; otherwise the crown rate
replaces the surface rate in crowned cells. Crowned cells are burned, so
the rate is also carried into the unburned band ahead of them, every cell
within :math:`2 + \lceil R_{max}\Delta t / \Delta x \rceil` cells taking
the largest crown rate of its burned or already reached neighbours when
that exceeds its surface rate, and the ratio to the surface rate is a
per-cell factor that the directional level-set stages multiply into the
rate they rebuild, as the acceleration factor is (until October 2026 the
crown rate reached the front only on the FARSITE path and the isotropic
level set; the directional path, the default, rebuilt the surface rate at
every stage). The crown fraction burned of
Scott and Reinhardt (2001), :math:`(R - R_{surface}) / (R_{crown} - R_{surface})`,
is written as ``fire_crown_fraction_burned``.

**Heat release.** Crowned cells consume the canopy load
:math:`C_{BD} \times` depth with the e-folding time
:math:`\mathrm{depth} / R_{crown}` (at least 1 s) and add its heat to the
sensible flux as a step mean, :math:`w_c h_c (1 - e^{-\Delta t/\tau}) /
\Delta t`, so the heat handed over is the heat of the canopy consumed; it
enters the atmosphere with the same exponential injection and e-folding
depth as the surface flux. ``fire_crown_active`` and
``fire_crown_load`` are checkpointed so crowned cells stay crowned across a
restart.

References
----------

- Albini, F. A. (1983). Potential spotting distance from wind-driven surface fires. USDA Forest Service Research Paper INT-309. https://firedoc.nist.gov/article/zXczXYQBWEcjUZEYBfF9
- Scott, J. H. (2006). An analytical model for estimating spotting distance from surface fires. In: Fuels Management: How to Measure Success, USDA Forest Service RMRS-P-41.
- Van Wagner, C. E. (1977). Conditions for the start and spread of crown fire. Canadian Journal of Forest Research, 7(1), 23-34.
- Cruz, M. G., Alexander, M. E. and Wakimoto, R. H. (2005). Development and testing of models for predicting crown fire rate of spread in conifer forest stands. Canadian Journal of Forest Research, 35(7), 1626-1639.
- Rothermel, R. C. (1991). Predicting behavior and size of crown fires in the northern Rocky Mountains. USDA Forest Service Research Paper INT-438.
- Scott, J. H. and Reinhardt, E. D. (2001). Assessing crown fire potential by linking models of surface and crown fire behavior. USDA Forest Service Research Paper RMRS-RP-29.
