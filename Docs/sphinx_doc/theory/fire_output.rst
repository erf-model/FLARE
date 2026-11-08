.. role:: cpp(code)
   :language: c++

.. _sec:FireOutput:

Fire Diagnostics and Output
===========================

Flame diagnostics
-----------------

Computed every fire step on burned cells (:math:`\phi < 0`) and zero
elsewhere:

- **Byram fireline intensity** :math:`I_B = h\, w_0\, R` [kW/m], with
  :math:`h` the heat content [kJ/kg], :math:`w_0` the fuel load the front
  consumes [kg/m²] (the cell's initial load) and :math:`R` the rate of
  spread, on burning cells; it is zero once the cell's remaining load falls
  below 1 % of :math:`w_0` (burned out). It drives ember lofting, crown-fire
  initiation and the flame diagnostics below. :math:`w_0` is the whole
  load, dead and live, so on heavy fuels the intensity is above Byram's
  :math:`I_R t_r R / 60` of the flaming front alone (BehavePlus): for
  Anderson 10 the 100-h class is 0.23 lb/ft² of the 0.55 lb/ft² load, so
  the crown criterion and ember launch trigger more readily there than
  BehavePlus would. Until October 2026 it was
  :math:`h\,(w_0 - w)\, R`, zero on the cell the front had just reached and
  largest far behind it; that form is kept as the plotfile variable
  ``fire_heat_release``.
- **Byram (1959) flame length** :math:`L = 0.0775\, I_B^{0.46}` [m], the SI form of
  :math:`0.45\, I_B^{0.46}` ft with :math:`I_B` in BTU/ft/s.
- **Flame temperature**, by :cpp:`erf.fire.flame_temp_method`:
  ``"byram_radiant"`` (default) :math:`T = T_a + 800\,(I_B/1000)^{0.25}`;
  ``"mcalpine_heat"``, the adiabatic estimate of McAlpine and Xanthopoulos
  (1989) from the heat content and moisture, clipped to
  :math:`T_a + 100 < T < 1500` K. :math:`T_a` is
  :cpp:`erf.fire.flame_temp_T_amb` [K]. A third method,
  ``"nelson_emissivity"``, was removed in November 2026: it put the fireline
  intensity [W/m] into the Stefan-Boltzmann law as a flux [W/m²].
- **Flame tilt** from vertical, :math:`\arctan(U / v_b)` in degrees, from
  the effective wind speed and the buoyancy velocity
  :math:`v_b = (g I_B / (\rho c_p T_a))^{1/3}` with
  :cpp:`erf.fire.flame_tilt_rho_air` and :cpp:`erf.fire.flame_tilt_T_amb` [K]. Off unless
  :cpp:`erf.fire.compute_flame_tilt` is set; a diagnostic only.

Fire plotfiles
--------------

The fire grid has its own plotfile stream, separate from the atmospheric
plotfiles, written to ``erf.fire.plot_file`` (default ``plt_fire_``) every
``erf.fire.plot_int`` steps or every ``erf.fire.plot_per`` seconds. Each is
a single-level AMReX plotfile on the fire grid, whose domain is the region
the fire grid covers (the refined region when the fire runs on a finer level),
with a small JSON sidecar giving the grid ratio, the atmospheric level of the
fire grid and the variable count. Variables appear in this fixed
order; the optional blocks are present only when their feature is on:

.. list-table::
   :widths: 30 12 58
   :header-rows: 1

   * - Variable
     - Unit
     - Present when
   * - ``fire_phi``
     - -
     - always
   * - ``fire_ros``
     - m/s
     - always
   * - ``fire_wind_eff_u``, ``fire_wind_eff_v``
     - m/s
     - always
   * - ``fire_wind_ref_u``, ``fire_wind_ref_v``
     - m/s
     - always
   * - ``fire_extract_z``
     - m
     - always
   * - ``fire_slope_x``, ``fire_slope_y``
     - -
     - always
   * - ``fire_fuel_mc_1hr``, ``fire_fuel_mc_10hr``, ``fire_fuel_mc_100hr``
     - fraction
     - always
   * - ``fire_surface_temp_K``, ``fire_surface_rh``
     - K, 0-1
     - always
   * - ``fire_heat_flux``
     - W/m²
     - always
   * - ``fire_fuel_load``
     - kg/m²
     - always
   * - ``fire_fireline_intensity``, ``fire_flame_length``
     - kW/m, m
     - always
   * - ``fire_arrival_time``
     - s
     - always
   * - ``fire_heat_release``
     - kW/m
     - always (:math:`h\,(w_0 - w)\, R`, the intensity's form before 2026-10)
   * - ``fire_fuel_mc_lh``, ``fire_fuel_mc_lw``
     - fraction
     - live moisture components present (BEHAVE path)
   * - ``fire_lofting_height``, ``fire_spot_brand_count``, ``fire_spot_max_dist``, ``fire_spot_active``
     - m, -, m, -
     - ``erf.fire.spotting.enable``
   * - ``fire_flame_tilt``
     - deg
     - ``erf.fire.compute_flame_tilt`` (crown fire on or off)
   * - ``fire_crown_active``, ``fire_crown_load``, ``fire_crown_fraction_burned``
     - -, kg/m², -
     - ``erf.fire.crown.enable``
   * - ``fire_flame_temp``
     - K
     - ``erf.fire.crown.enable``
   * - ``fire_fuel_model_code``
     - -
     - a spatial fuel map is loaded
   * - ``fire_ros_model_weight``
     - 0-1
     - ``erf.fire.ros_model = hybrid``
   * - ``fire_structure_height``
     - m
     - ``erf.fire.structures.enable`` or the hybrid ``structure`` selector
   * - ``fire_nonburnable``
     - 0/1
     - structures, ``fuel_map.nonburnable_codes``, ``firebreak.use_mask`` or ``suppression.enable`` configured
   * - ``fire_structure_id``, ``fire_heat_load``, ``fire_peak_intensity``, ``fire_ember_landings``
     - -, J/m², kW/m, -
     - ``erf.fire.exposure.enable`` (see below)
   * - ``fire_structure_state``, ``fire_structure_ignition_time``, ``fire_structure_rad_flux``
     - -, s, W/m²
     - ``erf.fire.structures.ignition.enable`` (:doc:`wui_structure_ignition`)
   * - ``fire_precip_mm_hr``
     - mm/hr
     - dynamic dead-fuel moisture with ``erf.fire.precip_source = atmosphere`` or a positive ``erf.fire.precip_rate_mm_hr``
   * - ``fire_suppression_mask``
     - -1/0/1
     - ``erf.fire.suppression.enable``; 1 suppressed this step, -1 a line cell that failed its hold test
   * - ``fire_ros_factor``
     - 0-1
     - ``erf.fire.suppression.enable``; the drop multiplier on the rate of spread, 1 outside drops
   * - ``fire_line_progress``
     - -
     - ``erf.fire.suppression.enable``; ordinal of the line a built cell belongs to, 0 elsewhere

Fire statistics CSV
-------------------

When :cpp:`erf.fire.write_fire_stats_csv` is true (default) one line per
fire step is appended to :cpp:`erf.fire.fire_stats_csv_file` (default
``fire_stats.csv``) with the columns

``step, time_s, burned_area_ha, perimeter_km, active_front_cells, head_ros_ms,
major_axis_m, minor_axis_m, heat_flux_max_Wm2, spot_fires_this_step,
max_spot_dist_m, edge_band_cells, edge_contact_time_s, precip_max_mm_hr``.

Burned area and perimeter come from the arrival-time field, the head rate is
the maximum rate of spread over burning cells, and the axes are those of the
burned region's bounding ellipse. The spotting columns are zero unless
spotting is on. The last two are the boundary guard of
:ref:`sec:FireGridEdge`: the number of burning cells within
:cpp:`erf.fire.boundary_guard_cells` of a non-periodic edge of the fire grid
on this step, and the time the fire first entered that band, :math:`-1`
until it does. Both are zero and :math:`-1` with
:cpp:`erf.fire.boundary_guard_action = "none"`.

Arrival-time probes
-------------------

:cpp:`erf.fire.probes = x1 y1 x2 y2 ...` [m] lists points at which the
arrival time is reported. Each probe prints one line,
``[FIRE PROBE] n x= y= cell=(i,j) arrival_time_s=``, on the fire step its
cell first burns, and nothing afterwards; a probe outside the domain is
reported once and ignored. The obstacle regression suite uses these to
tabulate arrival at building faces.

Structure exposure
------------------

With :cpp:`erf.fire.exposure.enable` (which needs
:cpp:`erf.fire.structures.enable`) the fire layer reports what each building
experienced. At initialisation the points of the structure heightmap above
:cpp:`erf.fire.structures.min_height` are grouped into connected footprints
(4-connected on the file's grid, numbered in scan order) and each fire cell
takes the id of its nearest point, the same sampling as the mask, so
``fire_structure_id`` is nonzero exactly where the mask holds a structure.
Three accumulators run on the fire grid: the heat load
:math:`\int Q\,dt` [J/m²] of every cell, the running maximum of the Byram
fireline intensity [kW/m], and the number of embers that landed on each
masked cell (brands the spotting model drops there instead of igniting).

Every :cpp:`erf.fire.exposure.interval` fire steps one line per structure is
appended to :cpp:`erf.fire.exposure.file` (default ``fire_exposure.csv``),
evaluated over the *wall band*: the burnable cells within
:cpp:`erf.fire.exposure.ring` fire cells of the footprint (a cell between
two footprints counts for both). The columns are

``time_s, structure_id, x_m, y_m, height_m, footprint_cells, wall_cells,
wall_burned_frac, t_first_s, t_last_s, residence_s, peak_intensity_kWm,
heat_load_mean_MJm2, heat_load_max_MJm2, embers``:

the footprint centroid and sampled height, the size of the footprint and
of its band, the fraction of the band burned, the earliest and latest
arrival of the front in the band and their difference (how long the front
spent passing the structure; both are -1 until it arrives), the largest
peak intensity in the band, the mean and largest heat load there in MJ/m²,
and the embers on the footprint. The first time a structure's band burns a
line ``[FIRE EXPOSURE] structure n x= y= reached_at_s=`` is printed, and
every row set ends with a summary of how many structures have been reached
and the largest peak intensity and heat load among them. A restarted run
appends to an existing file and writes the header only into a new one.

The four fields are added to the fire plotfile as ``fire_structure_id``,
``fire_heat_load``, ``fire_peak_intensity`` and ``fire_ember_landings``, and
the three accumulators are checkpointed as ``FireHeatLoad``,
``FirePeakIntensity`` and ``FireEmberLandings``. The heat load is the
fire's own release at the ground next to the wall, not a radiant flux onto
the wall; a view-factor model of the wall flux from the flame length and
tilt of the band cells is the natural next step. ``Exec/RegTests/FireExposure``
runs the obstacle scenario with and without immersed-forcing buildings and
with spotting, and ``FireRestart`` checks that the CSV resumes exactly.

With :cpp:`erf.fire.structures.ignition.enable` the structures ignite from
these accumulators, burn and load their neighbours (:ref:`sec:WUIStructureIgnition`).
The CSV then gains ``state, t_ignition_s, cause, structure_flux_Wm2,
incident_flux_max_Wm2`` (the columns are absent when the option is off), the
plotfile ``fire_structure_state``, ``fire_structure_ignition_time`` and
``fire_structure_rad_flux``, and the checkpoint ``FireStructureState``.

Models in use
-------------

Every run prints, once at start-up and before the first step, the model or
option it uses for each piece of fire physics, whatever
:cpp:`erf.fire.fire_debug` says: the rate-of-spread model and its options,
the front propagation and its directional options, the wind (interpolation,
height, sampling, wind adjustment factor, terrain factors), the fuel, burnout
and moisture models, the coupling to the atmosphere, the flame temperature,
and whether acceleration, spotting, crown fire, threshold ignition, smoke,
structures and suppression are on, with the ignition sources. Each line gives
the ``erf.fire`` key, the value the run uses, and whether the inputs set it or
it is the default; a value the code forced says so. An excerpt::

  [FIRE] Models (erf.fire key = value; set in the inputs, or the default):
  [FIRE]   Rate of spread                      ros_model = rothermel  (set)
  [FIRE]     Rothermel reaction velocity       reaction_velocity_formula = albini  (default)
  [FIRE]     midflame wind limit               use_wind_limit = true  (default)
  [FIRE]     wind limit form                   wind_limit = rothermel  (default)
  [FIRE]   Front propagation                   propagation_method = levelset  (set)

An option is listed only when it acts in the run: a model's own options only
when that model is in use (the Cheney-Gould pasture with ``ros_model =
cheney_gould``, the level-set gradient with ``propagation_method = levelset``),
the heat-injection options only when heat reaches the atmosphere, the wind
interpolation only when the wind is not prescribed.

.. _sec:FireFarsiteWinds:

FARSITE and FlamMap gridded winds
---------------------------------

With :cpp:`erf.fire.farsite_winds.interval` > 0 the fire layer writes its
reference wind as gridded winds for FARSITE and FlamMap. The layout is the
one WindNinja writes (``farsiteAtm.cpp`` and ``ascii_grid.cpp`` in
`firelab/windninja <https://github.com/firelab/windninja>`_), so the files
load in the same way as WindNinja output. In :cpp:`erf.fire.farsite_winds.dir`:

* ``<prefix>_<step>_vel.asc``: the wind speed, 2 decimals.
* ``<prefix>_<step>_ang.asc``: the direction the wind blows **from**, in
  degrees clockwise from north, 0 decimals. A wind toward the north is 180.
* ``<prefix>_<step>_vel.prj`` and ``_ang.prj``: copies of
  :cpp:`erf.fire.farsite_winds.prj_file`, when it is given.
* ``<prefix>.atm``: the atmosphere file that lists the grids by date and time:

.. code-block:: text

   WINDS
   ENGLISH
   12 30 1100 flare_000001_vel.asc flare_000001_ang.asc
   12 30 1101 flare_000061_vel.asc flare_000061_ang.asc

The grids are ESRI ASCII with a tab-separated header (``ncols``, ``nrows``,
``xllcorner``, ``yllcorner``, ``cellsize``, ``NODATA_value`` = -9999) and the
northern row first. Cells that no fire box covers hold ``NODATA_value``.

**Which wind.** The grids hold ``fire_wind_ref``: :math:`u` and :math:`v`
at the reference height above ground, before the wind adjustment factor and
FLARE's parametric terrain correction (:cpp:`erf.fire.use_terrain_wind`). The
ERF wind itself already includes the resolved terrain, as a WindNinja field
does. This is the open wind that FARSITE expects as input. With lagged or
passive coupling the wind is the one at the step's start time; with
synchronous coupling it is extracted after the dynamics, and is stamped at
the step's end.
The speed and direction are

.. math::

   S = c\,\sqrt{u^2 + v^2}, \qquad
   \theta_\mathrm{from} = \operatorname{atan2}(-u, -v)\ \text{in degrees, wrapped to } [0, 360),

where :math:`u` is the east component [m/s], :math:`v` the north component
[m/s], and :math:`c` = 3600/1609.344 (mph) or 3.6 (km/h). The fire grid's x
is taken as east and y as north. With :cpp:`erf.fire.farsite_winds.coarsen`
= :math:`n`, each output cell is the mean of :math:`u` and :math:`v` over
:math:`n \times n` fire cells, and the speed and direction come from that
mean vector. The speed of a mean vector is at most the mean speed. For a block
whose directions spread by 30 deg the written speed is about 2 % lower,
and for 90 deg about 16 % lower.

**Units fix the height.** In a ``.atm`` file, ``ENGLISH`` means mph at 20 ft
above the vegetation, and ``METRIC`` means km/h at 10 m. The run therefore
aborts unless the reference height matches the units: 6.096 m for
``english`` (:cpp:`erf.fire.wind_ref_ht` = 6.096, or
:cpp:`erf.fire.use_per_fuel_wind_ht`), or 10 m for ``metric``. The default
:cpp:`wind_ref_ht` of 6.1 m is within the 1 cm tolerance of 20 ft.

**Dates.** The ``.atm`` lines carry month, day and HHMM in UTC. They are
taken from the top-level :cpp:`start_datetime` plus the simulated time,
rounded to the nearest minute. The run aborts without :cpp:`start_datetime`,
also with ``init_type`` WRFInput or Metgrid, which otherwise take the date
from the file,
and an interval below 60 s also aborts.

**When.** A step whose start time is at or past the next multiple of the
interval writes the wind it extracted at that start time. A run that starts
at :math:`t` = 0 therefore writes at 0, :math:`\Delta`, :math:`2\Delta`, and
so on, each stamped with its step's start time, which lags the multiple by
less than one time step. A time step longer than the interval skips grids
(a warning is printed once). Two lines in the same minute would break
FARSITE's time order, so the run aborts when that happens: keep the time step
below the interval minus 30 s. The write times follow from the simulated time alone, so a restart
needs nothing from the checkpoint. A restart appends to the ``.atm`` file it
finds, after dropping the lines of steps at or after its own (a restart from
an older checkpoint); a file with the other unit line aborts the run. A run
that starts at :math:`t` = 0 begins a new file.
Every write gathers the output grid onto one rank (with ``coarsen`` > 1 the
block means are taken first, on the ranks that own the fire grid). For a
large fire grid, use ``coarsen`` to keep that grid small. A non-finite wind
is written as ``NODATA_value``.

**Map position.** :cpp:`erf.fire.farsite_winds.xllcorner` and
:cpp:`yllcorner` place the grids on the FARSITE landscape (the landscape's
lower-left corner, in its projected coordinates). Without them, the fire
grid's own ``prob_lo`` is written. The fire cells must be square.

**Other formats.** ``Exec/Tools/farsite_winds_convert.py`` reads the
``.atm`` file (or one ``_vel.asc``) and writes VTK (standard library only),
shapefile wind arrows (pyshp), Google Earth KMZ (pyproj and a projection),
a PDF map (matplotlib), or a geospatial PDF of the speed (GDAL):

.. code-block:: bash

   python3 Exec/Tools/farsite_winds_convert.py farsite_winds/flare.atm --to vtk shp kmz --stride 4

The regression case ``Exec/RegTests/FireFarsiteWinds`` checks the files
against a prescribed wind in closed form and against the fire plotfile.

Debug output
------------

:cpp:`erf.fire.fire_debug` prints, every step, the wind extraction range,
the rate-of-spread maximum and mean over burning cells, the number of
subcycles, the count of active fire cells, the flux and intensity maxima,
and, where relevant, the hybrid weight summary and the coupling tendency
check. The regression scripts parse these lines, so their formats are
stable.

Checkpoint and restart
----------------------

The atmospheric checkpoint carries the fire state: fields under the
``Level_<lev>`` directory of the level the fire grid refines, and the text
files ``FireState`` and ``FireSuppression`` at the top of the checkpoint.
The fire engine writes its own fields and the
coupling layer (:ref:`sec:FireCoupling`) writes the lagged flux buffers
``FireQAtmPrev`` and ``FireQLatAtmPrev`` that the next step injects. The
native engine writes:

* always: ``FirePhi``, ``FireArrivalTime``, ``FireHeatFlux``,
  ``FireFlameLength``, ``FireFuelLoad``, ``FireFuelMC``, ``FireDispAccum``,
  ``FireCrownFractionBurned`` and the text file ``FireState`` (the step and
  level-set subcycle counters, the heat-flux partition factor, the boundary
  guard's contact time and reach flag, the start of the next ignition window,
  the anchor level and the corner of the fire region, and the already
  reported probes and structures);
* with crown fire, ``FireCrownActive``, ``FireCrownLoad`` and
  ``FireCrownRosActive``; with the stick moisture model, ``FireStickMC``;
  with :cpp:`erf.fire.precip_source = atmosphere`, ``FirePrecipAccumPrev``;
  with exposure, ``FireHeatLoad``, ``FirePeakIntensity`` and
  ``FireEmberLandings``; with structure ignition, ``FireStructureState``;
  with spotting, ``FireAlbiniData``; with the temporal acceleration model,
  ``FireAccelState``;
* with suppression, the four fields ``FireSuppressionMask``,
  ``FireSuppressionFactor``, ``FireSuppressionProgress`` and
  ``FireSuppressionFailed`` with the action state in ``FireSuppression``
  (:ref:`sec:FireSuppression`).

On restart the fire engine is initialised from the inputs first, so a
spatial fuel map, firebreaks, a hybrid weight, a structure mask or a
suppression action file must still be available, and the checkpointed fields
are then read over the initial ones: the engine's first, then the lagged
fluxes. A field missing from an older checkpoint keeps its initial value, and
a checkpoint that holds the fire on another level or region than the run's
aborts. Diagnostics not in these lists are recomputed on the first step. A
fire run restarts on the level-0 boxes it was written with: a restart that
re-makes them (:cpp:`erf.regrid_level_0_on_restart`, or more ranks than the
checkpoint has level-0 boxes) stops with a message saying so.

References
----------

- Byram, G. M. (1959). Combustion of forest fuels. In: Forest Fire: Control and Use, McGraw-Hill, 61-89.
- McAlpine, R. S. and Xanthopoulos, G. (1989). Predicted vs. observed fire spread rates in Ponderosa pine fuel beds. Proceedings of the 10th Conference on Fire and Forest Meteorology.
