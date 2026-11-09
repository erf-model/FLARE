# A conductor span heated by a fire below it (conductors.thermal_model = ieee738). A 300 m span of a
# Drake-like ACSR conductor (28.1 mm, 1.628 kg/m, 1.5 m slack, 600 A) hangs between two points 30 m
# above flat ground in still air at 300 K. FLARE's fire engine releases a prescribed 2e5 W/m^2 over a
# 40 m disc under the middle of the span from t = 0 (erf.fire.prescribed_heat.*; nothing burns
# otherwise). The coupling is passive: the fire puts no heat into the air, so the air stays at rest
# and at the sounding's temperature, and the conductor's temperature follows its heat balance alone.
# Each burning fire cell radiates 0.3 of its heat release from 2 m above the ground. The conductor
# starts at its steady temperature without the fire (Joule heating against convection and radiation),
# heats under the fire, lengthens and sags. check_fire_heating.py recomputes the hottest node's
# irradiance from the fire cells and the node positions, steps its heat balance with its own IEEE 738
# terms, and checks the thermal strain and the span's sag (README.md).
#
# The command line moves the fire and the lines onto a refined level (CTest variants).

max_step = 200
stop_time = 100.0
erf.fixed_dt = 0.5

geometry.prob_lo     = 0.    0.    0.
geometry.prob_hi     = 1500. 1000. 300.
amr.n_cell           = 30    20    12
amr.max_grid_size    = 10       # box lengths divide by the fire grid ratio
amr.max_grid_size_z  = 64       # the fire module forbids a decomposition in z
amr.blocking_factor  = 2
geometry.is_periodic = 1 1 0
zlo.type = "surface_layer"      # the fire module needs the surface layer at the ground
zhi.type = "SlipWall"
erf.most.z0 = 0.1
erf.most.surf_temp = 300.0

amr.max_level = 0

erf.molec_diff_type = "ConstantAlpha"   # the surface layer needs a diffusive closure; all zero
erf.dynamic_viscosity = 0.0
erf.alpha_T = 0.0
erf.alpha_C = 0.0
erf.les_type = "None"

erf.init_type = "input_sounding"
erf.input_sounding_file = "input_sounding"   # still air, theta 300 K; surface pressure in hPa
erf.prob_name = "ABL"
erf.use_gravity = true

erf.sum_interval = -1
erf.check_int    = -1
erf.plot_int_1   = -1
erf.v = 0
amr.v = 0

# the fire: a prescribed heat disc and nothing burning
erf.fire.enable = true
erf.fire.grid_ratio = 5                 # 10 m fire cells on level 0
erf.fire.anchor_level = 0
erf.fire.propagation_method = "levelset"
erf.fire.ros_model = "rothermel"
erf.fire.fuel_model_id = 1
erf.fire.moisture_1hr = 0.08
erf.fire.moisture_10hr = 0.08
erf.fire.moisture_100hr = 0.10
erf.fire.moisture_dynamic = false
erf.fire.use_waf = false
erf.fire.use_terrain_wind = false
erf.fire.ignition_r = 0.0               # nothing burns: all the heat is the prescribed disc
erf.fire.coupling_type = "passive"      # the air keeps its temperature
erf.fire.prescribed_heat.flux = 2.0e5
erf.fire.prescribed_heat.center = 750.0 500.0
erf.fire.prescribed_heat.radius = 40.0
erf.fire.write_fire_stats_csv = false
erf.fire.plot_int = -1

# the span
erf.structures.model           = conductors
conductors.lines               = S1
conductors.S1.end_a            = 600. 500. 30.
conductors.S1.end_b            = 900. 500. 30.
conductors.S1.length           = 301.5
conductors.S1.diameter         = 0.0281
conductors.S1.mass_per_length  = 1.628
conductors.S1.axial_stiffness  = 3.0e7
conductors.S1.output_root      = S1
conductors.diagnostics_dir     = conductors
conductors.node_output_int     = 1
conductors.air_density         = 1.17
conductors.drag_on_flow        = false

# its temperature
conductors.thermal_model              = ieee738
conductors.thermal.radiative_fraction = 0.3
conductors.thermal.cutoff_radius      = 300.0       # m: every burning cell of the disc is well within it
conductors.thermal.source_height      = 2.0
conductors.S1.absorptivity            = 0.5
conductors.S1.emissivity              = 0.5
conductors.S1.heat_capacity           = 1310.0      # J/(m K): aluminium 1.116 kg/m x 955 + steel 0.512 kg/m x 476
conductors.S1.thermal_expansion       = 1.93e-5     # 1/K
conductors.S1.current                 = 600.0       # A
conductors.S1.resistance              = 7.283e-5 8.688e-5   # ohm/m at 25 and 75 C
conductors.S1.resistance_temperatures = 25.0 75.0
