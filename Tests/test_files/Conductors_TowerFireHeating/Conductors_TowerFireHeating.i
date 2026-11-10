# A lattice tower heated by a fire beside it (conductors.thermal.towers = true, with the conductors'
# thermal model conductors.thermal_model = ieee738). One line of two 150 m spans hangs from a
# suspension tower at (750, 500) between two fixed ends, in still air at 300 K over flat ground. The
# tower bends as a frame ERF generates from its type (8 panels of equal-leg angles below the cross-arm,
# crossed bracing, a cross-arm truss each side). FLARE's fire engine releases a prescribed 2e5 W/m^2
# over a 20 m disc under the tower, its centre 5 m along the line and 10 m along the cross-arm from
# the tower's axis (so no two mirrored members heat alike), from t = 0; nothing else burns.
# The coupling is passive, so the air keeps the sounding's temperature, and every member's steel
# temperature follows its EN 1993-1-2 heat balance alone: each burning fire cell radiates 0.3 of its
# heat release from 2 m above the ground, and the member exchanges heat with the air by convection and
# radiation. The members near the fire pass 100 C, where the steel's stiffness starts to fall, and the
# tower's frame is rebuilt at their temperatures. check_tower_heating.py recomputes every member's
# temperature from the frame the run wrote and the fire cells (README.md).
#
# The command line moves the fire and the tower onto a refined level (CTest variants).

max_step = 150
stop_time = 300.0
erf.fixed_dt = 2.0

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

# the fire: a prescribed heat disc under the tower and nothing burning
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
erf.fire.prescribed_heat.center = 755.0 510.0
erf.fire.prescribed_heat.radius = 20.0
erf.fire.write_fire_stats_csv = false
erf.fire.plot_int = -1

# the line and its tower
erf.structures.model           = conductors
conductors.lines               = S1
conductors.tower_types         = lattice
conductors.S1.end_a            = 600. 500. 30.
conductors.S1.towers           = 750. 500. 30.
conductors.S1.end_b            = 900. 500. 30.
conductors.S1.stringing_tension = 20000
conductors.S1.diameter         = 0.0281
conductors.S1.mass_per_length  = 1.628
conductors.S1.axial_stiffness  = 3.0e7
conductors.S1.tower_type       = lattice
conductors.S1.insulator_length = 2.5
conductors.S1.insulator_mass   = 60.
conductors.S1.segments         = 10
conductors.S1.output_root      = S1
conductors.diagnostics_dir     = conductors
conductors.air_density         = 1.17
conductors.drag_on_flow        = false

conductors.lattice.base_width    = 6
conductors.lattice.top_width     = 1.5
conductors.lattice.solidity      = 0.2
conductors.lattice.arm_length    = 12
conductors.lattice.arm_depth     = 1.2
conductors.lattice.frame_panels  = 8
conductors.lattice.leg_angle     = 0.15 0.012
conductors.lattice.brace_angle   = 0.10 0.008
conductors.lattice.damping_ratio = 0.02

# the conductors' temperature (needed for the towers'), as Conductors_FireHeating
conductors.thermal_model              = ieee738
conductors.thermal.radiative_fraction = 0.3
conductors.thermal.cutoff_radius      = 300.0
conductors.thermal.source_height      = 2.0
conductors.S1.absorptivity            = 0.5
conductors.S1.emissivity              = 0.5
conductors.S1.heat_capacity           = 1310.0
conductors.S1.thermal_expansion       = 1.93e-5
conductors.S1.current                 = 600.0
conductors.S1.resistance              = 7.283e-5 8.688e-5
conductors.S1.resistance_temperatures = 25.0 75.0

# the tower's steel temperature
conductors.thermal.towers                    = true
# not the defaults (0.7, 25, 0.01), so that the check shows each key is read
conductors.thermal.tower_emissivity          = 0.6
conductors.thermal.tower_convection          = 20.0    # W/(m^2 K)
conductors.thermal.tower_reduction_tolerance = 0.02
