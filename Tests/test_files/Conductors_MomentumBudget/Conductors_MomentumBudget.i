# The y-momentum budget of a conductor line in a crosswind (multi-level tests of the conductor model
# on ERF's structure coupling). A periodic box of uniform air (density 1, 300 K, no gravity, no
# Coriolis force, no diffusion, slip walls) blows across one line at 10 m/s. Nothing but the
# line's drag changes the air's y-momentum, so its change over the run equals the impulse the
# model reports putting into the air: the sum over the steps of the lines' level of the integrated
# source in conductors/total_load.dat (source_y, N) times the step. With a finer level the force
# must survive the averaging down: on the finest level (the default conductors.anchor_level = -1)
# or on level 0 with the finer level away from the drag. The line is thick (6 m) so that its
# impulse is about 1e-3 of the air's momentum, which fvolumesum's 6 digits resolve.
max_step  = 10
stop_time = 1.0e9
erf.fixed_dt = 0.5          # level 0: 10 m/s on 25 m cells, Courant 0.2

geometry.prob_lo     = 0.    0.    0.
geometry.prob_hi     = 600.  400.  150.
geometry.is_periodic = 1 1 0
amr.n_cell           = 24    16    6
amr.max_grid_size    = 64
amr.blocking_factor  = 2
amr.max_level        = 0
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.init_type = "Uniform"   # constant density and temperature
erf.prob_name = "ABL"
prob.rho_0 = 1.0
prob.T_0   = 300.0
prob.V_0   = 10.0           # the crosswind, along y
erf.molec_diff_type = "None"
erf.les_type        = "None"
erf.use_gravity     = false
erf.use_coriolis    = false

erf.v            = 0
amr.v            = 0
erf.sum_interval = -1
erf.check_int    = -1
erf.plot_file_1  = plt
erf.plot_int_1   = 10
erf.plot_vars_1  = density x_velocity y_velocity z_velocity conductor_fx conductor_fy conductor_fz

erf.structures.model    = conductors   # FLARE's conductor model (Structures/Conductors)
conductors.lines               = S1
conductors.S1.end_a            = 150. 200. 75.
conductors.S1.end_b            = 450. 200. 75.
conductors.S1.length           = 301.5
conductors.S1.diameter         = 6.0     # a test value: 300 N/m at 10 m/s
conductors.S1.mass_per_length  = 60.0    # heavier than the 28 kg/m of air it displaces
conductors.S1.axial_stiffness  = 3.0e9
conductors.S1.output_root      = S1
conductors.diagnostics_dir     = conductors
conductors.air_density         = 1.0
conductors.drag_on_flow        = true
conductors.epsilon             = 1.5     # the drag reaches 3 epsilon = 4.5 cells of the lines' level
