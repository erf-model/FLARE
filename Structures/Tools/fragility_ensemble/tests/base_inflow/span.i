# A deck for the fragility ensemble's tests: it is rewritten, never run. One span in a wind
# that enters through x-low as the profile of inflow_profile (8 m/s along +x at 10 m).
max_step = 10
erf.fixed_dt = 0.5          # s
geometry.prob_lo     = 0.    0.    0.
geometry.prob_hi     = 1500. 1000. 300.
amr.n_cell           = 30    20    12
geometry.is_periodic = 0 0 0
xlo.type           = "Inflow"
xlo.dirichlet_file = "inflow_profile"   # z u v w theta
xlo.KE             = 2.0
xhi.type           = "Outflow"
ylo.type           = "SlipWall"
yhi.type           = "SlipWall"
zlo.type           = "SlipWall"
zhi.type           = "SlipWall"
erf.init_type           = "input_sounding"
erf.input_sounding_file = "input_sounding"
erf.terrain_file_name   = "terrain.txt"
conductors.diagnostics_dir = conductors
FILE = lines.inputs
