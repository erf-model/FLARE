# A deck for the fragility ensemble's tests: it is rewritten, never run. One span in a
# prescribed wind of 15 m/s along +y with 1.5 m/s upwards.
max_step = 10
geometry.is_periodic = 0 1 0
xlo.type     = "Inflow"
xlo.velocity = 10.0 0. 0.
xhi.type     = "Outflow"
erf.init_type = "uniform"
erf.structures.model = conductors
conductors.lines = S1
conductors.S1.rated_breaking_strength = 13890.
conductors.prescribed_velocity = 0. 15. 1.5
