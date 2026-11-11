# A deck for the fragility ensemble's tests: it is rewritten, never run. One span in a
# crosswind v = 15 + 0.1 z m/s (16 m/s along +y at 10 m), periodic in x and y.
max_step = 10
geometry.is_periodic = 1 1 0
zlo.type = "SlipWall"
zhi.type = "SlipWall"
erf.init_type = "input_sounding"
erf.input_sounding_file = "input_sounding"
erf.structures.model = conductors
conductors.lines = S1
conductors.S1.rated_breaking_strength = 13890.
