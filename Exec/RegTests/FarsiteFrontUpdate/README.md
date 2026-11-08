# FarsiteFrontUpdate

The two front updates of the FARSITE path, `erf.fire.farsite.front_update`,
against the Richards (1990) rates they are built on.

A 50 m radius disc of short grass (Anderson fuel model 1, 5.5 % moisture) burns
for 10 minutes in a uniform 3 m/s westerly on flat ground, with 10 m fire cells
and one-way coupling. With `erf.fire.farsite.shape = ellipse` (the default
since 2026-10) the spread shape is the Richards (1990) ellipse of Anderson's
length-to-width ratio `L/W` with Alexander's (1985) head-to-back ratio
`HB = (L/W + sqrt((L/W)^2 - 1))^2`: the head runs at `R`, the back at `R / HB`
and the flanks at the semi-minor rate `a = (R + R/HB) / (2 L/W)`, and the exact
burned region is the disc swept by that ellipse. With `shape = rectangle` the
normal speed is `R` times `a cos + b |sin|` ahead of the wind and
`b |sin| - c cos` behind, with `a = 1`, `c = 0.2` and `b = 1.2 / (2 L/W)`: the
support function of the rectangle `[-c R, a R] x [-b R, b R]`, head, back and
flanks at `a R`, `c R` and `b R`, a head-to-back ratio of 5 at every wind.
Either region stays convex.

- `inputs_front_cell`: the default. The unburned cells next to the burned
  region are the front cells; each burns when the ellipse grown from its
  burned neighbours' arrival times reaches it (a Hopf-Lax update), so the front
  moves one row at a time.
- `inputs_rectangle`: the same update with the rectangle, the shape this path
  grew until 2026-10.
- `inputs_legacy`: the update before 2026-09. Every cell at or below
  `farsite.phi_threshold` stamped a target one cell ahead, and since phi was
  rebuilt as 0 on unburned cells, the first unburned row stamped too: two rows
  per cell of travel.

```bash
[MPIRUN="mpirun -np 2"] ./run_farsitefrontupdate.sh /path/to/erf_exec
```

`check_farsitefrontupdate.py` fits the head, back and flank positions over the
plotfiles and requires, for `front_cell` and `rectangle`, the head within 3 %
of `R`, the back and flanks within 10 % of their shape's rates, and the burned
area at least 0.93 (ellipse, whose union of discrete sources scallops) or
0.95 (rectangle) of its convex hull; `legacy` is reported.

Measured on 2026-11-08 (R = 0.971 m/s, L/W = 1.91, so HB = 12.6: the
ellipse's back runs at R / HB = 0.077 m/s and its flanks at b R = 0.274 m/s).
The deck pins `erf.fire.wind_below_first_cell = clamp` (the default), since
its 80 m atmosphere cells put the first centre at 40 m and the case wants the
sounding's 3 m/s itself as the midflame wind:

| update | head [m/s] | back [m/s] | flanks [m/s] | area at 600 s [ha] | area / hull |
|---|---|---|---|---|---|
| ellipse (Richards/Alexander) | 0.971 | 0.077 | 0.274 | | 1 |
| `front_cell` (ellipse, default) | 0.975 (+0.4 %) | 0.081 (+1.9 m over 540 s) | 0.276 (+1.0 m) | 24.3 | 0.945 or more |
| rectangle (Richards coefficients) | 0.971 | 0.194 | 0.305 | | 1 |
| `rectangle` | 0.975 (+0.4 %) | 0.191 (-1.8 m) | 0.308 (+1.9 m) | 34.6 | 0.973 or more |
| `legacy` (ellipse normal speed) | 1.828 (1.88 R) | 0.162 (2.09 c R) | 0.543 (1.98 b R) | 53.7 | 0.693 at 60 s |

The 2026-10-09 numbers had L/W = 5.00, because the Anderson length-to-breadth
formula was fed the 3 m/s wind as 6.7 (mi/h) in a correlation written for m/s
(Finney 1998, eq. 13). With the wind in m/s the ellipse is 2.6 times wider for
the same head rate: the back rose from 0.010 to 0.077 m/s, the flanks from
0.098 to 0.274 m/s, and the `front_cell` area at 600 s from 11.6 to 24.3 ha.

Before 2026-10 the rectangle was the only shape (`front_cell` 0.975 / 0.191 /
0.112 m/s, 19.2 ha; `legacy` 1.828 / 0.382 / 0.354 m/s, 50.1 ha).

`Tests/Unit/Fire/ERF_GTestFarsiteSpreadAccumulation.cpp` checks the same
rates cell by cell, and that the arrival times do not depend on the box
decomposition.
