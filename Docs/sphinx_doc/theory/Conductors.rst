.. _sec:Conductors:

Conductor lines in the wind
===========================

A conductor line (``Structures/Conductors``, inputs ``conductors.*``)
is an overhead power-line conductor or shield wire hanging between two
dead-end attachments, either as a single span or as a section of spans over
suspension towers, where it hangs from insulator strings or is clamped. Its
motion in the wind is computed by MoorDyn-C (:doc:`../CouplingToMoorDyn`)
while ERF supplies the wind. The questions it answers are those of a line
exposed to a fire wind: how far the spans blow out, how their clearance to
the ground and to the neighbouring lines changes, and what tension the
attachments carry, down to the transformers the lines end on.

The lines are a structure model on ERF's structure coupling
(:ref:`sec:StructureCoupling`), built into FLARE with
``ERF_ENABLE_STRUCTURES`` (:doc:`../CouplingToMoorDyn`) and selected with
``erf.structures.model = conductors``; their inputs are read under
``conductors.*``. An embedded-boundary terrain (``erf.terrain_type = EB``)
stops the run: ERF hands a structure model no mesh heights or cell volumes on
a mesh of constant :math:`\Delta z`.

Model
-----

MoorDyn integrates the lumped-mass dynamics of each span, one MoorDyn line
per span: a span is a chain of ``segments`` elastic segments of axial
stiffness ``EA`` (N) with internal damping, each node carrying its share of
the mass per unit length, the weight net of the air's buoyancy and the drag of
the air normal to the line, relative to the node (ERF writes the added-mass
coefficients Ca and CaAx and the axial drag coefficient CdAx as 0). The fluid
is ERF's air: MoorDyn's "water" density is ``conductors.air_density``, its
free surface and flat bottom are placed ``surface_offset`` above and below
ERF's ``z = 0``, beyond the domain's top and bottom, and the fluid velocity at every line node comes from ERF
through MoorDyn's external wave-kinematics interface instead of a wave model.
The dead ends are fixed points. The tower attachments are fixed points too,
unless the towers bend (`Moving towers`_), when they are MoorDyn coupled
points that ERF moves; a line on rigid towers has no coupled degree of freedom
and needs nothing back from ERF but the wind.

For a span of chord :math:`c` (m) and unstretched length :math:`L > c` (m)
the still-air shape is close to a parabola of mid-span sag
:math:`s = \sqrt{3 c (L - c) / 8}` (m), with the horizontal tension
:math:`H = w c^2 / (8 s)` (N) for a weight :math:`w` (N/m, net of the air's
buoyancy) per unstretched metre, as MoorDyn takes the mass. MoorDyn's stationary solver finds the elastic
catenary itself, and for every level span between fixed points the start-up
log prints MoorDyn's sag and end tension against those of the elastic
catenary (see `Verification`_). A steady crosswind of speed :math:`U` (m/s)
normal to the span loads it with :math:`q = \tfrac12 \rho C_d D U^2` (N/m) per
unit length, :math:`\rho` the ``air_density`` (kg/m^3), :math:`C_d` the line's
``drag_coefficient`` and :math:`D` its ``diameter`` (m), so the span swings
out of its vertical plane towards the quasi-static blowout angle
:math:`\phi_b = \arctan(q / w)` and its tension rises with the effective
weight :math:`\sqrt{w^2 + q^2}`. The diagnostics report the swing angle of
the middle node about the chord against exactly this reference.

Sections and insulator strings
------------------------------

A line with ``towers`` is a section: several spans hanging at suspension
towers from insulator strings, or clamped to them. On a real line the
conductor runs over many towers. At a suspension tower it
hangs from an insulator string, a chain of porcelain or glass discs a few
metres long that swings freely, and the spans on either side pull along the
line with nearly the same force, so the string only swings across the line
with the wind. At a dead-end (strain) tower the conductor is anchored. A
single span between two fixed points cannot hang from suspension strings:
nothing would balance its horizontal tension, some 13 kN for the test
conductor, and the strings would be dragged along the line until nearly
horizontal. The ``towers`` of a line are therefore the suspension points
between its two dead-ends: the line becomes a section of one span more than
the towers, each span with its own unstretched ``length``, written as one
MoorDyn system. With ``insulator_length`` the conductor hangs at each tower
from a string, a short MoorDyn line from a fixed point at the tower down to a
free point where the spans either side meet, with the string's mass, the
wind on its disc diameter, an axial stiffness of 1e7 N (a 2.5 m string gives
4 mm under a 15 kN conductor) and two segments; without it the conductor is
clamped to the towers, as a shield wire is at the tower peaks.

The span quantities (sag, offset, swing, tensions, clearance, drag) are
reported per span, measured from the chord between the span's attachment
points at the towers, so the sag and the blowout of a span on strings
include the drop and the swing of its strings: the blowout is the conductor's
displacement from the vertical plane of the tower attachments, the number
that sets its clearance. Each string's angle from the vertical, its angle
across the line (positive to the left of the direction from ``end_a`` to
``end_b``, the side a positive offset is on) and the tension at its top are
reported too.

In a steady crosswind each string of a long section swings across the line
to the angle :math:`\theta` from the vertical whose tangent is the transverse
load over the vertical load it carries, the wind span over the weight span:
with equal spans of length :math:`L` (m) either side,

.. math::

   \tan\theta = \frac{q L + q_i L_i \cos\theta / 2}{w L + W_i / 2},

where :math:`q_i = \rho C_{d,i} D_i U^2 / 2` (N/m) is the wind load per unit
length of the string of length :math:`L_i` (m, ``insulator_length``) and disc
diameter :math:`D_i` (m, ``insulator_diameter``), :math:`C_{d,i} = 1` the drag
coefficient ERF writes for a string, and :math:`W_i` (N) the string's weight
net of buoyancy: the moments about the string's top, the swung string seeing
the normal wind :math:`U\cos\theta`, so that its own drag :math:`q_i L_i
\cos^2\theta` acts at half its length, normal to it.

A string takes only the weight span: a span strung to the horizontal tension
:math:`H` pulls down on its end by :math:`w h / 2 + H \Delta z / h` (N), with
:math:`h` (m) its horizontal length and :math:`\Delta z` (m) the height of the
end above the other, so a tower in a dip, the spans either side rising away from
it, can be pulled up. Its string then carries none of the conductor and flips
over the cross-arm, the line slackening on either side; in practice such a
point needs a strain tower or a taller one. The start-up log warns, naming the
line and the tower, about every string that carries less than a tenth of the
weight of the half spans either side in still air. Against the real MoorDyn, a
section whose dead ends stand 80 m above its two towers, strung to 20 kN (each
span rising 80 m over 300 m pulls up 5.3 kN a side against 2.4 kN of weight),
flags both strings, and a level section neither, its strings carrying one
span's weight each.

Lockstep with ERF
-----------------

Each line is advanced in step with ERF on one level, with the wind sampled
where its nodes are at the start of the step. The lines live on their anchor
level (``anchor_level``, ``amr.max_level`` by default, which must exist from
the start and for the whole run, or the run stops).

With ``drag_on_flow`` the drag goes into the flow of the anchor level only,
and reaches the coarser levels by averaging down. A finer level over the drag
would remove it there when it is averaged down, so after every step ERF's
structure coupling checks the drag's support, a box per spread point grown by
its line's or tower type's reach (three ``epsilon`` cells of the anchor level,
at least one cell), against the next finer level, and stops the run where
they overlap. A finer level elsewhere is allowed. Without ``drag_on_flow``
nothing goes into the flow, the support is empty, and the lines may sit on a
level a finer level covers; they then take the wind of their own level, which
under the finer level is that level's flow averaged down. ERF does not
synchronise momentum fluxes at coarse/fine faces: the momentum along the drag
that crosses a face normal to it, mostly the pressure term, differs between
the levels. With a finer patch whose faces across the wind were 75 to 150 m
from a line, the air's momentum missed the lines' impulse by 23 % to 111 %,
against 2e-6 with faces only along the wind (``Conductors_MomentumBudget*``),
so a finer level that holds lines with ``drag_on_flow`` should reach well
beyond them across the wind.

Every step of that
level, the wind is handed to MoorDyn at the line nodes (MoorDyn lists its
points after the line nodes, the attachments and the strings' lower ends, and
then one entry at its own origin; these get no wind), and the line is
advanced by ERF's step in ``substeps``
MoorDyn calls; MoorDyn sub-steps internally with its own time step, set by
the Courant factor ``moordyn_cfl`` and bounded further by ``moordyn_dt`` when
that is given. MoorDyn-C keeps the wind it is handed until the next call,
whatever time it is stamped with, so the lines take the wind of the step's
start held over the step, half a step behind its middle (an event gust is
evaluated at the step's middle). MoorDyn's clock starts at zero when the lines are created, at
ERF's time zero in a fresh run, and ERF's clock and MoorDyn's (plus ERF's time
at its zero) must agree at the end of every step.

The wind handed to MoorDyn is ERF's velocity at the start of the step,
sampled with the actuator core at the points where the line is at that
moment, not where it hung: a span blown out by several metres sits in a
different part of the flow, which matters on a ridge or near a plume. Each
component is interpolated from its own staggered grid, bilinearly in the
horizontal and linearly in the physical height, so a field linear in
``x``, ``y`` and the height is handed over exactly. The fluid acceleration
MoorDyn also accepts enters only its fluid-inertia (Froude-Krylov) load,
which in air is about a thousandth of the line's own inertia, and is left at
zero. Before every
sampling the points are checked: a point that has left the domain, or whose
surrounding cells are not on the anchor level's grids (a refined patch that
does not cover the whole line), stops the run, naming the line, the point
and its position. ``prescribed_velocity`` replaces the sampled wind by a uniform
one, for testing.

Every step the values crossing between ERF and MoorDyn are checked for NaN and
infinity: the wind handed to each line and to the towers, the node positions
and coupled-point forces MoorDyn returns, the lines' pull on the towers in the
coupling iteration, and the forces and positions spread into the flow. A
non-finite value stops the run, naming the line (or the towers) and the point.

The Courant factor matters more than its name suggests. For a 300 m span of
795 kcmil (thousand circular mils) ACSR (aluminium conductor, steel
reinforced; the Drake conductor) with 20 segments, MoorDyn's own default of 0.5 (an internal
step of 11 ms) and 0.25 make the line integration diverge within half a
second in a 20 m/s wind; 0.2 runs but more than doubles the still-air sag;
0.15 and 0.1 agree to the millimetre. ERF
therefore writes 0.1 by default and refuses values above 0.15. A MoorDyn
step that still diverges aborts with ``MOORDYN_NAN_ERROR`` and the advice to
reduce the internal step.

Verification
------------

Against the real MoorDyn-C 2.7.1, a 300 m span of 795 kcmil ACSR with 1.5 m
of slack and 20 segments in air of density 1.2 kg/m^3 (the unit tests
``ConductorVerification``):

* in a steady crosswind the mean swing angle over 20 to 60 s is the
  quasi-static blowout angle of a level span, :math:`\arctan(q/w)` (on a
  span inclined at :math:`\beta` the weight's part normal to the chord is
  :math:`w\cos\beta`): 6.023 against 6.029
  degrees at 10 m/s, 22.893 against 22.903 at 20 m/s and 43.519 against
  43.548 at 30 m/s;
* released after a one-second gust into still air, the span swings with the
  first out-of-plane period of a cable, :math:`T = 2c/\sqrt{H/m}` (s),
  :math:`m` the ``mass_per_length`` (kg/m) (Irvine,
  *Cable Structures*, 1981), which, unlike the in-plane modes, does not depend
  on the cable's extensibility (Irvine's :math:`\lambda^2`); with
  :math:`H = m g c^2/8s` it grows as the square root of the sag :math:`s`: 6.658 s
  measured over five periods against 6.655 s;
* a section of three such spans over two suspension towers, hanging from
  2.5 m strings of 60 kg, swings its strings in a steady 20 m/s crosswind to
  22.382 degrees on average over 30 to 60 s, against 22.391 degrees from the
  wind span over the weight span (0.04 %), and at most 0.11 degrees along
  the line, the spans either side balancing;
* in still air the span hangs in the elastic catenary (Irvine, *Cable
  Structures*, 1981), with the weight :math:`w` per unstretched metre as
  MoorDyn's mass: the horizontal tension :math:`H` solves
  :math:`c = HL/EA + (2H/w)\,\sinh^{-1}(wL/2H)`, the sag is
  :math:`(H/w)\,(\sqrt{1 + (wL/2H)^2} - 1) + wL^2/8EA` and the end tension
  :math:`\sqrt{H^2 + (wL/2)^2}`. For the test span this gives a sag of
  13.5836 m, a horizontal tension of 13 251 N and an end tension of 13 468 N;
  MoorDyn gives 13.599 m and 13 430 N with 20 segments and converges to
  13.5838 m with 160, its stretched length 301.6339 m equal to the
  catenary's (the parabola estimate, 12.99 m, and the inextensible catenary,
  13.01 m, are both too low). The
  start-up log prints MoorDyn's sag and end tension against the elastic
  catenary for every level span between fixed points.
* the dead ends of that span carry its weight and the catenary's horizontal
  tension: the pulls MoorDyn gives on the two fixed points add up to
  4812.82 N downwards against the line's weight of 4812.96 N, and each
  pulls 13 234 N along the chord against the catenary's 13 251 N (0.13 %).

The blowout and swing period agree within 0.1 %, the string swing within
0.3 % and the still-air sag within 0.11 %; the tests allow 1 % (0.5 % for the
catenary). The verification tests skip on the stub library, which has no line
dynamics and hangs an elastic parabola.

Stringing
---------

A line is either given an unstretched ``length`` per span or strung to a
``stringing_tension``, the horizontal tension :math:`H` every span carries in
still air, as a line crew strings a section. A section strung so hangs its
strings plumb, the spans either side of every tower pulling alike, whatever
their lengths; equal sag ratios instead would leave a short span far slacker
than a long one and drag the string between them along the line. Once the
attachments stand on the terrain, each span's length is set from the points
the conductor hangs from (the bottoms of the strings at the towers), a chord
:math:`c` over a horizontal distance :math:`h`: the parabola of horizontal
tension :math:`H` under the weight :math:`w` per unstretched metre is

.. math::

   l = c + \frac{w^2 h^4}{24 H^2 c}

long, the tension along the chord :math:`H c / h` stretches the line by that
over :math:`EA`, and the unstretched length is
:math:`l / (1 + H c / (EA\,h))`. A short span may then be shorter than its
chord unstretched, held taut by its stretch, which MoorDyn handles like any
other. On the terrain test case MoorDyn's still-air solve gives every span
of the two lines strung to 20 kN end tensions of 19.8 to 20.9 kN, the
tension along the line being a little above the horizontal one on the
inclined spans.

Placement on terrain
--------------------

The ``z`` of each attachment, the ends and the towers, is its height above
the terrain surface at its ``(x, y)``: the surface height is read from the
mesh with the actuator core's ``terrain_heights`` (the ``k = 0`` node plane,
bilinear between the nodes) and added once at start-up, so on a flat mesh
``z`` is the absolute height, and a tower on a hill top holds its conductor
that much higher. Each span's ``length`` must exceed the distance between
the points the conductor hangs from where they stand (the attachments, or the
bottoms of the insulator strings at the towers): on a slope a dead end 10 m above a summit
and a tower 30 m above the hillside below it can be closer, or further
apart, than their heights above the ground suggest, so the slack is checked
once they are placed. The attachments must lie inside the domain.
``<diagnostics_dir>/ground.dat`` (columns ``line point x y ground z``) lists
each attachment of each line (point ``a``, the towers ``t1``, ``t2``, ... from
``end_a``, then ``b``) with the terrain height found under it (``ground``, m)
and its resulting absolute height (``z``, m), and each transformer (its name
in the ``line`` column, point ``transformer``) with the terrain height under
its centre (``ground``) and the height of its top (``z``).

With an immersed terrain (``erf.terrain_type = ImmersedForcing``) the mesh is
flat and the hills are solid cells inside it, so the mesh's bottom says nothing
about the ground. The lines then take the terrain's height from the surface the
immersed boundary is built from (``erf.terrain_file_name`` or the problem's
own terrain, at the nodes of the anchor level), bilinear between its nodes as on
a fitted mesh's bottom: the ends, the towers, the transformers and the ground
under every node for the clearance. The wind is sampled at the nodes' absolute
heights, which on the flat mesh are the mesh's own. Placed on an immersed ramp,
a section and its towers stand exactly where they stand on the same ramp as a
fitted mesh. A precursor run on the flat mesh can therefore seed a run over
immersed hills from its checkpoint. An anelastic run takes no acoustic substeps,
so its immersed forcing must act on the slow step (``erf.immersed_forcing_substep``
false, the default for anelastic runs), and the point-implicit form
(``erf.if_implicit_drag = true``) keeps it stable at the flow's step: explicitly,
the drag on a solid cell, :math:`C_d |u| / \Delta z`, about 60 s\ :sup:`-1` for 20
m/s on 16 m cells, is far beyond a 0.3 s step.

Drag on the flow
----------------

With ``drag_on_flow`` the lines act on the air as well. In each step of the
anchor level ERF advances the lines first, from the wind at the start of the
step; the air's drag on every line node at the end of that MoorDyn step,
reversed (and on every node of the lattice towers, see `Towers`_), is then
spread onto ERF's face-centred momentum sources with the actuator core's
Gaussian of width ``epsilon`` cells, discretely normalised so that the source
integrates back to the force exactly, and held constant as a source over the
same ERF step, from :math:`t` to :math:`t + \Delta t`. Each node carries its
share of the line's drag, so the force per unit length is resolved as long
as the node spacing stays below the kernel width. The width is
``conductors.epsilon`` unless a line or a tower type sets its own:
``conductors.<line>.epsilon`` for the nodes of that line, its insulator
strings' included, and ``conductors.<type>.epsilon`` for the drag nodes of
the towers of that type (a line's towers take their type's width, not the
line's). Every width is in units of the anchor level's :math:`\Delta x`, and
all lines and towers share that level. The width is bounded below by the mesh,
not by the body: a conductor a few centimetres thick is never resolved, so its
lines keep the narrowest width the mesh allows, about one cell or more, while a
tower wider than a cell can take a wider kernel that matches its plan size.
Below one cell the kernel's shape depends on where the node sits in its cell
(its second moment varies by about 9 % at 0.75 cells and by a factor of two at
0.5 cells; the integral stays exact), and a tower's drag gathered on a few
faces can overshoot the wind next to it in one explicit step, so the run warns
for any width below one cell. The widths are not part of the checkpoint, so a
restart may change them. A conductor's drag is
small next to the momentum of the wind it crosses (4 to 5 N/m in a 17 m/s
wind for 795 kcmil ACSR), which is why the switch is off by default; it is
there for dense bundles, many spans in a small domain, and to close the
momentum budget. The plot variables ``conductor_fx``, ``conductor_fy`` and
``conductor_fz`` hold the source averaged to the cell centres (N/m^3) on
the anchor level and are available only with the switch on. The source goes
into the anchor level only, so with ``drag_on_flow`` no finer level may cover
the drag's support (`Lockstep with ERF`_; ERF's structure coupling stops the
run where one does), and the anchor level's grids must hold
every line and tower node with the kernel's reach, 3 ``epsilon`` cells of its
line's or tower type's width, around it (the run stops, naming the node,
otherwise). On a terrain-following or
stretched mesh a node's cell is not its height over :math:`\Delta z`, so there
the grids must hold the reach across in every cell the heights within reach
may lie in, bounded by the mesh's lowest and highest ground and its smallest
and largest level spacing: a conservative check, which over steep terrain asks
for grids over more of the height than the kernel reaches. The kernel itself
finds its faces by their physical heights on every mesh.

Diagnostics
-----------

At the first step and then every ``diagnostics_int`` steps each span appends
a row to ``<output_root>.dat`` for a single span, or
``<output_root>_span<k>.dat`` for span ``k`` of a section (``k`` from 1 at
``end_a``; whitespace separated, like ERF's data logs): the time,
the middle node's position, its sag below the chord, its lateral offset from
the chord's vertical plane, the swing angle in degrees, the tension at the
first and last node, the largest tension on the line, the wind handed to
MoorDyn at the middle node, the smallest clearance of any node above the
terrain under it with that node's horizontal position, and the air's total
drag on the span (columns ``time mid_x mid_y mid_z mid_sag mid_offset
swing_deg tension_a tension_b max_tension mid_u mid_v mid_w min_clearance
min_clearance_x min_clearance_y drag_x drag_y drag_z``). All vectors are in
ERF's frame (z up), in m, m/s or N: ``drag_x``, ``drag_y``, ``drag_z`` are the
force of the air on the span, ``mid_sag`` is positive below the chord, and
``mid_offset`` and ``swing_deg`` are positive to the left of the direction
from ``end_a`` to ``end_b``. A line on strings also writes
``<output_root>_insulators.dat``: for each tower ``t<j>`` (``j`` from 1 at the
tower next to ``end_a``) the string's swing from the vertical and across the
line in degrees and the tension at its top (N). The clearance uses the
terrain surface under each node (on a fitted mesh the ``k = 0`` node plane
from ``terrain_heights``; with an immersed terrain the surface the immersed
boundary is built from), so a blown-out span on a slope is measured against
the ground it hangs over at that moment.

From ``stats_start`` on each span keeps running statistics, the mean, root
mean square, minimum and maximum of its swing angle (``swing_deg``), mid-span
offset (``mid_offset``, m), tensions at the first and last node
(``tension_a``, ``tension_b``, N), largest tension (``max_tension``, N),
minimum clearance (``min_clearance``, m) and the air's drag across the span
(``drag_normal``, N: its horizontal part normal to the span's horizontal chord,
positive along :math:`\hat z \times` the span's direction (towards ``end_b``);
the y component for a span along +x: a checkpoint holding ``drag_y``, the
ERF-frame y drag, in its place continues the span's statistics on a span along
+x and starts all of them afresh, saying so, on any other), in
``<output_root>_stats.csv`` (``<output_root>_span<k>_stats.csv`` in a
section), and the strings theirs in ``<output_root>_insulators_stats.csv``;
these are the numbers a turbulent-wind run is judged on. With
``node_output_int`` every node of every line (the ``node`` column counts
from 0: the spans' nodes from ``end_a``, then each string's from its top;
position, clearance, tension, wind and drag) is written to
``<output_root>_nodes.dat`` for plotting the line's shape.
``<diagnostics_dir>/total_load.dat`` holds, every ``diagnostics_int`` steps,
the air's drag on all lines and all lattice towers (``drag_x``, ``drag_y``,
``drag_z``, N, ERF frame), the force they put into the air (``force_on_air_x``
to ``force_on_air_z``: minus the drag with ``drag_on_flow``, zero without) and
the integral of the momentum source (``source_x`` to ``source_z``, N), which
equals that force: the comparison checks that the spreading keeps the force.
The force is the lines' drag at the end of the ERF step and the towers' members'
drag of its last coupling step, held over the next, so over a step the air takes
that drag times the step rather than the drag the lines and members felt
integrated over MoorDyn's calls and the coupling steps; over a run the two
differ by half a step times the change of the drag from the start to the end. The MoorDyn input file the line was built from is kept
next to it
(``<diagnostics_dir>/<name>.moordyn.txt``) for inspection or for running
MoorDyn on its own.

Clearance between lines
-----------------------

Every pair of lines is watched for how close their conductors come: the
closest approach of the two polylines through their span nodes (the
strings left out, since a string belongs to its own phase), found exactly
from the closest points of every pair of segments. Every
``diagnostics_int`` steps ``<diagnostics_dir>/separation.dat`` gets, for
each pair ``A-B``, that distance, the point midway between the closest
points and a clash flag, 1 when the distance is below
``flashover_distance``, the gap an arc is assumed to jump. The flashover
distance depends on the voltage and the line's insulation coordination; the
default of 1 m is a placeholder to be set for the line studied. From
``stats_start`` each pair keeps the mean, root mean square, minimum and
maximum of the distance and of the flag, whose mean is the fraction of the
time spent in a clash, in ``<diagnostics_dir>/separation_A-B_stats.csv``.
The start-up log prints how far apart each pair hangs in still air and
warns when a pair starts inside the flashover distance. Two identical
parallel lines are almost equally far apart along their whole length, so
where they come closest can move by a span's length on a millimetre's
difference;
the distance itself is well defined.

Transformers
------------

A transformer is a box on the terrain on which lines are dead-ended; ERF
reports the lines' load on it and how close any conductor comes to it. A
transformer (``conductors.transformers``, a block of its own per name) is
a box of ``size`` (length along ``x``, width along ``y``, height) standing on
the terrain: its base is at the terrain height under the centre of its
footprint, ``position``. Every line end, ``end_a`` or ``end_b``, whose
``(x, y)`` lies on a footprint (edges included) is dead-ended on that
transformer, at the end's own height above the terrain, which must clear the
box's top; the MoorDyn line is unchanged, its end still a fixed point. Two
transformers may not overlap under an end, and a transformer on which no line
ends is refused.

A dead end takes the pull of its line: the net force MoorDyn finds on the
line's end node, which the fixed point holds still, that is the end
segment's tension with the node's share of the weight and the drag. The
load on a transformer is the sum :math:`\mathbf{F}` of those pulls and
their moment about the centre of its base,

.. math::

   \mathbf{M} = \sum_e (\mathbf{x}_e - \mathbf{x}_b) \times \mathbf{F}_e ,

where :math:`\mathbf{x}_e` is the dead end and :math:`\mathbf{x}_b` the base
centre. The horizontal force :math:`F_h = |(F_x, F_y)|` and the overturning
moment :math:`M_h = |(M_x, M_y)|`, the one that tips the box, are checked
against ``allowable_force`` and ``allowable_moment``; the flag is up when
either is exceeded, and an allowable of 0 (the default) is not checked.
Lines pulling from different sides partly cancel, so the check is on the
vector sum, which is what the foundation and the bushings carry. In still air
two lines ending on a transformer from opposite sides balance; in a crosswind
they do not.

Every conductor is also watched for how close it comes to each box: the
closest approach of the polyline through its span nodes to the box, found by
a golden-section search along each segment (the distance from a point moving
along a segment to a box is convex). The lines ending on a transformer count
too: their ends clear its top by their standoff, and a conductor dropping
steeply from its bushing towards a valley comes closer. The clearance is
flagged against ``flashover_distance`` like a pair of lines.

Every ``diagnostics_int`` steps ``<diagnostics_dir>/transformers.dat`` gets,
for each transformer ``T``, the force components ``T_Fx``, ``T_Fy``,
``T_Fz``, ``T_Fh``, the moment ``T_Mx``, ``T_My``, ``T_Mh``, the allowable
flag ``T_over``, the clearance ``T_clearance`` and its flag ``T_clash``
(``T_F*`` is the force of the lines on the transformer, N, and ``T_M*`` its
moment about the centre of the base, N m, both in ERF's frame); from
``stats_start`` the horizontal force, the overturning moment, both flags and
the clearance are kept as running statistics in
``<diagnostics_dir>/transformer_T_stats.csv``. The start-up log prints, for
each transformer, the height of its base, the line ends on it, the still-air
load and the closest conductor.

Towers
------

The suspension towers of a line with a ``tower_type`` are lattice towers
loaded by the wind (``conductors.tower_types``, a block of its own per
type). A tower stands on the terrain under its suspension point: a square
lattice body tapering linearly from ``base_width`` at the ground to
``top_width`` at the cross-arm (with an optional ``peak`` above it at the top
width), and a lattice cross-arm of ``arm_length`` centred on the body at the
height the conductor hangs from, across the line (normal to the mean
horizontal direction of the spans either side), whose face is ``arm_depth``
deep out to the body's corners at the cross-arm's bottom and tapers linearly to
nothing at the tips: a truss of chords from the body's corners at the
cross-arm's top and bottom converging on each tip, the shape of the generated
frame's cross-arm (:doc:`TowerFrame`), with about half the area of a rectangle
of the same depth. Its
drag nodes stand at the centroids of their pieces of that face and carry each
piece's mean depth, so the drag and its moment on the face are the exact
integrals. The towers carry the wind's drag on their members and the pull of the
line hanging from them, and stand on a foundation of four footings; they are
rigid unless their type has a ``frequency`` (see `Moving towers`_) or a frame
(``frame_file`` or ``frame_panels``; see :doc:`TowerFrame`).

A circuit's phases and its shield wire hang from one row of towers: a line
with ``share_towers = <line>`` hangs from the towers of that line, which has
the ``tower_type``, each at its own point on them, its ``towers`` points. The
phases sit across the cross-arm and the shield wire on the ``peak`` above it;
a point off the tower (further from its axis than half the cross-arm, below
its base, or above the top of its peak) stops the run naming the line and the
tower. A shared tower stands
under the owning line's point, and the other lines' points are placed above
its base, not above the ground under each point, so that a cross-arm on a slope
stays level. The tower takes each line's pull at its own point, the
foundation all of them, and the lines that share a moving tower step with it
together.

Each member is cut into drag nodes at the middle of equal segments:
``segments`` up the body to the cross-arm, as many on the ``peak`` as fit at
that spacing (at least one), and four on the arm: two on each of its parts
outside the body, from the body's half width at the middle of the face's depth
at the body (half ``arm_depth`` below the cross-arm) to the tip, each at the
centroid of its piece of the tapering face with the piece's mean depth as its
width, so that only the body's nodes carry the wind where the arm passes through
it (``arm_outside_shaft``, the default, which needs ``arm_length`` longer than
the body is wide there; ``false`` is the earlier layout: the four along the
whole arm at the cross-arm's height, ``arm_depth`` deep all along, which counts
the part within the body twice). A node stands for a
length :math:`\ell` (m) of member with axis :math:`\mathbf{e}` (a unit
vector), and the wind :math:`\mathbf{U}` (m/s) loads it as a slender member
does, by the flow normal to it, relative to the node's own velocity
:math:`\mathbf{v}` (m/s, zero for a rigid tower):

.. math::

   \mathbf{F} = \tfrac{1}{2} \rho\, C_f\, b\, \ell\, |\mathbf{U}_n| \mathbf{U}_n ,
   \qquad \mathbf{U}_n = (\mathbf{U} - \mathbf{v}) - \big((\mathbf{U} - \mathbf{v})\cdot\mathbf{e}\big)\mathbf{e} ,

with :math:`\mathbf{F}` (N) the force of the air on the node, :math:`\rho`
the ``air_density``, :math:`b` (m) the face width times the ``solidity``
:math:`\sigma` (the members' area over the face's outline) and :math:`C_f`
the ``drag_coefficient`` or, by default, the force coefficient of a square
lattice tower of flat-sided members on the projected area of one face,
:math:`C_f = 4\sigma^2 - 5.9\sigma + 4` (ASCE 7, the load standard of the
American Society of Civil Engineers), which counts the windward and the
leeward face together: 2.98 at :math:`\sigma = 0.2`. Wind along a diagonal of
the square body meets more of its members than wind normal to a face, so the
body's and the peak's nodes take the factor :math:`1 + g \sin^2 2\theta`, with
:math:`\theta` the angle of :math:`\mathbf{U}_n` from a face's normal and
:math:`g = \min(0.75\sigma, 0.2)`: along the diagonal ASCE 7's
:math:`1 + 0.75\sigma` (1.15 at :math:`\sigma = 0.2`), which ASCE 7 caps at 1.2,
and in between the angular shape of IEC 60826's :math:`1 + 0.2 \sin^2 2\theta`
(``diagonal_wind_factor``, the default; ``false``: no factor). The cross-arm
takes no such factor. A rotor-less tower in
AeroDyn is the same law with :math:`b` the tower's diameter. The wind is
ERF's velocity sampled at the nodes with the lines' sampler at the start of
each step, or the prescribed velocity; a node just above sloping ground,
below the averaged bottom face of its cell, is read from the bottom cell. The
aerodynamics sit behind a narrow interface (node positions, axes and
velocities and the wind in, a force per node out), so another model of the
members' loads can stand in for the ERF drag.

The line pulls on the tower where it hangs from it: the net force MoorDyn
finds on the top node of the insulator string there, or, for a conductor
clamped to the tower, on the end nodes of the two spans meeting there. In
still air a tower of a level section carries the weight of one span (half of
each span either side) and its string, the spans' pulls along the line
cancelling; in a crosswind it also takes the wind on that length, the wind
span. Drag and pull together load the foundation: the base shear (the
horizontal force), the overturning moment about the centre of the base and
the downward load :math:`P` (N), the tower's ``weight`` included. The four
legs stand at :math:`(x_i, y_i) = (\pm s/2, \pm s/2)` in the frame of the line
and the cross-arm, :math:`s` (m) the ``leg_spacing`` (the base width by
default), and, as a rigid square, share :math:`P` equally and the moment
:math:`\mathbf{M}` (N m) of all loads about the centre of the base linearly,

.. math::

   N_i = \frac{P}{4} + \frac{M_y x_i - M_x y_i}{4 (s/2)^2} ,

:math:`N_i` (N) the axial force of leg :math:`i`, compression positive, so
that the reactions balance the loads. The largest compression,
:math:`\max(0, \max_i N_i)`, and the largest uplift,
:math:`\max(0, -\min_i N_i)`, reported as a positive force (0 when no leg
pulls), are flagged against ``allowable_compression`` and ``allowable_uplift``, a
footing's bearing and pull-out capacity (0, the default, is not checked). A
pull along the diagonal loads the corner leg :math:`\sqrt{2}` times more
than the same pull face-on.

Every ``diagnostics_int`` steps ``<diagnostics_dir>/towers.dat`` gets, for
each tower ``<line>_t<k>``, the drag ``drag_Fx``, ``drag_Fy``, ``drag_Fz``,
the line's pull ``line_Fx``, ``line_Fy``, ``line_Fz``, the ``shear``, the
``overturning`` moment, the ``vertical`` load, the ``max_compression`` and
``max_uplift`` of the legs and the allowable flag ``over``, at the time the
step starts (the flow the drag was found from and the line as it was then);
from ``stats_start`` each tower keeps the statistics of its horizontal drag
and pull, shear, overturning moment, leg compression and uplift and the flag
in ``<diagnostics_dir>/tower_<line>_t<k>_stats.csv`` (``k`` from 1 at the
tower next to ``end_a``). In ``towers.dat`` all vectors are in ERF's frame,
in N and N m: ``drag_F*`` is the force of the air on the tower's members,
``line_F*`` the force of the lines on the tower, ``vertical`` is positive
down and ``max_uplift`` is a positive tension. The start-up log prints
each tower's still-air pull and leg loads. With ``drag_on_flow`` the
towers' drag goes into the momentum sources with the lines', and
``total_load.dat`` counts both. A checkpoint carries the statistics and the
last node forces, so the drag a restart puts into the flow is the
checkpointed step's.

Verification (unit tests ``MemberDrag``, ``Tower``): a node takes only the
normal component of the relative wind (none along its axis, none moving with
the wind); in a uniform wind a 30 m tower of 6 m base, 1.5 m top width and
solidity 0.2 carries exactly :math:`p\, C_f \sigma (b_0 + b_t) z_a / 2` on its
body, :math:`p = \rho U^2 / 2` the dynamic pressure, :math:`b_0` and
:math:`b_t` the base and top widths and :math:`z_a` the cross-arm height (the
width is linear, so the segments' midpoint sum is exact), and
:math:`p\, C_f \sigma\, d\, [(s_0 - s_i) + (L_a/2 - s_0)/2]` on each part of its
cross-arm outside the body, :math:`d` the ``arm_depth``, :math:`L_a` the
``arm_length``, :math:`s_i` the body's half width at the middle of the face's
depth and :math:`s_0` its half width at the cross-arm's bottom, where the
taper starts (the closed integral of the face: for the test tower 3.15 m\ :sup:`2`
a side, 50.4 % of the rectangle the model took before); the base moment equals the hand
integral up to the midpoint sum's known :math:`z_a^3/(12 n^2)` on the body's
quadratic part, :math:`n` the ``segments``, the cross-arm's part exact; wind along the cross-arm loads only the body; in a log-law
wind ten segments give the drag and base moment of a 20 000-point integral
within 1 %. The legs' reactions are in equilibrium with the loads for any
orientation, face-on and diagonal pulls give the hand values, and each
allowable flags on its own. Against the real MoorDyn, the middle tower of four
300 m spans on 2.5 m strings carries 5403.6 N down in still air against the
weight span's 5400.1 N; in a steady 20 m/s crosswind it takes 2149 N across
the line against :math:`qL` plus the swung string's drag, 2154 N, and 5343 N
down against the weight span less the string drag's lift, 5350 N. (The string,
swung :math:`\theta` across the line, sees the normal wind
:math:`U\cos\theta`, along :math:`(\cos\theta, \sin\theta)` in the plane
across the line, so its drag :math:`q_i L_i \cos^2\theta` lifts as well as
pushes.) A tower next to a dead end differs: the span from the dead end
drops to the string's bottom, and the lower end takes less of its weight. On
the terrain test case the hilltop tower carries 20.9 kN of drag, 22.1 kN with
the rectangular cross-arm of earlier versions: the difference is the drag on
half the cross-arm's face outside the body, which the taper removes.

Moving towers
-------------

A tower type with a ``frequency`` (Hz, its first bending frequency on a rigid
foundation) bends under its loads, and its line hangs from a cross-arm that
moves: MoorDyn takes the cross-arms as coupled points, which ERF moves and on
which MoorDyn hands back the line's pull, the way OpenFAST couples MoorDyn to
a platform. The tower is a single mode, the same in :math:`x` and :math:`y`
(a square tower), on a foundation that can tilt and slide
(``foundation_rotational_stiffness`` :math:`k_r` in N m/rad and
``foundation_lateral_stiffness`` :math:`k_l` in N/m, rigid when 0, the
default), with the structural ``damping_ratio`` :math:`\zeta` (0.02 by
default). Its mass, ``weight`` over gravity, is spread by length of member,
:math:`m_i` (kg) at node :math:`i`: the shaft's over its drag nodes, and the
whole cross-arm's (its full ``arm_length``) at the cross-arm's height, shared by
the cross-arm's drag nodes, so where ``arm_outside_shaft`` puts those nodes does
not change the mode's mass or stiffness. A
node at height :math:`z` (m) above the base moves horizontally by
:math:`\psi(z)` times the cross-arm's displacement :math:`x_a` (m),

.. math::

   \psi(z) = a_b \left(\frac{z}{z_a}\right)^2 + a_r \frac{z}{z_a} + a_l ,
   \qquad a_b = \frac{1/K_b}{C},\; a_r = \frac{z_a^2/k_r}{C},\; a_l = \frac{1/k_l}{C},\;
   C = \frac{1}{K_b} + \frac{z_a^2}{k_r} + \frac{1}{k_l} ,

the cantilever's bending, the footing's tilt and its slide, each in
proportion to its compliance under a load at the cross-arm, :math:`z_a` (m)
the cross-arm's height above the base and :math:`C` (m/N) the total
compliance (a term with a stiffness of 0, a rigid foundation, is left out).
The generalized mass is :math:`M_g = \sum m_i \psi(\hat z_i)^2` (kg) and the stiffness
:math:`K_g = 1/C` (N/m), the strain energy of that shape exactly, with the
bending stiffness :math:`K_b = (2\pi f)^2 \sum m_i (\hat z_i/z_a)^4` (N/m) that
gives the type's frequency :math:`f` (Hz) on a rigid foundation, :math:`\hat z_i`
the height of node :math:`i`'s mass: its own height on the body, and
:math:`z_a` for the cross-arm's nodes, whose mass is the cross-arm's at its
height (their drag acts at their own height, :math:`\psi_i = \psi(z_i)` in
:math:`Q` below, and their inertia reaches the foundation from :math:`z_a`); a foundation
that gives lowers the frequency to :math:`\sqrt{K_g/M_g}/2\pi`. In each
horizontal direction

.. math::

   M_g \ddot x_a + 2 \zeta \omega M_g \dot x_a + K_g x_a = Q ,
   \qquad Q = \sum_i \psi_i F_i + \sum_a \psi_a F_a ,
   \qquad \omega^2 = K_g/M_g ,

with :math:`F_i` (N) the members' drag and :math:`F_a` (N) the pull of the
line hanging at attachment :math:`a`, weighted by the shape at the
attachment's height (:math:`\psi_a = 1` at the cross-arm, above 1 for a shield
wire on the peak). A load held over a step is
integrated exactly, so any step is stable. The tower does not twist, rise or
sink, and its weight does not add to the overturning as it leans
(:math:`P`-:math:`\Delta`).

A line whose towers move steps in coupling steps, at least ``substeps`` and
at least 20 over the period of its fastest tower; lines that share towers step
together. In each, every tower advances under its members' drag and the lines'
mean pull over the coupling step, both held over the step: a one-mode tower by
its exact solution for a held load, a frame tower by the implicit midpoint rule
(Newmark's method started from the acceleration the held loads give). Then
MoorDyn moves each line's points on the towers from where they were to where
the towers have taken them, at a constant velocity over the coupling step
(MoorDyn moves a coupled point linearly), in one call per internal step of
MoorDyn's, and the lines' pull on each tower is read after every call and
averaged over the step: MoorDyn's second-order Runge-Kutta step leaves the
forces at its midpoint, where it applies them over the whole step. A tower whose
attachment moves by :math:`\Delta\mathbf x` under a held pull
:math:`\bar{\mathbf F}` takes the work :math:`\bar{\mathbf F} \cdot
\Delta\mathbf x`; with :math:`\bar{\mathbf F}` the lines' mean pull over the
step, that is the work the lines give up, since the attachment moves at the
constant velocity :math:`\Delta\mathbf x / \Delta t_c`. The exchange then
neither makes nor loses energy, to the coupling tolerance, so it gains none from
motions the coupling step :math:`\Delta t_c` does not resolve, and the coupling
step need only resolve the towers' sway. The fastest motions at a
tower are not the tower's: an insulator string (axial stiffness
:math:`10^7` N over 2.5 m) with its own mass and the conductor's near its foot
rings at about 45 Hz against a frame's stiff cross-arm, far above what 20
coupling steps a period of a 4 Hz tower resolve. Exchanging the pull at the
step's end instead, with a frame stepping by Newmark's method to it, lets the
interface gain energy from that ring when the coupling step is near half its
period. ``Conductors_FrameTowers_CouplingStep_MoorDyn`` runs
``Conductors_FrameTowers`` (0.9 s steps, 76 coupling steps chosen) with 76 and
with 160 coupling steps, and the middle span's tensions must agree to three
digits. The
mean pull is not known before MoorDyn's calls, so the coupling step is
iterated: the lines' MoorDyn states and the towers' are kept at its start and
restored for each iteration (MoorDyn's in-memory serialization), and the
estimate of the mean pull, from the pull at the step's start, is updated with
Aitken's relaxation until it changes by less than :math:`10^{-4}` of the
largest pull (at most 50 iterations). A single exchange, the towers stepping
with the pull from the end of the last call, is only conditionally stable: a
short, nearly taut span pulls back in milliseconds, and a pull that lags acts
on the tower as negative damping, :math:`k \Delta t_c / 2` (N s/m) for a span
of stiffness :math:`k` (N/m) on a tower holding the start pull, which outgrows
the tower's own damping once :math:`k` exceeds a fraction of the tower's
stiffness (about a quarter at a damping ratio of 0.02 and 20 coupling steps a
period). The
members' drag is held from the start of each coupling step, a damping force a
step late that an undamped frame's modes above about 6 times its first
frequency would grow under: a frame's ``damping_ratio`` must be positive.
Iterated, the coupling is the implicit one; it takes about four iterations
where the towers carry strings and a clamped shield wire, and one once they are
still. ``<diagnostics_dir>/coupling.dat`` logs, every ``diagnostics_int``
steps, the most iterations a coupling step took since its last row and how many
did not converge (the run warns once and goes on with the last iterate; in a
group of lines sharing towers, 20 unconverged coupling steps in a row, and at
least two ERF steps' worth, stop it, the lines and towers diverging; a checkpoint
carries the counts of coupling.dat's next row). A
pull that is not finite stops the run instead, naming the line and the time. The members' drag is
found again in every coupling step from the wind of the ERF step and the
members' current velocity, so the wind damps the sway: a member moving with the
wind feels less of it. The foundation takes the loads
less the inertia of the nodes (each node's mass times its acceleration), so
a tower swaying freely still loads its footings; both are those at the
coupling step's end, the acceleration under the lines' mean pull over that
step. MoorDyn's pull on a coupled
point leaves out the inertia of the line's end node, a few kilograms against
the tower's tonnes.

``towers.dat`` gets the cross-arm's displacement ``arm_dx`` and ``arm_dy``
for each moving tower (on a frame tower, that of the cross-arm's centre on the
tower's axis, which a twist of the frame leaves in place), and its statistics the size of that displacement,
``arm_displacement``. The start-up log prints each moving tower's frequency
on its foundation, generalized mass and stiffness. A tower starts upright and
at rest, so a run that starts in a wind rings down from the sudden load, over
about :math:`1/(\zeta\omega)`: 4 s at 2 Hz and 2 %. A checkpoint carries each
tower's sway, and the line restarts with its cross-arms where the towers had
taken them. The tower model sits behind a narrow interface (the node loads
and the line's pull in, the motion of every node and of the cross-arm and the
nodes' inertia out), so a frame model of the lattice can stand in its place.

Verification (unit tests ``OneModeTower``, ``ConductorLine``,
``Conductors``): a load held over a step is integrated exactly however the
time is cut, and a step of a thousand periods lands on the static
deflection; released, the tower rings down with the damped period to
:math:`10^{-4}` and the logarithmic decrement of its damping ratio to
:math:`10^{-3}`; a load :math:`F` at the cross-arm deflects it by
:math:`F C` with the foundation's compliances in series, and a load
:math:`p_l (z/z_a)^2` (N/m) spread up the body by the hand value of its
midpoint sum; swaying freely, its base shear is
:math:`\omega^2 x_a \sum m_i \psi_i`. Swaying in a steady
20 m/s wind, a 30 m tower at 2 Hz with 0.5 % structural damping rings down at
0.9047 % against the 0.5 % plus the 0.4047 % of the linearised relative drag,
:math:`\sum \psi_i^2 \rho C_f b_i \ell_i U / (2 M_g \omega)`. Against both
MoorDyn libraries, the pull MoorDyn hands back on a coupled point is the one
the tower is reported to carry, a coupled step moves the cross-arm by its
velocity times the step, and a coupled point held still is a fixed one to
:math:`10^{-9}` m. In a 10 m/s crosswind the towers of a clamped section
settle 2.6 mm downwind, where their stiffness balances the 4023 N of drag on
the body (weighted by the shape) and the line's pull of 672 N; that pull is
within 2 % of the rigid towers'. A circuit's shield wire clamped 40 m from
its dead end to the peak of a tower at 2 Hz with 2 % damping, in the same wind,
settles at 4 mm when the coupling steps are iterated; with a single exchange per
coupling step the tower is still growing through 18 mm after 15 s. On the moving-towers test case the towers
bend at 1.67 Hz (2 Hz lowered by footings of :math:`10^9` N m/rad) and the
hilltop tower leans about 3 cm in its 25 m/s wind.

A tower type that gives ``conductors.<type>.frame_file`` bends as a frame
model instead of one mode: the SubDyn input file of its lattice, in tower-local
axes, gives every member's stiffness and mass (:ref:`sec:TowerFrame`, section
"Coupling to the conductor lines"). The drag nodes and the attachments are tied
to the frame's nodes, so the coupling above runs unchanged, and the footings'
loads come from the frame's support reactions. On the frame-towers test case
(the moving-towers lines on a 76-member frame with a first natural frequency of
4.22 Hz) the coupling converges in three or four iterations per coupling step
with the real library, the hilltop tower carries a base shear of 21.4 kN against
the one-mode tower's 21.7 kN, and its cross-arm moves 3.4 mm against the
one-mode tower's 26 mm, the frame being stiffer.

A tower type can instead have ERF generate its frame from its dimensions
(``conductors.<type>.frame_panels`` with the angles of its legs and its
bracing; :ref:`sec:TowerFrame`, section "Generated lattice towers"), and every
member of a generated frame, or of a frame file given with ``member_file``, is
checked against its strength in tension and compression (ASCE 10-15) every step;
``towers.dat`` carries each tower's largest utilisation and the member it is in.
``steel_temperature`` sets the steel's temperature, which lowers its stiffness
and strength (EN 1993-1-2). On the generated-towers test case (the moving-towers
lines on generated towers of 8 panels, 310 members and 11.1 t) the hilltop
tower's most utilised member, a leg at the base, reaches 0.12 of its design
strength at 20 C, 0.23 at 500 C and 0.66 at 650 C, while the tower's first
natural frequency falls from 3.31 Hz to 2.56 and 1.55 Hz and its cross-arm's
largest sway grows from 5 mm to 12 and 36 mm; the coupling converges in three
iterations per coupling step with the real library.

A network over hills
--------------------

``Structures/Exec/PowerLines`` puts the pieces together: ``make_case.py``
draws four hills of 87 to 99 m on a 3 km by 2 km domain, six transformers
(three on hilltops, three on flat ground) and the five lines of a minimum
spanning tree between them, each a section of three spans on insulator
strings with lattice towers placed for ground clearance, all strung to 20 kN, and a
neutral log-law inflow of 18 m/s at 30 m. A tower in a dip, which the spans
either side would pull up, is raised until the line weighs on it: four of the
ten stand 36 to 46 m tall rather than 30 m. The k-equation RANS
(Reynolds-averaged) flow spins up for 600 s (1.2 million cells), then the
lines run 120 s from its checkpoint with the real MoorDyn. The numbers below
come from one run of the committed decks against MoorDyn-C 2.7.1; no test
checks them. The wind 30 m above the ground speeds up to about
25 m/s over the hilltops and slows in their lee. The lines settle within
about 30 s, every span at 20 to 22 kN: the spans running across the wind blow
out up to 2.3 m at mid-span (swings of 12 to 22 degrees), the one running along
it barely moves (2 degrees), and the transformers' horizontal loads settle
between 10.8 kN (three lines from different sides) and 22.1 kN, under their
25 kN allowable. The lattice towers (6 m base, 1.5 m top, solidity 0.2, 12 m
cross-arm) carry 18.7 to 20.6 kN of wind drag on the hilltops, 10.5 to 15 kN on
the raised towers and about 7 kN on the 30 m towers on flat ground (with the
rectangular cross-arm of earlier versions; the tapering cross-arm lowers a
tower's drag by 5 to 7 %, as the regression decks show). With the
line's pull and a 60 kN tower weight on 6 m footings, the towers on the hills
and the raised ones lift a footing by 11 to 22 kN (the allowable is 50 kN),
their worst legs carrying 42 to 53 kN in compression against 29 kN on flat
ground. The towers bend at 1.44 to 1.67 Hz: 2 Hz on a rigid foundation,
lowered by footings of :math:`10^9` N m/rad the more the taller the tower
(1.67 Hz at 30 m, 1.44 Hz at 46 m, from the one-mode model above; the start-up
log prints each tower's frequency). They settle leaning 20 to 32 mm downwind,
10 mm on flat ground. With a
steady RANS wind the lines hold a steady blowout; their gust response needs a
turbulent inflow.

The same lines in a turbulent wind (``les/`` in the same directory): a periodic
precursor over flat land (Deardorff large-eddy simulation, LES, 16 m cells,
3072 by 1536 by 768 m,
:math:`z_0` = 0.1 m, :math:`u_*` about 1 m/s) spins up for 7200 s and writes
boundary planes every 1.5 s; a run over three immersed hills restarts from its
checkpoint and takes its inflow from the planes, with three circuits (three
phases and a shield wire each, 12 lines) on five shared towers that bend. The
numbers below come from one run of the committed decks; no test checks them.
Over 1500 s of statistics the wind at the conductors' middles peaks about 1.3 times
its mean (16.6 and 21.4 m/s on the line across the wind); that line swings 15
degrees on average and 22 at its peaks, the lines along the wind 3 and 6; the
conductors' tension moves by about 1 %, the strings taking up the swing; the
towers' sway and base shear peak at 1.4 to 2 times their means, up to 35 mm and
20 kN; and the footings' uplift varies most, a tower lifting 9 kN on average
reaching 33 kN in a gust. Each coupling step converges in at most four
iterations.

ASCE 74 design check
--------------------

This section describes the quasi-static wire wind load of ASCE Manual of
Practice 74 (Guidelines for Electrical Transmission Line Structural Loading)
that ERF can compute for every span (``ERF_ASCE74.H``), the design value a run's
lines are compared with. ``conductors.asce74_wind`` gives the 3-second gust
:math:`V` at 10 m over open country and ``conductors.asce74_exposure`` the
terrain exposure, B (suburban or wooded) or C (open country, the default).

The wind load per metre of conductor is

.. math::

   F = Q\, k_z V^2 G_w C_f d,

with :math:`Q = \rho/2` (``air_density``; ASCE 74's 0.613 for 1.226 kg/m\ :sup:`3`),
:math:`C_f` the line's ``drag_coefficient`` and :math:`d` its ``diameter``. The
exposure factor :math:`k_z = 2.01 (z/z_g)^{2/\alpha}` raises the gust with the wire's
height :math:`z`, so that :math:`k_z V^2` is the square of the 3-second gust at that
height. The wire's gust response factor

.. math::

   G_w = \frac{1 + 2.7 E \sqrt{B_w}}{k_v^2}, \qquad
   E = 4.9 \sqrt{\kappa} \left(\frac{33\ \mathrm{ft}}{z}\right)^{1/\alpha}, \qquad
   B_w = \frac{1}{1 + 0.8 L / L_s},

with :math:`k_v = 1.43`, the 3-second gust over the 10-minute mean at 33 ft over
open country, accounts for the gust not reaching the whole span :math:`L` at once
(:math:`G_w` falls from about 0.8 for a 100 m span to about 0.6 for 700 m). The
constants, from ASCE 74's third edition:

========  ========  ==================  =================  ==================
Exposure  alpha     z_g                 kappa              L_s
========  ========  ==================  =================  ==================
B         7.0       1200 ft (365.8 m)   0.010              170 ft (51.8 m)
C         9.5       900 ft (274.3 m)    0.005              220 ft (67.1 m)
========  ========  ==================  =================  ==================

Under :math:`F` and its net weight :math:`W` per metre the span swings about its
chord by :math:`\arctan(F/(W\cos\beta))`, with :math:`\beta` the chord's angle to
the horizontal: on an inclined chord the tension carries the weight's part along
it, and only :math:`W\cos\beta` opposes the wind (``asce74_inclined_spans``, the
default; ``false`` takes the span as level under the whole weight: the swing
:math:`\arctan(F/W)`, and the sag and tension under :math:`\sqrt{F^2 + W^2}`
without the upper end's share below; the still-air sag of the effective height
below is the inclined chord's either way). It hangs as the elastic catenary (above)
under :math:`\sqrt{F^2 + (W\cos\beta)^2}` across the chord, as a level span of the
chord's length: its sag in the swung plane and its mid-span's sideways blowout,
the sag times the sine of the swing. The tension is the upper end's: the two
ends' tensions differ by the weight times the rise, :math:`W c \sin\beta`, since
the wind does no work along the chord, and the level span's lies half way
between them, so the upper end carries the level span's plus
:math:`W c \sin\beta / 2`: against the exact inclined elastic catenary,
within -0.7 to +0.02 % for chords of 50 to 800 m up to 30 degrees, slacks up to
1 % and axial stiffnesses of :math:`10^7` N or more, falling to 1.7 % low at 3 %
slack in still air (the level span's alone is up to 14 % low at 10 to 30
degrees, more on the slacker spans). Weight and wind are both taken per
unstretched metre; the wind per stretched metre would raise the tension by
:math:`\varepsilon F^2/(F^2 + W^2)` for a strain :math:`\varepsilon`, under
0.1 % at strains below :math:`10^{-3}`. A span's height is the
wire's effective height: the mean of its two attachment points' heights above the
ground under them (the bottoms of the insulator strings at the towers), less a
third of its sag in still air, measured vertically (``asce74_wire_height =
effective``, the default; ``attachment``: the mean alone, the higher and more
conservative height). Its length is its chord.
At start-up ERF writes ``<diagnostics_dir>/asce74.csv``: per span its height,
chord (between the points the conductor hangs from, the bottoms of the
insulator strings at the towers, the chord its unstretched length spans),
unstretched length, :math:`k_z`, :math:`G_w`, the pressure
:math:`Q k_z V^2`, the load, the weight, the swing (degrees), the sag, the
blowout and the tension.

With these constants the product :math:`k_z G_w` reproduces the ASCE wire load
ratios that Peyrot (Wind loading: uncertainties and honesty suggest
simplification, ASCE Electrical Transmission and Substation Structures
Conference, 2009, Tables 1 and 2) gives for a Drake conductor at five
attachment heights and spans: to the table's two digits for a 100 m span 10 m
up (0.80 over open country, 0.64 suburban), and within 0.03 for spans of 350 to
720 m at 20 to 50 m, the paper not stating the height it took for a sagging
span (unit test ``ASCE74.WireLoadRatiosArePeyrots``).

``Structures/Exec/PowerLines/les/compare_asce74.py`` compares a run with it:
for every span, the LES wind at mid-span normal to the span gives its mean and
its peak 3-second average :math:`V_3` (the gust at that height, so that
:math:`k_z V^2 = V_3^2`); ASCE 74's predicted peak load is then
:math:`(\rho/2) C_f d V_3^2 G_w`, against the LES's peak span load (the span's drag
normal to it over its chord), and the swing and the tension under it against
the LES's peaks. The LES's own gust response factor, its peak span load over
:math:`(\rho/2) C_f d V_3^2`, measures how much of the point gust the whole span
feels, which is what :math:`G_w` models.

On the hills' LES of the previous section, with the middle phase of each
circuit clamped at fixed points where the towers stood (eight spans of 36 to
265 m, 20 to 38 m up, the check at their effective heights of 20 to 36 m), over 600 s of statistics: the 3-second gust at mid-span
is 1.2 to 1.5 times the mean wind across the span. On the spans the wind
crosses (16 m/s across the 213 m span, 9 to 11 m/s across two of about 250 m) the
whole span feels 0.74 to 0.83 of the point gust's load where :math:`G_w` gives
0.70, and 0.97 against 0.83 on a 46 m span. That compares :math:`G_w` with the
LES's own 3-second gust, 1.20 to 1.37 times the mean on these spans, where
:math:`G_w` assumes ASCE 74's :math:`k_v = 1.43` (and a 16 m LES resolves only
part of the 3-second gust: its smallest eddies are some 6 to 8 cells, 100 m and
more). Taken over the mean span drag (the LES's peak span load over its mean),
ASCE 74's :math:`k_v^2 G_w` (1.43 for these spans) is 10 to 13 % above the LES's
1.27 to 1.30, and 19 % above on the 46 m span. Taken over the load of the mean
wind, :math:`q(\bar V)`, which is ASCE 74's own basis (the mean drag also holds
the fluctuations' and the yaw's share), the LES's peak is 1.24 to 1.56 times it
on the three crossed spans, and ASCE 74 is 10 % and 15 % above on two of them
and 9 % short on the third (the 254 m span). The peak
tensions (the upper end's) agree within 2.1 %, and the LES's peak swing exceeds the quasi-static one
by up to 5 degrees, the span overshooting as it swings. On the spans the wind
runs nearly along (under 7 m/s across) the ratios scatter (ASCE 74 over the LES
0.71 to 1.07 on the peak over the mean drag, the LES's span factor 0.80 to 1.41)
and do not test :math:`G_w`. These numbers come from one run of one realization; no
test checks them.

Gusts from the RANS turbulence
------------------------------

A steady RANS wind gives the lines their mean load only. With any
``conductors.gust_type`` (``factor``, ``event`` or ``random``) every span
also gets a gust factor and a peak wind load from the mean wind and the
turbulent kinetic energy :math:`k` of the k-equation RANS closure where it hangs
(``Structures/Conductors/ERF_Gusts.H``). It needs
``erf.rans_type = kEqn`` on the conductors' anchor level and the flow's wind
(not ``prescribed_velocity``), and every span needs a horizontal extent; the run
stops otherwise, naming the input. With ``factor`` the lines' motion is not
changed: the factor is a diagnostic. ``event`` and ``random`` also add a gust to
the wind the lines and towers take (below).

The streamwise velocity fluctuation is :math:`\sigma_u = c \sqrt{k}`. In a
neutral surface layer :math:`\sigma_u` is about :math:`2.5\,u_*`, and the RANS
closure's equilibrium there is :math:`k = u_*^2 / C_{\mu 0}^2`, so the default
is :math:`c = 2.5\,C_{\mu 0}` (1.39 with the default ``erf.Cmu0`` of 0.5562);
``conductors.gust_sigma_factor`` sets :math:`c` instead. The RANS
:math:`k` (about :math:`3.2\,u_*^2` there) is well below a measured one (about
:math:`5.5\,u_*^2`), so :math:`c` must not be taken from a measured
:math:`\sigma_u/\sqrt{k}` (about 1.05). The lateral fluctuation is
:math:`\sigma_v = 0.8\,\sigma_u`. A span sees the wind's component normal to it:
with the mean horizontal wind :math:`U` at an angle :math:`\phi` to the span,
:math:`U_n = U \sin\phi`, and the fluctuation normal to the span is
:math:`\sigma_n = (\sigma_u^2 \sin^2\phi + \sigma_v^2 \cos^2\phi)^{1/2}`. With
the turbulence intensity :math:`I_n = \sigma_n / U_n` the span's gust response
factor takes ASCE 74's form with that intensity in place of the exposure's,

.. math::

   G = 1 + 2 g I_n \sqrt{B}, \qquad B = \frac{1}{1 + 0.8 L / L_s},

with :math:`L` the span's chord, :math:`g` the peak factor
(``gust_peak_factor``, default 2.7) and :math:`L_s` a length scale
(``gust_span_length_scale``, default ASCE 74's for ``asce74_exposure`` when it
is given, else exposure C's: 220 ft = 67.056 m for C, 170 ft = 51.816 m for B). ASCE 74 writes
:math:`G_w k_v^2 = 1 + 2.7 E \sqrt{B_w}` for wind normal to the span, its
:math:`E` being twice its exposure's turbulence intensity; for wind normal to
the span, :math:`g = 2.7` and the exposure's :math:`L_s`, the two are the same
when :math:`I_n = E/2` (unit test
``Gusts.WithTheExposuresIntensityItIsASCE74sFactor``, exposures B and C at
several heights and spans). The peak wind load per length is :math:`G` times
the mean one, :math:`(\rho/2) C_d d U_n^2`; the peak wind at a point is
:math:`U + g \sigma_u`. The form is linear in the fluctuation: on a span the
wind runs nearly along (:math:`U_n` below :math:`\sigma_n`) it misses the
load of the fluctuation itself and underestimates the peak load, and :math:`G`
grows as :math:`U_n` falls; the ``valid`` column of ``gusts.csv`` is 1 where
:math:`U_n \ge \sigma_n` and 0 on such spans. With :math:`U_n = 0` the
intensity, the factor and the loads are written as 0.

Every step that ends at or after ``stats_start`` takes, per span, the root-mean-square over
its nodes of the horizontal wind speed and of its component normal to the
span's chord in the horizontal plane (the mean load goes with the mean of the
square), and the mean of :math:`k = (\rho k)/\rho` over the nodes (a negative
:math:`k` from extrapolating below the first cell centre counts as 0), all
from the state at the step's start, where the nodes are then; the samples are
stamped with the step's start time, so the first one's time (``t_first``) can be
one step before ``stats_start``. Their running means per span are
``<output_root>_gusts_stats.csv`` of each line, and whenever the statistics
are written ``<diagnostics_dir>/gusts.csv`` gets one row per span from them:
its attachment height (the mean height above the ground of its two conductor
points, ``asce74.csv``'s height with ``asce74_wire_height = attachment``; the
nodes sampled sag below it), chord (between
those points, as in ``asce74.csv``), wind,
normal wind, :math:`k`, :math:`\sigma_u`, :math:`\sigma_n`, :math:`I_n`,
:math:`G`, the point gust, the mean load, the peak load and ``valid``. Nothing
is written before the first sample. The means are checkpointed with the other statistics;
a run restarted from a checkpoint written without gusts starts them at the
restart.

Travelling and random gusts
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``gust_type = event`` and ``random`` add a fluctuation :math:`u'` to the RANS
wind every line node and tower drag node takes, along the mean horizontal wind
of the span, insulator string or tower it acts on (none where that mean is
zero; the vertical wind is left alone); a negative gust larger than a low wind
reverses it. The gust acts on the lines and towers only: these types cannot be
used with ``drag_on_flow``, which would put the gust's drag into the RANS flow.
Every step :math:`\sigma_u = c \sqrt{k}` is taken from the mean :math:`k` over
each span's nodes, and over each tower's drag nodes weighted by each node's
share of the tower's fluctuating load, the face area it stands for (its drag
width times its length) times the horizontal wind on it (the drag's
fluctuation is :math:`2 q C_f A U u'`; the area alone in still air), so that
a tower's turbulence is its structure's and does not follow the number or
the layout of its drag nodes. An insulator string takes
the mean amplitude (``event``) or the mean gust (``random``) of the two spans
beside it; with independent random gusts on those spans its fluctuation is
the smaller, which matters little for a string a few metres long next to spans
of hundreds.

``event`` is one travelling 1 - cos gust. A plane front moving at the speed
:math:`c_g` (``gust_event_speed``) along the horizontal unit vector
:math:`\mathbf{e}` (``gust_event_direction``, degrees counterclockwise from
:math:`x`, default 0) crosses :math:`\mathbf{x}_0` (``gust_event_origin``,
default the domain's centre) at :math:`t_0` (``gust_event_time``), so it reaches
a point :math:`\mathbf{x}` at
:math:`t_a = t_0 + \mathbf{e} \cdot (\mathbf{x} - \mathbf{x}_0) / c_g` and lasts
:math:`D` there (``gust_event_duration``, default 10.5 s, the duration of IEC
61400-1's extreme operating gust, not its shape or amplitude):

.. math::

   u' = g \sigma_u \frac{1 - \cos 2 \pi s}{2}, \qquad s = \frac{t - t_a}{D},
   \qquad 0 \le s \le 1,

and 0 outside, with :math:`g` the peak factor (``gust_peak_factor``). It is
taken at the middle of each step, where the wind MoorDyn holds over the step
is valid; where the front has already passed a point when the run starts, the
gust acts there from the first step. The event is a point peak put on the whole
front at once: unlike the random gusts and the factor it is not reduced by
:math:`B`, so on a long span it is conservative. On a 400 m span crossed at
right angles with :math:`I = 0.15` its peak load is 1.97 times the mean, where
the factor gives 1.34.

``random`` gives every span and every tower an Ornstein-Uhlenbeck process
:math:`z` of unit variance and integral time scale :math:`T = L_u / U`,
:math:`U` the root-mean-square horizontal wind over its nodes (a tower's
weighted by its drag nodes' load shares), stepped exactly over each step
:math:`\Delta t`,

.. math::

   z' = z\, e^{-\Delta t / T} + \left(1 - e^{-2 \Delta t / T}\right)^{1/2} \xi,
   \qquad u' = \sigma_u \sqrt{B}\, z,

with :math:`\xi` standard normal. :math:`L_u` is ``gust_integral_length``, or
the integral length of IEC 61400-1's Kaimal spectrum, :math:`8.1 \Lambda_1`
with :math:`\Lambda_1 = 0.7 z` below 60 m and 42 m above, at the span's
attachment height or the tower's load-weighted mean drag-node height, where
its fluctuating load centres (about half the height of the test towers, whose
body tapers from 6 m to 1.5 m; ASCE 74 takes the exposure's turbulence at two
thirds of a structure's height, a proxy for the same centre on a prismatic
structure). :math:`B`, the part of
the point fluctuation's variance the whole structure feels, is the span's
:math:`B` above, and for a tower ASCE 74's
:math:`B_t = 1 / (1 + 0.375 h / L_s)`, :math:`h` its height (the cross-arm's
height plus the peak). A line that shares another's towers (``share_towers``) takes the
processes of that line's spans, and a line that names another in
``conductors.<line>.gust_with`` takes them too (a circuit on separate
lines), so the conductors of a circuit move together; otherwise each line's
gusts are independent of the others'. Each tower has its own process. :math:`\xi` comes from a hash of ``gust_seed``, the process
and the step count alone, so every rank draws the same numbers, another seed
gives another series, and a restart, which carries :math:`z`, continues the
same series. The first step draws :math:`z` from its stationary law, so the
series is stationary from the start. The model is the simplest with the right
variance and correlation time: it has no spectrum, no coherence along a span
beyond :math:`B`, and no lateral or vertical fluctuation. Its time scale is the
point's, while a span-averaged gust is more correlated in time (at low
frequencies the whole span feels the point fluctuation), so a span's gust has
too little energy at low frequencies and too much at high ones. A sampled
Ornstein-Uhlenbeck path has no derivative, so its peaks grow slowly as the step
shrinks: the mean largest value of a unit process over 600 s with
:math:`T = 9.5` s is 2.80, 2.91 and 3.01 at steps of 1, 0.3 and 0.03 s, so a
peak load moves by 1 to 2 % between RANS steps of 1 and 0.3 s. A tower's gust
and those of the spans beside it are independent: where the lines' pull is
about half the members' drag (8 against 18 kN on a hilltop tower at an
intensity of 0.15), the fluctuation of the tower's total load is 16 % below the
fully correlated one, and its peak about 6 % low.

Every ``diagnostics_int`` steps ``<diagnostics_dir>/gust_series.dat`` gets a
row, stamped with the step's start time, of the gust at each span's middle node
and at each tower's highest drag node (the event's at the step's middle, half
a step after the stamp). The unit tests ``Gusts.*`` check the
shape and arrival of the event, the normal numbers, the exact step, and the
variance and the autocorrelation :math:`e^{-\tau / T}` of a long series;
``Conductors.AnEventGustPassesEveryNodeAndTowerWithItsFront``,
``Conductors.RandomGustsFollowTheirProcessesAndACircuitSharesThem``,
``Conductors.ALineTakesAnotherLinesGustsAndAChangedSharingDrawsAfresh`` and
``Conductors.RandomGustsContinueAcrossARestart`` check the wind the lines and
towers take, and ``Conductors.TheGustFactorRunsOnLinesWithLatticeTowers`` the
factor on a section with lattice towers.

``Structures/Exec/PowerLines/les/rans`` runs the clamped network of the
ASCE 74 comparison in a k-equation RANS of the same hills: a terrain-fitted
mesh of the LES's domain and spacing (the k-equation's wall distance is
measured from the mesh's bottom, which an immersed terrain leaves flat), with
the LES precursor's mean profile as the inflow.
``Structures/Exec/PowerLines/les/compare_gusts.py`` compares every span with
the LES over 600 s of each. The LES's :math:`\sigma_u` at mid-span (its resolved
fluctuation and the precursor's subgrid energy at that height) over the RANS
:math:`\sqrt{k}` along the span gives :math:`c` = 0.76 to 1.59 across the eight
spans, median 1.09 (0.96 from the resolved part alone), where the default is 1.39.
With the default the random gusts' peak span load over the mean is 1.16 times
the LES's (median; 0.80 to 1.35) and the factor's :math:`G` 1.18 times (0.93 to
1.71); with :math:`c` = 1.10 (``gust_sigma_factor``) they are 1.07 (0.73 to
1.20) and 1.08 (0.84 to 1.48). The peak tensions agree with the LES's within
3 % with either; the peak swing is 0.78 to 2.0 times the LES's (median 1.39 with
the default; 2.0 on the 33 m span). The RANS mean load is 0.74 to 1.57 times the LES's, following the
mean wind across each span, so the peak loads themselves are 1.31 (default)
and 1.20 times the LES's (medians). On the two long spans of the line whose
chords lie at about 22 degrees to the wind (262 and 230 m) the random gusts'
peak over the mean falls below the LES's (0.80 and 0.89 with the default): they
have no lateral fluctuation, which carries much of the load normal to such a
span; on the same line's 33 m span, at the same angle, it is 1.25. The default is kept: these numbers come from one LES
realization of one network over one terrain, and no test checks them.

Component capacities and first failure
--------------------------------------

This section describes how the model compares the load on each structural
component of the lines and towers with that component's capacity, and how it
records the first time a capacity is exceeded. It is the first stage of the
line-failure work: it identifies which component is critical, when, and by
how much. By default the structure is not changed when a capacity is
exceeded. With ``conductors.failure_break`` the line lets go at the component,
and with ``conductors.tower_failure = release`` a tower over a footing's or
member's capacity lets go of its lines (`Release of a line at a capacity
exceedance`_).
The implementation is ``Structures/Conductors/ERF_Capacities.{H,cpp}``
(class ``CapacityMonitor``), called from ``Conductors``.

Demand-to-capacity ratio
~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`c` index the monitored components and :math:`t_n` the times at
which they are checked (defined under `Instants of the check`_). The
demand-to-capacity ratio of component :math:`c` at :math:`t_n` is

.. math::
   :label: cap-ratio

   r_c(t_n) = \frac{D_c(t_n)}{C_c(t_n)} ,

where :math:`D_c \ge 0` (N) is the demand on the component and
:math:`C_c > 0` (N) its capacity. The capacities are constant within a run. A
frame member's capacity is evaluated at each check from its steel temperature
(:doc:`TowerFrame`), which is fixed today; a coupling that heats the members
would make it vary. A component **exceeds** its capacity at
:math:`t_n` when :math:`D_c(t_n) > C_c(t_n)`. A demand equal to the capacity
is not an exceedance. The model records, for each component, the peak ratio
and the first exceedance time:

.. math::
   :label: cap-peak

   \hat r_c = \max_{n} \, r_c(t_n) , \qquad
   t_c^{*} = \min \{\, t_n : D_c(t_n) > C_c(t_n) \,\} .

The run's **first failure** is the component with the smallest
:math:`t_c^{*}`. When two components first exceed at the same instant, the one
listed first in ``capacities.csv`` is reported. The first exceedances of one
ERF step, those of every group of lines and towers, are logged in time order
once all groups have stepped, so the stop of ``conductors.failure_stop``
follows the earliest of them, whichever group stepped first.

Component demands and capacities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A component is monitored only when its capacity is given. The capacity of
each kind of component is a single input, so no input is added for footings
or members: the existing tower-type and frame inputs are reused. A run whose
tower types already give allowables, or whose frames give member design
data, is therefore monitored without further input; its flow, lines, towers
and existing logs are unchanged, and it writes ``capacities.csv`` and
``failures.dat`` besides them.

.. list-table::
   :header-rows: 1
   :widths: 16 44 40

   * - Kind
     - Demand :math:`D_c`
     - Capacity :math:`C_c` (input)
   * - conductor (one per span)
     - largest tension in the span, :math:`\max_s T_s`, over its nodes and
       its two ends
     - rated breaking strength, ``conductors.<line>.rated_breaking_strength``
   * - insulator (one per string)
     - tension at the top of the string, :math:`T_{\mathrm{top}}`
     - rated mechanical load, ``conductors.<line>.insulator_rated_load``
   * - hardware (one per dead end and per tower attachment)
     - magnitude of the line's pull on the support, :math:`|\mathbf{F}_a|`
     - rated load, ``conductors.<line>.hardware_rated_load``
   * - footing uplift (one per tower)
     - largest leg tension, Eq. :eq:`cap-legs`
     - ``conductors.<type>.allowable_uplift``
   * - footing compression (one per tower)
     - largest leg compression, Eq. :eq:`cap-legs`
     - ``conductors.<type>.allowable_compression``
   * - member (one per frame member with design data)
     - axial force on the governing side, Eq. :eq:`cap-member`
     - ASCE 10-15 design strength on that side

**Conductor.** :math:`T_s` (N) is the tension in segment :math:`s` of the
span as MoorDyn reports it. MoorDyn-C reports segment (mid-segment)
tensions, so the largest lies slightly below the tension at the attachment,
by about :math:`w\,\Delta s\,V / (2T)` with :math:`\Delta s` the segment
length and :math:`V` the vertical part of the end tension: 0.16 % for a
300 m span of 20 segments. The bundled stub includes the end node. The rated breaking strength (RBS) is the
conductor manufacturer's rating; for example, the 795 kcmil 26/7 ACSR
(aluminium conductor, steel reinforced) Drake conductor has an RBS of
140.1 kN.

**Insulator string.** :math:`T_{\mathrm{top}}` (N) is the tension in the
string's top segment, where it hangs from the cross-arm. The capacity is the
string's rated mechanical load: the specified mechanical load (SML) of a
composite insulator (IEC 61109), or the mechanical and electrical (M&E)
rating of its porcelain or glass cap-and-pin units (ANSI C29.2B). A string of
units in series has the rating of one unit, not the sum.

**Attachment hardware.** :math:`\mathbf{F}_a` (N) is the force of the line on
its support at attachment :math:`a`: the dead end ``end_a``, each tower, and
``end_b``. At a dead end it is the net force on the line's end node. At a
tower it is the net force on the top node of the insulator string, or on the
two span ends meeting there when the conductor is clamped. Both include the
node's share of the conductor's weight and drag with MoorDyn-C; the bundled
stub reports the end tension alone. On a tangent tower with the conductor
clamped, the force is the net of the two spans' pulls, mostly the weight
between the mid-spans. The capacity applies to the
whole load path: dead-end fittings, suspension clamps and shackles share one
rating per line.

**Footings.** The tower's foundation is four legs on a square of side
:math:`s` (the type's ``leg_spacing``, or ``base_width`` when it is not
given), with half-spacing :math:`a = s/2`. The
legs share the vertical load :math:`P` (N) and the overturning moment
:math:`(M_x, M_y)` (N m) about the base centre as a rigid cap would:

.. math::
   :label: cap-legs

   R_i = \frac{P}{4} + \frac{M_y\, x_i - M_x\, y_i}{4 a^2} , \quad i = 1,\dots,4 ,
   \qquad
   D_{\mathrm{up}} = \max\!\big(0, -\min_i R_i\big) , \quad
   D_{\mathrm{comp}} = \max\!\big(0, \max_i R_i\big) ,

where :math:`R_i` (N, compression positive) is the axial force in leg
:math:`i` at offset :math:`(x_i, y_i) = a(\pm \hat{\mathbf{e}}_{\parallel} \pm
\hat{\mathbf{e}}_{\perp})` from the base centre, with
:math:`\hat{\mathbf{e}}_{\parallel}` along the line and
:math:`\hat{\mathbf{e}}_{\perp}` along the cross-arm.
:math:`P = W - F_z` is the tower weight :math:`W` (the type's ``weight``)
plus the downward part of the members' drag and of the lines' pull, the
lines' pull being taken where each line hangs. For a moving tower the nodes'
inertial forces are included (see `Moving towers`_). For a frame tower, the
legs are the frame's support reactions instead (:doc:`TowerFrame`).

**Members.** For a frame member with tension :math:`N_t \ge 0`, compression
:math:`N_c \ge 0` and design strengths :math:`P_t` and :math:`P_c` in
tension and compression (ASCE 10-15; see :doc:`TowerFrame`), the governing
side is the one with the larger ratio:

.. math::
   :label: cap-member

   (D_c, C_c) =
   \begin{cases}
     (N_t, P_t) & \text{if } N_t / P_t \ge N_c / P_c , \\
     (N_c, P_c) & \text{otherwise,}
   \end{cases}

so :math:`r_c` equals the member's utilisation in
``tower_<name>_members.csv``.

Instants of the check
~~~~~~~~~~~~~~~~~~~~~

The check runs at the end of every coupling step between the lines and the
rest of the model:

* lines on fixed supports, and towers that stand still, are checked once at
  the end of each step of the anchor level, :math:`t_n = t + \Delta t`;
* lines on moving towers, with those towers, are checked at the end of each of
  the :math:`N` coupling steps of length :math:`h = \Delta t / N` that their
  group takes within the step (see `Moving towers`_), at
  :math:`t_n = t + k h`, :math:`k = 1, \dots, N`.

A tower that stands still is loaded, for the check, by the lines' pull at the
end of the step and by the members' drag held over the step. ``towers.dat``
reports the same tower at the step's start, so its values can differ from the
check by the change of the lines' pull over one step.

With ``conductors.capacity_check_interval`` the lines on fixed supports and
the standing towers are also checked inside the step (`Checks inside a
step`_), and with ``conductors.capacity_check_start`` no check is made before
that time (`When the checks start`_).

The lines step only on their anchor level, so the check does too. With
subcycling, a line on a refined anchor level is checked at every step of that
level. If the refinement ratio in time is :math:`r_t`, that is
:math:`r_t` checks per coarse step, stamped with the anchor level's times.
A call on another level records nothing. This follows the multi-level
requirement of the structure coupling (erf-model/ERF#5009, raised by
A. Almgren): the model acts on, and is reported at, its own level's steps.

Assumptions and limitations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Nominal values.** Demands are the computed loads and capacities are the
  given ratings, without load or strength factors. The RBS and the SML are
  failing loads, so an exceedance of them marks a failure, not a design
  limit. A design check against a factored limit, for example a fraction of
  the RBS, is made by entering that limit as the capacity.
* **Sampling.** Demands are sampled at the coupling instants only: once per
  ERF step for lines on fixed supports, once per coupling step for lines on
  moving towers. MoorDyn's internal steps between them are not sampled. For a
  300 m span the swing modes have periods near 6 s, so at
  :math:`\Delta t` = 0.1, 0.25 and 0.5 s a sinusoidal swing is sampled 60, 24
  and 12 times per period and its peak is underestimated by at most 0.14 %,
  0.9 % and 3.4 % of the fluctuating part, :math:`1 - \cos(\pi/N)`. The
  axial elastic mode (about 7 Hz) and impulsive loads, such as a release or
  ice shedding, are shorter than one step and are missed. A smaller ERF step
  samples more often, and so does ``conductors.capacity_check_interval``
  (`Checks inside a step`_).
* **Axial loads only.** Conductor, insulator and hardware demands are force
  magnitudes. Bending, torsion and fatigue of the hardware and the conductor
  at the clamp are not checked. Member checks are the axial checks of
  ASCE 10-15, without combined bending.
* **Rigid-cap footings.** Eq. :eq:`cap-legs` assumes a rigid cap on four
  identical legs. Foundation stiffness and soil-structure interaction enter
  only through the tower models, not the split.
* **Failure response.** By default an exceedance does not change the
  structure; ``conductors.failure_break`` and ``conductors.tower_failure``
  make the line let go (`Release of a line at a capacity exceedance`_). With
  ``conductors.failure_stop = true`` the run stops instead, with an abort after
  ``failures.dat`` and ``capacities.csv`` are written. The step's other logs,
  plotfile and checkpoint are not written; the last periodic checkpoint
  remains available. Exceedances restored from a checkpoint do not stop a
  restarted run; only a new one does.

Outputs and restart
~~~~~~~~~~~~~~~~~~~

``<diagnostics_dir>/capacities.csv`` holds one row per component: its kind
and name, the last demand and capacity, the peak ratio :math:`\hat r_c` with
its demand, capacity and time, and, once it has exceeded its capacity,
:math:`t_c^{*}` with the demand and capacity then. It is rewritten every
``diagnostics_int`` steps, at each first exceedance, at each checkpoint and on
a restart.
``<diagnostics_dir>/failures.dat`` appends one row,
``time kind component demand capacity ratio``, at each component's first
exceedance, and the run prints the same line, at most five a step (a frame
tower can have dozens of members over their strength at once). Component names are
``<line>`` or ``<line>_span<k>`` (conductor), ``<line>_string<j>``,
``<line>_end_a``, ``<line>_tower<j>``, ``<line>_end_b``,
``<tower>_footing_uplift``, ``<tower>_footing_compression`` and
``<tower>_m<id>``.

A checkpoint stores every component's record (17 significant digits) under
``<chk>/conductors/capacities``. A restart continues the peaks and keeps the
first exceedances, which are not reported again. ``failures.dat`` is written
again from the checkpoint's records, so it holds the first exceedances up to
the checkpoint, also in a new diagnostics directory. A restart without any
capacity removes ``capacities.csv`` of the run being restarted, as a run from
the start without capacities removes ``failures.dat`` and ``capacities.csv``
of an earlier run in the same directory. A component whose capacity is first given
on the restart starts there. A checkpoint record for a component no longer
monitored is dropped with a message. A capacity changed on the restart keeps
the component's record from before it: the earlier peak ratio and first
exceedance were found against the earlier capacity.

Verification of the checks
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Test
     - What is pinned
   * - ``Conductors.CapacitiesAreCheckedAfterEveryStepOfTheAnchorLevelOnly``
     - Two levels, lines on level 1. A level-0 call records nothing. The
       conductor demand equals the elastic catenary's still-air end tension
       to 0.5 % (stub −0.03 %, MoorDyn-C −0.28 %). An RBS of 0.99 times that tension is first exceeded at the
       end of the first level-1 step, and each first exceedance is logged
       once.
   * - ``Conductors.ACapacityAboveTheDemandIsNeverExceeded``
     - An RBS of 1.01 times the catenary tension: peak ratio
       :math:`1/1.01`, no exceedance.
   * - ``Conductors.AStandingTowersFootingsAndAttachmentsAreChecked``
     - Still air: leg compression :math:`(W - F_z)/4`, zero uplift; the
       tower hardware demand equals the line's pull there.
   * - ``Conductors.AStandingTowerIsCheckedUnderTheLinesPullAtTheStepsEnd``
     - A crosswind blows the section out in the first step: the tower's check
       after it equals the load the tower carries over the next step, not the
       still-air load it held during the first.
   * - ``Conductors.AMovingTowerIsCheckedAfterEveryCouplingStep``
     - Towers bending at 2 Hz, :math:`N = 10` coupling steps per ERF step: the
       first exceedance is at :math:`t = h = \Delta t / 10`.
   * - ``Conductors_Capacities_TwoLevel`` (CTest, stub)
     - ``Conductors_TwoLevel`` with the span on level 1 (0.125 s under level
       0's 0.25 s) and an RBS of 13 892 N: the first row of ``failures.dat`` is
       at 1.875 s, the first time the span's largest tension in ``S1.dat``
       passes the RBS, between two level-0 steps (a check at level 0's steps
       would give 2.0 s); the peak demand in ``capacities.csv`` equals the
       largest tension in ``S1.dat``.
   * - ``Conductors_Capacities_Restart`` (CTest, stub)
     - A checkpoint at step 5 (1.25 s) between the conductor's first
       exceedance (0.5 s) and the dead ends' (1.875 s): the restarted run's
       ``failures.dat`` equals the uninterrupted run's.
   * - ``Conductors_FailureStopAbort``, ``Conductors_FailureStopNeedsCapacityAbort`` (CTests)
     - The run stops at the first step naming the conductor; ``failure_stop``
       without a capacity stops at start-up.
   * - ``Conductors.FailureStopEndsTheRunAtTheFirstExceedance``,
       ``Conductors.FailureStopNeedsACapacity``
     - The abort names the component and time after ``failures.dat`` is
       written; ``failure_stop`` without a capacity is refused.
   * - ``Conductors.CapacityRecordsContinueAcrossARestart``
     - First exceedances and peaks survive a checkpoint; rows after it are
       dropped; a capacity added on the restart starts there.
   * - ``Capacities.*``
     - The monitor alone: strict exceedance, peak by ratio under a changing
       capacity, first failure on ties, bit-exact state round trip, the table
       and log formats, member governing side, refused inputs.

References
~~~~~~~~~~

* ASCE (2015). *Design of Latticed Steel Transmission Structures*,
  ASCE/SEI 10-15. American Society of Civil Engineers, Reston, VA.
* ASCE (2020). *Guidelines for Electrical Transmission Line Structural
  Loading*, ASCE Manuals and Reports on Engineering Practice No. 74, 4th ed.
  American Society of Civil Engineers, Reston, VA.
* ANSI C29.2B. *American National Standard for Insulators - Wet-Process
  Porcelain and Toughened Glass - Transmission Suspension Type.* American
  National Standards Institute, New York.
* IEC 61109:2008. *Insulators for overhead lines - Composite suspension and
  tension insulators for a.c. systems with a nominal voltage greater than
  1 000 V - Definitions, test methods and acceptance criteria.* International
  Electrotechnical Commission, Geneva.

Wind intensity beside the peak demand
-------------------------------------

This section describes a record a run can keep for a later fragility analysis:
beside each component's peak demand-to-capacity ratio (`Component capacities
and first failure`_), the wind at a fixed height near the component. One run
then gives, for every component, a pair of wind intensity and peak demand. It
is switched on with
``conductors.intensity_measure = wind`` and is off by default, in which case
nothing changes. The implementation is
``Structures/Conductors/ERF_Intensity.{H,cpp}`` (class ``IntensityRecorder``),
called from ``Conductors``.

Intensity points
~~~~~~~~~~~~~~~~

The wind is taken at **intensity points**, fixed for the run: one under the
middle of every span (the middle of the chord between its attachment points),
one at every attachment of every line (the two dead ends and each tower
attachment), and one at every tower's base. Each point sits at the height
:math:`h` (``conductors.intensity_height``, 10 m by default) above the ground
there, the terrain or the immersed surface under it. A span's conductor takes
the point under its middle, an insulator string and the hardware the point at
their attachment, and a tower's footings and members the point at its base.

Averaged wind
~~~~~~~~~~~~~

At every step of the lines' anchor level, from :math:`t_n` to
:math:`t_{n+1} = t_n + \Delta t_n`, the model samples the flow's horizontal
velocity :math:`(u, v)` at each point at :math:`t_n`. This is the wind the
step's loads are computed from. The intensity at point :math:`p` is the
horizontal speed averaged over the last :math:`T` seconds
(``conductors.intensity_averaging_time``, 3 s by default):

.. math::
   :label: im-average

   S_p(t_{n+1}) = \frac{\sum_{i \le n} w_i \, s_{p,i}}{\sum_{i \le n} w_i} , \qquad
   s_{p,i} = \sqrt{u_{p,i}^2 + v_{p,i}^2} , \qquad
   w_i = \left| \, [t_i, t_{i+1}] \cap [t_{n+1} - T, \, t_{n+1}] \, \right| ,

where :math:`s_{p,i}` is the sample of step :math:`i` and :math:`w_i` the
length of that step inside the averaging window. The oldest step that reaches
past the window's start counts with the part inside it. Before the run is
:math:`T` old the window holds what has been sampled. With :math:`T = 0`, or a
step longer than :math:`T`, the intensity is the last sample. The direction is
that of the averaged vector,

.. math::
   :label: im-direction

   \theta_p = \operatorname{atan2}\left( \sum_i w_i v_{p,i} , \ \sum_i w_i u_{p,i} \right) ,

in degrees from :math:`+x` towards :math:`+y`: the direction the wind blows
towards, in (−180°, 180°], and 0 when the averaged vector vanishes. With the
k-equation RANS (:ref:`RANS`) the turbulent kinetic energy :math:`k` at the
point is averaged with the same weights.

**What the average means.** :math:`S_p` is the :math:`T`-second average of
the wind the flow resolves, built from one sample per step. In a flow that
resolves gusts (a large-eddy simulation) it approaches the :math:`T`-second
gust at the point when the step is much shorter than :math:`T` (see
Limitations). In a RANS flow the velocity is already a mean, so :math:`S_p` is
the mean wind, and the gustiness is in :math:`k`. The record holds both and
applies no gust factor. The speed is averaged as a speed, so :math:`S_p` is
never below the length of the averaged vector; the two differ by a fraction of
a per cent in a gusty wind, and by 10 % when the wind turns through 90° inside
the window.

What is kept
~~~~~~~~~~~~

For each point the model keeps the run's largest intensity over the steps
whose window is full,

.. math::
   :label: im-max

   \hat S_p = \max_{n \,:\, t_n - t_0 \ge T} \, S_p(t_n) ,

with its time, direction and :math:`k`; :math:`t_0` is the time of the first
sample. An average over a window that is not yet full is a shorter average,
and is not compared: a run shorter than :math:`T` has no :math:`\hat S_p`. For each component :math:`c` with
point :math:`p(c)` it keeps the intensity of the step at which the peak ratio
:math:`\hat r_c` of Eq. :eq:`cap-peak` was reached, and the intensity of the
step of its first exceedance :math:`t_c^{*}`:

.. math::
   :label: im-kept

   S_c^{peak} = S_{p(c)}(t_{\hat n}) , \quad r_c(t_{\hat n}) = \hat r_c ; \qquad
   S_c^{*} = S_{p(c)}(t_c^{*}) .

**Why both the peak's wind and the run's largest.** A component's largest
demand need not occur at the largest wind: a span swinging back, a gust model,
a heated member or a neighbour's load can move it. :math:`S_c^{peak}` is the
wind that accompanied the peak demand; :math:`\hat S_p` is the hazard level
the run reached.

Outputs and restart
~~~~~~~~~~~~~~~~~~~

``<diagnostics_dir>/demand_table.csv`` holds one row per component, in
``capacities.csv``'s order: ``kind``, ``component``, the peak ratio with its
demand, capacity and time, the component's ``point``, the wind kept at the
peak (``peak_speed``, ``peak_direction``, ``peak_k``, and ``peak_wind_time``,
the end of the step it was sampled for), the point's largest
(``max_speed``, ``max_direction``, ``max_k``, ``max_time``), ``exceeded``, and
the first exceedance's time, demand and capacity with the wind kept then
(``first_speed``, ``first_direction``, ``first_k``, ``first_wind_time``).
``<diagnostics_dir>/intensity.csv`` holds one row per point: its name,
position (``x``, ``y`` and the elevation ``z`` it is sampled at), height above
the ground, the averaged wind now and the run's largest. Speeds are in m/s, directions in degrees, :math:`k` in
m\ :sup:`2`/s\ :sup:`2`; a field that does not exist is empty (:math:`k`
without the k-equation RANS, the first exceedance of a component that never
exceeded, the wind of a component never recorded, the largest of a run
shorter than :math:`T`). Both files are rewritten at the first step, every
``diagnostics_int`` steps, at each first exceedance, at each checkpoint and on
a restart, as ``capacities.csv`` is.
Point names are ``<line>`` or ``<line>_span<k>``, ``<line>_end_a``,
``<line>_tower<j>``, ``<line>_end_b`` and ``tower_<tower>``. A tower
attachment's point and its tower's point stand at nearly the same place and
read nearly the same wind. A line whose name equals one of these names stops
the run at start-up.

A checkpoint stores the samples inside the averaging window, every point's
largest and every component's kept winds (17 significant digits) under
``<chk>/conductors/intensity``, and a restart continues them: the restarted
run's tables equal the uninterrupted run's. The height and the averaging time
must be the checkpoint's, and the samples must still carry (or not carry)
:math:`k`, which changes with ``erf.rans_type`` on the anchor level and with
``conductors.prescribed_velocity``; otherwise the restart stops and names the
input. A checkpoint without the record (the measure switched on at the
restart) starts the averages and the maxima there, with a message: a
component whose restored peak is not passed again then has no wind beside it,
and the largest covers the time since the restart. A restart with the measure
off removes the two tables the earlier run left, with a message.

Limitations
~~~~~~~~~~~

* **The flow's wind only.** The gusts of ``conductors.gust_type`` are added
  to the wind the lines and towers take, not to the flow, so they are not in
  :math:`S_p`. With a RANS flow and a gust model the peak demand includes the
  gust and the intensity does not; :math:`k` at the point is recorded for a
  gust estimate made afterwards.
* **No exposure conversion.** :math:`S_p` is the model's wind at the point,
  over the terrain and roughness of the run. It is not converted to another
  exposure or height.
* **One sample per step.** :math:`S_p` is a sum over step-start samples, so it
  is a :math:`T`-second average only when the step is short against
  :math:`T`. For a synthetic wind (mean 15 m/s, standard deviation 3 m/s,
  integral time scale 5 s, :math:`T = 3` s) its root-mean-square difference
  from the true 3-second average is 2 % of that average's standard deviation
  at :math:`\Delta t = 0.125` s, 18 % at 0.9 s and 62 % at 3 s, where the
  largest value in 10 minutes is also biased high (by 0.4 m/s). A step of
  :math:`T/6` or less keeps the difference under about 17 %. The average also
  ends half a step early, since the last sample is taken at the step's start.
* **Below the first cell centre.** When :math:`h` is below the first cell
  centre, the wind and :math:`k` are extrapolated linearly from the first two
  cell centres, not taken from the wall law. For a logarithmic profile with a
  first centre at 12.5 m this reads 2.5 % high at 10 m.
* **Slow components.** The wind kept at a peak is that of the peak's step. A
  span swinging with a period of several seconds reaches its largest tension
  after the gust that caused it, so its :math:`S_c^{peak}` can be below that
  gust; for such components :math:`\hat S_p`, or an averaging time at least
  as long as the period, is the better intensity.
* **k from the RANS only.** A large-eddy simulation carries no :math:`k`
  column; its gustiness is in the resolved wind.
* **Start of the step.** The wind is sampled at the start of each step and
  stamped with its end, as the loads it produces are. On towers that move, the
  capacity checks inside one ERF step (`Instants of the check`_) all take the
  wind of that step.
* **A prescribed wind.** With ``conductors.prescribed_velocity`` the
  intensity is that velocity's horizontal speed at every point, without
  :math:`k`; ``intensity_height`` is then read and checked but changes no
  value.
* **Refined levels.** Every point must lie on the grids of the lines' anchor
  level with the cells around it; a refined patch that holds the lines but
  not the air at :math:`h` above the ground under them stops the run, naming
  the point.

Verification of the intensity record
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Test
     - What is pinned
   * - ``Intensity.TheAverageIsOverTheLastTSeconds``,
       ``Intensity.StepsOfDifferentLengthsCountWithTheirPartInTheWindow``,
       ``Intensity.AStepLongerThanTheWindowAndAZeroWindowKeepTheLastSample``
     - Eq. :eq:`im-average` by hand: equal steps before and after the window
       fills; steps of 2, 1 and 2 s with :math:`T = 3` s; the oldest step cut
       at the window's start; a step longer than :math:`T`; :math:`T = 0`.
   * - ``Intensity.TheSpeedIsAveragedAsASpeedAndTheDirectionIsTheMeanVectors``
     - A turning wind keeps its mean speed; Eq. :eq:`im-direction`'s
       convention at 0°, 45°, 180° and −90°, and the zero vector.
   * - ``Intensity.TheRunsMaximumKeepsItsTimeDirectionAndK``,
       ``Intensity.TheMaximumWaitsForAFullWindow``
     - Eq. :eq:`im-max` with its time, direction and :math:`k`, not the last
       value; averages over 1 s and 2 s of a 3 s window (20 and 11 m/s) are
       not the largest, which is 8 m/s at 3 s.
   * - ``Intensity.AComponentKeepsTheWindAtItsPeakNotTheLastWind``
     - Eq. :eq:`im-kept`: the wind of the peak's step while the wind goes on
       rising; a component without a point keeps none.
   * - ``Intensity.TheTablesCarryEveryField``,
       ``Intensity.BadInputsAreRefused``
     - Both tables' headers and rows, with empty fields; the two inputs'
       ranges and every refused call.
   * - ``Intensity.ARestartContinuesTheWindowAndTheMaxima``
     - A run split by a state round trip equals the straight run to the bit;
       a changed height, averaging time or :math:`k` is refused by name; a
       point that joins at a restart survives the next round trip.
   * - ``Conductors_Intensity_TwoLevel``, ``Conductors_Intensity_Height40``
       (CTests, stub)
     - ``Conductors_TwoLevel`` (crosswind :math:`v = 15 + 0.1 z` m/s, lines on
       level 1): every point reads :math:`15 + 0.1 h` m/s at 90° to 0.2 %, for
       :math:`h` = 10 m and for 40 m with :math:`T` = 1 s; the 2.5 s run has
       no largest with :math:`T` = 3 s and one stamped at 1 s or later with
       :math:`T` = 1 s; the table's peaks and first exceedances are
       ``capacities.csv``'s, and the conductor's peak demand is the largest
       tension in ``S1.dat``.
   * - ``Conductors_Intensity_PeakWind`` (CTest, stub)
     - ``Conductors_EventGust`` (two lines on lattice towers over hills,
       k-equation RANS): 12 of 15 components peak before the run's end; each
       kept wind is stamped with its component's peak time and differs from
       the point's last value; one point per span, attachment and tower, each
       attachment's at 10 m above the terrain of ``ground.dat`` (which
       differs by 53 m between attachments); every :math:`k` positive; no
       largest stamped before 3 s.
   * - ``Conductors_Intensity_Restart``,
       ``Conductors_Intensity_RestartWindowAbort`` (CTests, stub)
     - A checkpoint at 1.25 s with :math:`T` = 1 s (the window full and
       trimmed): the restarted run's two tables equal the uninterrupted
       run's; a restart with another averaging time stops, naming the input.
   * - ``Conductors_Intensity_*Abort`` (CTests)
     - An unknown measure, a height or averaging time out of range or given
       without the measure, the measure without a capacity, a point above the
       refined patch, and a line named like another's attachment each stop
       the run with the input's name.

An ensemble of runs for fragility curves
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A fragility curve needs many pairs of wind intensity and peak demand: runs at
several wind speeds (stripes), wind directions and, with random gusts, several
gust series. ``Structures/Tools/fragility_ensemble/fragility_ensemble.py``
(Python, standard library only) writes, runs and gathers such a set from one
base deck and a plan. It changes no model: every member is an ordinary run
with ``conductors.intensity_measure = wind``. The fitting of curves to the
gathered tables is a later step and is not part of the script.

**Members.** The plan (a JSON file) lists ``speeds`` (m/s), optionally
``directions`` (degrees, the convention of Eq. :eq:`im-direction`) and
``seeds``, and the members are all their combinations. Each member gets its
own directory with a deck that differs from the base deck only in the inputs
the manifest (``manifest.json``) lists for it. The base deck's ``FILE``
includes are written into the member's deck. The script reads one
``key = value`` per line and stops on any other line.

**The wind of a member.** Let :math:`(u, v, w)(z)` be the base deck's wind
profile and :math:`S_b` and :math:`\theta_b` its horizontal speed and
direction at the reference height :math:`z_r` (``reference_height``, 10 m by
default), found by linear interpolation between the profile's rows, as ERF
reads the file. Between a first row without wind (the ground of a
logarithmic profile) and the next row a straight line is not the profile, so
the script stops when :math:`z_r` lies there. :math:`z_r` is independent of
``conductors.intensity_height``. A member
with stripe :math:`S` and direction :math:`\theta` takes, at every height,

.. math::
   :label: ens-wind

   \begin{pmatrix} u' \\ v' \end{pmatrix} = \frac{S}{S_b}
   \begin{pmatrix} \cos\Delta & -\sin\Delta \\ \sin\Delta & \cos\Delta \end{pmatrix}
   \begin{pmatrix} u \\ v \end{pmatrix} , \qquad
   w' = \frac{S}{S_b} \, w , \qquad \Delta = \theta - \theta_b ,

so the profile keeps its shape and its turning with height, and reads
:math:`S` towards :math:`\theta` at :math:`z_r`. For a logarithmic profile
this is the same roughness length with the friction velocity multiplied by
:math:`S / S_b`. Without ``directions`` the base direction is kept
(:math:`\Delta = 0`). The stripe is a nominal value that orders the runs; the
intensity of a fragility analysis is the wind each run records at its
intensity points, which the terrain and the lines' surroundings change. The
largest wind of the scaled profile, at any height, must stay within 150 m/s;
it is written to the manifest as ``max_wind``.

The script applies Eq. :eq:`ens-wind` to the files the base deck takes its
wind from:

* with ``erf.init_type = input_sounding``, to the sounding
  (``erf.input_sounding_file``) and, when a face is an ``Inflow`` with a
  ``dirichlet_file``, to that profile, which is then the one :math:`S_b` and
  :math:`\theta_b` are read from;
* with ``conductors.prescribed_velocity``, to that velocity in the deck
  (there is no profile and no reference height; the flow is not changed).

The script stops, with the input's name, on a base deck that sets or drives
its wind in a way it does not scale: a constant ``<face>.velocity``,
``<face>.inflow_profile``, boundary planes, another ``erf.init_type``, a
restart (``amr.restart``), ``erf.abl_pressure_grad`` or ``erf.abl_geo_wind``.
Other inputs that hold a wind speed are left as they are, among them
``conductors.asce74_wind`` and the event gust's inputs.

**Boundaries.** In a direction that is not periodic, the face the member's
wind enters becomes the ``Inflow`` and the face it leaves the ``Outflow``,
each with the inputs the base deck gives its own inflow and outflow face
(the profile file, the inflow :math:`k`). A wind along a pair of faces
(:math:`|\cos|` or :math:`|\sin|` of :math:`\theta` below
:math:`10^{-12}`) makes both of them side faces, with the type of the base
deck's side faces or the plan's ``side_type``. An oblique wind therefore
enters through two faces and leaves through two; ERF keeps one inflow profile
per level, which both faces read. A lateral face of another type (``Open``,
for one) stops the script. An inflow turbulent kinetic
energy ``<face>.KE`` is multiplied by :math:`(S / S_b)^2`, which keeps the
inflow's turbulence intensity. ERF uses ``<face>.KE`` only together with
``<face>.density``; without it the inflow :math:`k` is zero-gradient and
follows the interior, and the script says so. Periodic directions need no
change, and a base deck without an inflow face cannot take a wind through its
walls: the script stops.

**Seeds.** A seed is ``conductors.gust_seed``, so ``seeds`` needs
``conductors.gust_type = random`` (`Gusts from the RANS turbulence`_). Members
that differ only in the seed share one flow and differ in the gusts the lines
and towers take.

**Other inputs.** ``inputs`` in the plan sets deck inputs for every member (a
capacity, a gust type, an averaging time), and a stripe may carry its own
(``{"speed": 40, "inputs": {"erf.fixed_dt": 0.4}}``). The time step, the
run's length and the averaging time are not changed by the script. A faster
stripe raises the advective Courant number with the wind aloft, not with the
wind at :math:`z_r`, and an oblique wind raises it again (by up to
:math:`\sqrt 2` on square cells), so the base deck's step must suit the
fastest member, or the fast stripes carry their own step. For each member the
manifest holds ``courant``, the largest :math:`(|u'|/\Delta x + |v'|/\Delta y)
\, \Delta t` of its profile on level 0 with ``erf.fixed_dt``, and the script
prints a note for a member above 1; a compressible deck's step is also bound
by the sound speed, which the number does not include. The inputs the
members' wind is derived from cannot be given in the plan: the wind and
restart inputs, ``geometry.is_periodic``, the lateral faces' inputs,
``conductors.intensity_measure``, ``conductors.diagnostics_dir`` and, with
``seeds``, ``conductors.gust_seed``. A file the deck names that stands beside
it (a terrain file, a member design file) must be listed in
``support_files``, which are copied to every member; otherwise the script
stops.

**Running.** ``fragility_ensemble.py run --plan plan.json --exe erf_exec``
runs the members one after another. Before each run it waits until the
machine's 1-minute load is below ``--load-start`` (6); during a run it stops
the run's processes when the load passes ``--load-pause`` (12) and continues
them below ``--load-resume`` (7). ``--ranks`` (4) may not exceed
``--max-ranks`` (4). A member that finished with return code 0 is skipped
when the command is given again; a member that failed or was interrupted runs
again from its start. ``--dry-run`` writes the directories and the manifest
and runs nothing. The manifest is written before the directories and holds a
SHA-256 of the base deck, the wind files and the support files: a directory
written from another plan or from files that have changed since is refused,
and so is a second driver on a directory in use (``run.lock``). The run's
processes are found as the launcher's descendants with ``ps``; when the
driver is interrupted or receives SIGTERM it ends the run with it, and that
member runs again next time.

**Shared spin-up.** The seed sets only ``conductors.gust_seed``, which the
flow does not see before the lines are added, so the members of one stripe
and direction share one flow up to the restart.
With the plan key ``spinup: {"max_step": n}`` the script writes one spin-up
run per distinct flow under ``spinup/``: the member's deck without any
``conductors.`` or ``erf.structures.`` input, run for n steps with a
checkpoint at step n (``erf.check_int``; ``erf.check_per`` and the
``amr.check_*`` keys, which ERF reads after the ``erf.`` ones, are removed). Members are grouped by a SHA-256 of that deck and of
their wind files. Each member restarts from the checkpoint
(``amr.restart``) with the structure model and runs to its own
``max_step``. The spin-up leaves the structure model out because a restart
keeps the model's peak demands, first exceedances, intensity record and
failure logs from the checkpoint; without them the lines are made afresh at
the restart time, so every member's records start on the spun-up flow. A
spin-up runs once, before its members; when it fails, or ends without its
checkpoint, its members are marked failed without a run.

**Gathering.** After the runs (or with ``collect``) the script writes
``demand_table_all.csv`` and ``intensity_all.csv``: the rows of every
complete member's ``demand_table.csv`` and ``intensity.csv`` behind the four
key columns ``member``, ``stripe``, ``nominal_direction`` and ``seed``. A
member is complete when its run finished and both tables can be read with
the header of the first complete member. ``members_status.csv`` has one row
per member with its status (``complete``, ``not_run``, ``failed``,
``no_demand_table``, ``no_intensity_table``, ``bad_demand_table``,
``bad_intensity_table``), return code and, for a failed run, the abort
message of its log. Members that are not complete are listed on the screen
and the command returns 1; a member goes into both tables or into neither.
The failure logs ``breaks.dat``, ``tower_failures.dat`` and ``failures.dat``
of the complete members are gathered the same way, and
``damage_table_all.csv`` holds one row per member: the largest peak ratio,
the number of components above 1, whether each log was written, and the
counts of failures, released lines, releases and lost towers.

**Limits.**

* The members' flows start from the scaled sounding, or with ``spinup``
  from a spin-up of it. A base deck that restarts from its own checkpoint
  (``amr.restart``) would start every member from that checkpoint's wind, so
  the script stops on it.
* An inflow through two faces has no regression test in ERF beyond the trial
  described in the tool's README; the corner column between an inflow and an
  outflow face takes the profile on one side and a zero gradient on the other.
  The profile is applied by mesh level index, so the terrain should be flat
  at every lateral face that can become an inflow; the script cannot check
  this.
* A direction within :math:`10^{-12}` of a face's plane counts as along it;
  any other direction, however close, makes that face an inflow with a small
  normal wind. Directions a few degrees off a face are best avoided.
* A wind along a periodic direction with side faces across it has no inflow;
  with surface drag nothing sustains it.
* The scaling keeps a neutral flow similar. Stratification is not scaled: a
  profile with an inversion has a bulk Richardson number that varies as
  :math:`1/S^2` between stripes.
* The run's length is the same for every stripe, so the number of
  flow-through times and of gust time scales (:math:`L_u / U`) in a run
  changes with the stripe.
* Members that differ in the seed only repeat the same flow solve and record
  the same intensity; their peak demands differ through the gusts. The data
  therefore suit an analysis stripe by stripe. The same seeds are used in
  every stripe, so the gust series of different stripes are not independent
  draws.
* The load gate reads the load of the machine the script runs on and stops
  processes started there; ranks on other nodes, or started by a launcher
  that does not keep them as its descendants, are not paused
  (``--no-load-gate``).
* At most 100 000 members.
* A run that stops by design at a first failure (``conductors.failure_stop``)
  returns an error code and is reported as failed; ensembles run with
  ``failure_stop`` off.
* Large-eddy members are not covered: their inflow needs a precursor or a
  perturbation method that has not been checked against one.

**Tests.** ``FragilityEnsemble_DryRun`` compares the member set, the decks,
the wind files and the manifest of a dry run with values worked out by hand
for three stored base decks (an inflow face with a profile, a periodic deck,
a prescribed wind), and requires a changed source file or a lost manifest to
be refused; ``FragilityEnsemble_Refused`` gives the plan's keys, the run
options and the base deck values out of range and requires a stop that names
each; ``FragilityEnsemble_Run`` uses a stand-in executable of several
processes to show that finished members are skipped, failed ones run again,
an earlier attempt's tables are not taken, the gate waits, stops every
process of the run and continues it, and a driver that is told to end leaves
no process behind;
``FragilityEnsemble_Collect`` gathers stored tables, among them a member
without a demand table, a failed one, one not run and one with another
header; ``FragilityEnsemble_Spinup`` compares the spin-up decks and the
members' restart paths with values written by hand, requires each spin-up to
run once, a failed or checkpoint-less one to hold back its members, a lost
checkpoint to be run again, an earlier attempt's checkpoint not to be taken,
the checkpoint name to follow ``erf.check_file`` and ``erf.file_name_digits``,
and the gathered failure logs and damage table to hold the stand-in's rows. No test runs ERF. The tool's README records two ensembles run by
hand: four members on ``Conductors_TwoLevel``, whose recorded wind equals the
planned stripe within 0.01 % and the planned direction within 0.001°, and an
oblique wind through two inflow faces on ``Conductors_RandomGusts``.

Fragility curves from an ensemble
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A fragility curve gives the probability that a component's demand exceeds its
capacity as a function of the wind speed. This subsection describes how
``Structures/Tools/fragility_fit/fragility_fit.py`` (Python, standard library
only) fits such curves to the table an ensemble gathers
(`An ensemble of runs for fragility curves`_). The fit is the multiple stripes
analysis of Baker (2015). It changes no model and reads no run directory, only
``demand_table_all.csv``.

**Data.** Member :math:`i` of the ensemble gives, for component :math:`c`,
the wind :math:`x_i` recorded at the component's intensity point (m/s) and
the peak ratio :math:`r_{i}` of demand to nominal capacity
(Eq. :eq:`cap-peak`). The wind is the point's largest averaged wind of the run
(``max_speed``, Eq. :eq:`im-max`) by default; ``intensity`` selects instead
the wind at the peak demand (``peak_speed``, Eq. :eq:`im-kept`) or the
member's nominal ``stripe``. A member that recorded no such wind is left out
of that component's curve and named, and so is a row without a peak ratio (a
component its run never recorded). The plan's keys and their ranges are in
the tool's ``README.md``.

**Outcome of a member.** With the nominal capacity, a member's outcome is the
model's own rule, exceedance when the demand is above the capacity:

.. math::
   :label: frag-outcome

   y_i = 1 \ \text{if } r_i > 1 , \qquad y_i = 0 \ \text{otherwise}.

A component kind may instead be given a capacity that is lognormal about its
nominal value :math:`C_n`, with mean :math:`m \, C_n` (``mean_factor``,
1 when left out) and coefficient of variation :math:`\delta` (``cov``). Then
:math:`\ln (C / C_n)` is normal with standard deviation :math:`\zeta` and
mean :math:`\ln m - \zeta^2 / 2`, and the outcome is the probability that the
capacity lies below the peak demand,

.. math::
   :label: frag-capacity

   y_i = \Phi\!\left( \frac{\ln r_i - \ln m + \zeta^2 / 2}{\zeta} \right) ,
   \qquad \zeta^2 = \ln \left( 1 + \delta^2 \right) ,

where :math:`\Phi` is the standard normal distribution function. The tool
has no default for :math:`m` and :math:`\delta`: they are the user's, with
their source. Two published examples: for the yield strength of tower steel,
Liu and Yan (2022, their Table 2) and Zhu et al. (2024, their Table 3) use a
lognormal distribution with a coefficient of variation of 0.07, the standard
(nominal) strength being its 0.05 quantile in the former; for the strength of
a conductor and of an insulator, Mazaheri and Alipour (2025) use lognormal
distributions with a coefficient of variation of 0.1 and do not state the
ratio of mean to nominal. None of the sources read for this tool gives
measured statistics for conductor, insulator, hardware or footing capacity
against a rated value.

**Fragility function and fit.** The curve is the lognormal distribution
function (Baker 2015, Eq. 1)

.. math::
   :label: frag-function

   P(\text{exceedance} \mid V = x) = \Phi\!\left( \frac{\ln (x / \theta)}{\beta} \right) ,

with median :math:`\theta` (m/s) and dispersion :math:`\beta`, the standard
deviation of the logarithm of the wind speed at which the component exceeds.
The estimates :math:`\hat\theta` and :math:`\hat\beta` maximise the
likelihood of the members' outcomes (Baker 2015, Eq. 11, with one observation
per member and its recorded wind),

.. math::
   :label: frag-likelihood

   \ln L(\theta, \beta) = \sum_i \left[ \, y_i \ln p_i + (1 - y_i) \ln (1 - p_i) \, \right] ,
   \qquad p_i = \Phi\!\left( \frac{\ln (x_i / \theta)}{\beta} \right) .

With the outcomes of Eq. :eq:`frag-outcome` this is the binomial likelihood
of the paper. With those of Eq. :eq:`frag-capacity` it is the limit of that
likelihood, per draw, when the capacity is drawn many times for every member.
Members with the same wind add up to the counts :math:`z_j` of :math:`n_j`
of the paper's stripes. The tool maximises Eq. :eq:`frag-likelihood` as the
probit regression :math:`\Phi(a + b \ln x)` the paper names as its
equivalent (:math:`\beta = 1/b`, :math:`\theta = e^{-a/b}`), by Fisher
scoring with step halving.

**Why maximum likelihood.** The observed fraction at a wind speed has a
variance that changes with the probability, so a least-squares fit of the
fractions weighs the stripes wrongly and gives biased parameters (Baker
2015). The likelihood also needs no repeated observations at one wind
speed, so every member enters with the wind it recorded, which the terrain
and the lines' surroundings make differ from the nominal stripe.

**Outcomes without a fit.** Equation :eq:`frag-likelihood` has no maximum
for some data, and for other data its maximum is of no use. The tool then
reports why in place of the two parameters:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - ``status``
     - Meaning
   * - ``fitted``
     - :math:`\hat\theta` and :math:`\hat\beta` are given.
   * - ``no_exceedance``
     - The outcomes add up to less than :math:`10^{-6}` exceedances.
       ``bracket_low`` is the largest wind of the data.
   * - ``all_exceed``
     - The outcomes add up to all but less than :math:`10^{-6}` members.
       ``bracket_high`` is the smallest wind of the data.
   * - ``separated``
     - No wind with a member that does not exceed lies above a wind with one
       that does. The likelihood then grows as :math:`\beta \to 0` and has
       no maximum. ``bracket_low`` is the largest wind at which no member
       exceeds and ``bracket_high`` the smallest at which all do; one of
       them is empty when the data have no such wind.
   * - ``one_intensity``
     - All members recorded one wind speed: one fraction cannot give two
       parameters.
   * - ``not_increasing``
     - The best curve does not rise with the wind (the same fraction at
       every wind, or a falling one).
   * - ``flat``
     - A best curve exists, but it rises by less than 0.05 between the
       smallest and the largest wind of the data: the data do not show the
       rise, and the median would be an extrapolation.
   * - ``not_converged``
     - The search did not reach a point where the likelihood's slope
       vanishes.
   * - ``no_intensity``
     - No member recorded the chosen wind.

The brackets say where the data change; they are not bounds on the median.
With two members at the largest wind and none exceeding, a median equal to
that wind would give this outcome one time in four.

``separated`` is the expected result for a component with a nominal capacity
in a set of RANS runs: the flow is the same in every member of a stripe, so
the demand rises with the stripe and crosses the capacity once. A dispersion
needs something that varies between members at one wind speed: gust seeds
(``conductors.gust_type = random``), a capacity distribution, or resolved
turbulence. With a capacity distribution the outcomes are never exactly 0 or
1, so ``separated`` arises only for the kinds left at their nominal capacity.

**Groups.** A group is a set of components (a line, a tower, the network),
named in the plan by a regular expression searched in the component names, a
list of kinds, or both. A member exceeds for the group when any of its
components does. With capacity distributions the capacities of the components
are taken as independent given the demands of the member's run, so

.. math::
   :label: frag-series

   y_i^{g} = 1 - \prod_{c \in g} \left( 1 - y_{i,c} \right) .

This is the form of the series rule Panteli et al. (2017, Eq. 4) use for
towers that fail independently of one another; here it is applied inside one
member, whose demands come from one run and are not taken as independent.
Every member must hold every component of the group, or the tool stops. A
group's components stand at different points, so the group's curve is fitted
against the member's nominal stripe; its median is therefore a nominal wind
speed and is not on the same scale as a component's, which is a recorded one.

**Spread of the estimates.** Baker (2015) assesses a way of collecting data
by simulating outcomes from an assumed fragility function, estimating the
function from them, and repeating this; the paper reports the bias and the
standard deviation of the estimates. The tool does the same with the fitted
function in place of an assumed one: it simulates an outcome for every member
from Eq. :eq:`frag-function` with :math:`\hat\theta` and :math:`\hat\beta`
at the member's wind, fits the simulated set, and repeats this ``resamples``
times (1000 by default). It reports the standard deviation of the two
estimates and, which the paper does not propose, their central interval of
probability ``confidence`` (0.9). Each curve has its own random sequence,
made from ``seed`` and the curve's name, so a curve's spread does not depend
on the other curves. Simulated sets without a fit are left out, the number
used is reported, and a note is printed when it is below 90 %: the spread
then describes only the sets that could be fitted.

The spread says how well this number of members at these wind speeds
determines the curve if the curve is right; it holds no model error. It is
given only when every outcome is 0 or 1. Outcomes that are probabilities
vary less between members than simulated 0 and 1 outcomes do, so the same
procedure would overstate their spread; those curves are reported without
one.

**Damage states.** With the plan key ``damage_states`` the tool also reads
the driver's ``damage_table_all.csv`` and fits one curve per listed state of
the line network. With :math:`r_i` the member's largest peak ratio,
:math:`L_i` its number of released lines and :math:`T_i` its number of lost
towers, member :math:`i` reaches

.. math::

   \mathrm{DS1}: r_i > 1, \quad
   \mathrm{DS2}: L_i \ge 1, \quad
   \mathrm{DS3}: T_i \ge 1, \quad
   \mathrm{DS4}: T_i \ge N.

**What the level is.** The states are cumulative:
:math:`P(\mathrm{DS}k \mid V)` is the probability of reaching level
:math:`k` or higher, which nests the counts. The member's level is the
highest of the four states it reaches, read from every record its run wrote
(a log not written is a state not reached), and its outcome is 1 for every
listed state up to that level and 0 above it, so a state's curve does not
depend on which other states are listed. A member that reaches a state
without a lower one is named in a note. The damage table must hold the same
members, stripes, directions and seeds as the demand table. A lost tower releases its spans, so DS3 implies DS2 in the runs.
The curves are fitted against the nominal stripe with the likelihood above.

**Why N has no default.** :math:`N`, the number of lost towers that makes a
cascade, is the plan's ``cascade_towers``, an integer in [2, 100 000]. No
source available to the tool fixes it, so it is the user's choice; DS4
without it stops the tool. The states are the runs' own outcomes at the
capacities of their decks, so ``capacity`` acts on the component and group
curves only.

**Through ERF.** ``Structures/Exec/CascadeStripes`` runs two stripes (10 and
25 m/s, two seeds each, a shared spin-up per stripe) with the line failures
on and a synthetic footing uplift allowable between the stripes' demands. The
CTest ``CascadeStripes_DamageStates`` requires every member to restart from
its stripe's checkpoint with no conductor state in it, no failure at 10 m/s,
at least 2 towers lost at 25 m/s, and the fit's four states 0 of 2 at
10 m/s and 2 of 2 at 25 m/s. The damage counts alone would not catch a
spin-up that kept the lines: its towers fail during the spin-up and the
members restart with the same counts.

**Outputs.** ``fragility_curves.csv`` has one row per component and per
group (and per direction with ``by_direction``): the wind the curve was
fitted against, the number of members, the sum of the outcomes, the status, :math:`\hat\theta`, :math:`\hat\beta`,
their intervals and standard deviations, the bracket, and the log-likelihood.
``fragility_points.csv`` has one row per curve and stripe for a plot: the
members, their mean wind, the sum and fraction of outcomes, and the fitted
probability at that mean wind.

**Limits.**

* Directions are pooled unless ``by_direction`` is set. A pooled curve is a
  lognormal fit to the mixture of the directions, each weighted by its number
  of members, and its dispersion then holds the difference between the
  directions as well as the scatter within one. The tool does not weigh directions by a
  wind rose; Zhu et al. (2024, Eq. 12) integrate a fragility per direction
  over the joint distribution of speed and direction.
* The recorded wind of a RANS run holds no gust, while the demand of a run
  with a gust model does (`Wind intensity beside the peak demand`_). The
  fitted median is then a mean wind speed.
* A lognormal capacity acts on the peak ratio of a run made with the nominal
  capacity. It is not right for a run in which a failure changes the loads on
  the other components (a released conductor).
* With a capacity distribution the fit estimates the fragility that includes
  the capacity's scatter. It is exact when that fragility is lognormal in the
  wind speed and the nearest lognormal function otherwise.
* Capacities of different components are independent in Eq.
  :eq:`frag-series`. Capacities that vary together (a common material lot, a
  common foundation soil) make a group less likely to exceed than the
  equation gives, so the group's curve is an upper bound on the probability.
* A lognormal function is assumed. The tool gives no test of that
  assumption beyond the points in ``fragility_points.csv``.
* No spread is given for curves fitted to probabilities. A spread from
  resampling the members would need its own source and is not in the tool.
* Damage states are fitted one by one and are not forced to keep their
  order. Two lognormal curves with different :math:`\beta` always cross at
  some wind, often outside the stripes; inside the stripes' range they can
  cross when the members are few.

**Tests.** ``FragilityFit_Published`` reproduces numbers printed in Baker
(2015): with a true :math:`\theta = 1` and :math:`\beta = 0.4`, 40
observations at each of three intensities 0.6, 1 and 1.5 give a standard
deviation of :math:`\hat\theta` of 0.056, and 20 give 0.078 (the paper's
p. 589; this tool gives 0.054 and 0.079 from 20 000 simulated sets, within
the sampling error of the paper's own simulation); 45 at each of 0.4, 0.8
and 1.2 give coefficients of variation of 0.06 and 0.20 (its Table 3; this
tool gives 0.057 and 0.19). Two columns of
that table are not reproduced. For 30 observations per stripe the paper
prints 0.08 and 0.26, and this tool and an independent estimator both give
0.070 and 0.23. For two stripes the paper prints 0.07 and 0.40; about one
simulated set in seven has no collapse at the lower stripe and therefore no
maximum, this tool leaves those sets out, and the paper does not say how it
treated them. ``FragilityFit_ByHand`` fits two stripes whose answer follows
from :math:`\Phi(0)` and :math:`\Phi(1)`, compares the optimum with a grid
search of the likelihood that shares no code with the optimiser (also for a
shallow curve that needs shortened steps), and checks the statuses of the
table above except ``not_converged``, which no test reaches, the tail of
:math:`\ln \Phi`, winds a hundred-thousandth apart, and Eq.
:eq:`frag-capacity`. ``FragilityFit_Table`` runs two stored tables through
the command (a separated conductor, a fitted hardware item, groups by Eq.
:eq:`frag-series`, directions, a spread that does not depend on the order of
the rows), ``FragilityFit_Refused`` gives every plan key a value out of
range and the tool tables it cannot use, and ``FragilityFit_DamageStates``
counts DS1 to DS4 by hand for twelve members, requires DS1 to equal the
curve of the group of every component, moves the cascade threshold, names a member that skips a state, requires a
state's counts not to change with the other states listed, and refuses a
damage table from other members. The tool's README records the fit of
one ensemble that was run.

**References.**

* Baker, J. W. (2015). Efficient analytical fragility function fitting using
  dynamic structural analysis. *Earthquake Spectra* 31(1), 579–599.
  doi:10.1193/021113EQS025M.
* Liu, C., and Yan, Z. (2022). Fragility analysis of wind-induced collapse
  of a transmission tower considering corrosion. *Buildings* 12, 1500.
  doi:10.3390/buildings12101500. Open access (CC BY 4.0).
* Mazaheri, P., and Alipour, A. (2025). Fragility functions for weather
  events to enhance resiliency of transmission lines: from experiments to
  statistics. *Electrical Transmission and Substation Structures 2025*,
  443–451. American Society of Civil Engineers.
* Panteli, M., Pickering, C., Wilkinson, S., Dawson, R., and Mancarella, P.
  (2017). Power system resilience to extreme weather: fragility modeling,
  probabilistic impact assessment, and adaptation measures. *IEEE
  Transactions on Power Systems* 32(5), 3747–3757.
  doi:10.1109/TPWRS.2016.2641463. Creative Commons Attribution 3.0.
* Zhu, C., Yang, Q., Huang, G., Zhang, X., and Wang, D. (2024). Fragility
  analysis and wind directionality-based failure probability evaluation of
  transmission tower under strong winds. *Journal of Wind Engineering and
  Industrial Aerodynamics* 246, 105668. doi:10.1016/j.jweia.2024.105668.

Ground contact
--------------

This section describes how the terrain holds up a conductor that sags or
falls to the ground. Without it a conductor has no ground: MoorDyn's bottom is
placed far below the domain, and a node that reaches the terrain leaves the
grids or the domain and stops the run. Ground contact is therefore needed
before a conductor can be released and followed to the ground (see
`Release of a line at a capacity exceedance`_), and it lets a run continue when a slack or heated span
sags onto the terrain. It is switched on with ``conductors.ground_contact``
and off by default, in which case every result is unchanged. The terrain
reaches MoorDyn through its 3-D seafloor; the record of which spans lie on the
ground is ``Structures/Conductors/ERF_GroundContact.{H,cpp}`` (class
``GroundContactLog``), called from ``Conductors``.

Contact law
~~~~~~~~~~~

MoorDyn-C pushes a node up when it lies below the seafloor. For node
:math:`i` of a line of diameter :math:`d` (m), at height :math:`z_i` (m) and
vertical velocity :math:`\dot z_i` (m/s), with segment length :math:`l_i` (m),
the contact force is vertical,

.. math::
   :label: gc-force

   B_i = \bigl[\, k_b \, (z_s - z_i) - c_b \, \dot z_i \,\bigr] \, d \, l_i
   \quad \text{when } z_i < z_s , \qquad B_i = 0 \ \text{otherwise},

where :math:`z_s(x_i, y_i)` (m) is the seafloor height under the node,
:math:`k_b` (Pa/m) the contact stiffness and :math:`c_b` (Pa s/m) the contact
damping. MoorDyn's contact acts on the node, which is the conductor's axis.
The seafloor is therefore the terrain :math:`h(x, y)` raised by the
conductor's radius,

.. math::
   :label: gc-seafloor

   z_s(x, y) = h(x, y) + \tfrac{1}{2} d ,

so that the conductor's surface, not its axis, rests on the terrain.

**Stiffness.** :math:`k_b` is MoorDyn's default, :math:`3 \times 10^6` Pa/m. A
conductor of mass :math:`\mu` (kg/m) at rest sinks below the seafloor by

.. math::
   :label: gc-rest

   \delta = \frac{\mu \, g}{k_b \, d} ,

0.19 mm for :math:`\mu` = 1.628 kg/m and :math:`d` = 28.1 mm, against a
radius of 14 mm. The run refuses a conductor with :math:`\delta > d/4`.

**Damping.** A node of mass :math:`\mu l_i` on the contact spring
:math:`k_b d l_i` oscillates at

.. math::
   :label: gc-freq

   \omega_g = \sqrt{\frac{k_b \, d}{\mu}} ,

independent of the segment length (227.6 rad/s for the conductor above). The
damping that makes this contact critical is

.. math::
   :label: gc-damping

   c_b = 2 \sqrt{\frac{k_b \, \mu}{d}} ,

26 367 Pa s/m for that conductor. The model writes this value for each line.
MoorDyn's default, :math:`3 \times 10^5` Pa s/m, is chosen for mooring chain.
On a conductor it is eleven times critical. With it, a 300 m span released
at one end 30 m above flat ground diverged as it reached the ground, while
with Eq. :eq:`gc-damping` it came to rest.

**Time step.** MoorDyn integrates Eq. :eq:`gc-force` explicitly with a
second-order Runge-Kutta scheme (explicit midpoint) at its internal step
:math:`\Delta t_M`. Its amplification factor on a linear oscillator is
:math:`G(z) = 1 + z + z^2/2`, with :math:`z = \lambda \Delta t_M` for each root
:math:`\lambda` of the damped node. The insulator strings share :math:`c_b`, so
their damping ratio differs from one. For line type :math:`j` with
:math:`a_j = c_b d_j / (2 \mu_j)`, the rate

.. math::
   :label: gc-rate

   R_j = 2 \max \Bigl( \omega_{g,j} , \; a_j + \sqrt{\max(a_j^2 - \omega_{g,j}^2, 0)} \Bigr)

is twice the largest root magnitude: :math:`2 \omega_{g,j}` for a critically
damped or under-damped type and close to :math:`2 c_b d_j / \mu_j` for a
strongly over-damped one. A critically damped linear contact is stable while
:math:`R \Delta t_M < 4`. The one-sided contact chatters earlier. Nine
conductors of 0.5 to 3 kg/m and 15 to 40 mm, each dropped 30 m onto flat
ground, landed cleanly in every case at :math:`R \Delta t_M` = 2. They sank 1 to
11 cm transiently at 2.5 and 3, some failed at 3.5, and all failed at 3.9. At
start-up each line must satisfy

.. math::
   :label: gc-step

   \max_j R_j \, \Delta t_M \le 1.5 ,

over its line types :math:`j` (the conductor and, with strings, the
insulator); otherwise the run stops, naming ``conductors.moordyn_dt`` and
``conductors.moordyn_cfl``. The conductor above, in 15 m segments at
``moordyn_cfl = 0.1``, gives 1.0; the default insulator string (60 kg over
2.5 m, 0.254 m), under-damped with a damping ratio of 0.78, gives less. A type
with a damping ratio below 0.1 would not be covered by the check; under the
conductor's :math:`c_b` it would need :math:`\mu_j / d_j > 5800` kg/m\ :sup:`2`.

The seafloor
~~~~~~~~~~~~

Each line's seafloor file, ``<diagnostics_dir>/<line>.seafloor.txt``, holds
the anchor level's terrain nodes under the line and one span length around it
(plus the string length when the line hangs from strings), clipped to the
domain. That is the farthest a span can swing after letting go of one end.
With ``conductors.failure_break`` the margin is the sum of the line's span
lengths, since several releases can leave more than one span hanging from one
support as a single conductor. MoorDyn interpolates bilinearly between
the nodes, as the model measures the terrain under a node, so Eq.
:eq:`gc-seafloor` holds node for node. The grid is the smallest range of
nodes that covers the line and that margin. A grid of more than
:math:`10^6` points is coarsened to every :math:`n`-th node. Nodes that the
anchor level's grids do not hold (a refined anchor level that does not cover
the whole area) take the lowest height sampled. The conductor cannot reach
them in any case, since a node off the anchor level's grids stops the run.
MoorDyn reads the file in single precision, which rounds heights near
:math:`10^4` m by up to 0.5 mm. It reads the grid's coordinates the same way:
at a UTM northing near :math:`4.5 \times 10^6` m their spacing in single
precision is 0.5 m. The run stops at start-up if two neighbouring grid
coordinates become equal; a domain placed near the origin avoids this. On a
periodic domain the grid is clipped at the domain's faces, not wrapped, so a
span that swings across a periodic face rests on the seafloor's edge values
there.

On the ground
~~~~~~~~~~~~~

A node is on the ground when its axis is at most one diameter above the
terrain, that is, its surface is within one radius of it. The number of span
:math:`k`'s nodes on the ground at time :math:`t_n` is

.. math::
   :label: gc-onground

   n_k(t_n) = \# \{\, i \in k : z_i(t_n) - h(x_i, y_i) \le d \,\} .

It is evaluated at the end of every anchor-level step. ``GroundContactLog``
keeps, for each span, the first time :math:`n_k > 0` and the largest
:math:`n_k` with the time it was first reached. A resting conductor lies at
:math:`z - h = d/2 - \delta`, well inside the criterion. A span's end nodes are
counted with it, so a node on the ground at a clamped tower counts in both
spans beside it.

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Input or constant
     - Value
     - Effect
   * - ``conductors.ground_contact``
     - ``false``
     - ``true``: the terrain is each line's seafloor; spans on the ground are
       logged
   * - :math:`k_b`
     - :math:`3 \times 10^6` Pa/m
     - contact stiffness, MoorDyn's default (``ConductorInputs::ground_stiffness``)
   * - :math:`c_b`
     - Eq. :eq:`gc-damping`
     - contact damping, critical for the conductor
   * - step limit
     - 1.5
     - right-hand side of Eq. :eq:`gc-step` (``ConductorInputs::max_ground_step``)
   * - grid limit
     - :math:`10^6` points
     - seafloor coarsening (``ConductorInputs::max_seafloor_points``)

Assumptions and limitations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Vertical contact.** MoorDyn's contact force is vertical, not normal to the
  terrain. On a slope of angle :math:`\theta` the vertical force has no
  along-slope component, so a node resting on a slope is in equilibrium and
  does not slide under gravity; its horizontal motion is not damped by the
  ground. The vertical offset :math:`d/2` of Eq. :eq:`gc-seafloor` lets the
  conductor's surface reach :math:`(d/2)(1 - \cos\theta)` into the slope,
  1.9 mm at 30° for :math:`d` = 28 mm.
* **No friction.** MoorDyn's seabed friction is a viscous approximation and is
  not used. A conductor on the ground moves along it only under its own
  tension and the flow's drag, resisted by nothing else. Wind loads near the ground are
  small, and the flow sampled there is ERF's.
* **The ground is the terrain.** Vegetation, buildings and other obstacles do
  not hold a conductor up. The terrain does not move during a run.
* **Start-up on the ground.** When a span already reaches the ground at the
  start, MoorDyn's stationary initial-shape solve stops early. The span then
  settles dynamically over the first seconds; in the tests it was at rest
  within 15 s.
* **Strings.** The insulator strings share the conductor's contact damping
  (one value per MoorDyn system) and rest on the same raised seafloor.
* **A node pressed below the terrain.** A conductor that lands at speed sinks
  below the terrain for an instant (Eq. :eq:`br-sink`). The flow is then
  sampled, and the drag spread, at :math:`d/2` above the terrain under that
  node.
* **Capacities.** Ground contact changes the loads the components carry and
  so their demands (see `Component capacities and first failure`_); it adds
  no capacity of its own.

Outputs and restart
~~~~~~~~~~~~~~~~~~~

``<diagnostics_dir>/ground_contact.csv`` holds one row per span:
``span,nodes,nodes_on_ground,most_on_ground,most_time,first_touch_time``. The
two times are empty for a span that has not touched the ground. It is
rewritten every ``diagnostics_int`` steps, at each checkpoint and on a
restart. The run prints a line at each span's first touch, at most five per
step. Span names are those of ``capacities.csv``: ``<line>`` or
``<line>_span<k>``.

A checkpoint stores each span's record (17 significant digits) under
``<chk>/conductors/ground_contact``. A restart rewrites the seafloor files
from the terrain and continues the records; first touches are not reported
again. A span the checkpoint does not hold starts at the restart, as do all
spans when the checkpoint was written without ground contact. A restart
with ``ground_contact = false`` keeps no records and removes the table; it
warns that a span resting on the ground will fall through it and stop the
run. A fresh run without ground contact removes the table and the seafloor
files an earlier run left.

Verification of ground contact
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Test
     - What is pinned
   * - ``Conductors.TheSeafloorHoldsTheTerrainNodesUnderTheLineAndOneSpanAround``
     - On the ramp :math:`h = 0.05x + 0.02y`: the seafloor file holds the
       smallest range of the anchor level's nodes covering the line and one
       span length around it (clipped to the domain), each at :math:`h + d/2`
       less ``surface_offset``. The
       MoorDyn input names the file, with :math:`k_b` and the :math:`c_b` of
       Eq. :eq:`gc-damping`.
   * - ``Conductors.ASlackSpanRestsOnTheTerrainWithGroundContactAndStopsTheRunWithout``
     - 330 m of conductor over a 300 m chord, ends 30 m above the ramp, in
       still air: every node on the ground rests at :math:`d/2` above the
       terrain where it lies (stub exactly, MoorDyn-C within 1 mm), and the count of Eq.
       :eq:`gc-onground` and its first time (the first step) are logged.
       Without contact the same run stops.
   * - ``Conductors.GroundContactRefusesAMoorDynStepTooLongForIt``
     - Five 60 m segments and ``moordyn_dt`` = 0.01 s violate Eq.
       :eq:`gc-step` (4.0 with MoorDyn-C, 4.6 with the stub): refused at
       start-up; so is a conductor of 200 kg/m at 28.1 mm, which would sink
       more than a quarter of its diameter.
   * - ``Conductors.GroundContactRecordsContinueAcrossARestart``
     - First touch, count and most survive a checkpoint and continue as in
       the uninterrupted run; a restart without contact keeps none; one with
       contact from a checkpoint written without it starts them there.
   * - ``Conductors.OnAMeshWithoutTerrainTheSpanRestsOnTheDomainsFloor``
     - Without terrain the seafloor is ``prob_lo[2]`` plus :math:`d/2`, and
       the span rests there, inside the domain.
   * - ``Conductors.TheContactRateIsTheFastestRootOfEachLineTypeOnTheGround``
     - Eq. :eq:`gc-rate` for a critically damped, an under-damped and an
       over-damped type.
   * - ``Conductors.AnImmersedTerrainPlacesEverythingOnItsSurface``
     - An immersed terrain's surface and the fitted mesh's nodes give the
       same seafloors.
   * - ``Conductors_GroundContact_Restart``, ``Conductors_GroundContact_Terrain_Restart`` (CTests, stub and MoorDyn-C)
     - Flat ground (``Conductors_PrescribedWind``, 330 m span) and the hills
       of ``Conductors_Terrain`` (both lines strung at 1.5 kN): the spans rest
       on the ground from the first step (observed: 13 of 21 nodes flat,
       four of five spans on the hills; these counts are not pinned), and the
       restarted run's ``ground_contact.csv`` and span logs equal the
       uninterrupted run's.
   * - ``Conductors_GroundContactOffAbort``, ``Conductors_GroundContactStepAbort`` (CTests)
     - The flat-ground span without contact stops with a node leaving the
       domain; the step check of Eq. :eq:`gc-step` in a full run.
   * - ``MoorDynInputWriter.TheSeafloorGridIsCoarsenedToItsLimitAndMissingHeightsTakeTheLowest``
     - The coarsening and the fill of nodes no grid holds.
   * - ``GroundContactLog.*``, ``MoorDynInputWriter.*Seafloor*``, ``MoorDynInputWriter.*GroundContact*``
     - The log alone (first touch once, the most and its time, refused
       counts and names, table, bit-exact round trip); the seafloor file
       format and the contact rows of the MoorDyn input.

References
~~~~~~~~~~

* Hall, M. and Goupee, A. (2015). Validation of a lumped-mass mooring line
  model with DeepCwind semisubmersible model test data. *Ocean Engineering*,
  104, 590-603.
* MoorDyn-C v2 documentation, input file options ``SeafloorFile``, ``kbot``
  and ``cbot``. https://moordyn.readthedocs.io

Release of a line at a capacity exceedance
------------------------------------------

This section describes what happens to a line when one of its components
exceeds its capacity and the run is to go on. The capacity check of
`Component capacities and first failure`_ only records an exceedance; with
``conductors.failure_stop`` it stops the run there. With
``conductors.failure_break = true`` the line instead lets go at that
component: a span's end leaves its support, or an insulator string leaves its
cross-arm, and the line is followed as it falls, swings and comes to rest on
the terrain. The option is off by default, in which case every result is
unchanged. It needs ``conductors.ground_contact = true`` (`Ground contact`_)
and at least one capacity on a line. The rule and the log are
``Structures/Conductors/ERF_LineBreaks.{H,cpp}`` (``decide_break``, class
``BreakLog``); the line makes the release in MoorDyn
(``ConductorLine::release_span_end``, ``release_string``).

Trigger and timing
~~~~~~~~~~~~~~~~~~

Let :math:`D_c(t_n)` (N) be the demand on component :math:`c` and
:math:`C_c` (N) its capacity, both as defined in `Component capacities and
first failure`_, checked at the end :math:`t_n` (s) of every step of the
lines' level. The component's first exceedance is at

.. math::
   :label: br-trigger

   t_c^{*} = \min \{\, t_n : D_c(t_n) > C_c \,\} .

A release is made once per component, at :math:`t_c^{*}`, after the step's
exceedances have been logged in time order. The line is whole during the step
that ends at :math:`t_c^{*}` and released from the next step on. In the rows
written at :math:`t_c^{*}` the tensions, positions and drag are therefore
those of the whole line (the span logs, ``failures.dat``, the tension in
``breaks.dat``), while the pull on the support that was let go is already
zero (``towers.dat``, ``transformers.dat``, the hardware's demand from the
next check on).

For lines on towers that bend, the check runs after every coupling step inside
the level's step (`Component capacities and first failure`_). With
``conductors.capacity_check_interval`` (`Checks inside a step`_) every line
is checked several times inside the step. An exceedance
found inside a step carries its own time in ``failures.dat``, and the release is made
at the end of the level's step that contains it. The tensions of Eq.
:eq:`br-end`, and the tension written to ``breaks.dat``, are those at the end
of that step.

What lets go
~~~~~~~~~~~~

Spans are numbered :math:`k = 1 \dots N` from ``end_a``, attachments
:math:`a = 0 \dots N` (0 at ``end_a``, :math:`N` at ``end_b``, the towers
between). :math:`T_k^a` and :math:`T_k^b` (N) are the tensions at span
:math:`k`'s ends towards ``end_a`` and ``end_b``.

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Component over its capacity
     - Release
   * - Conductor of span :math:`k`
     - Span :math:`k` leaves its support at the end with the higher tension,
       Eq. :eq:`br-end`
   * - Insulator string at tower :math:`j`
     - The string leaves the cross-arm; the conductor it carried goes with it
   * - Hardware at a dead end
     - The span there leaves the dead end
   * - Hardware at a tower with strings
     - The string leaves the cross-arm
   * - Hardware at a tower the spans are clamped to
     - Of the two spans meeting there, the one with the higher tension at the
       tower leaves it
   * - Tower footing or member
     - Nothing, unless ``conductors.tower_failure = release``: then the tower
       lets go of every line on it (`Tower failure`_)

The end of a span that lets go is

.. math::
   :label: br-end

   e_k = \begin{cases} b & \text{if } T_k^b > (1 + \epsilon) \, T_k^a , \\
   a & \text{otherwise,} \end{cases} \qquad \epsilon = 10^{-6} .

**Why a tolerance.** On a level span in still air the two end tensions are
equal up to the line solver's residual (one part in :math:`10^8` with
MoorDyn-C). Without :math:`\epsilon` that residual would choose the end. A
height difference of 1 m between the ends changes the tensions by the
conductor's weight over 1 m, about one part in :math:`10^3`, far above
:math:`\epsilon`.

An end already released is not a candidate: a span with one end free can only
let go at the other, and a component with nothing left to release does
nothing. A span released in a step is not released again in the same step by
another component's exceedance, since its tension has not fallen yet.

After the release
~~~~~~~~~~~~~~~~~

MoorDyn-C detaches the line's end from its point and puts it on a new free
point at the same place (``MoorDyn_BreakLine``). The end then moves under the
conductor's tension, weight, drag and the ground's contact force. Nothing
else in the line changes: the nodes, their numbering and the span's logs are
the same.

**Rest state.** A span released at one end comes to rest hanging from the
other. With the held support at height :math:`h` (m) above the ground under
it and the span longer than :math:`h`, the hanging length is :math:`h` less
the conductor's radius and the tension at the held end is its weight,

.. math::
   :label: br-rest

   T_{\mathrm{rest}} = w \left( h - \tfrac{1}{2} d \right) ,

where :math:`w` (N/m) is the conductor's weight per unit length less
buoyancy and :math:`d` (m) its diameter. The rest of the span lies on the
ground. Equation :eq:`br-rest` holds in still air. In a wind the ground,
which has no friction, lets the fallen length add its drag to the held end.

**What MoorDyn-C gives.** MoorDyn-C lumps the line at its nodes. At rest the
tension at the held end is the weight of the nodes that hang clear of the
ground,

.. math::
   :label: br-lumped

   T_{\mathrm{rest}}^{\mathrm{MD}} = w \, l \, n , \qquad
   w \left( h - \tfrac{1}{2} d \right) - w \, l \;\le\; T_{\mathrm{rest}}^{\mathrm{MD}} \;\le\; w \left( h - \tfrac{1}{2} d \right) ,

with :math:`l` (m) the segment length and :math:`n` the number of such nodes
below the support. It is below Eq. :eq:`br-rest` by up to one segment's
weight: 3 % for 1 m segments on a 30 m support, 50 % for 15 m segments. A
span that hangs free of the ground puts its whole weight :math:`w L` on its
support whatever the segment length. In still air a released span swings for
minutes before it rests, since only the air's drag damps it.

**Wind and drag on a fallen conductor.** The ground's contact is a spring and
a damper (Eq. :eq:`gc-force`), so a conductor that lands at speed sinks below
the seafloor of Eq. :eq:`gc-seafloor` for an instant. With the critical
damping of Eq. :eq:`gc-damping`, a node landing at speed :math:`v` (m/s)
sinks at most

.. math::
   :label: br-sink

   \delta_{\max} = \frac{v}{e \, \omega_g} ,

where :math:`\omega_g` is the contact's frequency of Eq. :eq:`gc-freq` and
:math:`e` is Euler's number: 39 mm at 24 m/s, the speed after a 30 m fall
without drag, and 34 mm at the 21 m/s the conductor above reaches in air
(:math:`\omega_g` = 227.6 rad/s; 34.0 to 34.5 mm in MoorDyn-C runs of a
conductor dropped flat). That is more than
the 14 mm the seafloor is raised by, so the node's axis passes below the
terrain. A node below the terrain is on no grid. For sampling the wind and the
air's state, and for spreading the drag, a node whose axis is below the
terrain is therefore taken at :math:`d/2` above it. The line's own positions
are not changed. This applies whenever ``ground_contact`` is on. A node at or
above the terrain, which includes every node at rest on it, is sampled where
it is; before, a node below the terrain stopped the run.

**Supports.** A support that was let go takes no load from that end. A tower
keeps the pull of what is still attached to it, which is no longer balanced
between its two sides; the footing and member checks see that load.

.. list-table::
   :header-rows: 1
   :widths: 34 16 50

   * - Input or constant
     - Default
     - Effect
   * - ``conductors.failure_break``
     - ``false``
     - ``true``: a line's component over its capacity lets go of the line at
       the end of that step; needs ``ground_contact = true`` and a line's
       capacity (or, with ``tower_failure = release``, a tower's); not with
       ``failure_stop``
   * - :math:`\epsilon`
     - :math:`10^{-6}`
     - tie tolerance of Eq. :eq:`br-end` (``break_tie_tolerance``)

Assumptions and limitations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Where a conductor breaks.** A conductor over its strength is released at a
  support, not cut in the span. A conductor parted in mid-span would leave two
  pieces, each hanging from its own support; here one support keeps the whole
  span.
* **One release per component.** A component is released at its first
  exceedance only. A later, larger load on a component already over its
  capacity releases nothing more. This includes a component whose release was
  not made because its span had let go at its other end in the same step, and
  a component already over its capacity in a checkpoint when a restart first
  sets ``failure_break``; the restart names each such component.
* **Towers.** A footing or member over its capacity is logged and releases
  nothing by default. With ``conductors.tower_failure = release`` the tower
  lets go of its lines (`Tower failure`_). Either way a release changes the
  loads on what remains, which may then exceed its own capacity.
* **Nominal capacities.** As for the capacity check, the capacities are the
  inputs' values with no strength factor and no scatter.
* **Timing.** The release is made at the end of the step that finds the
  exceedance, so the time of a release is resolved to one step of the lines'
  level.
* **Dynamic loads after a release.** The whip of a falling span is resolved by
  MoorDyn-C's internal step, but by default the demand is checked once per
  step of the lines' level: a peak shorter than that step may be missed (see
  the runs under `Verification of the release`_). Set
  ``conductors.capacity_check_interval`` to check inside the step (`Checks
  inside a step`_).
* **Sensitivity of a fall.** A falling conductor that lands and whips
  amplifies small differences. Two MoorDyn-C runs of the same case on 2 and
  on 3 MPI ranks, whose flows differ in the last digits, made the same
  releases at the same times; their span logs agreed to :math:`10^{-13}` until
  the span landed, to :math:`10^{-8}` for 2 s after, and differed by 10 %
  later. Read the motion after a landing as one realisation, and a release
  that follows from it (a cascade) as one that may differ between runs when
  its demand is close to the capacity. A restart on the same number of ranks
  reproduces the uninterrupted run.
* **Coarse segments.** The tension of a span that has come to rest on the
  ground is that of Eq. :eq:`br-lumped`, below Eq. :eq:`br-rest` by up to one
  segment's weight.
* **Refined levels.** The grids of the lines' level must cover the ground
  under the lines, with the reach of the sampler, since that is where a
  released conductor goes. A run stops, naming the node, when one leaves
  those grids or the domain.
* **No electrical model.** A conductor on the ground or swinging near another
  is only reported (``ground_contact.csv``, ``separation.dat``).
* **The stub.** The bundled MoorDyn stand-in has no dynamics: a released span
  is placed at once in the rest state of Eq. :eq:`br-rest`, so a stub run
  shows that state by construction and only the MoorDyn-C runs test it. A tie in Eq.
  :eq:`br-end` is the rule there, since its span carries the same tension at
  both ends. Release times agree between the stub and MoorDyn-C when the
  demand passes the capacity at the same step.

Checks inside a step
~~~~~~~~~~~~~~~~~~~~

By default the demands are compared with the capacities once per step
:math:`\Delta t` (s) of the lines' level, at its end. A load peak shorter than
the step can then pass unseen. With
``conductors.capacity_check_interval`` :math:`= \tau > 0` the lines advance
within the step in

.. math::
   :label: ci-pieces

   n = \left\lceil \frac{\Delta t}{\tau} \right\rceil , \qquad h = \frac{\Delta t}{n} \le \tau ,

equal pieces of length :math:`h`, and the components are checked after each,
at :math:`t + h, t + 2h, \dots, t + \Delta t`. The wind is the one sampled at
the start of the step, held over all pieces, as without the option. The
lines' own integration is MoorDyn's and keeps its internal step; only the
points at which ERF reads the loads are added. Lines on towers that bend take
at least :math:`n` coupling steps, each of which they already check.

The towers that stand still are checked with the lines, under the lines' pull
at each of these times. Releases and tower failures are still made at the end
of the step (Eqs. :eq:`br-trigger` and :eq:`tf-time`); the times in
``failures.dat`` and the peaks in ``capacities.csv`` are those of the pieces.

**Choosing the interval.** After a release the remaining support sees a peak
15 to 105 ms wide (`Verification of the release`_). An interval of 10 ms
places at least one check inside the narrowest of these; one of 125 ms can
miss more than half of its height. In MoorDyn-C runs of a released 300 m
span the checks saw, at the worst phase, 32 to 44 % of the peak at 0.25 s,
41 to 46 % at 0.125 s, 51 to 69 % at 0.05 s and 85 to 99 % at 0.01 s (the
peak was 12 ms wide on the finest line).

**Cost and effect on the answer.** Each piece of a line on fixed supports is
``conductors.substeps`` calls of MoorDyn and one check; a line on towers that
bend takes the larger of its coupling steps and :math:`n`. Advancing a line in
more calls changes its answer, since MoorDyn divides each call into its own
steps and ends it on a shorter one: by about :math:`10^{-7}` in the tension
of a line at rest; for a falling conductor by :math:`10^{-6}` in position and
up to a few per cent in the tension read at a check. The lines on fixed supports are advanced in pieces whenever any
capacity is checked, also one that is only on a tower that bends. With
:math:`\tau = 0`, the default, nothing changes.

**What stays per step.** The wind, the gusts, the conductors' temperatures,
the statistics, ``ground_contact.csv`` and the wind of
``conductors.intensity_measure`` are sampled once per step, as without the
option. A peak found inside a step keeps the averaged wind of that step's
end beside it in ``demand_table.csv``.

.. list-table::
   :header-rows: 1
   :widths: 38 14 48

   * - Input or constant
     - Default
     - Effect
   * - ``conductors.capacity_check_interval``
     - 0
     - :math:`\tau` (s), in [0, 3600]; 0: one check per step of the lines'
       level; needs a capacity to check
   * - most checks per step
     - 100 000
     - a smaller :math:`\tau` stops the run at its first step, naming the
       input (``ConductorInputs::max_capacity_checks``)

When the checks start
~~~~~~~~~~~~~~~~~~~~~

By default the components are checked from the first step. A run often
starts away from the state it is meant to study: the flow spins up, or lines
switched on at a restart start from their still-air shape and their towers,
if they bend, ring from rest. The loads of that start-up are then recorded as
peaks, and with ``failure_break`` or ``tower_failure`` they can release a
line. With ``conductors.capacity_check_start`` :math:`= t_s` a check stamped
with a time :math:`t < t_s` is not made:

.. math::
   :label: cs-start

   \text{checked at } t \iff t \ge t_s \, (1 - \delta) .

The slack :math:`\delta = 10^{-9}` (four units of the last place in single
precision; taken of 1 s when :math:`t_s < 1` s) covers the roundoff of a time
that is a sum of steps: ten steps of 0.1 s end at 0.9999999999999999 s, and
the check stamped there is the one meant by :math:`t_s = 1` s.

Before :math:`t_s` nothing is recorded: no demand, no peak, no exceedance,
no release and no tower failure, and the components' demand and peak fields
in ``capacities.csv`` are empty. The lines and towers themselves advance as
without the option. The first check is the first one stamped at or after
:math:`t_s`: the end of a step, or of a piece or coupling step inside it
(`Checks inside a step`_). The time is the run's, so a restart keeps it.
With ``conductors.intensity_measure = wind`` the largest averaged wind of
each point (Eq. :eq:`im-max`) is likewise looked for from :math:`t_s` on; the
averaging window is filled from the first step as before.

.. list-table::
   :header-rows: 1
   :widths: 38 14 48

   * - Input
     - Default
     - Effect
   * - ``conductors.capacity_check_start``
     - 0
     - :math:`t_s` (s), finite and :math:`\ge 0`; needs a capacity to check

Tower failure
~~~~~~~~~~~~~

With ``conductors.tower_failure = release`` a tower fails when one of its
checked components first exceeds its capacity: the footings' largest uplift
or compression against the tower type's allowable values, or a frame tower's
member against its design strength (`Component capacities and first
failure`_). With :math:`c \in t` the checked components of tower :math:`t`,
its failure time is

.. math::
   :label: tf-time

   t_t^{*} = \min_{c \in t} \; t_c^{*} ,

with :math:`t_c^{*}` from Eq. :eq:`br-trigger`. A tower fails once. At the end
of the step that finds :math:`t_t^{*}` every line still attached to it is let
go:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Line at the tower
     - Release
   * - Hangs from an insulator string
     - The string leaves the cross-arm
   * - Clamped to the tower
     - The span before the tower leaves it, then the span after it

Lines of other circuits that share the tower are released with it. Ends
released before are left as they are. A failed tower lets go of every end
still on it, also of a span that let go at its other end in the same step; a
span the tower lets go of is not released again in that step by its own
conductor or hardware. The tower itself is not removed: it
stays where it stands, keeps its wind load and carries no line afterwards.
Its footing and member checks go on and see no further load from the lines.

**Why this is a cascade mechanism.** A tower that lets go of a conductor
passes the conductor's weight to its neighbours. On a line with strings, the
next tower then carries its share of the longer span beside it. An example,
by the elastic catenary and by MoorDyn-C: spans of 300, 500 and 100 m with
attachments 60 m up, strung to pull equally along the line. Towers 1 and 2
carry the conductor of 397 m and 282 m. After tower 1 fails, the conductor
hangs 800 m from the dead end to tower 2, which carries that of 423 m at
rest: 1.44 times its vertical load before. Before it rests the load is
larger: in the MoorDyn-C run it reaches 21.5 kN against 7.3 kN at rest, as
the released conductor falls and swings. Whether the next tower fails depends
on its own capacity and on that dynamic load (`Checks inside a step`_).

The transfer is smaller where the released conductor reaches the ground,
which then carries part of it. With spans of 300, 400 and 200 m (unstretched
lengths 301.5, 402 and 201 m) on 30 m towers, the dead ends 30 m up, the 700 m span lies on the ground at rest and tower 2 carries 1.11
times its load before, against 1.47 times clear of the ground.

.. list-table::
   :header-rows: 1
   :widths: 34 16 50

   * - Input
     - Default
     - Effect
   * - ``conductors.tower_failure``
     - ``none``
     - ``release``: a tower whose footing or member exceeds its capacity lets
       go of every line on it at the end of that step; needs
       ``failure_break = true`` and a capacity on a tower

**Limits of the tower model.** The failed tower does not collapse, lose
members or shed load to other members: only its lines are released. A real
tower that fails falls with its conductors and pulls on its neighbours as it
goes; here the neighbours see only the released conductor. Which towers fail
together is decided by the step: all towers whose first exceedance falls in
one step of the lines' level fail at its end.

``<diagnostics_dir>/tower_failures.dat`` holds one row per failed tower:
``time tower component releases``, with the component whose exceedance failed
it and the number of releases made. Each release is also a row of
``breaks.dat``, naming that component. A checkpoint stores the failed towers
under ``<chk>/conductors/tower_failures``; a restart does not fail them
again. A component already over its capacity in a checkpoint does not fail
its tower when a restart first sets the option; the restart names it.

State at rest after a release
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After a span has let go, the rest of its section comes to a new state at
rest. On a line with insulator strings the string beside the lost span is no
longer pulled from both sides: it leans into the span that stays, that span
slackens, and its tension falls to a residual value. The tower then carries
that residual tension along the line, a load it does not see while the line is
whole. This state follows from the balance of forces alone, and
``Structures/Tools/broken_wire/broken_wire.py`` computes it without ERF or
MoorDyn, as the reference the model is compared with (`Verification of the
release`_). The section lies in one vertical plane, :math:`x` along the line
and :math:`z` up, in still air.

**One cable.** Every span and every string is a heavy elastic cable of
unstretched length :math:`L` (m), weight :math:`w` (N per unstretched metre,
net of buoyancy) and axial stiffness :math:`EA` (N). Let :math:`s` (m) be the
unstretched arc length from the cable's end :math:`a`, and
:math:`(H, V_a)` (N) the tension there as a vector pointing into the cable.
The piece between :math:`a` and :math:`s` is held by the tensions at its two
ends and by its weight :math:`w s`, so the tension vector and its size at
:math:`s` are

.. math::
   :label: bw-tension

   \mathbf{T}(s) = \left( H ,\; V_a + w s \right) , \qquad
   T(s) = \sqrt{H^2 + (V_a + w s)^2} .

The cable lies along its tension and each of its elements is stretched by the
factor :math:`1 + T / EA`:

.. math::
   :label: bw-shape

   \frac{dx}{ds} = H \left( \frac{1}{T} + \frac{1}{EA} \right) , \qquad
   \frac{dz}{ds} = (V_a + w s) \left( \frac{1}{T} + \frac{1}{EA} \right) .

Integrated over the cable, with :math:`V_b = V_a + w L`, the end :math:`b`
lies at

.. math::
   :label: bw-reach

   \Delta x = \frac{H L}{EA} + \frac{H}{w} \left[ \sinh^{-1} \frac{V_b}{|H|} - \sinh^{-1} \frac{V_a}{|H|} \right] , \qquad
   \Delta z = \frac{V_a L + \tfrac{1}{2} w L^2}{EA} + \frac{T(L) - T(0)}{w}

from the end :math:`a`. For a level span, :math:`V_a = -wL/2`, these are the
elastic catenary of `Verification`_. Given the two ends, Eq. :eq:`bw-reach`
is solved for :math:`(H, V_a)` by Newton's method. The cable pulls on what
holds its end :math:`a` with :math:`(H, V_a)` and on what holds its end
:math:`b` with :math:`-(H, V_b)`.

**The strings' feet.** The dead ends and the cross-arms are fixed. The lower
end of each string, where the spans either side meet, is free, and at rest
the pulls :math:`\mathbf{F}_c` (N) of the cables that end there add up to zero:

.. math::
   :label: bw-balance

   \sum_{c \,\in\, \text{foot } j} \mathbf{F}_c \;+\; \mathbf{P}_j = 0 ,
   \qquad j = 1 \dots N - 1 ,

where :math:`\mathbf{P}_j` (N) is a load hung on the foot (below). The two
coordinates of every foot are the unknowns. They are found by Newton's method
from the whole section's state, the pull of the element that lets go being
taken away in steps.

**What a release changes.** A span released at one end no longer joins its two
attachments. It hangs from the one it is still on and loads it with the
weight of its hanging length, straight down, since the ground has no
friction. Let :math:`h_r` (m) be the height of that attachment above the
level the conductor's axis rests at on the ground, half a diameter above the
ground. The hanging unstretched length :math:`l_h` (m) is the whole span if it
hangs clear, and otherwise stretches under its own weight to reach the
ground:

.. math::
   :label: bw-hang

   \mathbf{P} = (0, \, -w \, l_h) , \qquad
   l_h \left( 1 + \frac{w \, l_h}{2 \, EA} \right) = h_r \quad \text{if } L \left( 1 + \frac{w L}{2 EA} \right) > h_r , \qquad l_h = L \text{ otherwise.}

Where the attachment is a fixed point :math:`h` above the ground,
:math:`h_r = h - d/2`, and this is Eq. :eq:`br-rest` less the stretch
(0.24 mm of 30 m). A string
released from its cross-arm (`Tower failure`_) leaves its foot without a
support: the foot carries the string's weight, :math:`\mathbf{P} = (0, -W_i)`,
and the two spans meeting there hang as one.

**The lumped chain.** MoorDyn integrates a chain of :math:`n` straight links of
unstretched length :math:`l = L / n`, with the weight :math:`w l` at every
inner node and half of it at each end. Link :math:`i` then carries

.. math::
   :label: bw-lumped

   \mathbf{T}_i = \left( H ,\; V_a + w \, l \, (i - 1) \right) , \qquad i = 1 \dots n ,

lies along it with the length :math:`l \, (1 + T_i / EA)`, and the chain pulls
on its supports with :math:`(H, V_a - w l / 2)` and
:math:`-(H, V_a + w L - w l / 2)`. With ``segments`` the tool solves this
chain in place of Eqs. :eq:`bw-tension` to :eq:`bw-reach`, the strings in two
links as ERF writes them, so that MoorDyn's state at rest can be compared
link for link and the effect of the links' length read from the difference
between the two. A released span in links loads its support with the nodes
that hang clear of the ground, :math:`w \, l` each, and with the held end's own
node, :math:`w \, l / 2`. That load changes by a whole node's weight as a node
reaches the ground, so for a foot height at which a node just touches, the
chain as the tool treats it has no state at rest and the tool reports none;
MoorDyn's ground is a spring and carries such a node in part.

**Assumptions.** Still air; one vertical plane; fixed dead ends and
cross-arms (towers that do not bend); the strings' feet without mass; a
ground without friction, so that a fallen length pulls on nothing; a string
that stays in tension; every conductor still held at both ends, or left
without its string, clear of the ground (the tool stops, saying so, where one
would lie on it). The tool reports the state at rest only: the loads
while the released span falls are larger (`Verification of the release`_).

**Comparison with MoorDyn-C.** The unit tests ``BrokenWire.*`` release a span
or a string of a conductor line in still air and compare the towers' and the
dead end's pulls with the tool's for the same chain of 20 links per span. In
still air the section keeps swaying about its state at rest for minutes,
because the air's drag falls with the square of the speed, so the pulls
compared are means over 60 s. They agree within 1 N in 13 500 N for the whole
section and within 0.3 N after a release; the tests require
:math:`3 \times 10^{-4}` of the largest pull. ``Structures/Exec/BrokenWire``
makes the same comparison through ERF, the release made by a capacity
(``conductors.failure_break``):

.. list-table::
   :header-rows: 1
   :widths: 34 16 16 17 17

   * - Section, span 1 released
     - Span 2 keeps, model
     - Force balance
     - Tower 1 along the line, model (N)
     - Force balance (N)
   * - 3 spans of 300 m, 2.5 m strings, at tower 1
     - 70.7 %
     - 70.6 %
     - 9344.9
     - 9342.8
   * - 3 spans of 300 m, 4 m strings, at tower 1
     - 62.5 %
     - 62.4 %
     - 8260.0
     - 8258.0
   * - 300, 500 and 100 m, 2.5 m strings, at tower 1
     - 87.0 %
     - 87.0 %
     - 11748.3
     - 11746.9
   * - 3 spans of 300 m, 2.5 m strings, at its dead end
     - 71.2 %
     - 71.2 %
     - 9413.5
     - 9411.6

Every tower pull of these four runs, along the line and vertical, is within
0.02 % of tower 1's largest pull (means over the last 240 s of 600 s; MoorDyn-C
2.7.1). The longer the strings and the more equal the spans, the more tension
is lost.

**Loads while the span falls.** The same case records the largest loads after
the release, sampled every 0.01 s (`Checks inside a step`_). With a middle
span of 0.5 to 2 times the outer spans of 300 m they were 3.3 to 4.4 times the
loads at rest after the release in the strings and at the dead end, and 2.4
to 4.0 times in the towers' legs. On eight spans of 300 m, with every tower
letting go of its lines when a leg passes its bearing capacity
(``conductors.tower_failure = release``), all seven towers failed for
capacities up to 5 times the legs' load before the release, two or three at 8
times, none or one at 12 times and none from 16 times on; where two numbers
are given the larger is with checks every 0.01 s. These are single runs of
rigid towers that stay standing, with no structural damping beyond MoorDyn's:
a fall amplifies roundoff ("Sensitivity of a fall" above), and a tower that
yields would take up part of the load. They are measured, not compared with
a design guide's load factors, to which we had no access.

Outputs and restart
~~~~~~~~~~~~~~~~~~~

``<diagnostics_dir>/breaks.dat`` holds one row per release:
``time line released component tension``. ``released`` is
``span<k>_end_a``, ``span<k>_end_b`` or ``string<j>``; ``component`` is the
name in ``capacities.csv`` of the component whose exceedance released it;
``tension`` (N) is the tension at the released end, or at the string's top,
when it let go. The run prints a line per release.

A checkpoint stores the releases in order (17 significant digits) under
``<chk>/conductors/breaks``. MoorDyn-C's saved state holds one free point per
release and fits only a system with the same releases. On a restart each line
therefore makes its releases again, in order, before it loads its state. This
is done whether or not the restart sets ``failure_break``; without it no
further release is made. A restart from a checkpoint that holds releases
needs ``ground_contact = true`` and is refused without it. On a restart ``breaks.dat`` and ``tower_failures.dat`` are
written again from the checkpoint's records, so they hold the releases and
failed towers up to the checkpoint, also in a new diagnostics directory, and
only their headers when an option is first set at the restart. A run from the
start without ``failure_break`` removes ``breaks.dat`` of an earlier run in
the same directory, and one without ``tower_failure = release`` removes
``tower_failures.dat``.

Verification of the release
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Test
     - What is pinned
   * - ``DecideBreak.*``
     - The table above and Eq. :eq:`br-end`: the higher end, end a within
       :math:`\epsilon` and end b beyond it, the one end still attached,
       nothing when none is; refused indices.
   * - ``BreakLog.*``
     - The rows of ``breaks.dat``, the releases of one line in order, a
       bit-exact round trip of the checkpoint text, refused records.
   * - ``LineRelease.AReleasedSpanHangsFromItsOtherEndAndItsSupportTakesNothing``
     - Eqs. :eq:`br-rest` and :eq:`br-lumped`: 31 m of conductor in 1 m
       segments from a support 30 m above flat ground, released at its low
       end, in still air. :math:`w (h - d/2)` = 478.68 N: the stub gives it
       exactly; MoorDyn-C must lie within one segment's weight below it
       (478.9 N while the node at the foot still hovers above the ground,
       462.9 N once it rests on it). One more kinematics point; the released
       support takes nothing at once.
   * - ``LineRelease.ASpanShorterThanItsSupportsHeightHangsFree``
     - 21.35 m of conductor from a support 30 m up, released at its low end,
       hangs clear of the ground: the support carries :math:`w L` = 340.8 N
       (stub exactly; MoorDyn-C 340.8 N, within 2 %) and the free end stays
       clear of the ground (8.65 m above it at rest, where the stub places it).
   * - ``LineRelease.AStateAfterAReleaseFitsOnlyALineThatMakesItAgain``
     - A state saved at the release and one saved ten steps later both
       continue in a new line that makes the release again; without it the
       load is refused; a step redone from an in-memory state is identical.
   * - ``LineRelease.AClampedTowerKeepsOnlyTheSpanStillOnIt``, ``LineRelease.AReleasedStringDropsItsConductorAndItsTowerTakesNothing``
     - A tower's pull is that of what is still attached: one span's horizontal
       tension and half its weight, then nothing; a released string's tower
       takes nothing and the next tower takes more.
   * - ``LineRelease.ASpanOnStringsLeavesTheStringsLowerEnd``, ``LineRelease.AReleaseAtAMovingTowersCrossArmContinuesFromASavedState``
     - On a line with strings a span's end at a tower is released from the
       string's lower end, and the string stays on its tower. At a tower that
       moves, MoorDyn's force on the cross-arm is the pull of what is still on
       it, and a state saved after the tower has moved on continues in a new
       line.
   * - ``LineRelease.ReleasingTwiceOrOutsideTheLineIsRefused``
     - A second release of the same end or string, and an index outside the
       line, stop naming the line.
   * - ``Conductors.FailureBreakLetsGoAtTheHigherEndAndTheRunGoesOn``
     - An inclined span, end b 10 m higher: with MoorDyn-C the tension at b
       is higher by :math:`w \times` 10 m and b lets go, at the end of the
       step of the exceedance, once; the span comes to lie on the ground.
   * - ``Conductors.ASpanIsReleasedOnceInAStepAndADeadEndReleasesItsOwnEnd``
     - Three exceedances in one step release the span once; the hardware at a
       dead end releases that end.
   * - ``Conductors.ReleasesContinueAcrossARestart``
     - Checkpoints at the step of the release and four steps later restart
       into the uninterrupted run's positions, also without
       ``failure_break``; without ``ground_contact`` the restart is refused.
       A restart that first sets ``failure_break`` after the exceedance
       releases nothing.
   * - ``Conductors.StringsOverTheirLoadLetGoAndTheirTowersTakeNothing``, ``Conductors.LinesOnMovingTowersLetGoAtTheEndOfTheLevelsStep``
     - Through the whole model: both strings of a section let go and their
       standing towers carry no line load; on towers that bend, exceedances
       found inside the step are released at its end, one release per span.
   * - ``Conductors.FailureBreakNeedsGroundContactALineCapacityAndNoFailureStop``
     - The three refusals; off by default.
   * - ``Conductors_Break_TwoLevel`` (CTest, stub and MoorDyn-C)
     - A level 300 m span 30 m above flat ground in a sheared crosswind, the
       line on a refined level, strength 5 kN against 13.4 kN: released at
       0.125 s, the first step of the lines' level, with either library; the
       released end then carries next to nothing, the span reaches the ground
       (MoorDyn-C: at 2.125 s), and after 6 s the held end carries less than
       half the strength. The stub's held end must carry Eq. :eq:`br-rest`,
       478.71 N, which it places by construction. MoorDyn-C's span, in 15 m
       segments, still swings then (476 N in this run; Eq. :eq:`br-lumped`
       gives 241 N at rest).
   * - ``Conductors_Break_MidRun`` (CTest, stub)
     - The strength of ``Conductors_Capacities_TwoLevel``: the release is at
       1.875 s, between two steps of the coarse level.
   * - ``Conductors_Break_Restart``, ``Conductors_Break_Terrain_Restart`` (CTests, stub and MoorDyn-C), ``Conductors_Break_AfterRestart`` (stub)
     - The restarted runs' ``breaks.dat``, ``failures.dat``, ground records
       and span, tower and transformer logs equal the uninterrupted runs' to
       ten digits. On the hills both dead ends of a line let go at 0.9 s;
       with MoorDyn-C the swing then overloads the hardware at a tower, whose
       string lets go at 4.5 s, the step of the checkpoint.
   * - ``Conductors_Break_NeedsGroundContactAbort``, ``_WithFailureStopAbort``, ``_NeedsLineCapacityAbort`` (CTests)
     - The three refusals at start-up.
   * - ``Conductors.CapacityChecksInsideTheStepStampTheFirstExceedanceThere``
     - Eq. :eq:`ci-pieces`: a step of 0.25 s with :math:`\tau` = 0.1 s is
       checked at 0.25/3 s, where a strength below the still-air tension is
       first exceeded; the release is at 0.25 s. Lines on towers that bend
       take 25 coupling steps with :math:`\tau` = 0.01 s, 10 without.
   * - ``Conductors.CapacityChecksInsideTheStepSeeTheLoadAfterARelease``
     - A clamped section lets go of both dead ends: the largest pull recorded
       on a tower with a check every 0.01 s is at least that with one check
       per 0.25 s step, and more than 5 % above it with MoorDyn-C.
   * - ``Conductors.CapacityCheckIntervalIsRefusedOutOfRangeOrWithoutACapacity``
     - Values outside [0, 3600], no capacity to check, and more than the
       most checks per step stop the run, naming the input.
   * - ``Conductors.CapacityChecksStartAtTheirStartTime``
     - Eq. :eq:`cs-start`: with :math:`t_s` = 0.6 s and steps of 0.25 s a
       strength below the still-air tension is first exceeded at 0.75 s,
       not at 0.25 s; before it the component has no record and the line is
       whole. With pieces of 0.25/3 s the first check is at 0.6667 s. A
       negative start and a start without a capacity stop the run.
   * - ``TowerReleases.*``, ``TowerFailureLog.*``
     - The table of `Tower failure`_: the string, or both sides of a clamped
       line, less what was released before; the rows of
       ``tower_failures.dat``, a bit-exact round trip, a second failure of a
       tower refused.
   * - ``Conductors.ATowerOverItsFootingsLoadLetsGoOfItsLines``
     - Standing towers whose legs carry 23.7 kN against a bearing capacity
       of 23 kN fail at the end of the first step, Eq. :eq:`tf-time`: a
       clamped line leaves each on both sides, a line on strings loses its
       strings; each tower is logged once and carries no line afterwards.
   * - ``Conductors.AFailedTowerPassesItsLoadToTheNext``
     - Spans of 300, 500 and 100 m on strings, 60 m up, legs rated 1.8 kN:
       tower 1's legs carry 1.98 kN and tower 2's 1.52 kN, so only tower 1
       fails at the first step (0.25 s). With MoorDyn-C tower 2, whose legs
       then carry 2.08 kN at rest and more while the conductor falls, fails
       at 1.0 s; without the release it never passes its capacity. The
       stub's string load hardly changes, so there the second tower stands.
   * - ``Conductors.FailedTowersStayFailedAcrossARestart``, ``Conductors.TowerFailureNeedsFailureBreakAndATowersCapacity``
     - Failed towers and their releases come back from a checkpoint and the
       lines follow the uninterrupted run's; the refusals; off by default.
   * - ``Conductors_TowerFailure_Terrain_Restart`` (CTest, stub and MoorDyn-C), three abort CTests
     - On the hills all three towers fail at the first step under a bearing
       capacity of 10 kN; the restarted run's logs equal the straight run's
       to ten digits. The refusals at start-up.

**Observed, not pinned.** In the MoorDyn-C run of ``Conductors_Break_TwoLevel``
the falling span's largest tension after the release is 7.7 kN at 3 s, 57 % of
the tension before it, as the span lands and pulls taut on its remaining
support. In 27 further MoorDyn-C runs of this span alone (10, 20 and 40
segments; internal steps of 0.1, 0.2 and 0.5 ms; crosswinds of 0, 15 and
30 m/s) every run stayed finite, the peak was 41 to 59 % of the tension
before the release, 1.1 to 2.5 s after it, and 15 to 105 ms wide at half its
height. A check every 0.125 s therefore sees 41 to 67 % of such a peak at the
worst phase, and one every 0.9 s 10 to 27 %. On the hills the hardware at tower 1, rated 15 kN, sees 17.3 kN at
4.5 s after the dead ends let go at 0.9 s. These are single runs on coarse
lines (15 m and 16 m segments).

References
~~~~~~~~~~

* MoorDyn-C v2 documentation, C API (``MoorDyn_BreakLine``).
  https://moordyn.readthedocs.io

No published formula or data set is used in this section: Eq. :eq:`br-rest`
is a balance of forces at rest, and the rule of Eq. :eq:`br-end` is this
model's own.

Conductor temperature under a fire
----------------------------------

This section describes how the model gives each conductor a temperature that
follows the heat it gains from a fire, and how the conductor's sag follows that
temperature. It is the first stage of the fire-to-structures coupling. The
towers' steel temperature is the second (see `Tower steel temperature under a
fire`_). The implementation is
``Structures/Conductors/ERF_LineThermal.{H,cpp}`` (class ``LineThermal``),
called from ``Conductors`` once per step of the lines' anchor level. It is
selected with ``conductors.thermal_model = ieee738``; without it nothing
changes.

Heat balance
~~~~~~~~~~~~

Every node :math:`n` of every span (insulator strings excluded) carries a
temperature :math:`T_n` (°C). It follows the heat balance of IEEE Std 738-2023
(eq 3) per unit length of conductor,

.. math::
   :label: lt-balance

   m c_p \frac{dT_n}{dt} = q_J(T_n) + q_s + q_{f,n} - q_c(T_n, T_{a,n}, V_n, \phi_n) - q_r(T_n, T_{s,n}) ,

where :math:`m c_p` (J/(m K)) is the conductor's heat capacity per unit length
(``heat_capacity``), :math:`q_J` the Joule heating, :math:`q_s` the solar
heating, :math:`q_{f,n}` the radiation absorbed from the fire, :math:`q_c` the
convective loss and :math:`q_r` the radiative loss, all in W/m. :math:`T_{a,n}`
is the air temperature at the node, :math:`V_n` the wind speed there,
:math:`\phi_n` the angle between the wind and the conductor's axis, and
:math:`T_{s,n}` the temperature of the surroundings the conductor radiates to.
IEEE 738 takes :math:`T_s = T_a`, the ambient air. Here :math:`T_{s,n}` is
the node's air temperature at the first step and stays fixed: inside a hot
plume the air is a poor emitter, and radiating to it as to a black enclosure
would overstate the heating (by 873 W/m for a conductor at 100 °C in air at
500 °C, more than the fire's own radiation in the test case below).

.. list-table::
   :header-rows: 1
   :widths: 18 52 30

   * - Term
     - Form (SI)
     - IEEE Std 738-2023
   * - :math:`q_J`
     - :math:`I^2 R(T)`, :math:`R` linear through the two given resistances
     - eq 33
   * - :math:`q_s`
     - :math:`\alpha D E_s`, with :math:`E_s` the sun's flux normal to the conductor
     - eq 23, with :math:`E_s` given
   * - :math:`q_c`
     - the larger in magnitude of forced convection,
       :math:`K_\phi\,\max[(1.01 + 1.35 Re^{0.52}),\ 0.754 Re^{0.6}]\, k_f (T - T_a)`,
       and natural convection, :math:`3.645 \rho_f^{1/2} D^{3/4} |T - T_a|^{5/4}`
       with the sign of :math:`T - T_a`
     - eqs 6, 7, 9, 10, 19
   * - :math:`q_r`
     - :math:`17.8\, D\, \varepsilon \left[ \left(\frac{T + 273}{100}\right)^4 - \left(\frac{T_s + 273}{100}\right)^4 \right]`
     - eq 21, with :math:`T_s` for :math:`T_a`

Here :math:`D` is the diameter (m), :math:`\alpha` the absorptivity and
:math:`\varepsilon` the emissivity of the conductor's surface,
:math:`Re = D \rho_f V / \mu_f` the Reynolds number and
:math:`K_\phi = 1.194 - \cos\phi + 0.194\cos 2\phi + 0.368\sin 2\phi` the wind
direction factor. The air's viscosity :math:`\mu_f`, density :math:`\rho_f`
(at ``conductors.thermal.elevation``) and conductivity :math:`k_f` are taken at
the film temperature :math:`(T + T_a)/2` (IEEE Std 738-2023 eqs 12, 13, 15, 17).
IEEE Std 738-2023 (4.4.6) interpolates the resistance between its values at
25 and 75 °C up to 100 °C, and between 25 and 200 °C above it. A run whose
conductors pass 100 °C under a fire should therefore give
``resistance_temperatures = 25 200`` with the resistance at 200 °C.
IEEE 738 states the convection correlations for a conductor hotter than the
air. When the air is hotter, as inside a plume, the same correlations are used
with the sign of :math:`T - T_a`, so heat flows into the conductor.

Radiation from the fire
~~~~~~~~~~~~~~~~~~~~~~~

The fire reaches the model through ERF's fire view (``StructureStep::fire``,
:ref:`sec:StructureCoupling`). Each fire cell :math:`j` with a positive heat
release :math:`q''_j` (W/m\ :sup:`2`, ``FireEngine::heat_release``) is treated
as a point source of radiant power :math:`\chi_r Q_j`, with
:math:`Q_j = q''_j\,\Delta A` the cell's heat release and :math:`\chi_r` the
radiative fraction (``conductors.thermal.radiative_fraction``). The source sits
at the cell centre, ``conductors.thermal.source_height`` above the terrain. A
conductor element at node :math:`n` intercepts radiation over its projected
width :math:`D \sin\psi_{nj}`, where :math:`\psi_{nj}` is the angle between the
ray from the source and the conductor's axis. The absorbed heat is

.. math::
   :label: lt-fire

   q_{f,n} = \alpha D E_n , \qquad
   E_n = \sum_j \frac{\chi_r Q_j \sin\psi_{nj}}{4 \pi \max(R_{nj}^2,\ \Delta A/\pi)} ,

where :math:`R_{nj}` is the distance from the source to the node. Each source
radiates equally in every direction, so its flux falls with the square of the
distance: :math:`\chi_r Q_j` spread over a sphere of radius :math:`R_{nj}`. The floor
:math:`\Delta A/\pi` is the squared radius of a disc of one fire cell's area. It
stops the sum from growing without bound when a node comes within a cell of a
source, where the point-source model does not hold anyway. Sources farther
from the nodes than ``conductors.thermal.cutoff_radius`` are left out; a run
with a fire must give it, since every burning cell within it of the lines is
gathered to every rank at every step. In a periodic direction the distance is
the nearest image's. The conductor's axis at a node comes from its
neighbouring nodes, so a span that sags or swings changes its exposure.

The burning cells are gathered from every rank once per step and sorted by
position, so the sums do not depend on the box decomposition. On a
terrain-following mesh the terrain under a source is read from the lines'
anchor level, which must therefore cover every source within the cutoff.

Air temperature and wind
~~~~~~~~~~~~~~~~~~~~~~~~

The air temperature at each node comes from the state of the lines' anchor
level at the start of the step, through ERF's equation of state for dry air
(:math:`\rho` and :math:`\rho\theta`). The wind is the one the lines already
take at their nodes. With ``erf.fire.coupling_type`` lagged or synchronous the
fire heats the air, so a plume that rises through the span heats the
conductor by convection as well as by radiation. With passive coupling only
the radiation reaches it.

Time integration and timing
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Equation :eq:`lt-balance` is stepped with backward Euler over the anchor
level's step :math:`\Delta t`, with the air, the wind and the fire held at
their values at the start of the step:

.. math::
   :label: lt-be

   m c_p \left(T_n^{k+1} - T_n^{k}\right) = \Delta t \; G\!\left(T_n^{k+1}\right) ,

where :math:`G` is the right-hand side of :eq:`lt-balance`. The new
temperature is found by bracketing and bisection to :math:`10^{-12}` of its
value. Over the range of the inputs the equation has one root, although
:math:`m c_p (T - T^k) - \Delta t\, G(T)` need not be monotone at steps of
order :math:`10^4` s and high current, where the Joule term's growth with
:math:`T` can exceed the losses' near :math:`T_a`. Backward Euler is stable
for any step and tends to the steady temperature for long steps. It is first
order: against a fine-step solution its error is about
:math:`0.5\,\Delta t/\tau` times the temperature change, with :math:`\tau`
the conductor's time constant (about 0.5 K at :math:`\Delta t` = 1 s for
:math:`\tau \approx` 120 s in a hot plume; 0.04 K for :math:`\tau \approx`
650 s under radiation alone). On the first step of a run each node starts at
its steady temperature without the fire, :math:`G(T) = 0` with
:math:`q_f = 0`.

The fire is read at the start of the lines' step, and the fire steps after
the dycore of its own level. When the fire is on the lines' level or a finer
one, its fields therefore belong to the lines' step time exactly (the end of
its last step), and on the run's first step the lines see no fire yet. When
the fire is on a coarser level, its fields belong to the end of the coarse
step that contains the lines' step (:ref:`sec:StructureCoupling`).

Thermal sag
~~~~~~~~~~~

Each span's unstretched length follows the mean temperature of its conductor:

.. math::
   :label: lt-length

   L_k = L_{k,\mathrm{ref}} \left[ 1 + \alpha_T \left( \bar T_k - T_{k,\mathrm{ref}} \right) \right] ,

where :math:`\alpha_T` is the thermal expansion coefficient
(``thermal_expansion``), :math:`\bar T_k` the trapezoid mean of the span's node
temperatures, :math:`L_{k,\mathrm{ref}}` the span's input unstretched length
and :math:`T_{k,\mathrm{ref}}` the temperature it belongs to:
``conductors.thermal.reference_temperature`` when given (a stringing
temperature, say), otherwise the span's mean steady temperature at the first
step, so that the span starts at its input length. MoorDyn takes the new length
before the lines step (``MoorDyn_SetLineUnstretchedLength``), so the span
lengthens and sags, and its tension falls.

A span that sags comes closer to the fire, which raises its irradiance. With
the geometry taken from the step's start this feedback is stable at any step
as long as its loop gain, the rise in irradiance per metre of sag times the
temperature rise per unit irradiance times the sag per kelvin, stays below
one (about 0.14 in the test case); above one the line would sag into the
flames, which is a physical runaway, not a numerical one.

Outputs and restart
~~~~~~~~~~~~~~~~~~~

``<diagnostics_dir>/conductor_temperature.dat`` has one row per span and
output step: the time, the line and span, the mean temperature, the hottest
node's temperature and index, that node's air temperature, wind speed, wind
angle and fire irradiance, the largest irradiance on the span, the thermal
strain :math:`\alpha_T(\bar T_k - T_{k,\mathrm{ref}})` and the reference
temperature; the hottest node is numbered from 0 over all the line's nodes.
Rows follow ``conductors.diagnostics_int``. The checkpoint holds every node's
temperature and surroundings temperature and each span's mean, reference
temperature and reference length (``<chk>/conductors/thermal``). A restart
restores them and sets the spans' lengths again, since MoorDyn-C's saved state
holds the segment lengths but not the unstretched length. A checkpoint without
them starts the temperatures afresh at the restart, with the spans' lengths as
they are then; a run without the thermal model that restarts from a
checkpoint with them says so, and its spans return to their input lengths.

Limitations
~~~~~~~~~~~

* The point-source model places all of a cell's radiation at one point and
  radiates it equally in every direction. It does not model the flame's shape,
  its tilt in the wind, or smoke that absorbs radiation between the flame and
  the conductor. With ``source_height`` 0 half of each source's power goes
  into the ground: over the test case's disc, an upward (Lambertian) emitter
  of the same power gives 1.8, 2.7 and 3.2 times the irradiance at 5, 17 and
  30 m above it. Take ``source_height`` near half the flame height, and read
  the irradiance as uncertain by a factor of two to three near the fire.
* :math:`\chi_r` and ``source_height`` are inputs. :math:`\chi_r` ranges from
  about 0.15 to 0.45 with the fuel, and the radiative heating scales with it
  directly, so it has no default.
* The air's density in the convection terms comes from IEEE 738's elevation
  formula, not from ERF's sampled density.
* Convection along the conductor at low wind uses the wind angle factor as
  given (0.388 for a wind along the axis). At low wind the convection is the
  larger of the forced and natural values, as IEEE Std 738-2023 (eq 6)
  recommends, not a combination of the two.
* The air temperature uses the equation of state for dry air. Moisture
  changes it by well under 1 K.
* The conductor is one temperature across its section: no radial gradient
  between the aluminium strands and the steel core. Annealing and loss of
  strength above about 93 °C (aluminium) are not modelled; the temperatures are
  reported so that they can be checked against such limits.
* One expansion coefficient serves at every temperature. An ACSR conductor
  passes a knee point (about 70 to 120 °C) above which the aluminium carries no
  tension and the steel core's smaller expansion and stiffness govern, so above
  it a single :math:`\alpha_T` overstates the sag.
* MoorDyn-C keeps each segment's volume proportional to its length, so a
  lengthened span gains mass in proportion (about 1 % at +500 K); the bundled
  stub keeps the mass.
* IEEE Std 738-2023 (Table 2) states the forced convection fits for air at 15
  to 260 °C, conductors at 21 to 1004 °C and diameters of 0.01 to 150 mm. Plume
  air hotter than 260 °C, or air colder than 15 °C, lies outside that range,
  and so does air hotter than the conductor.
* Insulator strings are not heated. The towers' steel is, with
  ``conductors.thermal.towers`` (see `Tower steel temperature under a fire`_).

Verification of the thermal model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Test
     - What is pinned
   * - ``LineThermal.TheIEEE738TermsMatchTheIndependentCalculation``
     - The air's properties, :math:`K_\phi`, forced and natural convection and
       radiation at :math:`D` = 28.1 mm, :math:`V` = 0.61 m/s across,
       :math:`T` = 100 °C, :math:`T_a` = 40 °C, against an independent Python
       calculation (``check_fire_heating.py``), to :math:`10^{-12}`.
   * - ``LineThermal.TheTermsMatchIEEE7382023sSampleCalculations``
     - The Drake conductor of IEEE Std 738-2023 (4.6.1 to 4.6.3): the air's
       properties, natural and forced convection, radiation, solar heating,
       :math:`R` (100 °C), the 1025 A rating, 97.5 °C at 1000 A and the first two
       10 s steps after a step to 1200 A, to the standard's printed digits.
   * - ``LineThermal.TheSteadyAndImplicitTemperaturesMatchTheIndependentCalculation``
     - Steady temperatures with 1000 A, with the sun and with a fire, and
       backward-Euler steps, to :math:`10^{-9}` K; without heat sources the
       steady temperature is the air's.
   * - ``LineThermal.APointSourceUnderALongConductorGivesTheLineIntegral``
     - Over a straight conductor at height :math:`H` above one source,
       :math:`\int E\, dx = \chi_r Q/(4\pi) \cdot 2a/(H\sqrt{a^2 + H^2})`
       to :math:`10^{-5}`.
   * - ``LineThermal.TheSpansLengthFollowsItsMeanTemperatureAndSurvivesARestart``
     - On the stub, MoorDyn's unstretched length equals :eq:`lt-length`; a
       checkpoint restores the temperatures to the bit and sets the length again.
   * - ``Conductors_FireHeating`` and its two-level variants (CTests, stub)
     - A 2e5 W/m\ :sup:`2` disc of 40 m radius under a 300 m span in still air,
       passive coupling. An independent calculation of every span node: the air
       temperature from the sounding (:math:`10^{-4}` K), the irradiance from the
       fire cells and node positions (:math:`10^{-6}`), the temperatures with its
       own IEEE 738 terms (hottest node and span mean, :math:`10^{-3}` K), the
       reference temperature, the strain of every row and the final sag from the
       stub's elastic parabola (:math:`10^{-4}` m). One level;
       the lines on a level-1 patch with the fire on level 0; the fire on a
       level-1 patch with the lines on level 0. The first step sees no fire
       except when the fire is on the coarser level.
   * - ``Conductors_FireHeating_Restart`` (CTest, stub)
     - A checkpoint at 10 s: the temperature log and the span's log equal the
       uninterrupted run's.
   * - ``LineThermal.TheConductorRadiatesToItsSurroundingsNotToThePlume``
     - Plume air at 500 °C with surroundings at 26.7 °C: the heat gain and a
       step against the independent calculation.
   * - ``LineThermal.APeriodicSourceIsSeenAcrossTheSeam``,
       ``LineThermal.TheWindsAngleToTheSpanAndTheReferenceTemperatureAreUsed``
     - The nearest image across a periodic boundary; a wind along and across
       the span gives the steady temperatures of :math:`\phi` = 0 and 90°;
       ``reference_temperature`` sets the span's starting length.
   * - ``Conductors_Thermal_*Abort`` (CTests)
     - An unknown model, a missing property, a value out of range and a fire
       without ``radiative_fraction`` stop the run, naming the key.

References
~~~~~~~~~~

* IEEE (2023). *IEEE Standard for Calculating the Current-Temperature
  Relationship of Bare Overhead Conductors*, IEEE Std 738-2023 (revision of
  IEEE Std 738-2012). IEEE, New York. ISBN 978-1-5044-9820-3.
  doi:10.1109/IEEESTD.2023.10382442. The equation
  numbers above are this edition's.

Tower steel temperature under a fire
------------------------------------

This section describes how the model gives each member of a lattice tower a
steel temperature that follows the heat it gains from a fire, and how the
tower's frame softens as the steel heats. It is the second stage of the
fire-to-structures coupling, after `Conductor temperature under a fire`_. The
implementation is ``Structures/Conductors/ERF_TowerThermal.{H,cpp}`` (class
``TowerThermal``), called from ``Conductors`` once per step of the anchor level,
after the conductors' temperatures and before the lines and towers step. It is
selected with ``conductors.thermal.towers = true``, which needs
``conductors.thermal_model = ieee738`` (the fire's sources and the air
temperature come from the same pass) and at least one tower that bends as a
frame model (``frame_file`` or ``frame_panels``, see :ref:`sec:TowerFrame`). Every
such tower is heated. Without it nothing changes.

Heat balance of a member
~~~~~~~~~~~~~~~~~~~~~~~~

Each frame member :math:`m` carries one temperature :math:`\theta_m` (°C). It
follows the lumped heat balance of unprotected steel of EN 1993-1-2 (its
relation 4.25, as the JRC report of Vassart et al. (2014), section 2.5.2.2,
gives it), per unit length of member,

.. math::
   :label: tt-balance

   \rho_a A_m \frac{dH(\theta_m)}{dt} =
     P_{x,m} \left[ \alpha_c \left( T_{a,m} - \theta_m \right)
     + \varepsilon\, \sigma \left( T_{s,m}^4 - \theta_m^4 \right) \right]
     + \varepsilon\, W_m E_m ,

where :math:`\rho_a` = 7850 kg/m\ :sup:`3` is the steel's density,
:math:`A_m` (m\ :sup:`2`) the member's area, :math:`H` (J/kg) the steel's
enthalpy, :math:`P_{x,m}` (m) the perimeter that exchanges heat with the air,
:math:`\alpha_c` (W/(m\ :sup:`2` K)) the convection coefficient
(``conductors.thermal.tower_convection``), :math:`\varepsilon` the steel's
emissivity (``conductors.thermal.tower_emissivity``, 0.7 by default, also its
absorptivity; the JRC report, section 2.5.2.2), :math:`\sigma` the Stefan-Boltzmann constant, :math:`W_m` (m)
the width that intercepts the fire's radiation and :math:`E_m` (W/m\ :sup:`2`)
the fire's irradiance. The fourth powers are of absolute temperatures.
:math:`T_{a,m}` is ERF's air temperature along the member and
:math:`T_{s,m}` the temperature of the surroundings it radiates to: its air
temperature at the first step, held fixed, as for the conductors.

The enthalpy is the integral of the specific heat of carbon steel of
EN 1993-1-2 (its section 3.4.1.2, which the JRC report plots in Fig. 2.4.7),

.. math::
   :label: tt-enthalpy

   H(\theta) = \int_{20}^{\theta} c_a(s)\, ds , \qquad
   c_a(\theta) =
   \begin{cases}
     425 + 0.773\,\theta - 1.69 \times 10^{-3}\,\theta^2 + 2.22 \times 10^{-6}\,\theta^3 , & 20 \le \theta < 600 , \\
     666 + 13002 / (738 - \theta) , & 600 \le \theta < 735 , \\
     545 + 17820 / (\theta - 731) , & 735 \le \theta < 900 , \\
     650 , & 900 \le \theta ,
   \end{cases}

in J/(kg K), integrated in closed form. Below 20 °C the 20 °C value is used.
The peak of 5000 J/(kg K) at 735 °C is the steel's phase change.

**Section.** The section enters through :math:`A_m`, :math:`P_{x,m}` and
:math:`W_m`:

.. list-table::
   :header-rows: 1
   :widths: 30 23 22 25

   * - Member
     - :math:`A_m`
     - :math:`P_{x,m}`
     - :math:`W_m`
   * - Equal-leg angle :math:`b \times t` (its design data)
     - :math:`t(2b - t)`
     - :math:`4b`
     - :math:`[(2 + \sqrt 2)\, b + (2 - \sqrt 2)\, t]/\pi`
   * - Circular, diameter :math:`D`, wall :math:`t` (0: solid)
     - :math:`\pi [D^2 - (D - 2t)^2]/4`
     - :math:`\pi D`
     - :math:`D`
   * - Rectangular :math:`S_a \times S_b`, wall :math:`t` (0: solid)
     - :math:`S_a S_b - (S_a - 2t)(S_b - 2t)`
     - :math:`2(S_a + S_b)`
     - :math:`2(S_a + S_b)/\pi`

:math:`P_{x,m}` is EN 1993-1-2's box perimeter, the shadow factor
:math:`k_{sh}` times the section's perimeter; for an angle the box and the
section's own perimeter are both :math:`4b`, so :math:`k_{sh}` = 1.
:math:`W_m` is the mean projected width of the section's convex hull, its
perimeter over :math:`\pi` (Cauchy's formula); an angle's hull is its two outer
faces, its two toe edges and the chord between the toes. For a tube it is the
diameter, so a tube intercepts the fire as a conductor does.

**A hybrid of EN's convention and the geometry.** The heat exchanged with the
air uses EN 1993-1-2's box perimeter, as the standard does, while the fire's
beam is intercepted by the geometric width. For a 100 × 8 mm angle in air at
26.85 °C, the steady temperatures under 1e4, 5e4 and 2e5 W/m\ :sup:`2` are
89.4, 278.4 and 611.2 °C as coded; 138.2, 364.9 and 677.2 °C with
:math:`\alpha_c` = 10 W/(m\ :sup:`2` K); 91.0, 287.7 and 637.1 °C with the
hull's perimeter, :math:`(2 + \sqrt 2) b + (2 - \sqrt 2) t`, for the
radiation exchange; and 98.5, 306.3 and 653.2 °C with the box's width
:math:`4b/\pi` for the fire. EN's choices give the cooler steel. A member with an
arbitrary SubDyn section has no perimeter in the frame file; it needs its
design data as an angle (``member_file``), or the run stops.

**Why EN 1993-1-2.** It is the standard method for the temperature of
unprotected structural steel in a fire, and its specific heat and reduction
factors are the ones the member checks already use (`Moving towers`_). Its
heating method, its reduction factors (Table 3.1) and a worked example are
reproduced in the open-access JRC report of Vassart et al. (2014), section 2.5.2
and Table 2.4.1, against which the model is checked (see below). Its
section factor :math:`P_{x,m}/A_m` (260 m\ :sup:`-1` for a 100 × 8 mm angle)
sets how fast a member heats: thin bracing heats faster than the legs.

Fire radiation and air temperature
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fire's irradiance is that of the conductors, :eq:`lt-fire`: the point
sources at every burning fire cell within ``cutoff_radius``, at
``source_height`` above the terrain, with the radiative fraction
:math:`\chi_r`, projected on the member's axis (:math:`\sin\psi`) and with the
same floor on the squared distance. :math:`E_m` and :math:`T_{a,m}` are the
means over three Gauss-Legendre points along the member (at
:math:`1/2 \mp \sqrt{3/5}/2` and :math:`1/2` of its length, weights 5/18,
8/18 and 5/18). The points stand on the frame as built, in ERF's axes: the
tower's base plus the frame's coordinates turned so that its :math:`x` axis
runs along the line and its :math:`y` axis along the cross-arm. The air
temperature comes from ERF's equation of state at those points, so a plume
(lagged or synchronous coupling) heats the steel by convection as well.

Time integration and the frame
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Equation :eq:`tt-balance` is stepped with backward Euler in the enthalpy, with
the air and the fire held at their values at the start of the step:

.. math::
   :label: tt-be

   \rho_a A_m \left[ H\!\left(\theta_m^{k+1}\right) - H\!\left(\theta_m^{k}\right) \right]
     = \Delta t\; G_m\!\left(\theta_m^{k+1}\right) ,

where :math:`G_m` is the right-hand side of :eq:`tt-balance`. The left side
rises and :math:`G_m` falls with :math:`\theta`, so the root is unique; it is
bracketed and found by bisection to :math:`10^{-10}` K. Stepping the enthalpy
rather than the temperature keeps the energy of every step exact across the
peak of :math:`c_a`, at any step. Backward Euler tends to the steady
temperature for long steps, short of it by :math:`(\tau/\Delta t)` times the
remaining change when :math:`\Delta t \gg \tau` (by :math:`1/(1 + \Delta t/\tau)` of
it for the linearised balance), with :math:`\tau` the member's time constant
:math:`\rho_a A_m c_a / [P_{x,m} (\alpha_c + 4 \varepsilon \sigma \theta^3)]` (about 450 s for a
100 × 8 mm angle near ambient, 350 s at 250 °C). On the first step every member starts at its
air temperature, the steady state without the fire.

**Rebuilding the frame.** The frame's stiffness and the members' strengths
follow EN 1993-1-2's reduction factors :math:`k_E(\theta)` and
:math:`k_y(\theta)` (Table 3.1). The tower's frame is rebuilt at its members'
temperatures (``FrameTower::set_temperature``) when, for any member and for
:math:`k = k_y` or :math:`k_E`,

.. math::
   :label: tt-rebuild

   \frac{\left| k(\theta_m) - k(\theta_m^{f}) \right|}{\max\left[ k(\theta_m),\ k(\theta_m^{f}) \right]}
   > \delta_k ,

where :math:`\theta_m^{f}` is the temperature the frame was last built at and
:math:`\delta_k` is ``conductors.thermal.tower_reduction_tolerance``. The
change is relative, so the frame's stiffness lags the steel's by at most the
fraction :math:`\delta_k` at any temperature; an absolute change would let it lag
by 8 % of the remaining stiffness at 650 °C and 25 % above 900 °C for
:math:`\delta_k` = 0.01. The rebuild re-solves the frame's static equilibrium
under its weight, its first natural frequency and the Newmark factorisation,
and carries the frame on from the same position and velocity. Both factors are
1 up to 100 °C, so no rebuild happens below it. The members' strength checks
do not wait for a rebuild: they take the steel's temperature of the step
(``FrameTower::set_check_temperature``), under the forces of the frame as last
built. The coupling's step count follows the frame's lower frequency. A rebuild
costs a frame assembly, an eigenvalue solve and a static solve; between 500 and
700 °C, where the factors fall fastest, a member heating by a few kelvin per
second can rebuild its frame every step.

Outputs and restart
~~~~~~~~~~~~~~~~~~~

``<diagnostics_dir>/tower_temperature.dat`` has one row per heated tower and
output step: the time, the tower, the members' mean and largest temperature,
the member id of the largest, the largest member irradiance, the least
:math:`k_E` of the frame, the rebuilds so far and the frame's largest member
temperature :math:`\theta^f`. ``tower_<tower>_member_temperature.dat`` has
every member's temperature, in the frame's member order (the ids are in its
header). Rows follow ``conductors.diagnostics_int``. The checkpoint holds every
member's temperature and surroundings temperature and the rebuild count
(``<chk>/conductors/tower_thermal``); the frame's own temperatures come back
with the tower's state, which rebuilds it hot, and the members' checks take the
member temperatures again. A restart without ``conductors.thermal.towers`` from
a checkpoint whose frame was rebuilt hot stops, naming the key. From a
checkpoint of a run with the heating whose frames were not yet rebuilt, it
continues with a warning: the members stay at their type's
``steel_temperature``, and the earlier run's tower temperature logs are
removed. A restart with the heating from a checkpoint without it (written
before the heating was turned on) starts the member temperatures at the
restart, at the air temperature then, which also becomes their surroundings'
temperature, and writes a row at its first step. The ``tower_*`` settings
without ``conductors.thermal.towers`` are not read, and the run says so.

Limitations
~~~~~~~~~~~

* One temperature per member: no gradient along a member or across its
  section. Long legs near a fire are hotter at one end than the other; the
  Gauss mean of the irradiance sets the member's one temperature.
* The members' thermal expansion is not applied to the frame: heating lowers
  the stiffness and the strength but does not lengthen a member, so the forces
  from restrained expansion are missing. The members are checked against
  their reduced strengths, but creep at high temperature is not modelled.
* No member shades another, and the tower does not shade the conductors; the
  fire's radiation reaches every member as if it stood alone.
* The members stand where the frame was built; their displacement (of order
  centimetres) does not move them toward or away from the fire. The frame
  stands on a flat footing at the height of the tower's centre, so on a slope
  part of a leg may lie below the terrain; such a point takes the air at the
  terrain, and its irradiance is taken where it stands.
* :math:`\alpha_c` is a constant (25 W/(m\ :sup:`2` K) by default, the value
  for the standard fire; the JRC report, section 2.5.2.2). It does not follow the wind. For a member heated
  by the fire's radiation in cooler air it overstates the convective cooling
  in calm air, which gives cooler steel (see the hybrid paragraph above).
* The specific heat's piecewise formulas have not been checked against the text
  of EN 1993-1-2, which we do not have; the JRC report plots the curve (Fig.
  2.4.7, its peak of 5000 J/(kg K) near 735 °C), and its worked example with the
  varying specific heat (667 °C in about 17 minutes) agrees with the model to
  6 s.
* The default emissivity is the JRC report's 0.7 for steel; a galvanised or
  weathered surface may differ. Set ``tower_emissivity`` for the surface.
* The point-source model's limits (see `Conductor temperature under a fire`_)
  apply; flames that engulf a member heat it far faster than
  :eq:`lt-fire` gives. Near the fire the irradiance depends on the fire's cell
  size, through the sources' positions and the floor :math:`\Delta A/\pi`:
  in the test case, 5 m fire cells give the hottest member 1.35 times the
  irradiance of 10 m cells and heat it to 165 °C instead of 131 °C in 300 s.
* A member that reaches 1200 °C stops the run: EN 1993-1-2's steel keeps no
  strength or stiffness there, and the frame cannot be built.

Verification of the towers' heating
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Test
     - What is pinned
   * - ``TowerThermal.TheSpecificHeatIsEN1993s``,
       ``TowerThermal.TheEnthalpyIsTheSpecificHeatsIntegral``
     - :math:`c_a` at 20, 300, 600, 735 and 900 °C; :math:`H` at 300, 700, 800
       and 1100 °C against an independent numerical integral, to
       :math:`10^{-9}` relative, and across the peak against a Simpson sum.
   * - ``TowerThermal.TheJRCWorkedExampleReachesItsCriticalTemperature``
     - The unprotected beam of the JRC report (Vassart et al., 2014, section
       2.8.3.7): :math:`k_{sh} A_m/V` = 131.3 m\ :sup:`-1` in the ISO 834 fire,
       emissivity 0.7, :math:`\alpha_c` = 25 W/(m\ :sup:`2` K), 3 s steps, reaches
       its critical temperature of 667 °C in about 17 minutes; here in 17 min 6 s.
   * - ``TowerThermal.TheSectionsAreTheAnglesTheTubesAndTheBoxes``
     - :math:`A`, :math:`P_x` and :math:`W` of an angle, solid and hollow tubes
       and a box; an arbitrary section without angle data is refused.
   * - ``TowerThermal.TheSteadyAndImplicitTemperaturesMatchTheIndependentCalculation``
     - A 100 × 8 mm angle: steady temperatures with 1e4 and 5e4 W/m\ :sup:`2`
       (89.399 and 278.387 °C), with an emissivity of 0.35 and 10 W/(m\ :sup:`2` K),
       and in 400 °C air, backward-Euler steps, and a 60 s step from 700 °C
       across the peak, to :math:`10^{-9}` K, each step balancing its energy.
   * - ``TowerThermal.AFrameIsRebuiltWhenItsSteelHasSoftenedAndSurvivesARestart``
     - On a generated tower: no rebuild below 100 °C; a long step under a fire
       brings every member to its steady temperature and rebuilds the frame at
       them, which lowers its frequency; a warming below the tolerance leaves
       the frame but moves the members' checks to the new temperatures; the
       checkpoint restores the temperatures to the bit.
   * - ``Conductors_TowerFireHeating`` and its two-level variants (CTests, stub)
     - A 2e5 W/m\ :sup:`2` disc of 20 m radius under a generated 310-member
       tower for 300 s, passive coupling. An independent calculation
       (``check_tower_heating.py``) of every member from the frame file the run
       wrote and the fire cells, with its own enthalpy, steps and rebuild rule,
       with the emissivity, convection and tolerance away from their defaults:
       every member to :math:`10^{-3}` K at every step, and every row's mean,
       largest, hottest member, largest irradiance, least :math:`k_E`,
       rebuilds and the frame's largest temperature. One level; the tower on a
       level-1 patch; the fire on a level-1 patch.
   * - ``Conductors_TowerFireHeating_Restart`` (CTest, stub)
     - A checkpoint at 270 s, after the first rebuild: the temperature logs and
       the towers' log equal the uninterrupted run's.
   * - ``Conductors_TowerFireHeating_RestartUnheatedAbort``,
       ``Conductors_Thermal_Tower*Abort``, ``Conductors_Thermal_TowersNeed*Abort`` (CTests)
     - A restart without the heating from a hot frame; the heating without the
       thermal model or a frame tower, with a type's ``steel_temperature``, a
       value out of range, and a member past 1200 °C, each stop the run.

References
~~~~~~~~~~

* Vassart, O., Zhao, B., Cajot, L. G., Robert, F., Meyer, U. and Frangi, A.
  (2014). *Eurocodes: Background & Applications. Structural Fire Design. Worked
  examples*. M. Poljanšek, B. Nikolova, L. Sousa, S. Dimova and A. Pinto (eds.),
  JRC Science and Policy Report EUR 26698 EN. Publications Office of the European
  Union, Luxembourg. doi:10.2788/85432. Open access (reproduction authorised
  with acknowledgement):
  https://eurocodes.jrc.ec.europa.eu/sites/default/files/2022-06/2012_11_WS_fire.pdf

Restart
-------

A checkpoint carries the lines under ``<chk>/conductors``: each line's whole
MoorDyn state (node and free-point positions and velocities, internal forces
and the time integrator's state, through MoorDyn's own save), the running
statistics of every span, of each line's strings, of every pair of lines,
transformer and tower, the gusts' per-span means (started afresh when the
checkpoint has none), the random gusts' processes (drawn afresh at the
restart's first step when the checkpoint has none, or other ones: another
number, or lines taking each other's gusts differently), each moving tower's sway, the capacity records (see `Component capacities and first failure`_), the intensity records (see `Wind intensity beside the peak demand`_), the spans' ground-contact records (see `Ground contact`_), the releases of the lines in order and the failed towers (see `Release of a line at a capacity exceedance`_), the conductors' temperatures (see `Conductor temperature under a fire`_), the towers' member temperatures (see `Tower steel temperature under a fire`_), the step count and the time.
On a restart the line is created from the same inputs, initialised without the initial-shape solve and given that state, so
it continues blown out exactly where the checkpoint left it, on MoorDyn's
clock; the statistics go on accumulating, and the diagnostics continue on the
step count of the original run. With ``drag_on_flow`` the restored lines' drag
is spread into the momentum sources again at once, so a plotfile written at
the restart shows the source of the checkpointed step. The span and string
logs, ``<output_root>_nodes.dat``, ``total_load.dat``, ``separation.dat``,
``transformers.dat`` and ``coupling.dat`` are appended to after the rows a run wrote beyond the
checkpoint time are dropped, so a run that went on past its last checkpoint
and is restarted from it leaves no duplicated stretch; ``failures.dat``,
``breaks.dat`` and ``tower_failures.dat`` are written again from the
checkpoint's records; ``towers.dat``,
``gust_series.dat`` and the frame towers' ``tower_<tower>_frame.dat``, whose rows
carry the time a step starts at, lose the row at
the checkpoint time as well, since the restarted run writes it again.

The lines of a restart must be those of the run that wrote the checkpoint.
The run stops, naming the line, when the checkpoint holds a different number
of lines, another line name at a position of ``conductors.lines``, or a
line with another number of nodes (spans, segments, strings). It also stops,
naming the statistics set or the tower, when a pair's, transformer's or
tower's statistics or a moving tower's sway is missing or has another number
of quantities, when the checkpoint holds tower sway and no tower type moves
(a ``frequency`` or a frame), or the other way round, or when ``surface_offset`` differs
from the one the checkpoint records, the frame MoorDyn's saved state is in,
or when the checkpoint holds releases of the lines and ``ground_contact`` is
false, or a failed tower the run does not have.
The line nodes restart where the checkpoint left them, so a line's attachment
points, unstretched lengths and strings' length must be those it was built
with: a restart that moves an end or a tower, or changes a length, a stringing
tension or ``insulator_length``, stops, naming the line (checkpoints written
before this check carry no such record). A frame tower's restart with another
``steel_temperature`` stops too (its frame and state would not match), unless
``conductors.thermal.towers`` heats the towers, whose frames then carry the
fire's temperatures (see `Tower steel temperature under a fire`_). The
checkpoint records each tower's ``diagonal_wind_factor``, ``arm_outside_shaft``
and ``angle_principal_axes`` (on a generated frame), from its type, which set the
towers' drag nodes, drag and frame: a restart that changes one stops, naming the
key, and a checkpoint without the record, written with the earlier forms,
restarts only with ``diagonal_wind_factor`` and ``arm_outside_shaft`` set to
``false``, and ``angle_principal_axes`` too on a generated frame. The
checkpoint also carries the iterations and unconverged steps counted towards
``coupling.dat``'s next row. Other values that leave the counts unchanged (sizes,
stiffnesses, drag coefficients) are not compared: the restarted run uses the
values of its inputs. A checkpoint without conductor state, from a run without conductor
lines, starts the lines afresh from their still-air shape, with MoorDyn's
clock at zero at the restart time.

Checks at start-up
------------------

A run with conductor lines refuses to start with an anchor
level that is not a level of the run, with a ``surface_offset`` that does
not put MoorDyn's free surface above the domain top and its bottom below the
domain bottom, with a span whose length does not exceed the distance between
the points the conductor hangs from placed on the terrain, with a
length missing for a span of a section, with both ``length`` and
``stringing_tension`` on a line, with a ``mass_per_length`` no larger than
the mass of the air the conductor displaces, with insulator strings on a line
without towers or longer than a tower is high, with ``insulator_mass`` or
``insulator_diameter`` on a line without strings, with the attachment points
either side of a string's tower at the same ``(x, y)``, with two lines
writing to the same ``output_root``, with a ``tower_type`` that names no
tower type or is on a line without towers, with a ``share_towers`` line that
has no ``tower_type``, shares another line's towers or has another number of
towers, with a name used twice among the lines, transformers and tower types,
with an attachment or transformer outside the domain, with a shared point off
its tower, with a line that turns back on itself at a tower, with a line end
on two overlapping transformers or not above the top of the transformer it
ends on, with a transformer no line ends on, with transformers or tower
types but no lines, with the key ``conductors.spans`` (the lines are
listed in ``conductors.lines``), with a Real input that is not finite, or
with any other value outside its documented range. Every message names the
input key.

Tests on more than one level
----------------------------

``Conductors_MomentumBudget.i`` is a periodic box of uniform air (density 1,
no gravity, Coriolis force or diffusion, slip walls) blowing at 10 m/s across
one line, made 6 m thick so that its impulse is a thousandth of the air's
momentum. Nothing but the drag changes the air's y-momentum, so

.. math::

   \int_V v \, dV \Big|_{t} - \int_V v \, dV \Big|_{0}
   = \frac{1}{\rho_0} \sum_n S_y^n \, \Delta t ,

where :math:`V` is the box, :math:`v` the y-velocity integrated over the
composite (each place on its finest level), :math:`\rho_0 = 1` kg/m\ :sup:`3`
the density, :math:`S_y^n` the integrated source the model put into the air
in step :math:`n` of its level (``source_y`` in ``conductors/total_load.dat``)
and :math:`\Delta t` that level's step.

.. list-table::
   :header-rows: 1
   :widths: 40 40 20

   * - Test
     - Levels
     - Measured (stub, MoorDyn-C)
   * - ``Conductors_MomentumBudget``
     - one level
     - 5e-7, 3.7e-7
   * - ``Conductors_MomentumBudget_FinestLevel``
     - the lines on level 1, a slab over the whole box that holds them
     - 1.3e-5, 1.0e-3
   * - ``Conductors_MomentumBudget_FinestLevelPatch``
     - the lines on level 1, a patch over the whole width across the line
     - 1.9e-6, 5.3e-4
   * - ``Conductors_MomentumBudget_FinerLevelBelow``
     - the lines on level 0, level 1 below the drag's reach
     - 1.1e-6, 9.8e-5
   * - ``Conductors_Abort_FinerLevelCoversLines``
     - the lines on level 0, level 1 over them
     - the run stops
   * - ``Conductors_Abort_FinerLevelCoversDragTail``
     - the lines on level 0, level 1 over only the last 12.5 m of the drag kernel's reach
     - the run stops
   * - ``Conductors_OneWayUnderFinerLevel``
     - lines without ``drag_on_flow`` on level 0 under level 1
     - the flow of a run without the lines, bit for bit

MoorDyn-C's line swings and bobs where the stub's holds its quasi-static
shape, and the velocities it sends through the coarse/fine faces carry
momentum the levels do not synchronise, so its two-level budgets close less
tightly; a lost drag misses by tens of percent. With ERF's check switched off and level 1 over the whole box, the lines on
level 0 lost all of their drag to the averaging down: the air's momentum did
not change, against an impulse of :math:`-4.5 \times 10^5` N s. The unit test
``Conductors.EachLineAndTowerSpreadsItsDragWithItsOwnWidth`` checks that the
support holds every spread point with its own line's or tower type's reach,
and ``Conductors.WithoutDragOnFlowNothingIsPutIntoTheFlow`` that it is empty
without ``drag_on_flow``.
