.. _sec:ConductorInputs:

Conductor lines (MoorDyn lines)
===============================

Overhead conductors, shield wires and their insulator strings, their motion
in the wind computed by MoorDyn-C; see :ref:`sec:Conductors` for the model
and :doc:`CouplingToMoorDyn` for building FLARE with MoorDyn. The lines are
FLARE's structure model on ERF's structure coupling (:ref:`sec:StructureCoupling`):
build with ``-DERF_ENABLE_STRUCTURES=ON`` (CMake only) and select them with
``erf.structures.model = conductors``; a build without them stops on that
line, naming the models it has. ``amrex.fpe_trap_invalid``,
``amrex.fpe_trap_zero`` and ``amrex.fpe_trap_overflow`` may be on: they are
masked while MoorDyn initialises and steps the lines, whose initial-condition
solver overflows an intermediate value it then discards
(:doc:`CouplingToMoorDyn`), and apply to the rest of the run. Every Real input of a line, a transformer, a
tower type or the shared settings must be finite: a NaN or an infinity stops
the run, naming the key. Lines are named in ``conductors.lines``
(each entry names a whole line, a single span or a section); each has its own
``conductors.<line>.*`` block and is either a single span between two
dead-ends or, with ``towers``, a section of spans over suspension towers.
Transformers, named in ``conductors.transformers`` with blocks
``conductors.<transformer>.*``, are boxes on the terrain the lines end on;
tower types, named in ``conductors.tower_types`` with blocks
``conductors.<type>.*``, make a line's towers lattice towers loaded by the
wind, which bend under their loads when the type has a ``frequency`` or a frame
(``frame_file`` or ``frame_panels``).

+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| Parameter                                                 | Definition                                               | Acceptable Values    | Default                  |
+===========================================================+==========================================================+======================+==========================+
| **conductors.lines**                                      | names of the conductor lines; the feature is off when    | Strings              | none                     |
|                                                           | absent. conductors.spans is not an input: a run that     |                      |                          |
|                                                           | gives it stops, naming conductors.lines                  |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.end_a**                               | first dead-end attachment (m): x, y in ERF coordinates;  | 3 Reals inside the   | must be set              |
|                                                           | the point stands z less geometry.prob_lo[2] above the    | domain               |                          |
|                                                           | terrain surface at (x, y) (with prob_lo[2] = 0, z is its |                      |                          |
|                                                           | height above the terrain; on a flat mesh, its absolute   |                      |                          |
|                                                           | height)                                                  |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.end_b**                               | last dead-end attachment, as end_a                       | 3 Reals inside the   | must be set              |
|                                                           |                                                          | domain               |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.towers**                              | the suspension points between the ends, in order from    | 3 Reals per tower    | none (a single span)     |
|                                                           | end_a, three Reals each as end_a; the line is then a     | inside the domain;   |                          |
|                                                           | section of one span more than the towers                 | with insulator       |                          |
|                                                           |                                                          | strings, the         |                          |
|                                                           |                                                          | attachment points    |                          |
|                                                           |                                                          | either side of a     |                          |
|                                                           |                                                          | tower not at the     |                          |
|                                                           |                                                          | same x, y            |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.length**                              | unstretched length of each span (m), one per span; the   | Reals > the chord of | must be set, or          |
|                                                           | slack above the chord between the points the conductor   | each span between    | stringing_tension        |
|                                                           | hangs from (the attachments, or the bottoms of the       | those points on the  |                          |
|                                                           | insulator strings at the towers), placed on the terrain, | terrain              |                          |
|                                                           | sets the sag                                             |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.stringing_tension**                   | instead of length: the horizontal tension every span     | Real > 0             | none                     |
|                                                           | carries in still air (N); each span's length is set from |                      |                          |
|                                                           | its chord on the terrain, between the bottoms of the     |                      |                          |
|                                                           | strings at the towers                                    |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.diameter**                            | conductor diameter (m)                                   | Real > 0             | must be set              |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.mass_per_length**                     | mass per unstretched metre (kg/m), as MoorDyn takes it   | Real > air_density * | must be set              |
|                                                           |                                                          | pi * diameter^2 / 4  |                          |
|                                                           |                                                          | (the air the         |                          |
|                                                           |                                                          | conductor displaces) |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.axial_stiffness**                     | axial stiffness EA (N)                                   | Real > 0             | must be set              |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.drag_coefficient**                    | drag coefficient normal to the line, on the diameter (-) | Real >= 0            | 1.0                      |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.damping_ratio**                       | internal damping of the line (and its insulator strings) | Real in (0, 1]       | 0.5                      |
|                                                           | as a fraction of critical (-)                            |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.segments**                            | segments per span (count; nodes per span = segments + 1) | Integer >= 2         | 20                       |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.insulator_length**                    | length of the suspension insulator string the conductor  | Real >= 0, its       | 0                        |
|                                                           | hangs from at each tower (m); 0: the conductor is        | bottom above the     |                          |
|                                                           | clamped to the towers                                    | ground under each    |                          |
|                                                           |                                                          | tower; needs towers  |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.insulator_mass**                      | mass of one insulator string (kg)                        | Real > 0; refused    | must be set with         |
|                                                           |                                                          | when                 | insulator_length         |
|                                                           |                                                          | insulator_length is  |                          |
|                                                           |                                                          | 0                    |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.insulator_diameter**                  | disc diameter of the insulator strings, which the wind   | Real > 0; refused    | 0.254                    |
|                                                           | acts on (m)                                              | when                 |                          |
|                                                           |                                                          | insulator_length is  |                          |
|                                                           |                                                          | 0                    |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.rated_breaking_strength**             | the conductor's rated breaking strength (N), checked     | Real > 0             | none: not checked        |
|                                                           | against the largest tension of every span after every    |                      |                          |
|                                                           | coupling step (capacities.csv, failures.dat)             |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.insulator_rated_load**                | the rated mechanical load (specified mechanical load, or | Real > 0; needs      | none: not checked        |
|                                                           | the M&E rating of one unit) of each insulator string     | towers and           |                          |
|                                                           | (N), checked against the tension at the string's top     | insulator_length > 0 |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.hardware_rated_load**                 | the rated load of the attachment hardware (N), checked   | Real > 0             | none: not checked        |
|                                                           | against the magnitude of the line's pull at each dead    |                      |                          |
|                                                           | end and each tower                                       |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.absorptivity**                        | the conductor surface's absorptivity of the sun's and    | Real in [0, 1];      | none                     |
|                                                           | the fire's radiation (IEEE 738)                          | needed with          |                          |
|                                                           |                                                          | thermal_model =      |                          |
|                                                           |                                                          | ieee738              |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.emissivity**                          | the conductor surface's emissivity (IEEE 738)            | Real in [0, 1];      | none                     |
|                                                           |                                                          | needed with          |                          |
|                                                           |                                                          | thermal_model =      |                          |
|                                                           |                                                          | ieee738              |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.heat_capacity**                       | the conductor's heat capacity per unit length, m c_p,    | Real > 0; needed     | none                     |
|                                                           | summed over its materials (J/(m K))                      | with thermal_model = |                          |
|                                                           |                                                          | ieee738              |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.thermal_expansion**                   | the conductor's thermal expansion coefficient (1/K);     | Real in [0, 1e-4];   | none                     |
|                                                           | each span's unstretched length follows its mean          | needed with          |                          |
|                                                           | temperature                                              | thermal_model =      |                          |
|                                                           |                                                          | ieee738              |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.current**                             | the current the conductor carries (A), for its Joule     | Real in [0, 1e5]     | 0                        |
|                                                           | heating                                                  |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.resistance**                          | the conductor's AC resistance at the two                 | two Reals > 0,       | none                     |
|                                                           | resistance_temperatures (ohm/m); linear between and      | staying positive     |                          |
|                                                           | beyond them; not used without a current                  | from -150 to 3000 C; |                          |
|                                                           |                                                          | needed with a        |                          |
|                                                           |                                                          | current              |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.resistance_temperatures**             | the temperatures of the two resistance values (C); for   | two increasing Reals | 25 75                    |
|                                                           | conductors heated past 100 C give 25 200 (IEEE Std       | in [-50, 500]        |                          |
|                                                           | 738-2023, 4.4.6)                                         |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.tower_type**                          | the suspension towers of the line are lattice towers of  | one of the tower     | none                     |
|                                                           | this type, loaded by the wind on their members and by    | types; needs towers  |                          |
|                                                           | the line, on four footings (towers.dat); none: the       |                      |                          |
|                                                           | towers are points                                        |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.share_towers**                        | the line hangs from the towers of this line (which has   | a line of            | none: its own towers     |
|                                                           | the tower_type), each at its own towers point on them: a | conductors.lines     |                          |
|                                                           | circuit's phases across the cross-arm and its shield     | with a tower_type    |                          |
|                                                           | wire on the peak; the points stand above the tower's     | and as many towers   |                          |
|                                                           | base                                                     | that shares no other |                          |
|                                                           |                                                          | line's towers; no    |                          |
|                                                           |                                                          | tower_type of its    |                          |
|                                                           |                                                          | own                  |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.gust_with**                           | with gust_type = random the line takes this line's       | a line of            | none: its own, or those  |
|                                                           | random gusts, span by span, so that the conductors of a  | conductors.lines     | of the line whose towers |
|                                                           | circuit on separate lines move together                  | with as many spans   | it shares                |
|                                                           |                                                          | that takes no other  |                          |
|                                                           |                                                          | line's gusts; needs  |                          |
|                                                           |                                                          | gust_type = random;  |                          |
|                                                           |                                                          | not with             |                          |
|                                                           |                                                          | share_towers         |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.epsilon**                             | width of the Gaussian that spreads the drag of this      | Real > 0             | conductors.epsilon       |
|                                                           | line's nodes, its insulator strings' included, in cells  |                      |                          |
|                                                           | (dx); the towers it hangs from take their type's; given  |                      |                          |
|                                                           | without drag_on_flow it stops the run, and with          |                      |                          |
|                                                           | drag_on_flow = false given it is ignored with a warning  |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<line>.output_root**                         | prefix of the line's diagnostics files:                  | String, different    | <diagnostics_dir>/<line> |
|                                                           | <output_root>.dat for a single span,                     | for every line; no   |                          |
|                                                           | <output_root>_span<k>.dat per span of a section (k from  | log of the line may  |                          |
|                                                           | 1 at end_a), <output_root>_insulators.dat,               | take the name of a   |                          |
|                                                           | <output_root>_nodes.dat and the matching _stats.csv      | log the run writes   |                          |
|                                                           | files                                                    | in diagnostics_dir,  |                          |
|                                                           |                                                          | nor its statistics   |                          |
|                                                           |                                                          | the name or file of  |                          |
|                                                           |                                                          | the run's own        |                          |
|                                                           |                                                          | (towers,             |                          |
|                                                           |                                                          | transformers, pairs  |                          |
|                                                           |                                                          | of lines)            |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.diagnostics_dir**                            | directory of the MoorDyn input files                     | String               | conductors               |
|                                                           | (<line>.moordyn.txt), ground.dat, total_load.dat,        |                      |                          |
|                                                           | separation.dat, transformers.dat, towers.dat,            |                      |                          |
|                                                           | coupling.dat, asce74.csv, gusts.csv, gust_series.dat,    |                      |                          |
|                                                           | the generated frames (frame_<tower>.dat), the frame      |                      |                          |
|                                                           | towers' logs (tower_<tower>_frame.dat,                   |                      |                          |
|                                                           | tower_<tower>_members.csv), the statistics of pairs,     |                      |                          |
|                                                           | transformers and towers, and the lines' default          |                      |                          |
|                                                           | output_root                                              |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.diagnostics_int**                            | write a diagnostics row every this many steps of the     | Integer >= 1         | 1                        |
|                                                           | anchor level (steps), and at the first step              |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.anchor_level**                               | the level the lines are stepped with (-1:                | -1 .. amr.max_level  | -1                       |
|                                                           | amr.max_level); it must exist from the start and for the | (a level that exists)|                          |
|                                                           | whole run, or the run stops                              |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.air_density**                                | the fluid density MoorDyn applies the drag with (kg/m^3) | Real > 0             | 1.225                    |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.substeps**                                   | MoorDyn calls per ERF step for lines on fixed supports,  | Integer >= 1         | 1                        |
|                                                           | the wind held over the step; lines on moving towers take |                      |                          |
|                                                           | at least this many coupling steps (count)                |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.moordyn_dt**                                 | an upper bound on MoorDyn's internal step (s); 0: no     | Real >= 0            | 0                        |
|                                                           | bound beyond moordyn_cfl                                 |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.moordyn_cfl**                                | the Courant factor (-) that sets MoorDyn's internal      | Real in (0, 0.15]    | 0.1                      |
|                                                           | step; larger values make the line integration diverge or |                      |                          |
|                                                           | give a wrong sag                                         |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.moordyn_log_level**                          | MoorDyn's verbosity: 0 debug, 1 messages, 2 warnings, 3  | 0 .. 3               | 2                        |
|                                                           | errors only                                              |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.surface_offset**                             | MoorDyn's free surface sits this far above ERF's z = 0   | Real >               | 10000                    |
|                                                           | and its flat bottom this far below (m): the free surface | geometry.prob_hi[2], |                          |
|                                                           | must be above the domain top and the bottom below the    | with -surface_offset |                          |
|                                                           | domain bottom                                            | <                    |                          |
|                                                           |                                                          | geometry.prob_lo[2]  |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.prescribed_velocity**                        | a uniform wind at every line node and every tower drag   | 3 Reals              | none (the flow is        |
|                                                           | node (m/s) instead of the flow's, for testing            |                      | sampled at the nodes)    |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.stats_start**                                | time (s) from which the running statistics start: of     | Real >= 0            | 0                        |
|                                                           | every span (<output_root>_stats.csv, or                  |                      |                          |
|                                                           | <output_root>_span<k>_stats.csv in a section: mean, rms, |                      |                          |
|                                                           | minimum, maximum of the swing, mid-span offset, end and  |                      |                          |
|                                                           | maximum tensions, minimum clearance and the y component  |                      |                          |
|                                                           | of the air's drag), of the insulator strings, of the     |                      |                          |
|                                                           | closest approach of every pair of lines, of the          |                      |                          |
|                                                           | transformers' loads and clearances and of the towers'    |                      |                          |
|                                                           | loads. The conductors' times, here, in gust_event_time   |                      |                          |
|                                                           | and in their logs, are ERF's elapsed time since the run  |                      |                          |
|                                                           | began, which a restart continues (erf.start_time is not  |                      |                          |
|                                                           | added)                                                   |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.node_output_int**                            | write every node of every line (position, clearance,     | Integer >= 0         | 0                        |
|                                                           | tension, wind, drag) to <output_root>_nodes.dat, and     |                      |                          |
|                                                           | every frame tower's node displacements and member        |                      |                          |
|                                                           | utilisations to diagnostics_dir/tower_<tower>_frame.dat, |                      |                          |
|                                                           | every this many steps (steps); 0: never                  |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.drag_on_flow**                               | put the air's drag on the lines and the towers,          | Boolean; no finer    | false                    |
|                                                           | reversed, back into the flow as a momentum source on the | level may cover the  |                          |
|                                                           | anchor level (and allow the plot variables conductor_fx, | drag's support (the  |                          |
|                                                           | conductor_fy, conductor_fz)                              | run stops)           |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.epsilon**                                    | width of the Gaussian that spreads the drag onto the     | Real > 0             | 2.0                      |
|                                                           | momentum sources, in cells (dx), for every line and      |                      |                          |
|                                                           | tower type that sets no epsilon of its own (a warning    |                      |                          |
|                                                           | when every line and the type of every line's towers set  |                      |                          |
|                                                           | their own); given without drag_on_flow it stops the run, |                      |                          |
|                                                           | and with drag_on_flow = false given it is ignored with a |                      |                          |
|                                                           | warning                                                  |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.flashover_distance**                         | two lines whose conductors come closer than this are     | Real > 0             | 1.0                      |
|                                                           | flagged as clashing in separation.dat, and a conductor   |                      |                          |
|                                                           | this close to a transformer in transformers.dat (m)      |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.failure_stop**                               | stop the run (an abort naming the component) at the      | true or false; true  | false                    |
|                                                           | first coupling step where a component's demand exceeds   | needs a capacity: a  |                          |
|                                                           | its capacity, after failures.dat and capacities.csv are  | line's rated_breakin |                          |
|                                                           | written; exceedances restored from a checkpoint do not   | g_strength,          |                          |
|                                                           | stop it                                                  | insulator_rated_load |                          |
|                                                           |                                                          | or                   |                          |
|                                                           |                                                          | hardware_rated_load, |                          |
|                                                           |                                                          | a tower type's       |                          |
|                                                           |                                                          | allowables, or a     |                          |
|                                                           |                                                          | frame's member       |                          |
|                                                           |                                                          | design data          |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.thermal_model**                              | the conductors' temperature: ieee738 steps every span    | none or ieee738      | none                     |
|                                                           | node's IEEE 738 heat balance with the fire's radiation   |                      |                          |
|                                                           | and the air around it, and lengthens each span with its  |                      |                          |
|                                                           | mean temperature (conductor_temperature.dat; theory:     |                      |                          |
|                                                           | Conductor temperature under a fire)                      |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.thermal.radiative_fraction**                 | the fraction of the fire's heat release it radiates,     | Real in (0, 1];      | none                     |
|                                                           | chi_r                                                    | needed when the run  |                          |
|                                                           |                                                          | has a fire           |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.thermal.source_height**                      | the height above the terrain of the point sources that   | Real in [0, 1000]    | 0                        |
|                                                           | stand for the burning fire cells (m)                     |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.thermal.cutoff_radius**                      | sources farther than this from a node are left out of    | Real in [0, 1e6]; >  | 0                        |
|                                                           | its irradiance (m); only burning cells within it of the  | 0 needed when the    |                          |
|                                                           | lines are gathered every step                            | run has a fire       |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.thermal.elevation**                          | the lines' elevation above sea level (m), for the air's  | Real in [-500, 9000] | 0                        |
|                                                           | density                                                  |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.thermal.solar_flux**                         | the sun's flux on the conductors, normal to them (W/m^2) | Real in [0, 1500]    | 0                        |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.thermal.reference_temperature**              | the temperature at which the spans have their input      | Real in [-50, 300]   | none: each span's mean   |
|                                                           | unstretched lengths (C), a stringing temperature say     |                      | steady temperature at    |
|                                                           |                                                          |                      | the first step (it       |
|                                                           |                                                          |                      | starts at its input      |
|                                                           |                                                          |                      | length)                  |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.thermal.towers**                             | the towers' steel temperature: every member of every     | Boolean; needs       | false                    |
|                                                           | tower that bends as a frame follows EN 1993-1-2's heat   | thermal_model =      |                          |
|                                                           | balance of unprotected steel under the fire, and the     | ieee738, a frame     |                          |
|                                                           | frame is rebuilt as the steel softens                    | tower and no         |                          |
|                                                           | (tower_temperature.dat; theory: Tower steel temperature  | steel_temperature    |                          |
|                                                           |                                                          | other than 20        |                          |
|                                                           | under a fire)                                            |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.thermal.tower_emissivity**                   | the steel's emissivity, also its absorptivity of the     | Real in (0, 1]       | 0.7 (EN 1993-1-2 2.2)    |
|                                                           | fire's radiation                                         |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.thermal.tower_convection**                   | the convection coefficient between the members and the   | Real in [0, 200]     | 25 (EN 1991-1-2 3.2.1)   |
|                                                           | air (W/(m^2 K))                                          |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.thermal.tower_reduction_tolerance**          | rebuild a tower's frame when any member's k_y or k_E     | Real in (0, 0.2]     | 0.01                     |
|                                                           | (EN 1993-1-2 Table 3.1) has moved by more than this      |                      |                          |
|                                                           | fraction of the larger of its two values from the        |                      |                          |
|                                                           | frame's                                                  |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.ground_contact**                             | the terrain holds the lines up: each line's MoorDyn      | true or false;       | false                    |
|                                                           | system takes the anchor level's terrain under it, raised | MoorDyn's internal   |                          |
|                                                           | by the conductor's radius, as its seafloor               | step must pass the   |                          |
|                                                           | (<line>.seafloor.txt), with critical contact damping;    | contact's stability  |                          |
|                                                           | each span's first touch of the ground and the most nodes | check                |                          |
|                                                           | on it go to ground_contact.csv                           |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.asce74_wind**                                | the 3-second gust at 10 m over open country of an ASCE   | Real >= 0            | 0: no check              |
|                                                           | 74 design check (m/s): every span's quasi-static wind    |                      |                          |
|                                                           | load, swing, blowout and tension under it are written to |                      |                          |
|                                                           | diagnostics_dir/asce74.csv at start-up (see the ASCE 74  |                      |                          |
|                                                           | design check of the conductor theory)                    |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.asce74_exposure**                            | ASCE 74's terrain exposure of that check: B (suburban or | B or C; needs        | C                        |
|                                                           | wooded) or C (open country)                              | asce74_wind          |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.asce74_wire_height**                         | the wire height of that check: effective, the mean of    | effective or         | effective                |
|                                                           | the span's attachment heights above the ground less a    | attachment; needs    |                          |
|                                                           | third of its still-air sag (vertical), or attachment,    | asce74_wind          |                          |
|                                                           | the mean alone                                           |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.asce74_inclined_spans**                      | the check's swing about an inclined chord against the    | Boolean; needs       | true                     |
|                                                           | weight's part normal to it, W cos(beta), and the upper   | asce74_wind          |                          |
|                                                           | end's tension (false: the level span's swing, sag and    |                      |                          |
|                                                           | tension under the whole weight W; the effective height's |                      |                          |
|                                                           | still-air sag is the inclined chord's either way)        |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.gust_type**                                  | none; factor: every span's gust factor and peak wind     | none, factor, event  | none                     |
|                                                           | load from the mean wind and the k-equation RANS's k      | or random; not with  |                          |
|                                                           | along it, written to diagnostics_dir/gusts.csv; event:   | prescribed_velocity; |                          |
|                                                           | also one travelling 1 - cos gust, random: also a random  | every span needs a   |                          |
|                                                           | gust per span and tower, added to the wind the lines and | horizontal extent;   |                          |
|                                                           | towers take and written to                               | event and random not |                          |
|                                                           | diagnostics_dir/gust_series.dat (see Gusts from the RANS | with drag_on_flow    |                          |
|                                                           | turbulence in the conductor theory). Every type but none |                      |                          |
|                                                           | needs erf.rans_type = kEqn on the anchor level and the   |                      |                          |
|                                                           | flow's wind                                              |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.gust_sigma_factor**                          | c in sigma_u = c sqrt(k)                                 | Real > 0; needs a    | 2.5 erf.Cmu0 (1.39)      |
|                                                           |                                                          | gust_type            |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.gust_peak_factor**                           | g in the gust factor 1 + 2 g I_n sqrt(B), the point gust | Real > 0; needs a    | 2.7 (ASCE 74)            |
|                                                           | U + g sigma_u and the event's amplitude g sigma_u        | gust_type            |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.gust_span_length_scale**                     | L_s in the spans' B = 1/(1 + 0.8 L/L_s), L the span's    | Real > 0; needs a    | 67.056 (220 ft)          |
|                                                           | chord, and the towers' B_t = 1/(1 + 0.375 h/L_s) (m);    | gust_type            |                          |
|                                                           | the default is ASCE 74's for asce74_exposure when it is  |                      |                          |
|                                                           | given, else exposure C's (51.816 for B)                  |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.gust_event_time**                            | when the event's front crosses gust_event_origin (s,     | Real; needs and is   | must be set              |
|                                                           | ERF's time since the run began, as stats_start)          | needed by gust_type  |                          |
|                                                           |                                                          | = event              |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.gust_event_speed**                           | how fast the event's front moves (m/s)                   | Real > 0; needs and  | must be set              |
|                                                           |                                                          | is needed by         |                          |
|                                                           |                                                          | gust_type = event    |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.gust_event_direction**                       | the way the front moves, degrees counterclockwise from   | Real; needs          | 0                        |
|                                                           | +x                                                       | gust_type = event    |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.gust_event_origin**                          | a point x y the front crosses at gust_event_time (m)     | 2 Reals; needs       | the domain's centre      |
|                                                           |                                                          | gust_type = event    |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.gust_event_duration**                        | how long the event lasts at a point (s); the default is  | Real > 0; needs      | 10.5                     |
|                                                           | the duration of IEC 61400-1's extreme operating gust,    | gust_type = event    |                          |
|                                                           | not its shape or amplitude                               |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.gust_seed**                                  | picks the random gusts' series                           | Integer >= 0; needs  | 1                        |
|                                                           |                                                          | gust_type = random   |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.gust_integral_length**                       | the random gusts' integral length L_u (m); their time    | Real > 0; needs      | 8.1 Lambda_1 of IEC      |
|                                                           | scale is L_u / U                                         | gust_type = random   | 61400-1 at each span's   |
|                                                           |                                                          |                      | and tower's height       |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.tower_types**                                | names of the lattice tower types the lines' towers can   | Strings (not the     | none                     |
|                                                           | be, each with its own block; needs conductors.lines      | name of a line or    |                          |
|                                                           |                                                          | transformer)         |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.base_width**                          | the tower body's width at the ground (m)                 | Real > 0             | must be set              |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.top_width**                           | its width at the cross-arm (m), the body tapering        | Real in (0,          | must be set              |
|                                                           | linearly between                                         | base_width]          |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.solidity**                            | the members' area over the outline of a face, of the     | Real in (0, 1)       | must be set              |
|                                                           | body and the cross-arm (-)                               |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.arm_length**                          | the cross-arm's length across the line, centred on the   | Real > 0; with       | must be set              |
|                                                           | body (m)                                                 | arm_outside_shaft >  |                          |
|                                                           |                                                          | the body's width at  |                          |
|                                                           |                                                          | the middle of the    |                          |
|                                                           |                                                          | arm's face (at least |                          |
|                                                           |                                                          | top_width)           |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.arm_depth**                           | the depth of the cross-arm's face at the body (m); it    | Real >= 0 (0:        | top_width                |
|                                                           | tapers linearly to 0 at the tips                         | top_width); with     |                          |
|                                                           |                                                          | arm_outside_shaft    |                          |
|                                                           |                                                          | under twice the      |                          |
|                                                           |                                                          | cross-arm's height   |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.peak**                                | the body above the cross-arm, at the top width (m)       | Real >= 0            | 0                        |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.drag_coefficient**                    | the force coefficient on the members' area (-); 0: 4     | Real >= 0            | 0                        |
|                                                           | phi^2 - 5.9 phi + 4 from the solidity phi (square        |                      |                          |
|                                                           | lattice tower)                                           |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.segments**                            | drag nodes along the body up to the cross-arm (count)    | Integer >= 1         | 10                       |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.diagonal_wind_factor**                | the body's and the peak's drag times 1 + g sin^2(2       | Boolean              | true                     |
|                                                           | theta), theta the wind's angle from a face's normal and  |                      |                          |
|                                                           | g = min(0.75 solidity, 0.2): ASCE 7's factor along a     |                      |                          |
|                                                           | diagonal (false: none)                                   |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.arm_outside_shaft**                   | the cross-arm's drag nodes on its two parts outside the  | Boolean; true needs  | true                     |
|                                                           | body, from the body's half width at the middle of the    | arm_length > the     |                          |
|                                                           | face's depth at the body (arm_depth, or top_width when   | body's width there   |                          |
|                                                           | that is 0, below the cross-arm) to the tips, each at the | (> top_width)        |                          |
|                                                           | centroid of its piece of the tapering face (false, the   |                      |                          |
|                                                           | earlier layout: along its whole length at the cross-arm's|                      |                          |
|                                                           | height, the face's depth all along)                      |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.weight**                              | the tower's own weight on its foundation (N)             | Real >= 0            | 0                        |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.leg_spacing**                         | the spacing of the four legs' footings, a square (m)     | Real >= 0 (0:        | base_width               |
|                                                           |                                                          | base_width); not     |                          |
|                                                           |                                                          | with frame_file or   |                          |
|                                                           |                                                          | frame_panels         |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.allowable_uplift**                    | a footing's pull-out capacity (N): flagged in towers.dat | Real >= 0 (0: not    | 0                        |
|                                                           | above it, and the capacity of each tower's footing       | checked)             |                          |
|                                                           | uplift after every coupling step (capacities.csv,        |                      |                          |
|                                                           | failures.dat; a capacity for failure_stop)               |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.allowable_compression**               | a footing's bearing capacity (N): flagged above it, and  | Real >= 0 (0: not    | 0                        |
|                                                           | the capacity of each tower's footing compression after   | checked)             |                          |
|                                                           | every coupling step (capacities.csv, failures.dat; a     |                      |                          |
|                                                           | capacity for failure_stop)                               |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.frequency**                           | the tower's first bending frequency on a rigid           | Real >= 0; > 0 needs | 0                        |
|                                                           | foundation (Hz): the tower sways in that mode under its  | weight > 0; not with |                          |
|                                                           | loads and MoorDyn moves its cross-arm as a coupled       | frame_file or        |                          |
|                                                           | point; 0 without a frame: the tower stands still and its | frame_panels         |                          |
|                                                           | attachments are fixed MoorDyn points                     |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.damping_ratio**                       | the structural damping ratio (-, fraction of critical):  | Real in [0, 1); > 0  | 0.02                     |
|                                                           | of the one mode with frequency; with frame_file or       | with frame_file or   |                          |
|                                                           | frame_panels, the frame's Rayleigh damping at its first  | frame_panels; needs  |                          |
|                                                           | natural frequency and at ten times it. On a tower that   | frequency,           |                          |
|                                                           | stands still it stops the run, or with frequency = 0     | frame_file or        |                          |
|                                                           | given is ignored with a warning                          | frame_panels         |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.foundation_rotational_stiffness**     | the footing's resistance to tilting (N m/rad); on a      | Real >= 0 (0:        | 0                        |
|                                                           | tower that stands still as damping_ratio                 | rigid); > 0 needs    |                          |
|                                                           |                                                          | frequency            |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.foundation_lateral_stiffness**        | the footing's resistance to sliding (N/m); on a tower    | Real >= 0 (0:        | 0                        |
|                                                           | that stands still as damping_ratio                       | rigid); > 0 needs    |                          |
|                                                           |                                                          | frequency            |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.frame_file**                          | a SubDyn input file of the tower's frame in tower-local  | File name; not with  | none                     |
|                                                           | axes (origin at the base centre, x along the line, y     | frame_panels,        |                          |
|                                                           | along the cross-arm, z up; see the lattice tower frame   | frequency, weight or |                          |
|                                                           | model): the tower bends as that frame, which gives its   | the foundation       |                          |
|                                                           | stiffness, mass and footing loads; MoorDyn moves its     | stiffnesses          |                          |
|                                                           | cross-arm as a coupled point. Its four supports stand at |                      |                          |
|                                                           | z = 0, one per quadrant, its one interface joint is the  |                      |                          |
|                                                           | cross-arm's centre at the tower's cross-arm height       |                      |                          |
|                                                           | (within half the arm's face depth), every drag node and  |                      |                          |
|                                                           | line attachment of the tower lies within base_width of a |                      |                          |
|                                                           | frame node, and it has at most 6000 free degrees of      |                      |                          |
|                                                           | freedom                                                  |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.member_file**                         | the design data of every member of frame_file's frame    | File name; needs     | none: no checks          |
|                                                           | (role, yield strength, angle, net area, end conditions;  | frame_file           |                          |
|                                                           | see the lattice tower frame model): each member is       |                      |                          |
|                                                           | checked against its strength (ASCE 10-15) every step     |                      |                          |
|                                                           | and, as a capacity, after every coupling step            |                      |                          |
|                                                           | (capacities.csv, failures.dat; a capacity for            |                      |                          |
|                                                           | failure_stop)                                            |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.frame_panels**                        | generate the tower's frame from this type's dimensions,  | Integer in [0, 200], | 0: none                  |
|                                                           | with this many panels of the shaft below the cross-arm:  | and a frame of at    |                          |
|                                                           | a lattice of equal-leg angles (see the lattice tower     | most 6000 free       |                          |
|                                                           | frame model), one per type and cross-arm height, written | degrees of freedom   |                          |
|                                                           | to diagnostics_dir/frame_<tower>.dat as a SubDyn file;   | (about 120 panels    |                          |
|                                                           | its members are checked                                  | with crossed         |                          |
|                                                           |                                                          | bracing); not with   |                          |
|                                                           |                                                          | frame_file,          |                          |
|                                                           |                                                          | frequency, weight or |                          |
|                                                           |                                                          | the foundation       |                          |
|                                                           |                                                          | stiffnesses          |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.leg_angle**                           | the leg width and thickness (m) of the equal-leg angle   | 2 Reals, 0 <         | must be set with         |
|                                                           | of the generated frame's legs and cross-arm chords       | thickness < width;   | frame_panels             |
|                                                           |                                                          | needs frame_panels   |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.brace_angle**                         | the same for every other member of the generated frame   | 2 Reals, 0 <         | must be set with         |
|                                                           |                                                          | thickness < width;   | frame_panels             |
|                                                           |                                                          | needs frame_panels   |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.bracing**                             | the generated frame's face bracing: crossed (two         | crossed or single;   | crossed                  |
|                                                           | diagonals per panel face, bolted where they cross) or    | needs frame_panels   |                          |
|                                                           | single (one diagonal, alternating from panel to panel)   |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.yield_strength**                      | the yield strength of the generated frame's steel (Pa)   | Real >= 0 (0:        | 3.45e8                   |
|                                                           |                                                          | 3.45e8); needs       |                          |
|                                                           |                                                          | frame_panels         |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.angle_principal_axes**                | turn each generated angle (MSpin) so that its axis of    | Boolean; only with   | true                     |
|                                                           | symmetry, its major principal axis, lies as on a tower:  | frame_panels         |                          |
|                                                           | a leg's or a cross-arm chord's towards its corner, a     |                      |                          |
|                                                           | face member's at 45 degrees to its face (false: every    |                      |                          |
|                                                           | MSpin 0)                                                 |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.steel_temperature**                   | the steel's temperature (C): it lowers the frame's       | Real in (-273.15,    | 20                       |
|                                                           | stiffness and its members' yield strength by EN 1993-1-2 | 1200); not 20 needs  |                          |
|                                                           | Table 3.1                                                | frame_file or        |                          |
|                                                           |                                                          | frame_panels; a      |                          |
|                                                           |                                                          | restart cannot       |                          |
|                                                           |                                                          | change it; not with  |                          |
|                                                           |                                                          | thermal.towers       |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<type>.epsilon**                             | width of the Gaussian that spreads the drag of this      | Real > 0             | conductors.epsilon       |
|                                                           | type's towers' drag nodes, in cells (dx); given without  |                      |                          |
|                                                           | drag_on_flow it stops the run, and with drag_on_flow =   |                      |                          |
|                                                           | false given it is ignored with a warning                 |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.transformers**                               | names of the transformers the lines end on: every line   | Strings, each the    | none                     |
|                                                           | end whose x, y lies on a transformer's footprint is      | name of no line or   |                          |
|                                                           | dead-ended on it; needs conductors.lines                 | tower type; at least |                          |
|                                                           |                                                          | one line end on each |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<transformer>.position**                     | the centre of the transformer's footprint (m); the box   | 2 Reals inside the   | must be set              |
|                                                           | stands on the terrain under it                           | domain               |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<transformer>.size**                         | the box's length along x, width along y and height (m);  | 3 Reals > 0          | must be set              |
|                                                           | the line ends on it must be above its top                |                      |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<transformer>.allowable_force**              | the horizontal pull of the lines the transformer may     | Real >= 0 (0: not    | 0                        |
|                                                           | take (N); flagged in transformers.dat above it           | checked)             |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
| **conductors.<transformer>.allowable_moment**             | the overturning moment about the centre of its base it   | Real >= 0 (0: not    | 0                        |
|                                                           | may take (N m); flagged above it                         | checked)             |                          |
+-----------------------------------------------------------+----------------------------------------------------------+----------------------+--------------------------+
