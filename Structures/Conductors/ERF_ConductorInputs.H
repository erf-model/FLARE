#ifndef ERF_CONDUCTOR_INPUTS_H
#define ERF_CONDUCTOR_INPUTS_H

/**
 * Inputs of the conductor module (conductors.*): one LineInputs per conductor line, one
 * TransformerInputs per transformer, the lattice tower types, and the shared settings with their
 * range checks; also the still-air elastic catenary used to check MoorDyn's shape, and a helper
 * that finds non-finite values at the module's boundaries.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include <AMReX_REAL.H>

#include "ERF_Gusts.H"
#include "ERF_TowerInputs.H"

namespace erf_conductors {

/**
 * One conductor line from the conductors.<name>.* input block: one conductor (a phase or a
 * shield wire) dead-ended at end_a and end_b, stepped as one MoorDyn system. Without towers the line
 * is a single span between the two ends; with towers it is a section of several spans over
 * suspension towers, the conductor hanging at each tower from an insulator string (or clamped to
 * the tower when insulator_length is 0), so that the spans on either side of a tower balance their
 * pull along the line as on a real line. Keys marked "required" have no default.
 *
 *   conductors.lines            = S1
 *   conductors.S1.end_a         = 100. 500. 30.    (required; m; x, y in ERF coordinates, inside the domain;
 *                                                       z - geometry.prob_lo[2] is the height above the terrain
 *                                                       surface at x, y, so z is that height when prob_lo z is 0)
 *   conductors.S1.end_b         = 1000. 500. 30.   (required; m; as end_a)
 *   conductors.S1.towers        = 400. 500. 30.  700. 500. 30.   (m, the suspension points between the ends, in
 *                                                       order, each as end_a; default none: one span)
 *   conductors.S1.length        = 301.5 301.5 301.5   (m, unstretched, one per span; each must exceed the
 *                                                       distance between the points the conductor hangs from where
 *                                                       they stand on the terrain; required unless stringing_tension)
 *   conductors.S1.stringing_tension = 2.0e4        (N, > 0, instead of length: every span strung to this
 *                                                       horizontal tension in still air; default 0: none)
 *   conductors.S1.diameter      = 0.0281           (required; m, > 0)
 *   conductors.S1.mass_per_length = 1.628          (required; kg/m, more than the air the conductor displaces,
 *                                                       air_density * pi * diameter^2 / 4)
 *   conductors.S1.axial_stiffness = 3.0e7          (required; N, > 0; EA, Young's modulus E times the conductor's
 *                                                       cross-section area A)
 *   conductors.S1.drag_coefficient = 1.0           (>= 0, the normal drag coefficient on the diameter; default 1)
 *   conductors.S1.damping_ratio = 0.5              (in (0, 1], the internal damping of the conductor and of its
 *                                                       insulator strings as a fraction of critical; default 0.5)
 *   conductors.S1.segments      = 20               (>= 2, line segments per span; default 20)
 *   conductors.S1.insulator_length = 2.5           (m, >= 0 and below every tower's z: the suspension strings at
 *                                                       the towers; needs towers; default 0: clamped)
 *   conductors.S1.insulator_mass = 60.             (kg per string, > 0; required with insulator_length, refused
 *                                                       without it)
 *   conductors.S1.insulator_diameter = 0.254       (m, > 0, the disc diameter the wind acts on; default 0.254;
 *                                                       refused without insulator_length)
 *   conductors.S1.rated_breaking_strength = 1.401e5 (N, > 0; the conductor's rated breaking strength (RBS): every
 *                                                       span's largest tension is checked against it after every
 *                                                       coupling step (ERF_Capacities.H); default none: not checked)
 *   conductors.S1.insulator_rated_load = 1.11e5    (N, > 0; the rated mechanical load of each insulator string (its
 *                                                       specified mechanical load, SML), checked against the tension at
 *                                                       the string's top; needs insulator_length; default none)
 *   conductors.S1.hardware_rated_load = 1.5e5      (N, > 0; the rated load of the attachment hardware (dead-end
 *                                                       fittings, suspension clamps or shackles), checked against the
 *                                                       line's pull at each dead end and at each tower; default none)
 *   conductors.S1.output_root   = conductors/S1    (prefix of the line's diagnostics files, one per line;
 *                                                       default <diagnostics_dir>/<name>; no log of the line may
 *                                                       take the name of a log the run writes there itself, such
 *                                                       as total_load.dat, or towers.dat with lattice towers)
 *   conductors.S1.tower_type    = lattice          (one of conductors.tower_types: the suspension towers are
 *                                                       lattice towers of this type, loaded by the wind; needs towers;
 *                                                       default none: points only)
 *   conductors.S1.share_towers  = P2               (the line hangs from the towers of line P2, at its own points
 *                                                       on them: a circuit's phases and shield wire on one tower; default:
 *                                                       its own towers. P2 has a tower_type and shares no other line's
 *                                                       towers; this line has as many towers and no tower_type of its
 *                                                       own. Its towers' z is measured as end_a's but from the shared
 *                                                       tower's base instead of the terrain under the point, and each of
 *                                                       its tower points lies within arm_length/2 of the tower's axis
 *                                                       and between the base and the top, the cross-arm height + peak)
 *   conductors.S1.gust_with     = P2               (with gust_type = random the line takes the random gusts of line
 *                                                       P2's spans, span by span, so that the conductors of a circuit
 *                                                       on separate lines move together; default: its own, or those of
 *                                                       the line whose towers it shares. P2 has as many spans and takes
 *                                                       no other line's gusts; not with share_towers, which already
 *                                                       shares them)
 *   conductors.S1.epsilon       = 1.5              (> 0; the width of the kernel that spreads the drag of this line's
 *                                                       nodes, its insulator strings' included, into the flow, in units of
 *                                                       the anchor level's x spacing dx; default conductors.epsilon;
 *                                                       needs drag_on_flow = true as that key does, a warning with
 *                                                       drag_on_flow = false given. The towers the line hangs from take
 *                                                       their type's epsilon)
 *
 * Every Real input must be finite. Lines, transformers and tower types share the prefix
 * conductors.<name>, so no name may be used twice.
 */
struct LineInputs
{
    std::string name;                                //!< the line's name in conductors.lines
    std::array<amrex::Real,3> end_a{{0.0, 0.0, 0.0}};   //!< the first dead end (m; x, y in ERF coordinates, z as in the input block)
    std::array<amrex::Real,3> end_b{{0.0, 0.0, 0.0}};   //!< the last dead end (m; as end_a)
    std::vector<std::array<amrex::Real,3>> towers;   //!< the suspension points between the ends, in order (m; each as end_a)
    std::vector<amrex::Real> lengths;                //!< the unstretched length of each span (m), one per span
    amrex::Real stringing_tension{0.0};              //!< instead of lengths: the still-air horizontal tension of every span (N)
    amrex::Real diameter{0.0};                       //!< the conductor's diameter (m)
    amrex::Real mass_per_length{0.0};                //!< the conductor's mass per unit length (kg/m)
    amrex::Real axial_stiffness{0.0};                //!< EA (N)
    amrex::Real drag_coefficient{1.0};               //!< the normal drag coefficient on the diameter (-)
    amrex::Real damping_ratio{0.5};                  //!< the internal damping as a fraction of critical (-)
    int segments{20};                                //!< line segments per span
    amrex::Real insulator_length{0.0};               //!< the length of each insulator string (m); 0: clamped at the towers
    amrex::Real insulator_mass{0.0};                 //!< the mass of each insulator string (kg)
    amrex::Real insulator_diameter{0.254};           //!< the disc diameter of the insulator strings (m)
    bool insulator_mass_given{false};                //!< insulator_mass was given in the inputs
    bool insulator_diameter_given{false};            //!< insulator_diameter was given in the inputs
    amrex::Real rated_breaking_strength{0.0};        //!< the conductor's rated breaking strength (N); 0: not checked
    amrex::Real insulator_rated_load{0.0};           //!< the rated mechanical load of each insulator string (N); 0: not checked
    amrex::Real hardware_rated_load{0.0};            //!< the rated load of the dead-end and tower attachment hardware (N); 0: not checked
    bool rated_breaking_strength_given{false};       //!< rated_breaking_strength was given in the inputs
    bool insulator_rated_load_given{false};          //!< insulator_rated_load was given in the inputs
    bool hardware_rated_load_given{false};           //!< hardware_rated_load was given in the inputs
    std::string output_root;                         //!< the prefix of the line's diagnostics files
    std::string tower_type;                          //!< the type of the suspension towers; empty: none
    std::string share_towers;                        //!< the line whose towers this one hangs from; empty: its own
    std::string gust_with;                           //!< the line whose random gusts this one takes; empty: see gust_owner()
    amrex::Real epsilon{0.0};                        //!< the width of the kernel that spreads the line's drag (anchor dx)
    bool epsilon_given{false};                       //!< epsilon was given; else the line takes conductors.epsilon

    /** Segments of each insulator string. */
    static constexpr int insulator_segments = 2;
    /**
     * The axial stiffness of an insulator string (N): a string of porcelain or glass discs hardly
     * stretches, and 1e7 N lets a 2.5 m string give 4 mm under a 15 kN conductor while keeping
     * MoorDyn's internal step near that of the conductor.
     */
    static constexpr amrex::Real insulator_axial_stiffness = 1.0e7;
    /** The drag coefficient of an insulator string on its disc diameter (-). */
    static constexpr amrex::Real insulator_drag_coefficient = 1.0;

    /** The number of spans: one more than the towers. */
    int num_spans () const { return static_cast<int>(towers.size()) + 1; }
    /** The attachment points in order, end_a, the towers, end_b (k = 0 .. num_spans(); m). */
    const std::array<amrex::Real,3>& point (int k) const;
    std::array<amrex::Real,3>& point (int k);
    /** Whether the conductor hangs from insulator strings at the towers. */
    bool has_insulators () const { return insulator_length > 0.0 && !towers.empty(); }
    /**
     * Where the conductor hangs at attachment point k (k = 0 .. num_spans(): 0 is end_a, num_spans()
     * is end_b): the point itself, or the bottom of the insulator string at a tower (m).
     */
    std::array<amrex::Real,3> conductor_point (int k) const;
    /** The straight distance between the attachment points of span k (k from 0) (m). */
    amrex::Real chord (int k = 0) const;
    /**
     * The straight distance between the points the conductor of span k hangs from (conductor_point:
     * the bottoms of the insulator strings at the towers), the chord its unstretched length spans (m).
     */
    amrex::Real conductor_chord (int k = 0) const;
    /**
     * With stringing_tension, set every span's unstretched length so that it hangs in still air with
     * that horizontal tension H, as a line is strung. Between the points the conductor hangs from
     * (the bottoms of the insulator strings at the towers), a chord c over a horizontal distance h,
     * the parabola of horizontal tension H under the weight w per unit length is
     * l = c + w^2 h^4 / (24 H^2 c) long, the tension along the chord H c / h stretches the line by
     * that over EA, and the unstretched length is l / (1 + H c / (EA h)). All spans of a section
     * then pull alike on the towers between them and the strings hang plumb. Called once the
     * attachments stand on the terrain; w is net of the air's buoyancy (N/m). Nothing happens
     * without stringing_tension.
     */
    void lengths_from_stringing_tension (amrex::Real w);
    /**
     * The sag (m) at mid-span of span k (k from 0) for a parabola between its attachment points whose
     * length is the span's unstretched length; 0 without slack.
     */
    amrex::Real catenary_sag (int k = 0) const;
    /** The nodes of all of the line's MoorDyn lines: the spans', then the insulator strings'. */
    int num_line_nodes () const;
    /** The prefix of span k's diagnostics files (k from 0): output_root for a single span, output_root_span<k+1> in a section. */
    std::string span_root (int k) const;
    /**
     * The name under which span k's statistics are saved in a checkpoint (k from 0): name for a
     * single span, name_span<k+1> in a section.
     */
    std::string span_name (int k) const;
};

/**
 * One transformer of conductors.transformers: a box standing on the terrain where lines end.
 * Every line end (end_a or end_b) whose x, y lies on the box's footprint is dead-ended on it, at the
 * end's placed height, which must clear the box's top. The transformer takes the pull of those
 * lines and is checked against an allowable horizontal force and overturning moment; every
 * conductor is watched for how close it comes to the box. Every transformer needs a line end on it.
 *
 *   conductors.transformers        = T1 T2
 *   conductors.T1.position         = 400. 300.     (required; m, x, y of the footprint's centre, inside the domain;
 *                                                       the box stands on the terrain there)
 *   conductors.T1.size             = 8. 5. 6.      (required; m, > 0: the box's length along x, width along y and height)
 *   conductors.T1.allowable_force  = 2.0e4         (N, >= 0, the horizontal pull the transformer may take; default 0: not checked)
 *   conductors.T1.allowable_moment = 1.5e5         (N m, >= 0, the overturning moment about the centre of its base it
 *                                                       may take; default 0: not checked)
 */
struct TransformerInputs
{
    std::string name;                                   //!< the transformer's name in conductors.transformers
    std::array<amrex::Real,2> position{{0.0, 0.0}};     //!< x, y of the footprint's centre (m)
    std::array<amrex::Real,3> size{{0.0, 0.0, 0.0}};    //!< the length along x, width along y and height (m)
    amrex::Real allowable_force{0.0};                   //!< the allowable horizontal force (N); 0: not checked
    amrex::Real allowable_moment{0.0};                  //!< the allowable overturning moment (N m); 0: not checked

    /** Whether the point x, y (m) lies on the box's footprint (its edges included). */
    bool on_footprint (amrex::Real x, amrex::Real y) const
    {
        return x >= position[0] - amrex::Real(0.5) * size[0] && x <= position[0] + amrex::Real(0.5) * size[0] &&
               y >= position[1] - amrex::Real(0.5) * size[1] && y <= position[1] + amrex::Real(0.5) * size[1];
    }
};

/** The still-air shape of an elastic catenary between two points at the same height. */
struct Catenary
{
    amrex::Real sag{0.0};                //!< mid-span sag below the chord (m)
    amrex::Real horizontal_tension{0.0}; //!< H (N)
    amrex::Real end_tension{0.0};        //!< tension at the attachments (N)
    amrex::Real stretched_length{0.0};   //!< the line's length under its own weight (m)
    bool solved{false};                  //!< the end strain T / EA is below 50 % (false: no conductor this model describes)
};

/**
 * The elastic catenary of a span of unstretched length L, net weight w per unstretched metre (as MoorDyn's
 * mass per length) and axial stiffness EA hanging between two points at the same height a chord c apart
 * (Irvine, Cable Structures, 1981): the horizontal tension H solves c = H L / EA + (2H / w) asinh(w L / 2H),
 * which rises with H, so one bisection serves a slack span (L > c) and a taut one held by its stretch
 * (L <= c) alike. The sag is (H / w)(sqrt(1 + (w L / 2H)^2) - 1) + w L^2 / 8 EA, the end tension
 * sqrt(H^2 + (w L / 2)^2). `solved` is false when the end strain reaches 50 % (no conductor; the callers stop).
 * Used to check MoorDyn's still-air shape and for the ASCE 74 wire loads.
 * Units: chord, length (m); weight_per_length (N per unstretched m, net of the air's buoyancy); axial_stiffness (N).
 * Returns: all zeros when any argument is not positive, and when not `solved`.
 */
Catenary elastic_catenary (amrex::Real chord, amrex::Real length, amrex::Real weight_per_length,
                           amrex::Real axial_stiffness);

/**
 * Whether the horizontal chord from a to b runs along +x (its y extent within 1e-6 of its x extent):
 * there a span's ERF-frame y drag is its drag across the span, so the two statistics are one.
 */
bool along_plus_x (const std::array<amrex::Real,3>& a, const std::array<amrex::Real,3>& b);

/**
 * The first entry of v that is not finite, with v read as points of per_point components each: an
 * empty string when every entry is finite, else "<what> <p> is not finite (<its components>)", p the
 * point's 0-based index. Used at the conductor module's boundaries (the wind in, MoorDyn's outputs,
 * the forces spread into the flow).
 */
template <typename T>
std::string first_nonfinite (const std::vector<T>& v, std::size_t per_point, const std::string& what)
{
    const std::size_t n = std::max<std::size_t>(per_point, 1);
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (std::isfinite(static_cast<double>(v[i]))) { continue; }
        const std::size_t p = i / n;
        std::string values;
        for (std::size_t c = p * n; c < std::min(v.size(), (p + 1) * n); ++c) {
            values += (c > p * n ? ", " : "") + std::to_string(static_cast<double>(v[c]));
        }
        return what + " " + std::to_string(p) + " is not finite (" + values + ")";
    }
    return std::string();
}

/**
 * The conductors.* inputs: the named conductor lines (one LineInputs block each), the
 * transformers and tower types, and the settings they share.
 *
 *   conductors.lines            = S1 S2        (the conductor lines' names; absent or empty: the module is off)
 *   conductors.diagnostics_dir  = conductors   (where the MoorDyn input files and the diagnostics go)
 *   conductors.diagnostics_int  = 1            (>= 1; write a diagnostics row every this many steps)
 *   conductors.anchor_level     = -1           (0 .. amr.max_level, or -1: amr.max_level; the level the lines are
 *                                                   stepped on and the wind is sampled on, which must exist)
 *   conductors.air_density      = 1.225        (kg/m^3, > 0, the fluid density MoorDyn applies the drag with)
 *   conductors.substeps         = 1            (>= 1; MoorDyn calls per ERF step for lines on fixed supports, the
 *                                                   wind held over the step; lines on moving towers take at least
 *                                                   this many coupling steps)
 *   conductors.moordyn_dt       = 0            (s, >= 0, an upper bound on MoorDyn's internal step; 0: none)
 *   conductors.moordyn_cfl      = 0.1          (in (0, 0.15], MoorDyn's Courant factor, which sets its internal step)
 *   conductors.moordyn_log_level = 2           (0 debug, 1 messages, 2 warnings, 3 errors only)
 *   conductors.surface_offset   = 10000.       (m; z_md = z_erf - surface_offset, so MoorDyn's free surface sits
 *                                                   at ERF z = surface_offset; it must exceed geometry.prob_hi[2], and
 *                                                   -surface_offset must lie below geometry.prob_lo[2])
 *   conductors.prescribed_velocity = 0. 15. 0. (m/s, ERF frame; default none: the flow is sampled; a uniform
 *                                                   wind at every line node and tower member node instead of the flow,
 *                                                   for testing)
 *   conductors.stats_start      = 0.           (s, >= 0; the running statistics of every span, insulator string
 *                                                   set, line pair, transformer and tower take the steps that end at
 *                                                   or after this time)
 *   conductors.node_output_int  = 0            (>= 0; write every node of every line, insulator strings included,
 *                                                   to <output_root>_nodes.dat every this many steps; 0: never)
 *   conductors.drag_on_flow     = false        (put the force the lines and the towers' members exert on the air,
 *                                                   minus the air's drag on them, into the flow as a momentum source
 *                                                   on the anchor level, which must then be amr.max_level)
 *   conductors.epsilon          = 2.0          (> 0; the radius eps of the kernel exp(-r^2/eps^2) that spreads that
 *                                                   force, cut off at 3 eps, in units of the anchor level's x spacing dx;
 *                                                   needs drag_on_flow = true, a warning with drag_on_flow = false given;
 *                                                   the width of every line and tower type that sets no epsilon of its
 *                                                   own, conductors.<line>.epsilon or conductors.<type>.epsilon)
 *   conductors.flashover_distance = 1.0        (m, > 0; two lines closer than this are flagged as clashing, and a
 *                                                   conductor this close to a transformer)
 *   conductors.transformers     = T1 T2        (boxes on the terrain the lines end on, each with its own block;
 *                                                   needs lines; default none)
 *   conductors.tower_types      = lattice      (lattice tower types the lines' towers can be, each with its own
 *                                                   block; needs lines; default none)
 *   conductors.asce74_wind      = 40.          (m/s, >= 0; the 3-second gust at 10 m over open country a design
 *                                                   check by ASCE 74 takes (ERF_ASCE74.H): every span's quasi-static
 *                                                   wind load, swing, blowout and tension under it are written to
 *                                                   <diagnostics_dir>/asce74.csv at start-up; default 0: none)
 *   conductors.asce74_exposure  = C            (B or C, ASCE 74's terrain exposure of that check; default C;
 *                                                   needs asce74_wind)
 *   conductors.asce74_wire_height = effective  (effective or attachment, the wire height of that check: the mean
 *                                                   of the span's attachment heights less a third of its still-air sag,
 *                                                   or the attachment height alone; default effective; needs asce74_wind)
 *   conductors.asce74_inclined_spans = true    (the check's swing about an inclined chord against the weight's
 *                                                   part normal to it, W cos beta, and the upper end's tension;
 *                                                   default true, false: the level span's swing, sag and tension
 *                                                   under the whole weight W, the effective height's still-air
 *                                                   sag the inclined chord's either way; needs asce74_wind)
 *   conductors.intensity_measure = none        (none or wind. wind: beside every component checked against a capacity,
 *                                                   record the wind at conductors.intensity_height above the ground
 *                                                   where the component stands (a span's middle, a line's attachment,
 *                                                   a tower's base), averaged over conductors.intensity_averaging_time:
 *                                                   its value at the component's peak demand-to-capacity ratio and at
 *                                                   its first exceedance, and its largest value in the run. Writes
 *                                                   <diagnostics_dir>/demand_table.csv and intensity.csv; needs a
 *                                                   capacity to check. Default none)
 *   conductors.intensity_height = 10.          (m above the ground, in (0, 500]; needs intensity_measure = wind)
 *   conductors.intensity_averaging_time = 3.   (s, in [0, 600]; 0 keeps the last sample; needs intensity_measure = wind)
 *   conductors.failure_stop     = false        (stop the run, with an abort naming the component, at the first
 *                                                   coupling step where a component's demand exceeds its capacity:
 *                                                   a conductor's rated_breaking_strength, an insulator string's or the
 *                                                   hardware's rated load, a tower type's allowable_uplift or
 *                                                   allowable_compression, or a frame member's design strength. The
 *                                                   peaks go to <diagnostics_dir>/capacities.csv and each component's
 *                                                   first exceedance to <diagnostics_dir>/failures.dat with or without
 *                                                   it; default false; true needs at least one capacity)
 *   conductors.ground_contact   = false        (the terrain holds the lines up: each line's MoorDyn system takes the
 *                                                   anchor level's terrain under it as its seafloor
 *                                                   (<diagnostics_dir>/<line>.seafloor.txt), so a conductor that sags
 *                                                   or falls to the ground rests on it instead of passing through it;
 *                                                   each span's first touch of the ground and the most nodes on it go to
 *                                                   <diagnostics_dir>/ground_contact.csv; default false)
 *   conductors.gust_type        = factor       (none, factor, event or random; default none (ERF_Gusts.H). Every
 *                                                   type but none writes every span's gust factor and peak wind load
 *                                                   from the mean wind and the k-equation RANS's k along it to
 *                                                   <diagnostics_dir>/gusts.csv; event adds one travelling 1 - cos
 *                                                   gust to the wind the lines and towers take, random a random gust
 *                                                   per span and per tower, both written to
 *                                                   <diagnostics_dir>/gust_series.dat. Needs erf.rans_type = kEqn on the
 *                                                   anchor level and the flow's wind, not prescribed_velocity, and
 *                                                   every span with a horizontal extent; event and random cannot be
 *                                                   used with drag_on_flow)
 *   conductors.gust_sigma_factor = 1.39        (> 0; c in sigma_u = c sqrt(k); default 2.5 Cmu0, Cmu0 the RANS
 *                                                   constant erf.Cmu0; needs a gust_type)
 *   conductors.gust_peak_factor = 2.7          (> 0; g in the gust factor 1 + 2 g I_n sqrt(B), the peak wind
 *                                                   U + g sigma_u and the event's amplitude g sigma_u; default ASCE
 *                                                   74's 2.7; needs a gust_type)
 *   conductors.gust_span_length_scale = 67.056 (m, > 0; L_s in the spans' B = 1/(1 + 0.8 L/L_s) and the towers'
 *                                                   B_t = 1/(1 + 0.375 h/L_s); default ASCE 74's 220 ft of open
 *                                                   country, exposure C, or that of asce74_exposure when it is
 *                                                   given (51.816 m, 170 ft, for B); needs a gust_type)
 *   conductors.gust_event_time  = 600.         (s; when the event's front crosses gust_event_origin; needs and
 *                                                   is needed by gust_type = event)
 *   conductors.gust_event_speed = 12.          (m/s, > 0; how fast the front moves; needs and is needed by
 *                                                   gust_type = event)
 *   conductors.gust_event_direction = 0.       (degrees counterclockwise from +x, the way the front moves;
 *                                                   default 0; needs gust_type = event)
 *   conductors.gust_event_origin = 600. 500.   (m, x y; a point the front crosses at gust_event_time; default
 *                                                   the domain's centre; needs gust_type = event)
 *   conductors.gust_event_duration = 10.5      (s, > 0; how long the gust lasts at a point; default 10.5 s, the
 *                                                   duration of IEC 61400-1's extreme operating gust, not its shape
 *                                                   or amplitude; needs gust_type = event)
 *   conductors.gust_seed        = 1            (>= 0; picks the random gusts' series; default 1; needs
 *                                                   gust_type = random)
 *   conductors.gust_integral_length = 170.     (m, > 0; the random gusts' integral length L_u; default IEC
 *                                                   61400-1's 8.1 Lambda_1 at each span's attachment height and each
 *                                                   tower's mean drag-node height; needs gust_type = random)
 *
 * Every Real input must be finite.
 */
struct ConductorInputs
{
    bool active{false};                                 //!< conductors.lines lists at least one line
    std::vector<LineInputs> lines;                      //!< the lines, in the order of conductors.lines
    std::vector<TransformerInputs> transformers;        //!< the transformers, in the order of conductors.transformers
    std::vector<erf_towers::TowerType> tower_types;     //!< the tower types, in the order of conductors.tower_types
    std::string diagnostics_dir{"conductors"};
    int diagnostics_int{1};
    int anchor_level{-1};
    amrex::Real air_density{1.225};                     //!< (kg/m^3)
    int substeps{1};
    amrex::Real moordyn_dt{0.0};                        //!< (s)
    amrex::Real moordyn_cfl{0.1};                       //!< (-)
    int moordyn_log_level{2};
    /**
     * The largest Courant factor accepted. MoorDyn's own default, 0.5, and 0.25 make a 300 m ACSR
     * (aluminium conductor, steel reinforced) span diverge within a second in a 20 m/s wind; 0.2 runs
     * but more than doubles the still-air sag; 0.15 and 0.1 agree to the millimetre.
     */
    static constexpr amrex::Real max_moordyn_cfl = 0.15;
    amrex::Real surface_offset{10000.0};                //!< (m)
    bool has_prescribed_velocity{false};
    std::array<amrex::Real,3> prescribed_velocity{{0.0, 0.0, 0.0}};   //!< (m/s, ERF frame)
    amrex::Real stats_start{0.0};                       //!< (s)
    int node_output_int{0};
    bool drag_on_flow{false};
    amrex::Real epsilon{2.0};                           //!< (units of the anchor level's dx)
    amrex::Real flashover_distance{1.0};                //!< (m)
    bool failure_stop{false};                           //!< stop the run at the first capacity exceedance
    std::string intensity_measure{"none"};              //!< none or wind: the wind beside each component's peak demand
    amrex::Real intensity_height{10.0};                 //!< the height above the ground the wind is taken at (m)
    amrex::Real intensity_averaging_time{3.0};          //!< the time the wind is averaged over (s); 0: the last sample
    bool intensity_on () const { return intensity_measure == "wind"; }
    bool ground_contact{false};                         //!< the terrain holds the lines up (MoorDyn's 3-D seafloor)
    /**
     * The ground's contact stiffness per unit area of conductor on it (Pa/m), MoorDyn's own default:
     * a node of diameter d and length l pressed delta into the ground is pushed up with
     * ground_stiffness d l delta. A conductor of 1.6 kg/m and 28 mm rests 0.2 mm deep.
     */
    static constexpr amrex::Real ground_stiffness = 3.0e6;
    /**
     * The contact damping (Pa s/m) that makes a node of the line type of mass mu per metre (kg/m) and
     * diameter d (m) critically damped on the ground: 2 sqrt(ground_stiffness mu / d). MoorDyn's
     * default, 3e5 Pa s/m, is set for mooring chain; on a conductor it is eleven times critical, and the
     * explicit time integration diverges as the line lands.
     */
    static amrex::Real ground_damping (amrex::Real mu, amrex::Real d);
    /** The contact frequency sqrt(ground_stiffness d / mu) (rad/s) of a node of that line type on the ground. */
    static amrex::Real ground_frequency (amrex::Real mu, amrex::Real d);
    /**
     * The fastest rate (1/s) of a node of that line type on the ground under the contact damping cb
     * (Pa s/m): with omega its ground_frequency and a = cb d / (2 mu), the largest magnitude of the
     * roots of the damped oscillator, twice max(omega, a + sqrt(max(a^2 - omega^2, 0))). It is
     * 2 omega for a critically damped type (the conductor), 2 omega for an under-damped one, and
     * tends to 2 cb d / mu for an over-damped one (a light insulator string under the conductor's cb).
     */
    static amrex::Real ground_rate (amrex::Real mu, amrex::Real d, amrex::Real cb);
    /**
     * The rate times MoorDyn's internal step allowed. MoorDyn's second-order Runge-Kutta scheme
     * (explicit midpoint) is stable for a linear, critically damped oscillator while the rate times
     * the step is below 4, but the one-sided contact chatters earlier: 300 m spans of 0.5 to 3 kg/m
     * and 15 to 40 mm dropped onto flat ground landed cleanly in every case at 2, sank 1 to 11 cm
     * transiently at 2.5 and 3, and failed at 3.5 for some and at 3.9 for all. The check keeps 1.5.
     * An almost undamped type (damping ratio below 0.1) would not be covered, but needs
     * mu / d > 5800 kg/m^2 under the conductor's cb.
     */
    static constexpr amrex::Real max_ground_step = 1.5;
    /**
     * The most points of a line's seafloor file: past it the terrain is taken at every stride-th node
     * of the anchor level. MoorDyn holds the grid on every rank (8 MB at this size) and reads it at
     * start-up.
     */
    static constexpr long max_seafloor_points = 1000000;
    amrex::Real asce74_wind{0.0};                       //!< (m/s); 0: no ASCE 74 table
    std::string asce74_exposure;                        //!< B or C; empty: C
    std::string asce74_wire_height{"effective"};        //!< effective (attachment height less a third of the sag) or attachment
    bool has_asce74_wire_height{false};                 //!< asce74_wire_height was given (it needs asce74_wind)
    bool asce74_inclined_spans{true};                   //!< W cos(beta) across an inclined chord, the upper end's tension; false: level, W
    bool has_asce74_inclined_spans{false};              //!< asce74_inclined_spans was given (it needs asce74_wind)
    std::string gust_type{"none"};                      //!< none, factor, event or random
    bool has_gust_sigma_factor{false};
    amrex::Real gust_sigma_factor{0.0};                 //!< (-); not given: 2.5 Cmu0 of the k-equation RANS
    bool has_gust_peak_factor{false};
    amrex::Real gust_peak_factor{2.7};                  //!< (-)
    bool has_gust_span_length_scale{false};
    amrex::Real gust_span_length_scale{67.056};         //!< (m), 220 ft
    bool has_gust_event_time{false};
    amrex::Real gust_event_time{0.0};                   //!< (s)
    bool has_gust_event_speed{false};
    amrex::Real gust_event_speed{0.0};                  //!< (m/s)
    bool has_gust_event_direction{false};
    amrex::Real gust_event_direction{0.0};              //!< (degrees counterclockwise from +x)
    bool has_gust_event_origin{false};
    std::array<amrex::Real,2> gust_event_origin{{0.0, 0.0}};   //!< (m, x y); not given: the domain's centre
    bool has_gust_event_duration{false};
    amrex::Real gust_event_duration{10.5};              //!< (s)
    bool has_gust_seed{false};
    int gust_seed{1};
    bool has_gust_integral_length{false};
    amrex::Real gust_integral_length{0.0};              //!< (m); not given: IEC 61400-1's at each span's and tower's height
    /** The gust type (validate_settings() has checked that it parses; None when it does not). */
    GustType gust () const;
    /** Whether gust_type is anything but none. */
    bool gusts_on () const { return gust() != GustType::None; }
    /** Whether a gust is added to the wind the lines and towers take: gust_type event or random. */
    bool gusts_in_wind () const { return gust() == GustType::Event || gust() == GustType::Random; }

    /**
     * The line whose random gusts a line takes: its gust_with; else, for a line sharing towers, the towers' owner's
     * gust_with or the owner itself; else itself. That line takes no other's.
     */
    const LineInputs& gust_owner (const LineInputs& s) const;
    /** The line a line's towers belong to: the line whose towers it shares, or itself. */
    const LineInputs& tower_owner (const LineInputs& s) const;
    /** The tower type of a line's towers (those of the line it shares them with), or null when they have none. */
    const erf_towers::TowerType* tower_type (const LineInputs& s) const;
    /**
     * The width of the kernel that spreads the drag of a line's nodes (in units of the anchor level's x spacing):
     * the line's own epsilon when it gives one, else conductors.epsilon.
     */
    amrex::Real line_epsilon (const LineInputs& s) const { return s.epsilon_given ? s.epsilon : epsilon; }
    /** The same for the drag nodes of a tower of type t: the type's own epsilon when it gives one, else conductors.epsilon. */
    amrex::Real tower_epsilon (const erf_towers::TowerType& t) const { return t.epsilon_given ? t.epsilon : epsilon; }
    /**
     * The first line or tower type that gives an epsilon of its own, as its key conductors.<name>.epsilon (lines
     * first, in their order), or an empty string when none does.
     */
    static std::string first_own_epsilon (const ConductorInputs& in);
    /**
     * Whether every line, and the tower type of every line's towers, gives an epsilon of its own, so that no node
     * takes conductors.epsilon.
     */
    static bool every_width_own (const ConductorInputs& in);
    /**
     * The tower types that give an epsilon of their own while no line's towers are of that type, as their keys
     * conductors.<type>.epsilon: read, but taken by no node.
     */
    static std::vector<std::string> unused_type_widths (const ConductorInputs& in);
    /**
     * The widths below one cell that some node takes, as their keys (conductors.epsilon, then
     * conductors.<line>.epsilon, then conductors.<type>.epsilon of the types a line's towers are): under one
     * cell the kernel's shape depends on where a node sits in its cell, and a tower's drag gathered on a few faces can
     * overshoot the wind next to it in one explicit step.
     */
    static std::vector<std::string> narrow_widths (const ConductorInputs& in);
    /** Whether a line's towers move (their type has a frequency or a frame): MoorDyn then takes their cross-arms as coupled points. */
    bool towers_move (const LineInputs& s) const
    {
        const erf_towers::TowerType* t = tower_type(s);
        return t != nullptr && t->moves();
    }

    /** Read conductors.*; aborts with the offending key on a bad value. */
    static ConductorInputs read ();

    /**
     * The range checks of one line (every Real finite, then each key's range): an empty string when it
     * is acceptable, else the message to abort with, naming the key. With check_slack the lengths are
     * checked against the chords as well (validate_slack); read() leaves that to after the attachment
     * points are placed on the terrain.
     */
    static std::string validate_line (const LineInputs& s, bool check_slack = true);
    /**
     * Whether every span's unstretched length exceeds the distance between the points the conductor
     * hangs from (conductor_point: the attachment points, or the bottoms of the insulator strings at
     * the towers): an empty string, else the message. A line strung to a tension is not checked: a
     * short span is taut, held by its stretch. on_terrain says the points stand at their absolute
     * heights, whose chords differ from those of the heights above the terrain wherever the ground
     * under the two ends differs; it changes only the message.
     */
    static std::string validate_slack (const LineInputs& s, bool on_terrain = false);
    /** Whether every line has an output_root of its own: an empty string, else the message naming the second line. */
    static std::string validate_output_roots (const std::vector<LineInputs>& lines);
    /** Whether a line's tower_type names one of the tower types and the line has towers: an empty string, else the message. */
    static std::string validate_tower_type (const LineInputs& s, const std::vector<erf_towers::TowerType>& types);
    /**
     * Whether every line that shares towers names a line of the run that has a tower_type, shares
     * no other line's, and has as many towers, and sets no tower_type of its own: an empty string,
     * else the message.
     */
    static std::string validate_shared_towers (const std::vector<LineInputs>& lines);
    /** The range checks of one transformer: an empty string when it is acceptable, else the message to abort with. */
    static std::string validate_transformer (const TransformerInputs& t);
    /**
     * The range checks of the shared settings (every Real finite, then each key's range) and of each
     * line's mass_per_length against the air it displaces: an empty string when acceptable, else the
     * message, naming the key.
     */
    static std::string validate_settings (const ConductorInputs& in);
    /**
     * The solver settings the conductor module supports: an empty string when they hold, else the message.
     * @param max_level amr.max_level
     * @param anchor_level the level the lines are stepped on, already resolved (not -1)
     */
    static std::string validate_solver (int max_level, int anchor_level);
    /**
     * Whether MoorDyn's frame holds the whole domain: surface_offset must exceed prob_hi_z, so that
     * every node lies below MoorDyn's free surface at ERF z = surface_offset (MoorDyn applies no fluid
     * load above it), and -surface_offset must lie below prob_lo_z, so that no node reaches MoorDyn's
     * flat bottom at ERF z = -surface_offset. Units: m. Returns an empty string, else the message.
     */
    static std::string validate_frame (amrex::Real surface_offset, amrex::Real prob_lo_z, amrex::Real prob_hi_z);
    /** The level the lines are stepped on: requested, or max_level when requested is -1; Conductors::require_anchor_level
     *  checks that it exists. */
    static int resolve_anchor_level (int requested, int max_level) { return requested < 0 ? max_level : requested; }
};

} // namespace erf_conductors

#endif
