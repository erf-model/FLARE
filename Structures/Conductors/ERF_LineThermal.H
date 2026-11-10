#ifndef ERF_LINE_THERMAL_H
#define ERF_LINE_THERMAL_H

/**
 * @file ERF_LineThermal.H
 * @brief The conductors' temperature under a fire, and the sag it adds
 *
 * conductors.thermal_model = ieee738 gives every node of every span a temperature T [C] that follows
 * the heat balance of IEEE Std 738-2023 per unit length (eq 3),
 *
 *   m c_p dT/dt = q_J(T) + q_s + q_f - q_c(T, T_a, V) - q_r(T, T_a)          [W/m]
 *
 * with the Joule heating q_J of the line's current, the solar heating q_s, the radiation q_f the
 * conductor absorbs from the fire, and the convective and radiative losses q_c and q_r to the air
 * around the node. T_a is ERF's air temperature at the node and V its wind there, so a plume that the
 * fire heats (erf.fire.coupling_type lagged or synchronous) heats the conductor by convection as well.
 * The conductor radiates to surroundings at the air temperature its node had at the first step, not to
 * the plume (hot air is a poor emitter: treating it as a black enclosure would overstate the heating).
 * The fire radiates from every burning fire cell as a point source of chi_r q'' dA at
 * conductors.thermal.source_height above the ground, radiating equally in every direction. The
 * balance is stepped implicitly (backward Euler) with the anchor level's step.
 *
 * Each span's unstretched length follows the mean temperature of its conductor,
 * L = L_ref (1 + alpha_T (Tbar - T_ref)), where T_ref is its mean temperature at the first step and
 * L_ref its length then, so the span sags as it heats (MoorDyn sets the new length before the step).
 * The insulator strings are not heated.
 *
 * The terms, their equations and their limits are in Docs/sphinx_doc/theory/Conductors.rst
 * ("Conductor temperature under a fire"); the inputs in Docs/sphinx_doc/ConductorInputs.rst.
 */

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <array>
#include <memory>
#include <string>
#include <vector>

struct StructureFire;

namespace erf_conductors {

class ConductorLine;

/// The air's film properties of IEEE Std 738-2023 (SI): viscosity [kg/(m s)], density [kg/m^3], conductivity [W/(m K)].
struct AirFilm { double mu{0.0}, rho{0.0}, k{0.0}; };

/**
 * The air's properties at the film temperature T_film [C] (the mean of the conductor's and the air's) and
 * the elevation [m] above sea level: IEEE Std 738-2023 eqs 13, 15 and 17.
 */
AirFilm ieee738_air (double T_film, double elevation);

/** The wind direction factor K_angle of IEEE Std 738-2023 eq 19; phi [rad] is the angle between the wind and the conductor's axis. */
double ieee738_wind_angle_factor (double phi);

/**
 * The convective heat the conductor gives the air [W/m], IEEE Std 738-2023 eq 6: the larger in magnitude of forced
 * convection (the larger of eqs 9 and 10) and natural convection (eq 7), at the film temperature
 * (T + T_a)/2. Negative when the air is hotter than the conductor: the same correlations with the sign
 * of T - T_a (an extension of the standard, which states them for a conductor hotter than the air).
 *
 * @param[in] D          diameter [m]
 * @param[in] V          wind speed [m/s] (>= 0)
 * @param[in] phi        angle between the wind and the conductor's axis [rad]
 * @param[in] T, T_a     the conductor's and the air's temperature [C]
 * @param[in] elevation  above sea level [m]
 */
double ieee738_convection (double D, double V, double phi, double T, double T_a, double elevation);
/** Forced convection alone [W/m]: the larger of IEEE Std 738-2023 eqs 9 and 10 (signed as ieee738_convection). */
double ieee738_forced_convection (double D, double V, double phi, double T, double T_a, double elevation);
/** Natural convection alone [W/m]: IEEE Std 738-2023 eq 7 with |T - T_a|^1.25 and the sign of T - T_a. */
double ieee738_natural_convection (double D, double T, double T_a, double elevation);
/** The heat the conductor radiates to surroundings at T_s [W/m]: IEEE Std 738-2023 eq 21 (there with T_s = T_a). */
double ieee738_radiation (double D, double emissivity, double T, double T_s);

/** One conductor's thermal data: conductors.<line>.* */
struct LineThermalProps
{
    double diameter{0.0};            //!< D [m] (conductors.<line>.diameter)
    double absorptivity{0.0};        //!< alpha, of the sun's and the fire's radiation [-]
    double emissivity{0.0};          //!< epsilon [-]
    double heat_capacity{0.0};       //!< m c_p per unit length [J/(m K)]
    double thermal_expansion{0.0};   //!< alpha_T [1/K]
    double current{0.0};             //!< I [A]
    std::array<double,2> resistance{{0.0, 0.0}};               //!< R at the two temperatures below [ohm/m]
    std::array<double,2> resistance_temperatures{{25.0, 75.0}}; //!< [C]

    /** R(T) [ohm/m]: linear through the two points (IEEE Std 738-2023 eq 33). */
    double resistance_at (double T) const;
};

/** The settings every line shares: conductors.thermal.* */
struct ThermalSettings
{
    double radiative_fraction{0.0};  //!< chi_r: the fraction of the fire's heat release it radiates [-]
    bool   radiative_fraction_given{false};
    double source_height{0.0};       //!< height of the point sources above the ground [m]
    double cutoff_radius{0.0};       //!< sources farther than this from a node are left out [m]; 0: none
    double elevation{0.0};           //!< the lines' elevation above sea level [m], for the air's density
    double solar_flux{0.0};          //!< the sun's flux on the conductor, normal to it [W/m^2]
    double reference_temperature{0.0}; //!< the temperature at which the spans have their input lengths [C]
    bool   reference_temperature_given{false}; //!< not given: the spans' mean steady temperature at the first step
};

/**
 * The net heat a conductor gains per unit length [W/m]: q_J + q_s + q_f - q_c - q_r.
 *
 * @param[in] p       the conductor's thermal data
 * @param[in] s       the shared settings (elevation, solar flux)
 * @param[in] T, T_a  the conductor's and the air's temperature [C] (T_a for convection)
 * @param[in] T_s     the temperature of the surroundings the conductor radiates to [C]
 * @param[in] V, phi  the wind speed [m/s] and its angle to the conductor's axis [rad]
 * @param[in] E_fire  the fire's irradiance on the conductor [W/m^2] (fire_irradiance)
 */
double heat_gain (const LineThermalProps& p, const ThermalSettings& s, double T, double T_a, double T_s, double V, double phi,
                  double E_fire);

/**
 * The temperature [C] at which heat_gain() is zero (the steady state), with the surroundings at the air's
 * temperature; aborts, naming `where`, when none is found below 3000 C.
 */
double steady_temperature (const LineThermalProps& p, const ThermalSettings& s, double T_a, double V, double phi, double E_fire,
                           const std::string& where = std::string());

/**
 * One backward-Euler step of the balance: the T at which m c_p (T - T_old) = dt heat_gain(T), with the air,
 * the surroundings, the wind and the fire held at their values over the step. Aborts, naming `where`, when no
 * such T is found below 3000 C. The bracket is found by doubling steps from T_old and the root by bisection;
 * g(T) = m c_p (T - T_old) - dt heat_gain(T) need not be monotone (steps of 1e4 s at high current), but has
 * one root over the range of the inputs.
 */
double implicit_step (const LineThermalProps& p, const ThermalSettings& s, double T_old, double dt,
                      double T_a, double T_s, double V, double phi, double E_fire, const std::string& where = std::string());

/** A burning fire cell as a point source: where it radiates from [m] and its heat release q'' dA [W]. */
struct FireSource { double x{0.0}, y{0.0}, z{0.0}, Q{0.0}; };

/**
 * The fire cells whose heat release is positive, as point sources at their cell centres (z left at 0 for
 * the caller to set), the same list on every rank, sorted by y, x and Q. Only cells whose centre lies
 * horizontally within `reach` of the box [lo, hi] (the line nodes' bounds) are kept; reach < 0 keeps every
 * cell, and so does a direction whose period is positive (a periodic direction, where the distance wraps).
 * Collective: every rank calls it.
 */
std::vector<FireSource> gather_fire_sources (const StructureFire& fire, const std::array<double,3>& lo,
                                             const std::array<double,3>& hi, double reach,
                                             const std::array<double,2>& period = {{0.0, 0.0}});

/**
 * The fire's irradiance on a conductor at each node [W/m^2]: sum_j chi_r Q_j sin(psi_j) / (4 pi R_j^2),
 * where R_j is the distance from source j to the node and psi_j the angle between that ray and the
 * conductor's axis there (so the conductor intercepts D sin psi per unit length). R_j^2 is not taken
 * below r_min2 (the square of the radius of a disc of one fire cell's area); sources farther than
 * cutoff (> 0) are left out. In a direction with a positive period (a periodic domain's length) the
 * distance is the nearest image's.
 *
 * @param[in] nodes     3 coordinates per node [m]
 * @param[in] tangents  the conductor's unit axis at each node, 3 per node
 */
std::vector<double> fire_irradiance (const std::vector<double>& nodes, const std::vector<double>& tangents,
                                     const std::vector<FireSource>& sources, double radiative_fraction,
                                     double cutoff, double r_min2, const std::array<double,2>& period = {{0.0, 0.0}});

/** The unit axis of a conductor at each of a span's nodes (3 coordinates per node): central differences, one-sided at the ends. */
std::vector<double> span_tangents (const std::vector<double>& nodes);

/** The trapezoid mean of a span's node values (nodes equally spaced along the unstretched conductor). */
double span_mean (const std::vector<double>& v);

/**
 * The conductors' temperatures: created by Conductors when conductors.thermal_model = ieee738. It reads
 * its inputs, holds every span node's temperature, steps the balance, sets the spans' lengths, writes
 * <diagnostics_dir>/conductor_temperature.dat (per span and output step: the mean and hottest node's temperature,
 * that node's number (0-based over the line's nodes), air temperature, wind speed, angle to the axis and fire
 * irradiance, the fire's largest irradiance on the span, the thermal strain and the reference temperature) and its
 * checkpoint <chk>/conductors/thermal.
 */
class LineThermal
{
public:
    /** Whether conductors.thermal_model selects a thermal model (none, the default, or ieee738; anything else aborts). */
    static bool selected ();

    /**
     * Read conductors.thermal.* and conductors.<line>.* for every line (diameters from the line inputs) and
     * abort, naming the key, on a value out of range (check()) or a line named "thermal" (its keys would be
     * conductors.thermal.*).
     */
    LineThermal (const std::vector<std::string>& line_names, const std::vector<double>& diameters,
                 const std::string& diagnostics_dir);
    /** For tests: the settings and data given directly, checked as the constructor checks them. */
    LineThermal (const ThermalSettings& s, const std::vector<std::string>& line_names,
                 const std::vector<LineThermalProps>& props, const std::string& diagnostics_dir);

    /** The first problem with the settings or a line's data (the key and its range), or an empty string. */
    static std::string check (const ThermalSettings& s, const std::string& line, const LineThermalProps& p);

    const ThermalSettings& settings () const { return m_s; }
    const LineThermalProps& props (std::size_t line) const { return m_p[line]; }

    /**
     * Step every span node's temperature over [time, time + dt] and set each span's unstretched length to
     * its mean temperature (the lines step after this). On the first call of a fresh start each node starts
     * at its steady temperature without the fire, its surroundings take the air temperature it has then, and
     * the spans' reference temperatures (reference_temperature, or their mean steady temperature) and lengths
     * are taken then. The fire's irradiance at each node comes from `sources` (fire_irradiance, with this
     * class's radiative fraction and cutoff); the air and the wind are held over the step.
     *
     * @param[in,out] lines  the lines (their node positions are read, their spans' lengths set)
     * @param[in] wind       per line, the wind at its kinematics points (the nodes first, 3 per point) [m/s]
     * @param[in] T_air      per line, the air temperature at each node [C]
     * @param[in] sources    the fire's point sources (gather_fire_sources, z set); empty: no fire
     * @param[in] r_min2     the floor on the squared distance to a source [m^2] (fire_irradiance)
     * @param[in] period     the domain's length in x and y where periodic, else 0 (fire_irradiance)
     * @param[in] write_row  write the log's rows for this step (the conductors' diagnostics cadence)
     */
    void advance (double time, double dt, std::vector<std::unique_ptr<ConductorLine>>& lines,
                  const std::vector<std::vector<amrex::Real>>& wind, const std::vector<std::vector<double>>& T_air,
                  const std::vector<FireSource>& sources, double r_min2, const std::array<double,2>& period, bool write_row);

    /** Node m's temperature of line i [C] (0-based over the line's nodes; insulator-string nodes hold their line's air temperature). */
    double temperature (std::size_t i, unsigned m) const { return m_T[i][m]; }
    /** Span k's mean conductor temperature now, and its reference temperature [C]. */
    double span_temperature (std::size_t i, int k) const;
    double span_reference_temperature (std::size_t i, int k) const { return m_Tref[i][static_cast<std::size_t>(k)]; }
    /** Span k's reference unstretched length [m] (at span_reference_temperature). */
    double span_reference_length (std::size_t i, int k) const { return m_Lref[i][static_cast<std::size_t>(k)]; }
    /** The fire's irradiance [W/m^2] at node m of line i in the last step. */
    double node_irradiance (std::size_t i, unsigned m) const { return m_E[i][m]; }
    /** The temperature of the surroundings node m of line i radiates to [C]: its air temperature at the first step. */
    double surroundings_temperature (std::size_t i, unsigned m) const { return m_Ts[i][m]; }
    bool started () const { return m_started; }

    /**
     * Save the node temperatures and surroundings temperatures, and the spans' mean, reference temperatures and
     * reference lengths, as <dir>/thermal (I/O rank; dir = <chk>/conductors).
     */
    void write_checkpoint (const std::string& dir) const;
    /**
     * Read them back on a restart and set every span's length again (MoorDyn's saved state does not carry the
     * unstretched length).
     * Aborts when the checkpoint's lines or nodes do not match these lines'. A checkpoint without the file
     * (written without the thermal model) starts the temperatures afresh at the first step.
     */
    void restore (const std::string& dir, std::vector<std::unique_ptr<ConductorLine>>& lines);
    /** The log this class writes, for Conductors::log_files() (trimmed to the restart time on a restart). */
    std::string log_file () const { return m_dir + "/conductor_temperature.dat"; }

private:
    void set_lengths (std::vector<std::unique_ptr<ConductorLine>>& lines) const;
    void write_log (double time, const std::vector<std::unique_ptr<ConductorLine>>& lines, bool first) const;

    ThermalSettings m_s;
    std::vector<std::string> m_names;
    std::vector<LineThermalProps> m_p;
    std::string m_dir;
    bool m_started{false};
    bool m_restored{false};   //!< restore() ran: the log continues rather than starting afresh
    std::vector<std::vector<double>> m_T;     //!< per line, per node [C]
    std::vector<std::vector<double>> m_Ts;    //!< per line, per node: the surroundings' temperature [C]
    std::vector<std::vector<double>> m_E;     //!< per line, per node: the fire's irradiance in the last step [W/m^2]
    std::vector<std::vector<std::array<double,3>>> m_air; //!< per line, per node: the air temperature [C], wind speed [m/s] and angle [rad] of the last step
    std::vector<std::vector<double>> m_span_T; //!< per line, per span: the mean temperature now [C]
    std::vector<std::vector<double>> m_Tref;  //!< per line, per span [C]
    std::vector<std::vector<double>> m_Lref;  //!< per line, per span [m]
};

} // namespace erf_conductors

#endif
