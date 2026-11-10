#ifndef ERF_TOWER_THERMAL_H
#define ERF_TOWER_THERMAL_H

/**
 * @file ERF_TowerThermal.H
 * @brief The steel temperature of the lattice towers' members under a fire, and the frame it softens
 *
 * conductors.thermal.towers = true (with conductors.thermal_model = ieee738) gives every member of every
 * tower that bends as a frame model (frame_file or frame_panels) one temperature theta [C], which follows
 * the lumped heat balance of unprotected steel in EN 1993-1-2 (section 4.2.5.1) per unit length,
 *
 *   rho_a A dH/dt = P_x [alpha_c (T_a - theta) + eps sigma (T_s^4 - theta^4)] + eps W E          [W/m]
 *
 * where H(theta) is the steel's enthalpy from its specific heat c_a(theta) (EN 1993-1-2 eqs 3.2a-d),
 * rho_a = 7850 kg/m^3, A the member's area, P_x the perimeter that exchanges heat with the air (the box
 * perimeter of EN 1993-1-2, k_sh times the section's perimeter), alpha_c the convection coefficient, eps
 * the steel's emissivity, T_a ERF's air temperature along the member, T_s the air temperature the member
 * had at the first step (its surroundings, as for the conductors), W the width that intercepts the fire's
 * radiation (the mean projected width of the section's convex hull, its perimeter over pi) and E the
 * fire's irradiance of ERF_LineThermal.H (point sources at each burning fire cell, projected on the
 * member's axis). The balance is stepped implicitly in the enthalpy (backward Euler), so the spike of c_a
 * near 735 C is crossed without a step limit.
 *
 * A tower's frame is rebuilt at its members' temperatures (FrameTower::set_temperature, EN 1993-1-2's
 * k_y and k_E of Table 3.1) when any member's k_y or k_E has moved by more than the fraction
 * conductors.thermal.tower_reduction_tolerance from the frame's; the frame's stiffness, sag and first
 * frequency then follow the heated steel. The members' strength checks take their own temperatures every
 * step (FrameTower::set_check_temperature). The terms, their equations and their limits
 * are in Docs/sphinx_doc/theory/Conductors.rst ("Tower steel temperature under a fire"); the inputs in
 * Docs/sphinx_doc/ConductorInputs.rst.
 */

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace erf_towers { class FrameTower; struct FrameInputs; struct MemberDesign; }

namespace erf_conductors {

/** The density of structural steel [kg/m^3], EN 1993-1-2 3.2.2. */
inline constexpr double steel_density = 7850.0;

/**
 * The specific heat of carbon steel c_a [J/(kg K)] at theta [C], EN 1993-1-2 eqs 3.2a to 3.2d: a cubic up to
 * 600 C, the peak of 5000 at 735 C, 650 from 900 C. Below 20 C, the value at 20 C; above 1200 C, 650.
 */
double steel_specific_heat (double theta);

/** The steel's enthalpy per unit mass from 20 C [J/kg]: the integral of steel_specific_heat() from 20 C to theta, in closed form. */
double steel_enthalpy (double theta);

/**
 * A member's section, as the heat balance needs it: the area A [m^2], the perimeter P_x that exchanges heat with
 * the air [m] (EN 1993-1-2's box perimeter) and the width W that intercepts the fire's radiation [m].
 */
struct MemberSection { double area{0.0}, perimeter{0.0}, width{0.0}; };

/**
 * The section of frame member m (0-based, in.members order): from its design data when it is an equal-leg angle
 * b x t (designs, one row per member, may be null): A = t (2b - t), P_x = 4b, W = ((2 + sqrt 2) b + (2 - sqrt 2) t) / pi;
 * else from
 * its frame section: circular D, t (0: solid) A = pi (D^2 - (D - 2t)^2)/4, P_x = pi D, W = D; rectangular
 * Sa, Sb, t (0: solid) A = Sa Sb - (Sa - 2t)(Sb - 2t), P_x = 2 (Sa + Sb), W = P_x / pi. An arbitrary section
 * without an angle's design data has no perimeter: returns false and says why in err.
 */
bool member_section (const erf_towers::FrameInputs& in, const std::vector<erf_towers::MemberDesign>* designs, std::size_t m,
                     MemberSection& sec, std::string& err);

/** The settings of the towers' heating: conductors.thermal.tower_* */
struct TowerThermalSettings
{
    double emissivity{0.7};            //!< eps_m of the steel [-] (EN 1993-1-2 2.2(2)); also its absorptivity
    double convection{25.0};           //!< alpha_c [W/(m^2 K)] (EN 1991-1-2 3.2.1(2), the standard fire's)
    double reduction_tolerance{0.01};  //!< rebuild a frame when a member's k_y or k_E moves by more than this fraction [-]
};

/** The first problem with the settings (the key and its range), or an empty string. */
std::string check (const TowerThermalSettings& s);

/**
 * The net heat a member gains per unit length [W/m]: P_x [alpha_c (T_a - theta) + eps sigma (T_s^4 - theta^4)] +
 * eps W E, temperatures in C, E the fire's irradiance [W/m^2] (fire_irradiance() of ERF_LineThermal.H).
 */
double member_heat_gain (const MemberSection& sec, const TowerThermalSettings& s, double theta, double T_a, double T_s, double E);

/**
 * One backward-Euler step in the enthalpy: the theta at which rho_a A (H(theta) - H(theta_old)) = dt
 * member_heat_gain(theta), the air, the surroundings and the fire held over the step. The left side rises and the
 * gain falls with theta, so the root is unique; it is bracketed by doubling steps from theta_old (the trial points
 * held within [-200, 5000] C) and found by bisection to 1e-10 K. Aborts, naming where, when no root lies in that range.
 */
double member_implicit_step (const MemberSection& sec, const TowerThermalSettings& s, double theta_old, double dt,
                             double T_a, double T_s, double E, const std::string& where = std::string());

/** The temperature [C] at which member_heat_gain() is zero (the steady state). */
double member_steady_temperature (const MemberSection& sec, const TowerThermalSettings& s, double T_a, double T_s, double E,
                                  const std::string& where = std::string());

/**
 * The largest relative change of EN 1993-1-2's k_y or k_E between two sets of member temperatures (one per member):
 * |k(a) - k(b)| / max(k(a), k(b)), 0 where both are 0.
 */
double reduction_change (const std::vector<double>& a, const std::vector<double>& b);

/**
 * The 3 Gauss-Legendre points along a member from x_a to x_b (3 coordinates each), in the order
 * gauss_weights() gives their weights: at 1/2 - sqrt(3/5)/2, 1/2 and 1/2 + sqrt(3/5)/2 of its length.
 */
std::array<std::array<double,3>,3> gauss_points (const std::array<double,3>& x_a, const std::array<double,3>& x_b);
/** Their weights, 5/18, 8/18 and 5/18. */
const std::array<double,3>& gauss_weights ();

/**
 * The towers' member temperatures: created by Conductors when conductors.thermal.towers = true. It holds every
 * heated member's temperature, steps the balance, rebuilds the frames, writes
 * <diagnostics_dir>/tower_temperature.dat (per tower and output step: the members' mean and largest
 * temperature, the member id of the largest, the largest member irradiance, the least k_E of the frame, the
 * frame's rebuilds so far and its largest member temperature) and <diagnostics_dir>/tower_<tower>_member_temperature.dat
 * (every member's temperature), and its checkpoint <chk>/conductors/tower_thermal.
 */
class TowerThermal
{
public:
    /** Read conductors.thermal.towers (false by default; anything not a bool aborts). */
    static bool selected ();

    /**
     * Read conductors.thermal.tower_* and abort, naming the key, on a value out of range. names holds the heated
     * towers' names and sections each tower's member sections (frame order).
     */
    TowerThermal (const std::vector<std::string>& names, const std::vector<std::vector<MemberSection>>& sections,
                  const std::string& diagnostics_dir);
    /** For tests: the settings given directly, checked as the constructor checks them. */
    TowerThermal (const TowerThermalSettings& s, const std::vector<std::string>& names,
                  const std::vector<std::vector<MemberSection>>& sections, const std::string& diagnostics_dir);

    const TowerThermalSettings& settings () const { return m_s; }
    std::size_t num_towers () const { return m_names.size(); }

    /**
     * Step every member's temperature over [time, time + dt] and rebuild a tower's frame when its members' reduction
     * factors have moved more than the tolerance. On the first call of a fresh start each member starts at its air
     * temperature then (the steady state without the fire), which also becomes its surroundings' temperature.
     * Aborts, naming the tower and member, when a member reaches 1200 C (EN 1993-1-2's steel keeps no strength).
     *
     * @param[in,out] frames  per tower, its frame model (set_temperature() is called on it)
     * @param[in] T_air       per tower, per member: the air temperature along it [C] (Gauss mean)
     * @param[in] E           per tower, per member: the fire's irradiance [W/m^2] (Gauss mean; zeros without a fire)
     * @param[in] write_row   write the logs' rows for this step (the conductors' diagnostics cadence)
     */
    void advance (double time, double dt, const std::vector<erf_towers::FrameTower*>& frames,
                  const std::vector<std::vector<double>>& T_air, const std::vector<std::vector<double>>& E, bool write_row);

    /** Member m's temperature of tower t [C]. */
    double temperature (std::size_t t, std::size_t m) const { return m_T[t][m]; }
    const std::vector<double>& temperatures (std::size_t t) const { return m_T[t]; }
    /** The temperature of the surroundings member m of tower t radiates to [C]: its air temperature at the first step. */
    double surroundings_temperature (std::size_t t, std::size_t m) const { return m_Ts[t][m]; }
    /** How many times tower t's frame has been rebuilt at its members' temperatures. */
    long rebuilds (std::size_t t) const { return m_rebuilds[t]; }
    bool started () const { return m_started; }

    /** Save the member and surroundings temperatures and the rebuild counts as <dir>/tower_thermal (I/O rank). */
    void write_checkpoint (const std::string& dir) const;
    /**
     * Read them back on a restart. The frames' temperatures come back with the towers' own state (tower_motion).
     * Aborts when the checkpoint's towers or members do not match these. A checkpoint without the file (written
     * without the towers' heating) starts the temperatures afresh at the first step.
     */
    void restore (const std::string& dir);
    /** The logs this class writes, for Conductors::log_files() (trimmed to the restart time on a restart). */
    std::vector<std::string> log_files () const;

private:
    void write_logs (double time, const std::vector<erf_towers::FrameTower*>& frames, bool first) const;

    TowerThermalSettings m_s;
    std::vector<std::string> m_names;
    std::vector<std::vector<MemberSection>> m_sec;
    std::string m_dir;
    bool m_started{false};
    bool m_restored{false};
    std::vector<std::vector<double>> m_T;     //!< per tower, per member [C]
    std::vector<std::vector<double>> m_Ts;    //!< per tower, per member: the surroundings' temperature [C]
    std::vector<std::vector<double>> m_E;     //!< per tower, per member: the irradiance of the last step [W/m^2]
    std::vector<long> m_rebuilds;             //!< per tower
};

} // namespace erf_conductors

#endif
