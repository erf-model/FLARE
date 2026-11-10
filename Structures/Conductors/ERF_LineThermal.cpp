/**
 * @file ERF_LineThermal.cpp
 * @brief The conductors' temperature under a fire (IEEE Std 738-2023 heat balance) and the sag it adds
 */
#include "ERF_LineThermal.H"
#include "ERF_ConductorLine.H"

#include <ERF_DiagnosticsLog.H>
#include <ERF_StructureModel.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

using namespace amrex;

namespace erf_conductors {

namespace {

constexpr double pi = 3.14159265358979323846;
/// No conductor survives this temperature [C]: a solve that does not stop below it has gone wrong.
constexpr double T_ceiling = 3000.0;
/// The lowest air or conductor temperature the solver looks at [C].
constexpr double T_floor = -150.0;

double sign_of (double x) { return (x > 0.0) ? 1.0 : ((x < 0.0) ? -1.0 : 0.0); }

/**
 * The root of the increasing function g on [lo, hi] (g(lo) <= 0 <= g(hi)): bisection until the bracket is
 * 1e-12 of the temperature (at most 200 halvings; about 60 from a bracket of 1000 K).
 */
template <class G>
double bracketed_root (const G& g, double lo, double hi)
{
    for (int it = 0; it < 200; ++it) {
        const double m = 0.5 * (lo + hi);
        if (hi - lo <= 1.0e-12 * std::max(1.0, std::abs(m))) { return m; }
        const double gm = g(m);
        if (gm == 0.0) { return m; }
        if (gm < 0.0) { lo = m; } else { hi = m; }
    }
    return 0.5 * (lo + hi);
}

/**
 * The root of g (increasing in practice) near T0: bracket it by doubling steps away from T0, then solve.
 * what names the solve in the abort when no root lies in [T_floor, T_ceiling].
 */
template <class G>
double solve_temperature (const G& g, double T0, const std::string& what)
{
    const double g0 = g(T0);
    if (g0 == 0.0) { return T0; }
    double step = 1.0;
    if (g0 < 0.0) {
        double lo = T0;
        while (true) {
            const double hi = std::min(T0 + step, T_ceiling);
            const double ghi = g(hi);
            if (!std::isfinite(ghi)) { Abort(what + ": the heat balance is not finite at " + std::to_string(hi) + " C"); }
            if (ghi >= 0.0) { return bracketed_root(g, lo, hi); }
            if (hi >= T_ceiling) {
                Abort(what + ": the conductor would pass " + std::to_string(T_ceiling) + " C; the heat balance has no "
                      "solution below it (check heat_capacity, the current and the fire's radiative_fraction)");
            }
            lo = hi; step *= 2.0;
        }
    }
    double hi = T0;
    while (true) {
        const double lo = std::max(T0 - step, T_floor);
        const double glo = g(lo);
        if (!std::isfinite(glo)) { Abort(what + ": the heat balance is not finite at " + std::to_string(lo) + " C"); }
        if (glo <= 0.0) { return bracketed_root(g, lo, hi); }
        if (lo <= T_floor) { Abort(what + ": no solution of the heat balance above " + std::to_string(T_floor) + " C"); }
        hi = lo; step *= 2.0;
    }
}

bool in_range (double v, double lo, double hi) { return std::isfinite(v) && v >= lo && v <= hi; }

} // namespace

AirFilm ieee738_air (double T_film, double elevation)
{
    AirFilm a;
    a.mu  = 1.458e-6 * std::pow(T_film + 273.0, 1.5) / (T_film + 383.4);                                  // eq 13
    a.rho = (1.293 - 1.525e-4 * elevation + 6.379e-9 * elevation * elevation) / (1.0 + 0.00367 * T_film); // eq 15
    a.k   = 2.424e-2 + 7.477e-5 * T_film - 4.407e-9 * T_film * T_film;                                    // eq 17
    return a;
}

double ieee738_wind_angle_factor (double phi)
{
    return 1.194 - std::cos(phi) + 0.194 * std::cos(2.0 * phi) + 0.368 * std::sin(2.0 * phi);           // eq 19
}

double ieee738_forced_convection (double D, double V, double phi, double T, double T_a, double elevation)
{
    const AirFilm a = ieee738_air(0.5 * (T + T_a), elevation);
    const double Re = D * a.rho * V / a.mu;                                                                   // eq 11
    const double K = ieee738_wind_angle_factor(phi);
    const double q1 = K * (1.01 + 1.35 * std::pow(Re, 0.52)) * a.k * (T - T_a);                               // eq 9
    const double q2 = K * 0.754 * std::pow(Re, 0.6) * a.k * (T - T_a);                                        // eq 10
    return (std::abs(q1) >= std::abs(q2)) ? q1 : q2;
}

double ieee738_natural_convection (double D, double T, double T_a, double elevation)
{
    const AirFilm a = ieee738_air(0.5 * (T + T_a), elevation);
    return sign_of(T - T_a) * 3.645 * std::sqrt(a.rho) * std::pow(D, 0.75) * std::pow(std::abs(T - T_a), 1.25); // eq 7
}

double ieee738_convection (double D, double V, double phi, double T, double T_a, double elevation)
{
    const double qf = ieee738_forced_convection(D, V, phi, T, T_a, elevation);
    const double qn = ieee738_natural_convection(D, T, T_a, elevation);
    return (std::abs(qf) >= std::abs(qn)) ? qf : qn;
}

double ieee738_radiation (double D, double emissivity, double T, double T_s)
{
    const double a = (T + 273.0) / 100.0, b = (T_s + 273.0) / 100.0;
    return 17.8 * D * emissivity * (a * a * a * a - b * b * b * b);                                            // eq 21
}

double LineThermalProps::resistance_at (double T) const
{
    const double t0 = resistance_temperatures[0], t1 = resistance_temperatures[1];
    return resistance[0] + (resistance[1] - resistance[0]) / (t1 - t0) * (T - t0);                            // eq 33
}

double heat_gain (const LineThermalProps& p, const ThermalSettings& s, double T, double T_a, double T_s, double V, double phi,
                  double E_fire)
{
    const double qJ = (p.current > 0.0) ? p.current * p.current * p.resistance_at(T) : 0.0;
    const double qs = p.absorptivity * p.diameter * s.solar_flux;
    const double qf = p.absorptivity * p.diameter * E_fire;
    return qJ + qs + qf - ieee738_convection(p.diameter, V, phi, T, T_a, s.elevation)
           - ieee738_radiation(p.diameter, p.emissivity, T, T_s);
}

double steady_temperature (const LineThermalProps& p, const ThermalSettings& s, double T_a, double V, double phi, double E_fire,
                           const std::string& where)
{
    // g = -heat_gain has one sign change: the losses grow faster than the Joule heating above T_a
    return solve_temperature([&] (double T) { return -heat_gain(p, s, T, T_a, T_a, V, phi, E_fire); }, T_a,
                             "conductors.thermal_model" + (where.empty() ? std::string() : " (" + where + ")") +
                             ": the steady temperature");
}

double implicit_step (const LineThermalProps& p, const ThermalSettings& s, double T_old, double dt,
                      double T_a, double T_s, double V, double phi, double E_fire, const std::string& where)
{
    return solve_temperature([&] (double T) { return p.heat_capacity * (T - T_old) - dt * heat_gain(p, s, T, T_a, T_s, V, phi, E_fire); },
                             T_old, "conductors.thermal_model" + (where.empty() ? std::string() : " (" + where + ")") +
                             ": a temperature step");
}

std::vector<FireSource> gather_fire_sources (const StructureFire& fire, const std::array<double,3>& lo,
                                             const std::array<double,3>& hi, double reach, const std::array<double,2>& period)
{
    std::vector<double> mine;   // x, y, Q per source
    if (fire.heat_release != nullptr) {
        const MultiFab& q = *fire.heat_release;
        const auto plo = fire.geom.ProbLoArray();
        const auto dx = fire.geom.CellSizeArray();
        const double dA = static_cast<double>(dx[0]) * static_cast<double>(dx[1]);
        for (MFIter mfi(q); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            // the box on the host (ghost cells included, component 0): the fire's fields live in device memory on
            // a GPU build; dtoh_memcpy waits for the copy
            const FArrayBox& src = q[mfi];
            FArrayBox h(src.box(), 1, The_Pinned_Arena());
            Gpu::dtoh_memcpy(h.dataPtr(), src.dataPtr(0), h.nBytes());
            const auto a = h.const_array();
            const auto lo3 = lbound(bx), hi3 = ubound(bx);
            for (int k = lo3.z; k <= hi3.z; ++k) {
                for (int j = lo3.y; j <= hi3.y; ++j) {
                    for (int i = lo3.x; i <= hi3.x; ++i) {
                        const double v = static_cast<double>(a(i, j, k));
                        if (!(v > 0.0)) { continue; }
                        const double x = static_cast<double>(plo[0]) + (i + 0.5) * static_cast<double>(dx[0]);
                        const double y = static_cast<double>(plo[1]) + (j + 0.5) * static_cast<double>(dx[1]);
                        // the distance wraps in a periodic direction: no cell is left out there
                        const bool far_x = !(period[0] > 0.0) && (x < lo[0] - reach || x > hi[0] + reach);
                        const bool far_y = !(period[1] > 0.0) && (y < lo[1] - reach || y > hi[1] + reach);
                        if (reach >= 0.0 && (far_x || far_y)) { continue; }
                        mine.insert(mine.end(), {x, y, v * dA});
                    }
                }
            }
        }
    }
    // every rank gets every rank's sources
    const int nprocs = ParallelDescriptor::NProcs();
    if (mine.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        Abort("conductors.thermal_model: too many burning fire cells on one rank; set conductors.thermal.cutoff_radius smaller");
    }
    const int n = static_cast<int>(mine.size());
    std::vector<int> counts(static_cast<std::size_t>(nprocs), 0), disp(static_cast<std::size_t>(nprocs), 0);
    ParallelAllGather::AllGather(n, counts.data(), ParallelDescriptor::Communicator());
    Long sum = 0;
    for (int r = 0; r < nprocs; ++r) { sum += counts[static_cast<std::size_t>(r)]; }
    if (sum > static_cast<Long>(std::numeric_limits<int>::max())) {
        Abort("conductors.thermal_model: too many burning fire cells to gather; set conductors.thermal.cutoff_radius smaller");
    }
    int total = 0;
    for (int r = 0; r < nprocs; ++r) { disp[static_cast<std::size_t>(r)] = total; total += counts[static_cast<std::size_t>(r)]; }
    std::vector<double> all(static_cast<std::size_t>(total));
    const int root = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::Gatherv(mine.data(), n, all.data(), counts, disp, root);
    if (total > 0) { ParallelDescriptor::Bcast(all.data(), static_cast<std::size_t>(total), root); }
    std::vector<FireSource> src(static_cast<std::size_t>(total / 3));
    for (std::size_t j = 0; j < src.size(); ++j) { src[j] = FireSource{all[3*j], all[3*j+1], 0.0, all[3*j+2]}; }
    // one order whatever the boxes and ranks, so that the sums do not depend on the decomposition
    std::sort(src.begin(), src.end(), [] (const FireSource& a, const FireSource& b) {
        if (a.y != b.y) { return a.y < b.y; }
        if (a.x != b.x) { return a.x < b.x; }
        return a.Q < b.Q;
    });
    return src;
}

std::vector<double> fire_irradiance (const std::vector<double>& nodes, const std::vector<double>& tangents,
                                     const std::vector<FireSource>& sources, double radiative_fraction,
                                     double cutoff, double r_min2, const std::array<double,2>& period)
{
    // the nearest image's offset in a periodic direction
    auto image = [] (double d, double L) { return (L > 0.0) ? d - L * std::round(d / L) : d; };
    const std::size_t n = nodes.size() / 3;
    std::vector<double> E(n, 0.0);
    const double cut2 = (cutoff > 0.0) ? cutoff * cutoff : std::numeric_limits<double>::infinity();
    for (std::size_t m = 0; m < n; ++m) {
        const double tx = tangents[3*m], ty = tangents[3*m+1], tz = tangents[3*m+2];
        double sum = 0.0;
        for (const FireSource& s : sources) {
            const double rx = image(nodes[3*m] - s.x, period[0]), ry = image(nodes[3*m+1] - s.y, period[1]), rz = nodes[3*m+2] - s.z;
            const double r2 = rx * rx + ry * ry + rz * rz;
            if (r2 > cut2) { continue; }
            // sin of the angle between the ray and the axis: |t x r| / |r|
            double sin_psi = 1.0;
            if (r2 > 0.0) {
                const double cx = ty * rz - tz * ry, cy = tz * rx - tx * rz, cz = tx * ry - ty * rx;
                sin_psi = std::min(1.0, std::sqrt((cx * cx + cy * cy + cz * cz) / r2));
            }
            sum += s.Q * sin_psi / (4.0 * pi * std::max(r2, r_min2));
        }
        E[m] = radiative_fraction * sum;
    }
    return E;
}

std::vector<double> span_tangents (const std::vector<double>& nodes)
{
    const std::size_t n = nodes.size() / 3;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n >= 2, "span_tangents: a span has at least two nodes");
    std::vector<double> t(3 * n);
    for (std::size_t m = 0; m < n; ++m) {
        const std::size_t a = (m == 0) ? 0 : m - 1, b = (m + 1 == n) ? n - 1 : m + 1;
        double d[3], len2 = 0.0;
        for (int c = 0; c < 3; ++c) { d[c] = nodes[3*b+static_cast<std::size_t>(c)] - nodes[3*a+static_cast<std::size_t>(c)]; len2 += d[c] * d[c]; }
        const double len = std::sqrt(len2);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(len > 0.0, "span_tangents: two neighbouring nodes coincide");
        for (int c = 0; c < 3; ++c) { t[3*m+static_cast<std::size_t>(c)] = d[c] / len; }
    }
    return t;
}

double span_mean (const std::vector<double>& v)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(v.size() >= 2, "span_mean: a span has at least two nodes");
    double s = 0.5 * (v.front() + v.back());
    for (std::size_t m = 1; m + 1 < v.size(); ++m) { s += v[m]; }
    return s / static_cast<double>(v.size() - 1);
}

bool LineThermal::selected ()
{
    std::string model = "none";
    ParmParse("conductors").query("thermal_model", model);
    model = amrex::toLower(model);
    if (model != "none" && model != "ieee738") {
        Abort("conductors.thermal_model must be none or ieee738, got '" + model + "'");
    }
    return model == "ieee738";
}

std::string LineThermal::check (const ThermalSettings& s, const std::string& line, const LineThermalProps& p)
{
    const std::string t = "conductors.thermal.";
    if (s.radiative_fraction_given && !(std::isfinite(s.radiative_fraction) && s.radiative_fraction > 0.0 && s.radiative_fraction <= 1.0)) {
        return t + "radiative_fraction must be in (0, 1] (the fraction of the fire's heat release it radiates)";
    }
    if (!in_range(s.source_height, 0.0, 1.0e3)) { return t + "source_height must be in [0, 1000] m (above the ground)"; }
    if (!in_range(s.cutoff_radius, 0.0, 1.0e6)) { return t + "cutoff_radius must be in [0, 1e6] m (needed, > 0, when the run has a fire)"; }
    if (!in_range(s.elevation, -500.0, 9000.0)) { return t + "elevation must be in [-500, 9000] m above sea level"; }
    if (!in_range(s.solar_flux, 0.0, 1500.0)) { return t + "solar_flux must be in [0, 1500] W/m^2"; }
    if (s.reference_temperature_given && !in_range(s.reference_temperature, -50.0, 300.0)) {
        return t + "reference_temperature must be in [-50, 300] C (the temperature of the spans' input lengths)";
    }
    if (line.empty()) { return std::string(); }
    const std::string k = "conductors." + line + ".";
    if (!(std::isfinite(p.diameter) && p.diameter > 0.0)) { return k + "diameter must be positive (m)"; }
    if (!in_range(p.absorptivity, 0.0, 1.0)) { return k + "absorptivity must be in [0, 1] (needed with conductors.thermal_model = ieee738)"; }
    if (!in_range(p.emissivity, 0.0, 1.0)) { return k + "emissivity must be in [0, 1] (needed with conductors.thermal_model = ieee738)"; }
    if (!(std::isfinite(p.heat_capacity) && p.heat_capacity > 0.0)) {
        return k + "heat_capacity must be positive (J/(m K), the conductor's m c_p per unit length; needed with conductors.thermal_model = ieee738)";
    }
    if (!in_range(p.thermal_expansion, 0.0, 1.0e-4)) {
        return k + "thermal_expansion must be in [0, 1e-4] 1/K (needed with conductors.thermal_model = ieee738)";
    }
    if (!in_range(p.current, 0.0, 1.0e5)) { return k + "current must be in [0, 1e5] A"; }
    const auto& rt = p.resistance_temperatures;
    if (!(in_range(rt[0], -50.0, 500.0) && in_range(rt[1], -50.0, 500.0) && rt[1] > rt[0])) {
        return k + "resistance_temperatures must be two increasing values in [-50, 500] C";
    }
    if (p.current > 0.0) {
        if (!(std::isfinite(p.resistance[0]) && p.resistance[0] > 0.0 && std::isfinite(p.resistance[1]) && p.resistance[1] > 0.0)) {
            return k + "resistance must be two positive values (ohm/m at the two resistance_temperatures; needed with a current)";
        }
        // the line through them stays positive over every temperature the solver looks at
        if (!(p.resistance_at(T_floor) > 0.0 && p.resistance_at(T_ceiling) > 0.0)) {
            return k + "resistance: the line through the two values must stay positive from " + std::to_string(static_cast<int>(T_floor)) +
                   " to " + std::to_string(static_cast<int>(T_ceiling)) + " C";
        }
    }
    return std::string();
}

LineThermal::LineThermal (const ThermalSettings& s, const std::vector<std::string>& line_names,
                          const std::vector<LineThermalProps>& props, const std::string& diagnostics_dir)
    : m_s(s), m_names(line_names), m_p(props), m_dir(diagnostics_dir)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_names.size() == m_p.size(), "LineThermal: one set of thermal data per line");
    for (const auto& n : m_names) {
        if (n == "thermal") {
            Abort("conductors.lines: a line may not be named 'thermal' with conductors.thermal_model = ieee738 "
                  "(its keys would be the shared conductors.thermal.* settings)");
        }
    }
    std::string err = check(m_s, std::string(), LineThermalProps{});
    for (std::size_t i = 0; err.empty() && i < m_p.size(); ++i) { err = check(m_s, m_names[i], m_p[i]); }
    if (!err.empty()) { Abort(err); }
}

namespace {

ThermalSettings read_settings ()
{
    ThermalSettings s;
    ParmParse pt("conductors.thermal");
    s.radiative_fraction_given = pt.query("radiative_fraction", s.radiative_fraction) != 0;
    pt.query("source_height", s.source_height);
    pt.query("cutoff_radius", s.cutoff_radius);
    pt.query("elevation", s.elevation);
    pt.query("solar_flux", s.solar_flux);
    s.reference_temperature_given = pt.query("reference_temperature", s.reference_temperature) != 0;
    return s;
}

std::vector<LineThermalProps> read_props (const std::vector<std::string>& names, const std::vector<double>& diameters)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(names.size() == diameters.size(), "LineThermal: one diameter per line");
    std::vector<LineThermalProps> v(names.size());
    for (std::size_t i = 0; i < names.size(); ++i) {
        ParmParse pl("conductors." + names[i]);
        LineThermalProps& p = v[i];
        p.diameter = diameters[i];
        const std::string k = "conductors." + names[i] + ".";
        const char* required[] = {"absorptivity", "emissivity", "heat_capacity", "thermal_expansion"};
        double* dst[] = {&p.absorptivity, &p.emissivity, &p.heat_capacity, &p.thermal_expansion};
        for (int r = 0; r < 4; ++r) {
            if (pl.query(required[r], *dst[r]) == 0) {
                Abort(k + required[r] + " is needed with conductors.thermal_model = ieee738");
            }
        }
        pl.query("current", p.current);
        std::vector<double> rv;
        const bool r_given = pl.queryarr("resistance", rv) != 0;
        if (r_given) {
            if (rv.size() != 2) { Abort(k + "resistance needs two values (ohm/m at the two resistance_temperatures)"); }
            p.resistance = {{rv[0], rv[1]}};
        }
        std::vector<double> tv;
        if (pl.queryarr("resistance_temperatures", tv)) {
            if (tv.size() != 2) { Abort(k + "resistance_temperatures needs two values (C)"); }
            p.resistance_temperatures = {{tv[0], tv[1]}};
        }
        if (r_given && !(p.current > 0.0)) {
            Print() << "conductors: " << k << "resistance is not used: the line carries no current\n";
        }
    }
    return v;
}

} // namespace

LineThermal::LineThermal (const std::vector<std::string>& line_names, const std::vector<double>& diameters,
                          const std::string& diagnostics_dir)
    : LineThermal(read_settings(), line_names, read_props(line_names, diameters), diagnostics_dir)
{
    Print() << "conductors: thermal model IEEE 738 on " << m_names.size() << " line(s)"
            << (m_s.radiative_fraction_given ? ", fire radiation fraction " + std::to_string(m_s.radiative_fraction) : std::string())
            << ", sources " << m_s.source_height << " m above the ground"
            << (m_s.cutoff_radius > 0.0 ? ", cutoff " + std::to_string(m_s.cutoff_radius) + " m" : std::string())
            << (m_s.reference_temperature_given ? ", input lengths at " + std::to_string(m_s.reference_temperature) + " C"
                                                : std::string(", input lengths at the first step's steady temperature")) << "\n";
}

double LineThermal::span_temperature (std::size_t i, int k) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_started, "LineThermal::span_temperature: before the first step");
    return m_span_T.at(i).at(static_cast<std::size_t>(k));
}

void LineThermal::advance (double time, double dt, std::vector<std::unique_ptr<ConductorLine>>& lines,
                           const std::vector<std::vector<Real>>& wind, const std::vector<std::vector<double>>& T_air,
                           const std::vector<FireSource>& sources, double r_min2, const std::array<double,2>& period, bool write_row)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lines.size() == m_p.size() && wind.size() == lines.size() && T_air.size() == lines.size(),
                                     "LineThermal::advance: one set of data per line");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt > 0.0 && std::isfinite(dt), "LineThermal::advance: dt must be positive and finite");
    if (!sources.empty() && !m_s.radiative_fraction_given) {
        Abort("conductors.thermal.radiative_fraction is needed when the run has a fire (the fraction of its heat release it radiates)");
    }
    const bool fresh = !m_started;
    if (fresh) {
        m_T.assign(lines.size(), {});
        m_Ts.assign(lines.size(), {});
        m_Tref.assign(lines.size(), {});
        m_Lref.assign(lines.size(), {});
    }
    m_E.assign(lines.size(), {});
    m_air.assign(lines.size(), {});
    m_span_T.assign(lines.size(), {});
    for (std::size_t i = 0; i < lines.size(); ++i) {
        ConductorLine& line = *lines[i];
        const unsigned nn = line.num_nodes();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(T_air[i].size() == nn && wind[i].size() >= 3 * static_cast<std::size_t>(nn),
                                         "LineThermal::advance: the air temperature and wind need every node");
        const std::vector<Real> pos = line.node_positions();
        if (fresh) {
            m_T[i].assign(T_air[i].begin(), T_air[i].end());
            // the surroundings each node radiates to: the air around it before the fire
            m_Ts[i].assign(T_air[i].begin(), T_air[i].end());
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_T[i].size() == nn && m_Ts[i].size() == nn, "LineThermal::advance: the line's nodes changed");
        m_E[i].assign(nn, 0.0);
        m_air[i].assign(nn, {{0.0, 0.0, 0.0}});
        m_span_T[i].assign(static_cast<std::size_t>(line.num_spans()), 0.0);
        const LineThermalProps& p = m_p[i];
        for (int k = 0; k < line.num_spans(); ++k) {
            const unsigned n0 = line.span_first_node(k), ns = line.span_num_nodes(k);
            std::vector<double> x(3 * static_cast<std::size_t>(ns));
            for (unsigned m = 0; m < ns; ++m) {
                for (int c = 0; c < 3; ++c) { x[3*m+static_cast<std::size_t>(c)] = static_cast<double>(pos[3*(n0+m)+static_cast<std::size_t>(c)]); }
            }
            const std::vector<double> t = span_tangents(x);
            const std::vector<double> E = sources.empty() ? std::vector<double>(ns, 0.0)
                                        : fire_irradiance(x, t, sources, m_s.radiative_fraction, m_s.cutoff_radius, r_min2, period);
            std::vector<double> Ts(ns), T0(ns);
            for (unsigned m = 0; m < ns; ++m) {
                const unsigned node = n0 + m;
                const double u = static_cast<double>(wind[i][3*node]), v = static_cast<double>(wind[i][3*node+1]),
                             w = static_cast<double>(wind[i][3*node+2]);
                const double V = std::sqrt(u * u + v * v + w * w);
                // the angle between the wind and the axis; a calm takes the wind across the conductor
                double phi = 0.5 * pi;
                if (V > 0.0) {
                    const double c = std::abs(u * t[3*m] + v * t[3*m+1] + w * t[3*m+2]) / V;
                    phi = std::acos(std::min(1.0, c));
                }
                const double Ta = T_air[i][node];
                if (!std::isfinite(Ta) || !std::isfinite(V)) {
                    Abort("conductors." + m_names[i] + ": the air temperature or the wind is not finite at node " + std::to_string(node));
                }
                const std::string where = "line " + m_names[i] + ", node " + std::to_string(node);
                if (fresh) {
                    m_T[i][node] = steady_temperature(p, m_s, Ta, V, phi, 0.0, where);
                    T0[m] = m_T[i][node];
                }
                m_T[i][node] = implicit_step(p, m_s, m_T[i][node], dt, Ta, m_Ts[i][node], V, phi, E[m], where);
                m_E[i][node] = E[m];
                m_air[i][node] = {{Ta, V, phi}};
                Ts[m] = m_T[i][node];
            }
            m_span_T[i][static_cast<std::size_t>(k)] = span_mean(Ts);
            if (fresh) {
                // the reference: the temperature of the input lengths, or the span's mean steady temperature without the fire
                m_Tref[i].push_back(m_s.reference_temperature_given ? m_s.reference_temperature : span_mean(T0));
                m_Lref[i].push_back(line.span_unstretched_length(k));
            }
        }
        // the strings' nodes are not heated: they hold the air's temperature
        for (unsigned node = line.span_first_node(line.num_spans()); node < nn; ++node) { m_T[i][node] = T_air[i][node]; }
    }
    m_started = true;
    set_lengths(lines);
    if (fresh || write_row) { write_log(time + dt, lines, fresh && !m_restored); }
}

void LineThermal::set_lengths (std::vector<std::unique_ptr<ConductorLine>>& lines) const
{
    for (std::size_t i = 0; i < lines.size(); ++i) {
        for (int k = 0; k < lines[i]->num_spans(); ++k) {
            const double Tbar = m_span_T[i][static_cast<std::size_t>(k)];
            const double L = m_Lref[i][static_cast<std::size_t>(k)] * (1.0 + m_p[i].thermal_expansion * (Tbar - m_Tref[i][static_cast<std::size_t>(k)]));
            if (!(std::isfinite(L) && L > 0.0)) {
                Abort("conductors." + m_names[i] + ": span " + std::to_string(k) + "'s thermal length is not positive");
            }
            lines[i]->set_span_unstretched_length(k, L);
        }
    }
}

void LineThermal::write_log (double t, const std::vector<std::unique_ptr<ConductorLine>>& lines, bool truncate) const
{
    if (!ParallelDescriptor::IOProcessor()) { return; }
    std::ofstream out;
    // the header after truncating, or when a restart appends to a missing or empty file
    if (erf_actuator::open_log(out, log_file(), truncate)) {
        out << "# time line span T_mean_C T_max_C node_of_max T_air_C wind_m_s wind_angle_deg E_fire_W_m2 E_fire_max_W_m2 strain T_ref_C\n";
    }
    out << std::setprecision(10);
    for (std::size_t i = 0; i < lines.size(); ++i) {
        for (int k = 0; k < lines[i]->num_spans(); ++k) {
            const unsigned n0 = lines[i]->span_first_node(k), ns = lines[i]->span_num_nodes(k);
            double Tmax = -std::numeric_limits<double>::infinity(), Emax = 0.0;
            unsigned at = n0;
            for (unsigned m = n0; m < n0 + ns; ++m) {
                if (m_T[i][m] > Tmax) { Tmax = m_T[i][m]; at = m; }
                Emax = std::max(Emax, m_E[i][m]);
            }
            const double Tbar = m_span_T[i][static_cast<std::size_t>(k)];
            const double strain = m_p[i].thermal_expansion * (Tbar - m_Tref[i][static_cast<std::size_t>(k)]);
            const auto& air = m_air[i][at];
            out << t << " " << m_names[i] << " " << k << " " << Tbar << " " << Tmax << " " << at << " " << air[0] << " " << air[1]
                << " " << air[2] * 180.0 / pi << " " << m_E[i][at] << " " << Emax << " " << strain << " " << m_Tref[i][static_cast<std::size_t>(k)] << "\n";
        }
    }
}

void LineThermal::write_checkpoint (const std::string& dir) const
{
    if (!m_started || !ParallelDescriptor::IOProcessor()) { return; }
    std::ofstream out(dir + "/thermal", std::ios::trunc);
    if (!out) { Abort("cannot write the conductor checkpoint '" + dir + "/thermal'"); }
    out << std::setprecision(std::numeric_limits<double>::max_digits10) << "thermal_lines " << m_T.size() << "\n";
    for (std::size_t i = 0; i < m_T.size(); ++i) {
        out << m_names[i] << " " << m_T[i].size() << " " << m_Tref[i].size() << "\n";
        for (const auto* v : {&m_T[i], &m_Ts[i], &m_span_T[i], &m_Tref[i], &m_Lref[i]}) {
            for (const double x : *v) { out << x << " "; }
            out << "\n";
        }
    }
    if (!out) { Abort("cannot write the conductor checkpoint '" + dir + "/thermal'"); }
}

void LineThermal::restore (const std::string& dir, std::vector<std::unique_ptr<ConductorLine>>& lines)
{
    m_restored = true;
    if (!erf_actuator::file_exists_everywhere(dir + "/thermal")) {
        Print() << "conductors: the checkpoint '" << dir << "' holds no conductor temperatures; they start at the restart, "
                << "at their steady temperatures without the fire, and the spans' lengths are taken as they are then\n";
        return;
    }
    Vector<char> chars;
    ParallelDescriptor::ReadAndBcastFile(dir + "/thermal", chars);
    // ReadAndBcastFile ends the buffer with a NUL
    std::istringstream in{std::string(chars.dataPtr())};
    std::string tag;
    std::size_t nl = 0;
    if (!(in >> tag >> nl) || tag != "thermal_lines") { Abort("malformed conductor checkpoint '" + dir + "/thermal'"); }
    if (nl != lines.size()) {
        Abort("the conductor checkpoint '" + dir + "/thermal' holds " + std::to_string(nl) + " line(s), this run has " +
              std::to_string(lines.size()) + "; conductors.lines must match the run being restarted");
    }
    m_T.assign(nl, {}); m_Ts.assign(nl, {}); m_span_T.assign(nl, {}); m_Tref.assign(nl, {}); m_Lref.assign(nl, {});
    m_E.assign(nl, {}); m_air.assign(nl, {});
    for (std::size_t i = 0; i < nl; ++i) {
        std::string name;
        std::size_t nn = 0, nk = 0;
        if (!(in >> name >> nn >> nk)) { Abort("malformed conductor checkpoint '" + dir + "/thermal'"); }
        if (name != m_names[i] || nn != lines[i]->num_nodes() || nk != static_cast<std::size_t>(lines[i]->num_spans())) {
            Abort("the conductor checkpoint '" + dir + "/thermal': line " + name + " does not match line " + m_names[i] +
                  " of this run (its name, nodes or spans)");
        }
        auto read_n = [&] (std::vector<double>& v, std::size_t count) {
            v.resize(count);
            for (auto& x : v) { if (!(in >> x) || !std::isfinite(x)) { Abort("malformed conductor checkpoint '" + dir + "/thermal'"); } }
        };
        read_n(m_T[i], nn);
        read_n(m_Ts[i], nn);
        read_n(m_span_T[i], nk);
        read_n(m_Tref[i], nk);
        read_n(m_Lref[i], nk);
        m_E[i].assign(nn, 0.0);
        m_air[i].assign(nn, {{0.0, 0.0, 0.0}});
    }
    m_started = true;
    set_lengths(lines);
}

} // namespace erf_conductors
