/**
 * @file ERF_TowerThermal.cpp
 * @brief The steel temperature of the lattice towers' members under a fire (EN 1993-1-2) and the frame it softens
 */
#include "ERF_TowerThermal.H"

#include "ERF_FrameInputs.H"
#include "ERF_FrameTower.H"
#include "ERF_MemberChecks.H"

#include <ERF_DiagnosticsLog.H>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
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
constexpr double sigma_sb = 5.670374419e-8;   // W/(m^2 K^4)
constexpr double kelvin = 273.15;
/// EN 1993-1-2's steel keeps no strength or stiffness at this temperature [C].
constexpr double steel_limit = 1200.0;
/// The bracket search gives up beyond these [C]: a solve that does not stop inside them has gone wrong.
constexpr double theta_floor = -200.0;
constexpr double theta_ceiling = 5000.0;

bool in_range (double x, double lo, double hi) { return std::isfinite(x) && x >= lo && x <= hi; }

// the antiderivatives of the four pieces of c_a (EN 1993-1-2 section 3.4.1.2)
double F1 (double t) { return 425.0 * t + 0.5 * 7.73e-1 * t * t - 1.69e-3 / 3.0 * t * t * t + 2.22e-6 / 4.0 * t * t * t * t; }
double F2 (double t) { return 666.0 * t - 13002.0 * std::log(738.0 - t); }
double F3 (double t) { return 545.0 * t + 17820.0 * std::log(t - 731.0); }

/**
 * The root of the increasing function f, bracketed by doubling steps from x0 (upwards when f(x0) < 0, else
 * downwards), then bisection until the bracket is narrower than 1e-10. Aborts, naming where, when no bracket
 * is found in [theta_floor, theta_ceiling].
 */
template <typename F>
double increasing_root (F&& f, double x0, const std::string& where)
{
    const double f0 = f(x0);
    if (f0 == 0.0) { return x0; }
    double lo = x0, hi = x0;
    double step = 1.0;
    if (f0 < 0.0) {
        for (;;) {
            // the trial point stops at the ceiling; only a ceiling where f is still negative has no root below it
            hi = std::min(x0 + step, theta_ceiling);
            if (f(hi) >= 0.0) { break; }
            if (hi >= theta_ceiling) { Abort(where + ": no member temperature balances the heat below " + std::to_string(theta_ceiling) + " C"); }
            lo = hi;
            step *= 2.0;
        }
    } else {
        for (;;) {
            lo = std::max(x0 - step, theta_floor);
            if (f(lo) <= 0.0) { break; }
            if (lo <= theta_floor) { Abort(where + ": no member temperature balances the heat above " + std::to_string(theta_floor) + " C"); }
            hi = lo;
            step *= 2.0;
        }
    }
    while (hi - lo > 1.0e-10) {
        const double mid = 0.5 * (lo + hi);
        if (mid <= lo || mid >= hi) { break; }
        if (f(mid) < 0.0) { lo = mid; } else { hi = mid; }
    }
    return 0.5 * (lo + hi);
}

} // namespace

double steel_specific_heat (double theta)
{
    const double t = std::max(theta, 20.0);
    if (t < 600.0) { return 425.0 + 7.73e-1 * t - 1.69e-3 * t * t + 2.22e-6 * t * t * t; }
    if (t < 735.0) { return 666.0 + 13002.0 / (738.0 - t); }
    if (t < 900.0) { return 545.0 + 17820.0 / (t - 731.0); }
    return 650.0;
}

double steel_enthalpy (double theta)
{
    if (theta < 20.0) { return steel_specific_heat(20.0) * (theta - 20.0); }
    const double h600 = F1(600.0) - F1(20.0);
    if (theta < 600.0) { return F1(theta) - F1(20.0); }
    const double h735 = h600 + F2(735.0) - F2(600.0);
    if (theta < 735.0) { return h600 + F2(theta) - F2(600.0); }
    const double h900 = h735 + F3(900.0) - F3(735.0);
    if (theta < 900.0) { return h735 + F3(theta) - F3(735.0); }
    return h900 + 650.0 * (theta - 900.0);
}

bool member_section (const erf_towers::FrameInputs& in, const std::vector<erf_towers::MemberDesign>* designs, std::size_t m,
                     MemberSection& sec, std::string& err)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m < in.members.size(), "member_section: no such member");
    const erf_towers::FrameMember& mem = in.members[m];
    sec = MemberSection();
    if (designs != nullptr && m < designs->size() && (*designs)[m].b > 0.0 && (*designs)[m].t > 0.0) {
        const double b = (*designs)[m].b, t = (*designs)[m].t;
        sec.area = t * (2.0 * b - t);
        sec.perimeter = 4.0 * b;
        // the hull of the angle's outline: the two outer faces b, the toe edges t and the chord between the toes
        sec.width = ((2.0 + std::sqrt(2.0)) * b + (2.0 - std::sqrt(2.0)) * t) / pi;
        return true;
    }
    const erf_towers::FrameSection* s = in.section(mem.section, mem.shape);
    if (s == nullptr) { err = "member " + std::to_string(mem.id) + " has no section"; return false; }
    if (s->shape == erf_towers::SectionShape::Circular) {
        const double D = s->D, t = s->t;
        const double di = (t > 0.0) ? std::max(D - 2.0 * t, 0.0) : 0.0;
        sec.area = 0.25 * pi * (D * D - di * di);
        sec.perimeter = pi * D;
        sec.width = D;
    } else if (s->shape == erf_towers::SectionShape::Rectangular) {
        const double a = s->Sa, b = s->Sb, t = s->t;
        const double ai = (t > 0.0) ? std::max(a - 2.0 * t, 0.0) : 0.0, bi = (t > 0.0) ? std::max(b - 2.0 * t, 0.0) : 0.0;
        sec.area = a * b - ai * bi;
        sec.perimeter = 2.0 * (a + b);
        sec.width = sec.perimeter / pi;
    } else {
        err = "member " + std::to_string(mem.id) + " has an arbitrary section, whose perimeter the frame file does not give; "
              "give its design data as an equal-leg angle (member_file, b and t), or a circular or rectangular section";
        return false;
    }
    if (!(std::isfinite(sec.area) && sec.area > 0.0 && sec.perimeter > 0.0)) {
        err = "member " + std::to_string(mem.id) + " has a section without area or perimeter";
        return false;
    }
    return true;
}

std::string check (const TowerThermalSettings& s)
{
    const std::string t = "conductors.thermal.";
    if (!(std::isfinite(s.emissivity) && s.emissivity > 0.0 && s.emissivity <= 1.0)) {
        return t + "tower_emissivity must be in (0, 1] (the steel's emissivity, also its absorptivity)";
    }
    if (!in_range(s.convection, 0.0, 200.0)) { return t + "tower_convection must be in [0, 200] W/(m^2 K)"; }
    if (!(std::isfinite(s.reduction_tolerance) && s.reduction_tolerance > 0.0 && s.reduction_tolerance <= 0.2)) {
        return t + "tower_reduction_tolerance must be in (0, 0.2] (the relative change of k_y or k_E that rebuilds a frame)";
    }
    return std::string();
}

double member_heat_gain (const MemberSection& sec, const TowerThermalSettings& s, double theta, double T_a, double T_s, double E)
{
    const double th = theta + kelvin, ts = T_s + kelvin;
    const double exchange = s.convection * (T_a - theta) + s.emissivity * sigma_sb * (ts * ts * ts * ts - th * th * th * th);
    return sec.perimeter * exchange + s.emissivity * sec.width * E;
}

double member_implicit_step (const MemberSection& sec, const TowerThermalSettings& s, double theta_old, double dt,
                             double T_a, double T_s, double E, const std::string& where)
{
    const double mass = steel_density * sec.area;
    const double h0 = steel_enthalpy(theta_old);
    auto f = [&] (double th) { return mass * (steel_enthalpy(th) - h0) - dt * member_heat_gain(sec, s, th, T_a, T_s, E); };
    return increasing_root(f, theta_old, where.empty() ? std::string("conductors.thermal.towers: a member's backward-Euler step") : where);
}

double member_steady_temperature (const MemberSection& sec, const TowerThermalSettings& s, double T_a, double T_s, double E,
                                  const std::string& where)
{
    auto f = [&] (double th) { return -member_heat_gain(sec, s, th, T_a, T_s, E); };
    return increasing_root(f, T_a, where.empty() ? std::string("conductors.thermal.towers: a member's steady temperature") : where);
}

double reduction_change (const std::vector<double>& a, const std::vector<double>& b)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a.size() == b.size(), "reduction_change: one temperature per member in both");
    // the change relative to the larger of the two factors (0 when both are 0)
    auto rel = [] (double x, double y) { const double big = std::max(x, y); return big > 0.0 ? std::abs(x - y) / big : 0.0; };
    double d = 0.0;
    for (std::size_t m = 0; m < a.size(); ++m) {
        double ya = 0.0, ea = 0.0, yb = 0.0, eb = 0.0;
        erf_towers::steel_reduction(a[m], ya, ea);
        erf_towers::steel_reduction(b[m], yb, eb);
        d = std::max({d, rel(ya, yb), rel(ea, eb)});
    }
    return d;
}

std::array<std::array<double,3>,3> gauss_points (const std::array<double,3>& x_a, const std::array<double,3>& x_b)
{
    const double h = 0.5 * std::sqrt(0.6);
    const std::array<double,3> s{{0.5 - h, 0.5, 0.5 + h}};
    std::array<std::array<double,3>,3> p{};
    for (std::size_t k = 0; k < 3; ++k) {
        for (std::size_t d = 0; d < 3; ++d) { p[k][d] = x_a[d] + s[k] * (x_b[d] - x_a[d]); }
    }
    return p;
}

const std::array<double,3>& gauss_weights ()
{
    static const std::array<double,3> w{{5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0}};
    return w;
}

bool TowerThermal::selected ()
{
    bool on = false;
    ParmParse("conductors.thermal").query("towers", on);
    return on;
}

TowerThermal::TowerThermal (const TowerThermalSettings& s, const std::vector<std::string>& names,
                            const std::vector<std::vector<MemberSection>>& sections, const std::string& diagnostics_dir)
    : m_s(s), m_names(names), m_sec(sections), m_dir(diagnostics_dir), m_rebuilds(names.size(), 0)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_names.size() == m_sec.size(), "TowerThermal: one set of member sections per tower");
    const std::string err = check(m_s);
    if (!err.empty()) { Abort(err); }
}

namespace {

TowerThermalSettings read_tower_settings ()
{
    TowerThermalSettings s;
    ParmParse pt("conductors.thermal");
    pt.query("tower_emissivity", s.emissivity);
    pt.query("tower_convection", s.convection);
    pt.query("tower_reduction_tolerance", s.reduction_tolerance);
    return s;
}

} // namespace

TowerThermal::TowerThermal (const std::vector<std::string>& names, const std::vector<std::vector<MemberSection>>& sections,
                            const std::string& diagnostics_dir)
    : TowerThermal(read_tower_settings(), names, sections, diagnostics_dir)
{
    std::size_t members = 0;
    for (const auto& v : m_sec) { members += v.size(); }
    Print() << "conductors: tower steel heated (EN 1993-1-2) on " << m_names.size() << " frame tower(s), " << members
            << " member(s); emissivity " << m_s.emissivity << ", convection " << m_s.convection
            << " W/(m^2 K), frames rebuilt when k_y or k_E moves by " << m_s.reduction_tolerance << "\n";
}

void TowerThermal::advance (double time, double dt, const std::vector<erf_towers::FrameTower*>& frames,
                            const std::vector<std::vector<double>>& T_air, const std::vector<std::vector<double>>& E, bool write_row)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(frames.size() == m_names.size() && T_air.size() == m_names.size() && E.size() == m_names.size(),
                                     "TowerThermal::advance: one set of data per tower");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt > 0.0 && std::isfinite(dt), "TowerThermal::advance: dt must be positive and finite");
    const bool fresh = !m_started;
    if (fresh) {
        m_T.assign(m_names.size(), {});
        m_Ts.assign(m_names.size(), {});
    }
    m_E.assign(m_names.size(), {});
    for (std::size_t t = 0; t < m_names.size(); ++t) {
        const std::size_t nm = m_sec[t].size();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(frames[t] != nullptr && T_air[t].size() == nm && E[t].size() == nm &&
                                         frames[t]->temperature().size() == nm,
                                         "TowerThermal::advance: a frame, an air temperature and an irradiance per member are needed");
        if (fresh) {
            // each member starts at its air temperature, the steady state without the fire, which its surroundings keep
            m_T[t] = T_air[t];
            m_Ts[t] = T_air[t];
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_T[t].size() == nm, "TowerThermal::advance: the tower's members changed");
        const auto& ids = frames[t]->frame().inputs().members;
        for (std::size_t m = 0; m < nm; ++m) {
            if (!std::isfinite(T_air[t][m]) || !std::isfinite(E[t][m])) {
                Abort("conductors.thermal.towers: tower " + m_names[t] + ", member " + std::to_string(ids[m].id) +
                      ": the air temperature or the fire's irradiance is not finite");
            }
            m_T[t][m] = member_implicit_step(m_sec[t][m], m_s, m_T[t][m], dt, T_air[t][m], m_Ts[t][m], E[t][m],
                                             "conductors.thermal.towers: tower " + m_names[t] + ", member " + std::to_string(ids[m].id));
            if (!(m_T[t][m] < steel_limit)) {
                Abort("conductors.thermal.towers: tower " + m_names[t] + ", member " + std::to_string(ids[m].id) + " reached " +
                      std::to_string(m_T[t][m]) + " C at t = " + std::to_string(time + dt) +
                      " s; EN 1993-1-2's steel keeps no strength or stiffness at 1200 C");
            }
        }
        m_E[t] = E[t];
        // the members' checks at their own temperatures, every step
        frames[t]->set_check_temperature(m_T[t]);
        // the frame at the members' temperatures once their stiffness or strength has moved past the tolerance
        if (reduction_change(m_T[t], frames[t]->temperature()) > m_s.reduction_tolerance) {
            const std::string err = frames[t]->set_temperature(m_T[t]);
            if (!err.empty()) { Abort("conductors.thermal.towers: " + err); }
            ++m_rebuilds[t];
        }
    }
    m_started = true;
    if (fresh || write_row) { write_logs(time + dt, frames, fresh && !m_restored); }
}

std::vector<std::string> TowerThermal::log_files () const
{
    std::vector<std::string> f{m_dir + "/tower_temperature.dat"};
    for (const auto& n : m_names) { f.push_back(m_dir + "/tower_" + n + "_member_temperature.dat"); }
    return f;
}

void TowerThermal::write_logs (double t, const std::vector<erf_towers::FrameTower*>& frames, bool truncate) const
{
    if (!ParallelDescriptor::IOProcessor()) { return; }
    const std::vector<std::string> files = log_files();
    std::ofstream out;
    if (erf_actuator::open_log(out, files[0], truncate)) {
        out << "# time tower T_mean_C T_max_C member_of_max E_fire_max_W_m2 kE_min rebuilds T_frame_max_C\n";
    }
    out << std::setprecision(10);
    for (std::size_t tw = 0; tw < m_names.size(); ++tw) {
        const auto& ids = frames[tw]->frame().inputs().members;
        const auto& Tf = frames[tw]->temperature();
        double sum = 0.0, Tmax = -std::numeric_limits<double>::infinity(), Emax = 0.0, kmin = 1.0, Tfmax = -std::numeric_limits<double>::infinity();
        std::size_t at = 0;
        for (std::size_t m = 0; m < m_T[tw].size(); ++m) {
            sum += m_T[tw][m];
            if (m_T[tw][m] > Tmax) { Tmax = m_T[tw][m]; at = m; }
            Emax = std::max(Emax, m_E[tw][m]);
            double ky = 0.0, kE = 0.0;
            erf_towers::steel_reduction(Tf[m], ky, kE);
            kmin = std::min(kmin, kE);
            Tfmax = std::max(Tfmax, Tf[m]);
        }
        out << t << " " << m_names[tw] << " " << sum / static_cast<double>(std::max<std::size_t>(m_T[tw].size(), 1)) << " " << Tmax
            << " " << ids[at].id << " " << Emax << " " << kmin << " " << m_rebuilds[tw] << " " << Tfmax << "\n";
    }
    for (std::size_t tw = 0; tw < m_names.size(); ++tw) {
        std::ofstream mo;
        if (erf_actuator::open_log(mo, files[tw + 1], truncate)) {
            mo << "# time, then the temperature (C) of members";
            for (const auto& mem : frames[tw]->frame().inputs().members) { mo << " " << mem.id; }
            mo << "\n";
        }
        mo << std::setprecision(10) << t;
        for (const double x : m_T[tw]) { mo << " " << x; }
        mo << "\n";
    }
}

void TowerThermal::write_checkpoint (const std::string& dir) const
{
    if (!m_started || !ParallelDescriptor::IOProcessor()) { return; }
    std::ofstream out(dir + "/tower_thermal", std::ios::trunc);
    if (!out) { Abort("cannot write the conductor checkpoint '" + dir + "/tower_thermal'"); }
    out << std::setprecision(std::numeric_limits<double>::max_digits10) << "tower_thermal " << m_T.size() << "\n";
    for (std::size_t t = 0; t < m_T.size(); ++t) {
        out << m_names[t] << " " << m_T[t].size() << " " << m_rebuilds[t] << "\n";
        for (const auto* v : {&m_T[t], &m_Ts[t]}) {
            for (const double x : *v) { out << x << " "; }
            out << "\n";
        }
    }
    if (!out) { Abort("cannot write the conductor checkpoint '" + dir + "/tower_thermal'"); }
}

void TowerThermal::restore (const std::string& dir)
{
    m_restored = true;
    if (!erf_actuator::file_exists_everywhere(dir + "/tower_thermal")) {
        Print() << "conductors: the checkpoint '" << dir << "' holds no tower steel temperatures; they start at the restart, "
                << "at the air temperature\n";
        return;
    }
    Vector<char> chars;
    ParallelDescriptor::ReadAndBcastFile(dir + "/tower_thermal", chars);
    std::istringstream in{std::string(chars.dataPtr())};
    const std::string file = "the conductor checkpoint '" + dir + "/tower_thermal'";
    std::string tag;
    std::size_t nt = 0;
    if (!(in >> tag >> nt) || tag != "tower_thermal") { Abort("malformed " + file); }
    if (nt != m_names.size()) {
        Abort(file + " holds " + std::to_string(nt) + " heated tower(s), this run has " + std::to_string(m_names.size()) +
              "; the towers must match the run being restarted");
    }
    m_T.assign(nt, {});
    m_Ts.assign(nt, {});
    m_E.assign(nt, {});
    for (std::size_t t = 0; t < nt; ++t) {
        std::string name;
        std::size_t nm = 0;
        long nb = 0;
        if (!(in >> name >> nm >> nb) || nb < 0) { Abort("malformed " + file); }
        if (name != m_names[t] || nm != m_sec[t].size()) {
            Abort(file + ": tower " + name + " does not match tower " + m_names[t] + " of this run (its name or members)");
        }
        for (auto* v : {&m_T[t], &m_Ts[t]}) {
            v->resize(nm);
            for (auto& x : *v) { if (!(in >> x) || !std::isfinite(x)) { Abort("malformed " + file); } }
        }
        m_E[t].assign(nm, 0.0);
        m_rebuilds[t] = nb;
    }
    m_started = true;
}

} // namespace erf_conductors
