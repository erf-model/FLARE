#include "ERF_Intensity.H"
#include "ERF_Capacities.H"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>

#include <AMReX.H>

namespace erf_conductors {

namespace {
constexpr int k_digits = std::numeric_limits<double>::max_digits10;
constexpr double k_pi = 3.14159265358979323846;

void put (std::ostream& s, const IntensityValue& v)
{
    s << (v.set ? 1 : 0) << " " << v.speed << " " << v.direction << " " << v.k << " " << v.time;
}

bool get (std::istream& s, IntensityValue& v)
{
    int set = -1;
    if (!(s >> set >> v.speed >> v.direction >> v.k >> v.time) || (set != 0 && set != 1)) { return false; }
    v.set = (set == 1);
    return true;
}
} // namespace

std::string
IntensityRecorder::check (double height, double averaging_time)
{
    if (!(std::isfinite(height) && height > 0.0 && height <= 500.0)) {
        return "conductors.intensity_height must be in (0, 500] (m above the ground)";
    }
    if (!(std::isfinite(averaging_time) && averaging_time >= 0.0 && averaging_time <= 600.0)) {
        return "conductors.intensity_averaging_time must be in [0, 600] (s; 0 keeps the last sample)";
    }
    return "";
}

std::string
IntensityRecorder::configure (double height, double averaging_time, bool with_k)
{
    const std::string err = check(height, averaging_time);
    if (!err.empty()) { return err; }
    *this = IntensityRecorder();
    m_height = height;
    m_T = averaging_time;
    m_with_k = with_k;
    return "";
}

std::size_t
IntensityRecorder::add_point (const std::string& name, double x, double y, double z)
{
    if (name.empty() || name.find_first_of(" \t\n,") != std::string::npos) {
        amrex::Abort("IntensityRecorder::add_point: a point name must be one word without commas, got '" + name + "'");
    }
    for (const auto& p : m_p) {
        if (p.name == name) { amrex::Abort("IntensityRecorder::add_point: two intensity points are named '" + name + "'"); }
    }
    if (!(std::isfinite(x) && std::isfinite(y) && std::isfinite(z))) {
        amrex::Abort("IntensityRecorder::add_point: the position of '" + name + "' is not finite");
    }
    if (m_nsamples > 0) { amrex::Abort("IntensityRecorder::add_point: points are added before the first sample"); }
    IntensityPoint p;
    p.name = name;
    p.x = x;
    p.y = y;
    p.z = z;
    m_p.push_back(p);
    m_hist.emplace_back();
    return m_p.size() - 1;
}

void
IntensityRecorder::components (std::size_t n)
{
    m_link.assign(n, m_p.size());
    m_peak.assign(n, IntensityValue());
    m_first.assign(n, IntensityValue());
}

void
IntensityRecorder::link (std::size_t c, std::size_t p)
{
    if (c >= m_link.size()) { amrex::Abort("IntensityRecorder::link: no such component (components() first)"); }
    if (p >= m_p.size()) { amrex::Abort("IntensityRecorder::link: no such point"); }
    if (m_link[c] < m_p.size()) { amrex::Abort("IntensityRecorder::link: the component already has a point"); }
    m_link[c] = p;
}

std::size_t
IntensityRecorder::point_of (std::size_t c) const
{
    const std::size_t p = m_link.at(c);
    return (p < m_p.size()) ? p : m_p.size();
}

void
IntensityRecorder::average (std::size_t p, IntensityValue& out) const
{
    // newest first: sample i covers the dt_i before the ones newer than it; the window is the last T seconds
    double left = m_T, wsum = 0.0, s = 0.0, u = 0.0, v = 0.0, k = 0.0;
    const auto& h = m_hist[p];
    for (std::size_t i = 0; i < h.size(); ++i) {
        // T = 0: the newest sample alone, with its whole step as the weight
        const double w = (m_T > 0.0) ? std::min(m_dt[i], left) : m_dt[i];
        s += w * std::sqrt(h[i].u * h[i].u + h[i].v * h[i].v);
        u += w * h[i].u;
        v += w * h[i].v;
        k += w * h[i].k;
        wsum += w;
        left -= w;
        if (!(left > 0.0)) { break; }
    }
    out.speed = s / wsum;
    out.direction = (u == 0.0 && v == 0.0) ? 0.0 : std::atan2(v, u) * 180.0 / k_pi;
    out.k = m_with_k ? k / wsum : 0.0;
    out.time = m_time;
    out.set = true;
}

void
IntensityRecorder::sample (const std::vector<IntensitySample>& s, double dt, double t)
{
    if (s.size() != m_p.size()) { amrex::Abort("IntensityRecorder::sample: one sample per point"); }
    if (!(std::isfinite(dt) && dt > 0.0)) { amrex::Abort("IntensityRecorder::sample: the step must be positive and finite"); }
    if (!std::isfinite(t) || (m_has_time && !(t > m_time))) {
        amrex::Abort("IntensityRecorder::sample: the time must be finite and later than the last sample's");
    }
    for (std::size_t p = 0; p < s.size(); ++p) {
        if (!(std::isfinite(s[p].u) && std::isfinite(s[p].v))) {
            amrex::Abort("conductors: the wind at the intensity point " + m_p[p].name + " is not finite");
        }
        if (m_with_k && !(std::isfinite(s[p].k) && s[p].k >= 0.0)) {
            amrex::Abort("conductors: the turbulent kinetic energy at the intensity point " + m_p[p].name +
                         " is negative or not finite");
        }
    }
    m_time = t;
    m_has_time = true;
    ++m_nsamples;
    m_dt.push_front(dt);
    // keep the steps that reach into the window: drop the oldest while the newer ones already fill it
    double newer = 0.0;
    std::size_t keep = 0;
    for (std::size_t i = 0; i < m_dt.size(); ++i) {
        keep = i + 1;
        newer += m_dt[i];
        if (!(newer < m_T)) { break; }
    }
    while (m_dt.size() > keep) { m_dt.pop_back(); }
    // the window is full once the kept steps cover it (to round-off: ten steps of 0.3 s sum to just under 3 s)
    const bool full = !(newer < m_T * (1.0 - 1.0e-12));
    for (std::size_t p = 0; p < s.size(); ++p) {
        m_hist[p].push_front(s[p]);
        while (m_hist[p].size() > keep) { m_hist[p].pop_back(); }
        average(p, m_p[p].now);
        // the run's largest T-second average: an average over a window not yet full is a shorter one, and is not compared
        // (nor one stamped before the start time of the checks)
        if (full && !(t < m_max_start) && (!m_p[p].max.set || m_p[p].now.speed > m_p[p].max.speed)) { m_p[p].max = m_p[p].now; }
    }
}

void
IntensityRecorder::at_peak (std::size_t c)
{
    const std::size_t p = point_of(c);
    if (p < m_p.size() && m_p[p].now.set) { m_peak.at(c) = m_p[p].now; }
}

void
IntensityRecorder::at_first (std::size_t c)
{
    const std::size_t p = point_of(c);
    if (p < m_p.size() && m_p[p].now.set) { m_first.at(c) = m_p[p].now; }
}

void
IntensityRecorder::write_points (const std::string& fname) const
{
    std::ofstream out(fname, std::ios::trunc);
    if (!out) { amrex::Abort("cannot write '" + fname + "'"); }
    out << "point,x,y,z,height,speed,direction,k,time,max_speed,max_direction,max_k,max_time\n" << std::setprecision(10);
    for (const auto& p : m_p) {
        out << p.name << "," << p.x << "," << p.y << "," << p.z << "," << m_height << ",";
        auto value = [&] (const IntensityValue& v, bool last) {
            if (v.set) {
                out << v.speed << "," << v.direction << ",";
                if (m_with_k) { out << v.k; }
                out << "," << v.time;
            } else {
                out << ",,,";
            }
            out << (last ? "\n" : ",");
        };
        value(p.now, false);
        value(p.max, true);
    }
    if (!out) { amrex::Abort("cannot write '" + fname + "'"); }
}

void
IntensityRecorder::write_demand_table (const std::string& fname, const CapacityMonitor& cap) const
{
    if (cap.size() != m_link.size()) { amrex::Abort("IntensityRecorder::write_demand_table: not the monitor the links were made for"); }
    std::ofstream out(fname, std::ios::trunc);
    if (!out) { amrex::Abort("cannot write '" + fname + "'"); }
    out << "kind,component,peak_ratio,peak_demand,peak_capacity,peak_time,point,"
           "peak_speed,peak_direction,peak_k,peak_wind_time,max_speed,max_direction,max_k,max_time,"
           "exceeded,first_time,first_demand,first_capacity,first_speed,first_direction,first_k,first_wind_time\n"
        << std::setprecision(10);
    auto wind = [&] (const IntensityValue& v, bool last) {   // speed,direction,k,time
        if (v.set) {
            out << v.speed << "," << v.direction << ",";
            if (m_with_k) { out << v.k; }
            out << "," << v.time;
        } else {
            out << ",,,";
        }
        if (!last) { out << ","; }
    };
    for (std::size_t i = 0; i < cap.size(); ++i) {
        const CapacityComponent& c = cap.component(i);
        const std::size_t p = point_of(i);
        const bool has = (p < m_p.size());
        out << kind_name(c.kind) << "," << c.name << ",";
        if (c.seen) {
            out << c.peak_ratio << "," << c.peak_demand << "," << c.peak_capacity << "," << c.peak_time << ",";
        } else {
            out << ",,,,";
        }
        out << (has ? m_p[p].name : std::string()) << ",";
        wind(c.seen ? m_peak[i] : IntensityValue(), false);
        wind((has && m_p[p].max.set) ? m_p[p].max : IntensityValue(), false);
        out << (c.exceeded ? 1 : 0) << ",";
        if (c.exceeded) {
            out << c.first_time << "," << c.first_demand << "," << c.first_capacity << ",";
            wind(m_first[i], true);
        } else {
            out << ",,,,,,";
        }
        out << "\n";
    }
    if (!out) { amrex::Abort("cannot write '" + fname + "'"); }
}

std::string
IntensityRecorder::state (const CapacityMonitor& cap) const
{
    if (cap.size() != m_link.size()) { amrex::Abort("IntensityRecorder::state: not the monitor the links were made for"); }
    std::ostringstream s;
    s << std::setprecision(k_digits);
    s << "intensity " << m_height << " " << m_T << " " << (m_with_k ? 1 : 0) << " " << m_nsamples << " "
      << (m_has_time ? 1 : 0) << " " << m_time << "\n";
    s << "steps " << m_dt.size();
    for (const double dt : m_dt) { s << " " << dt; }
    s << "\n";
    for (std::size_t p = 0; p < m_p.size(); ++p) {
        s << "point " << m_p[p].name << " ";
        put(s, m_p[p].now);
        s << " ";
        put(s, m_p[p].max);
        s << " " << m_hist[p].size();
        for (const auto& h : m_hist[p]) { s << " " << h.u << " " << h.v << " " << h.k; }
        s << "\n";
    }
    for (std::size_t c = 0; c < m_link.size(); ++c) {
        if (point_of(c) >= m_p.size()) { continue; }
        s << "component " << cap.component(c).name << " ";
        put(s, m_peak[c]);
        s << " ";
        put(s, m_first[c]);
        s << "\n";
    }
    return s.str();
}

std::string
IntensityRecorder::restore (const std::string& text, const CapacityMonitor& cap, std::vector<std::string>& names_dropped)
{
    names_dropped.clear();
    if (cap.size() != m_link.size()) { amrex::Abort("IntensityRecorder::restore: not the monitor the links were made for"); }
    std::map<std::string, std::size_t> points, comps;
    for (std::size_t p = 0; p < m_p.size(); ++p) { points[m_p[p].name] = p; }
    for (std::size_t c = 0; c < cap.size(); ++c) { if (point_of(c) < m_p.size()) { comps[cap.component(c).name] = c; } }
    std::istringstream in(text);
    std::string line;
    int n = 0;
    bool head = false, steps = false, points_read = false;
    std::deque<double> dts;
    while (std::getline(in, line)) {
        ++n;
        std::istringstream ls(line);
        std::string key;
        if (!(ls >> key)) { continue; }
        const std::string bad = "malformed intensity record on line " + std::to_string(n);
        if (key == "intensity") {
            double height = 0.0, T = 0.0, time = 0.0;
            int with_k = -1, has_time = -1;
            long long nsamples = -1;
            if (head || !(ls >> height >> T >> with_k >> nsamples >> has_time >> time) || (with_k != 0 && with_k != 1) ||
                nsamples < 0 || (has_time != 0 && has_time != 1)) { return bad; }
            std::ostringstream m;
            m << std::setprecision(10);
            if (height != m_height) {
                m << "it was written with conductors.intensity_height = " << height << ", the run has " << m_height;
                return m.str();
            }
            if (T != m_T) {
                m << "it was written with conductors.intensity_averaging_time = " << T << ", the run has " << m_T;
                return m.str();
            }
            if ((with_k == 1) != m_with_k) {
                return std::string("it was written ") + (with_k == 1 ? "with" : "without") + " a turbulent kinetic energy and the run is " +
                       (m_with_k ? "with" : "without") + " one (erf.rans_type on the anchor level or conductors.prescribed_velocity changed)";
            }
            m_nsamples = static_cast<std::size_t>(nsamples);
            m_has_time = (has_time == 1);
            m_time = time;
            head = true;
        } else if (key == "steps") {
            // once, before any point: the points' samples are those of these steps, newest first
            if (!head || steps || points_read) { return bad; }
            long long count = -1;
            if (!(ls >> count) || count < 0) { return bad; }
            steps = true;
            dts.clear();
            for (long long i = 0; i < count; ++i) {
                double dt = 0.0;
                if (!(ls >> dt) || !(dt > 0.0)) { return bad; }
                dts.push_back(dt);
            }
            m_dt = dts;
        } else if (key == "point") {
            if (!head || !steps) { return bad; }
            points_read = true;
            std::string name;
            IntensityValue now, max;
            long long count = -1;
            // a point that joined at an earlier restart holds fewer samples than the window's steps: its newest ones
            if (!(ls >> name) || !get(ls, now) || !get(ls, max) || !(ls >> count) || count < 0 ||
                static_cast<std::size_t>(count) > dts.size()) { return bad; }
            std::deque<IntensitySample> h;
            for (long long i = 0; i < count; ++i) {
                IntensitySample smp;
                if (!(ls >> smp.u >> smp.v >> smp.k)) { return bad; }
                h.push_back(smp);
            }
            const auto it = points.find(name);
            if (it == points.end()) { names_dropped.push_back(name); continue; }
            m_p[it->second].now = now;
            m_p[it->second].max = max;
            m_hist[it->second] = h;
            points.erase(it);
        } else if (key == "component") {
            if (!head) { return bad; }
            std::string name;
            IntensityValue peak, first;
            if (!(ls >> name) || !get(ls, peak) || !get(ls, first)) { return bad; }
            const auto it = comps.find(name);
            if (it == comps.end()) { names_dropped.push_back(name); continue; }
            m_peak[it->second] = peak;
            m_first[it->second] = first;
        } else {
            return bad;
        }
    }
    if (!head) { return "it holds no intensity record"; }
    // a point the checkpoint did not hold keeps an empty history: its window fills from the next sample on
    return "";
}

} // namespace erf_conductors
