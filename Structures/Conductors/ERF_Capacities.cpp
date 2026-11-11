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
constexpr int k_kinds = 6;
} // namespace

const char*
kind_name (ComponentKind kind)
{
    switch (kind) {
        case ComponentKind::Conductor:          return "conductor";
        case ComponentKind::Insulator:          return "insulator";
        case ComponentKind::Hardware:           return "hardware";
        case ComponentKind::FootingUplift:      return "footing_uplift";
        case ComponentKind::FootingCompression: return "footing_compression";
        case ComponentKind::Member:             return "member";
    }
    return "unknown";
}

std::size_t
CapacityMonitor::add (ComponentKind kind, const std::string& name)
{
    if (name.empty() || name.find_first_of(" \t\n,") != std::string::npos) {
        amrex::Abort("CapacityMonitor::add: a component name must be one word without commas, got '" + name + "'");
    }
    for (const auto& c : m_c) {
        if (c.name == name) {
            amrex::Abort("conductors: two components checked against their capacities are named '" + name + "': rename a "
                         "line so that its name is not another line's component name (<line>_span<k>, <line>_string<j>, "
                         "<line>_end_a, <line>_end_b, <line>_tower<j>)");
        }
    }
    CapacityComponent c;
    c.kind = kind;
    c.name = name;
    m_c.push_back(c);
    return m_c.size() - 1;
}

bool
CapacityMonitor::record (std::size_t i, double demand, double capacity, double t)
{
    CapacityComponent& c = m_c.at(i);
    if (!(std::isfinite(demand) && demand >= 0.0)) {
        amrex::Abort(std::string("conductors: the demand on ") + kind_name(c.kind) + " " + c.name +
                     " is negative or not finite");
    }
    if (!(std::isfinite(capacity) && capacity > 0.0)) {
        amrex::Abort(std::string("conductors: the capacity of ") + kind_name(c.kind) + " " + c.name +
                     " is not positive and finite");
    }
    if (!std::isfinite(t)) { amrex::Abort("CapacityMonitor::record: the time is not finite"); }
    c.demand = demand;
    c.capacity = capacity;
    const double ratio = demand / capacity;
    if (!c.seen || ratio > c.peak_ratio) {
        c.peak_ratio = ratio;
        c.peak_demand = demand;
        c.peak_capacity = capacity;
        c.peak_time = t;
    }
    c.seen = true;
    if (!c.exceeded && demand > capacity) {
        c.exceeded = true;
        c.first_time = t;
        c.first_demand = demand;
        c.first_capacity = capacity;
        return true;
    }
    return false;
}

std::size_t
CapacityMonitor::first_failure () const
{
    std::size_t first = m_c.size();
    for (std::size_t i = 0; i < m_c.size(); ++i) {
        if (!m_c[i].exceeded) { continue; }
        if (first == m_c.size() || m_c[i].first_time < m_c[first].first_time) { first = i; }
    }
    return first;
}

void
CapacityMonitor::order_by_first_time (std::vector<std::size_t>& ids) const
{
    // (at one time: in the order they were added, which is capacities.csv's)
    std::sort(ids.begin(), ids.end(), [&] (std::size_t a, std::size_t b) {
        const double ta = m_c.at(a).first_time, tb = m_c.at(b).first_time;
        return ta < tb || (ta == tb && a < b);
    });
}

void
CapacityMonitor::write_table (const std::string& fname) const
{
    std::ofstream out(fname, std::ios::trunc);
    if (!out) { amrex::Abort("cannot write '" + fname + "'"); }
    out << "kind,component,demand,capacity,peak_ratio,peak_demand,peak_capacity,peak_time,"
           "exceeded,first_time,first_demand,first_capacity\n" << std::setprecision(10);
    for (const auto& c : m_c) {
        out << kind_name(c.kind) << "," << c.name << ",";
        if (c.seen) {
            out << c.demand << "," << c.capacity << "," << c.peak_ratio << "," << c.peak_demand << "," << c.peak_capacity << ","
                << c.peak_time << ",";
        } else {
            out << ",,,,,,";
        }
        out << (c.exceeded ? 1 : 0) << ",";
        if (c.exceeded) {
            out << c.first_time << "," << c.first_demand << "," << c.first_capacity;
        } else {
            out << ",,";
        }
        out << "\n";
    }
    if (!out) { amrex::Abort("cannot write '" + fname + "'"); }
}

const char*
CapacityMonitor::failure_log_header ()
{
    return "time kind component demand capacity ratio";
}

std::string
CapacityMonitor::failure_row (std::size_t i) const
{
    const CapacityComponent& c = m_c.at(i);
    std::ostringstream s;
    s << std::setprecision(10) << c.first_time << " " << kind_name(c.kind) << " " << c.name << " " << c.first_demand << " "
      << c.first_capacity << " " << c.first_demand / c.first_capacity;
    return s.str();
}

std::string
CapacityMonitor::state () const
{
    std::ostringstream s;
    s << std::setprecision(k_digits);
    for (const auto& c : m_c) {
        s << "component " << c.name << " " << static_cast<int>(c.kind) << " " << (c.seen ? 1 : 0) << " " << c.demand << " "
          << c.capacity << " " << c.peak_ratio << " " << c.peak_demand << " " << c.peak_capacity << " " << c.peak_time << " "
          << (c.exceeded ? 1 : 0) << " " << c.first_time << " " << c.first_demand << " " << c.first_capacity << "\n";
    }
    return s.str();
}

std::string
CapacityMonitor::restore (const std::string& text, std::vector<std::string>& names_dropped)
{
    names_dropped.clear();
    std::map<std::string, std::size_t> index;
    for (std::size_t i = 0; i < m_c.size(); ++i) { index[m_c[i].name] = i; }
    std::istringstream in(text);
    std::string line;
    int n = 0;
    while (std::getline(in, line)) {
        ++n;
        std::istringstream ls(line);
        std::string key;
        if (!(ls >> key)) { continue; }
        CapacityComponent r;
        int kind = -1, seen = 0, exceeded = 0;
        if (key != "component" ||
            !(ls >> r.name >> kind >> seen >> r.demand >> r.capacity >> r.peak_ratio >> r.peak_demand >> r.peak_capacity
                 >> r.peak_time >> exceeded >> r.first_time >> r.first_demand >> r.first_capacity) ||
            kind < 0 || kind >= k_kinds || (seen != 0 && seen != 1) || (exceeded != 0 && exceeded != 1)) {
            return "malformed capacity record on line " + std::to_string(n);
        }
        r.kind = static_cast<ComponentKind>(kind);
        r.seen = (seen == 1);
        r.exceeded = (exceeded == 1);
        const auto it = index.find(r.name);
        if (it == index.end() || m_c[it->second].kind != r.kind) {
            names_dropped.push_back(r.name);
            continue;
        }
        m_c[it->second] = r;
    }
    return "";
}

void
member_demand (double tension, double compression, double tension_capacity, double compression_capacity,
               double& demand, double& capacity)
{
    if (!(std::isfinite(tension) && std::isfinite(compression))) {
        amrex::Abort("member_demand: a member's axial force is not finite");
    }
    if (!(tension_capacity > 0.0 && compression_capacity > 0.0)) {
        amrex::Abort("member_demand: a member's design strengths must be positive (its steel is at 1200 C or above)");
    }
    if (tension / tension_capacity >= compression / compression_capacity) {
        demand = tension;
        capacity = tension_capacity;
    } else {
        demand = compression;
        capacity = compression_capacity;
    }
}

} // namespace erf_conductors
