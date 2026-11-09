#include "ERF_GroundContact.H"

#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>

#include <AMReX.H>

namespace erf_conductors {

namespace {
constexpr int k_digits = 17;   // a double's every bit, so a restart continues the records exactly
}

std::size_t
GroundContactLog::add (const std::string& name, int nodes)
{
    if (name.empty()) { amrex::Abort("GroundContactLog: a span needs a name"); }
    if (nodes < 2) { amrex::Abort("GroundContactLog: the span '" + name + "' needs at least two nodes"); }
    for (const auto& s : m_s) {
        if (s.name == name) { amrex::Abort("GroundContactLog: two spans are named '" + name + "'"); }
    }
    GroundSpan s;
    s.name = name;
    s.nodes = nodes;
    m_s.push_back(s);
    return m_s.size() - 1;
}

bool
GroundContactLog::record (std::size_t i, int on_ground, double t)
{
    GroundSpan& s = m_s.at(i);
    if (on_ground < 0 || on_ground > s.nodes) {
        amrex::Abort("GroundContactLog: " + std::to_string(on_ground) + " nodes of the span '" + s.name + "' on the ground, of " +
                     std::to_string(s.nodes));
    }
    s.now = on_ground;
    if (on_ground > s.most) {
        s.most = on_ground;
        s.most_time = t;
    }
    if (on_ground > 0 && !s.touched) {
        s.touched = true;
        s.first_time = t;
        return true;
    }
    return false;
}

void
GroundContactLog::write_table (const std::string& fname) const
{
    std::ofstream out(fname, std::ios::trunc);
    if (!out) { amrex::Abort("GroundContactLog: cannot write '" + fname + "'"); }
    out << "span,nodes,nodes_on_ground,most_on_ground,most_time,first_touch_time\n" << std::setprecision(10);
    for (const auto& s : m_s) {
        out << s.name << "," << s.nodes << "," << s.now << "," << s.most << ",";
        if (s.most > 0) { out << s.most_time; }
        out << ",";
        if (s.touched) { out << s.first_time; }
        out << "\n";
    }
}

std::string
GroundContactLog::state () const
{
    std::ostringstream o;
    o << std::setprecision(k_digits);
    for (const auto& s : m_s) {
        o << "span " << s.name << " " << s.nodes << " " << s.now << " " << s.most << " " << s.most_time << " "
          << (s.touched ? 1 : 0) << " " << s.first_time << "\n";
    }
    return o.str();
}

std::string
GroundContactLog::restore (const std::string& text, std::vector<std::string>& names_dropped)
{
    names_dropped.clear();
    std::map<std::string, std::size_t> index;
    for (std::size_t i = 0; i < m_s.size(); ++i) { index[m_s[i].name] = i; }
    std::istringstream in(text);
    std::string line;
    int n = 0;
    while (std::getline(in, line)) {
        ++n;
        std::istringstream ls(line);
        std::string key;
        if (!(ls >> key)) { continue; }
        GroundSpan r;
        int touched = -1;
        if (key != "span" || !(ls >> r.name >> r.nodes >> r.now >> r.most >> r.most_time >> touched >> r.first_time) ||
            (touched != 0 && touched != 1) || r.nodes < 2 || r.now < 0 || r.now > r.nodes || r.most < r.now || r.most > r.nodes) {
            return "malformed ground-contact record on line " + std::to_string(n);
        }
        r.touched = (touched == 1);
        const auto it = index.find(r.name);
        if (it == index.end() || m_s[it->second].nodes != r.nodes) {
            names_dropped.push_back(r.name);
            continue;
        }
        m_s[it->second] = r;
    }
    return "";
}

} // namespace erf_conductors
