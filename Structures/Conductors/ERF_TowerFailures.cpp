#include "ERF_TowerFailures.H"

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

#include <AMReX.H>
#include <AMReX_BLassert.H>

namespace erf_conductors {

std::vector<LineRelease> tower_releases (const TowerAttachment& a)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a.num_spans >= 2 && a.tower >= 0 && a.tower < a.num_spans - 1,
                                     "tower_releases: the tower is not one of the line's");
    std::vector<LineRelease> r;
    if (a.strings) {
        if (a.string_attached) { r.push_back(LineRelease{LineRelease::Kind::String, a.tower, 0}); }
    } else {
        if (a.before_attached) { r.push_back(LineRelease{LineRelease::Kind::SpanEnd, a.tower, 1}); }
        if (a.after_attached) { r.push_back(LineRelease{LineRelease::Kind::SpanEnd, a.tower + 1, 0}); }
    }
    return r;
}

namespace {
bool one_word (const std::string& s) { return !s.empty() && s.find_first_of(" \t\n\r") == std::string::npos; }
}

void TowerFailureLog::add (const TowerFailure& f)
{
    if (!one_word(f.tower) || !one_word(f.component)) {
        amrex::Abort("TowerFailureLog::add: the tower '" + f.tower + "' and the component '" + f.component + "' must be single words");
    }
    if (!std::isfinite(f.time) || f.releases < 0) {
        amrex::Abort("TowerFailureLog::add: the failure of tower " + f.tower + " is malformed (time or number of releases)");
    }
    if (!m_f.empty() && f.time < m_f.back().time) {
        amrex::Abort("TowerFailureLog::add: the failure of tower " + f.tower + " is earlier than the last one logged");
    }
    if (failed(f.tower)) { amrex::Abort("TowerFailureLog::add: tower " + f.tower + " has failed before"); }
    m_f.push_back(f);
}

bool TowerFailureLog::failed (const std::string& tower) const
{
    for (const TowerFailure& f : m_f) { if (f.tower == tower) { return true; } }
    return false;
}

const char* TowerFailureLog::log_header () { return "time tower component releases"; }

std::string TowerFailureLog::log_row (std::size_t i) const
{
    const TowerFailure& f = m_f.at(i);
    std::ostringstream out;
    out << std::setprecision(10) << f.time << " " << f.tower << " " << f.component << " " << f.releases;
    return out.str();
}

std::string TowerFailureLog::state () const
{
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const TowerFailure& f : m_f) { out << f.tower << " " << f.component << " " << f.releases << " " << f.time << "\n"; }
    return out.str();
}

std::string TowerFailureLog::restore (const std::string& text)
{
    std::vector<TowerFailure> read;
    std::istringstream in(text);
    std::string row;
    while (std::getline(in, row)) {
        std::istringstream ls(row);
        TowerFailure f;
        std::string extra;
        if (!(ls >> f.tower)) { continue; }   // a blank line
        if (!(ls >> f.component >> f.releases >> f.time) || (ls >> extra) || f.releases < 0 || !std::isfinite(f.time) ||
            (!read.empty() && f.time < read.back().time)) {
            return "malformed tower failure record '" + row + "'";
        }
        for (const TowerFailure& o : read) {
            if (o.tower == f.tower) { return "the tower failure record '" + row + "' names a tower that failed before"; }
        }
        read.push_back(f);
    }
    m_f = read;
    return std::string();
}

} // namespace erf_conductors
