#include "ERF_LineBreaks.H"

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

#include <AMReX.H>
#include <AMReX_BLassert.H>

namespace erf_conductors {

std::string LineRelease::name () const
{
    if (kind == Kind::String) { return "string" + std::to_string(index + 1); }
    return "span" + std::to_string(index + 1) + (end == 0 ? "_end_a" : "_end_b");
}

BreakDecision decide_break (ComponentKind kind, int index, const LineEnds& e)
{
    const int N = e.num_spans();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(N >= 1 && e.tension_b.size() == e.tension_a.size() && e.attached_a.size() == e.tension_a.size() &&
                                     e.attached_b.size() == e.tension_a.size(),
                                     "decide_break: one tension and one attached flag per span end");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(e.string_attached.size() == (e.strings ? static_cast<std::size_t>(N - 1) : std::size_t(0)),
                                     "decide_break: one flag per insulator string, none on a line without strings");
    BreakDecision d;
    auto span_end = [&d] (int k, int end) {
        d.any = true;
        d.release = LineRelease{LineRelease::Kind::SpanEnd, k, end};
    };
    // of two candidate span ends, the attached one with the higher tension; the first on a tie, which is a second
    // tension not higher than the first by more than break_tie_tolerance of it
    auto higher = [&] (int k0, int end0, int k1, int end1) {
        const auto s0 = static_cast<std::size_t>(k0), s1 = static_cast<std::size_t>(k1);
        const bool on0 = (end0 == 0 ? e.attached_a[s0] : e.attached_b[s0]) != 0;
        const bool on1 = (end1 == 0 ? e.attached_a[s1] : e.attached_b[s1]) != 0;
        const double t0 = (end0 == 0) ? e.tension_a[s0] : e.tension_b[s0];
        const double t1 = (end1 == 0) ? e.tension_a[s1] : e.tension_b[s1];
        if (on0 && (!on1 || !(t1 > t0 * (1.0 + break_tie_tolerance)))) { span_end(k0, end0); }
        else if (on1) { span_end(k1, end1); }
    };
    switch (kind) {
    case ComponentKind::Conductor:
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(index >= 0 && index < N, "decide_break: no such span");
        higher(index, 0, index, 1);
        break;
    case ComponentKind::Insulator:
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(e.strings && index >= 0 && index < N - 1, "decide_break: no such insulator string");
        if (e.string_attached[static_cast<std::size_t>(index)] != 0) {
            d.any = true;
            d.release = LineRelease{LineRelease::Kind::String, index, 0};
        }
        break;
    case ComponentKind::Hardware:
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(index >= 0 && index <= N, "decide_break: no such attachment");
        if (index == 0) {
            if (e.attached_a[0] != 0) { span_end(0, 0); }
        } else if (index == N) {
            if (e.attached_b[static_cast<std::size_t>(N - 1)] != 0) { span_end(N - 1, 1); }
        } else if (e.strings) {
            if (e.string_attached[static_cast<std::size_t>(index - 1)] != 0) {
                d.any = true;
                d.release = LineRelease{LineRelease::Kind::String, index - 1, 0};
            }
        } else {
            higher(index - 1, 1, index, 0);
        }
        break;
    case ComponentKind::FootingUplift:
    case ComponentKind::FootingCompression:
    case ComponentKind::Member:
        break;
    }
    return d;
}

namespace {
bool one_word (const std::string& s) { return !s.empty() && s.find_first_of(" \t\n\r") == std::string::npos; }
}

void BreakLog::add (const BreakRecord& r)
{
    if (!one_word(r.line) || !one_word(r.component)) {
        amrex::Abort("BreakLog::add: the line '" + r.line + "' and the component '" + r.component + "' must be single words");
    }
    if (!(std::isfinite(r.time) && std::isfinite(r.tension) && r.tension >= 0.0) || r.release.index < 0 ||
        (r.release.end != 0 && r.release.end != 1) || (r.release.kind == LineRelease::Kind::String && r.release.end != 0)) {
        amrex::Abort("BreakLog::add: the release of line " + r.line + " by " + r.component + " is malformed (time, tension, index or end)");
    }
    if (!m_r.empty() && r.time < m_r.back().time) {
        amrex::Abort("BreakLog::add: the release of line " + r.line + " by " + r.component + " is earlier than the last one logged");
    }
    for (const BreakRecord& o : m_r) {
        if (o.line == r.line && o.release == r.release) {
            amrex::Abort("BreakLog::add: line " + r.line + " was already released at " + r.release.name());
        }
    }
    m_r.push_back(r);
}

std::vector<LineRelease> BreakLog::releases_of (const std::string& line) const
{
    std::vector<LineRelease> out;
    for (const BreakRecord& r : m_r) { if (r.line == line) { out.push_back(r.release); } }
    return out;
}

const char* BreakLog::log_header () { return "time line released component tension"; }

std::string BreakLog::log_row (std::size_t i) const
{
    const BreakRecord& r = m_r.at(i);
    std::ostringstream out;
    out << std::setprecision(10) << r.time << " " << r.line << " " << r.release.name() << " " << r.component << " " << r.tension;
    return out.str();
}

std::string BreakLog::state () const
{
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const BreakRecord& r : m_r) {
        out << r.line << " " << (r.release.kind == LineRelease::Kind::String ? "string" : "span") << " " << r.release.index << " "
            << r.release.end << " " << r.component << " " << r.time << " " << r.tension << "\n";
    }
    return out.str();
}

std::string BreakLog::restore (const std::string& text)
{
    std::vector<BreakRecord> read;
    std::istringstream in(text);
    std::string row;
    while (std::getline(in, row)) {
        std::istringstream ls(row);
        BreakRecord r;
        std::string kind, extra;
        if (!(ls >> r.line)) { continue; }   // a blank line
        if (!(ls >> kind >> r.release.index >> r.release.end >> r.component >> r.time >> r.tension) || (ls >> extra) ||
            (kind != "span" && kind != "string") || r.release.index < 0 || (r.release.end != 0 && r.release.end != 1) ||
            (kind == "string" && r.release.end != 0) || !std::isfinite(r.time) || !(std::isfinite(r.tension) && r.tension >= 0.0) ||
            (!read.empty() && r.time < read.back().time)) {
            return "malformed release record '" + row + "'";
        }
        r.release.kind = (kind == "string") ? LineRelease::Kind::String : LineRelease::Kind::SpanEnd;
        for (const BreakRecord& o : read) {
            if (o.line == r.line && o.release == r.release) { return "the release record '" + row + "' is given twice"; }
        }
        read.push_back(r);
    }
    m_r = read;
    return std::string();
}

} // namespace erf_conductors
