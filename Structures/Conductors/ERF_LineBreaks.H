#ifndef ERF_LINE_BREAKS_H
#define ERF_LINE_BREAKS_H

/**
 * LineBreaks: what a conductor line lets go of when one of its components exceeds its capacity
 * (conductors.failure_break), and the log of the releases made.
 *
 * A release detaches one end of a span from its support, or an insulator string from its cross-arm;
 * the line goes on from the next step with that end free. decide_break() chooses the release for a
 * component's first exceedance from the line's layout and its tensions; BreakLog keeps the releases in
 * the order they were made, writes breaks.dat and carries them through a checkpoint, where a restart
 * needs them in order: MoorDyn's saved state fits only a system with the same releases made again.
 * Nothing here knows of MoorDyn: ConductorLine makes the release (ERF_ConductorLine.H).
 *
 * Conventions: spans k and strings j are numbered from 0 (string j hangs at tower j+1); a span's end 0
 * is the one towards the line's end_a, its end 1 the one towards end_b; tensions in N, times in s.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "ERF_Capacities.H"

namespace erf_conductors {

/** A support letting go of the line: a span's end, or an insulator string's top. */
struct LineRelease
{
    enum class Kind { SpanEnd, String };
    Kind kind{Kind::SpanEnd};
    int index{0};   //!< the span k, or the string j (both from 0)
    int end{0};     //!< for a span: 0 its end towards end_a, 1 its end towards end_b; 0 for a string
    bool operator== (const LineRelease& o) const { return kind == o.kind && index == o.index && end == o.end; }
    /** span<k+1>_end_a, span<k+1>_end_b or string<j+1>: how the logs name it. */
    std::string name () const;
};

/** A line's layout and state, as decide_break() needs them. */
struct LineEnds
{
    bool strings{false};                   //!< the line hangs from insulator strings at its towers
    std::vector<double> tension_a;         //!< per span: the tension at its end 0 (N)
    std::vector<double> tension_b;         //!< per span: the tension at its end 1 (N)
    std::vector<char> attached_a;          //!< per span: its end 0 is still on its support
    std::vector<char> attached_b;          //!< per span: its end 1 is still on its support
    std::vector<char> string_attached;     //!< per string: it still hangs from its cross-arm
    int num_spans () const { return static_cast<int>(tension_a.size()); }
};

/** What a component's exceedance releases; any is false when there is nothing left for it to release. */
struct BreakDecision
{
    bool any{false};
    LineRelease release;
};

/**
 * Two tensions closer than this fraction of the first are a tie for decide_break(): on a level span
 * in still air the two ends' tensions differ by the line solver's residual only (1e-8 of the tension
 * with MoorDyn-C), which must not choose the end.
 */
constexpr double break_tie_tolerance = 1.0e-6;

/**
 * The release for the first exceedance of a line's component:
 *  - a conductor (index = its span k): the span lets go at its end with the higher tension (end 0 on
 *    a tie, break_tie_tolerance), or at the one end still attached;
 *  - an insulator string (index = the string j): the string lets go of its cross-arm;
 *  - the hardware at attachment a (index = a, 0 at end_a to num_spans at end_b): at a dead end the
 *    span there lets go of it; at a tower with strings the string lets go of the cross-arm; at a
 *    tower the spans are clamped to, the span with the higher tension at the tower lets go of it (the
 *    span before the tower on a tie), or the one still attached.
 * Footings and members release nothing here (any = false; with conductors.tower_failure = release their
 * tower lets go of its lines, ERF_TowerFailures.H), as does a component whose ends were all released
 * before. Aborts on an index outside the line or on vectors of unequal length.
 */
BreakDecision decide_break (ComponentKind kind, int index, const LineEnds& e);

/** One release made. */
struct BreakRecord
{
    double time{0.0};        //!< the end of the step at which the line let go (s)
    std::string line;        //!< the line's name in conductors.lines
    LineRelease release;
    std::string component;   //!< the component whose exceedance released it (its name in capacities.csv)
    double tension{0.0};     //!< the tension at the released end, or at the string's top, when it let go (N)
};

class BreakLog
{
public:
    /** Add a release (the latest); its line and component names must hold no blanks and not be empty. */
    void add (const BreakRecord& r);
    std::size_t size () const { return m_r.size(); }
    bool empty () const { return m_r.empty(); }
    const BreakRecord& record (std::size_t i) const { return m_r.at(i); }
    /** The releases of one line, in the order they were made. */
    std::vector<LineRelease> releases_of (const std::string& line) const;

    /** The header of breaks.dat: one row per release. */
    static const char* log_header ();
    /** The row of release i for breaks.dat (whitespace separated, no newline). */
    std::string log_row (std::size_t i) const;

    /** The log as text: one record per release in order, every value at 17 significant digits. */
    std::string state () const;
    /** Replace the log by what state() wrote. Returns an empty string, or the message for a malformed record. */
    std::string restore (const std::string& text);

private:
    std::vector<BreakRecord> m_r;
};

} // namespace erf_conductors

#endif
