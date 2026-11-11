#ifndef ERF_TOWER_FAILURES_H
#define ERF_TOWER_FAILURES_H

/**
 * TowerFailures: what a tower lets go of when its footing or one of its members exceeds its capacity
 * (conductors.tower_failure = release), and the log of the towers that failed.
 *
 * A tower that fails releases every line attached to it: an insulator string leaves the cross-arm, and
 * a line clamped to the tower leaves it on both sides. The tower itself stays where it is and carries no
 * line afterwards. tower_releases() lists the releases for one line's attachment to the tower;
 * TowerFailureLog keeps the failed towers in order, writes tower_failures.dat and carries them through
 * a checkpoint. Nothing here knows of MoorDyn or of the towers' models: Conductors hands it the
 * attachments and ConductorLine makes the releases (ERF_LineBreaks.H, ERF_ConductorLine.H).
 *
 * Conventions: a line's towers j are numbered from 0 (tower j stands between spans j and j+1, both
 * from 0); times in s.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "ERF_LineBreaks.H"

namespace erf_conductors {

/** One line's attachment to a tower, as tower_releases() needs it. */
struct TowerAttachment
{
    bool strings{false};          //!< the line hangs from an insulator string at the tower
    int num_spans{2};             //!< the line's spans (its towers are num_spans - 1)
    int tower{0};                 //!< the tower's index j on the line (from 0)
    bool string_attached{true};   //!< with strings: the string still hangs from the cross-arm
    bool before_attached{true};   //!< without strings: span j's end 1 is still on the tower
    bool after_attached{true};    //!< without strings: span j+1's end 0 is still on the tower
};

/**
 * The releases that take the line off the tower, in the order to make them: the string, or span j's
 * end 1 and then span j+1's end 0; what was released before is left out, so the list may be empty.
 * Aborts when the tower is not one of the line's.
 */
std::vector<LineRelease> tower_releases (const TowerAttachment& a);

/** One tower's failure. */
struct TowerFailure
{
    double time{0.0};        //!< the end of the step at which its lines were released (s)
    std::string tower;       //!< the tower's name
    std::string component;   //!< the footing or member whose exceedance failed it (its name in capacities.csv)
    int releases{0};         //!< the releases made, over all its lines (0: nothing was left on it)
};

class TowerFailureLog
{
public:
    /** Add a failure (the latest). The tower must not have failed before; names are single words. */
    void add (const TowerFailure& f);
    std::size_t size () const { return m_f.size(); }
    bool empty () const { return m_f.empty(); }
    const TowerFailure& record (std::size_t i) const { return m_f.at(i); }
    /** Whether the tower has failed. */
    bool failed (const std::string& tower) const;

    /** The header of tower_failures.dat: one row per failed tower. */
    static const char* log_header ();
    /** The row of failure i (whitespace separated, no newline). */
    std::string log_row (std::size_t i) const;

    /** The log as text: one record per failure in order, the time at 17 significant digits. */
    std::string state () const;
    /** Replace the log by what state() wrote. Returns an empty string, or the message for a malformed record. */
    std::string restore (const std::string& text);

private:
    std::vector<TowerFailure> m_f;
};

} // namespace erf_conductors

#endif
