#ifndef ERF_GROUND_CONTACT_H
#define ERF_GROUND_CONTACT_H

/**
 * GroundContactLog: which spans of the conductor lines lie on the ground, from when, and how much of
 * each at most (conductors.ground_contact).
 *
 * A node is on the ground when its height above the terrain under it is at most the conductor's
 * diameter. Conductors counts each span's nodes on the ground at the end of every anchor-level step
 * and hands the count here; the log keeps the first time each span touched the ground and the most
 * nodes it had on it. It knows nothing of MoorDyn.
 */

#include <cstddef>
#include <string>
#include <vector>

namespace erf_conductors {

/** One span and what the run has seen of it on the ground. */
struct GroundSpan
{
    std::string name;            //!< unique over the log, e.g. S1_span2
    int nodes{0};                //!< the span's nodes (both ends included)
    int now{0};                  //!< its nodes on the ground at the last record
    int most{0};                 //!< the most nodes it has had on the ground
    double most_time{0.0};       //!< the first record with that many (s)
    bool touched{false};         //!< a node has been on the ground at some record
    double first_time{0.0};      //!< the first record at which one was (s)
};

class GroundContactLog
{
public:
    /** Add a span of nodes nodes (at least 2); its name must be new and not empty. Returns its index. */
    std::size_t add (const std::string& name, int nodes);

    /**
     * Record that span s has on_ground of its nodes on the ground at time t (0 to its node count,
     * else an abort naming it). Returns true when this is its first touch of the ground in the run.
     */
    bool record (std::size_t s, int on_ground, double t);

    std::size_t size () const { return m_s.size(); }
    bool empty () const { return m_s.empty(); }
    const GroundSpan& span (std::size_t s) const { return m_s.at(s); }

    /**
     * Write every span to fname (truncating it): one CSV row each with its name, node count, nodes
     * on the ground now, the most and when, and when it first touched the ground (empty when never).
     */
    void write_table (const std::string& fname) const;

    /** The records as text (one line per span, 17 significant digits), for a checkpoint. */
    std::string state () const;
    /**
     * Take the records of the spans in text (written by state()), matched by name. Spans text holds
     * that this log does not, or with another node count, are skipped and named in names_dropped.
     * Returns an error message for a malformed line, else an empty string.
     */
    std::string restore (const std::string& text, std::vector<std::string>& names_dropped);

private:
    std::vector<GroundSpan> m_s;
};

} // namespace erf_conductors

#endif
