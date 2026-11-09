#ifndef ERF_CAPACITIES_H
#define ERF_CAPACITIES_H

/**
 * CapacityMonitor: the demand on each structural component of the conductor lines and their towers
 * against its capacity, the peak demand-to-capacity ratio of each, and the first time each one's
 * demand exceeds its capacity (the first-failure log).
 *
 * Conventions: demands and capacities are nominal, in N (a member's are its axial force and its
 * design strength on the side that governs, tension or compression); a component exceeds its
 * capacity when its demand is larger than the capacity (ratio > 1). Times are ERF's (s). The monitor
 * knows nothing of MoorDyn or of the towers: Conductors hands it the demands after each coupling step.
 */

#include <cstddef>
#include <string>
#include <vector>

namespace erf_conductors {

/** What a monitored component is, which names its capacity in the inputs. */
enum class ComponentKind {
    Conductor,            //!< a span's conductor: its largest tension against conductors.<line>.rated_breaking_strength
    Insulator,            //!< a suspension string: the tension at its top against conductors.<line>.insulator_rated_load
    Hardware,             //!< a dead end or a tower attachment: the line's pull there against conductors.<line>.hardware_rated_load
    FootingUplift,        //!< a tower's footings: the largest leg uplift against conductors.<type>.allowable_uplift
    FootingCompression,   //!< a tower's footings: the largest leg compression against conductors.<type>.allowable_compression
    Member                //!< a frame tower's member: its axial force against its ASCE 10-15 design strength
};

/** The name of a kind as the logs write it: conductor, insulator, hardware, footing_uplift, footing_compression, member. */
const char* kind_name (ComponentKind kind);

/** One monitored component and what the run has seen of it. */
struct CapacityComponent
{
    ComponentKind kind{ComponentKind::Conductor};
    std::string name;               //!< unique over the monitor, e.g. S1_span2, S1_string1, S1_end_a, T3_footing_uplift, T3_m12
    double capacity{0.0};           //!< the capacity at the last record (N)
    double demand{0.0};             //!< the demand at the last record (N)
    bool seen{false};               //!< a demand has been recorded
    double peak_ratio{0.0};         //!< the largest demand / capacity recorded
    double peak_demand{0.0};        //!< the demand at that peak (N)
    double peak_capacity{0.0};      //!< the capacity at that peak (N)
    double peak_time{0.0};          //!< when it occurred (s)
    bool exceeded{false};           //!< the demand has exceeded the capacity at some record
    double first_time{0.0};         //!< the first record at which it did (s)
    double first_demand{0.0};       //!< the demand then (N)
    double first_capacity{0.0};     //!< the capacity then (N)
};

class CapacityMonitor
{
public:
    /** Add a component; its name must be new and not empty. Returns its index. */
    std::size_t add (ComponentKind kind, const std::string& name);

    /**
     * Record component c's demand and capacity at time t (both N; the capacity positive, the demand
     * >= 0, both finite, else an abort naming the component). Returns true when this is the
     * component's first exceedance: demand > capacity for the first time in the run.
     */
    bool record (std::size_t c, double demand, double capacity, double t);

    std::size_t size () const { return m_c.size(); }
    bool empty () const { return m_c.empty(); }
    const CapacityComponent& component (std::size_t c) const { return m_c.at(c); }
    const std::vector<CapacityComponent>& components () const { return m_c; }
    /**
     * The component that exceeded its capacity first in the run (the earliest first_time; on a tie,
     * the first added), or size() when none has.
     */
    std::size_t first_failure () const;
    /** Order component indices by their first exceedance time, keeping the given order among equal times. */
    void order_by_first_time (std::vector<std::size_t>& ids) const;

    /**
     * Write the table of every component to fname (truncating it): one CSV row each with its kind,
     * name, last demand and capacity, peak ratio with its demand, capacity and time, and whether and
     * when it first exceeded its capacity, with the demand and capacity then (empty fields when not).
     */
    void write_table (const std::string& fname) const;
    /** The header of the failure log: one row per component at its first exceedance. */
    static const char* failure_log_header ();
    /** The failure log's row for component c's first exceedance (whitespace separated, no newline). */
    std::string failure_row (std::size_t c) const;

    /**
     * The monitor's state as text: one record per component in order, the name first, every value
     * at 17 significant digits. restore() takes it back.
     */
    std::string state () const;
    /**
     * Take back what state() wrote for the components this monitor holds, matched by name and kind:
     * a component without a record starts afresh (a capacity given only on the restart), a record
     * without a component is dropped. Returns an empty string, or the message for a malformed record.
     * names_dropped receives the records without a component.
     */
    std::string restore (const std::string& text, std::vector<std::string>& names_dropped);

private:
    std::vector<CapacityComponent> m_c;
};

/**
 * The demand and capacity of a frame member's check on the side that governs: its tension and
 * tension capacity when tension / tension_capacity >= compression / compression_capacity, else its
 * compression and compression capacity (N). Capacities must be positive (else an abort): the steel's
 * strength reaches zero only at 1200 C (EN 1993-1-2), which the tower inputs refuse; a fire coupling
 * that lets the steel reach it must decide how a member without strength is reported first.
 */
void member_demand (double tension, double compression, double tension_capacity, double compression_capacity,
                    double& demand, double& capacity);

} // namespace erf_conductors

#endif
