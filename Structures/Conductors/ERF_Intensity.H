#ifndef ERF_INTENSITY_H
#define ERF_INTENSITY_H

/**
 * IntensityRecorder: the wind intensity measure beside each structural component's peak demand.
 *
 * An intensity point is a place on the ground plan (a span's middle, a line's attachment, a tower's base) at
 * which the run samples the horizontal wind at a fixed height above the ground. The recorder keeps, per point,
 * the wind speed averaged over the last T seconds (the averaging time: the T-second average of the wind the flow
 * resolves, the mean wind of a RANS flow), its largest value in the run once the window is full, with the time and
 * direction of that maximum, and, per component of a CapacityMonitor linked to a point, the averaged wind at
 * the component's peak demand-to-capacity ratio and at its first exceedance.
 *
 * The average at the end of a step is over the window (t - T, t], t the step's end. A sample taken for a step
 * of length dt covers the last dt of the window, the one before it the dt before that, and so on; the oldest
 * sample that reaches past the window's start counts with the part inside it:
 *   S = sum_i w_i s_i / sum_i w_i,   w_i = the length of sample i's step inside the window.
 * Before the run is T old the window holds what has been sampled. T = 0 keeps the last sample alone. The speed
 * s is the horizontal one, sqrt(u^2 + v^2), averaged as a speed; the direction is that of the averaged vector
 * (u, v), in degrees from +x towards +y (the direction the wind blows towards), in (-180, 180]; 0 when the
 * averaged vector is zero. The turbulent kinetic energy k, when the closure carries one, is averaged like s.
 *
 * SI units (m, s). The recorder knows nothing of ERF's mesh: Conductors hands it the samples.
 */

#include <cstddef>
#include <deque>
#include <string>
#include <vector>

namespace erf_conductors {

class CapacityMonitor;

/** One sample of the wind at an intensity point: the horizontal velocity (m/s) and k (m^2/s^2; unused without one). */
struct IntensitySample
{
    double u{0.0};
    double v{0.0};
    double k{0.0};
};

/** The averaged wind at a point: its speed (m/s), direction (degrees from +x towards +y) and k (m^2/s^2). */
struct IntensityValue
{
    double speed{0.0};
    double direction{0.0};
    double k{0.0};
    double time{0.0};      //!< the end of the step it belongs to (s)
    bool set{false};       //!< false until a value has been taken
};

/** An intensity point and what the run has seen at it. */
struct IntensityPoint
{
    std::string name;      //!< unique over the recorder, e.g. S1_span2, S1_end_a, S1_tower1, tower_T3
    double x{0.0};         //!< where it stands (m)
    double y{0.0};
    double z{0.0};         //!< the height it is sampled at: the ground there plus the recorder's height (m)
    IntensityValue now;    //!< the averaged wind at the last step
    IntensityValue max;    //!< the largest averaged speed of the run, with its direction, k and time
};

class IntensityRecorder
{
public:
    /**
     * Start afresh with the height above the ground the wind is sampled at (m) and the averaging time (s);
     * with_k says whether the samples carry a turbulent kinetic energy. Returns an empty string, or the
     * message naming the input that is out of range: conductors.intensity_height in (0, 500],
     * conductors.intensity_averaging_time in [0, 600].
     */
    std::string configure (double height, double averaging_time, bool with_k);
    /** The message configure() would return for these values, without changing the recorder. */
    static std::string check (double height, double averaging_time);

    double height () const { return m_height; }
    /**
     * The time (s) from which a point's largest value in the run is looked for: an average stamped earlier is not
     * compared (conductors.capacity_check_start, so that a run's start-up is in neither the peaks nor the maxima).
     * The averaging window is filled from the first sample whatever this is. Default 0.
     */
    void set_maximum_start (double t) { m_max_start = t; }
    double averaging_time () const { return m_T; }
    bool with_k () const { return m_with_k; }

    /** Add a point at (x, y), sampled at height z (m); its name must be new and one word without commas. Returns its index. */
    std::size_t add_point (const std::string& name, double x, double y, double z);
    std::size_t num_points () const { return m_p.size(); }
    const IntensityPoint& point (std::size_t p) const { return m_p.at(p); }

    /**
     * Component c of the capacity monitor (its index there) takes its intensity from point p. A component
     * can be linked once; components() must be called with the monitor's size before the first link.
     */
    void components (std::size_t n);
    void link (std::size_t c, std::size_t p);
    /** The point component c is linked to, or num_points() when it has none. */
    std::size_t point_of (std::size_t c) const;

    /**
     * Take one sample per point (in the points' order) for the step of length dt (> 0) that ends at t, and
     * update every point's averaged wind and, once the steps sampled cover the averaging time, its run maximum. Aborts on a sample that is not finite, a
     * negative k, a dt that is not positive or a time that does not advance.
     */
    void sample (const std::vector<IntensitySample>& s, double dt, double t);
    /** The number of steps sample() has taken. */
    std::size_t num_samples () const { return m_nsamples; }

    /** Component c's demand-to-capacity ratio reached a new peak at this step: keep its point's averaged wind. */
    void at_peak (std::size_t c);
    /** Component c exceeded its capacity for the first time at this step: keep its point's averaged wind. */
    void at_first (std::size_t c);
    const IntensityValue& peak_value (std::size_t c) const { return m_peak.at(c); }
    const IntensityValue& first_value (std::size_t c) const { return m_first.at(c); }

    /**
     * Write the table of the points to fname (truncating it): one CSV row each with its name, position (x, y and
     * the elevation z it is sampled at), sampling height above the ground, the averaged wind now and the run's maximum with its time.
     */
    void write_points (const std::string& fname) const;
    /**
     * Write the run's demand table to fname (truncating it): one CSV row per component of cap with its
     * kind and name, its peak ratio with the demand, capacity and time of the peak, its intensity point,
     * the averaged wind at the peak with the end of the step it was sampled for, the point's run maximum
     * with its time, and whether and when the component first exceeded its capacity with the averaged wind
     * then and its step's end. Fields that do not exist (a
     * component never recorded or never exceeded, one without a point, k without a closure) are empty.
     * cap must be the monitor the links were made for.
     */
    void write_demand_table (const std::string& fname, const CapacityMonitor& cap) const;

    /**
     * The recorder's state as text, every value at 17 significant digits: the configuration, each point
     * with its samples inside the window, and each linked component's kept winds by the component's name
     * in cap. restore() takes it back.
     */
    std::string state (const CapacityMonitor& cap) const;
    /**
     * Take back what state() wrote, for the points and components this recorder holds (matched by name): a
     * point or component without a record starts afresh, a record without one is dropped (names_dropped).
     * The height, the averaging time and with_k must be the checkpoint's, else the message says which
     * differs. Returns an empty string, or the message for a mismatch or a malformed record.
     */
    std::string restore (const std::string& text, const CapacityMonitor& cap, std::vector<std::string>& names_dropped);

private:
    double m_height{10.0};
    double m_T{3.0};
    bool m_with_k{false};
    double m_max_start{0.0};   //!< no run maximum from an average stamped before this (s)
    std::size_t m_nsamples{0};
    double m_time{0.0};                               //!< the end of the last sampled step (s)
    bool m_has_time{false};
    std::vector<IntensityPoint> m_p;
    std::deque<double> m_dt;                          //!< the step lengths inside the window, newest first
    std::vector<std::deque<IntensitySample>> m_hist;  //!< per point, the samples of those steps, newest first
    std::vector<std::size_t> m_link;                  //!< per component, its point (num_points() = none)
    std::vector<IntensityValue> m_peak;               //!< per component, the averaged wind at its peak ratio
    std::vector<IntensityValue> m_first;              //!< per component, the averaged wind at its first exceedance

    void average (std::size_t p, IntensityValue& out) const;
};

} // namespace erf_conductors

#endif
