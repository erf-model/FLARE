#ifndef ERF_FARSITE_WIND_OUTPUT_H
#define ERF_FARSITE_WIND_OUTPUT_H

/**
 * @file ERF_FarsiteWindOutput.H
 * @brief Gridded winds for FARSITE and FlamMap, in the layout WindNinja writes.
 *
 * With erf.fire.farsite_winds.interval > 0 the fire layer writes, every
 * interval seconds of simulated time,
 *
 *   <dir>/<prefix>_<step>_vel.asc   wind speed (mph or km/h, 2 decimals)
 *   <dir>/<prefix>_<step>_ang.asc   direction the wind blows from, degrees
 *                                   clockwise from north (0 decimals)
 *   <dir>/<prefix>_<step>_vel.prj, _ang.prj   copies of prj_file, when given
 *
 * and appends one line per pair to <dir>/<prefix>.atm:
 *
 *   WINDS
 *   ENGLISH                      (or METRIC)
 *   MM DD HHMM <vel file> <ang file>
 *
 * This matches farsiteAtm::writeAtmFile and ascii_grid<T>::write_Grid in
 * WindNinja (firelab/windninja, src/ninja): ENGLISH is mph at 20 ft above the
 * vegetation and METRIC is km/h at 10 m, so the fire layer's reference wind
 * height (erf.fire.wind_ref_ht, or 6.096 m with use_per_fuel_wind_ht) must be
 * that height. File names in the .atm carry no directory.
 *
 * The wind written is fire_wind_ref: u and v at the reference height above
 * ground, before the wind adjustment factor and the terrain corrections. x is
 * taken as east and y as north. The grids are written by the I/O rank after a
 * gather of the fire plane; uncovered cells hold NODATA_value -9999.
 *
 * The write times follow from the simulated time alone (no checkpoint state):
 * a step whose start time is at or past the next multiple of the interval
 * writes the wind it extracted at that start time. A restart continues the
 * .atm file, first dropping the lines of steps at or after its own (a restart
 * from an older checkpoint); a run that starts at t = 0 truncates it.
 *
 * With coarsen > 1 the u, v block means are taken on the fire grid's ranks
 * when its boxes coarsen evenly, so only the output grid is gathered; the
 * I/O rank holds one copy of the output grid (or the fire plane) per write.
 */

#include <AMReX_REAL.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <cstdio>
#include <ctime>
#include <string>

#include <ERF_FireParams.H>

namespace farsite_winds {

/// metres per second to miles per hour (1 mile = 1609.344 m)
constexpr double mph_per_mps = 3600.0 / 1609.344;
/// metres per second to kilometres per hour
constexpr double kph_per_mps = 3.6;
/// the NODATA_value of the written grids (WindNinja's default)
constexpr double nodata = -9999.0;

/// The wind height the .atm units imply: 20 ft (6.096 m) or 10 m.
inline double unit_height_m (bool english) { return english ? 6.096 : 10.0; }

/// Wind speed in the output unit from the east (u) and north (v) components [m/s].
inline double speed_out (double u, double v, bool english)
{
    return std::sqrt(u * u + v * v) * (english ? mph_per_mps : kph_per_mps);
}

/**
 * Direction the wind blows from, degrees clockwise from north, in [0, 360).
 * A wind toward the north (v > 0) is 180, one toward the east (u > 0) is 270.
 * This equals WindNinja's -atan2(u, -v) wrapped to [0, 360). Calm air is 0.
 */
inline double direction_from_deg (double u, double v)
{
    if (u == 0.0 && v == 0.0) { return 0.0; }
    const double pi = 3.14159265358979323846;
    double d = std::atan2(-u, -v) * 180.0 / pi;
    if (d < 0.0)     { d += 360.0; }
    if (d >= 360.0)  { d -= 360.0; }
    return d + 0.0;   // + 0.0 turns -0.0 (u = +0, v < 0) into 0
}

/// The .atm header: "WINDS" and the unit line.
inline std::string atm_header (bool english)
{
    return std::string("WINDS\n") + (english ? "ENGLISH\n" : "METRIC\n");
}

/// The ESRI ASCII header, tab separated with six decimals, as WindNinja writes it.
inline std::string asc_header (int ncols, int nrows, double xll, double yll, double cellsize)
{
    char buf[512];
    std::snprintf(buf, sizeof(buf),
                  "ncols\t%d\nnrows\t%d\nxllcorner\t%lf\nyllcorner\t%lf\ncellsize\t%lf\nNODATA_value\t%lf\n",
                  ncols, nrows, xll, yll, cellsize, nodata);
    return std::string(buf);
}

/// The file root of one write: <prefix>_<step, six digits>.
inline std::string grid_root (const std::string& prefix, int step)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "_%06d", step);
    return prefix + buf;
}

/**
 * The .atm line of one write: "MM DD HHMM <vel> <ang>" in UTC, the time
 * rounded to the nearest minute.
 *
 * @param[in] epoch_s  start_datetime as epoch seconds plus the simulated time
 */
inline std::string atm_line (double epoch_s, const std::string& vel_file, const std::string& ang_file)
{
    const std::time_t t = static_cast<std::time_t>(std::llround(epoch_s / 60.0) * 60);  // nearest minute
    std::tm tm_utc{};
#if defined(_WIN32)
    gmtime_s(&tm_utc, &t);
#else
    gmtime_r(&t, &tm_utc);
#endif
    char buf[32];
    std::strftime(buf, sizeof(buf), "%m %d %H%M", &tm_utc);
    return std::string(buf) + " " + vel_file + " " + ang_file + "\n";
}

/**
 * Slack, in intervals, for deciding that time t is on a multiple: 1e-9, or
 * four rounding units of t in amrex::Real when that is larger (single
 * precision: one float ulp at t = 86400 s is 8e-3 s).
 */
inline double interval_slack (double t, double interval)
{
    const double eps = static_cast<double>(std::numeric_limits<amrex::Real>::epsilon());
    return std::max(1.0e-9, 4.0 * eps * std::abs(t) / interval);
}

/// Index of the last interval multiple at or before time t.
inline long long interval_index (double t, double interval)
{
    return static_cast<long long>(std::floor(t / interval + interval_slack(t, interval)));
}

/// A direction rounded to whole degrees for the 0-decimal grid, 360 written as 0.
inline double direction_whole_deg (double d)
{
    const double r = std::round(d);
    return ((r >= 360.0) ? r - 360.0 : r) + 0.0;
}

/**
 * The step number of an .atm data line "MM DD HHMM <prefix>_<step>_vel.asc ...",
 * or -1 when the speed file is not named that way.
 */
inline int atm_line_step (const std::string& line, const std::string& prefix)
{
    std::istringstream ss(line);
    std::string mm, dd, hhmm, vel;
    if (!(ss >> mm >> dd >> hhmm >> vel)) { return -1; }
    const std::string head = prefix + "_";
    const std::string tail = "_vel.asc";
    if (vel.size() <= head.size() + tail.size() || vel.compare(0, head.size(), head) != 0
        || vel.compare(vel.size() - tail.size(), tail.size(), tail) != 0) { return -1; }
    const std::string digits = vel.substr(head.size(), vel.size() - head.size() - tail.size());
    if (digits.empty() || digits.find_first_not_of("0123456789") != std::string::npos) { return -1; }
    return std::stoi(digits);
}

} // namespace farsite_winds

/**
 * @brief Writes the reference wind as FARSITE gridded winds (see the file comment).
 */
class FarsiteWindWriter
{
public:
    /**
     * @param[in] params     erf.fire.farsite_winds.* (must be enabled)
     * @param[in] fire_geom  the fire grid geometry
     * @param[in] wind_ht    the reference wind height [m] (FireParams::sampled_wind_ht)
     *
     * Aborts, naming the key, when the wind height does not match the units,
     * when start_datetime is missing or malformed, when the fire cells are not
     * square, when coarsen does not divide the fire grid, or when prj_file
     * cannot be read. Creates the output directory.
     */
    FarsiteWindWriter (const FireParams::FarsiteWindsParams& params,
                       const amrex::Geometry& fire_geom,
                       amrex::Real wind_ht);

    /**
     * @brief Write the grids and the .atm line when this step is due.
     *
     * @param[in] wind_ref  u, v at the reference height on the fire grid [m/s]
     * @param[in] time      the step's start time [s], the time wind_ref was extracted at
     * @param[in] step      the fire step number (file names)
     */
    void maybe_write (const amrex::MultiFab& wind_ref, amrex::Real time, int step);

    /// The epoch minute an .atm line is stamped with (the time rounded to the nearest minute).
    static long long stamp_minute (double epoch_s) { return std::llround(epoch_s / 60.0); }

    /// Number of grid pairs written by this run (for tests and the log).
    int n_written () const { return m_n_written; }

private:
    void write_grids (const amrex::MultiFab& wind_ref, amrex::Real time, int step);

    FireParams::FarsiteWindsParams m_p;
    amrex::Geometry m_geom;
    bool        m_english   = true;
    double      m_start_epoch = 0.0;
    double      m_xll       = 0.0;
    double      m_yll       = 0.0;
    double      m_cellsize  = 0.0;
    int         m_ncols     = 0;
    int         m_nrows     = 0;
    std::string m_prj_text;
    long long   m_last_index = 0;
    bool        m_first_call = true;
    long long   m_last_minute = std::numeric_limits<long long>::min();  ///< stamp of the last .atm line
    double      m_prev_time   = -1.0;   ///< the previous call's time, for the dt > interval warning
    bool        m_warned_dt   = false;
    int         m_n_written  = 0;
};

#endif
