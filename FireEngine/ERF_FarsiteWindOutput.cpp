#include "ERF_FarsiteWindOutput.H"
#include "ERF_EpochTime.H"

#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFabUtil.H>

#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

using namespace amrex;

FarsiteWindWriter::FarsiteWindWriter (const FireParams::FarsiteWindsParams& params,
                                      const Geometry& fire_geom,
                                      Real wind_ht)
    : m_p(params), m_geom(fire_geom)
{
    AMREX_ALWAYS_ASSERT(m_p.enabled());
    m_english = (m_p.units == "english");

    // ENGLISH means 20 ft and METRIC 10 m in FARSITE's .atm (farsiteAtm.cpp in
    // WindNinja); a wind at another height would be read as if at that one.
    // wind_ref_ht's default 6.1 m is 20 ft to 4 mm, hence the 1 cm slack.
    const double h_unit = farsite_winds::unit_height_m(m_english);
    if (std::abs(static_cast<double>(wind_ht) - h_unit) > 0.01) {
        Abort("[FIRE] erf.fire.farsite_winds.units = " + m_p.units + " means a wind at "
              + std::to_string(h_unit) + " m, but the reference wind is sampled at "
              + std::to_string(wind_ht) + " m: set erf.fire.wind_ref_ht = "
              + (m_english ? "6.096 (or use_per_fuel_wind_ht)" : "10 without use_per_fuel_wind_ht")
              + ", or the other units");
    }

    // The .atm lines carry month, day and time: the run's start_datetime
    // (the top-level key ERF reads, YYYY-MM-DD HH:MM[:SS], UTC).
    {
        ParmParse pp_top;
        std::string start;
        if (!pp_top.query("start_datetime", start)) {
            Abort("[FIRE] erf.fire.farsite_winds.interval > 0 needs start_datetime"
                  " (YYYY-MM-DD HH:MM:SS, UTC) for the dates in the .atm file");
        }
        if (start.length() == 16) { start += ":00"; }
        const std::time_t e = (start.length() == 19) ? getEpochTime(start, "%Y-%m-%d %H:%M:%S") : -1;
        if (e == static_cast<std::time_t>(-1)) {
            Abort("[FIRE] start_datetime = \"" + start + "\" is not YYYY-MM-DD HH:MM[:SS]"
                  " (needed by erf.fire.farsite_winds)");
        }
        m_start_epoch = static_cast<double>(e);
    }

    // ESRI ASCII grids have one cell size
    const Real dx = m_geom.CellSize(0);
    const Real dy = m_geom.CellSize(1);
    if (std::abs(dx - dy) > Real(1.0e-6) * dx) {
        Abort("[FIRE] erf.fire.farsite_winds needs square fire cells (ESRI ASCII grids have one"
              " cellsize), got dx = " + std::to_string(dx) + " m, dy = " + std::to_string(dy) + " m");
    }
    const Box& dom = m_geom.Domain();
    const int nx = dom.length(0);
    const int ny = dom.length(1);
    if (nx % m_p.coarsen != 0 || ny % m_p.coarsen != 0) {
        Abort("[FIRE] erf.fire.farsite_winds.coarsen = " + std::to_string(m_p.coarsen)
              + " must divide the fire grid, " + std::to_string(nx) + " x " + std::to_string(ny) + " cells");
    }
    m_ncols    = nx / m_p.coarsen;
    m_nrows    = ny / m_p.coarsen;
    m_cellsize = static_cast<double>(dx) * m_p.coarsen;
    m_xll      = m_p.has_xllcorner ? static_cast<double>(m_p.xllcorner) : static_cast<double>(m_geom.ProbLo(0));
    m_yll      = m_p.has_yllcorner ? static_cast<double>(m_p.yllcorner) : static_cast<double>(m_geom.ProbLo(1));

    // The projection text, read once by the I/O rank (only it writes)
    if (!m_p.prj_file.empty()) {
        int ok = 1;
        if (ParallelDescriptor::IOProcessor()) {
            std::ifstream f(m_p.prj_file);
            ok = f.good() ? 1 : 0;
            if (ok) { std::stringstream ss; ss << f.rdbuf(); m_prj_text = ss.str(); }
        }
        ParallelDescriptor::Bcast(&ok, 1, ParallelDescriptor::IOProcessorNumber());
        if (ok == 0) {
            Abort("[FIRE] erf.fire.farsite_winds.prj_file = \"" + m_p.prj_file + "\" cannot be read");
        }
    }

    if (ParallelDescriptor::IOProcessor()) {
        if (!UtilCreateDirectory(m_p.dir, 0755)) { CreateDirectoryFailed(m_p.dir); }
    }
    ParallelDescriptor::Barrier();

    Print() << "[FIRE] FARSITE gridded winds every " << m_p.interval << " s to " << m_p.dir << "/"
            << m_p.prefix << ".atm (" << (m_english ? "ENGLISH, mph at 20 ft" : "METRIC, km/h at 10 m")
            << "), " << m_ncols << " x " << m_nrows << " cells of " << m_cellsize << " m, lower-left "
            << m_xll << ", " << m_yll << (m_p.prj_file.empty() ? ", no .prj" : "") << "\n";
}

void
FarsiteWindWriter::maybe_write (const MultiFab& wind_ref, Real time, int step)
{
    const long long k = farsite_winds::interval_index(time, m_p.interval);
    if (m_first_call) {
        m_first_call = false;
        // A start on a multiple (t = 0, or a restart there) writes now; a start
        // between two waits for the next one. The step that ends on a multiple
        // does not write: the step starting there does, with the wind of that time.
        const double slack = farsite_winds::interval_slack(time, m_p.interval) * m_p.interval;
        const bool on_multiple = std::abs(double(time) - double(k) * m_p.interval) <= slack
                                 || std::abs(time) <= Real(1.0e-12);
        m_last_index = on_multiple ? k - 1 : k;
    }
    // A step longer than the interval skips multiples: say so once
    if (!m_warned_dt && m_prev_time >= 0.0 && double(time) - m_prev_time > m_p.interval) {
        m_warned_dt = true;
        Print() << "[FIRE] WARNING: the time step (" << double(time) - m_prev_time << " s) exceeds"
                << " erf.fire.farsite_winds.interval = " << m_p.interval
                << " s: some grids are skipped, and each is stamped with its step's start time\n";
    }
    m_prev_time = double(time);
    if (k <= m_last_index) { return; }
    m_last_index = k;

    // FARSITE reads the .atm times in increasing order. A stamp is the step's
    // start time rounded to the minute, which lags the multiple by up to dt:
    // with dt > interval - 30 s two stamps can fall in the same minute.
    const long long minute = stamp_minute(m_start_epoch + double(time));
    if (minute <= m_last_minute) {
        Abort("[FIRE] erf.fire.farsite_winds: two .atm lines would carry the same minute (t = "
              + std::to_string(double(time)) + " s): raise erf.fire.farsite_winds.interval = "
              + std::to_string(m_p.interval) + " s, or keep the time step below interval - 30 s");
    }
    m_last_minute = minute;
    write_grids(wind_ref, time, step);
}

void
FarsiteWindWriter::write_grids (const MultiFab& wind_ref, Real time, int step)
{
    // With coarsen > 1 and fire boxes that coarsen evenly, average u and v on
    // the fire grid's own ranks (and device) first, so only the output grid
    // is gathered; otherwise gather the fire plane and average on the I/O rank.
    const int c = m_p.coarsen;
    const IntVect ratio(c, c, 1);
    const MultiFab* src = &wind_ref;
    MultiFab crse;
    int c_host = c;
    Box dom = m_geom.Domain();
    if (c > 1 && wind_ref.boxArray().coarsenable(ratio)) {
        BoxArray ba_c = wind_ref.boxArray();
        ba_c.coarsen(ratio);
        crse.define(ba_c, wind_ref.DistributionMap(), 2, 0);
        amrex::average_down(wind_ref, crse, 0, 2, ratio);
        src    = &crse;
        c_host = 1;
        dom.coarsen(ratio);
    }

    // Gather u, v onto one box owned by the I/O rank, in host-readable memory.
    // Cells no fire box covers keep the sentinel and are written as NODATA.
    constexpr Real sentinel = Real(1.0e30);
    BoxArray ba_one(dom);
    Vector<int> owner(1, ParallelDescriptor::IOProcessorNumber());
    DistributionMapping dm_one(owner);
    MultiFab one(ba_one, dm_one, 2, 0, MFInfo().SetArena(The_Pinned_Arena()));
    one.setVal(sentinel);
    one.ParallelCopy(*src, 0, 0, 2);
    Gpu::streamSynchronize();

    const std::string root  = farsite_winds::grid_root(m_p.prefix, step);
    const std::string vel_f = root + "_vel.asc";
    const std::string ang_f = root + "_ang.asc";

    if (ParallelDescriptor::IOProcessor()) {
        const auto a = one.const_array(0);
        const int i0 = dom.smallEnd(0);
        const int j0 = dom.smallEnd(1);
        const int kk = dom.smallEnd(2);

        // block averages of u and v (vector mean), NODATA when any cell is uncovered
        std::vector<double> spd(static_cast<std::size_t>(m_ncols) * m_nrows);
        std::vector<double> ang(spd.size());
        for (int J = 0; J < m_nrows; ++J) {
            for (int I = 0; I < m_ncols; ++I) {
                double su = 0.0, sv = 0.0;
                bool covered = true;
                for (int jj = 0; jj < c_host && covered; ++jj) {
                    for (int ii = 0; ii < c_host; ++ii) {
                        const Real u = a(i0 + I * c_host + ii, j0 + J * c_host + jj, kk, 0);
                        const Real v = a(i0 + I * c_host + ii, j0 + J * c_host + jj, kk, 1);
                        // uncovered (sentinel) or non-finite: NODATA
                        if (!(std::isfinite(u) && std::isfinite(v)) || u >= sentinel || v >= sentinel) { covered = false; break; }
                        su += u; sv += v;
                    }
                }
                const std::size_t n = static_cast<std::size_t>(J) * m_ncols + I;
                if (covered) {
                    su /= double(c_host * c_host); sv /= double(c_host * c_host);
                    spd[n] = farsite_winds::speed_out(su, sv, m_english);
                    // written with 0 decimals: round here so 359.6 is 0, not 360
                    ang[n] = farsite_winds::direction_whole_deg(farsite_winds::direction_from_deg(su, sv));
                } else {
                    spd[n] = farsite_winds::nodata;
                    ang[n] = farsite_winds::nodata;
                }
            }
        }

        // One grid: header, then rows from north to south, each value followed
        // by a tab (WindNinja's write_fpGridData), speed 2 decimals, angle 0.
        auto write_asc = [&](const std::string& name, const std::vector<double>& val, int decimals) {
            const std::string path = m_p.dir + "/" + name;
            std::FILE* f = std::fopen(path.c_str(), "w");
            if (f == nullptr) { FileOpenFailed(path); }
            const std::string hdr = farsite_winds::asc_header(m_ncols, m_nrows, m_xll, m_yll, m_cellsize);
            std::fputs(hdr.c_str(), f);
            for (int J = m_nrows - 1; J >= 0; --J) {
                for (int I = 0; I < m_ncols; ++I) {
                    std::fprintf(f, "%.*lf\t", decimals, val[static_cast<std::size_t>(J) * m_ncols + I]);
                }
                std::fputs("\n", f);
            }
            if (std::fclose(f) != 0) { FileOpenFailed(path); }
            if (!m_prj_text.empty()) {
                const std::string prj = m_p.dir + "/" + name.substr(0, name.size() - 4) + ".prj";
                std::ofstream pf(prj, std::ios::trunc);
                if (!pf.good()) { FileOpenFailed(prj); }
                pf << m_prj_text;
            }
        };
        write_asc(vel_f, spd, 2);
        write_asc(ang_f, ang, 0);

        // The .atm: a run that starts at t = 0 begins a new file; any other
        // start (a restart) appends to the file it finds, or begins one.
        const std::string atm = m_p.dir + "/" + m_p.prefix + ".atm";
        const bool fresh = (m_n_written == 0) && (std::abs(time) <= Real(1.0e-12) || !std::ifstream(atm).good());
        if (!fresh && m_n_written == 0) {
            // A restart from a checkpoint older than the file's last line would
            // append times that go backwards: keep only the lines of earlier steps.
            std::ifstream in(atm);
            std::string kept, line;
            for (int n = 0; std::getline(in, line); ++n) {
                if (n == 1 && line != (m_english ? "ENGLISH" : "METRIC")) {
                    Abort("[FIRE] " + atm + " is " + line + " but erf.fire.farsite_winds.units = "
                          + m_p.units + ": a restart keeps the units, or set another"
                          " erf.fire.farsite_winds.dir or prefix");
                }
                const int s_line = (n < 2) ? -1 : farsite_winds::atm_line_step(line, m_p.prefix);
                if (n < 2 || (s_line >= 0 && s_line < step)) { kept += line + "\n"; }
            }
            in.close();
            std::ofstream out(atm, std::ios::trunc);
            if (!out.good()) { FileOpenFailed(atm); }
            out << kept;
        }
        std::ofstream af(atm, fresh ? std::ios::trunc : std::ios::app);
        if (!af.good()) { FileOpenFailed(atm); }
        if (fresh) { af << farsite_winds::atm_header(m_english); }
        af << farsite_winds::atm_line(m_start_epoch + static_cast<double>(time), vel_f, ang_f);
        if (!af.good()) { FileOpenFailed(atm); }
    }
    ++m_n_written;
    Print() << "[FIRE] FARSITE winds t=" << time << " s: " << m_p.dir << "/" << vel_f << ", " << ang_f << "\n";
}
