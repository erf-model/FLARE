#include <gtest/gtest.h>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Math.H>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "ERF_FarsiteWindOutput.H"
#include "ERF_EpochTime.H"
#include "ERF_GTestTempDir.H"

/**
 * @file ERF_GTestFarsiteWindOutput.cpp
 * @brief erf.fire.farsite_winds: the WindNinja/FARSITE layout written by
 *        FarsiteWindWriter, checked against hand values: the direction
 *        convention, the unit factors, the .atm time stamps, the ESRI header,
 *        north-first row order, block averaging, the write cadence and the
 *        restart append.
 */

using namespace amrex;

namespace {

std::vector<std::string> read_lines (const std::filesystem::path& p)
{
    std::ifstream f(p);
    std::vector<std::string> out;
    for (std::string l; std::getline(f, l); ) { out.push_back(l); }
    return out;
}

std::vector<double> row_values (const std::string& line)
{
    std::istringstream ss(line);
    std::vector<double> v;
    for (double x; ss >> x; ) { v.push_back(x); }
    return v;
}

/// start_datetime set for one test and removed on exit
struct ScopedStart {
    explicit ScopedStart (const std::string& s) { ParmParse pp; pp.add("start_datetime", s); }
    ~ScopedStart () { ParmParse pp; pp.remove("start_datetime"); }
};

/// u = i (cell index east), v = 10 j (north). A free function: nvcc rejects an
/// extended __device__ lambda inside a TEST body (a private member function).
void fill_test_wind (MultiFab& wind)
{
    for (MFIter mfi(wind); mfi.isValid(); ++mfi) {
        const auto a = wind.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            a(i, j, k, 0) = Real(i);
            a(i, j, k, 1) = Real(10 * j);
        });
    }
}

} // namespace

// Wind toward north blows FROM the south (180), toward east FROM the west
// (270); WindNinja's -atan2(u, -v) wrapped to [0, 360).
TEST(FarsiteWinds, DirectionIsWhereTheWindComesFrom)
{
    using farsite_winds::direction_from_deg;
    EXPECT_NEAR(direction_from_deg( 0.0,  1.0), 180.0, 1e-12);
    EXPECT_NEAR(direction_from_deg( 1.0,  0.0), 270.0, 1e-12);
    EXPECT_NEAR(direction_from_deg( 0.0, -1.0),   0.0, 1e-12);
    EXPECT_NEAR(direction_from_deg(-1.0,  0.0),  90.0, 1e-12);
    EXPECT_NEAR(direction_from_deg( 3.0,  4.0), 180.0 + std::atan2(3.0, 4.0) * 180.0 / amrex::Math::pi<double>(), 1e-12);
    EXPECT_EQ  (direction_from_deg( 0.0,  0.0),   0.0);
    // wraps: from just west of north is near 360, never 360
    const double d = direction_from_deg(1.0e-12, -1.0);
    EXPECT_GE(d, 0.0);
    EXPECT_LT(d, 360.0);
}

// No "-0" and no "360" in the 0-decimal direction grid
TEST(FarsiteWinds, DirectionWrittenWithoutMinusZeroOr360)
{
    using farsite_winds::direction_from_deg;
    using farsite_winds::direction_whole_deg;
    EXPECT_FALSE(std::signbit(direction_from_deg(0.0, -1.0)));
    EXPECT_FALSE(std::signbit(direction_from_deg(-0.0, -1.0)));
    EXPECT_EQ(direction_whole_deg(359.6), 0.0);
    EXPECT_EQ(direction_whole_deg(359.4), 359.0);
    EXPECT_EQ(direction_whole_deg(0.4),   0.0);
    EXPECT_FALSE(std::signbit(direction_whole_deg(-0.0)));
    char buf[16];
    std::snprintf(buf, sizeof(buf), "%.0lf", direction_whole_deg(direction_from_deg(0.0, -1.0)));
    EXPECT_STREQ(buf, "0");
}

// The step of an .atm line, read back by the restart trim
TEST(FarsiteWinds, AtmLineStep)
{
    using farsite_winds::atm_line_step;
    EXPECT_EQ(atm_line_step("12 30 1101 flare_000061_vel.asc flare_000061_ang.asc", "flare"), 61);
    EXPECT_EQ(atm_line_step("12 30 1101 other_000061_vel.asc other_000061_ang.asc", "flare"), -1);
    EXPECT_EQ(atm_line_step("12 30 1101 flare_x_vel.asc a", "flare"), -1);
    EXPECT_EQ(atm_line_step("12 30", "flare"), -1);
}

// ENGLISH is mph (1 mile = 1609.344 m), METRIC is km/h
TEST(FarsiteWinds, SpeedUnits)
{
    EXPECT_NEAR(farsite_winds::speed_out(3.0, 4.0, true),  5.0 * 3600.0 / 1609.344, 1e-12);
    EXPECT_NEAR(farsite_winds::speed_out(3.0, 4.0, false), 18.0, 1e-12);
    EXPECT_DOUBLE_EQ(farsite_winds::unit_height_m(true),  6.096);
    EXPECT_DOUBLE_EQ(farsite_winds::unit_height_m(false), 10.0);
}

// "MM DD HHMM vel ang", UTC, rounded to the nearest minute
TEST(FarsiteWinds, AtmLineTimeStamp)
{
    const double e0 = static_cast<double>(getEpochTime("2021-12-30 11:00:00", "%Y-%m-%d %H:%M:%S"));
    EXPECT_EQ(farsite_winds::atm_line(e0, "a_vel.asc", "a_ang.asc"), "12 30 1100 a_vel.asc a_ang.asc\n");
    EXPECT_EQ(farsite_winds::atm_line(e0 + 3629.0, "v", "a"), "12 30 1200 v a\n");
    EXPECT_EQ(farsite_winds::atm_line(e0 + 3630.0, "v", "a"), "12 30 1201 v a\n");
    // across midnight and the year end
    EXPECT_EQ(farsite_winds::atm_line(e0 + 2.0 * 86400.0 - 11.0 * 3600.0, "v", "a"), "01 01 0000 v a\n");
    EXPECT_EQ(farsite_winds::atm_header(true),  "WINDS\nENGLISH\n");
    EXPECT_EQ(farsite_winds::atm_header(false), "WINDS\nMETRIC\n");
}

TEST(FarsiteWinds, AsciiHeaderAndNames)
{
    EXPECT_EQ(farsite_winds::asc_header(4, 2, 100.0, 200.5, 20.0),
              "ncols\t4\nnrows\t2\nxllcorner\t100.000000\nyllcorner\t200.500000\n"
              "cellsize\t20.000000\nNODATA_value\t-9999.000000\n");
    EXPECT_EQ(farsite_winds::grid_root("flare", 42), "flare_000042");
}

TEST(FarsiteWinds, IntervalIndex)
{
    EXPECT_EQ(farsite_winds::interval_index(0.0,   60.0), 0);
    EXPECT_EQ(farsite_winds::interval_index(59.99, 60.0), 0);
    EXPECT_EQ(farsite_winds::interval_index(60.0,  60.0), 1);
    // 0.1 summed 600 times lands a hair under 60
    double t = 0.0;
    for (int n = 0; n < 600; ++n) { t += 0.1; }
    EXPECT_EQ(farsite_winds::interval_index(t, 60.0), 1);
    // single precision: a time one float ulp under a multiple still counts
    const float tf = std::nextafter(86400.0f, 0.0f);
    if (sizeof(Real) == sizeof(float)) {
        EXPECT_EQ(farsite_winds::interval_index(tf, 60.0), 1440);
    }
}

// The writer end to end on an 8 x 4 fire grid of 10 m cells split in two
// boxes: u = i (cell index east), v = 10 j (north), coarsen 2. Output cell
// (I, J) averages fire cells 2I..2I+1, 2J..2J+1: u = 2I + 0.5, v = 10 (2J + 0.5).
// The first data row must be the northern one (J = 1).
TEST(FarsiteWinds, WriterGridsAtmAndCadence)
{
    ScopedStart start("2021-12-30 11:00");
    const auto dir = erf_gtest_temp_path("farsite_winds");

    Box dom(IntVect(0, 0, 0), IntVect(7, 3, 0));
    RealBox rb({100.0, 200.0, 0.0}, {180.0, 240.0, 1.0});
    Array<int, AMREX_SPACEDIM> per{0, 0, 0};
    Geometry geom(dom, rb, CoordSys::cartesian, per);
    BoxArray ba(dom);
    ba.maxSize(IntVect(4, 4, 1));
    DistributionMapping dm(ba);
    MultiFab wind(ba, dm, 2, 0);
    fill_test_wind(wind);

    FireParams::FarsiteWindsParams p;
    p.interval = 60.0;
    p.units    = "metric";
    p.dir      = dir.string();
    p.prefix   = "t";
    p.coarsen  = 2;
    {
        FarsiteWindWriter w(p, geom, 10.0);
        // t = 0 writes; 30 and 59.99 do not; the step starting at 60.5 writes
        // (the wind of its start); 119 does not; 121 writes
        int step = 0;
        for (Real t : std::vector<Real>{Real(0.0), Real(30.0), Real(59.99), Real(60.5), Real(119.0), Real(121.0)}) {
            w.maybe_write(wind, t, ++step);
        }
        EXPECT_EQ(w.n_written(), 3);
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor()) {
        const auto atm = read_lines(dir / "t.atm");
        ASSERT_EQ(atm.size(), 5u);
        EXPECT_EQ(atm[0], "WINDS");
        EXPECT_EQ(atm[1], "METRIC");
        EXPECT_EQ(atm[2], "12 30 1100 t_000001_vel.asc t_000001_ang.asc");
        EXPECT_EQ(atm[3], "12 30 1101 t_000004_vel.asc t_000004_ang.asc");
        EXPECT_EQ(atm[4], "12 30 1102 t_000006_vel.asc t_000006_ang.asc");

        const auto vel = read_lines(dir / "t_000001_vel.asc");
        const auto ang = read_lines(dir / "t_000001_ang.asc");
        ASSERT_EQ(vel.size(), 8u);
        EXPECT_EQ(vel[0], "ncols\t4");
        EXPECT_EQ(vel[1], "nrows\t2");
        EXPECT_EQ(vel[2], "xllcorner\t100.000000");
        EXPECT_EQ(vel[3], "yllcorner\t200.000000");
        EXPECT_EQ(vel[4], "cellsize\t20.000000");
        for (int r = 0; r < 2; ++r) {
            const int J = 1 - r;   // north row first
            const auto sv = row_values(vel[6 + r]);
            const auto sa = row_values(ang[6 + r]);
            ASSERT_EQ(sv.size(), 4u);
            ASSERT_EQ(sa.size(), 4u);
            for (int I = 0; I < 4; ++I) {
                const double u = 2.0 * I + 0.5;
                const double v = 10.0 * (2.0 * J + 0.5);
                EXPECT_NEAR(sv[I], std::sqrt(u * u + v * v) * 3.6, 0.005 + 1e-9) << "row " << r << " col " << I;
                EXPECT_NEAR(sa[I], std::round(farsite_winds::direction_from_deg(u, v)), 1e-9)
                    << "row " << r << " col " << I;
            }
        }
        // values end with a tab, as WindNinja writes them
        EXPECT_EQ(vel[6].back(), '\t');
        EXPECT_FALSE(std::filesystem::exists(dir / "t_000001_vel.prj"));
    }

    // A restart from the checkpoint of step 3 (t = 60, a multiple) drops the
    // lines of steps 4 and 6 written after it, writes step 4 again, and goes on
    {
        FarsiteWindWriter w(p, geom, 10.0);
        w.maybe_write(wind, Real(60.0), 4);
        w.maybe_write(wind, Real(90.0), 5);
        w.maybe_write(wind, Real(120.0), 6);
        EXPECT_EQ(w.n_written(), 2);
    }
    ParallelDescriptor::Barrier();
    if (ParallelDescriptor::IOProcessor()) {
        const auto atm = read_lines(dir / "t.atm");
        ASSERT_EQ(atm.size(), 5u);
        EXPECT_EQ(atm[2], "12 30 1100 t_000001_vel.asc t_000001_ang.asc");
        EXPECT_EQ(atm[3], "12 30 1101 t_000004_vel.asc t_000004_ang.asc");
        EXPECT_EQ(atm[4], "12 30 1102 t_000006_vel.asc t_000006_ang.asc");
    }

    // A fresh run (t = 0) starts the .atm over, and a prj_file is copied
    // beside each grid
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream(dir / "in.prj") << "PROJCS[\"test\"]";
    }
    ParallelDescriptor::Barrier();
    p.prj_file = (dir / "in.prj").string();
    {
        FarsiteWindWriter w(p, geom, 10.0);
        w.maybe_write(wind, 0.0, 1);
    }
    ParallelDescriptor::Barrier();
    if (ParallelDescriptor::IOProcessor()) {
        EXPECT_EQ(read_lines(dir / "t.atm").size(), 3u);
        const auto prj = read_lines(dir / "t_000001_ang.prj");
        ASSERT_EQ(prj.size(), 1u);
        EXPECT_EQ(prj[0], "PROJCS[\"test\"]");
        std::filesystem::remove_all(dir);
    }
}

// coarsen = 2 on fire boxes that coarsen evenly (4 x 4: averaged on the
// owning ranks before the gather) and on boxes that do not (3 + 3 + 2:
// averaged on the I/O rank after it) writes the same grids
TEST(FarsiteWinds, CoarsenOnDeviceEqualsOnHost)
{
    ScopedStart start("2021-12-30 11:00");
    const auto dir = erf_gtest_temp_path("farsite_coarsen");
    Box dom(IntVect(0, 0, 0), IntVect(7, 3, 0));
    RealBox rb({0.0, 0.0, 0.0}, {80.0, 40.0, 1.0});
    Array<int, AMREX_SPACEDIM> per{0, 0, 0};
    Geometry geom(dom, rb, CoordSys::cartesian, per);
    std::vector<std::string> text;
    for (int maxsz : {4, 3}) {
        BoxArray ba(dom);
        ba.maxSize(IntVect(maxsz, maxsz, 1));
        EXPECT_EQ(ba.coarsenable(IntVect(2, 2, 1)), maxsz == 4);
        DistributionMapping dm(ba);
        MultiFab wind(ba, dm, 2, 0);
        fill_test_wind(wind);
        FireParams::FarsiteWindsParams p;
        p.interval = 60.0;
        p.units    = "english";
        p.dir      = (dir / std::to_string(maxsz)).string();
        p.coarsen  = 2;
        FarsiteWindWriter w(p, geom, 6.096);
        w.maybe_write(wind, Real(0.0), 1);
        ParallelDescriptor::Barrier();
        if (ParallelDescriptor::IOProcessor()) {
            std::string all;
            for (const char* k : {"_vel.asc", "_ang.asc"}) {
                std::ifstream f(dir / std::to_string(maxsz) / (std::string("flare_000001") + k));
                std::stringstream ss; ss << f.rdbuf(); all += ss.str();
            }
            text.push_back(all);
        }
    }
    if (ParallelDescriptor::IOProcessor()) {
        ASSERT_EQ(text.size(), 2u);
        EXPECT_FALSE(text[0].empty());
        EXPECT_EQ(text[0], text[1]);
        std::filesystem::remove_all(dir);
    }
}
