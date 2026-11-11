// Unit tests of erf_conductors::IntensityRecorder (ERF_Intensity.H), the averaged wind kept beside each
// component's peak demand. Every expected number is worked out by hand in the test's comments.
//
// TheAverageIsOverTheLastTSeconds: equal steps, before and after the window is full.
// StepsOfDifferentLengthsCountWithTheirPartInTheWindow: time weights, the oldest sample cut at the window's start.
// AStepLongerThanTheWindowAndAZeroWindowKeepTheLastSample.
// TheSpeedIsAveragedAsASpeedAndTheDirectionIsTheMeanVectors: a turning wind keeps its speed; the direction's
//     convention (degrees from +x towards +y, the way the wind blows) and the zero vector.
// TheMaximumWaitsForAFullWindow: averages over a window not yet full are not the run's largest T-second average.
// TheRunsMaximumKeepsItsTimeDirectionAndK: the largest averaged speed, not the last, with k averaged alike.
// AComponentKeepsTheWindAtItsPeakNotTheLastWind: at_peak/at_first take the wind of that step; a component
//     without a point keeps none.
// TheTablesCarryEveryField: demand_table.csv and intensity.csv, with empty fields where nothing exists.
// ARestartContinuesTheWindowAndTheMaxima: state() and restore() in the middle of a run give the straight run's
//     state to the bit; a changed height, averaging time or closure is refused by name; a record without a point
//     is dropped, a point without a record starts afresh; a malformed record is refused.
// BadInputsAreRefused: the ranges of the two inputs, names, sample counts, steps, times and values.

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "ERF_Capacities.H"
#include "ERF_Intensity.H"
#include "ERF_GTestTempDir.H"
#include "ERF_GTestThrowOnAbort.H"

using erf_conductors::CapacityMonitor;
using erf_conductors::ComponentKind;
using erf_conductors::IntensityRecorder;
using erf_conductors::IntensitySample;

namespace {
IntensityRecorder one_point (double T, bool with_k = false)
{
    IntensityRecorder r;
    EXPECT_EQ(r.configure(10.0, T, with_k), "");
    r.add_point("P", 100.0, 200.0, 35.0);
    return r;
}
IntensitySample wind (double u, double v = 0.0, double k = 0.0) { IntensitySample s; s.u = u; s.v = v; s.k = k; return s; }
std::vector<std::string> lines_of (const std::string& f)
{
    std::ifstream in(f);
    std::vector<std::string> rows;
    for (std::string l; std::getline(in, l);) { rows.push_back(l); }
    return rows;
}
} // namespace

TEST(Intensity, TheAverageIsOverTheLastTSeconds)
{
    // T = 3 s, steps of 0.5 s: six samples fill the window. Speeds 1, 2, ..., 8 m/s.
    IntensityRecorder r = one_point(3.0);
    double t = 0.0;
    for (int n = 1; n <= 8; ++n) {
        t += 0.5;
        r.sample({wind(static_cast<double>(n))}, 0.5, t);
        if (n == 1) { EXPECT_DOUBLE_EQ(r.point(0).now.speed, 1.0); }
        if (n == 2) { EXPECT_DOUBLE_EQ(r.point(0).now.speed, 1.5) << "before the run is T old: what has been sampled"; }
        if (n == 6) { EXPECT_DOUBLE_EQ(r.point(0).now.speed, 3.5) << "(1 + ... + 6) / 6"; }
        if (n == 7) { EXPECT_DOUBLE_EQ(r.point(0).now.speed, 4.5) << "(2 + ... + 7) / 6: the first sample has left"; }
    }
    EXPECT_DOUBLE_EQ(r.point(0).now.speed, 5.5) << "(3 + ... + 8) / 6";
    EXPECT_DOUBLE_EQ(r.point(0).now.time, 4.0);
    EXPECT_EQ(r.num_samples(), 8u);
}

TEST(Intensity, StepsOfDifferentLengthsCountWithTheirPartInTheWindow)
{
    // T = 3 s. Steps (dt, speed): (2, 10), (1, 4), (2, 1).
    IntensityRecorder r = one_point(3.0);
    r.sample({wind(10.0)}, 2.0, 2.0);
    EXPECT_DOUBLE_EQ(r.point(0).now.speed, 10.0);
    r.sample({wind(4.0)}, 1.0, 3.0);
    // the last 3 s: 1 s at 4 and 2 s at 10
    EXPECT_DOUBLE_EQ(r.point(0).now.speed, (1.0 * 4.0 + 2.0 * 10.0) / 3.0);
    r.sample({wind(1.0)}, 2.0, 5.0);
    // the last 3 s: 2 s at 1 and 1 s at 4; the 10 has left the window
    EXPECT_DOUBLE_EQ(r.point(0).now.speed, (2.0 * 1.0 + 1.0 * 4.0) / 3.0);

    // the oldest sample cut at the window's start: (2, 10), (2, 1): 2 s at 1 and the last 1 s of the 10
    IntensityRecorder q = one_point(3.0);
    q.sample({wind(10.0)}, 2.0, 2.0);
    q.sample({wind(1.0)}, 2.0, 4.0);
    EXPECT_DOUBLE_EQ(q.point(0).now.speed, (2.0 * 1.0 + 1.0 * 10.0) / 3.0);
}

TEST(Intensity, AStepLongerThanTheWindowAndAZeroWindowKeepTheLastSample)
{
    IntensityRecorder r = one_point(3.0);
    r.sample({wind(7.0)}, 5.0, 5.0);
    r.sample({wind(2.0)}, 5.0, 10.0);
    EXPECT_DOUBLE_EQ(r.point(0).now.speed, 2.0) << "the step covers the whole window";
    IntensityRecorder z = one_point(0.0);
    z.sample({wind(7.0)}, 0.5, 0.5);
    z.sample({wind(2.0)}, 0.5, 1.0);
    EXPECT_DOUBLE_EQ(z.point(0).now.speed, 2.0) << "averaging time 0: the last sample";
    EXPECT_DOUBLE_EQ(z.point(0).max.speed, 7.0);
}

TEST(Intensity, TheSpeedIsAveragedAsASpeedAndTheDirectionIsTheMeanVectors)
{
    // a wind of 3 m/s turning from +x to +y: the mean speed stays 3 (the mean vector's length is 3 / sqrt 2)
    IntensityRecorder r = one_point(1.0);
    r.sample({wind(3.0, 0.0)}, 0.5, 0.5);
    EXPECT_DOUBLE_EQ(r.point(0).now.direction, 0.0);
    r.sample({wind(0.0, 3.0)}, 0.5, 1.0);
    EXPECT_DOUBLE_EQ(r.point(0).now.speed, 3.0);
    EXPECT_NEAR(r.point(0).now.direction, 45.0, 1e-12);
    // opposite winds: the speed is still 2, the mean vector is zero and its direction reads 0
    IntensityRecorder o = one_point(1.0);
    o.sample({wind(2.0, 0.0)}, 0.5, 0.5);
    o.sample({wind(-2.0, 0.0)}, 0.5, 1.0);
    EXPECT_DOUBLE_EQ(o.point(0).now.speed, 2.0);
    EXPECT_DOUBLE_EQ(o.point(0).now.direction, 0.0);
    // towards -x is 180, towards -y is -90
    IntensityRecorder w = one_point(0.0);
    w.sample({wind(-1.0, 0.0)}, 0.5, 0.5);
    EXPECT_NEAR(w.point(0).now.direction, 180.0, 1e-12);
    w.sample({wind(0.0, -1.0)}, 0.5, 1.0);
    EXPECT_NEAR(w.point(0).now.direction, -90.0, 1e-12);
}

TEST(Intensity, TheMaximumWaitsForAFullWindow)
{
    // T = 3 s, steps of 1 s: speeds 20, 2, 2, 5, 5, 5. The averages are 20, 11, 8, 3, 4, 5: the first two are over
    // 1 s and 2 s, shorter averages, and are not the run's largest 3-second average, which is 8 at t = 3 s
    IntensityRecorder r = one_point(3.0);
    const double s[6] = {20.0, 2.0, 2.0, 5.0, 5.0, 5.0};
    for (int n = 0; n < 6; ++n) {
        r.sample({wind(s[n])}, 1.0, 1.0 * (n + 1));
        if (n < 2) { EXPECT_FALSE(r.point(0).max.set) << "the window is not full after " << n + 1 << " s"; }
        if (n == 1) { EXPECT_DOUBLE_EQ(r.point(0).now.speed, 11.0) << "the value now is still the average of what was sampled"; }
    }
    EXPECT_DOUBLE_EQ(r.point(0).max.speed, 8.0);
    EXPECT_DOUBLE_EQ(r.point(0).max.time, 3.0);
    // ten steps of 0.3 s sum to just under 3 s in floating point: the window counts as full at the tenth
    IntensityRecorder q = one_point(3.0);
    for (int n = 1; n <= 10; ++n) {
        q.sample({wind(1.0)}, 0.3, 0.3 * n);
        EXPECT_EQ(q.point(0).max.set, n == 10) << n;
    }
}

TEST(Intensity, TheRunsMaximumKeepsItsTimeDirectionAndK)
{
    // T = 1 s, steps of 0.5 s: averages of pairs. Speeds 2, 6, 8, 4, 1 along +y, k = speed / 10.
    IntensityRecorder r = one_point(1.0, true);
    const double s[5] = {2.0, 6.0, 8.0, 4.0, 1.0};
    for (int n = 0; n < 5; ++n) { r.sample({wind(0.0, s[n], s[n] / 10.0)}, 0.5, 0.5 * (n + 1)); }
    // pair means: (2 over a half-full window), 4, 7, 6, 2.5: the maximum 7 at t = 1.5 s
    EXPECT_DOUBLE_EQ(r.point(0).now.speed, 2.5);
    EXPECT_DOUBLE_EQ(r.point(0).max.speed, 7.0);
    EXPECT_DOUBLE_EQ(r.point(0).max.time, 1.5);
    EXPECT_NEAR(r.point(0).max.direction, 90.0, 1e-12);
    EXPECT_DOUBLE_EQ(r.point(0).max.k, 0.7) << "k averaged over the same window: (0.6 + 0.8) / 2";
    EXPECT_DOUBLE_EQ(r.point(0).now.k, 0.25);
    // without a closure's k the value is 0 and not averaged
    IntensityRecorder n = one_point(1.0, false);
    n.sample({wind(1.0, 0.0, 5.0)}, 0.5, 0.5);
    EXPECT_DOUBLE_EQ(n.point(0).now.k, 0.0);
}

TEST(Intensity, AComponentKeepsTheWindAtItsPeakNotTheLastWind)
{
    CapacityMonitor m;
    const std::size_t a = m.add(ComponentKind::Conductor, "S1");
    const std::size_t b = m.add(ComponentKind::Member, "T1_m1");
    IntensityRecorder r;
    ASSERT_EQ(r.configure(10.0, 0.0, false), "");
    const std::size_t p = r.add_point("S1", 0.0, 0.0, 10.0);
    r.components(m.size());
    r.link(a, p);
    EXPECT_EQ(r.point_of(a), p);
    EXPECT_EQ(r.point_of(b), r.num_points()) << "no point";
    // the demand peaks at step 2 (ratio 1.2, the first exceedance) while the wind goes on rising
    const double demand[4] = {50.0, 120.0, 90.0, 110.0}, speed[4] = {5.0, 9.0, 11.0, 14.0};
    for (int n = 0; n < 4; ++n) {
        const double t = 0.5 * (n + 1);
        r.sample({wind(speed[n])}, 0.5, t);
        const double before = m.component(a).peak_ratio;
        const bool seen = m.component(a).seen;
        const bool first = m.record(a, demand[n], 100.0, t);
        if (!seen || m.component(a).peak_ratio > before) { r.at_peak(a); }
        if (first) { r.at_first(a); }
        r.at_peak(b);   // a component without a point: nothing kept
    }
    EXPECT_TRUE(r.peak_value(a).set);
    EXPECT_DOUBLE_EQ(r.peak_value(a).speed, 9.0) << "the wind at the peak's step, not the last (14) or the largest";
    EXPECT_DOUBLE_EQ(r.peak_value(a).time, 1.0);
    EXPECT_DOUBLE_EQ(r.first_value(a).speed, 9.0);
    EXPECT_DOUBLE_EQ(r.point(p).max.speed, 14.0);
    EXPECT_FALSE(r.peak_value(b).set);
    EXPECT_FALSE(r.first_value(b).set);
}

TEST(Intensity, TheTablesCarryEveryField)
{
    CapacityMonitor m;
    const std::size_t a = m.add(ComponentKind::Conductor, "S1");
    const std::size_t b = m.add(ComponentKind::Insulator, "S1_string1");
    const std::size_t c = m.add(ComponentKind::Member, "T1_m3");
    IntensityRecorder r;
    ASSERT_EQ(r.configure(10.0, 0.0, true), "");
    const std::size_t p = r.add_point("S1", 150.0, 20.0, 40.0);
    const std::size_t q = r.add_point("tower_T1", 300.0, 20.0, 42.0);
    r.components(m.size());
    r.link(a, p);
    r.link(b, q);
    // c has no point; b is never recorded
    r.sample({wind(0.0, 8.0, 0.5), wind(6.0, 0.0, 0.25)}, 0.5, 2.5);
    EXPECT_TRUE(m.record(a, 150.0, 100.0, 2.5));
    r.at_peak(a);
    r.at_first(a);
    m.record(c, 20.0, 100.0, 2.5);
    r.at_peak(c);
    r.sample({wind(0.0, 4.0, 0.125), wind(3.0, 0.0, 1.0)}, 0.5, 3.0);
    const auto dir = erf_gtest_temp_path("erf_gtest_intensity");
    std::filesystem::create_directories(dir);
    r.write_demand_table((dir / "demand_table.csv").string(), m);
    r.write_points((dir / "intensity.csv").string());
    const auto t = lines_of((dir / "demand_table.csv").string());
    const auto pts = lines_of((dir / "intensity.csv").string());
    std::filesystem::remove_all(dir);
    ASSERT_EQ(t.size(), 4u);
    EXPECT_EQ(t[0], "kind,component,peak_ratio,peak_demand,peak_capacity,peak_time,point,"
                    "peak_speed,peak_direction,peak_k,peak_wind_time,max_speed,max_direction,max_k,max_time,"
                    "exceeded,first_time,first_demand,first_capacity,first_speed,first_direction,first_k,first_wind_time");
    EXPECT_EQ(t[1], "conductor,S1,1.5,150,100,2.5,S1,8,90,0.5,2.5,8,90,0.5,2.5,1,2.5,150,100,8,90,0.5,2.5");
    EXPECT_EQ(t[2], "insulator,S1_string1,,,,,tower_T1,,,,,6,0,0.25,2.5,0,,,,,,,");
    EXPECT_EQ(t[3], "member,T1_m3,0.2,20,100,2.5,,,,,,,,,,0,,,,,,,");
    ASSERT_EQ(pts.size(), 3u);
    EXPECT_EQ(pts[0], "point,x,y,z,height,speed,direction,k,time,max_speed,max_direction,max_k,max_time");
    EXPECT_EQ(pts[1], "S1,150,20,40,10,4,90,0.125,3,8,90,0.5,2.5");
    EXPECT_EQ(pts[2], "tower_T1,300,20,42,10,3,0,1,3,6,0,0.25,2.5");
    // every row has the header's number of fields
    for (const auto& row : t) { EXPECT_EQ(std::count(row.begin(), row.end(), ','), 22); }
    for (const auto& row : pts) { EXPECT_EQ(std::count(row.begin(), row.end(), ','), 12); }
}

TEST(Intensity, ARestartContinuesTheWindowAndTheMaxima)
{
    auto make = [] (CapacityMonitor& m, IntensityRecorder& r, double height = 10.0, double T = 3.0, bool with_k = true) {
        m = CapacityMonitor();
        const std::size_t a = m.add(ComponentKind::Conductor, "S1");
        EXPECT_EQ(r.configure(height, T, with_k), "");
        const std::size_t p = r.add_point("S1", 1.0 / 3.0, 2.0, 30.0);
        r.add_point("S1_end_a", 0.0, 2.0, 31.0);
        r.components(m.size());
        r.link(a, p);
    };
    auto step = [] (CapacityMonitor& m, IntensityRecorder& r, int n) {
        // steps of unequal length and winds without short decimal forms
        const double dt = 0.4 + 0.1 * (n % 3), t = 0.1 + 0.55 * n + 0.1 * (n % 3);
        const double u = 5.0 + 3.0 * std::sin(0.9 * n), v = 2.0 * std::cos(1.3 * n);
        r.sample({wind(u, v, 0.1 * (n + 1) / 3.0), wind(0.5 * u, -v, 0.2)}, dt, t);
        const double before = m.component(0).peak_ratio;
        const bool seen = m.component(0).seen;
        const bool first = m.record(0, 80.0 + 30.0 * std::sin(0.7 * n), 100.0, t);
        if (!seen || m.component(0).peak_ratio > before) { r.at_peak(0); }
        if (first) { r.at_first(0); }
    };
    CapacityMonitor ms, mr;
    IntensityRecorder straight, first_leg, second_leg;
    make(ms, straight);
    for (int n = 1; n <= 14; ++n) { step(ms, straight, n); }
    make(mr, first_leg);
    for (int n = 1; n <= 8; ++n) { step(mr, first_leg, n); }
    const std::string chk = first_leg.state(mr), cap = mr.state();
    CapacityMonitor m2;
    make(m2, second_leg);
    std::vector<std::string> dropped;
    ASSERT_EQ(m2.restore(cap, dropped), "");
    ASSERT_EQ(second_leg.restore(chk, m2, dropped), "");
    EXPECT_TRUE(dropped.empty());
    EXPECT_EQ(second_leg.state(m2), chk) << "what was written comes back";
    for (int n = 9; n <= 14; ++n) { step(m2, second_leg, n); }
    EXPECT_EQ(second_leg.state(m2), straight.state(ms)) << "the restarted run is the straight run, to the bit";
    EXPECT_EQ(second_leg.num_samples(), 14u);
    EXPECT_GT(straight.point(0).max.speed, straight.point(0).now.speed) << "the case has a maximum before its end";

    // a changed input is refused by its name
    CapacityMonitor mc;
    IntensityRecorder other;
    make(mc, other, 20.0);
    EXPECT_NE(other.restore(chk, mc, dropped).find("conductors.intensity_height = 10, the run has 20"), std::string::npos);
    make(mc, other, 10.0, 1.0);
    EXPECT_NE(other.restore(chk, mc, dropped).find("conductors.intensity_averaging_time = 3, the run has 1"), std::string::npos);
    make(mc, other, 10.0, 3.0, false);
    EXPECT_NE(other.restore(chk, mc, dropped).find("with a turbulent kinetic energy and the run is without one (erf.rans_type"), std::string::npos);

    // a record without a point or component is dropped and named; a point without a record starts afresh
    CapacityMonitor mf;
    mf.add(ComponentKind::Conductor, "S9");
    IntensityRecorder fewer;
    ASSERT_EQ(fewer.configure(10.0, 3.0, true), "");
    fewer.add_point("S1_end_a", 0.0, 2.0, 31.0);
    const std::size_t fresh = fewer.add_point("S9", 5.0, 2.0, 30.0);
    fewer.components(mf.size());
    fewer.link(0, fresh);
    ASSERT_EQ(fewer.restore(chk, mf, dropped), "");
    ASSERT_EQ(dropped.size(), 2u);
    EXPECT_EQ(dropped[0], "S1");
    EXPECT_EQ(dropped[1], "S1") << "the point S1 and the component S1";
    EXPECT_TRUE(fewer.point(0).max.set);
    EXPECT_FALSE(fewer.point(fresh).now.set);
    fewer.sample({wind(1.0, 0.0, 0.1), wind(4.0, 0.0, 0.1)}, 0.5, 100.0);
    EXPECT_DOUBLE_EQ(fewer.point(fresh).now.speed, 4.0) << "its window holds its own first sample alone";
    // the point that joined holds fewer samples than the window's steps: the next checkpoint must still come back
    const std::string again = fewer.state(mf);
    CapacityMonitor mg;
    mg.add(ComponentKind::Conductor, "S9");
    IntensityRecorder back;
    ASSERT_EQ(back.configure(10.0, 3.0, true), "");
    back.add_point("S1_end_a", 0.0, 2.0, 31.0);
    back.add_point("S9", 5.0, 2.0, 30.0);
    back.components(mg.size());
    back.link(0, 1);
    ASSERT_EQ(back.restore(again, mg, dropped), "");
    EXPECT_EQ(back.state(mg), again);

    // malformed records
    make(mc, other);
    EXPECT_NE(other.restore("point S1 1 2 3", mc, dropped).find("malformed intensity record on line 1"), std::string::npos);
    EXPECT_NE(other.restore("", mc, dropped).find("holds no intensity record"), std::string::npos);
    EXPECT_NE(other.restore(chk + "what 1\n", mc, dropped).find("malformed intensity record"), std::string::npos);
    // the steps come once, before the points; counts are not negative; a point holds no more samples than steps
    make(mc, other);
    EXPECT_NE(other.restore(chk + "steps 1 0.5\n", mc, dropped).find("malformed intensity record"), std::string::npos);
    make(mc, other);
    EXPECT_NE(other.restore("intensity 10 3 1 -4 1 2\n", mc, dropped).find("malformed intensity record on line 1"), std::string::npos);
    make(mc, other);
    EXPECT_NE(other.restore("intensity 10 3 1 4 1 2\nsteps -1\n", mc, dropped).find("line 2"), std::string::npos);
    make(mc, other);
    EXPECT_NE(other.restore("intensity 10 3 1 4 1 2\nsteps 1 0.5\npoint S1 0 0 0 0 0 0 0 0 0 0 2 1 0 0 1 0 0\n", mc, dropped).find("line 3"),
              std::string::npos);
    make(mc, other);
    EXPECT_NE(other.restore("intensity 10 3 1 4 1 2\npoint S1 0 0 0 0 0 0 0 0 0 0 0\n", mc, dropped).find("line 2"), std::string::npos)
        << "a point before the steps";
}

TEST(Intensity, BadInputsAreRefused)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(IntensityRecorder::check(10.0, 3.0), "");
    EXPECT_EQ(IntensityRecorder::check(500.0, 600.0), "");
    EXPECT_EQ(IntensityRecorder::check(0.5, 0.0), "");
    for (const double h : {0.0, -1.0, 500.5, nan}) {
        EXPECT_NE(IntensityRecorder::check(h, 3.0).find("conductors.intensity_height must be in (0, 500]"), std::string::npos) << h;
    }
    for (const double T : {-0.1, 600.5, nan}) {
        EXPECT_NE(IntensityRecorder::check(10.0, T).find("conductors.intensity_averaging_time must be in [0, 600]"), std::string::npos) << T;
    }
    IntensityRecorder r;
    EXPECT_NE(r.configure(0.0, 3.0, false), "");
    ASSERT_EQ(r.configure(10.0, 3.0, true), "");
    EXPECT_FALSE(erf_gtest::abort_message([&] { r.add_point("", 0.0, 0.0, 0.0); }).empty());
    EXPECT_FALSE(erf_gtest::abort_message([&] { r.add_point("a b", 0.0, 0.0, 0.0); }).empty());
    EXPECT_FALSE(erf_gtest::abort_message([&] { r.add_point("a,b", 0.0, 0.0, 0.0); }).empty());
    EXPECT_FALSE(erf_gtest::abort_message([&] { r.add_point("P", nan, 0.0, 0.0); }).empty());
    r.add_point("P", 0.0, 0.0, 10.0);
    EXPECT_NE(erf_gtest::abort_message([&] { r.add_point("P", 1.0, 0.0, 10.0); }).find("two intensity points are named 'P'"), std::string::npos);
    EXPECT_FALSE(erf_gtest::abort_message([&] { r.link(0, 0); }).empty()) << "components() first";
    r.components(1);
    EXPECT_FALSE(erf_gtest::abort_message([&] { r.link(0, 3); }).empty());
    r.link(0, 0);
    EXPECT_FALSE(erf_gtest::abort_message([&] { r.link(0, 0); }).empty()) << "linked once";
    EXPECT_FALSE(erf_gtest::abort_message([&] { r.sample({}, 0.5, 0.5); }).empty()) << "one sample per point";
    EXPECT_FALSE(erf_gtest::abort_message([&] { r.sample({wind(1.0)}, 0.0, 0.5); }).empty());
    EXPECT_FALSE(erf_gtest::abort_message([&] { r.sample({wind(1.0)}, 0.5, nan); }).empty());
    EXPECT_NE(erf_gtest::abort_message([&] { r.sample({wind(nan)}, 0.5, 0.5); }).find("the wind at the intensity point P is not finite"), std::string::npos);
    EXPECT_NE(erf_gtest::abort_message([&] { r.sample({wind(1.0, 0.0, -1.0)}, 0.5, 0.5); }).find("turbulent kinetic energy at the intensity point P"), std::string::npos);
    r.sample({wind(1.0, 0.0, 0.1)}, 0.5, 0.5);
    EXPECT_FALSE(erf_gtest::abort_message([&] { r.sample({wind(1.0, 0.0, 0.1)}, 0.5, 0.5); }).empty()) << "the time must advance";
    EXPECT_FALSE(erf_gtest::abort_message([&] { r.add_point("Q", 0.0, 0.0, 10.0); }).empty()) << "no points after the first sample";
    CapacityMonitor two;
    two.add(ComponentKind::Conductor, "A");
    two.add(ComponentKind::Conductor, "B");
    EXPECT_FALSE(erf_gtest::abort_message([&] { (void)r.state(two); }).empty()) << "another monitor";
}
