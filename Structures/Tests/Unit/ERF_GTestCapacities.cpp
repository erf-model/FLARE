// Unit tests of erf_conductors::CapacityMonitor (ERF_Capacities.H), the demand-to-capacity record of the
// conductor lines' and towers' components.
//
// ADemandEqualToItsCapacityIsNotAnExceedance: the test is demand > capacity, and the first exceedance is
//     reported once, with its time, demand and capacity, later ones not again.
// ThePeakIsTheLargestRatioNotTheLargestDemand: with a capacity that changes between records (a member's strength at
//     its temperature), the peak follows demand / capacity.
// TheRunsFirstFailureIsTheEarliestExceedance: first_failure() over components exceeding at different times, the
//     first added on a tie, size() when none has; order_by_first_time sorts a step's exceedances, the first added first on a tie.
// TheStateComesBackByNameAndKind: state() and restore() round trip to the bit; a component without a record starts
//     afresh, a record without a component (or of another kind) is dropped and named, a malformed one refused.
// TheTableAndTheFailureRowCarryEveryField: capacities.csv's header and rows, empty fields before a record and
//     without an exceedance; the failure log's row.
// AMembersDemandIsOnTheSideThatGoverns: member_demand picks tension or compression by their ratios.
// BadNamesDemandsAndCapacitiesAreRefused: an empty or spaced name, a repeated one, a negative or NaN demand, a
//     capacity that is not positive.

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "ERF_Capacities.H"
#include "ERF_GTestTempDir.H"
#include "ERF_GTestThrowOnAbort.H"

using erf_conductors::CapacityMonitor;
using erf_conductors::ComponentKind;

TEST(Capacities, ADemandEqualToItsCapacityIsNotAnExceedance)
{
    CapacityMonitor m;
    const std::size_t c = m.add(ComponentKind::Conductor, "S1_span1");
    EXPECT_FALSE(m.record(c, 1.0e5, 1.4e5, 0.5));
    EXPECT_FALSE(m.record(c, 1.4e5, 1.4e5, 1.0)) << "demand == capacity is at the capacity, not over it";
    EXPECT_FALSE(m.component(c).exceeded);
    EXPECT_TRUE(m.record(c, 1.4e5 * (1.0 + 1.0e-12), 1.4e5, 1.5)) << "the first record over the capacity";
    EXPECT_FALSE(m.record(c, 2.0e5, 1.4e5, 2.0)) << "reported once";
    const auto& r = m.component(c);
    EXPECT_TRUE(r.exceeded);
    EXPECT_EQ(r.first_time, 1.5);
    EXPECT_EQ(r.first_capacity, 1.4e5);
    EXPECT_EQ(r.first_demand, 1.4e5 * (1.0 + 1.0e-12));
    EXPECT_EQ(r.peak_demand, 2.0e5);
    EXPECT_EQ(r.peak_time, 2.0);
    EXPECT_EQ(r.demand, 2.0e5);
}

TEST(Capacities, ThePeakIsTheLargestRatioNotTheLargestDemand)
{
    CapacityMonitor m;
    const std::size_t c = m.add(ComponentKind::Member, "T1_m7");
    m.record(c, 5.0e4, 1.0e5, 1.0);    // 0.5
    m.record(c, 6.0e4, 2.0e5, 2.0);    // 0.3: more demand, less ratio
    m.record(c, 4.0e4, 5.0e4, 3.0);    // 0.8: the steel weakened by heat
    m.record(c, 7.0e4, 1.0e5, 4.0);    // 0.7
    const auto& r = m.component(c);
    EXPECT_DOUBLE_EQ(r.peak_ratio, 0.8);
    EXPECT_EQ(r.peak_demand, 4.0e4);
    EXPECT_EQ(r.peak_capacity, 5.0e4);
    EXPECT_EQ(r.peak_time, 3.0);
    EXPECT_FALSE(r.exceeded);
    // the first record sets the peak even at a zero demand
    const std::size_t z = m.add(ComponentKind::FootingUplift, "T1_footing_uplift");
    m.record(z, 0.0, 1.0e5, 9.0);
    EXPECT_TRUE(m.component(z).seen);
    EXPECT_EQ(m.component(z).peak_time, 9.0);
}

TEST(Capacities, TheRunsFirstFailureIsTheEarliestExceedance)
{
    CapacityMonitor m;
    const std::size_t a = m.add(ComponentKind::Conductor, "a");
    const std::size_t b = m.add(ComponentKind::Hardware, "b");
    const std::size_t c = m.add(ComponentKind::Insulator, "c");
    EXPECT_EQ(m.first_failure(), m.size());
    m.record(a, 2.0, 1.0, 3.0);
    m.record(b, 2.0, 1.0, 2.0);
    EXPECT_EQ(m.first_failure(), b);
    m.record(c, 2.0, 1.0, 2.0);   // the same time as b: the first added wins
    EXPECT_EQ(m.first_failure(), b);
    // a step's exceedances in the order the groups stepped, logged in time order; at one time the component added
    // first (capacities.csv's order), whichever group stepped first
    std::vector<std::size_t> ids{a, c, b};
    m.order_by_first_time(ids);
    EXPECT_EQ(ids, (std::vector<std::size_t>{b, c, a}));
}

TEST(Capacities, TheStateComesBackByNameAndKind)
{
    CapacityMonitor m;
    const std::size_t a = m.add(ComponentKind::Conductor, "S1_span1");
    const std::size_t b = m.add(ComponentKind::FootingCompression, "T2_footing_compression");
    m.add(ComponentKind::Hardware, "S1_end_a");   // never recorded
    m.record(a, 1.0 / 3.0, 0.1, 0.1 + 0.2);
    m.record(a, 0.05, 0.1, 0.7);
    m.record(b, 2.0e5 / 7.0, 3.0e5, 1.0 / 7.0);
    const std::string text = m.state();

    // the restarted run: S1_end_a is gone, S1_end_b is new, T2's footing has changed its kind
    CapacityMonitor n;
    const std::size_t na = n.add(ComponentKind::Conductor, "S1_span1");
    const std::size_t nb = n.add(ComponentKind::FootingUplift, "T2_footing_compression");
    const std::size_t nc = n.add(ComponentKind::Hardware, "S1_end_b");
    std::vector<std::string> dropped;
    EXPECT_TRUE(n.restore(text, dropped).empty());
    ASSERT_EQ(dropped.size(), 2u);
    EXPECT_EQ(dropped[0], "T2_footing_compression");
    EXPECT_EQ(dropped[1], "S1_end_a");
    const auto& r = n.component(na);
    const auto& o = m.component(a);
    EXPECT_EQ(r.demand, o.demand);
    EXPECT_EQ(r.peak_ratio, o.peak_ratio);
    EXPECT_EQ(r.peak_time, o.peak_time) << "17 digits bring every double back to the bit";
    EXPECT_EQ(r.first_time, o.first_time);
    EXPECT_EQ(r.first_demand, o.first_demand);
    EXPECT_TRUE(r.exceeded);
    EXPECT_FALSE(n.component(nb).seen);
    EXPECT_FALSE(n.component(nc).seen);
    // a restored exceedance is not reported again
    EXPECT_FALSE(n.record(na, 1.0, 0.1, 1.0));

    EXPECT_FALSE(n.restore("component S1_span1 0 1 nan\n", dropped).empty());
    EXPECT_FALSE(n.restore("components 3\n", dropped).empty());
    EXPECT_FALSE(n.restore("component S1_span1 9 1 1 1 1 1 1 1 0 0 0 0\n", dropped).empty()) << "no such kind";
}

TEST(Capacities, TheTableAndTheFailureRowCarryEveryField)
{
    CapacityMonitor m;
    const std::size_t a = m.add(ComponentKind::Conductor, "S1");
    m.add(ComponentKind::Insulator, "S2_string1");
    const std::size_t c = m.add(ComponentKind::Member, "T1_m3");
    m.record(a, 150.0, 100.0, 2.5);
    m.record(c, 20.0, 100.0, 2.5);
    const auto dir = erf_gtest_temp_path("erf_gtest_capacities");
    std::filesystem::create_directories(dir);
    const std::string f = (dir / "capacities.csv").string();
    m.write_table(f);
    std::ifstream in(f);
    std::vector<std::string> rows;
    for (std::string l; std::getline(in, l);) { rows.push_back(l); }
    std::filesystem::remove_all(dir);
    ASSERT_EQ(rows.size(), 4u);
    EXPECT_EQ(rows[0], "kind,component,demand,capacity,peak_ratio,peak_demand,peak_capacity,peak_time,"
                       "exceeded,first_time,first_demand,first_capacity");
    EXPECT_EQ(rows[1], "conductor,S1,150,100,1.5,150,100,2.5,1,2.5,150,100");
    EXPECT_EQ(rows[2], "insulator,S2_string1,,,,,,,0,,,");
    EXPECT_EQ(rows[3], "member,T1_m3,20,100,0.2,20,100,2.5,0,,,");
    EXPECT_EQ(m.failure_row(a), "2.5 conductor S1 150 100 1.5");
    EXPECT_STREQ(CapacityMonitor::failure_log_header(), "time kind component demand capacity ratio");
    EXPECT_STREQ(erf_conductors::kind_name(ComponentKind::FootingUplift), "footing_uplift");
    EXPECT_STREQ(erf_conductors::kind_name(ComponentKind::FootingCompression), "footing_compression");
    EXPECT_STREQ(erf_conductors::kind_name(ComponentKind::Hardware), "hardware");
}

TEST(Capacities, AMembersDemandIsOnTheSideThatGoverns)
{
    double d = 0.0, c = 0.0;
    erf_conductors::member_demand(4.0e4, 3.0e4, 1.0e5, 5.0e4, d, c);   // 0.4 in tension, 0.6 in compression
    EXPECT_EQ(d, 3.0e4);
    EXPECT_EQ(c, 5.0e4);
    erf_conductors::member_demand(6.0e4, 1.0e4, 1.0e5, 5.0e4, d, c);   // 0.6 against 0.2
    EXPECT_EQ(d, 6.0e4);
    EXPECT_EQ(c, 1.0e5);
    erf_conductors::member_demand(0.0, 0.0, 1.0e5, 5.0e4, d, c);       // unloaded: tension, at 0
    EXPECT_EQ(d, 0.0);
    EXPECT_EQ(c, 1.0e5);
    EXPECT_FALSE(erf_gtest::abort_message([&] { erf_conductors::member_demand(1.0, 1.0, 0.0, 1.0, d, c); }).empty());
}

TEST(Capacities, BadNamesDemandsAndCapacitiesAreRefused)
{
    CapacityMonitor m;
    EXPECT_FALSE(erf_gtest::abort_message([&] { m.add(ComponentKind::Conductor, ""); }).empty());
    EXPECT_FALSE(erf_gtest::abort_message([&] { m.add(ComponentKind::Conductor, "a b"); }).empty());
    EXPECT_FALSE(erf_gtest::abort_message([&] { m.add(ComponentKind::Conductor, "a,b"); }).empty());
    const std::size_t c = m.add(ComponentKind::Conductor, "S1");
    EXPECT_NE(erf_gtest::abort_message([&] { m.add(ComponentKind::Hardware, "S1"); }).find("two components checked against their capacities are named 'S1'"), std::string::npos);
    const std::string neg = erf_gtest::abort_message([&] { m.record(c, -1.0, 1.0, 0.0); });
    EXPECT_NE(neg.find("conductor S1"), std::string::npos) << neg;
    EXPECT_FALSE(erf_gtest::abort_message([&] { m.record(c, std::numeric_limits<double>::quiet_NaN(), 1.0, 0.0); }).empty());
    EXPECT_FALSE(erf_gtest::abort_message([&] { m.record(c, 1.0, 0.0, 0.0); }).empty());
    EXPECT_FALSE(erf_gtest::abort_message([&] { m.record(c, 1.0, std::numeric_limits<double>::infinity(), 0.0); }).empty());
    EXPECT_FALSE(m.component(c).seen) << "a refused record leaves nothing behind";
}
