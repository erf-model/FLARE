// Unit tests of the towers' failures (conductors.tower_failure = release): which releases take a line off a failed
// tower (tower_releases) and the log of the failed towers (TowerFailureLog).
//
// TowerReleases.AStringOrBothSidesOfAClampedLineLeaveTheTower: with strings the tower's string; without, the span before
//     the tower at its end b and then the span after it at its end a; what was released before is left out.
// TowerReleases.ATowerOutsideTheLineIsRefused
// TowerFailureLog.RowsAndStateRoundTripInOrder: tower_failures.dat's rows, failed(), restore(state()) exact.
// TowerFailureLog.ASecondFailureOfATowerAndMalformedRecordsAreRefused

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "ERF_GTestThrowOnAbort.H"
#include "ERF_TowerFailures.H"

using erf_conductors::LineRelease;
using erf_conductors::TowerAttachment;
using erf_conductors::TowerFailure;
using erf_conductors::TowerFailureLog;
using erf_conductors::tower_releases;

TEST(TowerReleases, AStringOrBothSidesOfAClampedLineLeaveTheTower)
{
    TowerAttachment a;
    a.num_spans = 4;
    a.tower = 1;
    // clamped: span 2 (index 1) at its end b, then span 3 (index 2) at its end a
    auto r = tower_releases(a);
    ASSERT_EQ(r.size(), 2u);
    EXPECT_TRUE(r[0] == (LineRelease{LineRelease::Kind::SpanEnd, 1, 1})) << r[0].name();
    EXPECT_TRUE(r[1] == (LineRelease{LineRelease::Kind::SpanEnd, 2, 0})) << r[1].name();
    a.before_attached = false;
    r = tower_releases(a);
    ASSERT_EQ(r.size(), 1u);
    EXPECT_TRUE(r[0] == (LineRelease{LineRelease::Kind::SpanEnd, 2, 0})) << "only the side still on the tower";
    a.before_attached = true;
    a.after_attached = false;
    r = tower_releases(a);
    ASSERT_EQ(r.size(), 1u);
    EXPECT_TRUE(r[0] == (LineRelease{LineRelease::Kind::SpanEnd, 1, 1}));
    a.before_attached = false;
    EXPECT_TRUE(tower_releases(a).empty()) << "nothing left on the tower";
    // the string's flag does not count on a clamped line
    a.string_attached = false;
    a.before_attached = true;
    EXPECT_EQ(tower_releases(a).size(), 1u);
    // with strings: the string at this tower, whatever the spans' flags
    TowerAttachment s;
    s.strings = true;
    s.num_spans = 3;
    s.tower = 1;
    s.before_attached = false;
    s.after_attached = false;
    r = tower_releases(s);
    ASSERT_EQ(r.size(), 1u);
    EXPECT_TRUE(r[0] == (LineRelease{LineRelease::Kind::String, 1, 0})) << r[0].name();
    s.string_attached = false;
    EXPECT_TRUE(tower_releases(s).empty());
}

TEST(TowerReleases, ATowerOutsideTheLineIsRefused)
{
    TowerAttachment a;
    a.num_spans = 3;
    a.tower = 2;
    std::string msg = erf_gtest::abort_message([&] { tower_releases(a); });
    EXPECT_NE(msg.find("the tower is not one of the line's"), std::string::npos) << msg;
    a.tower = -1;
    msg = erf_gtest::abort_message([&] { tower_releases(a); });
    EXPECT_NE(msg.find("the tower is not one of the line's"), std::string::npos) << msg;
    a.tower = 0;
    a.num_spans = 1;
    msg = erf_gtest::abort_message([&] { tower_releases(a); });
    EXPECT_NE(msg.find("the tower is not one of the line's"), std::string::npos) << "a single span has no tower: " << msg;
}

TEST(TowerFailureLog, RowsAndStateRoundTripInOrder)
{
    TowerFailureLog g;
    EXPECT_TRUE(g.empty());
    EXPECT_STREQ(TowerFailureLog::log_header(), "time tower component releases");
    g.add(TowerFailure{0.1 + 0.2, "L1_t1", "L1_t1_footing_compression", 2});
    g.add(TowerFailure{4.5, "L1_t2", "L1_t2_m17", 0});
    ASSERT_EQ(g.size(), 2u);
    EXPECT_EQ(g.log_row(0), "0.3 L1_t1 L1_t1_footing_compression 2");
    EXPECT_EQ(g.log_row(1), "4.5 L1_t2 L1_t2_m17 0");
    EXPECT_TRUE(g.failed("L1_t1"));
    EXPECT_TRUE(g.failed("L1_t2"));
    EXPECT_FALSE(g.failed("L1_t3"));
    TowerFailureLog h;
    h.add(TowerFailure{0.1, "old", "old_footing_uplift", 1});
    EXPECT_EQ(h.restore(g.state()), "");
    ASSERT_EQ(h.size(), 2u) << "the log is replaced, not added to";
    EXPECT_FALSE(h.failed("old"));
    for (std::size_t i = 0; i < 2; ++i) {
        EXPECT_EQ(h.record(i).time, g.record(i).time) << "to the last bit";
        EXPECT_EQ(h.record(i).tower, g.record(i).tower);
        EXPECT_EQ(h.record(i).component, g.record(i).component);
        EXPECT_EQ(h.record(i).releases, g.record(i).releases);
    }
    EXPECT_EQ(h.restore(""), "");
    EXPECT_TRUE(h.empty());
}

TEST(TowerFailureLog, ASecondFailureOfATowerAndMalformedRecordsAreRefused)
{
    TowerFailureLog g;
    g.add(TowerFailure{1.0, "T1", "T1_footing_uplift", 2});
    std::string msg = erf_gtest::abort_message([&] { g.add(TowerFailure{2.0, "T1", "T1_m3", 0}); });
    EXPECT_NE(msg.find("tower T1 has failed before"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { g.add(TowerFailure{0.5, "T2", "T2_m3", 0}); });
    EXPECT_NE(msg.find("earlier than the last one logged"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { g.add(TowerFailure{2.0, "T 2", "T2_m3", 0}); });
    EXPECT_NE(msg.find("must be single words"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { g.add(TowerFailure{2.0, "T2", "T2_m3", -1}); });
    EXPECT_NE(msg.find("is malformed"), std::string::npos) << msg;
    EXPECT_EQ(g.size(), 1u);
    TowerFailureLog h;
    EXPECT_NE(h.restore("T1 T1_m3 2\n").find("malformed tower failure record"), std::string::npos) << "a field short";
    EXPECT_NE(h.restore("T1 T1_m3 -1 1.0\n").find("malformed tower failure record"), std::string::npos);
    EXPECT_NE(h.restore("T1 T1_m3 2 1.0 extra\n").find("malformed tower failure record"), std::string::npos);
    EXPECT_NE(h.restore("T1 T1_m3 2 2.0\nT2 T2_m3 2 1.0\n").find("malformed tower failure record"), std::string::npos) << "out of time order";
    EXPECT_NE(h.restore("T1 T1_m3 2 1.0\nT1 T1_m4 0 2.0\n").find("failed before"), std::string::npos);
    EXPECT_TRUE(h.empty()) << "a refused text leaves the log as it was";
}
