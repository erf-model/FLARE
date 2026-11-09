// Unit tests of GroundContactLog, the record of which conductor spans lie on the ground (conductors.ground_contact).
//
// AFirstTouchIsReportedOnceAndTheMostIsKept: record() returns true at a span's first record with a node on the
//     ground only; the first time, the count now and the most (with its first time) follow the records.
// CountsOutsideTheSpanAndDuplicateNamesAreRefused: a negative count, more nodes than the span has, a second span of
//     the same name, an empty name and a span of fewer than two nodes abort, naming the span.
// TheTableListsEverySpan: write_table() writes the header and one row per span, empty times for a span never on the ground.
// RecordsRoundTripAndUnmatchedSpansAreDropped: restore(state()) gives back every field; a span the log does not hold,
//     or held with another node count, is skipped and named; a malformed line is refused with its number.

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "ERF_GTestTempDir.H"
#include "ERF_GTestThrowOnAbort.H"
#include "ERF_GroundContact.H"

using erf_conductors::GroundContactLog;

TEST(GroundContactLog, AFirstTouchIsReportedOnceAndTheMostIsKept)
{
    GroundContactLog g;
    const auto s = g.add("S1_span1", 21);
    EXPECT_FALSE(g.record(s, 0, 0.5));
    EXPECT_FALSE(g.span(s).touched);
    EXPECT_TRUE(g.record(s, 3, 1.0)) << "the first touch";
    EXPECT_FALSE(g.record(s, 7, 1.5)) << "not again";
    EXPECT_FALSE(g.record(s, 7, 2.0));
    EXPECT_FALSE(g.record(s, 0, 2.5)) << "lifted off and back on is not a first touch";
    EXPECT_FALSE(g.record(s, 2, 3.0));
    const auto& r = g.span(s);
    EXPECT_TRUE(r.touched);
    EXPECT_EQ(r.first_time, 1.0);
    EXPECT_EQ(r.now, 2);
    EXPECT_EQ(r.most, 7);
    EXPECT_EQ(r.most_time, 1.5) << "the first record with the most";
}

TEST(GroundContactLog, CountsOutsideTheSpanAndDuplicateNamesAreRefused)
{
    GroundContactLog g;
    const auto s = g.add("S1_span1", 21);
    std::string msg = erf_gtest::abort_message([&] { g.record(s, -1, 0.5); });
    EXPECT_NE(msg.find("-1 nodes of the span 'S1_span1' on the ground, of 21"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { g.record(s, 22, 0.5); });
    EXPECT_NE(msg.find("22 nodes of the span 'S1_span1'"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { g.add("S1_span1", 21); });
    EXPECT_NE(msg.find("two spans are named 'S1_span1'"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { g.add("", 21); });
    EXPECT_NE(msg.find("a span needs a name"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { g.add("S1_span2", 1); });
    EXPECT_NE(msg.find("'S1_span2' needs at least two nodes"), std::string::npos) << msg;
    EXPECT_EQ(g.size(), 1u);
}

TEST(GroundContactLog, TheTableListsEverySpan)
{
    GroundContactLog g;
    const auto a = g.add("S1_span1", 21);
    g.add("S1_span2", 11);
    g.record(a, 4, 0.25);
    const auto dir = erf_gtest_temp_path("erf_gtest_ground_table");
    std::filesystem::create_directories(dir);
    const std::string f = (dir / "ground_contact.csv").string();
    g.write_table(f);
    std::ifstream in(f);
    std::vector<std::string> rows;
    for (std::string l; std::getline(in, l);) { rows.push_back(l); }
    std::filesystem::remove_all(dir);
    ASSERT_EQ(rows.size(), 3u);
    EXPECT_EQ(rows[0], "span,nodes,nodes_on_ground,most_on_ground,most_time,first_touch_time");
    EXPECT_EQ(rows[1], "S1_span1,21,4,4,0.25,0.25");
    EXPECT_EQ(rows[2], "S1_span2,11,0,0,,");
}

TEST(GroundContactLog, RecordsRoundTripAndUnmatchedSpansAreDropped)
{
    GroundContactLog g;
    const auto a = g.add("S1_span1", 21);
    const auto b = g.add("S1_span2", 11);
    g.record(a, 5, 0.1);
    g.record(a, 3, 0.30000000000000004);
    g.record(b, 0, 0.30000000000000004);
    GroundContactLog h;
    h.add("S1_span1", 21);
    h.add("S1_span2", 11);
    std::vector<std::string> dropped;
    EXPECT_EQ(h.restore(g.state(), dropped), "");
    EXPECT_TRUE(dropped.empty());
    for (std::size_t i = 0; i < 2; ++i) {
        EXPECT_EQ(h.span(i).now, g.span(i).now);
        EXPECT_EQ(h.span(i).most, g.span(i).most);
        EXPECT_EQ(h.span(i).most_time, g.span(i).most_time);
        EXPECT_EQ(h.span(i).touched, g.span(i).touched);
        EXPECT_EQ(h.span(i).first_time, g.span(i).first_time);
    }
    // a log whose second span has another node count, and without the first: both dropped
    GroundContactLog k;
    k.add("S1_span2", 21);
    EXPECT_EQ(k.restore(g.state(), dropped), "");
    ASSERT_EQ(dropped.size(), 2u);
    EXPECT_EQ(dropped[0], "S1_span1");
    EXPECT_EQ(dropped[1], "S1_span2");
    EXPECT_FALSE(k.span(0).touched);
    EXPECT_EQ(k.restore("span S1_span2 21 3 2 0.1 1 0.1\n", dropped), "malformed ground-contact record on line 1")
        << "fewer on the ground at most than now";
    EXPECT_EQ(k.restore("\nspan S1_span2 21 3\n", dropped), "malformed ground-contact record on line 2");
}
