// Unit tests of the conductor lines' releases (conductors.failure_break): the choice of what lets go
// (decide_break), the log of the releases (BreakLog), and a ConductorLine after a release, on the MoorDyn
// library the build links (MoorDyn-C or the stub).
//
// DecideBreak.AConductorLetsGoAtItsEndWithTheHigherTension: the higher of the two end tensions, end a on a tie,
//     the one end still attached, nothing when both were released.
// DecideBreak.TheHardwareReleasesWhatItHolds: a dead end its span; a tower with strings its string; a clamped
//     tower the span with the higher tension there, the span before it on a tie, the one still attached.
// DecideBreak.AStringLetsGoOnceAndFootingsAndMembersReleaseNothing
// DecideBreak.IndicesOutsideTheLineAreRefused
// BreakLog.RowsAndStateRoundTripInOrder: breaks.dat's rows, releases_of() per line in order, restore(state())
//     exact to the last bit.
// BreakLog.MalformedAndRepeatedRecordsAreRefused
// LineRelease.AReleasedSpanHangsFromItsOtherEndAndItsSupportTakesNothing: 31 m of conductor from a support 30 m above
//     flat ground to one 1 m above it, still air, released at the low end: one more kinematics point, that support
//     takes nothing at once, and the tension at the high end is the weight of the 30 m that hang, w h (closed
//     form; the stub exactly, MoorDyn-C within one 1 m segment's weight below it), the last metre on the ground.
// LineRelease.ASpanShorterThanItsSupportsHeightHangsFree: a released span shorter than its support's height hangs clear
//     of the ground and its support carries its whole weight, w L.
// LineRelease.AStateAfterAReleaseFitsOnlyALineThatMakesItAgain: save() after a release; a line created from it
//     with the release made again continues with the first one; without it the load is refused; state() and
//     set_state() redo a step after a release.
// LineRelease.AClampedTowerKeepsOnlyTheSpanStillOnIt: two spans clamped at a tower; with the first released at
//     the tower, the tower's pull is the second span's end force alone.
// LineRelease.AReleasedStringDropsItsConductorAndItsTowerTakesNothing: three spans on strings; the first string
//     let go: its tower takes nothing and the conductor there comes down.
// LineRelease.ASpanOnStringsLeavesTheStringsLowerEnd: on a line with strings a span's end at a tower is on the string's
//     free lower end, not on the cross-arm: released there, the string stays on its tower with the other span.
// LineRelease.AReleaseAtAMovingTowersCrossArmContinuesFromASavedState: spans clamped to a tower that moves (MoorDyn's
//     coupled point): after a release there the tower's coupled force is the remaining span's pull, and a state saved
//     after the tower has moved on continues in a new line.
// LineRelease.TheStubStillLoadsAStateSavedBeforeItKnewReleases: the stub's earlier saved-state format (stub only).
// LineRelease.ReleasingTwiceOrOutsideTheLineIsRefused
// LineRelease.AStateTakenBeforeAReleaseIsRefusedAfterIt: set_state() with a state from before a release aborts.
// LineRelease.TheLargestTensionOfASpanIncludesItsEndB: max_tension() is at least both end tensions on an inclined span.
// LineRelease.TheWindOnTheOtherSpansIsKeptAcrossARelease: a release leaves the wind held for the other nodes.

#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "ERF_ConductorInputs.H"
#include "ERF_ConductorLine.H"
#include "ERF_GTestTempDir.H"
#include "ERF_GTestThrowOnAbort.H"
#include "ERF_LineBreaks.H"
#include "ERF_MoorDynInputWriter.H"
#include "ERF_MoorDynSystem.H"
#include "ERF_TowerInputs.H"

using erf_conductors::BreakDecision;
using erf_conductors::BreakLog;
using erf_conductors::BreakRecord;
using erf_conductors::ComponentKind;
using erf_conductors::ConductorInputs;
using erf_conductors::ConductorLine;
using erf_conductors::LineEnds;
using erf_conductors::LineInputs;
using erf_conductors::LineRelease;
using erf_conductors::decide_break;

namespace {

LineEnds ends (int spans, bool strings)
{
    LineEnds e;
    e.strings = strings;
    e.tension_a.assign(static_cast<std::size_t>(spans), 1000.0);
    e.tension_b.assign(static_cast<std::size_t>(spans), 1000.0);
    e.attached_a.assign(static_cast<std::size_t>(spans), 1);
    e.attached_b.assign(static_cast<std::size_t>(spans), 1);
    if (strings) { e.string_attached.assign(static_cast<std::size_t>(spans - 1), 1); }
    return e;
}

const LineRelease::Kind span_end = LineRelease::Kind::SpanEnd;
const LineRelease::Kind string_top = LineRelease::Kind::String;

} // namespace

TEST(DecideBreak, AConductorLetsGoAtItsEndWithTheHigherTension)
{
    LineEnds e = ends(3, false);
    e.tension_a[1] = 1200.0; e.tension_b[1] = 1300.0;
    BreakDecision d = decide_break(ComponentKind::Conductor, 1, e);
    ASSERT_TRUE(d.any);
    EXPECT_TRUE(d.release == (LineRelease{span_end, 1, 1})) << d.release.name();
    e.tension_a[1] = 1300.5;
    d = decide_break(ComponentKind::Conductor, 1, e);
    EXPECT_TRUE(d.release == (LineRelease{span_end, 1, 0})) << d.release.name();
    e.tension_b[1] = 1300.5;
    d = decide_break(ComponentKind::Conductor, 1, e);
    EXPECT_TRUE(d.release == (LineRelease{span_end, 1, 0})) << "end a on a tie: " << d.release.name();
    // a tie is within one part in a million: a solver's residual does not choose the end
    e.tension_b[1] = 1300.5 * (1.0 + 0.5e-6);
    d = decide_break(ComponentKind::Conductor, 1, e);
    EXPECT_TRUE(d.release == (LineRelease{span_end, 1, 0})) << "half the tolerance is still a tie: " << d.release.name();
    e.tension_b[1] = 1300.5 * (1.0 + 2.0e-6);
    d = decide_break(ComponentKind::Conductor, 1, e);
    EXPECT_TRUE(d.release == (LineRelease{span_end, 1, 1})) << "twice the tolerance is higher: " << d.release.name();
    e.tension_b[1] = 1300.5;
    // the one end still attached, whatever the tensions
    e.attached_a[1] = 0;
    e.tension_a[1] = 9000.0;
    d = decide_break(ComponentKind::Conductor, 1, e);
    ASSERT_TRUE(d.any);
    EXPECT_TRUE(d.release == (LineRelease{span_end, 1, 1})) << d.release.name();
    e.attached_b[1] = 0;
    EXPECT_FALSE(decide_break(ComponentKind::Conductor, 1, e).any) << "nothing left to release";
    // the other spans are not touched by span 1's state
    d = decide_break(ComponentKind::Conductor, 0, e);
    EXPECT_TRUE(d.release == (LineRelease{span_end, 0, 0}));
}

TEST(DecideBreak, TheHardwareReleasesWhatItHolds)
{
    // the dead ends: the span there, at that end
    LineEnds e = ends(3, false);
    BreakDecision d = decide_break(ComponentKind::Hardware, 0, e);
    ASSERT_TRUE(d.any);
    EXPECT_TRUE(d.release == (LineRelease{span_end, 0, 0})) << d.release.name();
    d = decide_break(ComponentKind::Hardware, 3, e);
    EXPECT_TRUE(d.release == (LineRelease{span_end, 2, 1})) << d.release.name();
    e.attached_b[2] = 0;
    EXPECT_FALSE(decide_break(ComponentKind::Hardware, 3, e).any) << "already released there";
    // a clamped tower (attachment 2, between spans 1 and 2): the side with the higher tension at the tower
    e.tension_b[1] = 1500.0; e.tension_a[2] = 1400.0;
    d = decide_break(ComponentKind::Hardware, 2, e);
    EXPECT_TRUE(d.release == (LineRelease{span_end, 1, 1})) << d.release.name();
    e.tension_a[2] = 1600.0;
    d = decide_break(ComponentKind::Hardware, 2, e);
    EXPECT_TRUE(d.release == (LineRelease{span_end, 2, 0})) << d.release.name();
    e.tension_a[2] = 1500.0;
    d = decide_break(ComponentKind::Hardware, 2, e);
    EXPECT_TRUE(d.release == (LineRelease{span_end, 1, 1})) << "the span before the tower on a tie: " << d.release.name();
    // the tensions at the spans' far ends do not count
    e.tension_a[1] = 9000.0; e.tension_b[2] = 9500.0; e.tension_a[2] = 1400.0;
    d = decide_break(ComponentKind::Hardware, 2, e);
    EXPECT_TRUE(d.release == (LineRelease{span_end, 1, 1})) << d.release.name();
    e.attached_b[1] = 0;
    d = decide_break(ComponentKind::Hardware, 2, e);
    EXPECT_TRUE(d.release == (LineRelease{span_end, 2, 0})) << "the one still on the tower: " << d.release.name();
    e.attached_a[2] = 0;
    EXPECT_FALSE(decide_break(ComponentKind::Hardware, 2, e).any);
    // a tower with strings: its string, whatever the spans' tensions
    LineEnds s = ends(3, true);
    s.tension_b[0] = 5000.0;
    d = decide_break(ComponentKind::Hardware, 1, s);
    ASSERT_TRUE(d.any);
    EXPECT_TRUE(d.release == (LineRelease{string_top, 0, 0})) << d.release.name();
    d = decide_break(ComponentKind::Hardware, 2, s);
    EXPECT_TRUE(d.release == (LineRelease{string_top, 1, 0})) << d.release.name();
    s.string_attached[1] = 0;
    EXPECT_FALSE(decide_break(ComponentKind::Hardware, 2, s).any);
    // its dead ends still release their spans
    d = decide_break(ComponentKind::Hardware, 3, s);
    EXPECT_TRUE(d.release == (LineRelease{span_end, 2, 1})) << d.release.name();
}

TEST(DecideBreak, AStringLetsGoOnceAndFootingsAndMembersReleaseNothing)
{
    LineEnds s = ends(3, true);
    BreakDecision d = decide_break(ComponentKind::Insulator, 1, s);
    ASSERT_TRUE(d.any);
    EXPECT_TRUE(d.release == (LineRelease{string_top, 1, 0}));
    EXPECT_EQ(d.release.name(), "string2");
    s.string_attached[1] = 0;
    EXPECT_FALSE(decide_break(ComponentKind::Insulator, 1, s).any);
    for (const ComponentKind k : {ComponentKind::FootingUplift, ComponentKind::FootingCompression, ComponentKind::Member}) {
        EXPECT_FALSE(decide_break(k, 0, s).any) << erf_conductors::kind_name(k);
        EXPECT_FALSE(decide_break(k, 77, s).any) << "their index is not a place on the line";
    }
    EXPECT_EQ((LineRelease{span_end, 2, 0}).name(), "span3_end_a");
    EXPECT_EQ((LineRelease{span_end, 0, 1}).name(), "span1_end_b");
}

TEST(DecideBreak, IndicesOutsideTheLineAreRefused)
{
    const LineEnds e = ends(2, false);
    std::string msg = erf_gtest::abort_message([&] { decide_break(ComponentKind::Conductor, 2, e); });
    EXPECT_NE(msg.find("no such span"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { decide_break(ComponentKind::Hardware, 3, e); });
    EXPECT_NE(msg.find("no such attachment"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { decide_break(ComponentKind::Insulator, 0, e); });
    EXPECT_NE(msg.find("no such insulator string"), std::string::npos) << "a line without strings: " << msg;
    LineEnds bad = e;
    bad.attached_b.pop_back();
    msg = erf_gtest::abort_message([&] { decide_break(ComponentKind::Conductor, 0, bad); });
    EXPECT_NE(msg.find("one tension and one attached flag per span end"), std::string::npos) << msg;
    bad = e;
    bad.string_attached.push_back(1);
    msg = erf_gtest::abort_message([&] { decide_break(ComponentKind::Conductor, 0, bad); });
    EXPECT_NE(msg.find("one flag per insulator string"), std::string::npos) << msg;
}

TEST(BreakLog, RowsAndStateRoundTripInOrder)
{
    BreakLog g;
    EXPECT_TRUE(g.empty());
    EXPECT_STREQ(BreakLog::log_header(), "time line released component tension");
    g.add(BreakRecord{1.875, "L1", LineRelease{span_end, 1, 1}, "L1_span2", 13900.123456789012});
    g.add(BreakRecord{1.875, "L2", LineRelease{string_top, 0, 0}, "L2_tower1", 0.1 + 0.2});
    g.add(BreakRecord{2.5, "L1", LineRelease{span_end, 0, 0}, "L1_end_a", 7000.0});
    ASSERT_EQ(g.size(), 3u);
    EXPECT_EQ(g.log_row(0), "1.875 L1 span2_end_b L1_span2 13900.12346");
    EXPECT_EQ(g.log_row(2), "2.5 L1 span1_end_a L1_end_a 7000");
    const auto l1 = g.releases_of("L1");
    ASSERT_EQ(l1.size(), 2u);
    EXPECT_TRUE(l1[0] == (LineRelease{span_end, 1, 1})) << "in the order made";
    EXPECT_TRUE(l1[1] == (LineRelease{span_end, 0, 0}));
    EXPECT_EQ(g.releases_of("L2").size(), 1u);
    EXPECT_TRUE(g.releases_of("L3").empty());
    BreakLog h;
    h.add(BreakRecord{0.5, "old", LineRelease{span_end, 0, 0}, "old", 1.0});
    EXPECT_EQ(h.restore(g.state()), "");
    ASSERT_EQ(h.size(), 3u) << "the log is replaced, not added to";
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_EQ(h.record(i).time, g.record(i).time);
        EXPECT_EQ(h.record(i).line, g.record(i).line);
        EXPECT_TRUE(h.record(i).release == g.record(i).release);
        EXPECT_EQ(h.record(i).component, g.record(i).component);
        EXPECT_EQ(h.record(i).tension, g.record(i).tension) << "to the last bit";
    }
    EXPECT_EQ(h.state(), g.state());
    EXPECT_EQ(h.restore(""), "");
    EXPECT_TRUE(h.empty()) << "a checkpoint without releases";
}

TEST(BreakLog, MalformedAndRepeatedRecordsAreRefused)
{
    BreakLog g;
    g.add(BreakRecord{1.0, "L1", LineRelease{span_end, 0, 1}, "L1", 5.0});
    std::string msg = erf_gtest::abort_message([&] { g.add(BreakRecord{2.0, "L1", LineRelease{span_end, 0, 1}, "L1_end_b", 5.0}); });
    EXPECT_NE(msg.find("line L1 was already released at span1_end_b"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { g.add(BreakRecord{0.5, "L1", LineRelease{span_end, 0, 0}, "L1", 5.0}); });
    EXPECT_NE(msg.find("earlier than the last one logged"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { g.add(BreakRecord{2.0, "L 1", LineRelease{span_end, 0, 0}, "L1", 5.0}); });
    EXPECT_NE(msg.find("must be single words"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { g.add(BreakRecord{2.0, "L1", LineRelease{span_end, 0, 2}, "L1", 5.0}); });
    EXPECT_NE(msg.find("is malformed"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { g.add(BreakRecord{2.0, "L1", LineRelease{span_end, 0, 0}, "L1", -5.0}); });
    EXPECT_NE(msg.find("is malformed"), std::string::npos) << "a negative tension: " << msg;
    EXPECT_EQ(g.size(), 1u);
    BreakLog h;
    EXPECT_NE(h.restore("L1 span 0 1 L1 1.0\n").find("malformed release record"), std::string::npos) << "a field short";
    EXPECT_NE(h.restore("L1 rope 0 1 L1 1.0 5.0\n").find("malformed release record"), std::string::npos);
    EXPECT_NE(h.restore("L1 span 0 2 L1 1.0 5.0\n").find("malformed release record"), std::string::npos);
    EXPECT_NE(h.restore("L1 string 0 1 L1 1.0 5.0\n").find("malformed release record"), std::string::npos) << "a string has no end b";
    EXPECT_NE(h.restore("L1 span 0 1 L1 1.0 5.0 extra\n").find("malformed release record"), std::string::npos);
    EXPECT_NE(h.restore("L1 span 0 1 L1 2.0 5.0\nL1 span 0 0 L1 1.0 5.0\n").find("malformed release record"), std::string::npos)
        << "out of time order";
    EXPECT_NE(h.restore("L1 span 0 1 L1 1.0 5.0\nL1 span 0 1 L1_end_b 2.0 5.0\n").find("given twice"), std::string::npos);
    EXPECT_TRUE(h.empty()) << "a refused text leaves the log as it was";
}

namespace {

const std::filesystem::path& gtest_scratch_root ()
{
    struct Root {
        std::filesystem::path p;
        ~Root () { std::error_code ec; std::filesystem::remove_all(p, ec); }
    };
    static const Root root{[] {
        const std::filesystem::path p = erf_gtest_temp_path("erf_gtest_linebreaks");
        std::filesystem::create_directories(p);
        return p;
    }()};
    return root.p;
}

constexpr double pi = 3.14159265358979323846;
constexpr double diameter = 0.0281, mass = 1.628, rho = 1.2, gravity = 9.81;
constexpr double weight = (mass - rho * 0.25 * pi * diameter * diameter) * gravity;   // N/m, net of buoyancy
// the same computation repeated (a restored state): roundoff only
constexpr double rtol = (std::is_same<amrex::Real, float>::value) ? 1.0e-5 : 1.0e-9;

ConductorInputs settings ()
{
    ConductorInputs in;
    in.air_density = rho;
    in.ground_contact = true;
    in.diagnostics_dir = (gtest_scratch_root() / "lines").string();
    return in;
}

// spans of 300 m along x at y = 500, every attachment 30 m above flat ground at z = 0
LineInputs level_line (const std::string& name, int spans, double insulator = 0.0)
{
    LineInputs s;
    s.name = name;
    s.end_a = {{100.0, 500.0, 30.0}};
    s.end_b = {{amrex::Real(100.0 + 300.0 * spans), 500.0, 30.0}};
    for (int k = 1; k < spans; ++k) { s.towers.push_back({{amrex::Real(100.0 + 300.0 * k), 500.0, 30.0}}); }
    s.lengths.assign(static_cast<std::size_t>(spans), 301.5);
    s.diameter = diameter; s.mass_per_length = mass; s.axial_stiffness = 3.0e7;
    if (insulator > 0.0) { s.insulator_length = insulator; s.insulator_mass = 60.0; }
    s.output_root = (gtest_scratch_root() / "lines" / name).string();
    return s;
}

// the line over flat ground: the seafloor file MoorDyn reads holds the ground raised by the conductor's radius
std::unique_ptr<ConductorLine> make (const LineInputs& s, const ConductorInputs& in, const std::string& saved = "",
                                     const std::vector<LineRelease>& releases = {})
{
    const std::string file = in.diagnostics_dir + "/" + s.name + ".moordyn.txt";
    std::filesystem::create_directories(in.diagnostics_dir);
    std::ofstream sf(erf_conductors::seafloor_file(file), std::ios::trunc);
    const double z = 0.5 * diameter;
    sf << erf_conductors::seafloor_text({-1000.0, 3000.0}, {-1000.0, 2000.0}, {z, z, z, z}, static_cast<double>(in.surface_offset));
    sf.close();
    return std::make_unique<ConductorLine>(s, in, gravity, file, saved, std::vector<amrex::Real>{}, releases);
}

void still (ConductorLine& c, double t, double dt)
{
    c.set_wind(std::vector<amrex::Real>(3 * static_cast<std::size_t>(c.num_kinematics_points()), amrex::Real(0.0)), t + 0.5 * dt);
    c.step(t, dt);
    std::vector<amrex::Real> h(c.num_nodes(), amrex::Real(0.0));
    c.set_ground_under_nodes(h);
}

double norm (const std::array<amrex::Real,3>& f) { return std::sqrt(static_cast<double>(f[0]*f[0] + f[1]*f[1] + f[2]*f[2])); }

} // namespace

TEST(LineRelease, AReleasedSpanHangsFromItsOtherEndAndItsSupportTakesNothing)
{
    // 31 m of conductor in 31 segments from a support 30 m above flat ground down to one 1 m above it, almost under
    // the first: let go at the low end it is already close to where it comes to rest, so it hardly swings (a level
    // span let go at one end swings for minutes in still air, and its support's load with it)
    const ConductorInputs in = settings();
    LineInputs s = level_line("hang", 1);
    s.end_b = {{101.0, 500.0, 1.0}};
    s.lengths = {31.0};
    s.segments = 31;
    auto c = make(s, in);
    const unsigned nk = c->num_kinematics_points();
    const double dt = 0.1;
    double t = 0.0;
    for (int n = 0; n < 5; ++n) { still(*c, t, dt); t += dt; }
    ASSERT_GT(norm(c->end_force(1)), 1.0) << "the low support holds the line's end";
    EXPECT_TRUE(c->span_end_attached(0, 0));
    EXPECT_TRUE(c->span_end_attached(0, 1));
    EXPECT_TRUE(c->releases().empty());

    c->release_span_end(0, 1);
    EXPECT_EQ(c->num_kinematics_points(), nk + 1) << "MoorDyn's new free point takes the wind too";
    EXPECT_EQ(c->kinematics_points().size(), 3u * (nk + 1));
    EXPECT_TRUE(c->span_end_attached(0, 0));
    EXPECT_FALSE(c->span_end_attached(0, 1));
    ASSERT_EQ(c->releases().size(), 1u);
    EXPECT_TRUE(c->releases()[0] == (LineRelease{LineRelease::Kind::SpanEnd, 0, 1}));
    EXPECT_EQ(norm(c->end_force(1)), 0.0) << "at once: the support holds nothing";
    EXPECT_GT(norm(c->end_force(0)), 300.0) << "the high end still holds the line";
    EXPECT_EQ(c->num_nodes(), 32u) << "the line's nodes are the same";

    // 30 m hang from end a and the last metre lies on the ground: the tension at end a is the weight of the
    // hanging length, w h
    const double h = 30.0 - 0.5 * diameter;
    for (int n = 0; n < 60; ++n) { still(*c, t, dt); t += dt; }
    if (erf_moordyn::is_stub()) {
        // the stub places the hanging length exactly
        EXPECT_NEAR(c->tension_a(), weight * h, 1.0e-9 * weight * h);
        EXPECT_NEAR(norm(c->end_force(0)), weight * h, 1.0e-9 * weight * h);
    } else {
        // MoorDyn-C lumps the line at its nodes: the tension at end a is the weight of the nodes that hang clear of
        // the ground, w l n with l = 1 m here. The node at the foot of the hanging length hovers millimetres above
        // the ground for the first seconds (30 nodes hang, 478.9 N) and then rests on it (29 nodes, 462.9 N), so
        // the tension lies within one segment's weight below w h, whenever it is read; 1 % for the swing left
        EXPECT_LE(c->tension_a(), 1.01 * weight * h) << "at most the weight of the hanging 30 m, " << weight * h << " N";
        EXPECT_GE(c->tension_a(), 0.99 * weight * (h - 1.0)) << "at least that of 29 m";
        // the support's load: that tension and the end node's own half segment
        EXPECT_NEAR(norm(c->end_force(0)) - c->tension_a(), 0.5 * weight, 0.1 * weight);
    }
    EXPECT_EQ(norm(c->end_force(1)), 0.0);
    const auto a = c->node_position(0);
    EXPECT_NEAR(a[0], 100.0, 1.0e-6);
    EXPECT_NEAR(a[2], 30.0, 1.0e-6) << "end a stays on its support";
    const auto b = c->node_position(31);
    EXPECT_NEAR(b[2], 0.5 * diameter, 2.0e-3) << "the released end lies on the ground";
    const auto mid = c->node_position(15);
    EXPECT_NEAR(mid[0], 100.0, 0.5) << "the hanging length is plumb under end a";
    EXPECT_NEAR(mid[2], 15.0, 0.1);
    const int on = c->nodes_on_ground(0);
    EXPECT_GE(on, 1) << "the released end";
    EXPECT_LE(on, 2) << "and the node at the foot of the hanging length";
    std::cout << "[ measured ] the tension at end a at rest " << c->tension_a() << " N, the support's load " << norm(c->end_force(0))
              << " N, w h = " << weight * h << " N; nodes on the ground " << on << "\n";
}

TEST(LineRelease, ASpanShorterThanItsSupportsHeightHangsFree)
{
    // 21.35 m of conductor in 3 segments from a support 30 m above the ground: let go at its low end it hangs clear of
    // the ground, and its support carries its whole weight, w L (in double, 21.35 * 3 / 3 is above 21.35: the last node's arc length must be the length itself)
    const ConductorInputs in = settings();
    LineInputs s = level_line("free", 1);
    s.end_b = {{100.2, 500.0, 10.0}};   // almost under end a: let go, the span drops 1.35 m and hardly swings
    s.lengths = {21.35};
    s.segments = 3;
    auto c = make(s, in);
    const double dt = 0.1;
    double t = 0.0;
    for (int n = 0; n < 3; ++n) { still(*c, t, dt); t += dt; }
    c->release_span_end(0, 1);
    for (int n = 0; n < 60; ++n) { still(*c, t, dt); t += dt; }
    // (MoorDyn-C: the support's load is every node's weight, the end node's half segment included, so w L whatever
    // the segment length; 2 % for what is left of the drop's bounce)
    const double tol = erf_moordyn::is_stub() ? 1.0e-9 : 0.02;
    EXPECT_NEAR(norm(c->end_force(0)), weight * 21.35, tol * weight * 21.35) << "the whole span's weight";
    const auto b = c->node_position(3);
    if (erf_moordyn::is_stub()) {
        EXPECT_NEAR(b[2], 30.0 - 21.35, 1.0e-9) << "the free end hangs 8.65 m above the ground";
        EXPECT_NEAR(b[0], 100.0, 1.0e-9) << "plumb under its support";
    } else {
        // MoorDyn-C's three 7 m links still swing out of their slack loop after 6 s (the end 1.3 m short of plumb)
        EXPECT_GT(b[2], 30.0 - 21.35 - 0.01);
        EXPECT_LT(b[2], 12.0) << "clear of the ground, below where it was held";
    }
    EXPECT_EQ(c->nodes_on_ground(0), 0);
    std::cout << "[ measured ] the support's load " << norm(c->end_force(0)) << " N, w L = " << weight * 21.35 << " N; the free end at z = "
              << b[2] << " m\n";
}

TEST(LineRelease, AStateAfterAReleaseFitsOnlyALineThatMakesItAgain)
{
    const ConductorInputs in = settings();
    const LineInputs sa = level_line("after_a", 2);
    auto a = make(sa, in);
    const double dt = 0.1;
    double t = 0.0;
    for (int n = 0; n < 3; ++n) { still(*a, t, dt); t += dt; }
    a->release_span_end(1, 0);
    // saved at once, before any step of the released span, and again after ten
    const std::string state0 = in.diagnostics_dir + "/after_a.state0";
    a->save(state0);
    const double t0 = t;
    for (int n = 0; n < 10; ++n) { still(*a, t, dt); t += dt; }
    const std::string state1 = in.diagnostics_dir + "/after_a.state1";
    a->save(state1);
    const double t1 = t;
    // a step redone from a State taken after the release
    const ConductorLine::State st = a->state();
    still(*a, t, dt);
    const std::vector<amrex::Real> once = a->node_positions();
    a->set_state(st);
    still(*a, t, dt);
    const std::vector<amrex::Real> twice = a->node_positions();
    ASSERT_EQ(once.size(), twice.size());
    for (std::size_t i = 0; i < once.size(); ++i) { EXPECT_EQ(once[i], twice[i]) << "component " << i; }
    t += dt;
    for (int n = 0; n < 9; ++n) { still(*a, t, dt); t += dt; }

    const std::vector<LineRelease> made{LineRelease{LineRelease::Kind::SpanEnd, 1, 0}};
    for (const bool late : {false, true}) {
        LineInputs sb = level_line(late ? "after_c" : "after_b", 2);
        auto b = make(sb, in, late ? state1 : state0, made);
        EXPECT_FALSE(b->span_end_attached(1, 0));
        EXPECT_EQ(b->num_kinematics_points(), a->num_kinematics_points());
        ASSERT_EQ(b->releases().size(), 1u);
        // MoorDyn's clock comes back with the state: the line goes on at the saved time
        double tb = late ? t1 : t0;
        for (int n = 0; n < (late ? 10 : 20); ++n) { still(*b, tb, dt); tb += dt; }
        ASSERT_NEAR(tb, t, 1.0e-12);
        const std::vector<amrex::Real> pa = a->node_positions(), pb = b->node_positions();
        ASSERT_EQ(pa.size(), pb.size());
        double worst = 0.0;
        for (std::size_t i = 0; i < pa.size(); ++i) { worst = std::max(worst, std::abs(static_cast<double>(pa[i] - pb[i]))); }
        EXPECT_LE(worst, rtol * 700.0) << (late ? "saved ten steps after the release" : "saved at the release");
        EXPECT_NEAR(b->tension_a(0), a->tension_a(0), rtol * a->tension_a(0));
    }
    // without the release made again MoorDyn's state does not fit the line
    const std::string msg = erf_gtest::abort_message([&] { make(level_line("after_d", 2), in, state1); });
    EXPECT_NE(msg.find("MoorDyn_Load"), std::string::npos) << msg;
    // and a new line starts whole
    const std::string msg2 = erf_gtest::abort_message([&] { make(level_line("after_e", 2), in, "", made); });
    EXPECT_NE(msg2.find("releases are made again only for a saved state"), std::string::npos) << msg2;
}

TEST(LineRelease, AClampedTowerKeepsOnlyTheSpanStillOnIt)
{
    const ConductorInputs in = settings();
    auto c = make(level_line("clamped", 2), in);
    still(*c, 0.0, 0.1);
    // two equal level spans: their horizontal pulls on the tower cancel and it carries one span's weight
    const auto whole = c->tower_force(0);
    EXPECT_LT(std::abs(whole[0]), 20.0);
    EXPECT_NEAR(whole[2], -weight * 301.5, 0.02 * weight * 301.5);
    c->release_span_end(0, 1);
    const auto left = c->tower_force(0);
    // span 2 alone: its horizontal tension towards +x and half its weight
    EXPECT_GT(left[0], 10000.0) << "the horizontal pull of the span still on the tower";
    EXPECT_NEAR(left[2], -0.5 * weight * 301.5, 0.02 * weight * 301.5);
    EXPECT_GT(norm(c->end_force(0)), 10000.0) << "span 1 is still on its dead end";
    c->release_span_end(1, 0);
    EXPECT_EQ(norm(c->tower_force(0)), 0.0) << "both spans gone: the tower takes nothing";
    EXPECT_EQ(c->num_kinematics_points(), c->num_nodes() + 3u + 2u + 1u) << "the nodes, three attachments, two new points, the origin";
    for (int n = 0; n < 3; ++n) { still(*c, 0.1 * (n + 1), 0.1); }
    EXPECT_EQ(norm(c->tower_force(0)), 0.0);
}

TEST(LineRelease, AReleasedStringDropsItsConductorAndItsTowerTakesNothing)
{
    const ConductorInputs in = settings();
    auto c = make(level_line("strings", 3, 2.5), in);
    const double dt = 0.1;
    double t = 0.0;
    for (int n = 0; n < 3; ++n) { still(*c, t, dt); t += dt; }
    ASSERT_GT(-c->tower_force(0)[2], 0.9 * weight * 301.5) << "tower 1 carries a span's weight";
    const unsigned at_tower = c->span_first_node(1);   // the conductor's node under tower 1
    const double z_before = c->node_position(at_tower)[2];
    EXPECT_NEAR(z_before, 27.5, 0.1);
    EXPECT_TRUE(c->string_attached(0));
    const unsigned nk = c->num_kinematics_points();
    c->release_string(0);
    EXPECT_FALSE(c->string_attached(0));
    EXPECT_TRUE(c->string_attached(1));
    EXPECT_EQ(c->num_kinematics_points(), nk + 1);
    EXPECT_EQ(norm(c->tower_force(0)), 0.0) << "at once";
    const double other = -c->tower_force(1)[2];
    const int nsteps = erf_moordyn::is_stub() ? 3 : 100;
    for (int n = 0; n < nsteps; ++n) { still(*c, t, dt); t += dt; }
    EXPECT_EQ(norm(c->tower_force(0)), 0.0);
    // 603 m of conductor over 600 m between the dead end and tower 2 sag four times a span's 13 m and more
    EXPECT_LT(c->node_position(at_tower)[2], z_before - 10.0) << "the conductor under tower 1 comes down";
    if (erf_moordyn::is_stub()) {
        // the stub puts that point on the elastic parabola of the two spans together, between the dead end (30 m) and
        // string 2's lower end (27.5 m): 603 m stretched by its tension H = w C^2 / (8 S), with sag S where
        // C + 8 S^2 / (3 C) = 603 (1 + H / EA). Solved here by Newton's method; the point is the chord's middle less S
        const double C = std::hypot(600.0, 2.5), Lt = 603.0, EA = 3.0e7;
        double S = 26.0;
        for (int it = 0; it < 50; ++it) {
            const double f = C + 8.0 * S * S / (3.0 * C) - Lt * (1.0 + weight * C * C / (8.0 * S * EA));
            const double df = 16.0 * S / (3.0 * C) + Lt * weight * C * C / (8.0 * S * S * EA);
            S -= f / df;
        }
        EXPECT_NEAR(c->node_position(at_tower)[2], 28.75 - S, 1.0e-6) << "sag " << S << " m";
    }
    // tower 2 now carries half of the 600 m either side of the lost support, not half of 300 m; the stub's string
    // counts half of each span's chord, which hardly changes
    EXPECT_GT(-c->tower_force(1)[2], (erf_moordyn::is_stub() ? 1.0 : 1.2) * other) << "tower 2 takes more of the conductor's weight";
    std::cout << "[ measured ] the conductor under the released string: " << z_before << " m -> " << c->node_position(at_tower)[2]
              << " m; tower 2's load " << other << " -> " << -c->tower_force(1)[2] << " N\n";
}

TEST(LineRelease, ASpanOnStringsLeavesTheStringsLowerEnd)
{
    const ConductorInputs in = settings();
    auto c = make(level_line("on_strings", 3, 2.5), in);
    const double dt = 0.1;
    double t = 0.0;
    for (int n = 0; n < 3; ++n) { still(*c, t, dt); t += dt; }
    const double before = -c->tower_force(0)[2];
    ASSERT_GT(before, 0.9 * weight * 301.5);
    const unsigned nk = c->num_kinematics_points();
    // span 2's end at tower 1 hangs from string 1's lower end
    c->release_span_end(1, 0);
    EXPECT_FALSE(c->span_end_attached(1, 0));
    EXPECT_TRUE(c->string_attached(0)) << "the string stays on its cross-arm";
    EXPECT_EQ(c->num_kinematics_points(), nk + 1);
    EXPECT_NEAR(-c->tower_force(0)[2], before, 1.0e-9 * before) << "the string's pull, until the line steps";
    const int nsteps = erf_moordyn::is_stub() ? 3 : 100;
    for (int n = 0; n < nsteps; ++n) { still(*c, t, dt); t += dt; }
    // tower 1 keeps its string and span 1. The stub places the rest state: about half a span's weight less than
    // before. MoorDyn-C's span 1 swings for minutes with the string, and the tower's load with it
    const double after = -c->tower_force(0)[2];
    EXPECT_GT(after, 0.2 * before);
    if (erf_moordyn::is_stub()) { EXPECT_LT(after, 0.8 * before); }
    // span 2 hangs from string 2's lower end at tower 2 and reaches the ground
    EXPECT_GE(c->nodes_on_ground(1), 5);
    std::cout << "[ measured ] tower 1's load " << before << " -> " << after << " N; span 2's nodes on the ground "
              << c->nodes_on_ground(1) << "\n";
}

TEST(LineRelease, AReleaseAtAMovingTowersCrossArmContinuesFromASavedState)
{
    ConductorInputs in = settings();
    erf_towers::TowerType type;
    type.name = "lat";
    type.base_width = 6.0; type.top_width = 1.5; type.solidity = 0.2; type.arm_length = 12.0;
    type.weight = 9.0e4; type.frequency = 2.0;
    in.tower_types = {type};
    LineInputs sa = level_line("moving_a", 2);
    sa.tower_type = "lat";
    auto a = make(sa, in);
    ASSERT_TRUE(a->towers_move());
    auto wind = [] (ConductorLine& c, double t) {
        c.set_wind(std::vector<amrex::Real>(3 * static_cast<std::size_t>(c.num_kinematics_points()), amrex::Real(0.0)), t);
    };
    auto ground = [] (ConductorLine& c) { c.set_ground_under_nodes(std::vector<amrex::Real>(c.num_nodes(), amrex::Real(0.0))); };
    // the tower's cross-arm moves at 0.2 m/s across the line
    const double dt = 0.1;
    const std::vector<amrex::Real> velocity{0.0, 0.2, 0.0};
    std::vector<amrex::Real> at{0.0, 0.0, 0.0};
    double t = 0.0;
    auto move = [&] (ConductorLine& c) { wind(c, t + 0.5 * dt); c.step_coupled(t, dt, at, velocity); ground(c); };
    move(*a); t += dt; at[1] += static_cast<amrex::Real>(0.2 * dt);
    // span 1 lets go of the moving cross-arm
    a->release_span_end(0, 1);
    const auto left = a->tower_force(0);
    EXPECT_GT(left[0], 10000.0) << "span 2 alone pulls the tower along the line";
    // the tower moves on with span 2, then the state is saved
    move(*a); t += dt; at[1] += static_cast<amrex::Real>(0.2 * dt);
    {
        // MoorDyn's force on the coupled point is the pull of what is still on it
        const auto f = a->coupled_force(0), F = a->tower_force(0);
        const double mag = norm(F);
        EXPECT_GT(mag, 5000.0);
        for (int d = 0; d < 3; ++d) { EXPECT_NEAR(f[d], F[d], 1.0e-6 * mag) << "direction " << d; }
    }
    move(*a); t += dt; at[1] += static_cast<amrex::Real>(0.2 * dt);
    const std::string state = in.diagnostics_dir + "/moving_a.state";
    a->save(state);
    // a new line at the tower's displacement, with the release made again
    LineInputs sb = level_line("moving_b", 2);
    sb.tower_type = "lat";
    const std::string file = in.diagnostics_dir + "/moving_b.moordyn.txt";
    {
        std::ofstream sf(erf_conductors::seafloor_file(file), std::ios::trunc);
        const double z = 0.5 * diameter;
        sf << erf_conductors::seafloor_text({-1000.0, 3000.0}, {-1000.0, 2000.0}, {z, z, z, z}, static_cast<double>(in.surface_offset));
    }
    ConductorLine b(sb, in, gravity, file, state, at, {LineRelease{LineRelease::Kind::SpanEnd, 0, 1}});
    EXPECT_FALSE(b.span_end_attached(0, 1));
    for (int n = 0; n < 5; ++n) {
        move(*a);
        wind(b, t + 0.5 * dt); b.step_coupled(t, dt, at, velocity); ground(b);
        t += dt; at[1] += static_cast<amrex::Real>(0.2 * dt);
    }
    const std::vector<amrex::Real> pa = a->node_positions(), pb = b.node_positions();
    ASSERT_EQ(pa.size(), pb.size());
    double worst = 0.0;
    for (std::size_t i = 0; i < pa.size(); ++i) { worst = std::max(worst, std::abs(static_cast<double>(pa[i] - pb[i]))); }
    EXPECT_LE(worst, rtol * 700.0) << "the restored line follows the first";
    const auto fa = a->coupled_force(0), fb = b.coupled_force(0);
    for (int d = 0; d < 3; ++d) { EXPECT_NEAR(fb[d], fa[d], rtol * std::max(1.0, norm(fa))) << "direction " << d; }
}

TEST(LineRelease, TheStubStillLoadsAStateSavedBeforeItKnewReleases)
{
    if (!erf_moordyn::is_stub()) { GTEST_SKIP() << "the stub's own saved-state format"; }
    // the stub's state file before MoorDyn_BreakLine: the tag moordyn-stub-state and no drop direction per line
    const ConductorInputs in = settings();
    auto a = make(level_line("old_a", 1), in);
    for (int n = 0; n < 3; ++n) { still(*a, 0.1 * n, 0.1); }
    const std::string now = in.diagnostics_dir + "/old_a.state", old = in.diagnostics_dir + "/old_a.state1";
    a->save(now);
    {
        std::ifstream f(now);
        std::ofstream o(old);
        std::string tag, rest, row;
        f >> tag;
        ASSERT_EQ(tag, "moordyn-stub-state-2");
        std::getline(f, rest);
        o << "moordyn-stub-state" << rest << "\n";
        // the line's row: phi, phi_across, length, the two drop numbers, then the nodes
        std::getline(f, row);
        std::istringstream r(row);
        std::string v;
        for (int i = 0; r >> v; ++i) { if (i != 3 && i != 4) { o << v << " "; } }
        o << "\n";
        while (std::getline(f, row)) { o << row << "\n"; }
    }
    auto b = make(level_line("old_b", 1), in, old);
    const std::vector<amrex::Real> pa = a->node_positions(), pb = b->node_positions();
    ASSERT_EQ(pa.size(), pb.size());
    for (std::size_t i = 0; i < pa.size(); ++i) { EXPECT_NEAR(pb[i], pa[i], rtol * 700.0) << "component " << i; }
}

TEST(LineRelease, ReleasingTwiceOrOutsideTheLineIsRefused)
{
    const ConductorInputs in = settings();
    auto c = make(level_line("twice", 2), in);
    c->release_span_end(0, 1);
    std::string msg = erf_gtest::abort_message([&] { c->release_span_end(0, 1); });
    EXPECT_NE(msg.find("conductors.twice: span 1 was already released at its end b"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { c->release_span_end(2, 0); });
    EXPECT_NE(msg.find("conductors.twice: no span 3"), std::string::npos) << "counted from 1, as breaks.dat does: " << msg;
    msg = erf_gtest::abort_message([&] { c->release_span_end(0, 2); });
    EXPECT_NE(msg.find("conductors.twice: no span 1 with an end 2"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { c->release_string(0); });
    EXPECT_NE(msg.find("conductors.twice: no insulator string 1"), std::string::npos) << "a line without strings: " << msg;
    EXPECT_EQ(c->releases().size(), 1u);
    auto s = make(level_line("twice_s", 3, 2.5), in);
    s->release_string(1);
    msg = erf_gtest::abort_message([&] { s->release_string(1); });
    EXPECT_NE(msg.find("conductors.twice_s: insulator string 2 was already released"), std::string::npos) << msg;
}

TEST(LineRelease, AStateTakenBeforeAReleaseIsRefusedAfterIt)
{
    // the system has one more point after a release: a state from before it does not fit
    const ConductorInputs in = settings();
    auto c = make(level_line("stale_state", 2), in);
    still(*c, 0.0, 0.1);
    const auto before = c->state();
    c->release_span_end(0, 0);
    const std::string msg = erf_gtest::abort_message([&] { c->set_state(before); });
    EXPECT_FALSE(msg.empty()) << "a state taken before the release was accepted";
    // one taken after it is
    const auto after = c->state();
    still(*c, 0.1, 0.1);
    c->set_state(after);
    still(*c, 0.1, 0.1);
}

TEST(LineRelease, TheLargestTensionOfASpanIncludesItsEndB)
{
    // an inclined span whose end b is 40 m higher: its tension is largest at end b, the last node of the line,
    // which MoorDyn-C's own largest-tension call leaves out
    const ConductorInputs in = settings();
    LineInputs s = level_line("max_at_b", 1);
    s.end_b = {{400.0, 500.0, 70.0}};
    s.lengths.assign(1, 304.0);
    auto c = make(s, in);
    still(*c, 0.0, 0.1);
    EXPECT_GE(c->max_tension(0), c->tension_b(0));
    EXPECT_GE(c->max_tension(0), c->tension_a(0));
    if (!erf_moordyn::is_stub()) { EXPECT_GT(c->tension_b(0), c->tension_a(0) + 0.5 * weight * 40.0) << "the higher end carries more"; }
}

TEST(LineRelease, TheWindOnTheOtherSpansIsKeptAcrossARelease)
{
    // two spans clamped at a tower in a crosswind: releasing span 1 must leave span 2 blown out as it was
    const ConductorInputs in = settings();
    auto c = make(level_line("wind_kept", 2), in);
    auto blow = [&] (double t, double dt) {
        std::vector<amrex::Real> uvw(3 * static_cast<std::size_t>(c->num_kinematics_points()), amrex::Real(0.0));
        for (std::size_t p = 0; p < uvw.size() / 3; ++p) { uvw[3 * p + 1] = amrex::Real(20.0); }
        c->set_wind(uvw, t + 0.5 * dt);
        c->step(t, dt);
        std::vector<amrex::Real> h(c->num_nodes(), amrex::Real(0.0));
        c->set_ground_under_nodes(h);
    };
    const int nsteps = erf_moordyn::is_stub() ? 3 : 300;
    double t = 0.0;
    for (int n = 0; n < nsteps; ++n) { blow(t, 0.1); t += 0.1; }
    const double offset = static_cast<double>(c->mid_offset(1));
    ASSERT_GT(std::abs(offset), 1.0) << "span 2 is blown out";
    c->release_span_end(0, 0);
    // the wind MoorDyn holds for the nodes it had: unchanged by the release
    const unsigned mid = c->span_first_node(1) + (c->span_num_nodes(1) - 1) / 2;
    EXPECT_NEAR(static_cast<double>(c->wind_at_point(mid)[1]), 20.0, 1.0e-12);
    EXPECT_NEAR(static_cast<double>(c->mid_offset(1)), offset, 1.0e-9 * std::abs(offset)) << "until the line steps";
    // a step with the wind MoorDyn holds, none set again: span 2 stays blown out. The stub moves a span towards
    // the shape its wind gives, so after three steps it is still moving out, and with a wind zeroed at the release
    // it would move back; MoorDyn-C moves it little in 0.1 s
    c->step(t, 0.1);
    const double after = static_cast<double>(c->mid_offset(1));
    if (erf_moordyn::is_stub()) {
        EXPECT_GT(std::abs(after), std::abs(offset)) << "span 2 still moves out in the wind held across the release";
    } else {
        EXPECT_NEAR(after, offset, 0.2 * std::abs(offset)) << "span 2 after a step in the wind held across the release";
    }
}
