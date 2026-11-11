// Unit tests of a conductor section's state at rest after a release (conductors.failure_break,
// conductors.tower_failure), against the independent reference Structures/Tools/broken_wire/broken_wire.py: the
// balance of forces of heavy elastic cables at the strings' lower ends, with no code shared with ERF or
// MoorDyn. The reference numbers below are that tool's, for MoorDyn's lumped chain of 20 links per span (its
// "segments" = 20) and for continuous cables; Structures/Tools/broken_wire/tests/test_broken_wire.py
// (class Reference) recomputes them. The tests need MoorDyn-C: the stub has no dynamics and places a released
// span without moving the strings. In still air a section sways about its state at rest for minutes, so the
// pulls compared are means over 60 s; they must be within 3e-4 of the largest pull of the lumped chain's.
//
// BrokenWire.TheWholeSectionHangsAsTheForceBalanceSays: before any release, the towers' and the dead end's pulls.
// BrokenWire.ASpanCutAtATowerLeavesTheResidualTension: span 1 released at tower 1 on three sections (unequal
//     spans; equal spans on 2.5 m strings; equal spans on 4 m strings): at rest the string at tower 1 leans into
//     span 2, the tower carries the residual tension along the line, and both towers' and the far dead end's
//     pulls are the reference's.
// BrokenWire.ASpanCutAtItsDeadEndHangsFromTheStringAndLoadsIt: span 1 released at its dead end hangs from the
//     foot of the string at tower 1 and reaches the ground: tower 1 also carries the hanging length's weight.
// BrokenWire.AReleasedStringPassesItsSpansToTheNextTower: the string at tower 1 released: spans 1 and 2 hang as
//     one and tower 2 carries the reference's load.

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "ERF_ConductorInputs.H"
#include "ERF_ConductorLine.H"
#include "ERF_GTestTempDir.H"
#include "ERF_LineBreaks.H"
#include "ERF_MoorDynInputWriter.H"
#include "ERF_MoorDynSystem.H"

using erf_conductors::ConductorInputs;
using erf_conductors::ConductorLine;
using erf_conductors::LineInputs;

namespace {

const std::filesystem::path& gtest_scratch_root ()
{
    struct Root {
        std::filesystem::path p;
        ~Root () { std::error_code ec; std::filesystem::remove_all(p, ec); }
    };
    static const Root root{[] {
        const std::filesystem::path p = erf_gtest_temp_path("erf_gtest_brokenwire");
        std::filesystem::create_directories(p);
        return p;
    }()};
    return root.p;
}

constexpr double diameter = 0.0281, mass = 1.628, rho = 1.2, gravity = 9.81;

ConductorInputs settings ()
{
    ConductorInputs in;
    in.air_density = rho;
    in.ground_contact = true;
    in.diagnostics_dir = (gtest_scratch_root() / "lines").string();
    return in;
}

// a section along x at y = 500 over flat ground at z = 0: dead ends at x = 100 and 1000 m, two towers, 795 kcmil
// ACSR on strings of 60 kg. The drag coefficient is raised so that the section comes to rest within the test:
// in still air the drag takes no part in the state at rest
LineInputs section (const std::string& name, double tower2_x, const std::vector<amrex::Real>& lengths, double z_ends,
                    double insulator)
{
    LineInputs s;
    s.name = name;
    s.end_a = {{100.0, 500.0, amrex::Real(z_ends)}};
    s.end_b = {{1000.0, 500.0, amrex::Real(z_ends)}};
    s.towers.push_back({{400.0, 500.0, 60.0}});
    s.towers.push_back({{amrex::Real(tower2_x), 500.0, 60.0}});
    s.lengths = lengths;
    s.diameter = diameter; s.mass_per_length = mass; s.axial_stiffness = 3.0e7;
    s.drag_coefficient = 40.0;
    s.insulator_length = insulator; s.insulator_mass = 60.0;
    s.output_root = (gtest_scratch_root() / "lines" / name).string();
    return s;
}

std::unique_ptr<ConductorLine> make (const LineInputs& s, const ConductorInputs& in)
{
    const std::string file = in.diagnostics_dir + "/" + s.name + ".moordyn.txt";
    std::filesystem::create_directories(in.diagnostics_dir);
    std::ofstream sf(erf_conductors::seafloor_file(file), std::ios::trunc);
    const double z = 0.5 * diameter;
    sf << erf_conductors::seafloor_text({-1000.0, 3000.0}, {-1000.0, 2000.0}, {z, z, z, z}, static_cast<double>(in.surface_offset));
    sf.close();
    return std::make_unique<ConductorLine>(s, in, gravity, file, "", std::vector<amrex::Real>{}, std::vector<erf_conductors::LineRelease>{});
}

void still (ConductorLine& c, double t, double dt)
{
    c.set_wind(std::vector<amrex::Real>(3 * static_cast<std::size_t>(c.num_kinematics_points()), amrex::Real(0.0)), t + 0.5 * dt);
    c.step(t, dt);
    std::vector<amrex::Real> h(c.num_nodes(), amrex::Real(0.0));
    c.set_ground_under_nodes(h);
}

// the pulls the test reads: on the two towers and on dead end b, along the line (x) and vertical (z)
struct Pulls { double t1x, t1z, t2x, t2z, bx, bz; };

Pulls pulls (const ConductorLine& c)
{
    const auto a = c.tower_force(0), b = c.tower_force(1), e = c.end_force(1);
    return {static_cast<double>(a[0]), static_cast<double>(a[2]), static_cast<double>(b[0]), static_cast<double>(b[2]),
            static_cast<double>(e[0]), static_cast<double>(e[2])};
}

// advance in still air, at most to t_max, and return the time. In still air the section keeps swaying about its
// state at rest for minutes (the drag falls with the square of the speed), so the pulls compared are their means
// over 60 s, which a small sway about the state at rest does not move: 'mean' receives the last such mean, and the
// advance ends when two successive means differ by less than 'settled' N
double to_rest (ConductorLine& c, double t, double t_max, Pulls& mean, double settled = 3.0)
{
    const double dt = 0.1;
    const int window = 600;
    Pulls last{};
    bool have = false;
    while (t < t_max) {
        double s[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        for (int n = 0; n < window; ++n) {
            still(c, t, dt); t += dt;
            const Pulls p = pulls(c);
            s[0] += p.t1x; s[1] += p.t1z; s[2] += p.t2x; s[3] += p.t2z; s[4] += p.bx; s[5] += p.bz;
        }
        mean = {s[0] / window, s[1] / window, s[2] / window, s[3] / window, s[4] / window, s[5] / window};
        const double change = std::max({std::abs(mean.t1x - last.t1x), std::abs(mean.t1z - last.t1z), std::abs(mean.t2x - last.t2x),
                                        std::abs(mean.t2z - last.t2z), std::abs(mean.bx - last.bx), std::abs(mean.bz - last.bz)});
        if (have && change < settled) { break; }
        last = mean; have = true;
    }
    return t;
}

void expect_pulls (const Pulls& got, const Pulls& lumped, const Pulls& continuous, double tol_lumped, const std::string& what)
{
    // against the lumped chain MoorDyn integrates: tol_lumped of the largest pull; against continuous cables: the
    // difference the links' length makes, 1 % of the largest pull
    const double scale = std::max({std::abs(lumped.t1x), std::abs(lumped.t1z), std::abs(lumped.t2z), std::abs(lumped.bx)});
    const double a[6] = {got.t1x, got.t1z, got.t2x, got.t2z, got.bx, got.bz};
    const double l[6] = {lumped.t1x, lumped.t1z, lumped.t2x, lumped.t2z, lumped.bx, lumped.bz};
    const double c[6] = {continuous.t1x, continuous.t1z, continuous.t2x, continuous.t2z, continuous.bx, continuous.bz};
    const char* names[6] = {"tower 1 along the line", "tower 1 vertical", "tower 2 along the line", "tower 2 vertical",
                            "dead end b along the line", "dead end b vertical"};
    for (int i = 0; i < 6; ++i) {
        EXPECT_NEAR(a[i], l[i], tol_lumped * scale) << what << ": " << names[i] << " against the lumped chain";
        EXPECT_NEAR(a[i], c[i], 0.01 * scale) << what << ": " << names[i] << " against continuous cables";
    }
    std::cout << "[ measured ] " << what << ": tower 1 " << got.t1x << " " << got.t1z << ", tower 2 " << got.t2x << " " << got.t2z
              << ", dead end b " << got.bx << " " << got.bz << " N (reference, lumped: " << lumped.t1x << " " << lumped.t1z << ", "
              << lumped.t2x << " " << lumped.t2z << ", " << lumped.bx << " " << lumped.bz << ")\n";
}

// spans of 300, 500 and 100 m, slack in proportion to the cube of the span, every attachment 60 m up
const std::vector<amrex::Real> unequal{301.5, 506.944, 100.056};
// three spans of 300 m; the dead ends as high as the strings' feet
const std::vector<amrex::Real> equal{301.5, 301.5, 301.5};

} // namespace

TEST(BrokenWire, TheWholeSectionHangsAsTheForceBalanceSays)
{
    if (erf_moordyn::is_stub()) { GTEST_SKIP() << "needs the real MoorDyn-C: the stub hangs a parabola"; }
    const ConductorInputs in = settings();
    auto c = make(section("bw_whole", 900.0, unequal, 60.0, 2.5), in);
    Pulls got{};
    to_rest(*c, 0.0, 600.0, got);
    expect_pulls(got, {109.86, -6926.60, -23.17, -5093.22, -13530.64, -1137.39},
                 {109.85, -6926.47, -23.24, -5092.83, -13546.39, -1137.78}, 3.0e-4, "whole");
}

TEST(BrokenWire, ASpanCutAtATowerLeavesTheResidualTension)
{
    if (erf_moordyn::is_stub()) { GTEST_SKIP() << "needs the real MoorDyn-C: the stub does not move the strings"; }
    struct Case { const char* name; double tower2_x; const std::vector<amrex::Real>* lengths; double z_ends, insulator; Pulls lumped, continuous; double angle; };
    const std::vector<Case> cases{
        // span 2 keeps 11.77 kN of 13.55 kN (87 %); the string at tower 1 leans 69.6 degrees
        {"bw_cut_unequal", 900.0, &unequal, 60.0, 2.5, {11765.90, -4673.32, 23.08, -5096.71, -11788.99, -1093.96},
         {11779.34, -4673.38, 23.02, -5096.31, -11802.36, -1094.30}, 69.585},
        // equal spans: span 2 keeps 9.34 kN of 13.24 kN (71 %), span 3 10.95 kN; the strings lean 73.6 and 17.7 degrees
        {"bw_cut_equal", 700.0, &equal, 57.5, 2.5, {9342.85, -3047.20, 1609.45, -5350.80, -10952.30, -2402.14},
         {9354.33, -3047.29, 1610.02, -5350.72, -10964.35, -2402.13}, 73.582},
        // 4 m strings give more: span 2 keeps 8.26 kN (62 %)
        {"bw_cut_long_strings", 700.0, &equal, 56.0, 4.0, {8258.04, -3064.66, 1678.85, -5334.10, -9936.89, -2399.59},
         {8268.14, -3064.78, 1679.40, -5334.00, -9947.54, -2399.57}, 71.452}};
    const ConductorInputs in = settings();
    for (const auto& k : cases) {
        auto c = make(section(k.name, k.tower2_x, *k.lengths, k.z_ends, k.insulator), in);
        Pulls got{};
        double t = to_rest(*c, 0.0, 600.0, got);
        const double whole = -got.bx;
        // span 1's end at tower 1: the string's foot
        c->release_span_end(0, 1);
        t = to_rest(*c, t, t + 1200.0, got);
        expect_pulls(got, k.lumped, k.continuous, 3.0e-4, k.name);
        EXPECT_LT(-got.bx, 0.95 * whole) << k.name << ": the tension along the line has fallen";
        EXPECT_NEAR(static_cast<double>(c->insulator_swing(0)) * 180.0 / 3.14159265358979323846, k.angle, 0.2) << k.name << ": the string at tower 1 leans into span 2 (degrees)";
        EXPECT_NEAR(got.t1x, -got.t2x - got.bx, 3.0e-4 * got.t1x) << k.name << ": tower 1 carries what the other supports hold against";
    }
}

TEST(BrokenWire, ASpanCutAtItsDeadEndHangsFromTheStringAndLoadsIt)
{
    if (erf_moordyn::is_stub()) { GTEST_SKIP() << "needs the real MoorDyn-C: the stub does not move the strings"; }
    const ConductorInputs in = settings();
    auto c = make(section("bw_cut_dead_end", 900.0, unequal, 60.0, 2.5), in);
    Pulls got{};
    double t = to_rest(*c, 0.0, 600.0, got);
    c->release_span_end(0, 0);
    t = to_rest(*c, t, t + 1200.0, got);
    // the string's foot is 58.4 m up: three nodes of span 1 hang clear in 15 m links and load the foot with their
    // weight and the held end's node (842.3 N), against 59.0 m of a continuous conductor (941.2 N)
    const Pulls lumped{11801.89, -5512.27, 22.01, -5099.15, -11823.90, -1094.83};
    const double scale = 11801.89;
    EXPECT_NEAR(got.t1x, lumped.t1x, 3.0e-4 * scale);
    EXPECT_NEAR(got.t1z, lumped.t1z, 3.0e-4 * scale) << "tower 1 also carries the nodes of span 1 that hang clear";
    EXPECT_NEAR(got.t2z, lumped.t2z, 3.0e-4 * scale);
    EXPECT_NEAR(got.bx, lumped.bx, 3.0e-4 * scale);
    EXPECT_GE(c->nodes_on_ground(0), 10) << "the rest of span 1 lies on the ground";
    std::cout << "[ measured ] cut at the dead end: tower 1 " << got.t1x << " " << got.t1z << ", tower 2 " << got.t2x << " " << got.t2z
              << ", dead end b " << got.bx << " " << got.bz << " N; span 1's nodes on the ground " << c->nodes_on_ground(0) << "\n";
}

TEST(BrokenWire, AReleasedStringPassesItsSpansToTheNextTower)
{
    if (erf_moordyn::is_stub()) { GTEST_SKIP() << "needs the real MoorDyn-C: the stub's string load hardly changes"; }
    const ConductorInputs in = settings();
    auto c = make(section("bw_string", 900.0, unequal, 60.0, 2.5), in);
    Pulls got{};
    double t = to_rest(*c, 0.0, 600.0, got);
    const double before = -got.t2z;
    c->release_string(0);
    t = to_rest(*c, t, t + 1200.0, got);
    // spans 1 and 2 hang as one 800 m span, 5.7 m clear of the ground at its lowest; tower 2 carries 7333 N, 1.44
    // times its 5093 N before
    expect_pulls(got, {0.0, 0.0, -264.42, -7333.21, -25581.09, -1437.75}, {0.0, 0.0, -264.36, -7333.01, -25589.11, -1437.95}, 3.0e-4, "string 1 released");
    EXPECT_NEAR(-got.t2z / before, 1.44, 0.01);
    EXPECT_EQ(c->nodes_on_ground(0) + c->nodes_on_ground(1), 0) << "clear of the ground";
}
