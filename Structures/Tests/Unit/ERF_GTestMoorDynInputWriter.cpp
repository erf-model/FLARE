// Unit tests of erf_conductors::moordyn_input_text and write_moordyn_input, the MoorDyn-C input of one line.
//
// ASingleSpanLineIsWrittenInMoorDynsFrameWithAirAndExternalKinematics: heights lowered by the surface
//     offset (z_md = z_erf - surface_offset), the air density, gravity, external kinematics, a flat
//     bottom below everything, both dead ends Fixed, the line's properties and segments.
// MoorDynAcceptsTheFileAndTheEndsComeBackInERFsFrame: MoorDyn (the stub or the real library) puts the
//     end nodes on the attachment points, and to_erf_frame brings the positions back.
// ASectionIsWrittenWithItsInsulatorStringsAndFreePoints: Fixed points at every attachment point, a
//     Free point under each tower's string, the spans between them and the strings from the towers
//     down; MoorDyn hangs the conductor from the strings.
// WithGroundContactTheInputNamesTheSeafloorAndItsContact: the SeafloorFile, kbot and the conductor's critical cbot
//     rows with ground_contact only; a missing or blank-holding file name is refused; seafloor_file() swaps the ending.
// TheSeafloorGridIsCoarsenedToItsLimitAndMissingHeightsTakeTheLowest: seafloor_ticks takes the smallest stride that
//     fits, always ending on the last node; fill_missing_heights fills nodes no grid holds with the lowest height.
// TheSeafloorTextIsMoorDynsThreeDFormat: "nx ny", the ticks, then x y z rows with x the slower index and z lowered by
//     the surface offset; too few ticks or heights are refused.

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "ERF_GTestTempDir.H"
#include "ERF_GTestThrowOnAbort.H"
#include "ERF_ConductorInputs.H"
#include "ERF_MoorDynInputWriter.H"
#include "ERF_MoorDynSystem.H"

using erf_conductors::ConductorInputs;
using erf_conductors::LineInputs;

namespace {

// a scratch root drawn once per test process (ERF_GTestTempDir.H): the fixed names below it are this
// process's alone, so ctest -j and the shuffled rerun never share them; removed when the process exits
const std::filesystem::path& gtest_scratch_root ()
{
    struct Root {
        std::filesystem::path p;
        ~Root () { std::error_code ec; std::filesystem::remove_all(p, ec); }
    };
    static const Root root{[] {
        const std::filesystem::path p = erf_gtest_temp_path("erf_gtest_moordyninputwriter");
        std::filesystem::create_directories(p);
        return p;
    }()};
    return root.p;
}

// to_erf_frame returns Reals: positions of a few hundred metres carry a float's spacing, 3e-5 m at 500 m
constexpr double ptol = (std::is_same<amrex::Real, float>::value) ? 1.0e-4 : 1.0e-6;

LineInputs span ()
{
    LineInputs s;
    s.name = "S1";
    s.end_a = {{100.0, 500.0, 30.0}};
    s.end_b = {{400.0, 500.0, 40.0}};
    s.lengths = {301.5}; s.diameter = 0.0281; s.mass_per_length = 1.628; s.axial_stiffness = 3.0e7;
    s.drag_coefficient = 1.1; s.damping_ratio = 0.4; s.segments = 16;
    return s;
}

ConductorInputs settings ()
{
    ConductorInputs in;
    in.air_density = 1.15;
    in.surface_offset = 5000.0;
    in.moordyn_dt = 0.002;
    in.moordyn_cfl = 0.12;
    return in;
}

bool has_row (const std::string& text, const std::string& row) { return text.find(row) != std::string::npos; }

// the number at the start of the option row whose second token is key (a single-precision build writes 3e+06)
double row_numbers_before (const std::string& text, const std::string& key)
{
    std::istringstream in(text);
    for (std::string line; std::getline(in, line);) {
        std::istringstream ls(line);
        std::string v, k;
        if ((ls >> v >> k) && k == key) { return std::stod(v); }
    }
    return -1.0;
}

// the numbers of the row whose first token is name (the line types are written with as many digits
// as a Real holds, so they are compared as numbers)
std::vector<double> row_numbers (const std::string& text, const std::string& name)
{
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream ls(line);
        std::string first;
        if (!(ls >> first) || first != name) { continue; }
        std::vector<double> v;
        double x;
        while (ls >> x) { v.push_back(x); }
        return v;
    }
    return {};
}

void expect_numbers (const std::vector<double>& got, const std::vector<double>& want, const std::string& what)
{
    ASSERT_GE(got.size(), want.size()) << what;
    for (std::size_t i = 0; i < want.size(); ++i) {
        EXPECT_NEAR(got[i], want[i], 1.0e-6 * std::max(1.0, std::abs(want[i]))) << what << " column " << i;
    }
}

} // namespace

TEST(MoorDynInputWriter, ASingleSpanLineIsWrittenInMoorDynsFrameWithAirAndExternalKinematics)
{
    const std::string text = erf_conductors::moordyn_input_text(span(), settings(), 9.81);
    EXPECT_TRUE(has_row(text, "LINE TYPES"));
    EXPECT_TRUE(has_row(text, "POINT PROPERTIES"));
    EXPECT_TRUE(has_row(text, "LINES"));
    EXPECT_TRUE(has_row(text, "OPTIONS"));
    // the line type row: diameter, mass, EA, -damping ratio (MoorDyn's "-zeta" form), EI 0, Cd
    expect_numbers(row_numbers(text, "S1"), {0.0281, 1.628, 3.0e7, -0.4, 0.0, 1.1}, "line type");
    // the attachments, fixed, lowered by the surface offset
    EXPECT_TRUE(has_row(text, "1     Fixed     100   500   -4970")) << text;
    EXPECT_TRUE(has_row(text, "2     Fixed     400   500   -4960")) << text;
    // the line: type, attachments 1 and 2, length, segments
    EXPECT_TRUE(has_row(text, "1     S1      1        2         301.5   16")) << text;
    // the options
    EXPECT_TRUE(has_row(text, "0.002   dtM")) << text;
    EXPECT_TRUE(has_row(text, "0.12   CFL")) << text;
    EXPECT_TRUE(has_row(text, "9.81   g")) << text;
    EXPECT_TRUE(has_row(text, "1.15   WtrDnsty")) << text;
    EXPECT_TRUE(has_row(text, "10000   WtrDpth")) << text;
    EXPECT_TRUE(has_row(text, "1             WaveKin")) << text;
    EXPECT_TRUE(has_row(text, "0             ICgenDynamic")) << text;
    // no MoorDyn step when none is set: MoorDyn chooses from its CFL
    ConductorInputs no_dt = settings();
    no_dt.moordyn_dt = 0.0;
    EXPECT_FALSE(has_row(erf_conductors::moordyn_input_text(span(), no_dt, 9.81), "dtM"));
}

TEST(MoorDynInputWriter, MoorDynAcceptsTheFileAndTheEndsComeBackInERFsFrame)
{
    const auto dir = gtest_scratch_root() / "erf_gtest_moordyn_writer";
    const std::string fname = (dir / "S1.moordyn.txt").string();
    const LineInputs s = span();
    const ConductorInputs in = settings();
    erf_conductors::write_moordyn_input(fname, s, in, 9.81);
    ASSERT_TRUE(std::filesystem::exists(fname));

    std::string err;
    auto sys = erf_moordyn::MoorDynSystem::create(fname, "", MOORDYN_ERR_LEVEL, err);
    ASSERT_TRUE(sys) << err;
    EXPECT_EQ(sys->num_coupled_dof(), 0u);
    err = sys->init({}, {});
    ASSERT_TRUE(err.empty()) << err;
    ASSERT_EQ(sys->num_lines(), 1u);
    const unsigned nn = sys->line_num_nodes(1);
    EXPECT_EQ(nn, 17u);
    const auto a = erf_conductors::to_erf_frame(sys->line_node_position(1, 0), in.surface_offset);
    const auto b = erf_conductors::to_erf_frame(sys->line_node_position(1, nn - 1), in.surface_offset);
    for (int d = 0; d < 3; ++d) {
        EXPECT_NEAR(a[static_cast<std::size_t>(d)], s.end_a[static_cast<std::size_t>(d)], ptol) << "end a, dir " << d;
        EXPECT_NEAR(b[static_cast<std::size_t>(d)], s.end_b[static_cast<std::size_t>(d)], ptol) << "end b, dir " << d;
    }
    ASSERT_GT(sys->external_kinematics_init(err), 0u) << err;
}

TEST(MoorDynInputWriter, ASectionIsWrittenWithItsInsulatorStringsAndFreePoints)
{
    LineInputs s = span();
    s.end_b = {{1000.0, 500.0, 30.0}};
    s.towers = {{{400.0, 500.0, 30.0}}, {{700.0, 500.0, 30.0}}};
    s.lengths = {301.5, 301.5, 301.5};
    s.end_b[2] = 30.0;
    s.insulator_length = 2.5;
    s.insulator_mass = 60.0;
    s.segments = 20;
    ConductorInputs in = settings();
    const std::string text = erf_conductors::moordyn_input_text(s, in, 9.81);
    // the string's line type: disc diameter, mass per length, its fixed stiffness, the damping, its drag
    expect_numbers(row_numbers(text, "S1_insulator"), {0.254, 24.0, 1.0e7, -0.4, 0.0, 1.0}, "string line type");
    // fixed points at every attachment, free points 2.5 m under the towers
    EXPECT_TRUE(has_row(text, "1     Fixed     100   500   -4970")) << text;
    EXPECT_TRUE(has_row(text, "2     Fixed     400   500   -4970")) << text;
    EXPECT_TRUE(has_row(text, "4     Fixed     1000   500   -4970")) << text;
    EXPECT_TRUE(has_row(text, "5     Free      400   500   -4972.5")) << text;
    EXPECT_TRUE(has_row(text, "6     Free      700   500   -4972.5")) << text;
    // the spans hang from the free points; the strings run from the towers down to them
    EXPECT_TRUE(has_row(text, "1     S1      1        5         301.5   20")) << text;
    EXPECT_TRUE(has_row(text, "2     S1      5        6         301.5   20")) << text;
    EXPECT_TRUE(has_row(text, "3     S1      6        4         301.5   20")) << text;
    EXPECT_TRUE(has_row(text, "4     S1_insulator      2        5         2.5   2")) << text;
    EXPECT_TRUE(has_row(text, "5     S1_insulator      3        6         2.5   2")) << text;

    // MoorDyn hangs the line from the strings: the free points sit 2.5 m under the towers in still air
    in.moordyn_dt = 0.0;
    const auto dir = gtest_scratch_root() / "erf_gtest_moordyn_writer";
    std::filesystem::create_directories(dir);
    const std::string file = (dir / "section.txt").string();
    erf_conductors::write_moordyn_input(file, s, in, 9.81);
    std::string err;
    auto sys = erf_moordyn::MoorDynSystem::create(file, "", 3, err);
    ASSERT_TRUE(sys) << err;
    ASSERT_TRUE(sys->init({}, {}).empty());
    EXPECT_EQ(sys->num_lines(), 5u);
    for (unsigned p : {5u, 6u}) {
        const auto pos = erf_conductors::to_erf_frame(sys->point_position(p), in.surface_offset);
        EXPECT_NEAR(pos[2], 27.5, 0.05) << "free point " << p;
        EXPECT_NEAR(pos[1], 500.0, ptol) << "free point " << p;
    }
}

TEST(MoorDynInputWriter, WithGroundContactTheInputNamesTheSeafloorAndItsContact)
{
    ConductorInputs in = settings();
    EXPECT_FALSE(has_row(erf_conductors::moordyn_input_text(span(), in, 9.81, "S1.seafloor.txt"), "SeafloorFile"))
        << "no seafloor without ground_contact";
    in.ground_contact = true;
    const std::string text = erf_conductors::moordyn_input_text(span(), in, 9.81, "out/S1.seafloor.txt");
    EXPECT_TRUE(has_row(text, "out/S1.seafloor.txt   SeafloorFile")) << text;
    EXPECT_NEAR(row_numbers_before(text, "kbot"), 3.0e6, 1.0) << text;
    EXPECT_TRUE(has_row(text, "10000   WtrDpth")) << "the flat bottom stays, unused" << text;
    // cbot critical for the conductor: 2 sqrt(3e6 x 1.628 / 0.0281)
    const double cbot = row_numbers_before(text, "cbot");
    EXPECT_NEAR(cbot, 2.0 * std::sqrt(3.0e6 * 1.628 / 0.0281), 1.0e-6 * cbot);
    EXPECT_NEAR(static_cast<double>(ConductorInputs::ground_frequency(1.628, 0.0281)), std::sqrt(3.0e6 * 0.0281 / 1.628), 1.0e-3);
    std::string msg = erf_gtest::abort_message([&] { erf_conductors::moordyn_input_text(span(), in, 9.81); });
    EXPECT_NE(msg.find("conductors.S1: with ground_contact the seafloor file name must be given"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { erf_conductors::moordyn_input_text(span(), in, 9.81, "my run/S1.seafloor.txt"); });
    EXPECT_NE(msg.find("hold no blanks"), std::string::npos) << msg;
    EXPECT_EQ(erf_conductors::seafloor_file("conductors/S1.moordyn.txt"), "conductors/S1.seafloor.txt");
    EXPECT_EQ(erf_conductors::seafloor_file("S1.txt"), "S1.txt.seafloor.txt");
}

TEST(MoorDynInputWriter, TheSeafloorTextIsMoorDynsThreeDFormat)
{
    const std::vector<double> x{0.0, 100.0, 200.0}, y{50.0, 150.0};
    // z = 0.05 x + 0.02 y + 3 at (x[i], y[j]), x the slower index
    std::vector<double> z;
    for (const double xi : x) { for (const double yj : y) { z.push_back(0.05 * xi + 0.02 * yj + 3.0); } }
    const std::string text = erf_conductors::seafloor_text(x, y, z, 5000.0);
    EXPECT_EQ(text, "3 2\n0 100 200\n50 150\n"
                    "0 50 -4996\n0 150 -4994\n100 50 -4991\n100 150 -4989\n200 50 -4986\n200 150 -4984\n");
    std::string msg = erf_gtest::abort_message([&] { erf_conductors::seafloor_text({0.0}, y, {1.0, 2.0}, 0.0); });
    EXPECT_NE(msg.find("at least 2 x 2 ticks"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { erf_conductors::seafloor_text(x, y, {1.0}, 0.0); });
    EXPECT_NE(msg.find("one height per grid point"), std::string::npos) << msg;
}

TEST(MoorDynInputWriter, TheSeafloorGridIsCoarsenedToItsLimitAndMissingHeightsTakeTheLowest)
{
    std::vector<int> ii, jj;
    int stride = 0;
    // 11 x 9 nodes fit in 99 points
    erf_conductors::seafloor_ticks({{0, 1}}, {{10, 9}}, 99, ii, jj, stride);
    EXPECT_EQ(stride, 1);
    EXPECT_EQ(ii.size(), 11u);
    EXPECT_EQ(jj.size(), 9u);
    // at most 30 points: every 3rd node, (4 + 1) x (3 + 1) = 20 (stride 2 gives 6 x 5 = 30, which fits)
    erf_conductors::seafloor_ticks({{0, 1}}, {{10, 9}}, 30, ii, jj, stride);
    EXPECT_EQ(stride, 2);
    EXPECT_EQ(ii, (std::vector<int>{0, 2, 4, 6, 8, 10}));
    EXPECT_EQ(jj, (std::vector<int>{1, 3, 5, 7, 9}));
    erf_conductors::seafloor_ticks({{0, 1}}, {{10, 9}}, 29, ii, jj, stride);
    EXPECT_EQ(stride, 3);
    EXPECT_EQ(ii, (std::vector<int>{0, 3, 6, 9, 10})) << "the last node always included";
    EXPECT_EQ(jj, (std::vector<int>{1, 4, 7, 9}));
    EXPECT_LE(static_cast<long>(ii.size() * jj.size()), 29);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> h{5.0, nan, 2.0, nan};
    double lowest = 0.0;
    EXPECT_EQ(erf_conductors::fill_missing_heights(h, lowest), 2u);
    EXPECT_EQ(lowest, 2.0);
    EXPECT_EQ(h, (std::vector<double>{5.0, 2.0, 2.0, 2.0}));
    std::vector<double> none{nan, nan};
    EXPECT_EQ(erf_conductors::fill_missing_heights(none, lowest), 2u);
}
