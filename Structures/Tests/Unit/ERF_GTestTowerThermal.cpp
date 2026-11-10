// Unit tests of the towers' steel temperature under a fire (ERF_TowerThermal.H).
//
// The pinned values come from an independent Python calculation of the same EN 1993-1-2 balance (the steel's
// enthalpy integrated numerically, bisection for the roots; no shared code), the one check_tower_heating.py of
// Tests/test_files/Conductors_TowerFireHeating writes out.
//
// TheSpecificHeatIsEN1993s: c_a at 20, 300, 600 (the second formula), 735 (the peak of 5000), 900 C, and the 20 C
//     value below 20 C.
// TheEnthalpyIsTheSpecificHeatsIntegral: steel_enthalpy against the independent integral at 300, 700, 800 and
//     1100 C (across the peak), to 1e-9 relative, and against a Simpson sum here; continuous at the formula's breaks.
// TheSectionsAreTheAnglesTheTubesAndTheBoxes: member_section of an equal-leg angle (design data), a circular and a
//     rectangular section, solid and hollow; an arbitrary section without an angle's data is refused.
// TheSteadyAndImplicitTemperaturesMatchTheIndependentCalculation: a 100 x 8 angle in air at 26.85 C, steady with
//     1e4 and 5e4 W/m^2 (and with an emissivity of 0.35 and 10 W/(m^2 K)), steady in hot air (400 C) radiating to
//     20 C, two backward-Euler steps of 0.5 s, and a step of
//     60 s from 700 C across the peak of c_a, to 1e-9 K; each step balances its energy; an hour's cooling of a thin
//     angle from 600 C in air at -40 C (its root 600 K below, near the bracket's floor).
// TheJRCWorkedExampleReachesItsCriticalTemperature: the unprotected beam of the JRC report EUR 26698 EN (2014),
//     section 2.8.3.7, in the ISO 834 fire reaches 667 C in about 17 minutes, to half a minute.
// TheReductionChangeIsTheLargestOfKyAndKE: reduction_change of EN 1993-1-2 Table 3.1, relative to the larger factor.
// TheGaussPointsIntegrateACubic.
// EveryKeyIsRangeChecked: each out-of-range value is refused, naming its key.
// AFrameIsRebuiltWhenItsSteelHasSoftenedAndSurvivesARestart: on a generated tower, a small step leaves the frame;
//     a long step under a fire heats every member near its steady temperature and rebuilds the frame at the members'
//     temperatures; a second such step does not; a short warming below the tolerance leaves the frame but moves the
//     members' checks to the steel's temperatures; the checkpoint restores the temperatures and the rebuild count to
//     the bit; a checkpoint of other towers is refused.

#include <cmath>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "ERF_FrameTower.H"
#include "ERF_GTestTempDir.H"
#include "ERF_GTestThrowOnAbort.H"
#include "ERF_LatticeFrame.H"
#include "ERF_MemberChecks.H"
#include "ERF_TowerThermal.H"

using namespace erf_conductors;

namespace {

constexpr double pi = 3.14159265358979323846;

const std::filesystem::path& scratch ()
{
    struct Root {
        std::filesystem::path p;
        ~Root () { std::error_code ec; std::filesystem::remove_all(p, ec); }
    };
    static const Root root{[] {
        const std::filesystem::path p = erf_gtest_temp_path("erf_gtest_towerthermal");
        std::filesystem::create_directories(p);
        return p;
    }()};
    return root.p;
}

MemberSection angle_100x8 ()
{
    MemberSection s;
    s.area = 0.008 * (2.0 * 0.10 - 0.008);
    s.perimeter = 0.4;
    s.width = ((2.0 + std::sqrt(2.0)) * 0.10 + (2.0 - std::sqrt(2.0)) * 0.008) / pi;
    return s;
}

} // namespace

TEST(TowerThermal, TheSpecificHeatIsEN1993s)
{
    EXPECT_NEAR(steel_specific_heat(20.0), 425.0 + 0.773 * 20.0 - 1.69e-3 * 400.0 + 2.22e-6 * 8000.0, 1e-12);
    EXPECT_NEAR(steel_specific_heat(300.0), 564.74, 1e-9);
    EXPECT_NEAR(steel_specific_heat(600.0), 666.0 + 13002.0 / 138.0, 1e-9);
    EXPECT_NEAR(steel_specific_heat(735.0), 5000.0, 1e-9);
    EXPECT_NEAR(steel_specific_heat(900.0), 650.0, 1e-12);
    EXPECT_EQ(steel_specific_heat(-10.0), steel_specific_heat(20.0));
}

TEST(TowerThermal, TheEnthalpyIsTheSpecificHeatsIntegral)
{
    EXPECT_EQ(steel_enthalpy(20.0), 0.0);
    const double pins[][2] = {{300.0, 142920.31786666662}, {700.0, 419106.0750323175}, {800.0, 561600.8257005839},
                              {1100.0, 762063.8428885898}};
    for (const auto& p : pins) { EXPECT_NEAR(steel_enthalpy(p[0]), p[1], 1e-9 * p[1]) << p[0]; }
    // a Simpson sum here, split at the breaks
    auto simpson = [] (double a, double b) {
        const int n = 20000;
        const double h = (b - a) / n;
        double s = steel_specific_heat(a) + steel_specific_heat(b);
        for (int i = 1; i < n; ++i) { s += (i % 2 ? 4.0 : 2.0) * steel_specific_heat(a + i * h); }
        return s * h / 3.0;
    };
    EXPECT_NEAR(steel_enthalpy(760.0) - steel_enthalpy(650.0), simpson(650.0, 735.0) + simpson(735.0, 760.0), 1e-6);
    // continuous at the breaks; below 20 C the 20 C value
    for (const double b : {600.0, 735.0, 900.0}) { EXPECT_NEAR(steel_enthalpy(b - 1e-9), steel_enthalpy(b), 1e-4) << b; }
    EXPECT_NEAR(steel_enthalpy(0.0), -20.0 * steel_specific_heat(20.0), 1e-9);
}

TEST(TowerThermal, TheSectionsAreTheAnglesTheTubesAndTheBoxes)
{
    erf_towers::FrameInputs in;
    in.joints = {{1, {{0.0, 0.0, 0.0}}}, {2, {{0.0, 0.0, 5.0}}}};
    erf_towers::FrameSection arb; arb.id = 1; arb.shape = erf_towers::SectionShape::Arbitrary; arb.A = 1.0e-3;
    erf_towers::FrameSection tube; tube.id = 2; tube.shape = erf_towers::SectionShape::Circular; tube.D = 0.2; tube.t = 0.01;
    erf_towers::FrameSection rod = tube; rod.id = 3; rod.t = 0.0;
    erf_towers::FrameSection box; box.id = 4; box.shape = erf_towers::SectionShape::Rectangular; box.Sa = 0.3; box.Sb = 0.1; box.t = 0.01;
    in.sections = {arb, tube, rod, box};
    auto member = [] (int id, int sec, erf_towers::SectionShape sh) {
        erf_towers::FrameMember m; m.id = id; m.joint_a = 1; m.joint_b = 2; m.section = sec; m.shape = sh; return m;
    };
    in.members = {member(10, 1, erf_towers::SectionShape::Arbitrary), member(11, 2, erf_towers::SectionShape::Circular),
                  member(12, 3, erf_towers::SectionShape::Circular), member(13, 4, erf_towers::SectionShape::Rectangular),
                  member(14, 1, erf_towers::SectionShape::Arbitrary)};
    std::vector<erf_towers::MemberDesign> d(5);
    d[0].b = 0.10; d[0].t = 0.008;
    MemberSection s;
    std::string err;
    ASSERT_TRUE(member_section(in, &d, 0, s, err)) << err;
    EXPECT_NEAR(s.area, 0.001536, 1e-15);
    EXPECT_NEAR(s.perimeter, 0.4, 1e-15);
    EXPECT_NEAR(s.width, ((2.0 + std::sqrt(2.0)) * 0.1 + (2.0 - std::sqrt(2.0)) * 0.008) / pi, 1e-15);
    ASSERT_TRUE(member_section(in, &d, 1, s, err)) << err;
    EXPECT_NEAR(s.area, 0.25 * pi * (0.04 - 0.0324), 1e-15);
    EXPECT_NEAR(s.perimeter, 0.2 * pi, 1e-15);
    EXPECT_NEAR(s.width, 0.2, 1e-15);
    ASSERT_TRUE(member_section(in, &d, 2, s, err)) << err;
    EXPECT_NEAR(s.area, 0.25 * pi * 0.04, 1e-15);
    ASSERT_TRUE(member_section(in, &d, 3, s, err)) << err;
    EXPECT_NEAR(s.area, 0.03 - 0.28 * 0.08, 1e-15);
    EXPECT_NEAR(s.perimeter, 0.8, 1e-15);
    EXPECT_NEAR(s.width, 0.8 / pi, 1e-15);
    // an arbitrary section without an angle's data has no perimeter
    EXPECT_FALSE(member_section(in, &d, 4, s, err));
    EXPECT_NE(err.find("member 14 has an arbitrary section"), std::string::npos) << err;
    EXPECT_FALSE(member_section(in, nullptr, 0, s, err));
}

TEST(TowerThermal, TheSteadyAndImplicitTemperaturesMatchTheIndependentCalculation)
{
    const MemberSection a = angle_100x8();
    const TowerThermalSettings s;
    EXPECT_NEAR(member_steady_temperature(a, s, 26.85, 26.85, 0.0), 26.85, 1e-9);
    EXPECT_NEAR(member_steady_temperature(a, s, 26.85, 26.85, 1.0e4), 89.39851390916593, 1e-9);
    EXPECT_NEAR(member_steady_temperature(a, s, 26.85, 26.85, 5.0e4), 278.3874128640051, 1e-9);
    EXPECT_NEAR(member_steady_temperature(a, s, 400.0, 20.0, 0.0), 271.7523621444602, 1e-9);
    // the emissivity and the convection coefficient are the settings' (galvanised steel, a light wind)
    TowerThermalSettings g;
    g.emissivity = 0.35;
    g.convection = 10.0;
    EXPECT_NEAR(member_steady_temperature(a, g, 26.85, 26.85, 1.0e4), 100.59833027182998, 1e-9);
    double th = member_implicit_step(a, s, 26.85, 0.5, 26.85, 26.85, 5.0e4);
    EXPECT_NEAR(th, 27.20916514340992, 1e-9);
    th = member_implicit_step(a, s, th, 0.5, 26.85, 26.85, 5.0e4);
    EXPECT_NEAR(th, 27.56773978583844, 1e-9);
    // across the peak of c_a near 735 C in one step
    const double hot = member_implicit_step(a, s, 700.0, 60.0, 26.85, 26.85, 1.0e6);
    EXPECT_NEAR(hot, 882.0192116848236, 1e-8);
    // each step balances its energy
    EXPECT_NEAR(steel_density * a.area * (steel_enthalpy(hot) - steel_enthalpy(700.0)),
                60.0 * member_heat_gain(a, s, hot, 26.85, 26.85, 1.0e6), 1e-6);
    // a long step cooling a thin 30 x 3 angle from 600 C in air at -40 C: the root lies 600 K below, which the
    // bracket's trial points reach without stepping past the -200 C floor
    MemberSection thin;
    thin.area = 0.003 * (2.0 * 0.03 - 0.003);
    thin.perimeter = 0.12;
    thin.width = ((2.0 + std::sqrt(2.0)) * 0.03 + (2.0 - std::sqrt(2.0)) * 0.003) / pi;
    EXPECT_NEAR(member_implicit_step(thin, s, 599.99, 3600.0, -40.0, -40.0, 0.0), -1.1199022431394923, 1e-8);
    // a very long step lands on the steady state, short of it by backward Euler's (tau / dt) (T_ss - T_old), about
    // 400 s / 1e10 s x 250 K = 1e-5 K for this angle (tau 450 s near ambient, 350 s at 250 C)
    EXPECT_NEAR(member_implicit_step(a, s, 26.85, 1.0e10, 26.85, 26.85, 5.0e4), 278.3874128640051, 1e-4);
}

TEST(TowerThermal, TheJRCWorkedExampleReachesItsCriticalTemperature)
{
    // Vassart et al. (2014), Eurocodes: Background & Applications. Structural Fire Design. Worked examples, JRC report
    // EUR 26698 EN, doi:10.2788/85432, section 2.8.3.7: an unprotected steel beam with k_sh A_m/V = 0.706 x 186 =
    // 131.3 1/m in the ISO 834 fire, theta_g = 20 + 345 log10(8 t + 1) (t in minutes), with emissivity 0.7, alpha_c =
    // 25 W/(m^2 K) and steps of 3 s, reaches its critical temperature of 667 C in about 17 minutes with c_a varying
    // with temperature (16 min 30 s with c_a = 600 J/(kg K)). The member here is that section factor: area 0.01 m^2,
    // exchange perimeter 1.313 m; the gas is both its air and its surroundings, at the end of each step.
    MemberSection beam;
    beam.area = 0.01;
    beam.perimeter = 1.313;
    beam.width = 0.0;
    TowerThermalSettings s;
    s.emissivity = 0.7;
    s.convection = 25.0;
    auto gas = [] (double t) { return 20.0 + 345.0 * std::log10(8.0 * t / 60.0 + 1.0); };
    double th = 20.0, t = 0.0;
    while (th < 667.0 && t < 3600.0) {
        const double tg = gas(t + 3.0);
        th = member_implicit_step(beam, s, th, 3.0, tg, tg, 0.0);
        t += 3.0;
    }
    // the report gives the minute: 17 minutes, so within half a minute of it
    EXPECT_NEAR(t, 17.0 * 60.0, 30.0);
}

TEST(TowerThermal, TheReductionChangeIsTheLargestOfKyAndKE)
{
    EXPECT_EQ(reduction_change({20.0, 90.0}, {50.0, 20.0}), 0.0);
    EXPECT_NEAR(reduction_change({150.0, 20.0}, {20.0, 20.0}), 0.05, 1e-12);       // k_E 0.95 against 1
    EXPECT_NEAR(reduction_change({20.0, 450.0}, {20.0, 20.0}), 0.35, 1e-12);       // k_y 0.89, k_E 0.65 against 1
    // relative to the larger factor: k_y 0.78 -> 0.625 is 0.155/0.78, k_E 0.6 -> 0.455 is 0.145/0.6
    EXPECT_NEAR(reduction_change({550.0}, {500.0}), 0.145 / 0.6, 1e-12);
    // where the factors are small a few kelvin move them by more than the tolerance: k_E 0.0675 -> 0.06075 at 900 -> 930 C
    EXPECT_NEAR(reduction_change({930.0}, {900.0}), 0.1, 1e-12);
}

TEST(TowerThermal, TheGaussPointsIntegrateACubic)
{
    const auto p = gauss_points({{1.0, 2.0, 0.0}}, {{1.0, 2.0, 10.0}});
    const auto& w = gauss_weights();
    double I = 0.0, W = 0.0;
    for (std::size_t k = 0; k < 3; ++k) {
        EXPECT_EQ(p[k][0], 1.0);
        EXPECT_EQ(p[k][1], 2.0);
        const double z = p[k][2];
        I += w[k] * (z * z * z - 2.0 * z);
        W += w[k];
    }
    EXPECT_NEAR(W, 1.0, 1e-15);
    EXPECT_NEAR(I, (1.0e4 / 4.0 - 100.0) / 10.0, 1e-10);
    EXPECT_NEAR(p[1][2], 5.0, 1e-15);
}

TEST(TowerThermal, EveryKeyIsRangeChecked)
{
    EXPECT_TRUE(check(TowerThermalSettings{}).empty());
    TowerThermalSettings s;
    s.emissivity = 0.0;
    EXPECT_NE(check(s).find("tower_emissivity"), std::string::npos);
    s = TowerThermalSettings{}; s.emissivity = 1.01;
    EXPECT_NE(check(s).find("tower_emissivity"), std::string::npos);
    s = TowerThermalSettings{}; s.convection = -1.0;
    EXPECT_NE(check(s).find("tower_convection"), std::string::npos);
    s = TowerThermalSettings{}; s.convection = 201.0;
    EXPECT_NE(check(s).find("tower_convection"), std::string::npos);
    s = TowerThermalSettings{}; s.reduction_tolerance = 0.0;
    EXPECT_NE(check(s).find("tower_reduction_tolerance"), std::string::npos);
    s = TowerThermalSettings{}; s.reduction_tolerance = 0.3;
    EXPECT_NE(check(s).find("tower_reduction_tolerance"), std::string::npos);
    s = TowerThermalSettings{}; s.reduction_tolerance = 0.2; s.convection = 0.0; s.emissivity = 1.0;
    EXPECT_TRUE(check(s).empty());
    s.convection = std::nan("");
    EXPECT_NE(check(s).find("tower_convection"), std::string::npos);
    s = TowerThermalSettings{}; s.emissivity = 0.0;
    EXPECT_NE(erf_gtest::abort_message([&] { TowerThermal bad(s, {"t"}, {{angle_100x8()}}, scratch().string()); }).find("tower_emissivity"),
              std::string::npos);
}

TEST(TowerThermal, AFrameIsRebuiltWhenItsSteelHasSoftenedAndSurvivesARestart)
{
    // a generated 6-panel tower at (100, 200, 50), its cross-arm along y
    erf_towers::LatticeSpec spec;
    spec.base_width = 6.0; spec.top_width = 1.5; spec.arm_height = 30.0; spec.arm_length = 12.0; spec.arm_depth = 1.5; spec.peak = 3.0;
    spec.panels = 6; spec.leg_b = 0.15; spec.leg_t = 0.012; spec.brace_b = 0.10; spec.brace_t = 0.008;
    erf_towers::FrameInputs in;
    std::vector<erf_towers::MemberDesign> d;
    std::vector<int> load_joints;
    ASSERT_TRUE(erf_towers::lattice_frame(spec, in, d, &load_joints).empty());
    std::string err;
    std::shared_ptr<const erf_towers::Frame> frame(erf_towers::Frame::create(in, err));
    ASSERT_TRUE(frame) << err;
    auto designs = std::make_shared<const std::vector<erf_towers::MemberDesign>>(d);
    std::vector<std::size_t> links;
    for (const int id : load_joints) { links.push_back(frame->node_of_joint(id)); }
    erf_towers::TowerType type;
    type.name = "generated";
    type.base_width = 6.0; type.top_width = 1.5; type.solidity = 0.2; type.arm_length = 12.0; type.arm_depth = 1.5; type.peak = 3.0;
    type.frame_panels = 6; type.leg_angle = {0.15, 0.012}; type.brace_angle = {0.10, 0.008};
    erf_towers::Tower tw("t1", type, {{100.0, 200.0, 50.0}}, 30.0, {{0.0, 1.0, 0.0}});
    tw.add_attachment({{100.0, 200.0, 80.0}});
    erf_towers::FrameTower ft(tw, frame, 9.80665, designs, "generated", links);
    // a member's point in ERF's axes: the base plus the frame's offset (along = x, across = y here)
    const auto& j = in.joints.front().x;
    const auto pe = ft.erf_point(j);
    EXPECT_NEAR(pe[0], 100.0 + j[0], 1e-12);
    EXPECT_NEAR(pe[1], 200.0 + j[1], 1e-12);
    EXPECT_NEAR(pe[2], 50.0 + j[2], 1e-12);

    const std::size_t nm = in.members.size();
    std::vector<MemberSection> sec(nm);
    for (std::size_t m = 0; m < nm; ++m) { ASSERT_TRUE(member_section(in, designs.get(), m, sec[m], err)) << err; }
    const std::string dir = (scratch() / "diag").string();
    std::filesystem::create_directories(dir);
    TowerThermal tt(TowerThermalSettings{}, {"t1"}, {sec}, dir);
    std::vector<erf_towers::FrameTower*> frames{&ft};
    const std::vector<std::vector<double>> air{std::vector<double>(nm, 26.85)};
    const std::vector<std::vector<double>> cold{std::vector<double>(nm, 0.0)};
    std::vector<std::vector<double>> fire{std::vector<double>(nm, 0.0)};
    // a 5e4 W/m^2 fire on the lower half of the members
    for (std::size_t m = 0; m < nm; ++m) { if (m % 2 == 0) { fire[0][m] = 5.0e4; } }
    // a short step: nothing softens, the frame stays
    tt.advance(0.0, 1.0, frames, air, fire, true);
    EXPECT_EQ(tt.rebuilds(0), 0);
    EXPECT_EQ(ft.temperature(), std::vector<double>(nm, 20.0));
    EXPECT_EQ(tt.surroundings_temperature(0, 0), 26.85);
    // a long step: every heated member near its steady temperature, the frame rebuilt at them
    const double f_cold = static_cast<double>(ft.frequency());
    // (1e9 s: backward Euler's lag, (tau / dt) 250 K, is 1e-4 K)
    tt.advance(1.0, 1.0e9, frames, air, fire, true);
    EXPECT_EQ(tt.rebuilds(0), 1);
    EXPECT_EQ(ft.temperature(), tt.temperatures(0));
    EXPECT_LT(static_cast<double>(ft.frequency()), f_cold);
    for (std::size_t m = 0; m < nm; ++m) {
        EXPECT_NEAR(tt.temperature(0, m), member_steady_temperature(sec[m], tt.settings(), 26.85, 26.85, fire[0][m]), 1e-3) << m;
    }
    // the 100 x 8 bracing: the independent steady temperature
    for (std::size_t m = 0; m < nm; ++m) {
        if (d[m].b == 0.10 && fire[0][m] > 0.0) { EXPECT_NEAR(tt.temperature(0, m), 278.3874128640051, 1e-3); break; }
    }
    // nothing has moved since: no rebuild
    tt.advance(1.0e9, 1.0e9, frames, air, fire, true);
    EXPECT_EQ(tt.rebuilds(0), 1);
    // a little more heat for 10 s: the steel warms by a few kelvin, less than the tolerance, so the frame stays; the
    // members' checks take the steel's new temperatures (their E scales with k_E of the present temperature)
    const auto at_rebuild = ft.member_checks(tw);                 // at the steel's temperatures after the step above
    const std::vector<double> before = tt.temperatures(0);
    std::vector<std::vector<double>> more = fire;
    for (auto& e : more[0]) { e *= 2.0; }
    const std::vector<double> built = ft.temperature();
    tt.advance(2.0e9, 10.0, frames, air, more, true);
    EXPECT_EQ(tt.rebuilds(0), 1);
    EXPECT_EQ(ft.temperature(), built);
    const auto now = ft.member_checks(tw);
    bool warmer = false;
    for (std::size_t m = 0; m < nm; ++m) {
        double ky0 = 0.0, kE0 = 0.0, ky1 = 0.0, kE1 = 0.0;
        erf_towers::steel_reduction(before[m], ky0, kE0);
        erf_towers::steel_reduction(tt.temperature(0, m), ky1, kE1);
        EXPECT_NEAR(now[m].E, at_rebuild[m].E * kE1 / kE0, 1e-9 * at_rebuild[m].E) << m;
        warmer = warmer || (tt.temperature(0, m) > built[m] + 1.0 && kE1 < kE0);
    }
    EXPECT_TRUE(warmer) << "no member warmed enough to change its check";
    tt.advance(2.0e9 + 10.0, 1.0e9, frames, air, fire, true);
    // the checkpoint restores the temperatures and the count to the bit
    const std::string chk = (scratch() / "chk").string();
    std::filesystem::create_directories(chk);
    tt.write_checkpoint(chk);
    TowerThermal back(TowerThermalSettings{}, {"t1"}, {sec}, dir);
    back.restore(chk);
    EXPECT_TRUE(back.started());
    EXPECT_EQ(back.rebuilds(0), 1);
    for (std::size_t m = 0; m < nm; ++m) {
        EXPECT_EQ(back.temperature(0, m), tt.temperature(0, m));
        EXPECT_EQ(back.surroundings_temperature(0, m), 26.85);
    }
    // cooling without the fire brings the steel back to the air
    tt.advance(3.0e9, 1.0e10, frames, air, cold, true);
    EXPECT_EQ(tt.rebuilds(0), 2);
    // (backward Euler's lag of about (tau / dt) 250 K = 1e-5 K)
    for (std::size_t m = 0; m < nm; ++m) { EXPECT_NEAR(tt.temperature(0, m), 26.85, 1e-4); }
    // a checkpoint of another tower is refused
    TowerThermal other(TowerThermalSettings{}, {"t2"}, {sec}, dir);
    EXPECT_NE(erf_gtest::abort_message([&] { other.restore(chk); }).find("does not match tower t2"), std::string::npos);
    // a member past 1200 C stops the run
    std::vector<std::vector<double>> blaze{std::vector<double>(nm, 1.0e7)};
    TowerThermal melt(TowerThermalSettings{}, {"t1"}, {sec}, dir);
    EXPECT_NE(erf_gtest::abort_message([&] { melt.advance(0.0, 1.0e4, frames, air, blaze, false); }).find("keeps no strength"),
              std::string::npos);
}
