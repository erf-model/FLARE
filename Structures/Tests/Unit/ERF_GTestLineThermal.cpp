// Unit tests of the conductors' temperature under a fire (ERF_LineThermal.H).
//
// The pinned values come from Tests/test_files/Conductors_FireHeating/check_fire_heating.py, which writes
// the IEEE 738-2012 terms out independently of the C++ (Python, no shared code).
//
// TheIEEE738TermsMatchTheIndependentCalculation: the air's film properties, the wind angle factor, forced
//     convection (the larger of eqs 3a and 3b), natural convection and radiation at D = 28.1 mm,
//     V = 0.61 m/s across, T = 100 C, T_a = 40 C, to 1e-12.
// HotAirHeatsTheConductor: with the air hotter than the conductor the convection is negative (heat flows in).
// TheSteadyAndImplicitTemperaturesMatchTheIndependentCalculation: the steady temperature with 1000 A, with
//     the sun and with a fire, and backward-Euler steps of 0.5 s and 10 s, to 1e-9 K; with no heat source the
//     steady temperature is the air's.
// APointSourceUnderALongConductorGivesTheLineIntegral: summed over the nodes of a straight conductor at height H
//     over one source, the irradiance times the node spacing tends to chi Q / (4 pi) 2a / (H sqrt(a^2 + H^2)).
// AConductorAlongTheRayReceivesNothing: sin psi = 0; the cutoff leaves far sources out; the floor on R^2 holds.
// TheSpanHelpersAreTrapezoidAndCentralDifferences.
// FireSourcesAreTheBurningCellsInOneOrder: gather_fire_sources keeps the cells with q'' > 0 as Q = q'' dA at their
//     centres, within the reach of the nodes' box, sorted by y then x.
// APeriodicSourceIsSeenAcrossTheSeam: the nearest image in a periodic direction, for the distance and the cutoff.
// TheConductorRadiatesToItsSurroundingsNotToThePlume: heat_gain and a step with plume air at 500 C and surroundings
//     at 26.7 C against the independent calculation (radiating to the plume would add 873 W/m).
// EveryKeyIsRangeChecked: each out-of-range value is refused, naming its key.
// TheSpansLengthFollowsItsMeanTemperatureAndSurvivesARestart: on the stub, a span heated by a source under its
//     middle gets L = L_ref (1 + alpha_T (Tbar - T_ref)) set in MoorDyn; the checkpoint restores the temperatures
//     to the bit and sets the length again; a fire without radiative_fraction is refused.
// TheWindsAngleToTheSpanAndTheReferenceTemperatureAreUsed: a wind along and across the span gives the steady
//     temperatures of phi = 0 and 90 degrees; with reference_temperature the span starts at L (1 + alpha_T (Tbar - T_ref));
//     a line named thermal is refused.
// AStepInPlumeAirRadiatesToTheFirstStepsAir: after a first step in 30 C air, a step in 500 C air convects with 500 C and
//     radiates to 30 C.

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <ERF_StructureModel.H>

#include "ERF_ConductorInputs.H"
#include "ERF_ConductorLine.H"
#include "ERF_GTestTempDir.H"
#include "ERF_GTestThrowOnAbort.H"
#include "ERF_LineThermal.H"

using namespace erf_conductors;

namespace {

constexpr double pi = 3.14159265358979323846;
constexpr double D = 0.0281;

LineThermalProps drake ()
{
    LineThermalProps p;
    p.diameter = D; p.absorptivity = 0.5; p.emissivity = 0.5; p.heat_capacity = 1310.0; p.thermal_expansion = 1.93e-5;
    p.current = 1000.0; p.resistance = {{7.283e-5, 8.688e-5}}; p.resistance_temperatures = {{25.0, 75.0}};
    return p;
}

// three burning cells, (12, 3) and (2, 3) at 1e5 W/m^2 and (5, 14) at 2e5 W/m^2, and one negative value (not
// burning). A free function: nvcc does not take a device lambda inside a test body (a private member function).
void fill_three_burning_cells (amrex::MultiFab& q)
{
    q.setVal(0.0);
    for (amrex::MFIter mfi(q); mfi.isValid(); ++mfi) {
        auto a = q.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if ((i == 12 && j == 3) || (i == 2 && j == 3)) { a(i, j, k) = 1.0e5; }
            if (i == 5 && j == 14) { a(i, j, k) = 2.0e5; }
            if (i == 7 && j == 7) { a(i, j, k) = -1.0; }
        });
    }
}

const std::filesystem::path& scratch ()
{
    struct Root {
        std::filesystem::path p;
        ~Root () { std::error_code ec; std::filesystem::remove_all(p, ec); }
    };
    static const Root root{[] {
        const std::filesystem::path p = erf_gtest_temp_path("erf_gtest_linethermal");
        std::filesystem::create_directories(p);
        return p;
    }()};
    return root.p;
}

} // namespace

TEST(LineThermal, TheIEEE738TermsMatchTheIndependentCalculation)
{
    const AirFilm a = ieee738_air(70.0, 0.0);
    EXPECT_NEAR(a.mu, 2.042759261819309e-05, 1e-12 * 2.0e-5);
    EXPECT_NEAR(a.rho, 1.0287214575543002, 1e-12);
    EXPECT_NEAR(a.k, 0.0294523057, 1e-14);
    EXPECT_NEAR(ieee738_air(70.0, 1000.0).rho, 0.9124663855517542, 1e-12);
    EXPECT_NEAR(ieee738_wind_angle_factor(0.0), 0.388, 1e-12);
    EXPECT_NEAR(ieee738_wind_angle_factor(0.25 * pi), 0.8548932188134524, 1e-12);
    EXPECT_NEAR(ieee738_wind_angle_factor(0.5 * pi), 1.0, 1e-12);
    // forced: eq 3a (82.02 W/m) beats eq 3b (76.97 W/m) at Re = 863; natural 42.37 W/m
    EXPECT_NEAR(ieee738_forced_convection(D, 0.61, 0.5 * pi, 100.0, 40.0, 0.0), 82.02372268417102, 1e-10);
    EXPECT_NEAR(ieee738_natural_convection(D, 100.0, 40.0, 0.0), 42.370660576264946, 1e-10);
    EXPECT_NEAR(ieee738_convection(D, 0.61, 0.5 * pi, 100.0, 40.0, 0.0), 82.02372268417102, 1e-10);
    EXPECT_NEAR(ieee738_radiation(D, 0.5, 100.0, 40.0), 24.406167258312006, 1e-10);
    EXPECT_EQ(ieee738_convection(D, 0.61, 0.5 * pi, 40.0, 40.0, 0.0), 0.0);
}

TEST(LineThermal, HotAirHeatsTheConductor)
{
    EXPECT_NEAR(ieee738_convection(D, 2.0, pi / 3.0, 30.0, 200.0, 0.0), -398.6396136296722, 1e-9);
    EXPECT_LT(ieee738_natural_convection(D, 30.0, 200.0, 0.0), 0.0);
    EXPECT_LT(ieee738_radiation(D, 0.5, 30.0, 200.0), 0.0);
}

TEST(LineThermal, TheSteadyAndImplicitTemperaturesMatchTheIndependentCalculation)
{
    LineThermalProps p = drake();
    ThermalSettings s;
    EXPECT_NEAR(steady_temperature(p, s, 40.0, 0.61, 0.5 * pi, 0.0), 92.10441064276685, 1e-9);
    EXPECT_NEAR(steady_temperature(p, s, 40.0, 0.61, 0.5 * pi, 2.0e4), 236.69028687941528, 1e-9);
    EXPECT_NEAR(implicit_step(p, s, 60.0, 0.5, 40.0, 40.0, 0.61, 0.5 * pi, 2.0e4), 60.12570335109734, 1e-9);
    EXPECT_NEAR(implicit_step(p, s, 60.0, 10.0, 40.0, 40.0, 0.61, 0.5 * pi, 2.0e4), 62.48773151328221, 1e-9);
    ThermalSettings sun;
    sun.solar_flux = 1000.0;
    EXPECT_NEAR(steady_temperature(p, sun, 40.0, 0.61, 0.5 * pi, 0.0), 100.95064979068778, 1e-9);
    // nothing heats it: it stays at the air's temperature
    p.current = 0.0;
    EXPECT_NEAR(steady_temperature(p, s, 25.0, 0.0, 0.5 * pi, 0.0), 25.0, 1e-9);
    // a long step tends to the steady state: m c_p (T - T_old) / dt = 0.023 W/m short of it at dt = 1e7 s
    // (0.0096 K below), 2.3e-5 W/m at 1e10 s
    p = drake();
    EXPECT_NEAR(implicit_step(p, s, 60.0, 1.0e10, 40.0, 40.0, 0.61, 0.5 * pi, 2.0e4), 236.69028687941528, 1e-4);
}

TEST(LineThermal, APointSourceUnderALongConductorGivesTheLineIntegral)
{
    const double H = 15.0, a = 150.0, Q = 1.0e6, chi = 0.3;
    const int n = 3001;
    std::vector<double> x(3 * n), t(3 * n, 0.0);
    for (int m = 0; m < n; ++m) {
        x[3*m] = -a + 2.0 * a * m / (n - 1); x[3*m+1] = 0.0; x[3*m+2] = H;
        t[3*m] = 1.0;
    }
    const std::vector<FireSource> src{FireSource{0.0, 0.0, 0.0, Q}};
    const std::vector<double> E = fire_irradiance(x, t, src, chi, 0.0, 0.0);
    const double h = 2.0 * a / (n - 1);
    std::vector<double> Eh(E.size());
    for (std::size_t m = 0; m < E.size(); ++m) { Eh[m] = E[m]; }
    const double integral = span_mean(Eh) * 2.0 * a;
    const double exact = chi * Q / (4.0 * pi) * 2.0 * a / (H * std::sqrt(a * a + H * H));
    EXPECT_NEAR(integral, exact, 1e-5 * exact) << "node spacing " << h << " m";
    // directly above the source: chi Q / (4 pi H^2)
    EXPECT_NEAR(E[(n - 1) / 2], chi * Q / (4.0 * pi * H * H), 1e-12 * E[(n - 1) / 2]);
}

TEST(LineThermal, AConductorAlongTheRayReceivesNothing)
{
    // the conductor's axis points at the source: it intercepts nothing
    const std::vector<double> x{0.0, 0.0, 10.0}, up{0.0, 0.0, 1.0}, across{1.0, 0.0, 0.0};
    const std::vector<FireSource> src{FireSource{0.0, 0.0, 0.0, 1.0e6}};
    EXPECT_NEAR(fire_irradiance(x, up, src, 0.3, 0.0, 0.0)[0], 0.0, 1e-12);
    const double E = fire_irradiance(x, across, src, 0.3, 0.0, 0.0)[0];
    EXPECT_NEAR(E, 0.3e6 / (4.0 * pi * 100.0), 1e-9);
    // a cutoff below the distance leaves it out, above keeps it
    EXPECT_EQ(fire_irradiance(x, across, src, 0.3, 9.9, 0.0)[0], 0.0);
    EXPECT_EQ(fire_irradiance(x, across, src, 0.3, 10.1, 0.0)[0], E);
    // the floor on R^2
    EXPECT_NEAR(fire_irradiance(x, across, src, 0.3, 0.0, 400.0)[0], 0.3e6 / (4.0 * pi * 400.0), 1e-9);
}

TEST(LineThermal, TheSpanHelpersAreTrapezoidAndCentralDifferences)
{
    EXPECT_DOUBLE_EQ(span_mean({1.0, 2.0, 3.0, 10.0}), (0.5 + 2.0 + 3.0 + 5.0) / 3.0);
    const std::vector<double> x{0.0, 0.0, 0.0, 3.0, 4.0, 0.0, 6.0, 8.0, 0.0};
    const std::vector<double> t = span_tangents(x);
    for (int m = 0; m < 3; ++m) {
        EXPECT_NEAR(t[3*m], 0.6, 1e-15); EXPECT_NEAR(t[3*m+1], 0.8, 1e-15); EXPECT_EQ(t[3*m+2], 0.0);
    }
}

TEST(LineThermal, FireSourcesAreTheBurningCellsInOneOrder)
{
    const amrex::Box dom(amrex::IntVect(0, 0, 0), amrex::IntVect(15, 15, 0));
    const amrex::RealBox rb({0.0, 0.0, 0.0}, {160.0, 160.0, 10.0});
    const amrex::Geometry geom(dom, rb, 0, {0, 0, 0});
    amrex::BoxArray ba(dom);
    ba.maxSize(4);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab q(ba, dm, 1, 0);
    fill_three_burning_cells(q);
    const StructureFire fire{geom, 0, 0.0, &q, nullptr, [] () { return amrex::Real(-1.0); }};
    const auto all = gather_fire_sources(fire, {{0.0, 0.0, 0.0}}, {{1.0, 1.0, 1.0}}, -1.0);
    ASSERT_EQ(all.size(), 3u);
    EXPECT_EQ(all[0].x, 25.0);  EXPECT_EQ(all[0].y, 35.0); EXPECT_EQ(all[0].Q, 1.0e7);
    EXPECT_EQ(all[1].x, 125.0); EXPECT_EQ(all[1].y, 35.0);
    EXPECT_EQ(all[2].x, 55.0);  EXPECT_EQ(all[2].y, 145.0); EXPECT_EQ(all[2].Q, 2.0e7);
    // within 40 m of a node box at (0..10, 0..10): only (25, 35)
    const auto near = gather_fire_sources(fire, {{0.0, 0.0, 0.0}}, {{10.0, 10.0, 0.0}}, 40.0);
    ASSERT_EQ(near.size(), 1u);
    EXPECT_EQ(near[0].x, 25.0);
    // periodic in x: the reach does not cut in x, so (125, 35) comes back too
    const auto wrapped = gather_fire_sources(fire, {{0.0, 0.0, 0.0}}, {{10.0, 10.0, 0.0}}, 40.0, {{160.0, 0.0}});
    ASSERT_EQ(wrapped.size(), 2u);
    EXPECT_EQ(wrapped[1].x, 125.0);
}

TEST(LineThermal, APeriodicSourceIsSeenAcrossTheSeam)
{
    // a node at x = 995 m and a source at x = 5 m of a 1000 m periodic domain are 10 m apart
    const std::vector<double> x{995.0, 0.0, 0.0}, across{0.0, 1.0, 0.0};
    const std::vector<FireSource> src{FireSource{5.0, 0.0, 0.0, 1.0e6}};
    EXPECT_NEAR(fire_irradiance(x, across, src, 0.3, 0.0, 0.0, {{1000.0, 0.0}})[0], 0.3e6 / (4.0 * pi * 100.0), 1e-9);
    EXPECT_NEAR(fire_irradiance(x, across, src, 0.3, 0.0, 0.0)[0], 0.3e6 / (4.0 * pi * 990.0 * 990.0), 1e-12);
    // the cutoff applies to the nearest image
    EXPECT_GT(fire_irradiance(x, across, src, 0.3, 20.0, 0.0, {{1000.0, 0.0}})[0], 0.0);
}

TEST(LineThermal, TheConductorRadiatesToItsSurroundingsNotToThePlume)
{
    // a conductor at 100 C in plume air at 500 C: radiating to the plume would add 873 W/m
    const LineThermalProps p = drake();
    const ThermalSettings s;
    EXPECT_NEAR(heat_gain(p, s, 100.0, 500.0, 500.0, 2.0, 0.5 * pi, 0.0), 1946.7676676931278, 1e-9);
    EXPECT_NEAR(heat_gain(p, s, 100.0, 500.0, 26.7, 2.0, 0.5 * pi, 0.0), 1074.0204485781014, 1e-9);
    EXPECT_NEAR(implicit_step(p, s, 60.0, 1.0, 500.0, 26.7, 2.0, 0.5 * pi, 0.0), 60.89869479016551, 1e-9);
}

TEST(LineThermal, EveryKeyIsRangeChecked)
{
    const ThermalSettings ok;
    EXPECT_EQ(LineThermal::check(ok, "S1", drake()), "");
    auto with = [] (auto f) { ThermalSettings s; f(s); return LineThermal::check(s, "S1", drake()); };
    auto prop = [] (auto f) { LineThermalProps p = drake(); f(p); return LineThermal::check(ThermalSettings{}, "S1", p); };
    EXPECT_NE(with([] (ThermalSettings& s) { s.radiative_fraction_given = true; s.radiative_fraction = 0.0; }).find("radiative_fraction"), std::string::npos);
    EXPECT_NE(with([] (ThermalSettings& s) { s.radiative_fraction_given = true; s.radiative_fraction = 1.5; }).find("radiative_fraction"), std::string::npos);
    EXPECT_NE(with([] (ThermalSettings& s) { s.source_height = -1.0; }).find("source_height"), std::string::npos);
    EXPECT_NE(with([] (ThermalSettings& s) { s.cutoff_radius = -1.0; }).find("cutoff_radius"), std::string::npos);
    EXPECT_NE(with([] (ThermalSettings& s) { s.elevation = 1.0e4; }).find("elevation"), std::string::npos);
    EXPECT_NE(with([] (ThermalSettings& s) { s.solar_flux = -1.0; }).find("solar_flux"), std::string::npos);
    EXPECT_NE(prop([] (LineThermalProps& p) { p.absorptivity = 1.1; }).find("conductors.S1.absorptivity"), std::string::npos);
    EXPECT_NE(prop([] (LineThermalProps& p) { p.emissivity = -0.1; }).find("conductors.S1.emissivity"), std::string::npos);
    EXPECT_NE(prop([] (LineThermalProps& p) { p.heat_capacity = 0.0; }).find("conductors.S1.heat_capacity"), std::string::npos);
    EXPECT_NE(prop([] (LineThermalProps& p) { p.thermal_expansion = 1.0e-3; }).find("conductors.S1.thermal_expansion"), std::string::npos);
    EXPECT_NE(prop([] (LineThermalProps& p) { p.current = -1.0; }).find("conductors.S1.current"), std::string::npos);
    EXPECT_NE(prop([] (LineThermalProps& p) { p.resistance = {{0.0, 1.0e-4}}; }).find("conductors.S1.resistance"), std::string::npos);
    EXPECT_NE(prop([] (LineThermalProps& p) { p.resistance_temperatures = {{75.0, 25.0}}; }).find("resistance_temperatures"), std::string::npos);
    // resistance_temperatures is checked with no current too
    EXPECT_NE(prop([] (LineThermalProps& p) { p.current = 0.0; p.resistance_temperatures = {{75.0, 25.0}}; }).find("resistance_temperatures"), std::string::npos);
    // a resistance line that turns negative inside the solver's range
    EXPECT_NE(prop([] (LineThermalProps& p) { p.resistance = {{1.0e-5, 9.0e-5}}; }).find("must stay positive"), std::string::npos);
    EXPECT_NE(with([] (ThermalSettings& s) { s.reference_temperature_given = true; s.reference_temperature = 400.0; }).find("reference_temperature"), std::string::npos);
    // no current: the resistance is not needed
    EXPECT_EQ(prop([] (LineThermalProps& p) { p.current = 0.0; p.resistance = {{0.0, 0.0}}; }), "");
    EXPECT_NE(prop([] (LineThermalProps& p) { p.heat_capacity = std::numeric_limits<double>::quiet_NaN(); }).find("heat_capacity"), std::string::npos);
}

TEST(LineThermal, TheSpansLengthFollowsItsMeanTemperatureAndSurvivesARestart)
{
    ConductorInputs in;
    in.air_density = 1.2;
    in.diagnostics_dir = (scratch() / "span").string();
    std::filesystem::create_directories(in.diagnostics_dir);
    LineInputs s;
    s.name = "S1";
    s.end_a = {{100.0, 500.0, 30.0}};
    s.end_b = {{400.0, 500.0, 30.0}};
    s.lengths = {301.5}; s.diameter = D; s.mass_per_length = 1.628; s.axial_stiffness = 3.0e7;
    s.output_root = in.diagnostics_dir + "/S1";
    std::vector<std::unique_ptr<ConductorLine>> lines;
    lines.push_back(std::make_unique<ConductorLine>(s, in, 9.81, in.diagnostics_dir + "/S1.moordyn.txt"));
    const unsigned nn = lines[0]->num_nodes();
    const std::vector<std::vector<amrex::Real>> wind{std::vector<amrex::Real>(3 * lines[0]->num_kinematics_points(), 0.0)};
    const std::vector<std::vector<double>> air{std::vector<double>(nn, 30.0)};
    // a source of 1 GW under the middle: about 5 K there over the four steps, a tenth of that on the span's mean
    const std::vector<FireSource> fire{FireSource{250.0, 500.0, 2.0, 1.0e9}};

    ThermalSettings ts;
    {
        LineThermal no_fraction(ts, {"S1"}, {drake()}, in.diagnostics_dir);
        auto lines_copy = std::move(lines);
        const std::string msg = erf_gtest::abort_message([&] { no_fraction.advance(0.0, 1.0, lines_copy, wind, air, fire, 0.0, {{0.0, 0.0}}, true); });
        EXPECT_NE(msg.find("radiative_fraction is needed"), std::string::npos) << msg;
        lines = std::move(lines_copy);
    }
    ts.radiative_fraction = 0.3;
    ts.radiative_fraction_given = true;
    LineThermal th(ts, {"S1"}, {drake()}, in.diagnostics_dir);
    const double L0 = lines[0]->span_unstretched_length(0);
    for (int n = 0; n < 4; ++n) { th.advance(n * 1.0, 1.0, lines, wind, air, fire, 0.0, {{0.0, 0.0}}, true); }
    EXPECT_EQ(th.span_reference_length(0, 0), L0);
    EXPECT_NEAR(th.span_reference_temperature(0, 0), steady_temperature(drake(), ts, 30.0, 0.0, 0.5 * pi, 0.0), 1e-9);
    const double Tbar = th.span_temperature(0, 0);
    EXPECT_GT(Tbar, th.span_reference_temperature(0, 0) + 0.1) << "the fire heats the span";
    EXPECT_GT(th.temperature(0, nn / 2), th.temperature(0, 1)) << "the middle, over the source, is hottest";
    EXPECT_DOUBLE_EQ(lines[0]->span_unstretched_length(0), L0 * (1.0 + 1.93e-5 * (Tbar - th.span_reference_temperature(0, 0))));

    // the checkpoint: the temperatures to the bit, and the length set again
    const std::string chk = (scratch() / "chk").string();
    std::filesystem::create_directories(chk);
    th.write_checkpoint(chk);
    lines[0]->set_span_unstretched_length(0, L0);
    LineThermal back(ts, {"S1"}, {drake()}, in.diagnostics_dir);
    back.restore(chk, lines);
    ASSERT_TRUE(back.started());
    for (unsigned m = 0; m < nn; ++m) {
        EXPECT_EQ(back.temperature(0, m), th.temperature(0, m)) << "node " << m;
        EXPECT_EQ(back.surroundings_temperature(0, m), 30.0) << "node " << m;
    }
    EXPECT_EQ(back.span_temperature(0, 0), Tbar);
    EXPECT_EQ(lines[0]->span_unstretched_length(0), L0 * (1.0 + 1.93e-5 * (Tbar - th.span_reference_temperature(0, 0))));
    // a checkpoint of other lines is refused
    LineThermal other(ts, {"S2"}, {drake()}, in.diagnostics_dir);
    EXPECT_NE(erf_gtest::abort_message([&] { other.restore(chk, lines); }).find("does not match line S2"), std::string::npos);
}

TEST(LineThermal, TheWindsAngleToTheSpanAndTheReferenceTemperatureAreUsed)
{
    ConductorInputs in;
    in.air_density = 1.2;
    in.diagnostics_dir = (scratch() / "wind").string();
    std::filesystem::create_directories(in.diagnostics_dir);
    LineInputs s;
    s.name = "W1";
    s.end_a = {{100.0, 500.0, 30.0}};
    s.end_b = {{400.0, 500.0, 30.0}};
    s.lengths = {301.5}; s.diameter = D; s.mass_per_length = 1.628; s.axial_stiffness = 3.0e7;
    s.output_root = in.diagnostics_dir + "/W1";
    std::vector<std::unique_ptr<ConductorLine>> lines;
    lines.push_back(std::make_unique<ConductorLine>(s, in, 9.81, in.diagnostics_dir + "/W1.moordyn.txt"));
    const unsigned nn = lines[0]->num_nodes(), nk = lines[0]->num_kinematics_points();
    const std::vector<std::vector<double>> air{std::vector<double>(nn, 30.0)};
    // the span runs along x: a 3 m/s wind along x is at 0 to its middle, along y at 90 degrees
    // (pins from check_fire_heating.py: 83.63022336977076 C along, 51.69250070242272 C across)
    const double along = 83.63022336977076, across = 51.69250070242272;
    for (int d = 0; d < 2; ++d) {
        std::vector<std::vector<amrex::Real>> wind{std::vector<amrex::Real>(3 * nk, 0.0)};
        for (unsigned m = 0; m < nk; ++m) { wind[0][3*m+d] = 3.0; }
        ThermalSettings ts;
        LineThermal th(ts, {"W1"}, {drake()}, in.diagnostics_dir);
        th.advance(0.0, 1.0e-6, lines, wind, air, {}, 0.0, {{0.0, 0.0}}, false);
        const double want = (d == 0) ? along : across;
        // a step of 1e-6 s from the steady state stays there
        EXPECT_NEAR(th.temperature(0, nn / 2), want, 1e-6) << (d == 0 ? "along" : "across");
    }
    // reference_temperature: the input length belongs to 20 C, so the span starts lengthened by its steady
    // temperature's excess over it
    ThermalSettings ts;
    ts.reference_temperature = 20.0;
    ts.reference_temperature_given = true;
    LineThermal th(ts, {"W1"}, {drake()}, in.diagnostics_dir);
    const double L0 = lines[0]->span_unstretched_length(0);
    std::vector<std::vector<amrex::Real>> calm{std::vector<amrex::Real>(3 * nk, 0.0)};
    th.advance(0.0, 1.0, lines, calm, air, {}, 0.0, {{0.0, 0.0}}, false);
    EXPECT_EQ(th.span_reference_temperature(0, 0), 20.0);
    EXPECT_DOUBLE_EQ(lines[0]->span_unstretched_length(0), L0 * (1.0 + 1.93e-5 * (th.span_temperature(0, 0) - 20.0)));
    EXPECT_GT(lines[0]->span_unstretched_length(0), L0 + 0.1) << "about 90 K above 20 C: 0.5 m longer";
    // a line named thermal is refused
    EXPECT_NE(erf_gtest::abort_message([&] { LineThermal bad(ts, {"thermal"}, {drake()}, in.diagnostics_dir); }).find("named 'thermal'"),
              std::string::npos);
}

TEST(LineThermal, AStepInPlumeAirRadiatesToTheFirstStepsAir)
{
    ConductorInputs in;
    in.air_density = 1.2;
    in.diagnostics_dir = (scratch() / "plume").string();
    std::filesystem::create_directories(in.diagnostics_dir);
    LineInputs s;
    s.name = "P1";
    s.end_a = {{100.0, 500.0, 30.0}};
    s.end_b = {{400.0, 500.0, 30.0}};
    s.lengths = {301.5}; s.diameter = D; s.mass_per_length = 1.628; s.axial_stiffness = 3.0e7;
    s.output_root = in.diagnostics_dir + "/P1";
    std::vector<std::unique_ptr<ConductorLine>> lines;
    lines.push_back(std::make_unique<ConductorLine>(s, in, 9.81, in.diagnostics_dir + "/P1.moordyn.txt"));
    const unsigned nn = lines[0]->num_nodes(), nk = lines[0]->num_kinematics_points();
    const std::vector<std::vector<amrex::Real>> calm{std::vector<amrex::Real>(3 * nk, 0.0)};
    ThermalSettings ts;
    LineThermal th(ts, {"P1"}, {drake()}, in.diagnostics_dir);
    th.advance(0.0, 1.0, lines, calm, {std::vector<double>(nn, 30.0)}, {}, 0.0, {{0.0, 0.0}}, false);
    const double T1 = th.temperature(0, nn / 2);
    // the plume arrives: the air at the node is 500 C, the surroundings stay at the first step's 30 C
    th.advance(1.0, 1.0, lines, calm, {std::vector<double>(nn, 500.0)}, {}, 0.0, {{0.0, 0.0}}, false);
    EXPECT_EQ(th.surroundings_temperature(0, nn / 2), 30.0);
    EXPECT_NEAR(th.temperature(0, nn / 2), implicit_step(drake(), ts, T1, 1.0, 500.0, 30.0, 0.0, 0.5 * pi, 0.0), 1e-12);
    EXPECT_GT(std::abs(th.temperature(0, nn / 2) - implicit_step(drake(), ts, T1, 1.0, 500.0, 500.0, 0.0, 0.5 * pi, 0.0)), 0.1)
        << "radiating to the plume air would differ";
}
