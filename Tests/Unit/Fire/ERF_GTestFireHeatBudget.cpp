#include <gtest/gtest.h>
#include <cfenv>
#include <cmath>
#include <AMReX_REAL.H>

#include "ERF_FireHeatFlux.H"
#include "ERF_FireDiagnostics.H"

/**
 * @file ERF_GTestFireHeatBudget.cpp
 * @brief The heat a burning cell hands the atmosphere over a step is the heat
 *        of the fuel the step consumes: Q = h w (1 - e^{-dt/tau}) / dt, so the
 *        sum of Q dt over the steps is h (w_0 - w_N) exactly, where the
 *        instantaneous h w / tau at the start of the step overshot it by
 *        (dt/tau) / (1 - e^{-dt/tau}). Byram's intensity is h w_0 R on a cell
 *        the front has reached and still burns, not h (w_0 - w) R, which was
 *        zero on arrival and largest on burned-out cells; that form is kept as
 *        the heat-release diagnostic.
 */

using namespace amrex;

namespace {
constexpr double REL = (sizeof(Real) == 8) ? 1.0e-10 : 1.0e-4;
}

TEST(FireHeatBudget, StepMeanFluxIntegratesToTheFuelBurned)
{
    const Real h = 1.8e7_rt, tau = 60.0_rt, w0 = 1.2_rt;
    for (Real dt : {6.0_rt, 60.0_rt, 180.0_rt}) {
        Real w = w0, E = 0.0_rt;
        for (int n = 0; n < 20; ++n) {
            E += compute_heat_flux_cell(-1.0_rt, w, h, tau, dt) * dt;
            w  = deplete_fuel_load(-1.0_rt, w, tau, dt);
        }
        EXPECT_NEAR(E, h * (w0 - w), REL * h * w0) << "dt / tau = " << dt / tau;
    }
}

TEST(FireHeatBudget, TheInstantaneousFormOvershoots)
{
    const Real h = 1.8e7_rt, tau = 60.0_rt, w0 = 1.2_rt;
    const Real q_inst = compute_heat_flux_cell(-1.0_rt, w0, h, tau, 0.0_rt);
    const Real q_mean = compute_heat_flux_cell(-1.0_rt, w0, h, tau, tau);
    EXPECT_NEAR(q_inst, w0 * h / tau, REL * w0 * h / tau) << "dt <= 0: the instantaneous power";
    EXPECT_NEAR(q_inst / q_mean, 1.0 / (1.0 - std::exp(-1.0)), 1.0e3 * REL) << "1.58x at dt = tau";
    // the step mean tends to the instantaneous power as dt -> 0
    EXPECT_NEAR(compute_heat_flux_cell(-1.0_rt, w0, h, tau, 1.0e-3_rt), q_inst, 1.0e-4 * q_inst);
}

// A non-positive step gives the instantaneous power without evaluating the
// step mean at a positive exponent: e^{-dt/tau} with dt = -1e6 s overflowed
// (raised FE_OVERFLOW, an abort under the FPE traps) before the select
// discarded it (found by Copilot's review of hgopalan/ERF#501). An optimised
// build drops the unselected branch, so the flag shows the old code's
// overflow only in an unoptimised (Debug) build; the power is checked in both.
TEST(FireHeatBudget, ANegativeStepGivesTheInstantaneousPowerWithoutOverflow)
{
    volatile Real dt = -1.0e6_rt;   // volatile: evaluated at run time, not folded
    volatile Real tau = 2.0_rt;
    const Real h = 1.8e7_rt, w = 0.5_rt;
    std::feclearexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO);
    volatile Real q = compute_heat_flux_cell(-1.0_rt, w, h, tau, dt);   // computed before the flag test
    EXPECT_FALSE(std::fetestexcept(FE_OVERFLOW)) << "the step mean was evaluated at exp(+5e5)";
    EXPECT_FALSE(std::fetestexcept(FE_INVALID | FE_DIVBYZERO));
    EXPECT_NEAR(q, w * h / tau, REL * w * h / tau);
}

// A step below the 1e-30 floor of the divide is still a step: the rate is
// 1/tau, not dt/(1e-30 tau) (0.05 instead of 0.5 at dt = 1e-31 s, tau = 2 s)
TEST(FireHeatBudget, ATinyStepGivesTheInstantaneousPower)
{
    volatile Real dt = 1.0e-31_rt;
    const Real h = 1.8e7_rt, w = 0.5_rt, tau = 2.0_rt;
    EXPECT_NEAR(compute_heat_flux_cell(-1.0_rt, w, h, tau, dt), w * h / tau, REL * w * h / tau);
}

TEST(FireHeatBudget, UnburnedAndExhaustedCellsGiveNothing)
{
    EXPECT_EQ(compute_heat_flux_cell( 0.5_rt, 1.0_rt, 1.8e7_rt, 60.0_rt, 10.0_rt), 0.0_rt);
    EXPECT_EQ(compute_heat_flux_cell(-0.5_rt, 0.0_rt, 1.8e7_rt, 60.0_rt, 10.0_rt), 0.0_rt);
    EXPECT_EQ(compute_heat_flux_cell(-0.5_rt, 1.0_rt, 1.8e7_rt,  0.0_rt, 10.0_rt), 0.0_rt);
    EXPECT_EQ(deplete_fuel_load(0.5_rt, 1.0_rt, 60.0_rt, 10.0_rt), 1.0_rt);
}

TEST(ByramIntensity, AFreshCellCarriesTheFullLoad)
{
    const Real h = 18000.0_rt, w0 = 2.0_rt, R = 0.5_rt;
    EXPECT_NEAR(compute_fireline_intensity_kW_per_m(-1.0_rt, R, w0, w0, h), h * w0 * R, REL * h * w0 * R)
        << "nothing consumed yet: I_B = h w0 R";
    EXPECT_NEAR(compute_fireline_intensity_kW_per_m(-1.0_rt, R, w0, 0.5_rt * w0, h), h * w0 * R, REL * h * w0 * R)
        << "half consumed: still the load the front consumes";
    EXPECT_EQ(compute_fireline_intensity_kW_per_m(-1.0_rt, R, w0, 0.005_rt * w0, h), 0.0_rt) << "burned out";
    EXPECT_EQ(compute_fireline_intensity_kW_per_m( 1.0_rt, R, w0, w0, h), 0.0_rt) << "unburned";
    EXPECT_EQ(compute_fireline_intensity_kW_per_m(-1.0_rt, 0.0_rt, w0, w0, h), 0.0_rt) << "no spread";
    EXPECT_NEAR(compute_flame_length_m(h * w0 * R), 0.0775 * std::pow(h * w0 * R, 0.46), 1.0e3 * REL);
}

TEST(ByramIntensity, HeatReleaseIsTheConsumedLoad)
{
    const Real h = 18000.0_rt, w0 = 2.0_rt, R = 0.5_rt;
    EXPECT_EQ(compute_heat_release_kW_per_m(-1.0_rt, R, w0, w0, h), 0.0_rt) << "nothing consumed yet";
    EXPECT_NEAR(compute_heat_release_kW_per_m(-1.0_rt, R, w0, 0.5_rt * w0, h), 0.5 * h * w0 * R, REL * h * w0 * R);
    EXPECT_NEAR(compute_heat_release_kW_per_m(-1.0_rt, R, w0, 0.0_rt, h), h * w0 * R, REL * h * w0 * R) << "burned out";
    EXPECT_EQ(compute_heat_release_kW_per_m( 1.0_rt, R, w0, 0.0_rt, h), 0.0_rt);
}

// Motivation: the residence-time floor is Anderson's (1969) flame residence time
// as Rothermel (1983, GTR INT-143, appendix C) writes it, t_r = 384 / sigma min
// (sigma in 1/ft). Until November 2026 it was 301 rho_p / sigma s, with no source:
// 2.75 s for fuel model 1 (sigma 3500 1/ft, rho_p 32 lb/ft3) instead of 6.58 s.
// Pins: FM1 sigma 3500 -> 6.583 s; FM10 dead-weighted sigma 636.8 -> 36.18 s;
// the 1 s floor below sigma > 23040.
TEST(FireHeatBudget, ResidenceTimeIsAndersonsFlameResidence)
{
    EXPECT_NEAR(compute_residence_time_s(Real(3500.0)), Real(23040.0 / 3500.0), Real(1.0e-4));
    EXPECT_NEAR(compute_residence_time_s(Real(3500.0)), Real(6.583), Real(1.0e-3));
    EXPECT_NEAR(compute_residence_time_s(Real(636.8)), Real(36.18), Real(1.0e-2));
    EXPECT_NEAR(compute_residence_time_s(Real(30000.0)), Real(1.0), Real(1.0e-6));
    EXPECT_NEAR(compute_residence_time_s(Real(0.0)), Real(1.0), Real(1.0e-6));
}

// Motivation: the residence time takes the characteristic (surface-area-weighted)
// SAV of the bed, as Rothermel (1983, appendix C) does; his Table C-1 gives the
// times for the 13 NFFL models. The load-weighted dead mean used until November
// 2026 gave FM4 24.1 s and FM10 36.1 s, 1.8x and 2.8x the table.
// Pins (Table C-1, min): FM1 0.110, FM2 0.138, FM4 0.221, FM10 0.218.
TEST(FireHeatBudget, ResidenceTimeMatchesRothermelTableC1)
{
    const struct { int fm; double t_min; } table[] = {{1, 0.110}, {2, 0.138}, {4, 0.221}, {10, 0.218}};
    for (const auto& r : table) {
        const FuelModelParams fp = get_fuel_params(r.fm, FUEL_SET_ANDERSON13);
        const Real tau = compute_residence_time_s(characteristic_sav(fp));
        EXPECT_NEAR(tau / Real(60.0), Real(r.t_min), Real(0.0015)) << "fuel model " << r.fm;
    }
}
