#include <gtest/gtest.h>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <cmath>
#include <string>

#include "ERF_FireParams.H"
#include "ERF_AlbiniSpotting.H"
#include "ERF_FireWindAdjust.H"

/**
 * @file ERF_GTestEmberDrift.cpp
 * @brief Firebrands drift on the wind at their height: the reference wind of
 *        wind_ref_ht scaled by the neutral log profile ln(z/z0)/ln(z_ref/z0)
 *        (before October 2026 they drifted on the WAF-reduced midflame wind
 *        throughout a fall of 50-300 m). The landing distance is the forward
 *        Euler sum of that profile over the fall, replicated here; the cap is
 *        the source cell's own fuel on a fuel map; and the launch list is
 *        gathered and sorted, so the landings do not depend on the box
 *        decomposition.
 */

using namespace amrex;

namespace {

constexpr int  NX = 100, NY = 20;
constexpr Real DX = 10.0;
constexpr Real Z_REF = 6.1, Z0 = 0.1, V_T = 0.5, I_B = 2.0e4, U = 2.0;
constexpr int  N_TRAJ = 20;
constexpr double TOLP = (sizeof(Real) == 8) ? 1.0e-12 : 1.0e-5;

/// Launch block: cells 2-4 by 8-11, intensity I_B (above the 100 kW/m threshold, so every cell launches)
bool launches (int i, int j) { return i >= 2 && i <= 4 && j >= 8 && j <= 11; }

struct SpotCase
{
    BoxArray ba; DistributionMapping dm; Geometry geom;
    MultiFab phi, data, wind, intensity, fuel;
    SpotCase (int max_x, int max_y, int per_x = 0)
    {
        Box domain(IntVect(0, 0, 0), IntVect(NX - 1, NY - 1, 0));
        ba = BoxArray(domain);
        ba.maxSize(IntVect(max_x, max_y, 1));
        dm = DistributionMapping(ba);
        geom = Geometry(domain, RealBox(0.0, 0.0, 0.0, NX * DX, NY * DX, 1.0), CoordSys::cartesian, {per_x, 0, 0});
        phi.define(ba, dm, 1, 1); data.define(ba, dm, 4, 0); wind.define(ba, dm, 2, 0);
        intensity.define(ba, dm, 1, 0); fuel.define(ba, dm, 1, 0);
        phi.setVal(1.0_rt); data.setVal(0.0_rt); wind.setVal(U, 0, 1); wind.setVal(0.0_rt, 1, 1);
        intensity.setVal(0.0_rt); fuel.setVal(13.0_rt);
        for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
            auto p = phi.array(mfi);
            auto I = intensity.array(mfi);
            const Box& bx = mfi.validbox();
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j)
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i)
                    if (launches(i, j)) { p(i, j, 0) = -1.0_rt; I(i, j, 0) = I_B; }
        }
    }
    FireParams::SpottingParams params () const
    {
        FireParams::SpottingParams sp;
        sp.enable = true; sp.I_B_min = 100.0; sp.P_base = 1.0; sp.terminal_velocity = V_T;
        sp.n_traj_steps = N_TRAJ; sp.spot_radius = 10.0; sp.P_catch = 1.0; sp.random_seed = 7;
        sp.launch_from = "burned";
        return sp;
    }
    /// phi gathered to one box
    MultiFab gathered () const
    {
        BoxArray one(geom.Domain());
        MultiFab out(one, DistributionMapping(one), 1, 0);
        out.ParallelCopy(phi, 0, 0, 1);
        return out;
    }
};

/// The forward Euler drift of a brand lofted to H_z on flat ground
Real euler_drift (Real H_z)
{
    const Real dz = H_z / N_TRAJ, dt = dz / V_T;
    Real z = H_z, d = 0.0;
    for (int n = 0; n < N_TRAJ; ++n) {
        d += albini_wind_profile_factor(z, Z_REF, Z0) * U * dt;
        z -= dz;
        if (z <= 0.0_rt) { break; }
    }
    return d;
}

/// Range of x of the burned cells beyond the launch block (the landings)
void landing_x_range (const MultiFab& phi_one, Real& xmin, Real& xmax, int& n)
{
    xmin = 1.0e30; xmax = -1.0e30; n = 0;
    for (MFIter mfi(phi_one); mfi.isValid(); ++mfi) {
        auto p = phi_one.const_array(mfi);
        for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i)
                if (p(i, j, 0) < 0.0_rt && !launches(i, j)) {
                    const Real x = (i + 0.5_rt) * DX;
                    xmin = amrex::min(xmin, x); xmax = amrex::max(xmax, x); ++n;
                    EXPECT_GE(j, 6); EXPECT_LE(j, 13);   // no drift across the wind
                }
    }
}

} // namespace

TEST(EmberDrift, ProfileFactor)
{
    EXPECT_NEAR(albini_wind_profile_factor(Z_REF, Z_REF, Z0), 1.0, TOLP);
    EXPECT_NEAR(albini_wind_profile_factor(Z0, Z_REF, Z0), 0.0, TOLP);
    EXPECT_NEAR(albini_wind_profile_factor(0.0, Z_REF, Z0), 0.0, TOLP) << "floored at z0";
    EXPECT_NEAR(albini_wind_profile_factor(61.0, Z_REF, Z0), std::log(610.0) / std::log(61.0), TOLP);
    EXPECT_NEAR(albini_wind_profile_factor(200.0, Z_REF, Z0), 1.85, 0.01) << "about twice the 6.1 m wind at 200 m";
}

/// A roughness below 1e-6 m is accepted by the range checks (z0 > 0) and is
/// used as given, as in the ERF sampler, not floored at 1e-6 m (0.4 % at
/// z0 = 1e-7 m between the drift and the 10 m wind of the same column).
TEST(EmberDrift, SmallRoughnessIsNotFloored)
{
    const Real z0 = 1.0e-7;
    const double tol = (sizeof(Real) == 8) ? 1.0e-12 : 1.0e-5;
    EXPECT_NEAR(albini_wind_profile_factor(10.0, Z_REF, z0), std::log(10.0 / 1.0e-7) / std::log(6.1 / 1.0e-7), tol);
    EXPECT_NEAR(neutral_wind_height_factor(10.0, Z_REF, z0), std::log(10.0 / 1.0e-7) / std::log(6.1 / 1.0e-7), tol);
    EXPECT_NEAR(neutral_wind_height_factor(6.096, 10.0, z0), std::log(6.096 / 1.0e-7) / std::log(10.0 / 1.0e-7), tol);
}

/// INT-309's lofting height for an Anderson model at the case's wind: the
/// 2 m/s reference wind at 6.1 m is 2.24 m/s at 10 m, below the 5 m/s floor of
/// the Table 4 fits, so U10 = 5 m/s.
Real case_height (int anderson)
{
    Real A = 0.0, B = 0.0;
    EXPECT_TRUE(albini_thermal_fit(anderson, A, B));
    return 0.173 * std::sqrt(I_B * A * std::pow(5.0, B));
}

TEST(EmberDrift, LoftingHeightIsAlbinisEquationTwo)
{
    // INT-309 example 1: short grass, 2000 kW/m, U10 = 5 m/s. The worked
    // example uses B = -1.465 and prints E = 30 660 kJ/m, H = 30 m.
    EXPECT_NEAR(albini_thermal_strength(2000.0, 5.0, 162.0, -1.465), 30660.0, 15.0);
    EXPECT_NEAR(albini_lofting_height(2000.0, 5.0, 162.0, -1.465), 30.29, 0.01);
    // Example 2: mature chaparral (Table 4: 83.9, -1.04), 50 000 kW/m, 20 m/s:
    // f(U) = 3.72 s, E = 186 000 kJ/m, H = 75 m.
    Real A = 0.0, B = 0.0;
    ASSERT_TRUE(albini_thermal_fit(4, A, B));
    EXPECT_NEAR(A, 83.9, TOLP); EXPECT_NEAR(B, -1.040, TOLP);
    EXPECT_NEAR(albini_thermal_strength(50000.0, 20.0, A, B) / 50000.0, 3.72, 0.005);
    EXPECT_NEAR(albini_lofting_height(50000.0, 20.0, A, B), 74.6, 0.1);
    // Table 4's short grass (B = -1.456) at example 1
    ASSERT_TRUE(albini_thermal_fit(1, A, B));
    EXPECT_NEAR(albini_lofting_height(2000.0, 5.0, A, B), 30.5, 0.05);
    // below the threshold, nothing is lofted
    EXPECT_EQ(albini_lofting_height(99.0, 5.0, A, B, 100.0), 0.0);
}

TEST(EmberDrift, ThermalFitRangeAndFallback)
{
    Real A = 0.0, B = 0.0;
    // the wind is held to the 5-30 m/s range of the fits
    ASSERT_TRUE(albini_thermal_fit(1, A, B));
    EXPECT_EQ(albini_thermal_strength(1000.0, 0.0, A, B), albini_thermal_strength(1000.0, 5.0, A, B));
    EXPECT_EQ(albini_thermal_strength(1000.0, 45.0, A, B), albini_thermal_strength(1000.0, 30.0, A, B));
    EXPECT_GT(albini_thermal_strength(1000.0, 6.0, A, B), albini_thermal_strength(1000.0, 12.0, A, B));
    // model 8 and codes outside 1-13 have no fit
    EXPECT_FALSE(albini_thermal_fit(8, A, B));
    EXPECT_FALSE(albini_thermal_fit(0, A, B));
    EXPECT_FALSE(albini_thermal_fit(102, A, B));
    // a cell's own Anderson code when Table 4 has it, else the stand-in
    EXPECT_EQ(albini_thermal_model(4, "13", 0), 4);
    EXPECT_EQ(albini_thermal_model(8, "13", 9), 9);
    EXPECT_EQ(albini_thermal_model(8, "13", 0), 0);
    EXPECT_EQ(albini_thermal_model(1000, "13", 1), 1);
    EXPECT_EQ(albini_thermal_model(4, "40", 1), 1) << "under the Scott-Burgan table a code is not Anderson's";
}

TEST(EmberDrift, LandingDistanceIsTheLogProfileIntegral)
{
    // 20 000 kW/m in FM13 lofts a brand to 0.173 (2e4 71.8 5^-0.835)^(1/2)
    // = 106 m; at 0.5 m/s it falls 212 s, drifting 424 m on the uniform 2 m/s
    // wind and further on the log profile. FM13's cap (1000 m) does not bind.
    const Real H_z = case_height(13);
    EXPECT_NEAR(H_z, 105.9, 0.1);
    const Real D = euler_drift(H_z);
    EXPECT_GT(D, 1.2 * U * H_z / V_T) << "the profile carries the brand further than the reference wind";
    SpotCase r(NX, NY);
    const FireParams::SpottingParams sp = r.params();
    compute_albini_spotting(r.phi, r.data, r.wind, r.intensity, r.geom, sp, 1,
                            nullptr, nullptr, 0.5_rt, nullptr, "13", 13, Z_REF, Z0);
    MultiFab one = r.gathered();
    Real xmin, xmax; int n;
    landing_x_range(one, xmin, xmax, n);
    EXPECT_GE(n, 4) << "the brands landed and stamped cells";
    // the launch columns span x = 25..45; landings within spot_radius + a cell
    EXPECT_GE(xmin, 25.0 + D - sp.spot_radius - DX);
    EXPECT_LE(xmax, 45.0 + D + sp.spot_radius + DX);
    EXPECT_LT(xmax - xmin, 20.0 + 2.0 * sp.spot_radius + 2.0 * DX);
}

TEST(EmberDrift, TheCapIsTheSourceCellsFuel)
{
    // on a fuel map the launch cells are FM1 (cap 200 m) while the uniform
    // fuel is FM13: FM1's own fit lofts to 96 m, and the drift of more than
    // 500 m is cut to 200 m
    EXPECT_GT(euler_drift(case_height(1)), 500.0);
    SpotCase r(NX, NY);
    r.fuel.setVal(1.0_rt);
    const FireParams::SpottingParams sp = r.params();
    compute_albini_spotting(r.phi, r.data, r.wind, r.intensity, r.geom, sp, 1,
                            nullptr, nullptr, 0.5_rt, &r.fuel, "13", 13, Z_REF, Z0);
    MultiFab one = r.gathered();
    Real xmin, xmax; int n;
    landing_x_range(one, xmin, xmax, n);
    EXPECT_GE(n, 4);
    EXPECT_GE(xmin, 25.0 + 200.0 - sp.spot_radius - DX);
    EXPECT_LE(xmax, 45.0 + 200.0 + sp.spot_radius + DX);
}

TEST(EmberDrift, LandingsDoNotDependOnTheDecomposition)
{
    SpotCase a(NX, NY), b(25, 10);
    const FireParams::SpottingParams sp = a.params();
    compute_albini_spotting(a.phi, a.data, a.wind, a.intensity, a.geom, sp, 3,
                            nullptr, nullptr, 0.5_rt, nullptr, "13", 13, Z_REF, Z0);
    compute_albini_spotting(b.phi, b.data, b.wind, b.intensity, b.geom, sp, 3,
                            nullptr, nullptr, 0.5_rt, nullptr, "13", 13, Z_REF, Z0);
    MultiFab pa = a.gathered(), pb = b.gathered();
    Real diff = 0.0;
    int n = 0;
    for (MFIter mfi(pa); mfi.isValid(); ++mfi) {
        auto x = pa.const_array(mfi);
        auto y = pb.const_array(mfi);
        for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i) {
                diff = amrex::max(diff, std::abs(x(i, j, 0) - y(i, j, 0)));
                if (x(i, j, 0) < 0.0_rt && !launches(i, j)) { ++n; }
            }
    }
    EXPECT_GE(n, 4);
    EXPECT_EQ(diff, 0.0_rt) << "one box and eight boxes stamp the same landings";
}

namespace {
/// Largest lofting height the spotting call wrote (data component 0)
Real max_height (const MultiFab& data) { return data.max(0); }
}

TEST(EmberDrift, LoftingHeightUsesTheTenMetreWind)
{
    // 8 m/s at 6.1 m is 8 ln(10/0.1)/ln(6.1/0.1) = 8.962 m/s at 10 m, inside
    // the 5-30 m/s range of the fits, so the 10 m conversion sets the height:
    // FM13 at 20 000 kW/m gives 0.173 (2e4 71.8 8.962^-0.835)^(1/2) = 83.0 m
    // (87.0 m if the 6.1 m wind were used unconverted).
    const Real U8 = 8.0;
    const Real U10 = U8 * std::log(10.0 / Z0) / std::log(Z_REF / Z0);
    EXPECT_NEAR(U10, 8.962, 0.001);
    const Real H_expect = 0.173 * std::sqrt(I_B * 71.8 * std::pow(U10, -0.835));
    EXPECT_NEAR(H_expect, 83.0, 0.05);
    SpotCase r(NX, NY);
    r.wind.setVal(U8, 0, 1);
    FireParams::SpottingParams sp = r.params();
    sp.terminal_velocity = 8.0;   // the brands land inside the 1 km domain
    compute_albini_spotting(r.phi, r.data, r.wind, r.intensity, r.geom, sp, 1,
                            nullptr, nullptr, 0.5_rt, nullptr, "13", 13, Z_REF, Z0);
    EXPECT_NEAR(max_height(r.data), H_expect, 1.0e-6 * H_expect);
}

TEST(EmberDrift, BurningStructureUsesTheStandInFit)
{
    // A burning structure on an FM13 cell launches with thermal_fuel_model's
    // fit (FM1 here), not the cell's: Table 4 has no structure fit.
    SpotCase r(NX, NY);
    r.phi.setVal(1.0_rt);
    r.intensity.setVal(0.0_rt);
    MultiFab extra(r.ba, r.dm, 1, 0);
    extra.setVal(0.0_rt);
    for (MFIter mfi(extra); mfi.isValid(); ++mfi) {
        auto e = extra.array(mfi);
        const Box& bx = mfi.validbox();
        if (bx.contains(IntVect(3, 10, 0))) { e(3, 10, 0) = I_B; }
    }
    FireParams::SpottingParams sp = r.params();
    sp.thermal_fuel_model = 1;
    compute_albini_spotting(r.phi, r.data, r.wind, r.intensity, r.geom, sp, 1,
                            nullptr, nullptr, 0.5_rt, &r.fuel, "13", 13, Z_REF, Z0,
                            false, nullptr, nullptr, nullptr, true, &extra);
    EXPECT_NEAR(max_height(r.data), case_height(1), 1.0e-6 * case_height(1));
    EXPECT_GT(case_height(13), 1.05 * case_height(1)) << "the two fits differ, so the test can tell them apart";
}

TEST(EmberDrift, PerCellRoughnessReplacesTheScalar)
{
    // erf.fire.wind_sample_z0_source = surface_layer: the drift and the 10 m
    // wind of the lofting take the roughness per cell. A uniform field of 0.5 m
    // with the scalar argument left at 0.1 m must give the landings and lofting
    // heights of the scalar 0.5 m, and differ from those of the scalar 0.1 m.
    // 8 m/s, so the 10 m wind of the lofting is above the fit's 5 m/s floor
    // and depends on z0; brands that fall at 4 m/s stay inside the domain
    SpotCase f(25, 10), s(25, 10), r(25, 10);
    for (SpotCase* c : {&f, &s, &r}) { c->wind.setVal(8.0_rt, 0, 1); }
    FireParams::SpottingParams sp = f.params();
    sp.terminal_velocity = 4.0;
    MultiFab z0c(f.ba, f.dm, 1, 0);
    z0c.setVal(0.5_rt);
    compute_albini_spotting(f.phi, f.data, f.wind, f.intensity, f.geom, sp, 3,
                            nullptr, nullptr, 0.5_rt, nullptr, "13", 13, Z_REF, Z0,
                            false, nullptr, nullptr, nullptr, true, nullptr, nullptr,
                            nullptr, nullptr, &z0c);
    compute_albini_spotting(s.phi, s.data, s.wind, s.intensity, s.geom, sp, 3,
                            nullptr, nullptr, 0.5_rt, nullptr, "13", 13, Z_REF, 0.5_rt);
    compute_albini_spotting(r.phi, r.data, r.wind, r.intensity, r.geom, sp, 3,
                            nullptr, nullptr, 0.5_rt, nullptr, "13", 13, Z_REF, Z0);
    MultiFab pf = f.gathered(), ps = s.gathered(), pr = r.gathered();
    Real d_same = 0.0, d_scalar = 0.0;
    for (MFIter mfi(pf); mfi.isValid(); ++mfi) {
        auto a = pf.const_array(mfi);
        auto b = ps.const_array(mfi);
        auto c = pr.const_array(mfi);
        for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i) {
                d_same   = amrex::max(d_same,   std::abs(a(i, j, 0) - b(i, j, 0)));
                d_scalar = amrex::max(d_scalar, std::abs(a(i, j, 0) - c(i, j, 0)));
            }
    }
    EXPECT_EQ(d_same, 0.0_rt) << "a uniform 0.5 m field lands brands where the scalar 0.5 m does";
    EXPECT_GT(d_scalar, 0.0_rt) << "and not where the scalar 0.1 m does";
    EXPECT_NE(f.data.max(0), r.data.max(0)) << "the lofting height reads the 10 m wind";
    MultiFab dd(f.ba, f.dm, f.data.nComp(), 0);
    MultiFab::Copy(dd, f.data, 0, 0, dd.nComp(), 0);
    MultiFab::Subtract(dd, s.data, 0, 0, dd.nComp(), 0);
    for (int n = 0; n < dd.nComp(); ++n) {
        EXPECT_EQ(dd.norm0(n), 0.0_rt) << "spotting data component " << n << " as with the scalar 0.5 m";
    }
}

namespace {
/// Ground z = z_flat everywhere except (one low cell) or z = -slope x
MultiFab make_ground (const SpotCase& c, Real z_flat, Real slope, bool low_corner)
{
    MultiFab z(c.ba, c.dm, 1, 0);
    for (MFIter mfi(z); mfi.isValid(); ++mfi) {
        auto a = z.array(mfi);
        const Box& bx = mfi.validbox();
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j)
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                a(i, j, 0) = z_flat - slope * (i + 0.5_rt) * DX;
                if (low_corner && i == NX - 1 && j == 0) { a(i, j, 0) = 0.0_rt; }
            }
    }
    return z;
}

/// Smallest and largest landing distance the launch cells recorded (data component 2)
void launch_distances (const MultiFab& data, Real& dmin, Real& dmax, int& n)
{
    dmin = 1.0e30; dmax = -1.0e30; n = 0;
    for (MFIter mfi(data); mfi.isValid(); ++mfi) {
        auto d = data.const_array(mfi);
        const Box& bx = mfi.validbox();
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j)
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i)
                if (launches(i, j) && d(i, j, 0, 1) > 0.0_rt) {
                    dmin = amrex::min(dmin, d(i, j, 0, 2));
                    dmax = amrex::max(dmax, d(i, j, 0, 2));
                    ++n;
                }
    }
    ParallelDescriptor::ReduceRealMin(dmin);
    ParallelDescriptor::ReduceRealMax(dmax);
    ParallelDescriptor::ReduceIntSum(n);
}

/// The landing distance on the ground z = -slope x (slope * U * f < V_T), by
/// Simpson's rule in the height above the ground h (RK4 for dx/dh, which depends
/// on h only): dx/dh = U f(h) / (-V_T + slope U f(h))
Real slope_drift_exact (Real H_z, Real slope)
{
    const int n = 20000;
    const Real dh = -H_z / n;
    auto g = [=] (Real h) {
        const Real f = albini_wind_profile_factor(h, Z_REF, Z0);
        return U * f / (-V_T + slope * U * f);
    };
    Real h = H_z, x = 0.0;
    for (int k = 0; k < n; ++k) {
        const Real k1 = g(h), k2 = g(h + 0.5_rt * dh), k4 = g(h + dh);
        x += dh * (k1 + 4.0_rt * k2 + k4) / 6.0_rt;
        h += dh;
    }
    return x;
}
}

TEST(EmberDrift, OnAPlateauTheBrandFlightIsTheFlatOne)
{
    // Ground 500 m everywhere but one cell at 0 m far from the flight: the
    // brands fall H_z over n_traj_steps steps, as on flat ground (until Nov
    // 2026 the step was the fall to that lowest cell over n_traj_steps, 30 m
    // here, and the brands landed far beyond the flat-ground distance)
    const Real H_z = case_height(13);
    const Real D = euler_drift(H_z);
    SpotCase r(25, 10);
    MultiFab terrain = make_ground(r, 500.0_rt, 0.0_rt, true);
    const FireParams::SpottingParams sp = r.params();
    compute_albini_spotting(r.phi, r.data, r.wind, r.intensity, r.geom, sp, 1,
                            nullptr, nullptr, 0.5_rt, nullptr, "13", 13, Z_REF, Z0,
                            false, &terrain);
    Real dmin, dmax; int n;
    launch_distances(r.data, dmin, dmax, n);
    const Real tol = (sizeof(Real) == 8 ? 1.0e-9 : 1.0e-4) * D;
    EXPECT_EQ(n, 12) << "every launch cell's brand landed";
    EXPECT_NEAR(dmin, D, tol);
    EXPECT_NEAR(dmax, D, tol);
}

TEST(EmberDrift, OnASlopeTheBrandLandsWhereItMeetsTheGround)
{
    // Uniform slopes along the wind: the brand meets the ground inside its last
    // step and is landed there, by linear interpolation of its height above the
    // ground over the step. Forward Euler with n_traj_steps = 20 then lands
    // 1.7 % (downslope) and 1.9 % (upslope) beyond the Simpson integral, inside
    // the bounds of 2.2 % and 2.5 %; the end of the step lands 2.7 % and 3.9 %
    // too far. (The bounds suit this case only: the crossing's place in the last
    // step sets the end-of-step error. FrozenDrift.OnASlopeUniformWindLandsAt-
    // TheClosedForm pins the landing exactly.) Downslope, the brand also
    // steps on past the n_traj_steps of its loft (n_fall steps).
    const Real H_z = case_height(13);
    for (const Real slope : {0.05_rt, -0.05_rt}) {
        const Real E = slope_drift_exact(H_z, slope);
        SpotCase r(25, 10);
        MultiFab terrain = make_ground(r, 0.0_rt, slope, false);
        const FireParams::SpottingParams sp = r.params();
        compute_albini_spotting(r.phi, r.data, r.wind, r.intensity, r.geom, sp, 1,
                                nullptr, nullptr, 0.5_rt, nullptr, "13", 13, Z_REF, Z0,
                                false, &terrain);
        Real dmin, dmax; int n;
        launch_distances(r.data, dmin, dmax, n);
        EXPECT_EQ(n, 12) << "slope " << slope;
        const Real bound = (slope > 0.0_rt ? 0.022_rt : 0.025_rt) * E;
        EXPECT_GT(dmin, E) << "slope " << slope << ": forward Euler overshoots";
        EXPECT_LT(dmax, E + bound) << "slope " << slope << ": exact " << E;
    }
}

namespace {
/// Launch from column i_src (rows 8-11) only, on the uniform wind (u, 0)
void launch_column (SpotCase& c, int i_src, Real u)
{
    c.wind.setVal(u, 0, 1);
    c.phi.setVal(1.0_rt);
    c.intensity.setVal(0.0_rt);
    for (MFIter mfi(c.phi); mfi.isValid(); ++mfi) {
        auto p = c.phi.array(mfi);
        auto I = c.intensity.array(mfi);
        const Box& bx = mfi.validbox();
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j)
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i)
                if (i == i_src && j >= 8 && j <= 11) { p(i, j, 0) = -1.0_rt; I(i, j, 0) = I_B; }
    }
}

/// Set the x wind of column i alone
void set_column_u (SpotCase& c, int i_col, Real u)
{
    for (MFIter mfi(c.wind); mfi.isValid(); ++mfi) {
        auto w = c.wind.array(mfi);
        const Box& bx = mfi.validbox();
        if (i_col < bx.smallEnd(0) || i_col > bx.bigEnd(0)) { continue; }
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) { w(i_col, j, 0, 0) = u; }
    }
}

/// The x step of the last forward Euler step on flat ground, per unit u
Real euler_last_step (Real H_z)
{
    const Real dz = H_z / N_TRAJ, dt = dz / V_T;
    return albini_wind_profile_factor(dz, Z_REF, Z0) * dt;
}

/// Smallest and largest landing distance recorded in column i_src, and the count
void column_distances (const MultiFab& data, int i_src, Real& dmin, Real& dmax, int& n)
{
    dmin = 1.0e30; dmax = -1.0e30; n = 0;
    for (MFIter mfi(data); mfi.isValid(); ++mfi) {
        auto d = data.const_array(mfi);
        const Box& bx = mfi.validbox();
        for (int j = amrex::max(8, bx.smallEnd(1)); j <= amrex::min(11, bx.bigEnd(1)); ++j)
            if (i_src >= bx.smallEnd(0) && i_src <= bx.bigEnd(0) && d(i_src, j, 0, 1) > 0.0_rt) {
                dmin = amrex::min(dmin, d(i_src, j, 0, 2));
                dmax = amrex::max(dmax, d(i_src, j, 0, 2));
                ++n;
            }
    }
    ParallelDescriptor::ReduceRealMin(dmin);
    ParallelDescriptor::ReduceRealMax(dmax);
    ParallelDescriptor::ReduceIntSum(n);
}
} // namespace

TEST(EmberDrift, ABrandInTheOuterHalfCellLands)
{
    // Brands from the first column (centre x = 5 m) drift 0.5 % of the 2 m/s
    // flight, about 3.1 m, towards the x = 0 face: they land inside the domain,
    // in the outer half of the edge cell (near x = 1.9 m), where the wind is
    // the edge cell's. (Until Nov 2026 a brand past the outermost cell centre
    // left the domain.) Column 99 blows three times as hard: without a periodic
    // x the stencil must not reach it, and with one it must (the seam blends
    // columns 99 and 0), so the brands drift further.
    const Real u = -0.01_rt;
    const Real D = euler_drift(case_height(13)) * (-u / U);
    ASSERT_LT(D, 0.5_rt * DX);
    {
        SpotCase r(25, 10);
        launch_column(r, 0, u);
        set_column_u(r, NX - 1, 3.0_rt * u);
        const FireParams::SpottingParams sp = r.params();
        compute_albini_spotting(r.phi, r.data, r.wind, r.intensity, r.geom, sp, 1,
                                nullptr, nullptr, 0.5_rt, nullptr, "13", 13, Z_REF, Z0);
        Real dmin, dmax; int n;
        column_distances(r.data, 0, dmin, dmax, n);
        EXPECT_EQ(n, 4) << "every brand of the edge column landed";
        EXPECT_NEAR(dmin, D, TOLP * 1.0e3 * D);
        EXPECT_NEAR(dmax, D, TOLP * 1.0e3 * D);
    }
    SpotCase r(25, 10, 1);
    launch_column(r, 0, u);
    set_column_u(r, NX - 1, 3.0_rt * u);
    const FireParams::SpottingParams sp = r.params();
    compute_albini_spotting(r.phi, r.data, r.wind, r.intensity, r.geom, sp, 1,
                            nullptr, nullptr, 0.5_rt, nullptr, "13", 13, Z_REF, Z0);
    Real dmin, dmax; int n;
    column_distances(r.data, 0, dmin, dmax, n);
    EXPECT_EQ(n, 4);
    EXPECT_GT(dmin, 1.05_rt * D) << "periodic x: the seam stencil reads column 99";
    EXPECT_LT(dmax, 3.0_rt * D);
}

TEST(EmberDrift, ABrandThatLandsPastAFaceIsLost)
{
    // Brands from the first column drift 5 m plus half their last step towards
    // the x = 0 face: they start the last step inside the domain and meet the
    // ground half a step outside it, so they are lost. With a periodic x they
    // land on the other side.
    const Real H = case_height(13);
    const Real per_u = euler_drift(H) / U;              // drift per unit wind
    const Real target = 0.5_rt * DX + 0.5_rt * euler_last_step(H) * (0.5_rt * DX / per_u);
    const Real u = -target / per_u;
    {
        SpotCase r(25, 10);
        launch_column(r, 0, u);
        const FireParams::SpottingParams sp = r.params();
        compute_albini_spotting(r.phi, r.data, r.wind, r.intensity, r.geom, sp, 1,
                                nullptr, nullptr, 0.5_rt, nullptr, "13", 13, Z_REF, Z0);
        Real dmin, dmax; int n;
        column_distances(r.data, 0, dmin, dmax, n);
        EXPECT_EQ(n, 0) << "no periodic x: the brands land outside the domain";
    }
    SpotCase r(25, 10, 1);
    launch_column(r, 0, u);
    const FireParams::SpottingParams sp = r.params();
    compute_albini_spotting(r.phi, r.data, r.wind, r.intensity, r.geom, sp, 1,
                            nullptr, nullptr, 0.5_rt, nullptr, "13", 13, Z_REF, Z0);
    Real dmin, dmax; int n;
    column_distances(r.data, 0, dmin, dmax, n);
    EXPECT_EQ(n, 4) << "periodic x: the brands land across the seam";
    EXPECT_NEAR(dmin, target, TOLP * 1.0e3 * target);
}

TEST(EmberDrift, TheSpotDiscReachesAcrossASeam)
{
    // Brands from the second column (centre x = 15 m) drift 18 m and land at
    // x = -3 m, wrapped to 997 m. The 10 m disc reaches the first column
    // (centre x = 5 m, 8 m away across the seam) as well as column 99.
    const Real u = -18.0_rt / (euler_drift(case_height(13)) / U);
    SpotCase r(25, 10, 1);
    launch_column(r, 1, u);
    const FireParams::SpottingParams sp = r.params();
    compute_albini_spotting(r.phi, r.data, r.wind, r.intensity, r.geom, sp, 1,
                            nullptr, nullptr, 0.5_rt, nullptr, "13", 13, Z_REF, Z0);
    MultiFab one = r.gathered();
    int n0 = 0, n99 = 0;
    for (MFIter mfi(one); mfi.isValid(); ++mfi) {
        auto p = one.const_array(mfi);
        for (int j = 8; j <= 11; ++j) {
            if (p(0, j, 0) < 0.0_rt) { ++n0; }
            if (p(NX - 1, j, 0) < 0.0_rt) { ++n99; }
        }
    }
    ParallelDescriptor::ReduceIntSum(n0);
    ParallelDescriptor::ReduceIntSum(n99);
    EXPECT_EQ(n99, 4) << "the disc stamps column 99";
    EXPECT_EQ(n0, 4) << "the disc reaches column 0 across the seam";
}

TEST(EmberDrift, APeriodicSeamWrapsTheBrand)
{
    // Brands from x = 25-45 m drift about 627 m towards -x. Without a periodic
    // x they leave the domain at x = 0; with it they cross the seam, fly on
    // the wind on the other side and land at x = 398-418 m (one 1000 m period
    // on), with the landing distance of the unwrapped flight. (Until Nov 2026
    // they left the domain at the seam.)
    const Real D = euler_drift(case_height(13));
    {
        SpotCase r(25, 10, 0);
        r.wind.setVal(-U, 0, 1);
        const FireParams::SpottingParams sp = r.params();
        compute_albini_spotting(r.phi, r.data, r.wind, r.intensity, r.geom, sp, 1,
                                nullptr, nullptr, 0.5_rt, nullptr, "13", 13, Z_REF, Z0);
        Real dmin, dmax; int n;
        launch_distances(r.data, dmin, dmax, n);
        EXPECT_EQ(n, 0) << "no periodic x: the brands leave the domain";
    }
    SpotCase r(25, 10, 1);
    r.wind.setVal(-U, 0, 1);
    const FireParams::SpottingParams sp = r.params();
    compute_albini_spotting(r.phi, r.data, r.wind, r.intensity, r.geom, sp, 1,
                            nullptr, nullptr, 0.5_rt, nullptr, "13", 13, Z_REF, Z0);
    Real dmin, dmax; int n;
    launch_distances(r.data, dmin, dmax, n);
    const Real tol = (sizeof(Real) == 8 ? 1.0e-9 : 1.0e-4) * D;
    EXPECT_EQ(n, 12) << "periodic x: every brand landed";
    EXPECT_NEAR(dmin, D, tol);
    EXPECT_NEAR(dmax, D, tol);
    MultiFab one = r.gathered();
    Real xmin, xmax; int nb;
    landing_x_range(one, xmin, xmax, nb);
    const Real L = NX * DX;
    EXPECT_GE(nb, 4) << "the landings stamped cells";
    EXPECT_GE(xmin, 25.0 + L - D - sp.spot_radius - DX);
    EXPECT_LE(xmax, 45.0 + L - D + sp.spot_radius + DX);
}
