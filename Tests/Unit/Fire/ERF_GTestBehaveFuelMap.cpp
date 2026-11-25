#include <gtest/gtest.h>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <cmath>
#include <vector>

#include "ERF_FuelModels.H"
#include "ERF_BehaveModel.H"
#include "ERF_FireParams.H"   // ERF_CheneyGouldModel.H reads FireParams without including it
#include "ERF_DirectionalRos.H"

/**
 * @file ERF_GTestBehaveFuelMap.cpp
 * @brief BEHAVE on a two-fuel Scott-Burgan map (SH5 and GR2) against a
 *        closed-form multi-class Rothermel rate of spread.
 *
 * Two defects are pinned. (1) The BEHAVE kernels spread with the uniform fuel
 * on every cell of a fuel map: erf.fire.rothermel_per_fuel now gives them a
 * per-slot state table on the isotropic, per-cell moisture and level-set
 * paths. (2) compute_behave_state used the dead-category SAV as the fuel's
 * characteristic SAV and the category-mean SAVs in W' and M'_f; SH5, a fuel
 * with a coarse dead and a fine live load, spread 0.65 m/s against 1.76 m/s.
 *
 * The reference values come from an independent evaluation of Rothermel
 * (1972) with Albini's (1976) multi-class weighting as Andrews (2018,
 * GTR-RMRS-371) sets it out, written in Python apart from this code:
 *   sigma = sum_i f_i sigma_i over dead and live, f_i the category area shares;
 *   W' = sum_dead w exp(-138/sigma) / sum_live w exp(-500/sigma), per class;
 *   M'_f weighted by w exp(-138/sigma), per class;
 *   net loads with the size-class weights g_ij; heat sink
 *   rho_b sum_i f_i sum_j f_ij exp(-138/sigma_ij) (250 + 1116 M_ij);
 *   wind limited to 0.9 I_R ft/min.
 * Conditions: M_1h 3 %, M_10h 4 %, M_100h 5 %, live herbaceous and woody 60 %
 * (two thirds of the herbaceous load transferred to dead over the 0.30-1.20
 * window), midflame wind 6.5 or 3.0 m/s, flat ground. SB40 loads are
 * Scott and Burgan (2005), 2000/43560 lb/ft2 per ton/acre.
 */

using namespace amrex;

namespace {

constexpr double REL = (sizeof(Real) == 8) ? 1.0e-9 : 1.0e-4;

constexpr Real M1 = 0.03, M10 = 0.04, M100 = 0.05, ML = 0.60;
constexpr Real LO = 0.30, HI = 1.20;

constexpr int SH5 = 145, GR2 = 102;

// Closed-form multi-class rates [m/s] (see the file comment)
constexpr double SH5_R0   = 0.01926557908634948;
constexpr double SH5_U3   = 0.7129771792986564;
constexpr double SH5_U6p5 = 1.7584449850336146;
constexpr double GR2_R0   = 0.009369088043270534;
constexpr double GR2_U3   = 0.3522333821108386;
constexpr double GR2_U6p5 = 0.9489605773204224;   // the 0.9 I_R wind limit binds

BehaveState state (int code)
{
    return compute_behave_state(behave_fuel_params(code, FUEL_SET_SCOTT_BURGAN40),
                                M1, M10, M100, ML, ML, LO, HI);
}

void expect_rel (double got, double want, const char* what)
{
    EXPECT_NEAR(got, want, REL * std::abs(want)) << what;
}

/// Untransferred fuel parameters for every slot, Scott-Burgan set, no custom codes
std::vector<FuelModelParams> slot_params ()
{
    std::vector<FuelModelParams> t(FUEL_SLOT_COUNT);
    for (int s = 0; s < FUEL_SLOT_COUNT; ++s) {
        t[s] = (s >= FUEL_SLOT_CUSTOM_BASE)
            ? get_fuel_params(0, FUEL_SET_SCOTT_BURGAN40, Real(-1.0))   // zero load, as the layer's table
            : behave_fuel_params(fuel_code_from_slot(s), FUEL_SET_SCOTT_BURGAN40);
    }
    return t;
}

template <class T>
Gpu::DeviceVector<T> to_device (const std::vector<T>& h)
{
    Gpu::DeviceVector<T> d(h.size());
    Gpu::copy(Gpu::hostToDevice, h.begin(), h.end(), d.begin());
    Gpu::streamSynchronize();
    return d;
}

/// One row of cells; codes, moistures and a uniform wind along +x
struct Row
{
    BoxArray ba; DistributionMapping dm;
    explicit Row (int n) { ba = BoxArray(Box(IntVect(0, 0, 0), IntVect(n - 1, 0, 0))); dm = DistributionMapping(ba); }
};

void set_codes (MultiFab& fuel, const std::vector<int>& codes)
{
    for (MFIter mfi(fuel); mfi.isValid(); ++mfi) {
        auto f = fuel.array(mfi);
        for (int i = 0; i < static_cast<int>(codes.size()); ++i) { f(i, 0, 0) = Real(codes[i]); }
    }
}

void set_wind (MultiFab& wind, Real U)
{
    wind.setVal(0.0_rt);
    wind.setVal(U, 0, 1);
}

void set_moistures (MultiFab& mc)
{
    mc.setVal(M1, 0, 1); mc.setVal(M10, 1, 1); mc.setVal(M100, 2, 1);
    mc.setVal(ML, 3, 1); mc.setVal(ML, 4, 1);
}

Real cell (const MultiFab& mf, int i)
{
    Real v = -1.0;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) { v = mf.const_array(mfi)(i, 0, 0); }
    return v;
}

} // namespace

// The state of each fuel alone against the closed form: the characteristic
// SAV and the per-class W' and M'_f. SH5 was 0.650 m/s at 6.5 m/s before.
TEST(BehaveFuelMap, StateMatchesClosedFormMultiClassRothermel)
{
    const BehaveState sh5 = state(SH5), gr2 = state(GR2);
    expect_rel(sh5.r_0, SH5_R0, "SH5 no-wind rate");
    expect_rel(gr2.r_0, GR2_R0, "GR2 no-wind rate");
    expect_rel(behave_ros_cell(3.0, 0.0, 0.0, 0.0, sh5), SH5_U3,   "SH5 at 3 m/s");
    expect_rel(behave_ros_cell(6.5, 0.0, 0.0, 0.0, sh5), SH5_U6p5, "SH5 at 6.5 m/s");
    expect_rel(behave_ros_cell(3.0, 0.0, 0.0, 0.0, gr2), GR2_U3,   "GR2 at 3 m/s");
    expect_rel(behave_ros_cell(6.5, 0.0, 0.0, 0.0, gr2), GR2_U6p5, "GR2 at 6.5 m/s");
}

// The table: slot 0 and zero-load slots give no spread, the others match the
// state of their own code.
TEST(BehaveFuelMap, TableHoldsEachSlotsOwnState)
{
    const auto fp = slot_params();
    const auto tbl = build_fuel_behave_table(fp.data(), M1, M10, M100, ML, ML, LO, HI,
                                             true, fire_wind_limit::rothermel, behave_net_load::weighted);
    ASSERT_EQ(static_cast<int>(tbl.size()), FUEL_SLOT_COUNT);
    EXPECT_EQ(tbl[0].r_0, 0.0);
    EXPECT_EQ(tbl[FUEL_SLOT_CUSTOM_BASE].r_0, 0.0);
    expect_rel(tbl[fuel_slot(SH5)].r_0, SH5_R0, "SH5 slot");
    expect_rel(tbl[fuel_slot(GR2)].r_0, GR2_R0, "GR2 slot");
    for (const auto& bs : tbl) { EXPECT_TRUE(std::isfinite(bs.r_0)); }
}

// The isotropic kernel, uniform moisture: each cell spreads with its own fuel.
// Before, every cell took bs_default (here GR2, the deck's fuel_model_id).
// Code 0 does not spread; a code outside the table keeps the uniform fuel.
TEST(BehaveFuelMap, IsotropicKernelReadsTheCellsFuel)
{
    Row g(4);
    MultiFab ros(g.ba, g.dm, 1, 0), wind(g.ba, g.dm, 2, 0), slopes(g.ba, g.dm, 2, 0), fuel(g.ba, g.dm, 1, 0);
    set_codes(fuel, {SH5, GR2, 0, 999});
    set_wind(wind, 6.5_rt);
    slopes.setVal(0.0_rt);

    const auto h_fp = slot_params();
    const auto d_fp = to_device(h_fp);
    const auto d_bs = to_device(build_fuel_behave_table(h_fp.data(), M1, M10, M100, ML, ML, LO, HI,
                                                        true, fire_wind_limit::rothermel,
                                                        behave_net_load::weighted));
    const FuelModelParams fp_uniform = behave_fuel_params(GR2, FUEL_SET_SCOTT_BURGAN40);
    const BehaveState bs_uniform = state(GR2);

    fill_behave_ros(ros, wind, slopes, fp_uniform, bs_uniform, LO, HI, nullptr, false, true,
                    fire_wind_limit::rothermel, behave_net_load::weighted,
                    &fuel, d_bs.data(), d_fp.data(), FUEL_SLOT_COUNT, FUEL_SET_SCOTT_BURGAN40);
    expect_rel(cell(ros, 0), SH5_U6p5, "SH5 cell");
    expect_rel(cell(ros, 1), GR2_U6p5, "GR2 cell");
    EXPECT_EQ(cell(ros, 2), 0.0);
    expect_rel(cell(ros, 3), GR2_U6p5, "unknown code keeps the uniform fuel");

    // without the tables the map is ignored, as before
    fill_behave_ros(ros, wind, slopes, fp_uniform, bs_uniform, LO, HI);
    expect_rel(cell(ros, 0), GR2_U6p5, "no table: uniform fuel");
}

// The per-cell moisture kernel (moisture_dynamic) starts from the cell's own
// untransferred fuel parameters.
TEST(BehaveFuelMap, PerCellMoistureKernelReadsTheCellsFuel)
{
    Row g(4);
    MultiFab ros(g.ba, g.dm, 1, 0), wind(g.ba, g.dm, 2, 0), slopes(g.ba, g.dm, 2, 0),
             fuel(g.ba, g.dm, 1, 0), mc(g.ba, g.dm, 5, 0);
    set_codes(fuel, {SH5, GR2, 0, 999});
    set_wind(wind, 3.0_rt);
    slopes.setVal(0.0_rt);
    set_moistures(mc);

    const auto h_fp = slot_params();
    const auto d_fp = to_device(h_fp);
    const auto d_bs = to_device(build_fuel_behave_table(h_fp.data(), M1, M10, M100, ML, ML, LO, HI,
                                                        true, fire_wind_limit::rothermel,
                                                        behave_net_load::weighted));
    const FuelModelParams fp_uniform = behave_fuel_params(GR2, FUEL_SET_SCOTT_BURGAN40);

    fill_behave_ros(ros, wind, slopes, fp_uniform, state(GR2), LO, HI, &mc, true, true,
                    fire_wind_limit::rothermel, behave_net_load::weighted,
                    &fuel, d_bs.data(), d_fp.data(), FUEL_SLOT_COUNT, FUEL_SET_SCOTT_BURGAN40);
    expect_rel(cell(ros, 0), SH5_U3, "SH5 cell");
    expect_rel(cell(ros, 1), GR2_U3, "GR2 cell");
    EXPECT_EQ(cell(ros, 2), 0.0);
    expect_rel(cell(ros, 3), GR2_U3, "unknown code keeps the uniform fuel");
}

// The level-set path: fill_directional_ros with a BEHAVE state that carries
// the table. The front normal points down the wind, so the normal rate is the
// head rate of each cell's own fuel. Before, the state's uniform bs was used.
TEST(BehaveFuelMap, DirectionalKernelReadsTheCellsFuel)
{
    const int n = 4;
    const Box domain(IntVect(0, 0, 0), IntVect(n - 1, 0, 0));
    RealBox rb({0.0_rt, 0.0_rt, 0.0_rt}, {Real(25 * n), 25.0_rt, 25.0_rt});
    Geometry geom(domain, rb, 0, {0, 0, 0});
    Row g(n);
    MultiFab R(g.ba, g.dm, 1, 0), phi(g.ba, g.dm, 1, 1), wind(g.ba, g.dm, 2, 0),
             slopes(g.ba, g.dm, 2, 0), fuel(g.ba, g.dm, 1, 0);
    set_codes(fuel, {SH5, GR2, 0, 999});
    set_wind(wind, 6.5_rt);
    slopes.setVal(0.0_rt);
    // phi increases along +x (into the unburned fuel), ghost cells included
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        auto p = phi.array(mfi);
        const Box gbx = mfi.growntilebox(1);
        amrex::LoopOnCpu(gbx, [&] (int i, int j, int k) { p(i, j, k) = Real(25 * i); });
    }

    const auto h_fp = slot_params();
    const auto d_bs = to_device(build_fuel_behave_table(h_fp.data(), M1, M10, M100, ML, ML, LO, HI,
                                                        true, fire_wind_limit::rothermel,
                                                        behave_net_load::weighted));
    DirectionalRosState st;
    st.model = DIRECTIONAL_ROS_BEHAVE;
    st.bs = state(GR2);
    st.bs_tbl = d_bs.data();
    st.bs_tbl_size = FUEL_SLOT_COUNT;

    fill_directional_ros(R, phi, wind, slopes, geom, st, nullptr, false,
                         DIRECTIONAL_SHAPE_PROJECTION, DIRECTIONAL_ELLIPSE_LW_MODEL, 8.0_rt,
                         DIRECTIONAL_WIND_COUPLING_PROJECTION, &fuel, nullptr, 0,
                         FUEL_SET_SCOTT_BURGAN40);
    expect_rel(cell(R, 0), SH5_U6p5, "SH5 cell");
    expect_rel(cell(R, 1), GR2_U6p5, "GR2 cell");
    EXPECT_EQ(cell(R, 2), 0.0);
    expect_rel(cell(R, 3), GR2_U6p5, "unknown code keeps the uniform fuel");

    // a state without the table keeps the uniform fuel
    st.bs_tbl = nullptr; st.bs_tbl_size = 0;
    fill_directional_ros(R, phi, wind, slopes, geom, st, nullptr, false,
                         DIRECTIONAL_SHAPE_PROJECTION, DIRECTIONAL_ELLIPSE_LW_MODEL, 8.0_rt,
                         DIRECTIONAL_WIND_COUPLING_PROJECTION, &fuel, nullptr, 0,
                         FUEL_SET_SCOTT_BURGAN40);
    expect_rel(cell(R, 0), GR2_U6p5, "no table: uniform fuel");
}
