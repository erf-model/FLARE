#include <gtest/gtest.h>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Loop.H>
#include <algorithm>
#include <cmath>
#include <functional>

#include "ERF_FireParams.H"
#include "ERF_FrozenWind.H"
#include "ERF_AlbiniSpotting.H"

/**
 * @file ERF_GTestFrozenDrift.cpp
 * @brief Firebrands drifting in the frozen GUST wind
 *        (erf.fire.spotting.drift_wind = frozen).
 *
 * The brand drifts on u, v of the 3D field at its height above the ground;
 * the fall stays at the terminal velocity v_t. In a field whose wind grows
 * linearly with the height above the ground, u = a h, the forward Euler drift
 * over N sub-steps of a brand lofted to H is closed form: the sub-steps
 * sample h_n = H - n dz, n = 0..N-1, with dz = H / N and dt = H / (v_t N), so
 *
 *     D = a dt sum_n (H - n dz) = a dt dz N (N + 1) / 2,
 *
 * which tends to the exact a H^2 / (2 v_t) as N grows. Linear interpolation
 * in height reproduces u = a h exactly between the first and the top cell
 * centres, so the run must land the brand at D to round-off.
 *
 * With drift_wind = frozen_uvw the brand also moves with the field's w. In a
 * uniform u = U and w = W on flat ground it falls (v_t - W) dt per step of
 * the v_t fall, dt = H / (v_t N), so it lands after
 *
 *     n = ceil(N v_t / (v_t - W))   steps,   D = U n dt,
 *
 * if n <= 10 N, and is dropped as aloft otherwise (or if W >= v_t).
 */

using namespace amrex;

namespace {

constexpr int  NX = 100, NY = 20;          // fire grid, 10 m cells
constexpr Real DX = 10.0;
constexpr Real Z_REF = 6.1, Z0 = 0.1, V_T = 0.5, I_B = 2.0e4, U_REF = 2.0;
constexpr int  N_TRAJ = 20;
constexpr Real SHEAR = 0.02;               // a [1/s]: u = a h
constexpr Real GDX = 50.0, GDZ = 5.0;      // GUST cells
constexpr double TOLP = (sizeof(Real) == 8) ? 1.0e-12 : 1.0e-5;
constexpr double TOLD = (sizeof(Real) == 8) ? 1.0e-9 : 1.0e-4;   // relative, on a landing distance

bool launches (int i, int j) { return i >= 2 && i <= 4 && j >= 8 && j <= 11; }

/// GUST columns over the fire region: ground(x, y), centres ground + (k + 1/2) GDZ,
/// u = a h and v = 0 (h the height above the column ground)
/// (and, given wh, the vertical wind w = wh(h) as component 4; u = uh(h) if given)
void make_columns (FrozenWindHost& h, int nz, const std::function<Real(Real, Real)>& ground,
                   Real a = SHEAR, Real b = 0.0, const std::function<Real(Real)>& wh = nullptr,
                   const std::function<Real(Real)>& uh = nullptr)
{
    const int gnx = static_cast<int>(NX * DX / GDX), gny = static_cast<int>(NY * DX / GDX);
    const Box gdom(IntVect(0, 0, 0), IntVect(gnx - 1, gny - 1, nz - 1));
    const int nc = wh ? 5 : 4;
    h.data.assign(static_cast<std::size_t>(nc * gdom.numPts()), 0.0);
    const Dim3 lo = lbound(gdom), hi = ubound(gdom);
    Array4<Real> out(h.data.data(), lo, Dim3{hi.x + 1, hi.y + 1, hi.z + 1}, nc);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < gny; ++j)
            for (int i = 0; i < gnx; ++i) {
                const Real zg = ground((i + 0.5) * GDX, (j + 0.5) * GDX);
                const Real hk = (k + 0.5) * GDZ;
                out(i, j, k, 0) = uh ? uh(hk) : a * hk;
                out(i, j, k, 1) = b * hk;
                out(i, j, k, 2) = zg + hk;
                out(i, j, k, 3) = zg;
                if (wh) { out(i, j, k, 4) = wh(hk); }
            }
    FrozenColumnsView& c = h.view;
    c.a = Array4<const Real>(h.data.data(), lo, Dim3{hi.x + 1, hi.y + 1, hi.z + 1}, nc);
    c.gilo = 0; c.gihi = gnx - 1; c.gjlo = 0; c.gjhi = gny - 1; c.kz0 = 0; c.kz1 = nz - 1;
    c.gx0 = 0.0; c.gy0 = 0.0; c.gdx = GDX; c.gdy = GDX;
}

Real flat (Real, Real) { return 0.0; }

FrozenDrift drift_of (const FrozenWindHost& h, bool vertical = false)
{
    FrozenDrift d;
    d.cols = h.view; d.bilinear = true; d.below_first_log = false; d.z0_log = Z0;
    d.vertical = vertical;
    return d;
}

// Largest |a - b| over the first nc components. a and b may sit on different
// box layouts or rank maps: both are copied onto one shared layout first
// (two DistributionMappings built on the same BoxArray can map boxes to
// different ranks, so a direct MultiFab::Subtract is valid on one rank only).
Real max_diff (const MultiFab& a, const MultiFab& b, int nc)
{
    BoxArray one(a.boxArray().minimalBox());
    const DistributionMapping dm(one);
    MultiFab pa(one, dm, nc, 0), pb(one, dm, nc, 0);
    pa.ParallelCopy(a, 0, 0, nc);
    pb.ParallelCopy(b, 0, 0, nc);
    MultiFab::Subtract(pa, pb, 0, 0, nc, 0);
    return pa.norminf(0, nc, IntVect(0));
}

struct SpotCase
{
    BoxArray ba; DistributionMapping dm; Geometry geom;
    MultiFab phi, data, wind, intensity;
    SpotCase (int max_x, int max_y)
    {
        Box domain(IntVect(0, 0, 0), IntVect(NX - 1, NY - 1, 0));
        ba = BoxArray(domain);
        ba.maxSize(IntVect(max_x, max_y, 1));
        dm = DistributionMapping(ba);
        geom = Geometry(domain, RealBox(0.0, 0.0, 0.0, NX * DX, NY * DX, 1.0), CoordSys::cartesian, {0, 0, 0});
        phi.define(ba, dm, 1, 1); data.define(ba, dm, 4, 0); wind.define(ba, dm, 2, 0);
        intensity.define(ba, dm, 1, 0);
        phi.setVal(1.0_rt); data.setVal(0.0_rt); wind.setVal(U_REF, 0, 1); wind.setVal(0.0_rt, 1, 1);
        intensity.setVal(0.0_rt);
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
        sp.launch_from = "burned"; sp.drift_wind = "frozen";
        return sp;
    }
    int n_above_top = -1;   ///< brands that flew above GUST's top centre, from the last run
    int n_aloft     = -1;   ///< brands dropped as still aloft, from the last run
    void run (const FrozenDrift* d, const MultiFab* terrain = nullptr)
    {
        compute_albini_spotting(phi, data, wind, intensity, geom, params(), 1,
                                nullptr, nullptr, 0.5_rt, nullptr, "13", 13, Z_REF, Z0,
                                false, terrain, nullptr, nullptr, true, nullptr, d, &n_above_top,
                                &n_aloft);
    }
    int launched () const { return static_cast<int>(data.sum(1) + 0.5); }
};

/// FM13 at 20 000 kW/m with the 2 m/s reference wind (U10 below the 5 m/s floor of the fits)
Real lofting_height ()
{
    Real A = 0.0, B = 0.0;
    EXPECT_TRUE(albini_thermal_fit(13, A, B));
    return 0.173 * std::sqrt(I_B * A * std::pow(5.0, B));
}

/// Forward Euler drift of a brand lofted to H on flat ground in u = a min(h, h_top)
Real euler_shear_drift (Real H, Real h_top)
{
    const Real dt = (H / V_T) / N_TRAJ, dz = V_T * dt;
    Real d = 0.0;
    for (int n = 0; n < N_TRAJ; ++n) { d += SHEAR * std::min(H - n * dz, h_top) * dt; }
    return d;
}

} // namespace

TEST(FrozenDrift, PointWindFollowsTheTerrain)
{
    // A sloping, bumpy ground: each column is read at its own ground + h, so
    // u = a h wherever the point is (sampling at the absolute height h would
    // give a (h - ground) instead).
    auto ground = [](Real x, Real y) { return 0.05 * x + 20.0 * std::sin(y / 60.0); };
    FrozenWindHost h;
    make_columns(h, 40, ground, SHEAR, -0.5 * SHEAR);
    for (Real x : {60.0, 333.3, 712.0}) {
        for (Real y : {30.0, 101.0, 170.0}) {
            for (Real hh : {7.5, 41.2, 150.0}) {
                Real u, v, zs, gs, ex;
                frozen_point_wind(h.view, x, y, hh, true, false, Z0, u, v, zs, gs, ex);
                EXPECT_NEAR(u, SHEAR * hh, TOLP * 10.0) << "x " << x << " y " << y << " h " << hh;
                EXPECT_NEAR(v, -0.5 * SHEAR * hh, TOLP * 10.0);
                EXPECT_NEAR(zs - gs, hh, TOLP * 1000.0) << "sampled at the column ground + h";
                EXPECT_LT(ex, 0.0) << "below the top centre";
            }
        }
    }
}

TEST(FrozenDrift, PointWindBelowTheFirstCentreAndAboveTheTop)
{
    FrozenWindHost h;
    make_columns(h, 10, flat);     // centres 2.5 .. 47.5 m
    Real u, v, zs, gs, ex;
    // below the first centre: its wind as it is, or log-scaled from it
    frozen_point_wind(h.view, 300.0, 100.0, 1.0, true, false, Z0, u, v, zs, gs, ex);
    EXPECT_NEAR(u, SHEAR * 2.5, TOLP);
    frozen_point_wind(h.view, 300.0, 100.0, 1.0, true, true, Z0, u, v, zs, gs, ex);
    EXPECT_NEAR(u, SHEAR * 2.5 * std::log(1.0 / Z0) / std::log(2.5 / Z0), TOLP);
    // above the top centre: the top centre's wind, and the excess
    frozen_point_wind(h.view, 300.0, 100.0, 80.0, true, false, Z0, u, v, zs, gs, ex);
    EXPECT_NEAR(u, SHEAR * 47.5, TOLP);
    EXPECT_NEAR(ex, 80.0 - 47.5, TOLP * 100.0);
}

TEST(FrozenDrift, LandingDistanceIsTheShearClosedForm)
{
    const Real H = lofting_height();
    EXPECT_NEAR(H, 105.9, 0.1);
    const Real dz = H / N_TRAJ, dt = H / (V_T * N_TRAJ);
    const Real D = SHEAR * dt * dz * N_TRAJ * (N_TRAJ + 1) / 2.0;
    EXPECT_NEAR(D, euler_shear_drift(H, 1.0e30), TOLD * D);
    EXPECT_NEAR(D, SHEAR * H * H / (2.0 * V_T) * (N_TRAJ + 1.0) / N_TRAJ, TOLD * D)
        << "the Euler sum is the exact a H^2 / (2 v_t) times (N + 1) / N";
    EXPECT_GT(dz, 2.5) << "every sub-step samples above the first centre (2.5 m), where u = a h is exact";

    FrozenWindHost h;
    make_columns(h, 40, flat);     // top centre 197.5 m, above H
    const FrozenDrift d = drift_of(h);
    SpotCase r(NX, NY);
    r.run(&d);
    EXPECT_NEAR(r.data.max(0), H, 1.0e-6 * H) << "the lofting height still takes the reference wind";
    // every source cell lands its brands D downwind (component 2: largest landing distance)
    EXPECT_NEAR(r.data.max(2), D, TOLD * D);
    EXPECT_GT(D, 200.0);

    // the profile drift on the same reference wind lands the brands elsewhere
    SpotCase p(NX, NY);
    p.run(nullptr);
    EXPECT_GT(std::abs(p.data.max(2) - D), 100.0) << "the test tells the two drifts apart";
}

TEST(FrozenDrift, OnRaisedGroundTheBrandReadsItsHeightAboveTheGround)
{
    // The same case on a plateau 300 m up (the fire's terrain and GUST's
    // ground both at 300 m): the brand reads the wind at its height above the
    // ground, so it lands at the same closed-form distance. Reading it at its
    // absolute height (300 m + h) would put every sub-step above GUST's top.
    const Real H = lofting_height();
    const Real D = euler_shear_drift(H, 1.0e30);
    FrozenWindHost h;
    make_columns(h, 40, [](Real, Real) { return Real(300.0); });
    const FrozenDrift d = drift_of(h);
    SpotCase r(NX, NY);
    MultiFab terrain(r.ba, r.dm, 1, 0);
    terrain.setVal(300.0_rt);
    r.run(&d, &terrain);
    EXPECT_NEAR(r.data.max(2), D, TOLD * D);
}

TEST(FrozenDrift, AboveTheTopTheBrandTakesTheTopWind)
{
    // GUST top centre at 47.5 m, the brand lofted to 106 m: above 47.5 m it
    // drifts on the top centre's wind
    const Real H = lofting_height();
    FrozenWindHost h;
    make_columns(h, 10, flat);
    const FrozenDrift d = drift_of(h);
    SpotCase r(NX, NY);
    r.run(&d);
    const Real D = euler_shear_drift(H, 47.5);
    EXPECT_NEAR(r.data.max(2), D, TOLD * D);
    EXPECT_LT(D, euler_shear_drift(H, 1.0e30) - 50.0) << "the cap at the top centre matters here";
    // every brand starts above 47.5 m, so every one is counted and reported
    const auto launched = static_cast<int>(r.data.sum(1) + 0.5);
    EXPECT_GT(launched, 0);
    EXPECT_EQ(r.n_above_top, launched);

    // with the top centre above H no brand is counted
    FrozenWindHost tall;
    make_columns(tall, 40, flat);
    const FrozenDrift dt = drift_of(tall);
    SpotCase q(NX, NY);
    q.run(&dt);
    EXPECT_EQ(q.n_above_top, 0);
}

TEST(FrozenDrift, GatherCopiesTheColumnsUnderTheFireRegion)
{
    // A GUST mesh of 40 x 24 x 6 cells of 25 x 25 m, split into 8 x 8 column
    // boxes; the cell-centred fields hold their own indices and the node
    // heights a sloping ground, so every gathered value is known
    const Box gdom(IntVect(0, 0, 0), IntVect(39, 23, 5));
    const Geometry ggeom(gdom, RealBox(0.0, 0.0, 0.0, 1000.0, 600.0, 60.0), CoordSys::cartesian, {0, 0, 0});
    BoxArray ba(gdom);
    ba.maxSize(IntVect(16, 12, 6));          // not the 8 x 8 of the columns, so the copy re-maps
    const DistributionMapping dm(ba);
    MultiFab u(ba, dm, 1, 0), v(ba, dm, 1, 0), w(ba, dm, 1, 0), z(ba, dm, 1, 0);
    MultiFab znd(amrex::convert(ba, IndexType::TheNodeType()), dm, 1, 0);
    auto zn = [](int i, int j, int k) { return Real(0.3 * i + 0.2 * j + 10.0 * k); };
    auto zg = [](int i, int j) { return Real(0.3 * (i + 0.5) + 0.2 * (j + 0.5)); };
    auto uu = [](int i, int j, int k) { return Real(i + 100 * j + 10000 * k); };
    for (MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto au = u.array(mfi); auto av = v.array(mfi); auto aw = w.array(mfi); auto az = z.array(mfi);
        const Box& bx = mfi.validbox();
        LoopOnCpu(bx, [&] (int i, int j, int k) {
            au(i, j, k) = uu(i, j, k); av(i, j, k) = -2.0 * uu(i, j, k); aw(i, j, k) = 3.0 * uu(i, j, k);
            az(i, j, k) = zg(i, j) + 10.0 * (k + 0.5);
        });
        auto an = znd.array(mfi);
        LoopOnCpu(amrex::surroundingNodes(bx), [&] (int i, int j, int k) { an(i, j, k) = zn(i, j, k); });
    }
    FrozenWind fw;
    fw.define_from_fields(ggeom, "gather test", u, v, w, z, znd, 8);

    auto check = [&] (const Geometry& fgeom, const Box& expect, bool with_w = false) {
        FrozenWindHost h;
        fw.gather_columns(h, fgeom, with_w);
        const int nc = with_w ? 5 : 4;
        const FrozenColumnsView& c = h.view;
        const Dim3 lo = amrex::lbound(c.a), hi = amrex::ubound(c.a);
        EXPECT_EQ(lo.x, expect.smallEnd(0)); EXPECT_EQ(hi.x, expect.bigEnd(0));
        EXPECT_EQ(lo.y, expect.smallEnd(1)); EXPECT_EQ(hi.y, expect.bigEnd(1));
        EXPECT_EQ(lo.z, 0); EXPECT_EQ(hi.z, 5);
        EXPECT_EQ(c.gilo, 0); EXPECT_EQ(c.gihi, 39); EXPECT_EQ(c.gjlo, 0); EXPECT_EQ(c.gjhi, 23);
        EXPECT_EQ(c.kz0, 0); EXPECT_EQ(c.kz1, 5);
        EXPECT_EQ(h.data.size(), static_cast<std::size_t>(nc * expect.numPts()));
        EXPECT_EQ(c.a.nComp(), nc);
        int bad = 0;
        LoopOnCpu(expect, [&] (int i, int j, int k) {
            bad += (c.a(i, j, k, 0) != uu(i, j, k)) + (c.a(i, j, k, 1) != -2.0 * uu(i, j, k))
                 + (std::abs(c.a(i, j, k, 2) - (zg(i, j) + 10.0 * (k + 0.5))) > TOLP * 100.0)
                 + (std::abs(c.a(i, j, k, 3) - zg(i, j)) > TOLP * 100.0)
                 + (with_w && c.a(i, j, k, 4) != 3.0 * uu(i, j, k));
        });
        EXPECT_EQ(bad, 0) << "gathered u, v, z, ground or w differs from the GUST field";
    };
    // a fire region inside: x 200..500 m, y 100..300 m -> columns 6..21, 2..13 (one cell
    // of margin below, two above, as sample_reference_wind reads)
    const Box fdom(IntVect(0, 0, 0), IntVect(29, 19, 0));
    check(Geometry(fdom, RealBox(200.0, 100.0, 0.0, 500.0, 300.0, 1.0), CoordSys::cartesian, {0, 0, 0}),
          Box(IntVect(6, 2, 0), IntVect(21, 13, 5)));
    // the fire region is GUST's whole extent: the margin is clipped at GUST's edges
    check(Geometry(fdom, RealBox(0.0, 0.0, 0.0, 1000.0, 600.0, 1.0), CoordSys::cartesian, {0, 0, 0}), gdom);
    // with w (drift_wind = frozen_uvw) as component 4
    check(Geometry(fdom, RealBox(200.0, 100.0, 0.0, 500.0, 300.0, 1.0), CoordSys::cartesian, {0, 0, 0}),
          Box(IntVect(6, 2, 0), IntVect(21, 13, 5)), true);
}

TEST(FrozenDrift, LandingsDoNotDependOnTheDecomposition)
{
    FrozenWindHost h;
    make_columns(h, 40, flat, SHEAR, 0.1 * SHEAR);   // a cross-wind as well
    const FrozenDrift d = drift_of(h);
    SpotCase a(NX, NY), b(25, 10);
    a.run(&d);
    b.run(&d);
    EXPECT_EQ(max_diff(a.phi, b.phi, 1), 0.0_rt);
    EXPECT_GT(a.data.max(2), 200.0);
}

TEST(FrozenDrift, PointWindSamplesTheVerticalWind)
{
    // w = c h + d over a sloping, bumpy ground: each column is read at its own
    // ground + h, linearly between centres, so w is exact there; below the
    // first centre it is clamped or log-scaled with u and v; above the top
    // centre it takes the top centre's w
    const Real c = -0.004, d0 = 0.3;
    auto ground = [](Real x, Real y) { return 0.05 * x + 20.0 * std::sin(y / 60.0); };
    FrozenWindHost h;
    make_columns(h, 10, ground, SHEAR, 0.0, [=](Real hh) { return c * hh + d0; });   // centres 2.5 .. 47.5 m
    for (Real x : {60.0, 333.3, 712.0}) {
        for (Real hh : {7.5, 21.2, 44.0}) {
            Real u, v, zs, gs, ex, w = -99.0;
            frozen_point_wind(h.view, x, 101.0, hh, true, false, Z0, u, v, zs, gs, ex, &w);
            EXPECT_NEAR(w, c * hh + d0, TOLP * 10.0) << "x " << x << " h " << hh;
            EXPECT_NEAR(u, SHEAR * hh, TOLP * 10.0) << "u unchanged by asking for w";
        }
    }
    Real u, v, zs, gs, ex, w = -99.0;
    const Real w1 = c * 2.5 + d0;
    frozen_point_wind(h.view, 300.0, 100.0, 1.0, true, false, Z0, u, v, zs, gs, ex, &w);
    EXPECT_NEAR(w, w1, TOLP);
    frozen_point_wind(h.view, 300.0, 100.0, 1.0, true, true, Z0, u, v, zs, gs, ex, &w);
    EXPECT_NEAR(w, w1 * std::log(1.0 / Z0) / std::log(2.5 / Z0), TOLP);
    EXPECT_NEAR(w / u, w1 / (SHEAR * 2.5), TOLP * 10.0) << "the log rule keeps the wind's direction";
    frozen_point_wind(h.view, 300.0, 100.0, 80.0, false, false, Z0, u, v, zs, gs, ex, &w);
    EXPECT_NEAR(w, c * 47.5 + d0, TOLP);
}

namespace {
constexpr Real U_DRIFT = 0.5;   ///< uniform u of the frozen_uvw cases [m/s]

/// The frozen_uvw run in u = U, w = W on flat ground
void run_uniform (SpotCase& r, FrozenWindHost& h, Real W, Real U = U_DRIFT)
{
    make_columns(h, 40, flat, 0.0, 0.0, [=](Real) { return W; }, [=](Real) { return U; });
    const FrozenDrift d = drift_of(h, true);
    r.run(&d);
}

/// Steps to land with fall (v_t - W) dt per step of dt = H / (v_t N)
int steps_to_land (Real W)
{
    return static_cast<int>(std::ceil(N_TRAJ * V_T / (V_T - W) - 1.0e-9));
}
} // namespace

TEST(FrozenDrift, VerticalWindChangesTheFallClosedForm)
{
    const Real H  = lofting_height();
    const Real dt = H / (V_T * N_TRAJ);
    // W = 0, an updraft of half v_t, a downdraft of v_t, an updraft of 0.85 v_t
    // (134 steps, past N = 20 but within 10 N)
    for (Real W : {0.0_rt, 0.5_rt * V_T, -V_T, 0.85_rt * V_T}) {
        const int  n = steps_to_land(W);
        const Real D = U_DRIFT * n * dt;
        SpotCase r(NX, NY);
        FrozenWindHost h;
        run_uniform(r, h, W);
        EXPECT_GT(r.launched(), 0);
        EXPECT_EQ(r.n_aloft, 0) << "W " << W << ": every brand lands within 10 N steps";
        EXPECT_NEAR(r.data.max(2), D, TOLD * D) << "W " << W << ": n = " << n << " steps";
    }
    EXPECT_EQ(steps_to_land(0.5 * V_T), 2 * N_TRAJ);
    EXPECT_EQ(steps_to_land(-V_T), N_TRAJ / 2);
    EXPECT_GT(steps_to_land(0.85 * V_T), N_TRAJ);
    EXPECT_LE(steps_to_land(0.85 * V_T), albini_aloft_step_factor * N_TRAJ);
}

TEST(FrozenDrift, BrandsStillAloftAreDroppedAndCounted)
{
    // W = 0.95 v_t needs 400 steps (> 10 N = 200); W = v_t never lands. With
    // u = 0 the brands stay over their source, so none leaves the fire grid.
    for (Real W : {0.95_rt * V_T, V_T}) {
        EXPECT_TRUE(W >= V_T || steps_to_land(W) > albini_aloft_step_factor * N_TRAJ);
        SpotCase r(NX, NY);
        FrozenWindHost h;
        run_uniform(r, h, W, 0.0);
        EXPECT_EQ(r.n_aloft, 12) << "W " << W << ": the brand of every source cell is still aloft";
        EXPECT_EQ(r.data.max(2), 0.0_rt) << "no landing distance";
        // no spot fire: only the launching cells burn
        int burning = 0;
        for (MFIter mfi(r.phi); mfi.isValid(); ++mfi) {
            const auto p = r.phi.const_array(mfi);
            LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { burning += (p(i, j, k) < 0.0); });
        }
        ParallelDescriptor::ReduceIntSum(burning);
        EXPECT_EQ(burning, 12) << "W " << W;
    }
}

constexpr int NSLOPES = 400;  // grounds swept by VerticalWindOfZeroLandsAsFrozen

TEST(FrozenDrift, VerticalWindOfZeroLandsAsFrozen)
{
    // w = 0 everywhere: frozen_uvw lands every brand where frozen does, bit
    // for bit, over sloping, bumpy grounds (fire terrain and GUST ground
    // alike). A brand that passes within 1e-3 dz of a slope before its last
    // step is about one landing in a thousand, so many slopes are swept:
    // a landing tolerance applied with w = 0 makes some of them land a step
    // early.
    int n_cases = 0;
    Real longest = 0.0;
    for (int m = 0; m < NSLOPES; ++m) {
        const Real slope = 0.01 + 0.19 * m / NSLOPES;
        const Real phase = 0.37 * m;
        auto ground = [=](Real x, Real y) { return slope * x + 20.0 * std::sin(y / 60.0 + phase); };
        FrozenWindHost h;
        make_columns(h, 40, ground, SHEAR, 0.1 * SHEAR, [](Real) { return Real(0.0); });
        SpotCase a(NX, NY), b(NX, NY);
        // each case's terrain on its own layout (the two rank maps may differ)
        auto terrain_of = [&] (const SpotCase& c) {
            MultiFab t(c.ba, c.dm, 1, 0);
            for (MFIter mfi(t); mfi.isValid(); ++mfi) {
                auto ta = t.array(mfi);
                LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { ta(i, j, k) = ground((i + 0.5) * DX, (j + 0.5) * DX); });
            }
            return t;
        };
        MultiFab ta = terrain_of(a), tb = terrain_of(b);
        const FrozenDrift d2 = drift_of(h, false), d3 = drift_of(h, true);
        a.run(&d2, &ta);
        b.run(&d3, &tb);
        ASSERT_EQ(b.n_aloft, 0) << "slope " << slope;
        longest = std::max(longest, a.data.max(2));
        ASSERT_EQ(max_diff(a.data, b.data, 4), 0.0_rt) << "slope " << slope << ", phase " << phase;
        ASSERT_EQ(max_diff(a.phi, b.phi, 1), 0.0_rt) << "slope " << slope << ", phase " << phase;
        ++n_cases;
    }
    EXPECT_EQ(n_cases, NSLOPES);
    EXPECT_GT(longest, 100.0);
}
