#include <gtest/gtest.h>
#include <cmath>
#include <AMReX_REAL.H>
#include <ERF_Constants.H>

/// Round-off tolerance of the build precision: 1e-12 in double, 1e-5 in single.
static constexpr double TOL = (sizeof(amrex::Real) == 8) ? 1e-12 : 1e-5;

#include "ERF_LevelSetEllipse.H"

/**
 * @file ERF_GTestLevelSetEllipse.cpp
 * @brief The spread ellipse of the level-set path: rates at the head, flanks
 *        and back, the length-to-width and head-to-back ratios, the calm-wind
 *        circle and the cap on the ratio.
 */

TEST(LevelSetEllipse, HeadFlankBackRates)
{
    const amrex::Real R = 0.5;          // head rate [m/s]
    const amrex::Real U = 3.0;          // midflame wind [m/s]
    const SpreadEllipse e = spread_ellipse(R, U, 0.0, 8.0);
    const amrex::Real LB = anderson_LW_ratio(U);
    const amrex::Real s  = std::sqrt(LB * LB - 1.0);
    const amrex::Real HB = (LB + s) / (LB - s);
    EXPECT_NEAR(e.LB, LB, TOL);
    EXPECT_NEAR(e.HB, HB, (TOL * 1e3));
    // head, flank and back from the support function
    EXPECT_NEAR(ellipse_normal_speed(e, 1.0, 0.0), R, TOL);
    EXPECT_NEAR(ellipse_normal_speed(e, 0.0, 1.0), e.a, TOL);
    EXPECT_NEAR(ellipse_normal_speed(e, -1.0, 0.0), R / HB, TOL);
    // length over width of the envelope after unit time is LB
    const amrex::Real length = ellipse_normal_speed(e, 1.0, 0.0) + ellipse_normal_speed(e, -1.0, 0.0);
    const amrex::Real width  = 2.0 * ellipse_normal_speed(e, 0.0, 1.0);
    EXPECT_NEAR(length / width, LB, (TOL * 1e3));
    // the normal speed never goes negative around the ellipse
    for (int n = 0; n <= 360; n += 5) {
        const amrex::Real th = n * PI / 180.0;
        EXPECT_GE(ellipse_normal_speed(e, std::cos(th), std::sin(th)), 0.0);
    }
}

TEST(LevelSetEllipse, CalmWindIsCircle)
{
    const SpreadEllipse e = spread_ellipse(0.5, 0.0, 0.0, 8.0);
    EXPECT_NEAR(e.LB, 1.0, TOL);
    EXPECT_NEAR(e.a, 0.5, TOL);
    EXPECT_NEAR(e.c, 0.0, TOL);
    for (int n = 0; n < 360; n += 30) {
        const amrex::Real th = n * PI / 180.0;
        EXPECT_NEAR(ellipse_normal_speed(e, std::cos(th), std::sin(th)), 0.5, TOL);
    }
}

TEST(LevelSetEllipse, FixedRatioAndCap)
{
    const SpreadEllipse f = spread_ellipse(1.0, 10.0, 3.0, 8.0);   // fixed LB = 3 regardless of the wind
    EXPECT_NEAR(f.LB, 3.0, TOL);
    const SpreadEllipse c = spread_ellipse(1.0, 30.0, 0.0, 4.0);   // Anderson would give 8 (capped) at 30 m/s
    EXPECT_NEAR(c.LB, 4.0, TOL);
    const SpreadEllipse d = spread_ellipse(1.0, 30.0, 0.0, 8.0);   // Anderson's own cap
    EXPECT_NEAR(d.LB, 8.0, TOL);
}

// Motivation: anderson_LW_ratio takes the midflame wind in m/s (Finney 1998,
// RMRS-RP-4, eq. 13). Until November 2026 every caller passed mi/h, which gave
// L/W 1.70 for the BehavePlus case below and reached the cap of 8 at 3.8 m/s.
// Pins: BehavePlus (Andrews 2018, RMRS-GTR-371) gives L/W 1.24 for a 2.52 mi/h
// (1.1265 m/s) midflame wind; the BehavePlus mi/h form 0.936 e^{0.1147 U} +
// 0.461 e^{-0.0692 U} - 0.397 gives 1.3439 at 1.5 m/s (3.355 mi/h) and reaches
// 8 at 8.49 m/s.
TEST(LevelSetEllipse, AndersonRatioTakesMetresPerSecond)
{
    EXPECT_NEAR(anderson_LW_ratio(1.12654), 1.24, 5.0e-3);
    EXPECT_NEAR(anderson_LW_ratio(1.5), 1.34385, 1.0e-3);
    EXPECT_NEAR(anderson_LW_ratio(0.0), 1.0, TOL);
    EXPECT_LT(anderson_LW_ratio(8.4), 8.0);
    EXPECT_NEAR(anderson_LW_ratio(8.6), 8.0, TOL);
    // spread_ellipse applies the same ratio to the m/s wind it is given
    EXPECT_NEAR(spread_ellipse(1.0, 1.5, 0.0, 8.0).LB, anderson_LW_ratio(1.5), TOL);
}
