#ifndef ERF_FIRE_WIND_ADJUST_H
#define ERF_FIRE_WIND_ADJUST_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>
#include <cmath>

/**
 * @file ERF_FireWindAdjust.H
 * @brief What the native fire engine does to the wind it samples
 *
 * The native engine's corrections to the reference wind that the shared
 * sampler (ERF_FireWindExtract.H) hands it:
 * 1. Wind Adjustment Factor (WAF) for forest canopy effects
 * 2. FARSITE terrain wind corrections (ridge, shelter, valley, deflection)
 *
 * References:
 * - Rothermel (1972): USDA Forest Service Research Paper INT-115
 * - Dyer (1974): Chapter 2 in Modelling of Atmospheric Flow Fields
 * - Andrews (1986): BEHAVE: Fire Behavior Prediction and Fuel Modeling System
 * - FARSITE documentation (Finney 1998)
 */

/**
 * @brief Ratio of the neutral log-profile winds at two heights, U(z_to) / U(z_from).
 *
 * ln(z_to / z0) / ln(z_from / z0) with roughness length z0 [m]; both heights are
 * floored just above z0. The empirical spread models are fitted to the wind at a
 * fixed height (the WAF to the 20 ft wind, Cruz et al. 2005, FBP and Cheney et al.
 * 1998 to the 10 m open wind), and this brings the wind sampled at
 * erf.fire.wind_ref_ht to that height.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real neutral_wind_height_factor(amrex::Real z_to, amrex::Real z_from, amrex::Real z0) noexcept
{
    const amrex::Real z0s = amrex::max(z0, amrex::Real(1.0e-6));
    const amrex::Real zt  = amrex::max(z_to,   z0s * amrex::Real(1.000001));
    const amrex::Real zf  = amrex::max(z_from, z0s * amrex::Real(1.000001));
    return std::log(zt / z0s) / std::log(zf / z0s);
}

/**
 * @brief Compute unsheltered Wind Adjustment Factor (Andrews/Albini & Baughman 1979)
 *
 * @param[in] delta_ft Fuel bed depth [ft]
 * @return WAF, the midflame wind over the 20 ft (6.1 m) open wind, for fuels in
 *         the open (0.36 for a 1 ft grass bed, 0.46 for a 3 ft slash bed)
 *
 * Formula: WAF = 1.83 / ln((20 + 0.36*h) / (0.13*h)) (Andrews 2012, RMRS-GTR-266, eq. 1)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_waf_unsheltered(amrex::Real delta_ft) noexcept
{
    delta_ft = amrex::max(static_cast<amrex::Real>(delta_ft), static_cast<amrex::Real>(0.1));
    amrex::Real numerator = 20.0 + 0.36 * delta_ft;
    amrex::Real denominator = 0.13 * delta_ft;
    return 1.83 / std::log(numerator / denominator);
}

/**
 * @brief Compute sheltered Wind Adjustment Factor (FARSITE formula)
 *
 * @param[in] f_c      Canopy closure [fraction, 0–1]
 * @param[in] h_c_ft   Canopy height [ft]
 * @param[in] delta_ft Fuel bed depth [ft]
 * @return WAF
 *
 * Formula: WAF = 0.555 / (sqrt(f_c*h_c) * ln((20+0.36*h_c)/(0.13*h_c)))
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_waf_sheltered(amrex::Real f_c, amrex::Real h_c_ft,
                                   amrex::Real delta_ft) noexcept
{
    f_c = amrex::max(static_cast<amrex::Real>(f_c), static_cast<amrex::Real>(0.01));
    h_c_ft = amrex::max(static_cast<amrex::Real>(h_c_ft), static_cast<amrex::Real>(0.1));
    delta_ft = amrex::max(static_cast<amrex::Real>(delta_ft), static_cast<amrex::Real>(0.1));

    amrex::Real ln_term = std::log((20.0 + 0.36 * h_c_ft) / (0.13 * h_c_ft));
    amrex::Real denom = std::sqrt(f_c * h_c_ft) * ln_term;
    denom = amrex::max(static_cast<amrex::Real>(denom), static_cast<amrex::Real>(1.0e-6));

    return 0.555 / denom;
}

/**
 * @brief Apply canopy Beer-Lambert extinction to WAF
 *
 * @param[in] waf_open   WAF in open (f_c = 0)
 * @param[in] f_c        Canopy closure [0–1]
 * @param[in] alpha_c    Canopy extinction coefficient
 * @return WAF_sheltered
 *
 * Formula: WAF_sheltered = waf_open * exp(-alpha_c * f_c)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_waf_canopy_exp(amrex::Real waf_open, amrex::Real f_c,
                                    amrex::Real alpha_c) noexcept
{
    return waf_open * std::exp(-alpha_c * f_c);
}

/**
 * @brief Compute terrain curvature on fire grid
 *
 * Curvature κ = d²z/dx² + d²z/dy² (positive for convex peaks).
 * Used in FARSITE ridge/valley detection.
 *
 * @param[out] curvature       1-component MultiFab on fire grid
 * @param[in]  fire_slopes     2-component (dz/dx, dz/dy) on fire grid
 * @param[in]  geom_fire       Fire grid geometry
 */
void compute_terrain_curvature(
    amrex::MultiFab& curvature,
    const amrex::MultiFab& fire_slopes,
    const amrex::Geometry& geom_fire);

/**
 * @brief Apply FARSITE terrain wind corrections
 *
 * Modifies wind speed and direction based on terrain:
 * - Ridge speed-up: windward ridge faces get k_ridge factor
 * - Sheltering: lee-side slopes experience wind reduction
 * - Valley channeling: wind is deflected into valleys
 * - Directional deflection: wind can turn up to ±45° on slopes
 *
 * @param[in,out] fire_wind_eff  u, v on fire grid (modified in-place)
 * @param[in]     fire_slopes    dz/dx, dz/dy
 * @param[in]     fire_curvature Terrain curvature
 * @param[in]     k_ridge        Ridge speed-up factor (≈ 1.5)
 * @param[in]     k_shelter      Lee-side shelter factor (≈ 0.6)
 * @param[in]     k_valley       Valley channeling factor (≈ 0.8)
 * @param[in]     k_deflect      Wind deflection factor (≈ 0.3)
 *
 * Algorithm (per cell, as implemented in ERF_FireWindAdjust.cpp):
 * 1. Compute slope magnitude |∇z| and wind magnitude |U|
 * 2. Compute wind-upslope alignment: cos_upslope = (∇z · U) / (|∇z| |U|)
 * 3. Classify terrain (no curvature test for ridge and shelter, since a
 *    uniform planar slope has zero curvature):
 *    - Windward slope (|∇z| > 0.05, cos_upslope > 0): factor = k_ridge
 *    - Lee slope (|∇z| > 0.05, cos_upslope < -0.5): factor = k_shelter
 *    - Valley (otherwise, curv < -0.01): factor = k_valley + (1-k_valley)*(1-min(|∇z|,1))
 *    - Flat/other: factor = 1
 * 4. Clamp: factor = clamp(factor, 0.1, 3.0)
 * 5. Apply speed scaling: U *= factor
 * 6. For valleys (curv < -0.01 and |∇z| > 0.05): rotate the scaled wind by
 *    k_deflect * asin(sin θ), sin θ = (∇z × U_scaled)_z / (|∇z| |U_scaled|)
 *    the sine of the wind-slope angle, towards the slope aspect
 */
void apply_farsite_terrain_wind(
    amrex::MultiFab& fire_wind_eff,
    const amrex::MultiFab& fire_slopes,
    const amrex::MultiFab& fire_curvature,
    amrex::Real k_ridge, amrex::Real k_shelter,
    amrex::Real k_valley, amrex::Real k_deflect);

#endif
