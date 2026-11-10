/**
 * @file ERF_BurnCurve.H
 * @brief The burn curve of an ignited structure: t-squared growth, plateau, linear decay.
 *
 * Kept apart from ERF_StructureIgnition.H so that ERF_FireParams.H (which that
 * header includes) can call make_burn_curve for its start-up check.
 */
#ifndef ERF_BURN_CURVE_H
#define ERF_BURN_CURVE_H

#include <AMReX_REAL.H>
#include <AMReX_Math.H>
#include <AMReX_GpuQualifiers.H>

namespace erf_structure_ignition {

/**
 * @brief Heat release per unit floor area of a burning structure [W/m2].
 *
 * Growth q = q_peak (tau / t_growth)^2 until the peak, a plateau at q_peak
 * until the share (1 - decay_fraction) of the fuel load has been released, then
 * a linear decay to zero that releases the remaining decay_fraction. Zero before
 * ignition and after burnout.
 */
struct BurnCurve
{
    amrex::Real q_peak    = 0.0;   ///< Peak flux [W/m2]
    amrex::Real t_growth  = 0.0;   ///< Time from ignition to the peak [s]
    amrex::Real t_plateau = 0.0;   ///< Duration at the peak [s]
    amrex::Real t_decay   = 0.0;   ///< Duration of the linear decay [s]

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real duration () const noexcept { return t_growth + t_plateau + t_decay; }

    /// Flux at time tau since ignition [W/m2]; non-negative.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real flux_at (amrex::Real tau) const noexcept
    {
        if (tau < 0.0) { return 0.0; }
        // ratio from a floored growth time before the test (t_growth = 0 is a
        // legal instant-growth curve; the unselected tau/0 is speculated)
        const amrex::Real r = tau / amrex::max(t_growth, amrex::Real(1.0e-30));
        if (tau < t_growth) {
            return q_peak * r * r;
        }
        tau -= t_growth;
        if (tau < t_plateau) { return q_peak; }
        tau -= t_plateau;
        // ratio from a floored decay time before the test (decay_fraction = 0 gives
        // t_decay = 0, a legal curve; the unselected tau/0 is speculated)
        const amrex::Real d = tau / amrex::max(t_decay, amrex::Real(1.0e-30));
        if (tau < t_decay) { return q_peak * (1.0 - d); }
        return 0.0;
    }
};

/**
 * @brief Build the burn curve from the peak flux, the fuel load per unit area
 *        the time to the peak and the share of the load the decay releases.
 *
 * The growth phase releases q_peak t_growth / 3, the decay decay_fraction of the
 * load; the plateau takes what is left. Returns false when an input is out of
 * range (decay_fraction outside [0, 1)) or when the growth phase alone would
 * release more than (1 - decay_fraction) of the load, when no such curve exists.
 * FireParams::FireParams() calls this at start-up for the same decision, so
 * the start-up check and the curve built at the first ignition cannot
 * disagree (not even by one ulp).
 */
inline bool make_burn_curve (amrex::Real q_peak, amrex::Real fuel_load_J_m2,
                             amrex::Real t_growth, amrex::Real decay_fraction, BurnCurve& c)
{
    if (q_peak <= 0.0 || fuel_load_J_m2 <= 0.0 || t_growth < 0.0) { return false; }
    if (!(decay_fraction >= 0.0 && decay_fraction < 1.0)) { return false; }
    const amrex::Real e_growth = q_peak * t_growth / 3.0;
    const amrex::Real e_decay  = decay_fraction * fuel_load_J_m2;
    if (e_growth > fuel_load_J_m2 - e_decay) { return false; }
    c.q_peak    = q_peak;
    c.t_growth  = t_growth;
    c.t_plateau = (fuel_load_J_m2 - e_decay - e_growth) / q_peak;
    c.t_decay   = 2.0 * e_decay / q_peak;
    return true;
}

} // namespace erf_structure_ignition

#endif // ERF_BURN_CURVE_H
