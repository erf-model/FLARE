/**
 * @file ERF_FlareDustModel.cpp
 * @brief FLARE's dust layer behind ERF's dust coupling (ERF_DustModel.H)
 *
 * DustLayer does the work; this class makes it a DustModel. It holds what
 * ERF itself held for the layer until the layer left ERF's source tree: the
 * layer, the fire-dust coupling with its erf.fire_dust_* inputs and checks,
 * the fire fields copied onto the dust grid each step, and the layer's
 * checkpoint fields. The calls are made at the points ERF made them, in the
 * same order, so results do not change.
 */

#include <ERF.H>
#include <ERF_DustModel.H>
#include <ERF_DustModelRegistry.H>
#include <ERF_DustLayer.H>
#include <ERF_DustParams.H>
#include <ERF_DustAtmCoupling.H>
#ifdef ERF_USE_FIRE
#include <ERF_DustFireLofting.H>
#include <ERF_FireDustCoupling.H>
#include <ERF_FireEngine.H>
#include <ERF_FireGrid.H>
#endif

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include <cmath>
#include <memory>
#include <string>

using namespace amrex;

namespace {

class FlareDustModel : public DustModel
{
public:
    const char* name () const override { return "flare"; }

    void initialize (const DustInit& init) override
    {
        DustParams dust_params;  // reads all erf.dust.* from ParmParse
        m_layer = std::make_unique<DustLayer>();
        m_layer->initialize(init.erf, init.surface_layer, init.z_phys_nd, dust_params);
        if (init.restart_chkfile.empty()) {
            m_layer->remove_outputs_for_fresh_start();   // stale CSVs of an earlier run are not continued
            m_layer->write_diag_header();
        }
#ifdef ERF_USE_FIRE
        if (init.fire) { initialize_fire_dust_coupling(*init.fire, dust_params); }
#endif

        // On a restart the fields initialize() just built from the inputs are
        // replaced by the checkpointed ones, and the layer's counters resume.
        if (!init.restart_chkfile.empty()) {
            read_checkpoint(init.restart_chkfile, init.step, init.time);
        } else {
            amrex::Print() << "[DUST] Zeroed dust_scalar_comp=" << m_layer->get_dust_scalar_comp()
                           << " in vars_new/old after all initialization.\n";
        }
    }

    void add_sources (int lev, MultiFab& cc_src, const MultiFab& state,
                      const MultiFab& detJ_cc, const Geometry& geom, Real dt) override
    {
        if (lev != 0) { return; }
        // the cell thickness is detJ dz (detJ_cc is 1 on a flat grid and is
        // always allocated); the centre-to-centre spacing of z_phys_cc that
        // the kernels used until October 2026 created mass on stretched grids
        m_layer->apply_to_cc_source(cc_src, detJ_cc, geom);
        m_layer->apply_settling_to_cc_source(cc_src, state, detJ_cc, geom, dt);
        m_layer->apply_deposition_bc(cc_src, state, detJ_cc, geom, dt);
    }

    void advance (const DustStep& step) override
    {
#ifdef ERF_USE_FIRE
        if (step.fire && m_fire_dust_coupling.enabled) { copy_fire_to_dust_grid(*step.fire); }
#endif
        const int nz = step.geom.Domain().length(2);
        // the return fields of the step just taken (the surface dust density,
        // the surface latent flux), once; the slow-RHS hook extracted them at
        // every RK stage and only the last survived
        m_layer->extract_atm_return_fields(step.cons, step.surface_q1_flux, step.geom);
        m_layer->advance(step.dt, m_layer->get_params(), step.surface_layer,
                         &step.xvel, &step.yvel, &step.zvel, step.z_phys_cc, &step.geom, nz,
                         &step.cons);
    }

    void write_output (int step, double time, bool is_final) override
    {
        m_layer->write_output(step, time, is_final);
    }

    void write_checkpoint (const std::string& checkpointname) const override
    {
        // Everything on the dust grid that carries from one step to the next (see
        // DustLayer::checkpoint_fields), the counters behind it, and the
        // super-particles.
        for (auto& nf : m_layer->checkpoint_fields()) {
            VisMF::Write(*nf.second, MultiFabFileFullPrefix(0, checkpointname, "Level_", nf.first));
        }
        m_layer->write_checkpoint_state(checkpointname);
#ifdef ERF_USE_PARTICLES
        m_layer->checkpoint_particles(checkpointname);
#endif
        if (ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "[DUST] Dust state written to checkpoint " << checkpointname << "\n";
        }
    }

private:
    /**
     * Restore the dust layer from the checkpoint, after initialize() has
     * allocated the fields from the inputs and rasters, so anything the
     * checkpoint lacks (an older checkpoint, or a field this build does not
     * have) keeps its initialized value.
     */
    void read_checkpoint (const std::string& restart_chkfile, int step, double time)
    {
        const std::string probe = restart_chkfile + "/Level_0/DustDepositionRate_H";
        if (!amrex::FileExists(probe)) {
            amrex::Print() << "[DUST] No dust state found in checkpoint; starting the dust"
                           << " layer from its initial values at the atmosphere's step and time.\n";
            // the dust clock (blast windows, road start times, PHREEQC and MSHA
            // intervals) follows the atmosphere's, not zero
            m_layer->read_checkpoint_state(restart_chkfile, step, time);
            // a fresh dust leg: stale CSVs of an earlier dust attempt are not continued
            m_layer->remove_outputs_for_fresh_start();
            m_layer->write_diag_header();
            return;
        }
        amrex::Print() << "[DUST] Restoring dust state from checkpoint " << restart_chkfile << "\n";
        // the layout aborts first: VisMF::Read dies on a mismatched BoxArray or
        // component count with no mention of dust
        m_layer->check_checkpoint_layout(restart_chkfile);

        // The ghost cells come back from the file as they were, with no boundary
        // fill afterwards. No dust kernel reads a ghost cell of a checkpointed
        // field (the only stencils, the terrain slopes and the FARSITE wind
        // factors, work on fields built at start-up); if one is added, FillBoundary
        // it after this read and in advance() alike, or the restarted run differs
        // from the uninterrupted one at box edges.
        for (auto& nf : m_layer->checkpoint_fields()) {
            const std::string header = restart_chkfile + "/Level_0/" + nf.first + "_H";
            if (!amrex::FileExists(header)) {
                amrex::Print() << "[DUST] Checkpoint has no " << nf.first
                               << "; keeping the initialized values.\n";
                continue;
            }
            VisMF::Read(*nf.second, amrex::MultiFabFileFullPrefix(0, restart_chkfile, "Level_", nf.first));
        }
        m_layer->read_checkpoint_state(restart_chkfile, step, time);
        m_layer->trim_outputs_after_restart(step);   // rows past the restart step go
#ifdef ERF_USE_PARTICLES
        m_layer->restart_particles(restart_chkfile);
#endif
    }

#ifdef ERF_USE_FIRE
    /// The erf.fire_dust_* inputs, their checks, and the coupling handed to the layer
    void initialize_fire_dust_coupling (const DustFire& fire, const DustParams& dust_params)
    {
        amrex::ParmParse pp("erf");
        pp.query("fire_dust_coupling",        m_fire_dust_coupling.enabled);
        pp.query("fire_dust_crust_reduction", m_fire_dust_coupling.post_fire_crust_reduction);
        pp.query("fire_dust_wind_to_dust",   m_fire_dust_coupling.fire_wind_to_dust);
        pp.query("fire_dust_wind_z0",        m_fire_dust_coupling.fire_wind_z0);
        // the datum of the wind handed over is the height the fire samples
        // it at (erf.fire.wind_ref_ht, or 6.096 m with
        // erf.fire.use_per_fuel_wind_ht): the dust-side height follows
        // it, and a different value aborts
        {
            amrex::Real wind_ref_ht = fire.engine.reference_wind_height();
            m_fire_dust_coupling.fire_wind_zref = wind_ref_ht;
            // an engine that hands over no reference wind cannot feed the dust u*
            if (m_fire_dust_coupling.enabled && m_fire_dust_coupling.fire_wind_to_dust &&
                fire.engine.reference_wind() == nullptr) {
                amrex::Print() << "[FIRE-DUST] WARNING: the fire engine (" << fire.engine.name()
                               << ") hands over no reference wind; erf.fire_dust_wind_to_dust is switched off\n";
                m_fire_dust_coupling.fire_wind_to_dust = false;
            }
            if (m_fire_dust_coupling.enabled && m_fire_dust_coupling.fire_wind_to_dust) {
            amrex::Real zref_given = wind_ref_ht;
            if (pp.query("fire_dust_wind_zref", zref_given) &&
                std::abs(zref_given - wind_ref_ht) > 1.0e-6 * wind_ref_ht) {
                amrex::Abort("[FIRE-DUST] erf.fire_dust_wind_zref = " + std::to_string(zref_given)
                             + " m must match " + std::to_string(wind_ref_ht)
                             + " m, the height of the wind the fire hands over (fire_wind_ref: erf.fire.wind_ref_ht,"
                             " or 6.096 m with erf.fire.use_per_fuel_wind_ht); leave it unset to follow it");
            }
            amrex::ParmParse pp_most("erf.most");
            amrex::Real most_z0 = -1.0;
            if (pp_most.query("z0", most_z0) && most_z0 > 0.0 &&
                std::abs(most_z0 - m_fire_dust_coupling.fire_wind_z0) > 1.0e-6 * most_z0) {
                amrex::Print() << "[FIRE-DUST] WARNING: erf.fire_dust_wind_z0 = " << m_fire_dust_coupling.fire_wind_z0
                               << " m differs from erf.most.z0 = " << most_z0
                               << " m; inside the fire perimeter the fire wind u* (a neutral log law on"
                                  " fire_dust_wind_z0) is compared with the surface layer's u*\n";
            }
            } else {   // the keys are read only when the fire wind feeds the dust
                for (const char* k : {"fire_dust_wind_zref", "fire_dust_wind_z0"}) {
                    if (pp.contains(k)) {
                        amrex::Print() << "[FIRE-DUST] WARNING: erf." << k
                                       << " is read only with erf.fire_dust_coupling = true and"
                                          " erf.fire_dust_wind_to_dust = true\n";
                    }
                }
            }
        }
        pp.query("fire_dust_lofting_enabled",     m_fire_dust_coupling.fire_lofting_enabled);
        pp.query("fire_dust_lofting_k_loft",      m_fire_dust_coupling.lofting_k_loft);
        pp.query("fire_dust_lofting_Q_threshold", m_fire_dust_coupling.lofting_Q_threshold);
        pp.query("fire_dust_lofting_Q_ref",       m_fire_dust_coupling.lofting_Q_ref);
        // Ranges the coupling kernels divide by or take as fractions;
        // the keys are read only with erf.fire_dust_coupling = true.
        {
            auto& c = m_fire_dust_coupling;
            auto fail = [](const std::string& m) { amrex::Abort("[FIRE-DUST] " + m); };
            if (!(c.post_fire_crust_reduction >= 0.0 && c.post_fire_crust_reduction <= 1.0)) {
                fail("erf.fire_dust_crust_reduction must be in [0, 1], got " + std::to_string(c.post_fire_crust_reduction));
            }
            if (!(c.fire_wind_z0 > 0.0)) {
                fail("erf.fire_dust_wind_z0 must be > 0 m, got " + std::to_string(c.fire_wind_z0));
            }
            if (!(c.fire_wind_zref > c.fire_wind_z0)) {
                fail("erf.fire_dust_wind_zref (" + std::to_string(c.fire_wind_zref) + " m) must exceed"
                     " erf.fire_dust_wind_z0 (" + std::to_string(c.fire_wind_z0) + " m); the fire wind u* is"
                     " kappa |U| / ln(zref/z0)");
            }
            if (!(c.lofting_k_loft >= 0.0)) {
                fail("erf.fire_dust_lofting_k_loft must be >= 0, got " + std::to_string(c.lofting_k_loft));
            }
            if (!(c.lofting_Q_threshold >= 0.0)) {
                fail("erf.fire_dust_lofting_Q_threshold must be >= 0 W/m2, got " + std::to_string(c.lofting_Q_threshold));
            }
            // (erf.fire_dust_lofting_Q_ref <= 0 is the documented off switch:
            // apply_fire_lofting_to_emission returns before any division)
            if (!c.enabled) {
                for (const char* k : {"fire_dust_crust_reduction", "fire_dust_wind_to_dust", "fire_dust_wind_z0",
                                      "fire_dust_wind_zref", "fire_dust_lofting_enabled", "fire_dust_lofting_k_loft",
                                      "fire_dust_lofting_Q_threshold", "fire_dust_lofting_Q_ref"}) {
                    if (pp.contains(k)) {
                        fail(std::string("erf.") + k + " needs erf.fire_dust_coupling = true");
                    }
                }
            }
        }
        m_fire_dust_coupling.debug       = dust_params.dust_debug;
        // DustLayer::advance applies the burned-area crust reduction from its
        // baseline each step; without this registration the reduction was
        // applied here every step on top of the previous one.
        m_layer->set_fire_dust_coupling(&m_fire_dust_coupling);

        // Validate grid_ratio matching for fire-dust coupling
        if (m_fire_dust_coupling.enabled) {
            int fire_gr = fire.grid.C;   // erf.fire.grid_ratio
            int dust_gr = dust_params.grid_ratio;
            std::string msg = "[FIRE-DUST] fire_dust_coupling=true requires "
                "dust.grid_ratio == fire.grid_ratio for proper sub-grid mapping. "
                "Set erf.dust.grid_ratio = " + std::to_string(fire_gr) +
                " (currently " + std::to_string(dust_gr) + ").";
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fire_gr == dust_gr, msg.c_str());
        }
    }

    /// The fire fields the dust step reads, copied onto the dust grid
    void copy_fire_to_dust_grid (const DustFire& fire)
    {
        const FireEngine& fire_engine = fire.engine;
        // Every fire field the dust step reads is copied onto the dust
        // BoxArray by index overlap, with the DUST periodicity: the dust
        // domain is the atmosphere's refined by the grid ratio, and with the
        // atmosphere's periodicity the copy treated dust cells one atmosphere
        // domain apart as periodic images and landed the fire wind and heat
        // in the wrong place on more than one rank. The grid ratios are
        // asserted equal at start-up, so the kernels read cell for cell.
        const amrex::Periodicity& per = m_layer->get_dust_geom().periodicity();
        const amrex::BoxArray&            dust_ba = m_layer->get_dust_ba();
        const amrex::DistributionMapping& dust_dm = m_layer->get_dust_dm();

        // Burned area for the crust reduction (applied inside DustLayer::advance
        // from the crust baseline). The three dust-grid scratch fields are
        // defined once (the dust BoxArray never changes) and refilled every step
        if (fire_engine.level_set() != nullptr) {
            if (!m_fire_phi_scratch_for_crust.ok()) {
                m_fire_phi_scratch_for_crust.define(dust_ba, dust_dm, 1, 0);
            }
            m_fire_phi_scratch_for_crust.setVal(1.0_rt);   // unburned where no fire cell lands
            m_fire_phi_scratch_for_crust.ParallelCopy(
                *fire_engine.level_set(), 0, 0, 1, 0, 0, per);
            m_fire_dust_coupling.fire_phi_scratch = &m_fire_phi_scratch_for_crust;
        } else {
            m_fire_dust_coupling.fire_phi_scratch = nullptr;
        }

        // Fire wind raises dust u*: the wind at FireParams::sampled_wind_ht()
        // (fire_wind_ref), the datum erf.fire_dust_wind_zref names
        if (m_fire_dust_coupling.fire_wind_to_dust &&
            fire_engine.reference_wind() != nullptr &&
            fire_engine.level_set() != nullptr) {
            amrex::MultiFab& fire_wind_scratch = m_fire_wind_scratch_for_ustar;   // (u, v), 1 ghost cell
            if (!fire_wind_scratch.ok()) { fire_wind_scratch.define(dust_ba, dust_dm, 2, 1); }
            fire_wind_scratch.setVal(0.0_rt);
            fire_wind_scratch.ParallelCopy(
                *fire_engine.reference_wind(), 0, 0, 2, 0, 1, per);
            fire_wind_scratch.FillBoundary(per);

            const int C = amrex::max(fire.grid.C / m_layer->get_grid_ratio(), 1);

            // Into the fire-u* field, which advance() maxes into u* after it
            // has filled u* from the surface layer.
            amrex::MultiFab* ustar_in = m_layer->get_ustar_fire_mut();
            if (ustar_in) {
                ustar_in->setVal(0.0_rt);
                m_fire_dust_coupling.apply_fire_wind_to_dust_ustar(
                    *ustar_in,
                    fire_wind_scratch,
                    m_fire_phi_scratch_for_crust,   // inside the fire perimeter only
                    m_layer->get_dust_geom(),
                    m_fire_dust_coupling.fire_wind_z0,
                    m_fire_dust_coupling.fire_wind_zref,
                    C);
                ustar_in->FillBoundary(per);
            }
        }

        // Fire heat flux for the convective lofting, applied inside
        // DustLayer::advance before the budget and the particles read the flux
        if (m_fire_dust_coupling.fire_lofting_enabled &&
            fire_engine.heat_release()) {
            if (!m_fire_heat_scratch_for_lofting.ok()) {
                m_fire_heat_scratch_for_lofting.define(dust_ba, dust_dm, 1, 1);
            }
            m_fire_heat_scratch_for_lofting.setVal(0.0_rt);
            m_fire_heat_scratch_for_lofting.ParallelCopy(
                *fire_engine.heat_release(), 0, 0, 1, 0, 1, per);
            m_fire_heat_scratch_for_lofting.FillBoundary(per);
            m_layer->set_fire_lofting(&m_fire_heat_scratch_for_lofting,
                                      m_fire_dust_coupling.lofting_k_loft,
                                      m_fire_dust_coupling.lofting_Q_threshold,
                                      m_fire_dust_coupling.lofting_Q_ref);
        }
    }

    FireDustCoupling m_fire_dust_coupling;
    amrex::MultiFab  m_fire_heat_scratch_for_lofting;   ///< fire heat flux on the dust grid, for the lofting inside DustLayer::advance
    amrex::MultiFab  m_fire_phi_scratch_for_crust;      ///< fire level set on the dust grid, for the burned-area crust reduction and the fire-wind mask
    amrex::MultiFab  m_fire_wind_scratch_for_ustar;     ///< fire reference wind (u, v) on the dust grid, for the fire-wind u*
#endif

    std::unique_ptr<DustLayer> m_layer;
};

} // namespace

ERF_REGISTER_DUST_MODEL(FlareDustModel, "flare")
