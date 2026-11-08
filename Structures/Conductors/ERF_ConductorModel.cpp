/**
 * @file ERF_ConductorModel.cpp
 * @brief The conductor lines and their towers as a structure model on ERF's structure coupling
 *
 * erf.structures.model = conductors selects it. It wraps Conductors (one MoorDyn system per line,
 * the towers, the gusts, the drag on the flow) behind ERF's StructureModel interface
 * (Source/MovingBodies/Coupling/ERF_StructureModel.H), so FLARE adds the lines to a run without
 * editing ERF: what the fork did in ERF's own files (the constructor, InitData, Advance, the slow
 * RHS, the checkpoint and the plotfile) happens in the calls below. The inputs are read under
 * conductors.* (Conductors::create).
 */
#include <ERF_StructureModelRegistry.H>

#include "ERF_Conductors.H"

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <memory>
#include <string>
#include <vector>

using namespace amrex;

namespace {

const char* const k_plot_names[3] = {"conductor_fx", "conductor_fy", "conductor_fz"};

class ConductorModel : public StructureModel
{
public:
    explicit ConductorModel (int max_level)
        : m_c(Conductors::create(max_level))
    {
        if (!m_c) {
            Abort("erf.structures.model = conductors needs the lines: list them in conductors.lines");
        }
        {
            // ERF hands a structure model no mesh heights or cell volumes on a mesh of constant dz, so an
            // embedded boundary's offset and volume fractions would not reach the lines' sampling and spreading
            std::string terrain;
            ParmParse("erf").query("terrain_type", terrain);
            if (amrex::toLower(terrain) == "eb") {
                Abort("erf.structures.model = conductors does not support erf.terrain_type = EB; use a "
                      "terrain-following mesh or ImmersedForcing for the ground under the lines");
            }
        }
    }

    const char* name () const override { return "conductors"; }

    int level (int finest_level) const override
    {
        m_c->require_anchor_level(finest_level);
        return m_c->anchor_level();
    }

    std::vector<RealBox> support () const override { return m_c->support(); }

    // ERF builds the surface only with erf.terrain_type = ImmersedForcing, where the lines stand on it
    bool wants_immersed_terrain_surface () const override { return true; }

    void initialize (const StructureInit& init) override
    {
        if (!init.unit_map_factors) {
            Abort("erf.structures.model = conductors: the lines sample and spread with the actuator core, "
                  "which needs every map factor of the anchor level " + std::to_string(init.lev) + " to be 1");
        }
        // the anchor level's closure: gusts on the lines come from the k-equation RANS's k
        m_c->set_closure(init.keqn_rans, init.Cmu0);
        if (!m_c->ground_set()) {
            if (init.immersed_terrain_surface != nullptr) {
                // an immersed terrain stands on a flat mesh: the lines take the terrain's height from the
                // surface the immersed boundary was built from, at the nodes of the anchor level
                m_c->set_ground_surface(*init.immersed_terrain_surface, init.geom);
            }
            // places the attachment points, creates the lines, and on a restart continues them from the
            // checkpoint; with drag_on_flow the force on the air is then spread again as it was there
            m_c->set_ground(init.z_phys_nd, init.geom, init.restart_chkfile);
            m_c->restore_sources(init.lev, init.xvel, init.z_phys_nd, init.detJ_cc, init.geom);
        }
    }

    void advance (const StructureStep& step) override
    {
        if (!step.unit_map_factors) {
            Abort("erf.structures.model = conductors: a map factor of the anchor level " + std::to_string(step.lev)
                  + " is not 1 after a regrid; the actuator core needs map factors of 1");
        }
        m_c->advance(step.lev, step.time, static_cast<double>(step.dt), step.xvel, step.yvel, step.zvel,
                     step.z_phys_nd, step.detJ_cc, step.geom, &step.cons);
    }

    void add_momentum_sources (int lev, MultiFab& xmom_src, MultiFab& ymom_src, MultiFab& zmom_src) const override
    {
        m_c->add_momentum_sources(lev, xmom_src, ymom_src, zmom_src);
    }

    void write_checkpoint (const std::string& checkpointname) const override
    {
        m_c->write_checkpoint(checkpointname);
    }

    bool supports_restart () const override { return true; }

    // the force on the air (N/m^3), only when drag_on_flow puts it into the flow
    std::vector<std::string> plot_var_names () const override
    {
        if (!m_c->drag_on_flow()) { return {}; }
        return {k_plot_names[0], k_plot_names[1], k_plot_names[2]};
    }

    void fill_plot_var (const std::string& name, int lev, MultiFab& dst, int comp) const override
    {
        int c = 0;
        while (c < 3 && name != k_plot_names[c]) { ++c; }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(c < 3, "ConductorModel::fill_plot_var: not a conductor plot variable");
        MultiFab src(dst.boxArray(), dst.DistributionMap(), 3, 0);
        m_c->cell_sources(lev, src, 0);   // zero on every level but the anchor level
        MultiFab::Copy(dst, src, c, comp, 1, 0);
    }

private:
    std::unique_ptr<Conductors> m_c;
};

} // namespace

ERF_REGISTER_STRUCTURE_MODEL(ConductorModel, "conductors")
