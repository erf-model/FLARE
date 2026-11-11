// Unit tests of Conductors, the manager ERF holds for the conductor lines, transformers and towers.
//
// AttachmentsArePlacedAboveTheTerrainUnderEachEnd: the attachment points stand at their height above
//     the k = 0 node plane of z_phys_nd (bilinear between the nodes), as ground.dat records.
// OnAUniformMeshTheHeightsAreAbsolute: without terrain the given heights are absolute.
// LinesStepOnTheAnchorLevelOnlyAndLogEveryStep: other levels return at once; one row per step; an anchor
//     level that does not exist (yet, or any more) stops the run, naming conductors.anchor_level.
// TheFlowIsSampledAtTheLinesCurrentPosition: without a prescribed velocity, MoorDyn gets ERF's
//     velocity at the line's current position, not where it hung.
// LinesRunWithTheFloatingPointTrapsOn: amrex.fpe_trap_invalid, _zero and _overflow on (in the inputs and on
//     the thread) no longer stop the run at create(); MoorDyn starts and steps under them, and they stay on.
// ClearanceIsTheHeightAboveTheTerrainUnderEachNode: the clearance against the terrain under each node.
// TheSeafloorHoldsTheTerrainNodesUnderTheLineAndOneSpanAround: with ground_contact the seafloor file holds the anchor
//     level's terrain nodes within one span length of the line (clipped to the domain), raised by the conductor's radius
//     and lowered by surface_offset, and
//     the MoorDyn input names it with kbot = 3e6 Pa/m and the conductor's critical cbot = 2 sqrt(kbot mu / d).
// ASlackSpanRestsOnTheTerrainWithGroundContactAndStopsTheRunWithout: a span sagging 58 m from ends 30 m above a ramp
//     rests on it (its axis one radius above the terrain, MoorDyn-C within 1 mm), its nodes on the ground counted and logged from the
//     first step; without contact the same span reaches into the terrain and the run stops.
// GroundContactRefusesAMoorDynStepTooLongForIt: a MoorDyn internal step past half the measured stability limit of the
//     contact is refused at set_ground, naming moordyn_dt and moordyn_cfl; so is a conductor that would sink more than a
//     quarter of its diameter into the ground at rest.
// GroundContactRecordsContinueAcrossARestart: the spans' first touch and most nodes on the ground come back from the
//     checkpoint, are not reported again, and continue as the uninterrupted run's do; a restart without ground_contact
//     keeps none, and one with it from a checkpoint written without it starts them at the restart.
// OnAMeshWithoutTerrainTheSpanRestsOnTheDomainsFloor: without terrain the seafloor is prob_lo[2] plus the radius and the
//     span rests there, inside the domain.
// TheContactRateIsTheFastestRootOfEachLineTypeOnTheGround: ground_rate is 2 omega for a critically damped or
//     under-damped type and 2 (a + sqrt(a^2 - omega^2)) for an over-damped one.
// AStringReachingBelowTheGroundIsRefused: placed, a tower's string may not end below the ground.
// TheGridsMustReachAsFarAsTheDragDoes: the coverage reach is a cell, and with drag_on_flow 3 epsilon cells, for the
//     epsilon of the node's line or tower type.
// TheDragNeedsTheGridsAsFarAsItReaches: with a prescribed wind, the spreading itself refuses a patch the drag outreaches.
// TheSpansDragStatisticIsAcrossTheSpan: drag_normal of a span along +y in a wind along +x is minus its
//     x drag (the horizontal drag along z x the span's direction), not its y drag.
// TheSpreadDragIntegratesToMinusTheDragOnTheLines: the momentum source integrates to minus the drag.
// WithoutDragOnFlowNothingIsPutIntoTheFlow: drag_on_flow off leaves the sources and the support empty.
// TheASCE74CheckTakesTheWireAtItsEffectiveHeightOnItsInclinedChord: asce74.csv's height is the mean attachment
//     height less a third of the still-air sag, and its swing that about the inclined chord against W cos beta;
//     asce74_wire_height = attachment and asce74_inclined_spans = false give the attachment height and W.
// ARestartCannotMoveALinesPointsOrChangeItsLength: the checkpoint records the placed points and lengths,
//     and a restart whose inputs move an end or change a length stops, naming the line.
// ARestartContinuesTheLinesTheirSourcesStatisticsAndLogs: a restart continues where the checkpoint left off.
// TrimmingALogKeepsTheHeaderAndTheRowsUpToTheCheckpoint: trim_log_after.
// ASectionIsPlacedOnTheTerrainAtEveryTowerAndLogsEachSpan: a section's towers and per-span logs.
// TheClosestApproachOfTwoLinesIsFlaggedAgainstTheFlashoverDistance: the exact separation and the flag.
// ARestartContinuesASectionItsStringsAndTheSeparationOfTheLines: a section's restart; an older checkpoint's
//     drag_y statistics continue on spans along +x; a restart that changes insulator_length stops.
// TransformersTakeThePullOfTheLinesEndingOnThemAndContinueAcrossARestart: transformer loads,
//     allowables, clearances, and their restart.
// AStringingTensionSetsTheLengthsFromTheChordsOnTheTerrain: strung lengths from the placed chords.
// LatticeTowersStandAtTheSuspensionPointsAndCarryTheWindsDrag: tower placement, member drag, the
//     line's pull and the footing checks.
// TheTowersDragGoesIntoTheFlowWithTheLinesAndSurvivesARestart: the towers' drag in the sources.
// EachLineAndTowerSpreadsItsDragWithItsOwnWidth: a line's nodes take conductors.<line>.epsilon, its towers'
//     drag nodes their type's conductors.<type>.epsilon and a line without its own conductors.epsilon: the
//     sources are those of spreading every node with that width, and differ from one width for all; a line's
//     coverage check takes its own reach.
// EachTowerIsCoveredWithItsTypesReach: the coverage checks of the sampled towers and of the spreading take each tower's
//     reach from its type's width and each line's from its own, and name the tower or line that outreaches the grids.
// MovingTowersSettleWhereTheirStiffnessBalancesTheWindAndTheLine: bending towers at rest.
// MovingTowersContinueAcrossARestart: bending towers' restart, coupling.dat's counts across a checkpoint
//     between two of its rows included.
// ACircuitHangsFromOneRowOfTowersEachLineAtItsOwnPoint: shared towers, points on the tower's base.
// ACircuitOnBendingTowersMovesEveryLinesPoint: bending shared towers move every line's point.
// AShortTautSpanOnBendingTowersStaysStable: the iterated coupling keeps a stiff span stable.
// AnImmersedTerrainPlacesEverythingOnItsSurface: an immersed terrain on a flat mesh; with ground_contact its seafloors
//     (the surface's branch) equal the fitted mesh's (the nodes' branch).
// SetGroundRefusesASurfaceOffsetBelowTheDomainTop: conductors.surface_offset must hold the domain.
// AnAttachmentOutsideTheDomainIsRefusedNamingItsKey: the abort names end_a, end_b or the tower.
// ANonFiniteCouplingPullIsRefusedNotConverged: coupling_converged on NaN pulls; coupling_diverged at 20 in a row.
// ARestartChecksTheTowerSwayAndTheSurfaceOffset: restart_mismatch, and tower_forms_mismatch (the tower forms); a
//     tower_forms_version 1 checkpoint (the rectangular cross-arm's nodes) cannot restart the outside layout.
// GustsComeFromTheRANSkAlongEachSpan: with gust_type = factor, per span the root-mean-square wind and normal wind
//     over its nodes and the mean k = (rho k)/rho, from stats_start, sampled where the nodes are at each step's start;
//     gusts.csv's columns for a line across the wind and one at 45 degrees to it; nothing before stats_start; the
//     set_closure and conserved-state aborts.
// GustsSwitchedOnAtARestartStartAfresh: a checkpoint written without gusts restarts with gust_type = factor, the
//     gust statistics starting at the restart.
// AnEventGustPassesEveryNodeAndTowerWithItsFront: with gust_type = event, the wind at every line node and tower drag
//     node is the flow's plus g sigma_u times the 1 - cos shape at the node's phase at the step's middle, sigma_u from
//     the mean k over the span's nodes or the tower's load-weighted one (the face area each drag node stands for
//     times the wind on it), which an unweighted mean over the tower's nodes misses; the origin defaults to the domain's centre;
//     gust_series.dat.
// RandomGustsFollowTheirProcessesAndACircuitSharesThem: with gust_type = random, every span's gust is
//     sigma_u sqrt(B) z, z the Ornstein-Uhlenbeck process stepped with T = L_u / U, the phases and the shield wire on
//     shared towers taking the owner's spans' processes; a string the mean of the spans beside it; each tower
//     sigma_u sqrt(B_t) z of its own, its k, U and the height of L_u the load-weighted (area x wind) means over its drag nodes;
//     gust_series.dat's columns.
// RandomGustsContinueAcrossARestart: a restart continues the random gusts where the checkpoint left them, and a
//     checkpoint written without gusts draws them afresh from the stationary law at the restart's first step.
// TheGustFactorRunsOnLinesWithLatticeTowers: gust_type = factor on a section with lattice towers, whose drag nodes
//     it does not sample, writes gusts.csv and no gust series.
// CapacitiesAreCheckedAfterEveryStepOfTheAnchorLevelOnly: with two levels and the lines on level 1, a level-0 call
//     records nothing and each component's first exceedance is stamped with the end of the first level-1 step; the
//     still-air conductor tension (the demand) matches the elastic catenary; failures.dat holds each once.
// ACapacityAboveTheDemandIsNeverExceeded: an RBS 1 % over the still-air tension: peak ratio 1/1.01, no row.
// FailureStopEndsTheRunAtTheFirstExceedance: the abort names the component and time after failures.dat is written.
// FailureStopNeedsACapacity: failure_stop without any capacity is refused at set_ground.
// CapacityRecordsContinueAcrossARestart: peaks and first exceedances come back from the checkpoint and are not reported
//     again (nor stop a run restarted with failure_stop), a capacity added at the restart starts there (the trimming of rows after the checkpoint is the CTest
//     Conductors_Capacities_Restart's).
// FailureBreakLetsGoAtTheHigherEndAndTheRunGoesOn: an inclined span whose conductor is over its strength at the first
//     step lets go at its higher end (the higher tension) at that step's end, once; breaks.dat and failures.dat hold one
//     row each, the run goes on, the span comes to hang from its lower end with most of it on the ground.
// FailureBreakNeedsGroundContactALineCapacityAndNoFailureStop: the three refusals, each naming its input.
// ASpanIsReleasedOnceInAStepAndADeadEndReleasesItsOwnEnd: the conductor's and both dead ends' exceedances in one step
//     release the span once; with the hardware's capacity alone, end a's exceedance (the first logged) releases end a.
// ReleasesContinueAcrossARestart: a checkpoint written at the release's step and one written later both restart
//     into the uninterrupted run's positions, with the release made again, logged once, also with failure_break off;
//     without ground_contact the restart is refused; a restart that first sets failure_break after the exceedance
//     releases nothing.
// StringsOverTheirLoadLetGoAndTheirTowersTakeNothing: a section on strings over two standing towers, the strings'
//     rated load below what they carry: both let go at the first step and the towers carry no line load.
// LinesOnMovingTowersLetGoAtTheEndOfTheLevelsStep: on towers that bend the exceedances are found after the first
//     coupling step inside the level's step; the releases are made at that step's end, one per span.
// AStandingTowersFootingsAndAttachmentsAreChecked: in still air each standing tower's legs carry (W - F_z) / 4, which
//     exceeds a lower allowable_compression, no uplift; the hardware at a tower carries the line's pull there.
// AStandingTowerIsCheckedUnderTheLinesPullAtTheStepsEnd: a tower that stands still is checked under the lines' pull at
//     the end of the step, the load it carries over the next step, not the start-of-step load towers.dat reports.
// AMovingTowerIsCheckedAfterEveryCouplingStep: lines on towers that bend are checked with their towers at the end of
//     every coupling step, t = h = dt / N, not only at the end of the ERF step.
// AStepsFirstExceedancesAreLoggedInTimeOrderAcrossGroups: with two groups on bending towers, the step's first
//     exceedances are written in time order and failure_stop names the earliest, not the first group's.
// ACircuitsStringsHardwareAndSharedFootingsAreChecked: a circuit on shared bending towers: strings, hardware and the
//     shared towers' footings (once each, all four lines' load) are first exceeded at t = h, their demands the end state.
// AFrameTowersMembersAreCheckedAgainstTheirStrength: generated frame towers: one component per member, its ratio the
//     member's ASCE 10-15 utilisation.
// ALineTakesAnotherLinesGustsAndAChangedSharingDrawsAfresh: a line naming another in gust_with takes its random
//     gust, each along its span's mean horizontal wind in a wind that veers with the height; a restart whose lines
//     share differently, with as many processes, draws them afresh.

#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_RealBox.H>

#include <gtest/gtest.h>

#include "ERF_GTestTempDir.H"
#include "ERF_ActuatorSpreading.H"
#include "ERF_Conductors.H"
#include "ERF_Constants.H"
#include "ERF_Gusts.H"
#include "ERF_IndexDefines.H"
#include "ERF_GTestThrowOnAbort.H"
#include "ERF_MoorDynSystem.H"

namespace {

// a scratch root drawn once per test process (ERF_GTestTempDir.H): the fixed names below it are this
// process's alone, so ctest -j and the shuffled rerun never share them; removed when the process exits
const std::filesystem::path& gtest_scratch_root ()
{
    struct Root {
        std::filesystem::path p;
        ~Root () { std::error_code ec; std::filesystem::remove_all(p, ec); }
    };
    static const Root root{[] {
        const std::filesystem::path p = erf_gtest_temp_path("erf_gtest_conductors");
        std::filesystem::create_directories(p);
        return p;
    }()};
    return root.p;
}

using amrex::Real;

// the roundoff of the sampling, the terrain read and the source sums: about 1e-7 relative in single precision
constexpr Real roundoff = (std::is_same<Real, float>::value) ? Real(1.0e-5) : Real(1.0e-9);
// positions and heights of a few hundred metres carry a Real's spacing there: 3e-5 m at 500 m in single precision
constexpr Real postol = (std::is_same<Real, float>::value) ? Real(1.0e-4) : Real(1.0e-6);
// the same quantity computed twice the same way: a few units of the last place
constexpr Real tight = (std::is_same<Real, float>::value) ? Real(1.0e-6) : Real(1.0e-12);

// a 1200 m x 1000 m x 400 m box of 12 x 10 x 8 cells; terrain-following nodes over a ramp
// h = slope_x x + slope_y y, so the bilinear surface under any point is exact
struct Mesh {
    int nx = 12, ny = 10, nz = 8;
    Real Lx = 1200.0, Ly = 1000.0, H = 400.0;
    Real slope_x = 0.05, slope_y = 0.02;
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    std::unique_ptr<amrex::MultiFab> znd;
    amrex::MultiFab u, v, w;

    Real h (Real x, Real y) const { return slope_x * x + slope_y * y; }

    // the crosswind v = v0 + vy (y - 500) + vz z, u = w = 0: linear, so the sampler reproduces it exactly
    Real v0 = 10.0, vy = 0.02, vz = 0.05;
    Real v_at (Real y, Real z) const { return v0 + vy * (y - 500.0) + vz * z; }
    void fill_crosswind ()
    {
        const Real dx = Lx / nx, dy = Ly / ny, dz = H / nz;
        u.setVal(0.0); w.setVal(0.0);
        for (amrex::MFIter mfi(v); mfi.isValid(); ++mfi) {
            auto va = v.array(mfi);
            amrex::LoopOnCpu(mfi.growntilebox(), [&](int i, int j, int k) {
                amrex::ignore_unused(i, dx);
                va(i,j,k) = v_at(j * dy, (k + 0.5) * dz);
            });
        }
    }

    explicit Mesh (bool terrain)
    {
        const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(nx-1, ny-1, nz-1));
        const amrex::RealBox rb({AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(Lx, Ly, H)});
        const std::array<int,3> periodic{{0, 0, 0}};
        geom = amrex::Geometry(domain, &rb, 0, periodic.data());
        ba = amrex::BoxArray(domain);
        ba.maxSize(amrex::IntVect(4, 5, 1024));   // several boxes: the terrain read has one owner per point
        dm = amrex::DistributionMapping(ba);
        u.define(amrex::convert(ba, amrex::IntVect(1,0,0)), dm, 1, 1); u.setVal(0.0);
        v.define(amrex::convert(ba, amrex::IntVect(0,1,0)), dm, 1, 1); v.setVal(0.0);
        w.define(amrex::convert(ba, amrex::IntVect(0,0,1)), dm, 1, 1); w.setVal(0.0);
        if (terrain) {
            znd = std::make_unique<amrex::MultiFab>(amrex::convert(ba, amrex::IntVect(1,1,1)), dm, 1, 1);
            const Real dx = Lx / nx, dy = Ly / ny;
            for (amrex::MFIter mfi(*znd); mfi.isValid(); ++mfi) {
                auto za = znd->array(mfi);
                amrex::LoopOnCpu(mfi.growntilebox(), [&](int i, int j, int k) {
                    const Real hs = h(i * dx, j * dy);
                    za(i,j,k) = hs + (H - hs) * Real(k) / nz;
                });
            }
        }
    }
};

std::string scratch (const std::string& tag)
{
    const auto dir = gtest_scratch_root() / ("erf_gtest_conductors_" + tag);
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    return dir.string();
}

// one span, 300 m along x at y = 500, ends 30 m above the local surface
void set_inputs (const std::string& dir, const std::string& name, bool prescribed = true,
                 const std::array<Real,3>& a = {{300.0, 500.0, 30.0}}, const std::array<Real,3>& b = {{600.0, 500.0, 30.0}})
{
    amrex::ParmParse pp("conductors");
    pp.add("lines", name);
    pp.add("diagnostics_dir", dir);
    pp.add("air_density", 1.2);
    // ParmParse is global to the test binary: drop what another test left before setting this one's
    for (const char* key : {"prescribed_velocity", "drag_on_flow", "epsilon", "node_output_int", "stats_start", "flashover_distance",
                            "ground_contact", "moordyn_dt"}) { pp.remove(key); }
    if (prescribed) { pp.addarr("prescribed_velocity", std::vector<Real>{0.0, 10.0, 0.0}); }
    amrex::ParmParse ps("conductors." + name);
    ps.addarr("end_a", std::vector<Real>{a[0], a[1], a[2]});
    ps.addarr("end_b", std::vector<Real>{b[0], b[1], b[2]});
    ps.add("length", 301.5);
    ps.add("diameter", 0.0281);
    ps.add("mass_per_length", 1.628);
    ps.add("axial_stiffness", 3.0e7);
}

} // namespace

TEST(Conductors, AttachmentsArePlacedAboveTheTerrainUnderEachEnd)
{
    const std::string dir = scratch("terrain");
    set_inputs(dir, "Tterrain");
    Mesh m(true);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(m.znd.get(), m.geom);
    ASSERT_EQ(c->lines().size(), 1u);
    const auto& span = *c->lines().front();
    const unsigned last = span.num_nodes() - 1;
    const auto a = span.node_position(0);
    const auto b = span.node_position(last);
    // the ramp rises 15 m from x = 300 to x = 600: the two ends sit at different absolute heights
    EXPECT_NEAR(a[2], m.h(300.0, 500.0) + 30.0, postol);
    EXPECT_NEAR(b[2], m.h(600.0, 500.0) + 30.0, postol);
    EXPECT_NEAR(a[0], 300.0, postol);
    EXPECT_NEAR(b[0], 600.0, postol);
    EXPECT_GT(b[2] - a[2], 10.0);
    // ground.dat records the surface and the absolute height of each end
    std::ifstream g(dir + "/ground.dat");
    ASSERT_TRUE(g.good());
    std::string header, name, end;
    Real x, y, ground, z;
    std::getline(g, header);
    ASSERT_TRUE(static_cast<bool>(g >> name >> end >> x >> y >> ground >> z));
    EXPECT_EQ(name, "Tterrain"); EXPECT_EQ(end, "a");
    EXPECT_NEAR(ground, m.h(300.0, 500.0), postol);
    EXPECT_NEAR(z, m.h(300.0, 500.0) + 30.0, postol);
    ASSERT_TRUE(static_cast<bool>(g >> name >> end >> x >> y >> ground >> z));
    EXPECT_EQ(end, "b");
    EXPECT_NEAR(ground, m.h(600.0, 500.0), postol);
}

TEST(Conductors, OnAUniformMeshTheHeightsAreAbsolute)
{
    const std::string dir = scratch("flat");
    set_inputs(dir, "Tflat");
    Mesh m(false);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(nullptr, m.geom);
    const auto& span = *c->lines().front();
    EXPECT_NEAR(span.node_position(0)[2], 30.0, postol);
    EXPECT_NEAR(span.node_position(span.num_nodes() - 1)[2], 30.0, postol);
}

TEST(Conductors, LinesStepOnTheAnchorLevelOnlyAndLogEveryStep)
{
    const std::string dir = scratch("advance");
    set_inputs(dir, "Tadv");
    Mesh m(false);
    auto c = Conductors::create(1);   // two levels: the anchor is the finest, level 1
    ASSERT_TRUE(c);
    EXPECT_EQ(c->anchor_level(), 1);
    {
        // amr.max_level = 1 with no level 1 yet: the lines would stand on flat ground and never step
        const std::string msg = erf_gtest::abort_message([&] { c->require_anchor_level(0); });
        EXPECT_NE(msg.find("conductors.anchor_level"), std::string::npos) << msg;
        EXPECT_NE(msg.find("finest level is 0"), std::string::npos) << msg;
        EXPECT_TRUE(erf_gtest::abort_message([&] { c->require_anchor_level(1); }).empty());
    }
    c->set_ground(nullptr, m.geom);
    const auto& span = *c->lines().front();
    const Real off0 = span.mid_offset();
    c->advance(0, 0.0, 0.1, m.u, m.v, m.w, nullptr, nullptr, m.geom);   // level 0: nothing happens
    EXPECT_DOUBLE_EQ(span.mid_offset(), off0);
    EXPECT_FALSE(std::filesystem::exists(dir + "/Tadv.dat"));
    for (int s = 0; s < 3; ++s) { c->advance(1, 0.1 * s, 0.1, m.u, m.v, m.w, nullptr, nullptr, m.geom); }
    EXPECT_GT(span.mid_offset(), off0) << "the prescribed +y wind must move the span towards +y";
    std::ifstream f(dir + "/Tadv.dat");
    ASSERT_TRUE(f.good());
    std::string line;
    int rows = 0;
    while (std::getline(f, line)) { if (!line.empty() && line.rfind("time", 0) != 0) { ++rows; } }
    EXPECT_EQ(rows, 4) << "the initial row and one per step";
}

TEST(Conductors, AStringReachingBelowTheGroundIsRefused)
{
    // on a flat mesh whose domain starts 5 m up, z is the absolute height: a tower at z = 6 m with a 4 m
    // string passes the input check (6 > 4), but the string's bottom stands 3 m below the ground at 5 m
    const std::string dir = scratch("string_ground");
    set_inputs(dir, "Tsg", true, {{300.0, 500.0, 10.0}}, {{600.0, 500.0, 10.0}});
    amrex::ParmParse ps("conductors.Tsg");
    ps.addarr("towers", std::vector<Real>{450.0, 500.0, 6.0});
    ps.addarr("length", std::vector<Real>{151.0, 151.0});
    ps.add("insulator_length", 4.0);
    ps.add("insulator_mass", 50.0);
    Mesh m(false);
    const amrex::RealBox rb({AMREX_D_DECL(Real(0.0), Real(0.0), Real(5.0))}, {AMREX_D_DECL(m.Lx, m.Ly, m.H + Real(5.0))});
    const std::array<int,3> periodic{{0, 0, 0}};
    const amrex::Geometry raised(m.geom.Domain(), &rb, 0, periodic.data());
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    const std::string msg = erf_gtest::abort_message([&] { c->set_ground(nullptr, raised); });
    EXPECT_NE(msg.find("conductors.Tsg.towers"), std::string::npos) << msg;
    EXPECT_NE(msg.find("the insulator string's bottom lies 3.000000 m below the ground"), std::string::npos) << msg;
}

TEST(Conductors, TheGridsMustReachAsFarAsTheDragDoes)
{
    // the sampler needs the cells around a node; with drag_on_flow the drag reaches 3 epsilon cells
    const std::string dir = scratch("reach");
    set_inputs(dir, "Treach");
    Mesh m(false);
    const Real cell = std::max(static_cast<Real>(m.geom.CellSize(0)), static_cast<Real>(m.geom.CellSize(1)));
    {
        auto c = Conductors::create(0);
        EXPECT_EQ(c->coverage_reach(m.geom, c->inputs().epsilon), cell) << "without drag_on_flow the sampler's cell only";
    }
    amrex::ParmParse pp("conductors");
    pp.add("drag_on_flow", true);
    pp.add("epsilon", 2.0);
    auto c = Conductors::create(0);
    EXPECT_EQ(c->coverage_reach(m.geom, c->inputs().epsilon), std::max(cell, Real(6.0) * static_cast<Real>(m.geom.CellSize(0))));
    // a line or tower type with a width of its own reaches 3 of its widths
    EXPECT_EQ(c->coverage_reach(m.geom, 0.5), std::max(cell, Real(1.5) * static_cast<Real>(m.geom.CellSize(0))));
    EXPECT_EQ(c->coverage_reach(m.geom, 3.0), std::max(cell, Real(9.0) * static_cast<Real>(m.geom.CellSize(0))));
}

TEST(Conductors, TheDragNeedsTheGridsAsFarAsItReaches)
{
    // a prescribed wind samples nothing, so only the spreading can see that the level's patch (x 0 to 1000 m)
    // stops 400 m past the span's end: beyond epsilon (200 m), inside the drag's reach of 3 epsilon cells (600 m)
    const std::string dir = scratch("drag_reach");
    set_inputs(dir, "Tdr");
    amrex::ParmParse pp("conductors");
    pp.add("drag_on_flow", true);
    pp.add("epsilon", 2.0);
    Mesh m(false);
    const amrex::BoxArray ba(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(9, m.ny - 1, m.nz - 1)));
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab u(amrex::convert(ba, amrex::IntVect(1,0,0)), dm, 1, 1), v(amrex::convert(ba, amrex::IntVect(0,1,0)), dm, 1, 1),
                    w(amrex::convert(ba, amrex::IntVect(0,0,1)), dm, 1, 1);
    u.setVal(0.0); v.setVal(0.0); w.setVal(0.0);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(nullptr, m.geom);
    std::string msg = erf_gtest::abort_message([&] { c->advance(0, 0.0, 0.2, u, v, w, nullptr, nullptr, m.geom); });
    EXPECT_NE(msg.find("a node of line Tdr, is not covered, with its drag's reach of 3 epsilon cells"), std::string::npos) << msg;
    struct Cleanup { ~Cleanup () { amrex::ParmParse("conductors.Tdr").remove("epsilon"); } } cleanup;
    // the line's own width, 1 cell, reaches 300 m past its end, inside the patch
    amrex::ParmParse("conductors.Tdr").add("epsilon", 1.0);
    auto d = Conductors::create(0);
    ASSERT_TRUE(d);
    d->set_ground(nullptr, m.geom);
    msg = erf_gtest::abort_message([&] { d->advance(0, 0.0, 0.2, u, v, w, nullptr, nullptr, m.geom); });
    EXPECT_TRUE(msg.empty()) << msg;
    amrex::ParmParse("conductors.Tdr").remove("epsilon");
    // and a width of 1.5 cells, 450 m, outreaches it again
    amrex::ParmParse("conductors.Tdr").add("epsilon", 1.5);
    auto e = Conductors::create(0);
    ASSERT_TRUE(e);
    e->set_ground(nullptr, m.geom);
    msg = erf_gtest::abort_message([&] { e->advance(0, 0.0, 0.2, u, v, w, nullptr, nullptr, m.geom); });
    EXPECT_NE(msg.find("a node of line Tdr, is not covered"), std::string::npos) << msg;
}

TEST(Conductors, TheSpansDragStatisticIsAcrossTheSpan)
{
    const std::string dir = scratch("dragnormal");
    set_inputs(dir, "Tdn", false, {{500.0, 300.0, 30.0}}, {{500.0, 600.0, 30.0}});
    amrex::ParmParse("conductors").addarr("prescribed_velocity", std::vector<Real>{10.0, 0.0, 0.0});
    Mesh m(false);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(nullptr, m.geom);
    c->advance(0, 0.0, 0.2, m.u, m.v, m.w, nullptr, nullptr, m.geom);
    const auto D = c->lines().front()->span_drag(0);
    ASSERT_GT(D[0], 1.0) << "the wind across the span drags it along +x";
    std::ifstream st(dir + "/Tdn_stats.csv");
    ASSERT_TRUE(st.good());
    std::string row;
    double mean = 0.0;
    bool found = false;
    while (std::getline(st, row)) {
        std::stringstream ss(row);
        std::string cell;
        std::vector<std::string> cells;
        while (std::getline(ss, cell, ',')) { cells.push_back(cell); }
        if (cells.size() >= 5 && cells[3] == "drag_normal") { mean = std::stod(cells[4]); found = true; }
    }
    ASSERT_TRUE(found);
    EXPECT_NEAR(mean, -static_cast<double>(D[0]), 1.0e-6 * std::abs(static_cast<double>(D[0]))) << "one sample: the drag across";
}

TEST(Conductors, TheFlowIsSampledAtTheLinesCurrentPosition)
{
    const std::string dir = scratch("sampled");
    set_inputs(dir, "Tsampled", false);
    amrex::ParmParse pp("conductors");
    ASSERT_FALSE(pp.contains("prescribed_velocity")) << "set_inputs must have removed the prescribed velocity";
    Mesh m(false);
    m.fill_crosswind();
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(nullptr, m.geom);
    const auto& span = *c->lines().front();
    Real moved = 0.0;
    for (int s = 0; s < 20; ++s) {
        const std::vector<Real> where = span.kinematics_points();   // where the line is when the wind is sampled
        c->advance(0, 0.25 * s, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom);
        for (unsigned p = 0; p < span.num_kinematics_points(); ++p) {
            // the flow at the line nodes; no wind at MoorDyn's fixed entries after them
            const auto uvw = span.wind_at_point(p);
            const Real expect = (p < span.num_nodes()) ? m.v_at(where[3*p+1], where[3*p+2]) : Real(0.0);
            ASSERT_NEAR(uvw[0], 0.0, roundoff) << "step " << s << " point " << p;
            ASSERT_NEAR(uvw[1], expect, roundoff * std::max(Real(1.0), expect)) << "step " << s << " point " << p;
            ASSERT_NEAR(uvw[2], 0.0, roundoff) << "step " << s << " point " << p;
        }
        moved = std::max(moved, span.mid_offset());
    }
    // the line really moved, so sampling at the initial positions would have handed a different wind
    EXPECT_GT(moved, 0.5) << "the span must blow out in the sampled crosswind";
    EXPECT_GT(m.vy * moved, 1.0e-3) << "the field must change measurably over the distance the line moved";
}

TEST(Conductors, LinesRunWithTheFloatingPointTrapsOn)
{
    // the traps as given in the inputs (create() refused them before MoorDyn's calls were masked) and as
    // AMReX sets them on the thread; both restored when the test ends, a failed assertion included
    struct TrapsOn {
        const std::vector<std::string> keys{"fpe_trap_invalid", "fpe_trap_zero", "fpe_trap_overflow"};
        std::vector<int> had, was;
        amrex::FPExcept prev;
        TrapsOn ()
        {
            amrex::ParmParse pa("amrex");
            for (const auto& k : keys) {
                bool b = false;
                had.push_back(pa.contains(k.c_str()) ? 1 : 0);
                pa.query(k.c_str(), b);
                was.push_back(b ? 1 : 0);
                pa.add(k.c_str(), true);
            }
            prev = amrex::enableFPExcept(amrex::FPExcept::all);
        }
        ~TrapsOn ()
        {
            std::feclearexcept(FE_ALL_EXCEPT);
            amrex::setFPExcept(prev);
            amrex::ParmParse pa("amrex");
            for (std::size_t i = 0; i < keys.size(); ++i) {
                if (had[i] != 0) { pa.add(keys[i].c_str(), was[i] != 0); } else { pa.remove(keys[i].c_str()); }
            }
        }
        TrapsOn (const TrapsOn&) = delete;
        TrapsOn& operator= (const TrapsOn&) = delete;
    };
    const TrapsOn traps;
    const bool trapping = (amrex::getFPExcept() == amrex::FPExcept::all);
    {
        const std::string dir = scratch("traps");
        set_inputs(dir, "Ttraps", false);
        Mesh m(false);
        m.fill_crosswind();
        auto c = Conductors::create(0);
        ASSERT_TRUE(c);
        c->set_ground(nullptr, m.geom);   // MoorDyn_Init and its initial-condition solve
        const auto& span = *c->lines().front();
        for (int s = 0; s < 8; ++s) { c->advance(0, 0.25 * s, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom); }
        EXPECT_GT(span.mid_offset(), 0.1) << "the span must blow out in the sampled crosswind";
        if (trapping) { EXPECT_EQ(amrex::getFPExcept(), amrex::FPExcept::all) << "the traps must still be on"; }
    }
}

TEST(Conductors, GustsComeFromTheRANSkAlongEachSpan)
{
    const std::string dir = scratch("gusts");
    set_inputs(dir, "Tgust", false);
    amrex::ParmParse pp("conductors");
    // a second line at 45 degrees to the wind and 10 m higher: the normal component and the second line's height
    pp.remove("lines");
    pp.addarr("lines", std::vector<std::string>{"Tgust", "Tgust2"});
    {
        amrex::ParmParse ps("conductors.Tgust2");
        ps.addarr("end_a", std::vector<Real>{300.0, 300.0, 40.0});
        ps.addarr("end_b", std::vector<Real>{513.0, 513.0, 40.0});
        ps.add("length", 301.5);
        ps.add("diameter", 0.0281);
        ps.add("mass_per_length", 1.628);
        ps.add("axial_stiffness", 3.0e7);
    }
    pp.add("gust_type", std::string("factor"));
    // the steps end at 0.25, 0.5, ...: the first is not sampled
    pp.add("stats_start", 0.5);
    Mesh m(false);
    m.fill_crosswind();
    // rho and rho k with k linear in the height, which the cell sampler reproduces exactly
    const Real rho0 = 1.2, k0 = 0.5, kz = 0.01;
    amrex::MultiFab cons(m.ba, m.dm, RhoKE_comp + 1, 1);
    cons.setVal(0.0);
    const Real dz = m.H / m.nz;
    for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
        auto ca = cons.array(mfi);
        amrex::LoopOnCpu(mfi.growntilebox(), [&](int i, int j, int k) {
            ca(i,j,k,Rho_comp) = rho0;
            ca(i,j,k,RhoKE_comp) = rho0 * (k0 + kz * (k + 0.5) * dz);
        });
    }
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(nullptr, m.geom);
    {
        using erf_gtest::abort_message;
        EXPECT_NE(abort_message([&] { c->advance(0, 0.0, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom, &cons); })
                      .find("needs set_closure() first"), std::string::npos);
        EXPECT_NE(abort_message([&] { c->set_closure(false, 0.5562); }).find("needs the k-equation RANS"), std::string::npos);
    }
    c->set_closure(true, 0.5562);
    // Cmu0 arrives as amrex::Real: a float in single precision
    EXPECT_NEAR(c->gust_sigma_factor(), 2.5 * 0.5562, 1e-6);
    EXPECT_NE(erf_gtest::abort_message([&] { c->advance(0, 0.0, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom); })
                  .find("needs the conserved state"), std::string::npos);
    ASSERT_EQ(c->lines().size(), 2u);
    // the normals of the two chords' horizontal projections: +y for the line along x, (-1, 1)/sqrt(2) for the other
    const double ny[2] = {1.0, 1.0 / std::sqrt(2.0)};
    // the expected means: every sampled step, the root-mean-square wind and normal wind over each span's nodes and
    // the mean k, where the nodes are when the wind is sampled (u = w = 0, so the wind is |v|)
    double wind[2] = {0.0, 0.0}, normal[2] = {0.0, 0.0}, kmean[2] = {0.0, 0.0};
    const int steps = 8;
    for (int s = 0; s < steps; ++s) {
        if (s == 1) {
            EXPECT_FALSE(std::filesystem::exists(dir + "/gusts.csv")) << "no table before the first sample (stats_start)";
            EXPECT_FALSE(std::filesystem::exists(dir + "/Tgust_gusts_stats.csv")) << "no statistics before stats_start";
        }
        if (s > 0) {
            for (int l = 0; l < 2; ++l) {
                const auto& span = *c->lines()[static_cast<std::size_t>(l)];
                const std::vector<Real> where = span.kinematics_points();
                double w2 = 0.0, kk = 0.0;
                for (unsigned n = 0; n < span.num_nodes(); ++n) {
                    const double v = static_cast<double>(m.v_at(where[3*n+1], where[3*n+2]));
                    w2 += v * v;
                    kk += static_cast<double>(k0 + kz * where[3*n+2]);
                }
                wind[l] += std::sqrt(w2 / span.num_nodes());
                normal[l] += ny[l] * std::sqrt(w2 / span.num_nodes());
                kmean[l] += kk / span.num_nodes();
            }
        }
        c->advance(0, 0.25 * s, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom, &cons);
    }
    {
        // the running statistics count the sampled steps only
        std::ifstream st(dir + "/Tgust_gusts_stats.csv");
        std::string h, r;
        ASSERT_TRUE(std::getline(st, h) && std::getline(st, r));
        EXPECT_EQ(r.substr(0, r.find(',')), std::to_string(steps - 1));
    }
    std::ifstream f(dir + "/gusts.csv");
    ASSERT_TRUE(f.good());
    std::string header, row;
    std::getline(f, header);
    EXPECT_EQ(header, "line,span,height,chord,wind,normal_wind,k,sigma,normal_sigma,intensity,gust_response,gust_wind,mean_load,"
                      "peak_load,valid");
    const double tol = (std::is_same<Real, float>::value) ? 1e-4 : 1e-8;
    const char* names[2] = {"Tgust", "Tgust2"};
    const double height[2] = {30.0, 40.0}, chord[2] = {300.0, 213.0 * std::sqrt(2.0)};
    for (int l = 0; l < 2; ++l) {
        ASSERT_TRUE(std::getline(f, row)) << "one row per span";
        std::vector<std::string> col;
        std::stringstream ss(row);
        for (std::string x; std::getline(ss, x, ',');) { col.push_back(x); }
        ASSERT_EQ(col.size(), 15u);
        EXPECT_EQ(col[0], names[l]);
        EXPECT_EQ(col[1], "1");
        const double U = wind[l] / (steps - 1), Un = normal[l] / (steps - 1), k = kmean[l] / (steps - 1);
        EXPECT_NEAR(std::stod(col[2]), height[l], tol * height[l]) << "the attachment height above the flat ground";
        EXPECT_NEAR(std::stod(col[3]), chord[l], tol * chord[l]);
        EXPECT_NEAR(std::stod(col[4]), U, tol * U);
        EXPECT_NEAR(std::stod(col[5]), Un, tol * Un);
        EXPECT_NEAR(std::stod(col[6]), k, tol * k);
        // the formulas are pinned by the Gusts unit tests; here the columns must carry them
        const auto G = erf_conductors::span_gust(U, Un, k, chord[l], 0.0281, 1.0, 1.2, 2.5 * 0.5562, 2.7, 67.056);
        EXPECT_NEAR(std::stod(col[7]), G.sigma, tol * G.sigma);
        EXPECT_NEAR(std::stod(col[8]), G.normal_sigma, tol * G.normal_sigma);
        EXPECT_NEAR(std::stod(col[9]), G.intensity, tol * G.intensity);
        EXPECT_NEAR(std::stod(col[10]), G.gust_response, tol * G.gust_response);
        EXPECT_NEAR(std::stod(col[11]), G.gust_wind, tol * G.gust_wind);
        EXPECT_NEAR(std::stod(col[12]), 0.5 * 1.2 * 1.0 * 0.0281 * Un * Un, tol * G.mean_load);
        EXPECT_NEAR(std::stod(col[13]), G.peak_load, tol * G.peak_load);
        EXPECT_EQ(col[14], G.linear_valid ? "1" : "0");
        EXPECT_TRUE(G.linear_valid) << "about 10 and 7 m/s normal to the spans against sigma_n of about 1.2 m/s";
    }
    EXPECT_FALSE(std::getline(f, row)) << "one row per span";
}

TEST(Conductors, GustsSwitchedOnAtARestartStartAfresh)
{
    const std::string dir = scratch("gust_restart");
    set_inputs(dir, "Tfresh", false);
    Mesh m(false);
    m.fill_crosswind();
    amrex::MultiFab cons(m.ba, m.dm, RhoKE_comp + 1, 1);
    cons.setVal(0.0);
    for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
        auto ca = cons.array(mfi);
        amrex::LoopOnCpu(mfi.growntilebox(), [&](int i, int j, int k) {
            ca(i,j,k,Rho_comp) = 1.2;
            ca(i,j,k,RhoKE_comp) = 1.2 * 0.8;
        });
    }
    const double dt = 0.25;
    // a run without gusts writes a checkpoint after three steps
    auto a = Conductors::create(0);
    ASSERT_TRUE(a);
    a->set_closure(true, 0.5562);
    a->set_ground(nullptr, m.geom);
    for (int s = 0; s < 3; ++s) { a->advance(0, dt * s, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom, &cons); }
    const std::string chk = dir + "/chk00003";
    std::filesystem::create_directories(chk);
    a->write_checkpoint(chk);
    EXPECT_FALSE(std::filesystem::exists(dir + "/gusts.csv"));
    // the restart turns gusts on: the lines continue, the gust statistics start at the restart
    amrex::ParmParse("conductors").add("gust_type", std::string("factor"));
    auto b = Conductors::create(0);
    ASSERT_TRUE(b);
    b->set_closure(true, 0.5562);
    b->set_ground(nullptr, m.geom, chk);
    EXPECT_TRUE(b->restored());
    for (int s = 3; s < 6; ++s) { b->advance(0, dt * s, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom, &cons); }
    std::ifstream st(dir + "/Tfresh_gusts_stats.csv");
    std::string h, r;
    ASSERT_TRUE(std::getline(st, h) && std::getline(st, r));
    std::vector<std::string> col;
    std::stringstream ss(r);
    for (std::string x; std::getline(ss, x, ',');) { col.push_back(x); }
    ASSERT_GE(col.size(), 3u);
    EXPECT_EQ(col[0], "3") << "the three steps after the restart";
    EXPECT_DOUBLE_EQ(std::stod(col[1]), 0.75) << "the first sample is the restart step's start";
    EXPECT_TRUE(std::filesystem::exists(dir + "/gusts.csv"));
}

TEST(Conductors, ClearanceIsTheHeightAboveTheTerrainUnderEachNode)
{
    const std::string dir = scratch("clearance");
    set_inputs(dir, "Tclear");
    amrex::ParmParse("conductors").add("node_output_int", 1);
    Mesh m(true);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(m.znd.get(), m.geom);
    for (int s = 0; s < 3; ++s) { c->advance(0, 0.2 * s, 0.2, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); }
    const auto& span = *c->lines().front();
    Real lowest = 1.0e30;
    unsigned lowest_node = 0;
    for (unsigned n = 0; n < span.num_nodes(); ++n) {
        const auto p = span.node_position(n);
        const Real expect = p[2] - m.h(p[0], p[1]);   // the ramp: the bilinear surface is exact
        EXPECT_NEAR(span.clearance(n), expect, roundoff * std::max(Real(1.0), std::abs(expect))) << "node " << n;
        if (expect < lowest) { lowest = expect; lowest_node = n; }
    }
    unsigned node = 0;
    EXPECT_NEAR(span.min_clearance(node), lowest, roundoff * std::max(Real(1.0), std::abs(lowest)));
    EXPECT_EQ(node, lowest_node);
    EXPECT_LT(lowest, 30.0) << "the sagging middle sits lower above the ground than the attachments";
    // the node file: one row per node per write (the initial write and three steps)
    std::ifstream f(dir + "/Tclear_nodes.dat");
    ASSERT_TRUE(f.good());
    std::string line;
    int rows = 0;
    while (std::getline(f, line)) { if (!line.empty() && line.rfind("time", 0) != 0) { ++rows; } }
    EXPECT_EQ(rows, 4 * static_cast<int>(span.num_nodes()));
    // the statistics file lists every quantity
    std::ifstream st(dir + "/Tclear_stats.csv");
    ASSERT_TRUE(st.good());
    std::string all((std::istreambuf_iterator<char>(st)), std::istreambuf_iterator<char>());
    for (const char* q : {"swing_deg", "mid_offset", "tension_a", "tension_b", "max_tension", "min_clearance", "drag_normal"}) {
        EXPECT_NE(all.find(q), std::string::npos) << q;
    }
}

namespace {
// the rows of a seafloor file: its ticks and its heights (MoorDyn's frame), x the slower index
struct SeafloorFile {
    std::vector<double> x, y, z;
    explicit SeafloorFile (const std::string& fname)
    {
        std::ifstream in(fname);
        std::size_t nx = 0, ny = 0;
        in >> nx >> ny;
        x.resize(nx); y.resize(ny);
        for (auto& v : x) { in >> v; }
        for (auto& v : y) { in >> v; }
        for (std::size_t k = 0; k < nx * ny; ++k) {
            double px = 0.0, py = 0.0, pz = 0.0;
            in >> px >> py >> pz;
            EXPECT_DOUBLE_EQ(px, x[k / ny]) << "row " << k;
            EXPECT_DOUBLE_EQ(py, y[k % ny]) << "row " << k;
            z.push_back(pz);
        }
        EXPECT_TRUE(static_cast<bool>(in)) << fname;
    }
};
} // namespace

TEST(Conductors, TheSeafloorHoldsTheTerrainNodesUnderTheLineAndOneSpanAround)
{
    // the 300 m span along x at y = 500 on the ramp h = 0.05 x + 0.02 y (100 m cells): one span length,
    // 301.5 m, around the line covers x from -1.5 (the domain's edge, 0) to 901.5 (node 1000) and y
    // from 198.5 (node 100) to 801.5 (node 900)
    const std::string dir = scratch("seafloor");
    set_inputs(dir, "Tsea");
    amrex::ParmParse("conductors").add("ground_contact", true);
    Mesh m(true);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(m.znd.get(), m.geom);
    amrex::ParmParse("conductors").remove("ground_contact");
    const SeafloorFile sf(dir + "/Tsea.seafloor.txt");
    ASSERT_EQ(sf.x.size(), 11u);
    ASSERT_EQ(sf.y.size(), 9u);
    for (std::size_t i = 0; i < sf.x.size(); ++i) { EXPECT_DOUBLE_EQ(sf.x[i], 100.0 * i); }
    for (std::size_t j = 0; j < sf.y.size(); ++j) { EXPECT_DOUBLE_EQ(sf.y[j], 100.0 * (j + 1)); }
    // the terrain at each node raised by the conductor's radius, lowered by surface_offset into MoorDyn's frame
    for (std::size_t k = 0; k < sf.z.size(); ++k) {
        const double expect = static_cast<double>(m.h(static_cast<Real>(sf.x[k / sf.y.size()]), static_cast<Real>(sf.y[k % sf.y.size()]))) +
                              0.5 * 0.0281 - 10000.0;
        EXPECT_NEAR(sf.z[k], expect, 1.0e-6) << "point " << k;
    }
    // the MoorDyn input names it, with the contact stiffness and the conductor's critical damping
    std::ifstream md(dir + "/Tsea.moordyn.txt");
    const std::string text((std::istreambuf_iterator<char>(md)), std::istreambuf_iterator<char>());
    EXPECT_NE(text.find(dir + "/Tsea.seafloor.txt   SeafloorFile"), std::string::npos) << text;
    // the numbers before kbot and cbot (a single-precision build writes 3e+06): 3e6 and 2 sqrt(3e6 x 1.628 / 0.0281) = 26367 Pa s/m
    auto number_before = [&] (const std::string& key) {
        const std::size_t at = text.find(key);
        EXPECT_NE(at, std::string::npos) << key;
        const std::size_t from = text.rfind('\n', at) + 1;
        return std::stod(text.substr(from, at - from));
    };
    EXPECT_NEAR(number_before("   kbot"), 3.0e6, 1.0);
    EXPECT_NEAR(number_before("   cbot"), 2.0 * std::sqrt(3.0e6 * 1.628 / 0.0281), 0.1);
}

TEST(Conductors, ASlackSpanRestsOnTheTerrainWithGroundContactAndStopsTheRunWithout)
{
    // 330 m of conductor over a 300 m chord sags about 58 m from ends 30 m above the ramp. Without contact its
    // middle hangs into the terrain, where no grid holds the wind, and the run stops; with contact the middle
    // two thirds of it must come to rest on the ramp, not below it
    for (const bool contact : {true, false}) {
        const std::string tag = contact ? "on" : "off";
        const std::string dir = scratch("ground_" + tag);
        const std::string name = "Tgr" + tag;
        set_inputs(dir, name, false);
        amrex::ParmParse ps("conductors." + name);
        ps.remove("length");
        ps.add("length", 330.0);
        if (contact) { amrex::ParmParse("conductors").add("ground_contact", true); }
        Mesh m(true);
        auto c = Conductors::create(0);
        ASSERT_TRUE(c);
        c->set_ground(m.znd.get(), m.geom);
        // MoorDyn-C settles the line dynamically (its initial-shape solve stops early with the span on the ground)
        auto run = [&] { for (int s = 0; s < 30; ++s) { c->advance(0, 0.5 * s, 0.5, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); } };
        if (!contact) {
            const std::string msg = erf_gtest::abort_message(run);
            EXPECT_NE(msg.find("was sampled by 0 boxes"), std::string::npos) << "a node below the terrain: " << msg;
            EXPECT_FALSE(std::filesystem::exists(dir + "/ground_contact.csv"));
            continue;
        }
        run();
        amrex::ParmParse("conductors").remove("ground_contact");
        const auto& span = *c->lines().front();
        unsigned low = 0;
        const Real lowest = span.min_clearance(low);
        // the conductor's axis rests one radius above the terrain: the stub places it there, MoorDyn-C presses it
        // mu g / (kbot d) = 0.19 mm deeper and reads the seafloor in float (0.5 mm at 1e4 m)
        const Real tol = erf_moordyn::is_stub() ? postol : Real(1.0e-3);
        EXPECT_NEAR(lowest, Real(0.5 * 0.0281), tol) << "node " << low << " rests on the terrain";
        const int on = span.nodes_on_ground(0);
        EXPECT_GE(on, 5) << "a long stretch of the span lies on the ground";
        EXPECT_LT(on, static_cast<int>(span.span_num_nodes(0))) << "the ends stay up on their attachments";
        int counted = 0;
        for (unsigned n = 0; n < span.span_num_nodes(0); ++n) {
            if (span.clearance(n) <= Real(0.0281)) {
                ++counted;
                // every node on the ground rests on the ramp where it lies, not only the lowest
                EXPECT_NEAR(span.clearance(n), Real(0.5 * 0.0281), tol) << "node " << n;
            }
        }
        EXPECT_EQ(on, counted);
        // the log: the span touched the ground at the end of the first step and has as many nodes on it now
        std::ifstream f(dir + "/ground_contact.csv");
        std::string header, row;
        ASSERT_TRUE(std::getline(f, header) && std::getline(f, row));
        EXPECT_EQ(header, "span,nodes,nodes_on_ground,most_on_ground,most_time,first_touch_time");
        std::istringstream r(row);
        std::vector<std::string> fields;
        for (std::string v; std::getline(r, v, ',');) { fields.push_back(v); }
        ASSERT_EQ(fields.size(), 6u) << row;
        EXPECT_EQ(fields[0], name);
        EXPECT_EQ(std::stoi(fields[1]), static_cast<int>(span.span_num_nodes(0)));
        EXPECT_EQ(std::stoi(fields[2]), on);
        EXPECT_GE(std::stoi(fields[3]), on);
        EXPECT_DOUBLE_EQ(std::stod(fields[5]), 0.5);
    }
}

TEST(Conductors, GroundContactRefusesAMoorDynStepTooLongForIt)
{
    // five 60 m segments and moordyn_dt = 0.01 s: MoorDyn-C's Courant step (8.8 ms at CFL 0.1) or the stub's 10 ms,
    // against twice the contact frequency sqrt(3e6 x 0.0281 / 1.628) = 227.6 rad/s: 4.0 and 4.6, above 2
    const std::string dir = scratch("ground_step");
    set_inputs(dir, "Tgs");
    amrex::ParmParse ps("conductors.Tgs");
    ps.add("segments", 5);
    amrex::ParmParse pp("conductors");
    pp.add("ground_contact", true);
    pp.add("moordyn_dt", 0.01);
    Mesh m(true);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    std::string msg = erf_gtest::abort_message([&] { c->set_ground(m.znd.get(), m.geom); });
    pp.remove("ground_contact");
    pp.remove("moordyn_dt");
    EXPECT_NE(msg.find("conductors.Tgs: with ground_contact MoorDyn's internal step"), std::string::npos) << msg;
    EXPECT_NE(msg.find("lower conductors.moordyn_dt or conductors.moordyn_cfl"), std::string::npos) << msg;
    // a conductor too heavy for its diameter would sink into the ground: 200 kg/m at 28.1 mm sinks
    // mu g / (kbot d) = 23 mm, more than a quarter of its diameter (7 mm)
    const std::string dir2 = scratch("ground_sink");
    set_inputs(dir2, "Tgk");
    amrex::ParmParse pk("conductors.Tgk");
    pk.remove("mass_per_length");
    pk.add("mass_per_length", 200.0);
    pp.add("ground_contact", true);
    auto k = Conductors::create(0);
    ASSERT_TRUE(k);
    msg = erf_gtest::abort_message([&] { k->set_ground(m.znd.get(), m.geom); });
    pp.remove("ground_contact");
    EXPECT_NE(msg.find("conductors.Tgk: with ground_contact the conductor would sink"), std::string::npos) << msg;
}

TEST(Conductors, GroundContactRecordsContinueAcrossARestart)
{
    // the slack span of the test above on the ramp, checkpointed after four steps of 0.5 s; the first run goes on
    // four more, and the restarted run must carry the same records to the same end
    const std::string dir = scratch("ground_restart");
    set_inputs(dir, "Tgrr", false);
    amrex::ParmParse ps("conductors.Tgrr");
    ps.remove("length");
    ps.add("length", 330.0);
    amrex::ParmParse pp("conductors");
    pp.add("ground_contact", true);
    Mesh m(true);
    const double dt = 0.5;
    auto a = Conductors::create(0);
    a->set_ground(m.znd.get(), m.geom);
    for (int s = 0; s < 4; ++s) { a->advance(0, dt * s, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); }
    const std::string chk = dir + "/chk00004";
    std::filesystem::create_directories(chk);
    a->write_checkpoint(chk);
    ASSERT_TRUE(std::filesystem::exists(chk + "/conductors/ground_contact"));
    const erf_conductors::GroundSpan at_chk = a->ground_contacts().span(0);
    ASSERT_TRUE(at_chk.touched);
    for (int s = 4; s < 8; ++s) { a->advance(0, dt * s, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); }
    auto b = Conductors::create(0);
    b->set_ground(m.znd.get(), m.geom, chk);
    ASSERT_EQ(b->ground_contacts().size(), 1u);
    const auto& r = b->ground_contacts().span(0);
    EXPECT_TRUE(r.touched);
    EXPECT_EQ(r.first_time, at_chk.first_time) << "restored, not the restart's time";
    EXPECT_EQ(r.now, at_chk.now);
    EXPECT_EQ(r.most, at_chk.most);
    EXPECT_EQ(r.most_time, at_chk.most_time);
    for (int s = 4; s < 8; ++s) { b->advance(0, dt * s, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); }
    const auto& ra = a->ground_contacts().span(0);
    const auto& rb = b->ground_contacts().span(0);
    EXPECT_EQ(rb.first_time, ra.first_time);
    EXPECT_EQ(rb.now, ra.now);
    EXPECT_EQ(rb.most, ra.most);
    EXPECT_EQ(rb.most_time, ra.most_time);
    // a restart without ground_contact keeps no records
    pp.remove("ground_contact");
    auto o = Conductors::create(0);
    ASSERT_TRUE(std::filesystem::exists(dir + "/ground_contact.csv"));
    o->set_ground(m.znd.get(), m.geom, chk);
    EXPECT_TRUE(o->ground_contacts().empty());
    EXPECT_FALSE(std::filesystem::exists(dir + "/ground_contact.csv")) << "the earlier run's table is not this run's";
    // a checkpoint written without ground_contact restarted with it: the records start at the restart
    const std::string chk2 = dir + "/chk_off";
    std::filesystem::create_directories(chk2);
    o->write_checkpoint(chk2);
    EXPECT_FALSE(std::filesystem::exists(chk2 + "/conductors/ground_contact"));
    pp.add("ground_contact", true);
    auto n = Conductors::create(0);
    n->set_ground(m.znd.get(), m.geom, chk2);
    ASSERT_EQ(n->ground_contacts().size(), 1u);
    EXPECT_FALSE(n->ground_contacts().span(0).touched);
    n->advance(0, dt * 4, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom);
    EXPECT_EQ(n->ground_contacts().span(0).first_time, dt * 5) << "first touched at the end of the restart's first step";
    pp.remove("ground_contact");
}

TEST(Conductors, OnAMeshWithoutTerrainTheSpanRestsOnTheDomainsFloor)
{
    // the slack span over a flat mesh without terrain: the ground is z = prob_lo[2] = 0, and the conductor's axis rests
    // one radius above it (MoorDyn-C presses 0.19 mm deeper); without the radius it would leave the domain
    const std::string dir = scratch("ground_flat");
    set_inputs(dir, "Tgf", false);
    amrex::ParmParse ps("conductors.Tgf");
    ps.remove("length");
    ps.add("length", 330.0);
    amrex::ParmParse("conductors").add("ground_contact", true);
    Mesh m(false);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(nullptr, m.geom);
    for (int s = 0; s < 30; ++s) { c->advance(0, 0.5 * s, 0.5, m.u, m.v, m.w, nullptr, nullptr, m.geom); }
    amrex::ParmParse("conductors").remove("ground_contact");
    const SeafloorFile sf(dir + "/Tgf.seafloor.txt");
    for (const double z : sf.z) { EXPECT_NEAR(z, 0.5 * 0.0281 - 10000.0, 1.0e-6); }
    const auto& span = *c->lines().front();
    unsigned low = 0;
    const Real tol = erf_moordyn::is_stub() ? postol : Real(1.0e-3);
    EXPECT_NEAR(span.min_clearance(low), Real(0.5 * 0.0281), tol);
    EXPECT_GE(span.nodes_on_ground(0), 5);
}

TEST(Conductors, TheContactRateIsTheFastestRootOfEachLineTypeOnTheGround)
{
    using erf_conductors::ConductorInputs;
    const Real mu = 1.628, d = 0.0281;
    const Real cb = ConductorInputs::ground_damping(mu, d);
    const Real w = ConductorInputs::ground_frequency(mu, d);
    EXPECT_NEAR(w, std::sqrt(3.0e6 * 0.0281 / 1.628), 1.0e-3);
    // critical: 2 omega
    EXPECT_NEAR(ConductorInputs::ground_rate(mu, d, cb), 2 * w, 1.0e-3 * w);
    // the default string (60 kg over 2.5 m, 0.254 m) under the conductor's cb is under-damped (ratio 0.78): 2 omega_j
    const Real wi = ConductorInputs::ground_frequency(24.0, 0.254);
    EXPECT_NEAR(ConductorInputs::ground_rate(24.0, 0.254, cb), 2 * wi, 1.0e-3 * wi);
    // a light string, 20 kg over 2.5 m, is over-damped: 2 (a + sqrt(a^2 - omega^2)), a = cb d / (2 mu)
    const Real wl = ConductorInputs::ground_frequency(8.0, 0.254);
    const Real a = cb * Real(0.254) / Real(16.0);
    ASSERT_GT(a, wl);
    EXPECT_NEAR(ConductorInputs::ground_rate(8.0, 0.254, cb), 2 * (a + std::sqrt(a * a - wl * wl)), 1.0e-3 * a);
    EXPECT_GT(ConductorInputs::ground_rate(8.0, 0.254, cb), cb * Real(0.254) / Real(8.0)) << "more than cb d / mu";
}

TEST(Conductors, TheSpreadDragIntegratesToMinusTheDragOnTheLines)
{
    const std::string dir = scratch("drag");
    set_inputs(dir, "Tdrag");   // a prescribed 10 m/s crosswind along +y
    amrex::ParmParse("conductors").add("drag_on_flow", true);
    Mesh m(true);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    EXPECT_TRUE(c->drag_on_flow());
    c->set_ground(m.znd.get(), m.geom);
    // a terrain-following mesh: the face volumes carry detJ, the cell's physical height over its nominal
    // height, which shrinks where the ramp rises (the columns span H - h)
    amrex::MultiFab detJ(m.ba, m.dm, 1, 1);
    for (amrex::MFIter mfi(detJ); mfi.isValid(); ++mfi) {
        const auto za = m.znd->const_array(mfi);
        auto ja = detJ.array(mfi);
        const Real dz = m.H / m.nz;
        amrex::LoopOnCpu(mfi.growntilebox(), [&](int i, int j, int k) {
            ja(i,j,k) = Real(0.25) * ((za(i,j,k+1) - za(i,j,k)) + (za(i+1,j,k+1) - za(i+1,j,k)) +
                                      (za(i,j+1,k+1) - za(i,j+1,k)) + (za(i+1,j+1,k+1) - za(i+1,j+1,k))) / dz;
        });
    }
    EXPECT_LT(detJ.min(0), 0.99 * detJ.max(0)) << "the Jacobian must vary over the ramp";
    for (int s = 0; s < 2; ++s) { c->advance(0, 0.2 * s, 0.2, m.u, m.v, m.w, m.znd.get(), &detJ, m.geom); }
    const auto& span = *c->lines().front();
    const auto drag = span.total_drag();
    EXPECT_GT(drag[1], 1.0) << "a crosswind along +y must drag the line along +y";
    const auto& integral = c->source_integral();
    for (int d = 0; d < 3; ++d) {
        EXPECT_NEAR(integral[d], -drag[d], Real(0.1) * roundoff * std::max(Real(1.0), std::abs(drag[d]))) << "component " << d;
    }
    // the sources ERF adds integrate to the same force
    amrex::MultiFab sx(amrex::convert(m.ba, amrex::IntVect(1,0,0)), m.dm, 1, 0);
    amrex::MultiFab sy(amrex::convert(m.ba, amrex::IntVect(0,1,0)), m.dm, 1, 0);
    amrex::MultiFab sz(amrex::convert(m.ba, amrex::IntVect(0,0,1)), m.dm, 1, 0);
    sx.setVal(0.0); sy.setVal(0.0); sz.setVal(0.0);
    c->add_momentum_sources(0, sx, sy, sz);
    EXPECT_NEAR(erf_actuator::integrate_source(1, sy, &detJ, m.geom), -drag[1], Real(0.1) * roundoff * std::abs(drag[1]));
    c->add_momentum_sources(1, sx, sy, sz);   // not the anchor level: nothing is added
    EXPECT_NEAR(erf_actuator::integrate_source(1, sy, &detJ, m.geom), -drag[1], Real(0.1) * roundoff * std::abs(drag[1]));
    // the cell-centred plot field
    amrex::MultiFab cells(m.ba, m.dm, 3, 0);
    c->cell_sources(0, cells, 0);
    EXPECT_LT(cells.min(1), 0.0) << "the air is pushed against the wind (-y) around the line";
    c->cell_sources(1, cells, 0);
    EXPECT_DOUBLE_EQ(cells.norm0(1), 0.0);
}

TEST(Conductors, WithoutDragOnFlowNothingIsPutIntoTheFlow)
{
    const std::string dir = scratch("nodrag");
    set_inputs(dir, "Tnodrag");
    Mesh m(false);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    EXPECT_FALSE(c->drag_on_flow());
    c->set_ground(nullptr, m.geom);
    c->advance(0, 0.0, 0.2, m.u, m.v, m.w, nullptr, nullptr, m.geom);
    EXPECT_GT(c->lines().front()->total_drag()[1], 1.0) << "the line still feels the drag";
    for (int d = 0; d < 3; ++d) { EXPECT_DOUBLE_EQ(c->source_integral()[d], 0.0); }
    amrex::MultiFab sy(amrex::convert(m.ba, amrex::IntVect(0,1,0)), m.dm, 1, 0);
    amrex::MultiFab sx(amrex::convert(m.ba, amrex::IntVect(1,0,0)), m.dm, 1, 0);
    amrex::MultiFab sz(amrex::convert(m.ba, amrex::IntVect(0,0,1)), m.dm, 1, 0);
    sx.setVal(0.0); sy.setVal(0.0); sz.setVal(0.0);
    c->add_momentum_sources(0, sx, sy, sz);
    EXPECT_DOUBLE_EQ(sy.norm0(), 0.0);
    // nothing reaches the flow, so there is no support: one-way lines may sit under a finer level
    EXPECT_TRUE(c->support().empty());
}

namespace {
std::string slurp (const std::string& fname)
{
    std::ifstream in(fname);
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}
} // namespace

// asce74.csv's wire: at the effective height, the mean attachment height less a third of the still-air sag
// (vertical), and swinging about an inclined chord against W cos beta; the attachment height and the whole
// weight with asce74_wire_height = attachment and asce74_inclined_spans = false
TEST(Conductors, TheASCE74CheckTakesTheWireAtItsEffectiveHeightOnItsInclinedChord)
{
    const std::string dir = scratch("asce74");
    // ends 30 m and 130 m above flat ground, 300 m apart: a chord inclined by asin(100 / 316.2), 18.4 degrees, steep
    // enough that the sag's vertical part (1 / cos beta = 1.054) and W cos beta stand well clear of the tolerances
    set_inputs(dir, "Tas", true, {{300.0, 500.0, 30.0}}, {{600.0, 500.0, 130.0}});
    amrex::ParmParse("conductors.Tas").add("length", 317.8);
    amrex::ParmParse pp("conductors");
    pp.add("asce74_wind", 40.0);
    Mesh m(false);
    auto row = [&] () {
        auto c = Conductors::create(0);
        c->set_ground(nullptr, m.geom);
        std::istringstream in(slurp(dir + "/asce74.csv"));
        std::string hdr, line;
        std::getline(in, hdr);
        std::getline(in, line);
        std::map<std::string, double> r;
        std::istringstream h(hdr), l(line);
        std::string name, value;
        while (std::getline(h, name, ',') && std::getline(l, value, ',')) { r[name] = std::atof(value.c_str()); }
        return r;
    };
    auto r = row();
    const double c = std::hypot(300.0, 100.0), beta = std::asin(100.0 / c), L = 317.8, W = r["weight"], EA = 3.0e7;
    // an independent still-air sag: the parabola of a span stretched by its mean tension H/EA (to well under 1 %)
    double sag = 13.0;
    for (int it = 0; it < 50; ++it) {
        const double H = W * std::cos(beta) * c * c / (8.0 * sag);
        sag = std::sqrt(3.0 * c * (L * (1.0 + H / EA) - c) / 8.0);
    }
    const double sag_v = sag / std::cos(beta);
    EXPECT_NEAR(r["height"], 80.0 - sag_v / 3.0, 0.01 * sag_v / 3.0) << "mean attachment height 80 m less a third of the sag";
    // the chord in Real: beta to a few units of its precision
    const double swing_tol = std::is_same<Real, float>::value ? 1.0e-5 : 1.0e-8;
    EXPECT_NEAR(std::tan(r["swing_deg"] * 3.14159265358979323846 / 180.0), r["load"] / (W * std::cos(beta)), swing_tol);
    // the whole weight against the swing, the effective height still from the inclined chord's vertical sag
    const double h_effective = r["height"];
    pp.add("asce74_inclined_spans", false);
    r = row();
    EXPECT_EQ(r["height"], h_effective) << "the still-air sag does not depend on the swing's weight";
    EXPECT_NEAR(std::tan(r["swing_deg"] * 3.14159265358979323846 / 180.0), r["load"] / W, swing_tol);
    // the old forms
    pp.add("asce74_wire_height", std::string("attachment"));
    r = row();
    EXPECT_NEAR(r["height"], 80.0, 1e-9);
    EXPECT_NEAR(std::tan(r["swing_deg"] * 3.14159265358979323846 / 180.0), r["load"] / W, swing_tol);
    for (const char* k : {"asce74_wind", "asce74_wire_height", "asce74_inclined_spans"}) { pp.remove(k); }
}

TEST(Conductors, ARestartCannotMoveALinesPointsOrChangeItsLength)
{
    const std::string dir = scratch("restart_moved");
    set_inputs(dir, "Tmv");
    Mesh m(true);
    auto a = Conductors::create(0);
    ASSERT_TRUE(a);
    a->set_ground(m.znd.get(), m.geom);
    for (int step = 0; step < 2; ++step) { a->advance(0, 0.2 * step, 0.2, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); }
    const std::string chk = dir + "/chk00002";
    std::filesystem::create_directories(chk);
    a->write_checkpoint(chk);
    EXPECT_NE(slurp(chk + "/conductors/state").find("geometry Tmv 1 "), std::string::npos);
    amrex::ParmParse ps("conductors.Tmv");
    {
        // the same inputs restart
        auto b = Conductors::create(0);
        EXPECT_TRUE(erf_gtest::abort_message([&] { b->set_ground(m.znd.get(), m.geom, chk); }).empty());
    }
    {
        // end_b moved 1 m: the end segment would be stretched by 1 m at the first step
        ps.addarr("end_b", std::vector<Real>{599.0, 500.0, 30.0});
        auto b = Conductors::create(0);
        const std::string msg = erf_gtest::abort_message([&] { b->set_ground(m.znd.get(), m.geom, chk); });
        EXPECT_NE(msg.find("line Tmv's attachment point coordinate"), std::string::npos) << msg;
        EXPECT_NE(msg.find("cannot move its points"), std::string::npos) << msg;
        ps.addarr("end_b", std::vector<Real>{600.0, 500.0, 30.0});
    }
    {
        ps.add("length", 302.0);
        auto b = Conductors::create(0);
        const std::string msg = erf_gtest::abort_message([&] { b->set_ground(m.znd.get(), m.geom, chk); });
        EXPECT_NE(msg.find("line Tmv's unstretched length"), std::string::npos) << msg;
    }
}

TEST(Conductors, ARestartContinuesTheLinesTheirSourcesStatisticsAndLogs)
{
    const std::string dir = scratch("restart");
    set_inputs(dir, "Trst");   // a prescribed 10 m/s crosswind along +y
    amrex::ParmParse pp("conductors");
    pp.add("drag_on_flow", true);
    pp.add("node_output_int", 2);
    Mesh m(true);
    // the sources need only agree between the two runs here, so a unit Jacobian will do
    amrex::MultiFab detJ(m.ba, m.dm, 1, 1);
    detJ.setVal(1.0);
    const double dt = 0.2;

    // the run that writes the checkpoint after four steps and goes on for three more
    auto a = Conductors::create(0);
    ASSERT_TRUE(a);
    a->set_ground(m.znd.get(), m.geom);
    EXPECT_FALSE(a->restored());
    int step = 0;
    for (; step < 4; ++step) { a->advance(0, dt * step, dt, m.u, m.v, m.w, m.znd.get(), &detJ, m.geom); }
    const std::string chk = dir + "/chk00004";
    std::filesystem::create_directories(chk);
    a->write_checkpoint(chk);
    // the state file records the frame MoorDyn's saved state is in
    EXPECT_NE(slurp(chk + "/conductors/state").find("surface_offset = 10000"), std::string::npos);
    amrex::MultiFab src_at_chk(m.ba, m.dm, 3, 0);
    a->cell_sources(0, src_at_chk, 0);
    ASSERT_GT(src_at_chk.norm0(1), 0.0) << "the lines must push on the air";
    const std::vector<Real> pos_at_chk = a->lines().front()->node_positions();
    for (; step < 7; ++step) { a->advance(0, dt * step, dt, m.u, m.v, m.w, m.znd.get(), &detJ, m.geom); }
    const std::vector<Real> pos_end = a->lines().front()->node_positions();
    const std::string log_a = slurp(dir + "/Trst.dat");
    const std::string nodes_a = slurp(dir + "/Trst_nodes.dat");
    const std::string stats_a = slurp(dir + "/Trst_stats.csv");
    const std::string load_a = slurp(dir + "/total_load.dat");
    ASSERT_FALSE(log_a.empty());

    // the restarted run, in the same directory: it continues from the checkpoint
    auto b = Conductors::create(0);
    ASSERT_TRUE(b);
    b->set_ground(m.znd.get(), m.geom, chk);
    ASSERT_TRUE(b->restored());
    const auto& span_b = *b->lines().front();
    const std::vector<Real> pos_b = span_b.node_positions();
    ASSERT_EQ(pos_b.size(), pos_at_chk.size());
    for (std::size_t i = 0; i < pos_b.size(); ++i) {
        ASSERT_NEAR(pos_b[i], pos_at_chk[i], roundoff * std::max(Real(1.0), std::abs(pos_at_chk[i]))) << "component " << i;
    }
    EXPECT_GT(span_b.mid_offset(), 0.1) << "restored blown out, not hanging still";
    // the logs end at the checkpoint until the restarted run writes again
    EXPECT_LT(slurp(dir + "/Trst.dat").size(), log_a.size());
    // the momentum source is that of the step the checkpoint was written at
    b->restore_sources(0, m.u, m.znd.get(), &detJ, m.geom);
    amrex::MultiFab src_b(m.ba, m.dm, 3, 0);
    b->cell_sources(0, src_b, 0);
    amrex::MultiFab::Subtract(src_b, src_at_chk, 0, 0, 3, 0);
    EXPECT_LE(src_b.norm0(1), tight * src_at_chk.norm0(1));
    // the same three steps: the lines, the logs and the statistics come out as in the run without the restart
    for (step = 4; step < 7; ++step) { b->advance(0, dt * step, dt, m.u, m.v, m.w, m.znd.get(), &detJ, m.geom); }
    const std::vector<Real> pos_b_end = span_b.node_positions();
    for (std::size_t i = 0; i < pos_end.size(); ++i) {
        EXPECT_NEAR(pos_b_end[i], pos_end[i], roundoff * std::max(Real(1.0), std::abs(pos_end[i]))) << "component " << i;
    }
    EXPECT_EQ(slurp(dir + "/Trst.dat"), log_a);
    EXPECT_EQ(slurp(dir + "/Trst_nodes.dat"), nodes_a);
    EXPECT_EQ(slurp(dir + "/Trst_stats.csv"), stats_a);
    EXPECT_EQ(slurp(dir + "/total_load.dat"), load_a);

    // a checkpoint without conductor state (a precursor's, say): the lines start afresh, from
    // still air, at the restart time; MoorDyn's clock starts at zero there
    const std::string bare = dir + "/chk_bare";
    std::filesystem::create_directories(bare);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(m.znd.get(), m.geom, bare);
    EXPECT_FALSE(c->restored());
    EXPECT_NEAR(c->lines().front()->mid_offset(), 0.0, 1.0e-3);
    for (step = 7; step < 9; ++step) { c->advance(0, dt * step, dt, m.u, m.v, m.w, m.znd.get(), &detJ, m.geom); }
    EXPECT_GT(c->lines().front()->mid_offset(), 0.01) << "the span started at the restart time moves";
    // and a restart from a checkpoint of that run continues it on the same clock
    const std::string chk_c = dir + "/chk00009";
    std::filesystem::create_directories(chk_c);
    c->write_checkpoint(chk_c);
    auto d = Conductors::create(0);
    ASSERT_TRUE(d);
    d->set_ground(m.znd.get(), m.geom, chk_c);
    ASSERT_TRUE(d->restored());
    c->advance(0, dt * step, dt, m.u, m.v, m.w, m.znd.get(), &detJ, m.geom);
    d->advance(0, dt * step, dt, m.u, m.v, m.w, m.znd.get(), &detJ, m.geom);
    const std::vector<Real> pc = c->lines().front()->node_positions();
    const std::vector<Real> pd = d->lines().front()->node_positions();
    for (std::size_t i = 0; i < pc.size(); ++i) {
        EXPECT_NEAR(pd[i], pc[i], roundoff * std::max(Real(1.0), std::abs(pc[i]))) << "component " << i;
    }
}

TEST(Conductors, TrimmingALogKeepsTheHeaderAndTheRowsUpToTheCheckpoint)
{
    const std::string dir = scratch("trim");
    const std::string f = dir + "/log.dat";
    {
        std::ofstream out(f);
        out << "time a b\n# a comment\n0.5 1 2\n1 3 4\n1.0000000001 9 9\n1.5 5 6\n2 7 8\n";
    }
    erf_conductors::trim_log_after(f, 1.0);
    EXPECT_EQ(slurp(f), "time a b\n# a comment\n0.5 1 2\n1 3 4\n1.0000000001 9 9\n")
        << "a row printed to ten digits at the checkpoint time stays; later rows go";
    erf_conductors::trim_log_after(f, 0.0);
    EXPECT_EQ(slurp(f), "time a b\n# a comment\n");
    erf_conductors::trim_log_after(dir + "/missing.dat", 1.0);
    EXPECT_FALSE(std::filesystem::exists(dir + "/missing.dat"));
    // a log whose rows carry the time a step starts at loses the row at the checkpoint time too
    {
        std::ofstream out(f, std::ios::trunc);
        out << "time a\n0.5 1\n1 3\n1.0000000001 9\n1.5 5\n";
    }
    erf_conductors::trim_log_after(f, 1.0, true);
    EXPECT_EQ(slurp(f), "time a\n0.5 1\n");
}

TEST(Conductors, ASectionIsPlacedOnTheTerrainAtEveryTowerAndLogsEachSpan)
{
    const std::string dir = scratch("section");
    set_inputs(dir, "Tsec", true, {{100.0, 500.0, 30.0}}, {{1000.0, 500.0, 30.0}});
    amrex::ParmParse ps("conductors.Tsec");
    ps.remove("length");
    ps.addarr("length", std::vector<Real>{301.5, 301.5, 301.5});
    ps.addarr("towers", std::vector<Real>{400.0, 500.0, 30.0, 700.0, 500.0, 30.0});
    ps.add("insulator_length", 2.5);
    ps.add("insulator_mass", 60.0);
    Mesh m(true);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(m.znd.get(), m.geom);
    const auto& line = *c->lines().front();
    ASSERT_EQ(line.num_spans(), 3);
    EXPECT_EQ(line.inputs().lengths, (std::vector<Real>{301.5, 301.5, 301.5})) << "given lengths are kept";
    // the towers stand on the ramp, and the conductor hangs 2.5 m under each
    for (int k = 1; k <= 2; ++k) {
        const Real x = 100.0 + 300.0 * k;
        const auto p = line.node_position(line.span_first_node(k));
        EXPECT_NEAR(p[0], x, 0.05) << "tower " << k;
        EXPECT_NEAR(p[2], m.h(x, 500.0) + 30.0 - 2.5, 0.05) << "tower " << k;
    }
    std::ifstream g(dir + "/ground.dat");
    std::string header, name, label;
    Real x, y, ground, z;
    std::getline(g, header);
    std::vector<std::string> labels;
    while (g >> name >> label >> x >> y >> ground >> z) {
        labels.push_back(label);
        EXPECT_NEAR(ground, m.h(x, y), postol) << label;
        EXPECT_NEAR(z, m.h(x, y) + 30.0, postol) << label;
    }
    EXPECT_EQ(labels, (std::vector<std::string>{"a", "t1", "t2", "b"}));
    for (int s = 0; s < 3; ++s) { c->advance(0, 0.2 * s, 0.2, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); }
    // a log per span, one for the strings, and their statistics
    for (const char* f : {"/Tsec_span1.dat", "/Tsec_span2.dat", "/Tsec_span3.dat", "/Tsec_insulators.dat",
                          "/Tsec_span2_stats.csv", "/Tsec_insulators_stats.csv"}) {
        EXPECT_TRUE(std::filesystem::exists(dir + f)) << f;
    }
    std::ifstream ins(dir + "/Tsec_insulators.dat");
    std::getline(ins, header);
    EXPECT_EQ(header, "time t1_swing_deg t1_across_deg t1_tension t2_swing_deg t2_across_deg t2_tension");
    // the strings' statistics: a sample per step, the strings swinging with the +y wind
    std::ifstream st(dir + "/Tsec_insulators_stats.csv");
    std::getline(st, header);
    std::string row;
    std::getline(st, row);
    EXPECT_EQ(row.substr(0, row.find(',')), "3") << row;
    EXPECT_NE(row.find(",t1_swing_deg,"), std::string::npos) << row;
    const Real mean_swing = std::stod(row.substr(row.find("t1_swing_deg,") + 13));
    EXPECT_GT(mean_swing, 0.1) << row;
    EXPECT_FALSE(std::filesystem::exists(dir + "/separation.dat")) << "one line: no pairs to watch";
    ps.remove("towers");
    ps.remove("insulator_length");
    ps.remove("insulator_mass");
}

TEST(Conductors, TheClosestApproachOfTwoLinesIsFlaggedAgainstTheFlashoverDistance)
{
    const std::string dir = scratch("clash");
    // two parallel spans 6 m apart across a +y wind; the upwind one is lighter, swings further and
    // closes on the other
    set_inputs(dir, "Pa", true, {{300.0, 500.0, 30.0}}, {{600.0, 500.0, 30.0}});
    set_inputs(dir, "Pb", true, {{300.0, 506.0, 30.0}}, {{600.0, 506.0, 30.0}});
    amrex::ParmParse pp("conductors");
    pp.addarr("lines", std::vector<std::string>{"Pa", "Pb"});
    pp.add("flashover_distance", 5.5);
    amrex::ParmParse("conductors.Pa").add("mass_per_length", 1.0);
    amrex::ParmParse("conductors.Pb").add("mass_per_length", 3.0);
    Mesh m(false);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(nullptr, m.geom);
    ASSERT_EQ(c->pairs().size(), 1u);
    EXPECT_NEAR(c->separations()[0].distance, 6.0, 0.1) << "at rest the spans hang side by side";
    Real closest = 1.0e30;
    for (int s = 0; s < 40; ++s) {
        c->advance(0, 0.25 * s, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom);
        // the distance is the exact minimum between the two conductors at their current positions
        const auto exact = erf_conductors::closest_polylines(c->lines()[0]->conductor_path(), c->lines()[1]->conductor_path());
        ASSERT_NEAR(c->separations()[0].distance, exact.distance, tight * exact.distance) << "step " << s;
        closest = std::min(closest, c->separations()[0].distance);
    }
    EXPECT_LT(closest, 5.5) << "the lighter line must close on the heavier one";
    // separation.dat: the pair's columns, the flag off at rest and on once inside the flashover distance
    std::ifstream f(dir + "/separation.dat");
    std::string header;
    std::getline(f, header);
    EXPECT_EQ(header, "time Pa-Pb_distance Pa-Pb_x Pa-Pb_y Pa-Pb_z Pa-Pb_clash");
    std::vector<std::array<Real,6>> rows;
    std::array<Real,6> r{};
    while (f >> r[0] >> r[1] >> r[2] >> r[3] >> r[4] >> r[5]) { rows.push_back(r); }
    ASSERT_EQ(rows.size(), 41u) << "the initial row and one per step";
    for (const auto& row : rows) { EXPECT_EQ(row[5], (row[1] < 5.5) ? 1.0 : 0.0) << "t = " << row[0]; }
    EXPECT_EQ(rows.front()[5], 0.0);
    EXPECT_EQ(rows.back()[5], 1.0);
    EXPECT_GT(rows.back()[3], 500.0);
    EXPECT_LT(rows.back()[3], 506.0) << "the closest point lies between the lines";
    // the statistics: the fraction of samples in a clash is the mean of the flag
    const std::string stats = slurp(dir + "/separation_Pa-Pb_stats.csv");
    EXPECT_NE(stats.find(",clash,"), std::string::npos) << stats;
    pp.addarr("lines", std::vector<std::string>{});
    amrex::ParmParse("conductors.Pa").remove("mass_per_length");
    amrex::ParmParse("conductors.Pb").remove("mass_per_length");
}

TEST(Conductors, ARestartContinuesASectionItsStringsAndTheSeparationOfTheLines)
{
    const std::string dir = scratch("restart_circuit");
    // a section on strings and a single span beside its middle span
    set_inputs(dir, "Rs", true, {{100.0, 500.0, 30.0}}, {{1000.0, 500.0, 30.0}});
    // Rp runs along -x, so an older checkpoint's drag_y (its ERF-frame y drag) is not its drag across it
    set_inputs(dir, "Rp", true, {{700.0, 506.0, 30.0}}, {{400.0, 506.0, 30.0}});
    amrex::ParmParse pp("conductors");
    pp.addarr("lines", std::vector<std::string>{"Rs", "Rp"});
    amrex::ParmParse ps("conductors.Rs");
    ps.remove("length");
    ps.addarr("length", std::vector<Real>{301.5, 301.5, 301.5});
    ps.addarr("towers", std::vector<Real>{400.0, 500.0, 30.0, 700.0, 500.0, 30.0});
    ps.add("insulator_length", 2.5);
    ps.add("insulator_mass", 60.0);
    Mesh m(false);
    const double dt = 0.25;
    auto a = Conductors::create(0);
    ASSERT_TRUE(a);
    a->set_ground(nullptr, m.geom);
    int step = 0;
    for (; step < 4; ++step) { a->advance(0, dt * step, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom); }
    const std::string chk = dir + "/chk00004";
    std::filesystem::create_directories(chk);
    a->write_checkpoint(chk);
    for (; step < 7; ++step) { a->advance(0, dt * step, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom); }
    const auto files = std::vector<std::string>{"/separation.dat", "/separation_Rs-Rp_stats.csv", "/Rs_span2.dat",
                                                "/Rs_insulators.dat", "/Rs_insulators_stats.csv", "/Rp.dat"};
    std::vector<std::string> before;
    for (const auto& f : files) { before.push_back(slurp(dir + f)); ASSERT_FALSE(before.back().empty()) << f; }
    const Real sep = a->separations()[0].distance;

    auto b = Conductors::create(0);
    ASSERT_TRUE(b);
    b->set_ground(nullptr, m.geom, chk);
    ASSERT_TRUE(b->restored());
    EXPECT_LT(slurp(dir + "/separation.dat").size(), before[0].size()) << "the separation log ends at the checkpoint";
    for (step = 4; step < 7; ++step) { b->advance(0, dt * step, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom); }
    for (std::size_t i = 0; i < files.size(); ++i) { EXPECT_EQ(slurp(dir + files[i]), before[i]) << files[i]; }
    EXPECT_NEAR(b->separations()[0].distance, sep, roundoff * sep);
    EXPECT_NEAR(b->lines()[0]->insulator_swing(0), a->lines()[0]->insulator_swing(0), roundoff);
    {
        // a checkpoint listing drag_y in drag_normal's place: Rs's spans run along +x, where the two are one, so
        // their statistics continue; Rp runs along -x, so its start afresh
        int renamed = 0;
        for (const auto& e : std::filesystem::directory_iterator(chk + "/conductors")) {
            const std::string f = e.path().string();
            if (f.size() < 10 || f.compare(f.size() - 10, 10, "_stats.dat") != 0) { continue; }
            std::string text = slurp(f);
            const auto at = text.find("drag_normal ");
            if (at == std::string::npos) { continue; }
            text.replace(at, 12, "drag_y ");
            std::ofstream(f, std::ios::trunc) << text;
            ++renamed;
        }
        EXPECT_EQ(renamed, 4) << "Rs's three spans and Rp's";
        auto c = Conductors::create(0);
        EXPECT_TRUE(erf_gtest::abort_message([&] { c->set_ground(nullptr, m.geom, chk); }).empty());
        // its checkpoint still counts the samples taken before (afresh, it would count none)
        const std::string chk2 = dir + "/chk_renamed";
        std::filesystem::create_directories(chk2);
        c->write_checkpoint(chk2);
        int continued = 0, afresh = 0;
        for (const auto& e : std::filesystem::directory_iterator(chk2 + "/conductors")) {
            const std::string text = slurp(e.path().string());
            if (text.find("drag_normal ") == std::string::npos) { continue; }
            const bool along_minus_x = e.path().filename().string().rfind("Rp", 0) == 0;
            EXPECT_EQ(text.find("count = 0\n") != std::string::npos, along_minus_x) << e.path();
            ++(along_minus_x ? afresh : continued);
        }
        EXPECT_EQ(continued, 3);
        EXPECT_EQ(afresh, 1);
    }
    {
        // the explicit lengths stay, the strings grow 0.5 m: the nodes would restart where shorter strings held them
        ps.remove("insulator_length");
        ps.add("insulator_length", 3.0);
        auto c = Conductors::create(0);
        const std::string msg = erf_gtest::abort_message([&] { c->set_ground(nullptr, m.geom, chk); });
        EXPECT_NE(msg.find("line Rs's insulator_length was 2.5 and is 3"), std::string::npos) << msg;
    }
    pp.addarr("lines", std::vector<std::string>{});
    ps.remove("towers");
    ps.remove("insulator_length");
    ps.remove("insulator_mass");
}

TEST(Conductors, TransformersTakeThePullOfTheLinesEndingOnThemAndContinueAcrossARestart)
{
    const std::string dir = scratch("transformers");
    // a span from T1 to T2 on the ramp, 12 m above the ground at its ends (6 m over the tops), and
    // a second line dead-ended on open ground 2 m beside T2 and 3 m higher than its top; a +y wind
    set_inputs(dir, "Lt", true, {{302.0, 500.0, 12.0}}, {{598.0, 500.0, 12.0}});
    set_inputs(dir, "Lo", true, {{606.0, 500.0, 9.0}}, {{606.0, 800.0, 9.0}});
    amrex::ParmParse("conductors.Lt").remove("length");
    amrex::ParmParse("conductors.Lt").add("length", 297.5);
    amrex::ParmParse pp("conductors");
    pp.addarr("lines", std::vector<std::string>{"Lt", "Lo"});
    pp.addarr("transformers", std::vector<std::string>{"T1", "T2"});
    for (const char* t : {"T1", "T2"}) {
        amrex::ParmParse pt(std::string("conductors.") + t);
        pt.addarr("position", std::vector<Real>{t[1] == '1' ? Real(300.0) : Real(600.0), Real(500.0)});
        pt.addarr("size", std::vector<Real>{8.0, 5.0, 6.0});
    }
    amrex::ParmParse("conductors.T1").add("allowable_force", 1.0);   // any pull is over it
    Mesh m(true);
    const double dt = 0.25;
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(m.znd.get(), m.geom);
    ASSERT_EQ(c->transformers().size(), 2u);
    const auto& T1 = c->transformers()[0];
    const auto& T2 = c->transformers()[1];
    // the boxes stand on the ramp under their centres; Lt ends on both, Lo on neither
    EXPECT_NEAR(T1.base()[2], m.h(300.0, 500.0), roundoff * 1000);
    EXPECT_NEAR(T2.box_hi()[2], m.h(600.0, 500.0) + 6.0, roundoff * 1000);
    ASSERT_EQ(T1.ends().size(), 1u);
    ASSERT_EQ(T2.ends().size(), 1u);
    EXPECT_EQ(T1.ends()[0].end, 0);
    EXPECT_EQ(T2.ends()[0].end, 1);
    std::ifstream g(dir + "/ground.dat");
    std::string row, name, label;
    int boxes = 0;
    while (std::getline(g, row)) {
        std::istringstream ls(row);
        Real x, y, ground, top;
        if (ls >> name >> label >> x >> y >> ground >> top && label == "transformer") {
            ++boxes;
            EXPECT_NEAR(ground, m.h(x, y), postol) << name;
            EXPECT_NEAR(top, m.h(x, y) + 6.0, postol) << name;
        }
    }
    EXPECT_EQ(boxes, 2);

    auto check = [&](const Conductors& cc, const char* when) {
        const auto& span = *cc.lines()[0];
        const auto& L1 = cc.transformer_loads()[0];
        const auto& L2 = cc.transformer_loads()[1];
        const auto fa = span.end_force(0);
        const auto fb = span.end_force(1);
        for (int d = 0; d < 3; ++d) {
            EXPECT_EQ(L1.force[d], fa[d]) << when << " T1 carries Lt's end_a alone, component " << d;
            EXPECT_EQ(L2.force[d], fb[d]) << when << " T2 carries Lt's end_b alone, component " << d;
        }
        // the moment about the base centre of the pull at the end
        const auto& e = span.inputs().end_a;
        const auto b = cc.transformers()[0].base();
        const Real rx = e[0] - b[0], ry = e[1] - b[1], rz = e[2] - b[2];
        EXPECT_NEAR(L1.moment[0], ry * fa[2] - rz * fa[1], roundoff * 1.0e3 * std::abs(L1.moment[1])) << when;
        EXPECT_NEAR(L1.moment[1], rz * fa[0] - rx * fa[2], roundoff * 1.0e3 * std::abs(L1.moment[1])) << when;
        EXPECT_GT(L1.horizontal_force, 0.0);
        EXPECT_TRUE(L1.over_allowable) << when;
        EXPECT_FALSE(L2.over_allowable) << when << " no allowable on T2";
        // the closest conductor to each box: the exact distance over every line
        for (std::size_t t = 0; t < 2; ++t) {
            const auto& tr = cc.transformers()[t];
            Real best = 1.0e30;
            for (const auto& sp : cc.lines()) {
                best = std::min(best, erf_conductors::closest_polyline_box(sp->conductor_path(), tr.box_lo(), tr.box_hi()).distance);
            }
            EXPECT_EQ(cc.transformer_clearances()[t].distance, best) << when << " " << tr.name();
        }
    };
    check(*c, "at rest");
    // Lo's end is about 3.9 m from T2's top edge, nearer than Lt's 6 m standoff; only Lt comes near T1
    EXPECT_GT(c->transformer_clearances()[1].distance, 3.0);
    EXPECT_LT(c->transformer_clearances()[1].distance, 4.5);
    EXPECT_LE(c->transformer_clearances()[0].distance, 6.0 + 0.2);

    int step = 0;
    for (; step < 3; ++step) { c->advance(0, dt * step, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); }
    check(*c, "blown out");
    const std::string chk = dir + "/chk00003";
    std::filesystem::create_directories(chk);
    c->write_checkpoint(chk);
    for (; step < 6; ++step) { c->advance(0, dt * step, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); }
    const std::string log = slurp(dir + "/transformers.dat");
    const std::string stats1 = slurp(dir + "/transformer_T1_stats.csv");
    const std::string stats2 = slurp(dir + "/transformer_T2_stats.csv");
    // the log: a row at the start and one per step, T1 over its allowable throughout
    std::istringstream lg(log);
    std::string header;
    std::getline(lg, header);
    EXPECT_EQ(header.rfind("time T1_Fx T1_Fy T1_Fz T1_Fh T1_Mx T1_My T1_Mh T1_over T1_clearance T1_clash T2_Fx ", 0), 0u) << header;
    int rows = 0;
    while (std::getline(lg, row)) {
        std::istringstream ls(row);
        std::vector<Real> v;
        Real x;
        while (ls >> x) { v.push_back(x); }
        ASSERT_EQ(v.size(), 21u) << row;
        EXPECT_EQ(v[8], 1.0) << "T1_over";
        EXPECT_EQ(v[18], 0.0) << "T2_over";
        EXPECT_EQ(v[20], v[19] < 1.0 ? 1.0 : 0.0) << "T2_clash against the 1 m flashover distance";
        ++rows;
    }
    EXPECT_EQ(rows, 7);
    EXPECT_NE(stats1.find(",overturning_moment,"), std::string::npos) << stats1;

    // restarted from the checkpoint: the same log, statistics and loads
    auto r = Conductors::create(0);
    ASSERT_TRUE(r);
    r->set_ground(m.znd.get(), m.geom, chk);
    ASSERT_TRUE(r->restored());
    for (step = 3; step < 6; ++step) { r->advance(0, dt * step, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); }
    EXPECT_EQ(slurp(dir + "/transformers.dat"), log);
    EXPECT_EQ(slurp(dir + "/transformer_T1_stats.csv"), stats1);
    EXPECT_EQ(slurp(dir + "/transformer_T2_stats.csv"), stats2);
    for (int d = 0; d < 3; ++d) { EXPECT_EQ(r->transformer_loads()[0].force[d], c->transformer_loads()[0].force[d]); }

    pp.addarr("lines", std::vector<std::string>{});
    pp.addarr("transformers", std::vector<std::string>{});
    amrex::ParmParse("conductors.T1").remove("allowable_force");
}

TEST(Conductors, AStringingTensionSetsTheLengthsFromTheChordsOnTheTerrain)
{
    const std::string dir = scratch("stringing");
    // a section up the ramp, clamped at the towers, with spans of 300, 150 and 450 m strung to 20 kN
    set_inputs(dir, "Tstr", true, {{100.0, 500.0, 30.0}}, {{1000.0, 500.0, 30.0}});
    amrex::ParmParse ps("conductors.Tstr");
    ps.remove("length");
    ps.add("stringing_tension", 2.0e4);
    ps.addarr("towers", std::vector<Real>{400.0, 500.0, 30.0, 550.0, 500.0, 30.0});
    Mesh m(true);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(m.znd.get(), m.geom);
    const auto& line = *c->lines().front();
    const auto& s = line.inputs();
    ASSERT_EQ(s.lengths.size(), 3u);
    const std::array<Real,3> h{{300.0, 150.0, 450.0}};
    for (int k = 0; k < 3; ++k) {
        // the chord between the attachments on the ramp, which rises 5 m per 100 m
        const Real c_on_terrain = h[k] * std::sqrt(Real(1.0) + m.slope_x * m.slope_x);
        EXPECT_NEAR(s.chord(k), c_on_terrain, roundoff * 1.0e4) << k;
    }
    // every span pulls on its towers with about the stringing tension along the line: the real
    // MoorDyn within the parabola's few per cent, the stub (no stretch) within its sag error
    const Real tol = erf_moordyn::is_stub() ? Real(0.25) : Real(0.04);
    for (int k = 0; k < 3; ++k) {
        const Real T = line.tension_a(k);
        EXPECT_GT(T / 2.0e4, 1.0 - tol) << "span " << k;
        EXPECT_LT(T / 2.0e4, 1.0 + 2.0 * tol) << "span " << k;
    }
    const auto fa = line.end_force(0);
    EXPECT_NEAR(std::sqrt(fa[0] * fa[0] + fa[1] * fa[1]) / 2.0e4, 1.0, tol) << "the dead end takes the stringing tension";
    ps.remove("stringing_tension");
    ps.remove("towers");
}

namespace {
// a section of three spans along x on the ramp with lattice towers at x = 400 and 700
void set_towered_section (const std::string& dir, const std::string& name)
{
    set_inputs(dir, name, true, {{100.0, 500.0, 30.0}}, {{1000.0, 500.0, 30.0}});
    amrex::ParmParse ps("conductors." + name);
    ps.remove("length");
    ps.addarr("length", std::vector<Real>{301.5, 301.5, 301.5});
    ps.addarr("towers", std::vector<Real>{400.0, 500.0, 30.0, 700.0, 500.0, 30.0});
    ps.add("tower_type", std::string("lat"));
    amrex::ParmParse pp("conductors");
    pp.addarr("tower_types", std::vector<std::string>{"lat"});
    amrex::ParmParse pt("conductors.lat");
    for (const char* k : {"base_width", "top_width", "solidity", "arm_length", "arm_depth"}) { pt.remove(k); }
    pt.add("base_width", 6.0);
    pt.add("top_width", 1.5);
    pt.add("solidity", 0.2);
    pt.add("arm_length", 12.0);
    pt.add("arm_depth", 1.2);
    for (const char* k : {"weight", "allowable_uplift", "allowable_compression"}) { pt.remove(k); }
    pt.add("weight", 9.0e4);
}
void clear_towered_section (const std::string& name)
{
    amrex::ParmParse pt("conductors.lat");
    for (const char* k : {"frequency", "damping_ratio"}) { pt.remove(k); }
    amrex::ParmParse ps("conductors." + name);
    ps.remove("towers"); ps.remove("tower_type");
    amrex::ParmParse pp("conductors");
    pp.addarr("tower_types", std::vector<std::string>{});
    pp.remove("drag_on_flow");
    pp.addarr("lines", std::vector<std::string>{});
}
} // namespace

TEST(Conductors, LatticeTowersStandAtTheSuspensionPointsAndCarryTheWindsDrag)
{
    const std::string dir = scratch("towers");
    set_towered_section(dir, "Tw");
    Mesh m(true);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(m.znd.get(), m.geom);
    ASSERT_EQ(c->towers().size(), 2u);
    for (int k = 0; k < 2; ++k) {
        const auto& tw = c->towers()[static_cast<std::size_t>(k)];
        const Real x = 400.0 + 300.0 * k;
        EXPECT_EQ(tw.name(), "Tw_t" + std::to_string(k + 1));
        EXPECT_NEAR(tw.base()[0], x, (std::is_same<Real, float>::value) ? Real(1.0e-4) : Real(1.0e-9));
        EXPECT_NEAR(tw.base()[2], m.h(x, 500.0), roundoff * 1000) << "the base stands on the ramp";
        EXPECT_NEAR(tw.arm_height(), 30.0, roundoff * 1000) << "the cross-arm at the conductor's height";
        EXPECT_NEAR(std::abs(tw.across()[1]), 1.0, tight) << "the cross-arm across a line along x";
    }
    // the prescribed +y wind runs along the cross-arms: only the bodies carry it; the body's drag by
    // hand, q Cf phi (b0 + b1)/2 H (N), with Cf = 4 phi^2 - 5.9 phi + 4 for the solidity phi = 0.2
    const Real U = 10.0, q = 0.5 * 1.2 * U * U, cf = 4.0 * 0.04 - 5.9 * 0.2 + 4.0;
    const Real body = q * cf * 0.2 * 0.5 * (6.0 + 1.5) * 30.0;
    // the line's pull on each tower at the start of the last step, which the towers carry over it
    std::vector<std::array<Real,3>> pull;
    for (int s = 0; s < 3; ++s) {
        pull.clear();
        for (int t = 0; t < 2; ++t) { pull.push_back(c->lines()[0]->tower_force(t)); }
        c->advance(0, 0.25 * s, 0.25, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom);
    }
    for (const auto& tw : c->towers()) {
        const auto F = tw.total_force();
        EXPECT_NEAR(F[1], body, 1.0e-6 * body) << tw.name();
        EXPECT_NEAR(F[0], 0.0, roundoff * body) << tw.name();
        // M_x = -z F_y: the body's q cf phi int w(z) z dz, short by the midpoint sum's H^3 / (12 n^2) on z^2
        const Real H = 30.0, n = 10.0;
        const Real mom = q * cf * 0.2 * (6.0 * H * H / 2.0 + (1.5 - 6.0) / H * (H * H * H / 3.0 - H * H * H / (12.0 * n * n)));
        EXPECT_NEAR(tw.base_moment()[0], -mom, 1.0e-6 * mom) << tw.name();
    }
    // each tower carries its line's pull where the line hangs from it, and its footings the lot
    for (std::size_t t = 0; t < 2; ++t) {
        const auto& tw = c->towers()[t];
        const auto& F = pull[t];
        for (int d = 0; d < 3; ++d) { EXPECT_EQ(tw.line_force()[d], F[d]) << tw.name() << " " << d; }
        EXPECT_LT(F[2], 0.0) << "the clamped conductor weighs on the tower";
        const auto L = tw.foundation();
        EXPECT_NEAR(L.vertical, 9.0e4 - F[2] - tw.total_force()[2], roundoff * L.vertical);
        EXPECT_NEAR(L.legs[0] + L.legs[1] + L.legs[2] + L.legs[3], L.vertical, roundoff * L.vertical);
        EXPECT_NEAR(L.shear, std::hypot(F[0] + tw.total_force()[0], F[1] + tw.total_force()[1]), roundoff * L.shear);
    }
    // towers.dat: a row at the start of every step with every tower's drag, line pull and foundation load
    std::ifstream f(dir + "/towers.dat");
    std::string header, row;
    std::getline(f, header);
    EXPECT_EQ(header.rfind("time Tw_t1_drag_Fx Tw_t1_drag_Fy Tw_t1_drag_Fz Tw_t1_line_Fx Tw_t1_line_Fy Tw_t1_line_Fz Tw_t1_shear "
                           "Tw_t1_overturning Tw_t1_vertical Tw_t1_max_compression Tw_t1_max_uplift Tw_t1_over Tw_t2_drag_Fx", 0), 0u)
        << header;
    int rows = 0;
    while (std::getline(f, row)) {
        std::istringstream ls(row);
        std::vector<Real> v;
        Real x;
        while (ls >> x) { v.push_back(x); }
        ASSERT_EQ(v.size(), 25u);
        EXPECT_NEAR(v[2], body, 1.0e-6 * body);
        EXPECT_EQ(v[12], 0.0) << "no allowable given";
        ++rows;
    }
    EXPECT_EQ(rows, 3);
    for (const char* column : {",drag_h,", ",line_h,", ",shear,", ",overturning,", ",max_compression,", ",max_uplift,", ",over_allowable,"}) {
        EXPECT_NE(slurp(dir + "/tower_Tw_t1_stats.csv").find(column), std::string::npos) << column;
    }
    clear_towered_section("Tw");
}

TEST(Conductors, TheTowersDragGoesIntoTheFlowWithTheLinesAndSurvivesARestart)
{
    const std::string dir = scratch("towers_flow");
    set_towered_section(dir, "Tf");
    amrex::ParmParse pp("conductors");
    pp.add("drag_on_flow", true);
    pp.remove("prescribed_velocity");
    Mesh m(true);
    m.fill_crosswind();
    const double dt = 0.25;
    auto a = Conductors::create(0);
    ASSERT_TRUE(a);
    a->set_ground(m.znd.get(), m.geom);
    int step = 0;
    for (; step < 3; ++step) { a->advance(0, dt * step, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); }
    // the source integrates to minus the drag on the lines and on the towers
    std::array<Real,3> drag{{0.0, 0.0, 0.0}};
    for (const auto& s : a->lines()) { const auto f = s->total_drag(); for (int d = 0; d < 3; ++d) { drag[d] += f[d]; } }
    Real towers_y = 0.0;
    for (const auto& t : a->towers()) { const auto f = t.total_force(); towers_y += f[1]; for (int d = 0; d < 3; ++d) { drag[d] += f[d]; } }
    EXPECT_GT(towers_y, 0.0);
    for (int d = 0; d < 3; ++d) { EXPECT_NEAR(a->source_integral()[d], -drag[d], 1.0e-6 * std::abs(drag[1]) + roundoff) << d; }
    const std::string chk = dir + "/chk00003";
    std::filesystem::create_directories(chk);
    a->write_checkpoint(chk);
    const auto src = a->source_integral();
    for (; step < 6; ++step) { a->advance(0, dt * step, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); }
    const std::string log = slurp(dir + "/towers.dat");
    const std::string stats = slurp(dir + "/tower_Tf_t2_stats.csv");

    auto b = Conductors::create(0);
    ASSERT_TRUE(b);
    b->set_ground(m.znd.get(), m.geom, chk);
    ASSERT_TRUE(b->restored());
    // the restored towers' drag is the checkpointed step's, and so is the source spread from it
    b->restore_sources(0, m.u, m.znd.get(), nullptr, m.geom);
    for (int d = 0; d < 3; ++d) { EXPECT_EQ(b->source_integral()[d], src[d]) << d; }
    for (step = 3; step < 6; ++step) { b->advance(0, dt * step, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); }
    EXPECT_EQ(slurp(dir + "/towers.dat"), log);
    EXPECT_EQ(slurp(dir + "/tower_Tf_t2_stats.csv"), stats);
    // the checkpoint records the tower type's forms: a restart that changes one is refused, naming the key
    amrex::ParmParse("conductors.lat").add("arm_outside_shaft", false);
    std::string msg = erf_gtest::abort_message([&] {
        auto c = Conductors::create(0);
        c->set_ground(m.znd.get(), m.geom, chk);
    });
    EXPECT_NE(msg.find("conductors.lat.arm_outside_shaft = false"), std::string::npos) << msg;
    amrex::ParmParse("conductors.lat").remove("arm_outside_shaft");
    // a checkpoint that records the forms but not a tower's, and an older one without the record, written with the
    // earlier forms: both refused, the older naming the keys to set to false
    auto restart_with = [&] (const std::string& state) {
        std::ofstream(chk + "/conductors/state", std::ios::trunc) << state;
        return erf_gtest::abort_message([&] {
            auto c = Conductors::create(0);
            c->set_ground(m.znd.get(), m.geom, chk);
        });
    };
    const std::string state = slurp(chk + "/conductors/state");
    ASSERT_NE(state.find("tower_forms_version 2\n"), std::string::npos) << state;
    std::string without, older, rectangle;
    {
        std::istringstream lines(state);
        for (std::string l; std::getline(lines, l);) {
            if (l.rfind("tower_forms ", 0) != 0) { without += l + "\n"; }
            if (l.rfind("tower_forms", 0) != 0) { older += l + "\n"; }
            rectangle += (l == "tower_forms_version 2" ? std::string("tower_forms_version 1") : l) + "\n";
        }
    }
    msg = restart_with(without);
    EXPECT_NE(msg.find("holds no forms of tower Tf_t"), std::string::npos) << msg;
    // version 1 recorded the outside layout with the rectangular cross-arm: its node forces are not these nodes'
    msg = restart_with(rectangle);
    EXPECT_NE(msg.find("before the cross-arm's face tapered to the tips"), std::string::npos) << msg;
    EXPECT_NE(msg.find("conductors.lat.arm_outside_shaft = true"), std::string::npos) << msg;
    msg = restart_with(older);
    EXPECT_NE(msg.find("written before the tower forms were recorded"), std::string::npos) << msg;
    EXPECT_NE(msg.find("conductors.lat.diagonal_wind_factor"), std::string::npos) << msg;
    EXPECT_EQ(msg.find("angle_principal_axes"), std::string::npos) << "a tower without a generated frame: " << msg;
    // with those forms set to false the older checkpoint restarts
    amrex::ParmParse lat("conductors.lat");
    lat.add("diagonal_wind_factor", false);
    lat.add("arm_outside_shaft", false);
    EXPECT_TRUE(restart_with(older).empty()) << restart_with(older);
    lat.remove("diagonal_wind_factor");
    lat.remove("arm_outside_shaft");
    std::ofstream(chk + "/conductors/state", std::ios::trunc) << state;
    clear_towered_section("Tf");
}

namespace {
// the section's towers bend at 2 Hz with the damping ratio given
void make_towers_move (Real damping)
{
    amrex::ParmParse pt("conductors.lat");
    for (const char* k : {"frequency", "damping_ratio"}) { pt.remove(k); }
    pt.add("frequency", 2.0);
    pt.add("damping_ratio", damping);
}
} // namespace

TEST(Conductors, EachLineAndTowerSpreadsItsDragWithItsOwnWidth)
{
    // a section along x at y = 500 over two lattice towers, with its own width (1.2 cells); its towers' type with
    // another (0.8 cells); a single span along x at y = 250 without one, which takes conductors.epsilon (1.6 cells).
    // A prescribed crosswind: the drag of every node is along +y, the towers' bodies too
    const std::string dir = scratch("widths");
    set_towered_section(dir, "Tew");
    amrex::ParmParse pp("conductors");
    pp.add("drag_on_flow", true);
    pp.add("epsilon", 1.6);
    pp.addarr("lines", std::vector<std::string>{"Tew", "Tew2"});
    amrex::ParmParse("conductors.Tew").add("epsilon", 1.2);
    amrex::ParmParse("conductors.lat").add("epsilon", 0.8);
    {
        amrex::ParmParse ps("conductors.Tew2");
        ps.addarr("end_a", std::vector<Real>{300.0, 250.0, 30.0});
        ps.addarr("end_b", std::vector<Real>{600.0, 250.0, 30.0});
        ps.add("length", 301.5);
        ps.add("diameter", 0.0281);
        ps.add("mass_per_length", 1.628);
        ps.add("axial_stiffness", 3.0e7);
    }
    // ParmParse is global to the test binary: leave no width behind, whatever the test's outcome
    struct Cleanup {
        ~Cleanup ()
        {
            amrex::ParmParse("conductors.Tew").remove("epsilon");
            amrex::ParmParse("conductors.lat").remove("epsilon");
            amrex::ParmParse("conductors").remove("epsilon");
            clear_towered_section("Tew");
        }
    } cleanup;
    Mesh m(true);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(m.znd.get(), m.geom);
    ASSERT_EQ(c->lines().size(), 2u);
    ASSERT_EQ(c->towers().size(), 2u);
    c->advance(0, 0.0, 0.25, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom);

    // every node's position and force on the air, and the width it should be spread with (cells)
    std::vector<Real> pos, force, cells;
    for (const auto& span : c->lines()) {
        const std::vector<Real> p = span->node_positions();
        pos.insert(pos.end(), p.begin(), p.end());
        const Real e = (span->name() == "Tew") ? Real(1.2) : Real(1.6);
        for (unsigned n = 0; n < span->num_nodes(); ++n) {
            const auto f = span->node_drag(n);
            force.insert(force.end(), {-f[0], -f[1], -f[2]});
            cells.push_back(e);
        }
    }
    std::size_t tower_nodes = 0;
    for (const auto& t : c->towers()) {
        const auto nodes = t.current_nodes();
        for (std::size_t n = 0; n < nodes.size(); ++n) {
            pos.insert(pos.end(), {nodes[n].pos[0], nodes[n].pos[1], nodes[n].pos[2]});
            force.insert(force.end(), {-t.loads()[3*n], -t.loads()[3*n+1], -t.loads()[3*n+2]});
            cells.push_back(Real(0.8));
            ++tower_nodes;
        }
    }
    ASSERT_GT(tower_nodes, 0u);
    const Real dx = m.geom.CellSize(0);
    std::vector<Real> eps(cells.size());
    for (std::size_t p = 0; p < cells.size(); ++p) { eps[p] = cells[p] * dx; }
    auto sources = [&] (const std::vector<Real>& widths) {
        std::array<amrex::MultiFab,3> s{{amrex::MultiFab(amrex::convert(m.ba, amrex::IntVect(1,0,0)), m.dm, 1, 0),
                                         amrex::MultiFab(amrex::convert(m.ba, amrex::IntVect(0,1,0)), m.dm, 1, 0),
                                         amrex::MultiFab(amrex::convert(m.ba, amrex::IntVect(0,0,1)), m.dm, 1, 0)}};
        erf_actuator::spread_forces(pos, force, widths, m.znd.get(), nullptr, m.geom, s[0], s[1], s[2]);
        return s;
    };
    const auto want = sources(eps);
    amrex::MultiFab sx(amrex::convert(m.ba, amrex::IntVect(1,0,0)), m.dm, 1, 0);
    amrex::MultiFab sy(amrex::convert(m.ba, amrex::IntVect(0,1,0)), m.dm, 1, 0);
    amrex::MultiFab sz(amrex::convert(m.ba, amrex::IntVect(0,0,1)), m.dm, 1, 0);
    sx.setVal(0.0); sy.setVal(0.0); sz.setVal(0.0);
    c->add_momentum_sources(0, sx, sy, sz);
    const Real scale = want[1].norm0();
    ASSERT_GT(scale, 0.0) << "the crosswind's drag must reach the flow";
    // the same widths, points and forces in the same order: the same sources to the last place
    amrex::MultiFab diff(sy.boxArray(), sy.DistributionMap(), 1, 0);
    amrex::MultiFab::LinComb(diff, Real(1.0), sy, 0, Real(-1.0), want[1], 0, 0, 1, 0);
    EXPECT_LE(diff.norm0(), tight * scale) << "the y source differs from spreading each node with its own width";
    for (const auto& [got, ref] : {std::make_pair(&sx, &want[0]), std::make_pair(&sz, &want[2])}) {
        amrex::MultiFab d(got->boxArray(), got->DistributionMap(), 1, 0);
        amrex::MultiFab::LinComb(d, Real(1.0), *got, 0, Real(-1.0), *ref, 0, 0, 1, 0);
        EXPECT_LE(d.norm0(), tight * std::max(scale, ref->norm0()));
    }
    // and not those of one width for all: each choice a mutant would make moves the source by much more
    for (const Real one : {Real(0.8), Real(1.2), Real(1.6)}) {
        const auto other = sources(std::vector<Real>(eps.size(), one * dx));
        amrex::MultiFab d(sy.boxArray(), sy.DistributionMap(), 1, 0);
        amrex::MultiFab::LinComb(d, Real(1.0), sy, 0, Real(-1.0), other[1], 0, 0, 1, 0);
        EXPECT_GT(d.norm0(), Real(0.01) * scale) << "one width of " << one << " cells for every node";
    }
    // the swapped mapping, the line's width on the towers and the type's on the line, differs too
    {
        std::vector<Real> swapped(eps);
        for (std::size_t p = 0; p < swapped.size(); ++p) {
            if (cells[p] == Real(1.2)) { swapped[p] = Real(0.8) * dx; } else if (cells[p] == Real(0.8)) { swapped[p] = Real(1.2) * dx; }
        }
        const auto other = sources(swapped);
        amrex::MultiFab d(sy.boxArray(), sy.DistributionMap(), 1, 0);
        amrex::MultiFab::LinComb(d, Real(1.0), sy, 0, Real(-1.0), other[1], 0, 0, 1, 0);
        EXPECT_GT(d.norm0(), Real(0.01) * scale) << "the line's and its towers' widths swapped";
    }
    // the drag's support, which ERF's structure coupling checks against the finer levels: one box per spread
    // point, in the same order, the point grown by its own line's or tower type's reach (3 epsilon cells, at
    // least a cell). A support built with one width for all, or without the kernel's reach, would let a finer
    // level overlap the drag unseen and remove it when averaging down
    {
        const auto& sup = c->support();
        ASSERT_EQ(sup.size(), cells.size()) << "one box per spread point";
        const Real cell = std::max(m.geom.CellSize(0), m.geom.CellSize(1));
        for (std::size_t p = 0; p < cells.size(); ++p) {
            const Real r = std::max(cell, Real(3.0) * cells[p] * dx);
            for (int d = 0; d < 3; ++d) {
                EXPECT_NEAR(sup[p].lo(d), pos[3*p+d] - r, tight * r) << "point " << p << " dir " << d;
                EXPECT_NEAR(sup[p].hi(d), pos[3*p+d] + r, tight * r) << "point " << p << " dir " << d;
            }
        }
    }
    // the source still integrates to minus the drag on the lines and the towers
    std::array<Real,3> drag{{0.0, 0.0, 0.0}};
    for (const auto& s : c->lines()) { const auto f = s->total_drag(); for (int d = 0; d < 3; ++d) { drag[d] += f[d]; } }
    for (const auto& t : c->towers()) { const auto f = t.total_force(); for (int d = 0; d < 3; ++d) { drag[d] += f[d]; } }
    for (int d = 0; d < 3; ++d) { EXPECT_NEAR(c->source_integral()[d], -drag[d], 1.0e-6 * std::abs(drag[1]) + roundoff) << d; }
}

TEST(Conductors, EachTowerIsCoveredWithItsTypesReach)
{
    // a patch of the level from y = 100 to 900 m (cells 1 to 8 of 100 m) under a section along x at y = 500 over two
    // lattice towers: a width of 1 cell reaches 300 m, inside it; 1.5 cells reach 450 m, beyond it. The line and its towers
    // each with their own width; the check names the line or the tower that outreaches the patch. The flow is sampled
    // (tower_wind's check) for the towers, prescribed (spread_drag's check) for the line
    const std::string dir = scratch("tower_reach");
    set_towered_section(dir, "Ttr");
    amrex::ParmParse pp("conductors");
    pp.add("drag_on_flow", true);
    struct Cleanup {
        ~Cleanup ()
        {
            amrex::ParmParse("conductors.Ttr").remove("epsilon");
            amrex::ParmParse("conductors.lat").remove("epsilon");
            clear_towered_section("Ttr");
        }
    } cleanup;
    Mesh m(false);
    const amrex::BoxArray ba(amrex::Box(amrex::IntVect(0, 1, 0), amrex::IntVect(m.nx - 1, 8, m.nz - 1)));
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab u(amrex::convert(ba, amrex::IntVect(1,0,0)), dm, 1, 1), v(amrex::convert(ba, amrex::IntVect(0,1,0)), dm, 1, 1),
                    w(amrex::convert(ba, amrex::IntVect(0,0,1)), dm, 1, 1);
    u.setVal(0.0); v.setVal(10.0); w.setVal(0.0);
    auto run = [&] (Real line, Real tower, bool sampled) {
        amrex::ParmParse("conductors.Ttr").remove("epsilon");
        amrex::ParmParse("conductors.lat").remove("epsilon");
        amrex::ParmParse("conductors.Ttr").add("epsilon", line);
        amrex::ParmParse("conductors.lat").add("epsilon", tower);
        pp.remove("prescribed_velocity");
        if (!sampled) { pp.addarr("prescribed_velocity", std::vector<Real>{0.0, 10.0, 0.0}); }
        return erf_gtest::abort_message([&] {
            auto c = Conductors::create(0);
            c->set_ground(nullptr, m.geom);
            c->advance(0, 0.0, 0.2, u, v, w, nullptr, nullptr, m.geom);
        });
    };
    for (const bool sampled : {true, false}) {
        SCOPED_TRACE(sampled ? "the flow sampled" : "a prescribed wind");
        std::string msg = run(1.0, 1.0, sampled);
        EXPECT_TRUE(msg.empty()) << msg;
        // the towers' width outreaches the patch, the line's does not
        msg = run(1.0, 1.5, sampled);
        EXPECT_NE(msg.find("a node of tower Ttr_t1, is not covered"), std::string::npos) << msg;
        // and the other way round (with the flow sampled the line's own sampling check, in wind_at, stops it first)
        msg = run(1.5, 1.0, sampled);
        EXPECT_NE(msg.find(sampled ? "conductors.Ttr: the point" : "a node of line Ttr, is not covered"), std::string::npos) << msg;
        EXPECT_EQ(msg.find("tower Ttr_t"), std::string::npos) << "the towers fit: " << msg;
    }
}

TEST(Conductors, MovingTowersSettleWhereTheirStiffnessBalancesTheWindAndTheLine)
{
    const double dt = 0.25;
    // the same section with towers that stand still, for the line's pull on rigid towers
    const std::string rdir = scratch("still_towers");
    set_towered_section(rdir, "Ms");
    Mesh m(true);
    std::vector<std::array<Real,3>> rigid_pull;
    {
        auto c = Conductors::create(0);
        ASSERT_TRUE(c);
        c->set_ground(m.znd.get(), m.geom);
        for (int s = 0; s < 40; ++s) { c->advance(0, dt * s, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); }
        for (const auto& tw : c->towers()) { rigid_pull.push_back(tw.line_force()); }
        EXPECT_EQ(c->tower_models()[0], nullptr);
    }
    clear_towered_section("Ms");
    const std::string dir = scratch("moving_towers");
    set_towered_section(dir, "Mv");
    make_towers_move(0.3);   // heavily damped, so that they settle in a few periods
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(m.znd.get(), m.geom);
    ASSERT_TRUE(c->lines()[0]->towers_move());
    for (int s = 0; s < 40; ++s) { c->advance(0, dt * s, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); }
    for (std::size_t t = 0; t < 2; ++t) {
        const auto& tw = c->towers()[t];
        const auto* model = dynamic_cast<const erf_towers::OneModeTower*>(c->tower_models()[t].get());
        ASSERT_NE(model, nullptr);
        // at rest, the stiffness carries the members' drag by the mode shape and the line's pull
        std::array<double,2> Q{{tw.line_force()[0], tw.line_force()[1]}};
        for (std::size_t i = 0; i < tw.nodes().size(); ++i) {
            for (int d = 0; d < 2; ++d) { Q[d] += model->mode_shape(i) * tw.loads()[3*i+d]; }
        }
        // (the real MoorDyn's conductors still swing slowly after 10 s, lightly damped, and the tower
        // follows their pull, hence the 2 % tolerance below)
        // the cross-arm's displacement at its height (its drag nodes, half its face lower, move a little less)
        std::array<Real,3> x{};
        ASSERT_TRUE(model->arm_centre_displacement(x));
        EXPECT_GT(x[1], 1.0e-3) << tw.name() << " leans with the +y wind";
        EXPECT_NEAR(x[1] * model->stiffness() / Q[1], 1.0, 0.02) << tw.name();
        EXPECT_NEAR(x[0] * model->stiffness(), Q[0], 0.02 * Q[1]) << tw.name();
        EXPECT_LT(std::abs(model->v()[1]), 0.02 * 2.0 * 3.14159265 * model->frequency() * x[1]) << tw.name() << " has settled";
        RecordProperty(tw.name() + "_lean_m", std::to_string(x[1]));
        // a few millimetres of lean barely change the line's pull from the rigid towers'
        const auto& Fr = rigid_pull[t];
        for (int d = 0; d < 3; ++d) {
            EXPECT_NEAR(tw.line_force()[d], Fr[d], 0.02 * std::abs(Fr[2])) << tw.name() << " dir " << d;
        }
        // the line hangs from the cross-arm where the tower has taken it
        const auto p = c->lines()[0]->node_position(c->lines()[0]->span_first_node(static_cast<int>(t) + 1));
        const Real ptol = (std::is_same<Real, float>::value) ? Real(1.0e-4) : Real(1.0e-6);   // a Real's spacing at 500 m
        EXPECT_NEAR(p[0], 400.0 + 300.0 * t + x[0], ptol);
        EXPECT_NEAR(p[1], 500.0 + x[1], ptol);
    }
    // towers.dat carries each moving tower's cross-arm displacement, and the statistics its size
    std::ifstream f(dir + "/towers.dat");
    std::string header, row, last;
    std::getline(f, header);
    EXPECT_NE(header.find("Mv_t1_over Mv_t1_arm_dx Mv_t1_arm_dy Mv_t2_drag_Fx"), std::string::npos) << header;
    while (std::getline(f, row)) { last = row; }
    std::istringstream ls(last);
    std::vector<Real> v;
    Real x;
    while (ls >> x) { v.push_back(x); }
    ASSERT_EQ(v.size(), 29u);
    EXPECT_GT(v[14], 1.0e-3) << "the row's cross-arm displacement across";
    EXPECT_NE(slurp(dir + "/tower_Mv_t1_stats.csv").find(",arm_displacement"), std::string::npos);
    clear_towered_section("Mv");
}

TEST(Conductors, MovingTowersContinueAcrossARestart)
{
    const std::string dir = scratch("moving_restart");
    set_towered_section(dir, "Mr");
    make_towers_move(0.02);   // still swaying at the checkpoint
    // coupling.dat every 2 steps, so the checkpoint falls between two rows
    amrex::ParmParse("conductors").add("diagnostics_int", 2);
    Mesh m(true);
    const double dt = 0.25;
    auto a = Conductors::create(0);
    ASSERT_TRUE(a);
    a->set_ground(m.znd.get(), m.geom);
    int step = 0;
    for (; step < 3; ++step) { a->advance(0, dt * step, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); }
    const std::string chk = dir + "/chk00003";
    std::filesystem::create_directories(chk);
    a->write_checkpoint(chk);
    // the iterations counted towards coupling.dat's next row: a moving tower's step takes at least one
    auto logged = [&] (const std::string& d) {
        const std::string state = slurp(d + "/conductors/state");
        const auto at = state.find("coupling_logged = ");
        return (at == std::string::npos) ? std::string() : state.substr(at, state.find('\n', at) - at);
    };
    const std::string counted = logged(chk);
    ASSERT_FALSE(counted.empty());
    EXPECT_NE(counted.compare(0, 20, "coupling_logged = 0 "), 0) << "the checkpoint must fall between two rows: " << counted;
    for (; step < 6; ++step) { a->advance(0, dt * step, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); }
    const std::string coupling = slurp(dir + "/coupling.dat");
    const auto* ma = dynamic_cast<const erf_towers::OneModeTower*>(a->tower_models()[1].get());
    ASSERT_NE(ma, nullptr);
    EXPECT_GT(std::abs(ma->v()[1]), 1.0e-4) << "the test needs the towers moving";
    const std::string towers = slurp(dir + "/towers.dat");
    const std::string stats = slurp(dir + "/tower_Mr_t2_stats.csv");
    const std::string span = slurp(dir + "/Mr_span2.dat");

    auto b = Conductors::create(0);
    ASSERT_TRUE(b);
    b->set_ground(m.znd.get(), m.geom, chk);
    ASSERT_TRUE(b->restored());
    {
        // the restart takes the count up where the checkpoint left it (coupling.dat's max hides a lost count when
        // the steps after the checkpoint take as many iterations)
        const std::string again = dir + "/chk00003_again";
        std::filesystem::create_directories(again);
        b->write_checkpoint(again);
        EXPECT_EQ(logged(again), counted);
    }
    for (step = 3; step < 6; ++step) { b->advance(0, dt * step, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom); }
    const auto* mb = dynamic_cast<const erf_towers::OneModeTower*>(b->tower_models()[1].get());
    ASSERT_NE(mb, nullptr);
    EXPECT_EQ(mb->q()[1], ma->q()[1]);
    EXPECT_EQ(mb->v()[1], ma->v()[1]);
    EXPECT_EQ(slurp(dir + "/towers.dat"), towers);
    EXPECT_EQ(slurp(dir + "/tower_Mr_t2_stats.csv"), stats);
    EXPECT_EQ(slurp(dir + "/Mr_span2.dat"), span);
    EXPECT_FALSE(coupling.empty());
    EXPECT_EQ(slurp(dir + "/coupling.dat"), coupling);
    clear_towered_section("Mr");
}

namespace {
// a circuit along x on the ramp: the middle phase C2 owns the lattice towers at x = 400 and 700,
// the outer phases C1 and C3 hang 5.5 m to either side of it on the cross-arm (1.5 m at the dead
// ends) and the shield wire CS 7 m above the cross-arm on the peak; the phases on strings
void set_circuit (const std::string& dir, Real damping = 0.0)
{
    amrex::ParmParse pp("conductors");
    pp.addarr("lines", std::vector<std::string>{"C1", "C2", "C3", "CS"});
    pp.add("diagnostics_dir", dir);
    pp.add("air_density", 1.2);
    for (const char* key : {"prescribed_velocity", "drag_on_flow", "epsilon", "node_output_int", "stats_start", "flashover_distance",
                            "ground_contact", "moordyn_dt"}) { pp.remove(key); }
    pp.addarr("prescribed_velocity", std::vector<Real>{0.0, 10.0, 0.0});
    pp.addarr("tower_types", std::vector<std::string>{"lat"});
    amrex::ParmParse pt("conductors.lat");
    for (const char* k : {"base_width", "top_width", "solidity", "arm_length", "arm_depth", "peak", "weight", "allowable_uplift",
                          "allowable_compression", "frequency", "damping_ratio"}) { pt.remove(k); }
    pt.add("base_width", 6.0); pt.add("top_width", 1.5); pt.add("solidity", 0.2); pt.add("arm_length", 12.0);
    pt.add("arm_depth", 1.2); pt.add("peak", 8.0); pt.add("weight", 9.0e4);
    if (damping > 0.0) { pt.add("frequency", 2.0); pt.add("damping_ratio", damping); }
    const struct { const char* name; Real tower_y, end_y, z_tower, z_end; bool phase; } lines[] = {
        {"C1", 494.5, 498.5, 30.0, 30.0, true}, {"C2", 500.0, 500.0, 30.0, 30.0, true},
        {"C3", 505.5, 501.5, 30.0, 30.0, true}, {"CS", 500.0, 500.0, 37.0, 32.0, false}};
    for (const auto& l : lines) {
        amrex::ParmParse ps(std::string("conductors.") + l.name);
        for (const char* k : {"end_a", "end_b", "towers", "length", "stringing_tension", "diameter", "mass_per_length", "axial_stiffness",
                              "tower_type", "share_towers", "insulator_length", "insulator_mass", "output_root"}) { ps.remove(k); }
        ps.addarr("end_a", std::vector<Real>{100.0, l.end_y, l.z_end});
        ps.addarr("end_b", std::vector<Real>{1000.0, l.end_y, l.z_end});
        ps.addarr("towers", std::vector<Real>{400.0, l.tower_y, l.z_tower, 700.0, l.tower_y, l.z_tower});
        ps.add("stringing_tension", l.phase ? 20000.0 : 10000.0);
        ps.add("diameter", l.phase ? 0.0281 : 0.0111);
        ps.add("mass_per_length", l.phase ? 1.628 : 0.406);
        ps.add("axial_stiffness", l.phase ? 3.0e7 : 9.7e6);
        if (std::string(l.name) == "C2") { ps.add("tower_type", std::string("lat")); }
        else { ps.add("share_towers", std::string("C2")); }
        if (l.phase) { ps.add("insulator_length", 2.5); ps.add("insulator_mass", 60.0); }
        ps.add("output_root", dir + "/" + l.name);
    }
}
void clear_circuit ()
{
    amrex::ParmParse pp("conductors");
    pp.addarr("lines", std::vector<std::string>{});
    pp.addarr("tower_types", std::vector<std::string>{});
    amrex::ParmParse pt("conductors.lat");
    for (const char* k : {"peak", "frequency", "damping_ratio"}) { pt.remove(k); }
    for (const char* n : {"C1", "C2", "C3", "CS"}) {
        amrex::ParmParse ps(std::string("conductors.") + n);
        for (const char* k : {"tower_type", "share_towers", "stringing_tension", "insulator_length", "insulator_mass"}) { ps.remove(k); }
    }
}
} // namespace

TEST(Conductors, ACircuitHangsFromOneRowOfTowersEachLineAtItsOwnPoint)
{
    const std::string dir = scratch("circuit_towers");
    set_circuit(dir);
    Mesh m(true);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(m.znd.get(), m.geom);
    ASSERT_EQ(c->towers().size(), 2u) << "one row of towers for the four lines";
    for (std::size_t t = 0; t < 2; ++t) {
        const auto& tw = c->towers()[t];
        EXPECT_EQ(tw.name(), "C2_t" + std::to_string(t + 1));
        ASSERT_EQ(c->tower_lines()[t].size(), 4u);
        ASSERT_EQ(tw.attachments().size(), 4u);
        // the owner first, at the centre of the cross-arm, then the others in the order of conductors.lines
        EXPECT_EQ(c->lines()[c->tower_lines()[t][0].first]->name(), "C2");
        EXPECT_EQ(c->lines()[c->tower_lines()[t][1].first]->name(), "C1");
        EXPECT_EQ(c->lines()[c->tower_lines()[t][2].first]->name(), "C3");
        EXPECT_EQ(c->lines()[c->tower_lines()[t][3].first]->name(), "CS");
        // every line's point stands on the tower's base: the ramp rises 0.11 m across the 5.5 m to C1,
        // which the cross-arm does not follow
        for (std::size_t a = 0; a < 4; ++a) {
            const auto [line, j] = c->tower_lines()[t][a];
            const auto& p = c->lines()[line]->inputs().point(j + 1);
            const Real above = (c->lines()[line]->name() == "CS") ? Real(37.0) : Real(30.0);
            EXPECT_NEAR(p[2], tw.base()[2] + above, roundoff * 1000) << c->lines()[line]->name();
            EXPECT_EQ(tw.attachments()[a], p);
        }
        const auto& c1 = c->lines()[0]->inputs().point(static_cast<int>(t) + 1);
        EXPECT_GT(std::abs(c1[2] - (m.h(c1[0], c1[1]) + 30.0)), 0.1) << "the test must see the level cross-arm";
    }
    std::vector<std::array<Real,3>> pull(2, {{0.0, 0.0, 0.0}});
    for (int s = 0; s < 3; ++s) {
        for (std::size_t t = 0; t < 2; ++t) {
            pull[t] = {{0.0, 0.0, 0.0}};
            for (const auto& [line, j] : c->tower_lines()[t]) {
                const auto f = c->lines()[line]->tower_force(j);
                for (int d = 0; d < 3; ++d) { pull[t][d] += f[d]; }
            }
        }
        c->advance(0, 0.25 * s, 0.25, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom);
    }
    for (std::size_t t = 0; t < 2; ++t) {
        const auto& tw = c->towers()[t];
        for (int d = 0; d < 3; ++d) { EXPECT_NEAR(tw.line_force()[d], pull[t][d], roundoff * 1.0e5) << "the four lines' pull, dir " << d; }
        EXPECT_LT(pull[t][2], -3.0 * 4000.0) << "three phases and a shield wire weigh on the tower";
        const auto L = tw.foundation();
        EXPECT_NEAR(L.vertical, 9.0e4 - pull[t][2] - tw.total_force()[2], roundoff * L.vertical);
    }
    clear_circuit();
}

TEST(Conductors, ACircuitOnBendingTowersMovesEveryLinesPoint)
{
    const std::string dir = scratch("circuit_moving");
    set_circuit(dir, 0.3);
    amrex::ParmParse("conductors").add("diagnostics_int", 4);
    Mesh m(true);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(m.znd.get(), m.geom);
    for (const auto& span : c->lines()) { EXPECT_TRUE(span->towers_move()) << span->name() << " moves with the towers it shares"; }
    const double dt = 0.25;
    int most = 0, unconverged = 0;
    std::vector<int> per_step;
    for (int s = 0; s < 40; ++s) {
        c->advance(0, dt * s, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom);
        most = std::max(most, c->coupling_iterations());
        per_step.push_back(c->coupling_iterations());
        unconverged += c->coupling_unconverged();
    }
    EXPECT_GE(most, 2) << "the coupling steps are iterated while the towers move";
    EXPECT_EQ(unconverged, 0);
    // coupling.dat, every 4 steps, holds the most iterations over the 4 steps since its last row
    {
        std::istringstream log(slurp(dir + "/coupling.dat"));
        std::string row;
        std::getline(log, row);
        std::size_t k = 0;
        while (std::getline(log, row)) {
            std::istringstream rs(row);
            double t = 0.0;
            int it = -1, bad = -1;
            ASSERT_TRUE(static_cast<bool>(rs >> t >> it >> bad)) << row;
            ASSERT_LT(4 * k + 3, per_step.size());
            EXPECT_EQ(it, *std::max_element(per_step.begin() + static_cast<long>(4 * k), per_step.begin() + static_cast<long>(4 * k + 4)))
                << "row " << k;
            ++k;
        }
        EXPECT_EQ(k, 10u);
    }
    for (std::size_t t = 0; t < 2; ++t) {
        const auto& tw = c->towers()[t];
        const auto* model = dynamic_cast<const erf_towers::OneModeTower*>(c->tower_models()[t].get());
        ASSERT_NE(model, nullptr);
        // at rest the stiffness carries the drag by the shape and every line's pull by the shape where it hangs
        double Q = 0.0;
        for (std::size_t a = 0; a < tw.line_forces().size(); ++a) { Q += model->attachment_shape(a) * tw.line_forces()[a][1]; }
        for (std::size_t i = 0; i < tw.nodes().size(); ++i) { Q += model->mode_shape(i) * tw.loads()[3*i+1]; }
        EXPECT_GT(model->q()[1], 1.0e-3) << tw.name() << " leans with the +y wind";
        EXPECT_NEAR(model->q()[1] * model->stiffness() / Q, 1.0, 0.02) << tw.name();
        // each line's point has moved with the tower by the shape at its height
        const Real ptol = (std::is_same<Real, float>::value) ? Real(1.0e-4) : Real(1.0e-6);
        for (std::size_t a = 0; a < tw.attachments().size(); ++a) {
            const auto [line, j] = c->tower_lines()[t][a];
            const erf_conductors::ConductorLine& span = *c->lines()[line];
            // the point the line hangs from: the top of its string, or the first node of the span after the tower
            const unsigned node = span.num_insulators() > 0
                ? span.span_first_node(span.num_spans()) + static_cast<unsigned>((erf_conductors::LineInputs::insulator_segments + 1) * j)
                : span.span_first_node(j + 1);
            const auto p = span.node_position(node);
            const auto x = model->attachment_displacement(a);
            EXPECT_NEAR(p[1], tw.attachments()[a][1] + x[1], ptol) << span.name() << " on " << tw.name();
        }
    }
    clear_circuit();
}

TEST(Conductors, AShortTautSpanOnBendingTowersStaysStable)
{
    // the shield wire's first span is 40 m long, clamped to the peak and strung tight: its pull
    // changes by ~2.4e5 N/m of the peak's travel, 1.5 times the tower's own stiffness once weighted
    // by the shape there; with the pull lagging a coupling step behind, the exchange pumps energy in
    const std::string dir = scratch("short_span");
    set_circuit(dir, 0.02);
    for (const char* n : {"C1", "C2", "C3", "CS"}) {
        amrex::ParmParse ps(std::string("conductors.") + n);
        std::vector<Real> t;
        ps.getarr("towers", t);
        t[0] = 140.0;
        ps.remove("towers");
        ps.addarr("towers", t);
    }
    Mesh m(true);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(m.znd.get(), m.geom);
    const double dt = 0.25;
    Real peak = 0.0, last = 0.0;
    int unconverged = 0;
    for (int s = 0; s < 60; ++s) {
        c->advance(0, dt * s, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom);
        unconverged += c->coupling_unconverged();
        const auto x = c->towers()[0].arm_displacement();
        last = std::hypot(x[0], x[1]);
        peak = std::max(peak, last);
    }
    EXPECT_EQ(unconverged, 0);
    EXPECT_LT(peak, Real(0.05)) << "the first tower's cross-arm stays within a few centimetres";
    EXPECT_LT(last, peak) << "and settles";
    // about 4 mm, where the stiffness holds the drag and the pulls; without the iteration of each
    // coupling step the exchange pumps energy into the cross-arm and this bound fails
    EXPECT_LT(last, Real(0.006));
    clear_circuit();
}

TEST(Conductors, AnImmersedTerrainPlacesEverythingOnItsSurface)
{
    // the same section and towers on the ramp, once on the fitted mesh and once on a flat mesh with
    // the ramp given as the immersed terrain's surface
    // the ramp given as the immersed terrain's surface; with ground_contact both write the same seafloors
    std::vector<std::array<Real,3>> fitted;
    std::vector<Real> fitted_clearance;
    std::vector<std::vector<double>> fitted_seafloor;
    for (int pass = 0; pass < 2; ++pass) {
        const std::string dir = scratch(pass == 0 ? "ib_fitted" : "ib_surface");
        set_circuit(dir);
        amrex::ParmParse("conductors").add("ground_contact", true);
        Mesh m(pass == 0);
        auto c = Conductors::create(0);
        ASSERT_TRUE(c);
        if (pass == 1) {
            const amrex::Box nodes = amrex::surroundingNodes(m.geom.Domain());
            amrex::FArrayBox h(amrex::makeSlab(amrex::grow(nodes, 3), 2, 0), 1);
            const Real dx = m.Lx / m.nx, dy = m.Ly / m.ny;
            const auto a = h.array();
            amrex::LoopOnCpu(h.box(), [&](int i, int j, int k) { a(i,j,k) = m.h(i * dx, j * dy); });
            c->set_ground_surface(h, m.geom);
        }
        c->set_ground(pass == 0 ? m.znd.get() : nullptr, m.geom);
        std::vector<std::array<Real,3>> placed;
        std::vector<Real> clear;
        for (const auto& span : c->lines()) {
            for (int k = 0; k <= span->num_spans(); ++k) { placed.push_back(span->inputs().point(k)); }
            for (unsigned i = 0; i < span->num_nodes(); ++i) { clear.push_back(span->clearance(i)); }
        }
        for (const auto& tw : c->towers()) { placed.push_back(tw.base()); }
        amrex::ParmParse("conductors").remove("ground_contact");
        std::vector<std::vector<double>> seafloor;
        for (const auto& span : c->lines()) { seafloor.push_back(SeafloorFile(dir + "/" + span->name() + ".seafloor.txt").z); }
        if (pass == 0) { fitted = placed; fitted_clearance = clear; fitted_seafloor = seafloor; continue; }
        ASSERT_EQ(seafloor.size(), fitted_seafloor.size());
        for (std::size_t l = 0; l < seafloor.size(); ++l) {
            ASSERT_EQ(seafloor[l].size(), fitted_seafloor[l].size()) << "line " << l;
            for (std::size_t q = 0; q < seafloor[l].size(); ++q) {
                EXPECT_NEAR(seafloor[l][q], fitted_seafloor[l][q], 1.0e-6) << "line " << l << " point " << q;
            }
        }
        ASSERT_EQ(placed.size(), fitted.size());
        for (std::size_t p = 0; p < placed.size(); ++p) {
            for (int d = 0; d < 3; ++d) { EXPECT_NEAR(placed[p][d], fitted[p][d], roundoff * 1000) << "point " << p << " dir " << d; }
        }
        ASSERT_EQ(clear.size(), fitted_clearance.size());
        for (std::size_t i = 0; i < clear.size(); ++i) { EXPECT_NEAR(clear[i], fitted_clearance[i], roundoff * 1000) << "node " << i; }
        // the flat mesh's own bottom is z = 0: without the surface the hills would not be seen
        EXPECT_GT(fitted[0][2], 30.0 + 1.0);
    }
    clear_circuit();
}

TEST(Conductors, SetGroundRefusesASurfaceOffsetBelowTheDomainTop)
{
    const std::string dir = scratch("frame");
    set_inputs(dir, "Tframe");
    amrex::ParmParse pp("conductors");
    pp.add("surface_offset", 300.0);   // below the 400 m domain top
    Mesh m(false);
    auto c = Conductors::create(0);
    pp.remove("surface_offset");
    ASSERT_TRUE(c);
    const std::string msg = erf_gtest::abort_message([&] { c->set_ground(nullptr, m.geom); });
    EXPECT_NE(msg.find("conductors.surface_offset"), std::string::npos) << msg;
    EXPECT_NE(msg.find("geometry.prob_hi[2]"), std::string::npos) << msg;
    EXPECT_FALSE(c->ground_set()) << "the check comes before anything is placed";
}

TEST(Conductors, AnAttachmentOutsideTheDomainIsRefusedNamingItsKey)
{
    const std::string dir = scratch("outside");
    // the domain ends at x = 1200 m
    set_inputs(dir, "Tout", true, {{300.0, 500.0, 30.0}}, {{1300.0, 500.0, 30.0}});
    Mesh m(false);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    const std::string msg = erf_gtest::abort_message([&] { c->set_ground(nullptr, m.geom); });
    EXPECT_NE(msg.find("conductors.Tout.end_b at (1300"), std::string::npos) << msg;
    EXPECT_NE(msg.find("lies outside the domain horizontally"), std::string::npos) << msg;
}

TEST(Conductors, ANonFiniteCouplingPullIsRefusedNotConverged)
{
    const std::vector<double> F0{1000.0, 0.0, -2000.0}, F{1000.0, 0.0, -2000.0};
    bool converged = false;
    // the tolerance is 1e-4 of the largest pull at the step's start, 2000 N: 0.2 N
    EXPECT_TRUE(erf_conductors::coupling_converged(F0, F, {1000.1, 0.0, -2000.0}, 1.0e-4, converged).empty());
    EXPECT_TRUE(converged) << "0.1 N within 0.2 N";
    EXPECT_TRUE(erf_conductors::coupling_converged(F0, F, {1001.0, 0.0, -2000.0}, 1.0e-4, converged).empty());
    EXPECT_FALSE(converged) << "1 N beyond 0.2 N";
    const std::string err = erf_conductors::coupling_converged(F0, F, {std::numeric_limits<double>::quiet_NaN(), 0.0, -2000.0},
                                                               1.0e-4, converged);
    EXPECT_NE(err.find("is not finite"), std::string::npos) << err;
    EXPECT_FALSE(converged) << "a pull that is not finite never counts as converged";
    // unconverged coupling steps in a row stop the run at the limit; a converged one starts the count again
    int in_a_row = 0;
    for (int n = 1; n < 20; ++n) { EXPECT_FALSE(erf_conductors::coupling_diverged(false, in_a_row, 20)) << n; }
    EXPECT_FALSE(erf_conductors::coupling_diverged(true, in_a_row, 20));
    EXPECT_EQ(in_a_row, 0);
    for (int n = 1; n < 20; ++n) { EXPECT_FALSE(erf_conductors::coupling_diverged(false, in_a_row, 20)) << n; }
    EXPECT_TRUE(erf_conductors::coupling_diverged(false, in_a_row, 20)) << "the 20th in a row";
}

TEST(Conductors, ARestartChecksTheTowerSwayAndTheSurfaceOffset)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_TRUE(erf_conductors::restart_mismatch(false, false, 10000.0, 10000.0).empty());
    EXPECT_TRUE(erf_conductors::restart_mismatch(true, true, 10000.0, 10000.0).empty());
    EXPECT_TRUE(erf_conductors::restart_mismatch(false, false, nan, 5000.0).empty()) << "a checkpoint without the record";
    EXPECT_NE(erf_conductors::restart_mismatch(false, true, nan, 10000.0).find("holds no tower sway"), std::string::npos);
    EXPECT_NE(erf_conductors::restart_mismatch(true, false, nan, 10000.0).find("holds tower sway, but no tower type"), std::string::npos);
    EXPECT_NE(erf_conductors::restart_mismatch(false, false, 10000.0, 5000.0).find("conductors.surface_offset"), std::string::npos);
    // the tower forms: the same continue; a change, or an older checkpoint (no record: every form 0) with a form on, is
    // refused naming the keys
    using erf_conductors::tower_forms_mismatch;
    EXPECT_TRUE(tower_forms_mismatch("lattice", true, {{1, 1, 0}}, {{1, 1, 0}}).empty());
    EXPECT_TRUE(tower_forms_mismatch("lattice", false, {{0, 0, 0}}, {{0, 0, 0}}).empty()) << "an older checkpoint, old forms";
    const std::string changed = tower_forms_mismatch("lattice", true, {{1, 1, 1}}, {{1, 0, 1}});
    EXPECT_NE(changed.find("conductors.lattice.arm_outside_shaft = false"), std::string::npos) << changed;
    EXPECT_EQ(changed.find("diagonal_wind_factor"), std::string::npos) << changed;
    const std::string older = tower_forms_mismatch("lattice", false, {{0, 0, 0}}, {{1, 1, 1}});
    for (const char* k : {"diagonal_wind_factor", "arm_outside_shaft", "angle_principal_axes"}) {
        EXPECT_NE(older.find(std::string("conductors.lattice.") + k), std::string::npos) << older;
    }
    EXPECT_NE(older.find("to false"), std::string::npos) << older;
    EXPECT_EQ(older.find("= true"), std::string::npos) << "the keys alone: " << older;
}

namespace {
// rho and rho k = rho (k0 + kz z) in every cell of a flat mesh: k linear in the height, which the cell sampler reproduces
void fill_k (amrex::MultiFab& cons, const Mesh& m, Real k0, Real kz)
{
    cons.setVal(0.0);
    const Real dz = m.H / m.nz;
    for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
        auto ca = cons.array(mfi);
        amrex::LoopOnCpu(mfi.growntilebox(), [&](int i, int j, int k) {
            amrex::ignore_unused(i, j);
            ca(i,j,k,Rho_comp) = 1.2;
            ca(i,j,k,RhoKE_comp) = Real(1.2) * (k0 + kz * (k + Real(0.5)) * dz);
        });
    }
}

// the mean of k = k0 + kz z over the nodes [n0, n0 + nn) of the points pos (3 per point)
double mean_k (const std::vector<Real>& pos, unsigned n0, unsigned nn, Real k0, Real kz)
{
    double k = 0.0;
    for (unsigned n = n0; n < n0 + nn; ++n) { k += static_cast<double>(k0 + kz * pos[3*n+2]); }
    return k / nn;
}

// the rows of a whitespace-separated log after its header, as numbers
std::vector<std::vector<double>> read_rows (const std::string& fname, std::string& header)
{
    std::ifstream f(fname);
    std::getline(f, header);
    std::vector<std::vector<double>> rows;
    for (std::string line; std::getline(f, line);) {
        std::istringstream ls(line);
        rows.emplace_back();
        for (double x; ls >> x;) { rows.back().push_back(x); }
    }
    return rows;
}
} // namespace

TEST(Conductors, AnEventGustPassesEveryNodeAndTowerWithItsFront)
{
    const std::string dir = scratch("gust_event");
    set_towered_section(dir, "Tev");   // three spans along x from 100 to 1000 m, towers at 400 and 700 m
    amrex::ParmParse pp("conductors");
    pp.remove("prescribed_velocity");
    pp.add("gust_type", std::string("event"));
    pp.add("gust_event_time", 1.0);
    pp.add("gust_event_speed", 150.0);
    pp.add("gust_event_duration", 2.0);
    Mesh m(false);
    m.fill_crosswind();
    const Real k0 = 0.5, kz = 0.01;
    amrex::MultiFab cons(m.ba, m.dm, RhoKE_comp + 1, 1);
    fill_k(cons, m, k0, kz);
    {
        // without an origin the front crosses the domain's centre, moving along +x
        auto c0 = Conductors::create(0);
        ASSERT_TRUE(c0);
        c0->set_closure(true, 0.5562);
        c0->set_ground(nullptr, m.geom);
        EXPECT_NEAR(c0->gust_event().x0, 600.0, 1e-9);
        EXPECT_NEAR(c0->gust_event().y0, 500.0, 1e-9);
        EXPECT_NEAR(c0->gust_event().ex, 1.0, 1e-15);
        EXPECT_NEAR(c0->gust_event().ey, 0.0, 1e-15);
    }
    {
        // a direction in degrees counterclockwise from +x: 90 moves the front along +y
        pp.add("gust_event_direction", 90.0);
        auto c90 = Conductors::create(0);
        ASSERT_TRUE(c90);
        c90->set_closure(true, 0.5562);
        c90->set_ground(nullptr, m.geom);
        EXPECT_NEAR(c90->gust_event().ex, 0.0, 1e-6);
        EXPECT_NEAR(c90->gust_event().ey, 1.0, 1e-12);
        pp.add("gust_event_direction", 0.0);
    }
    // the front crosses the line's first end at t = 1 s and its last at 7 s; it lasts 2 s at a point
    pp.addarr("gust_event_origin", std::vector<Real>{100.0, 500.0});
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_closure(true, 0.5562);
    c->set_ground(nullptr, m.geom);
    EXPECT_NE(erf_gtest::abort_message([&] { c->advance(0, 0.0, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom); })
                  .find("gust_type = event needs the conserved state"), std::string::npos);
    const auto& e = c->gust_event();
    EXPECT_NEAR(e.x0, 100.0, 1e-9);
    EXPECT_EQ(e.time, 1.0);
    EXPECT_EQ(e.speed, 150.0);
    EXPECT_EQ(e.duration, 2.0);
    const double sf = c->gust_sigma_factor(), g = 2.7, dt = 0.25;
    const auto& line = *c->lines().front();
    ASSERT_EQ(c->towers().size(), 2u);
    int during = 0, before = 0;
    std::vector<double> mid_gust;   // the gust at span 1's middle node, step by step
    const int steps = 40;
    for (int s = 0; s < steps; ++s) {
        const double tm = dt * s + 0.5 * dt;   // the wind is held at the step's middle
        const std::vector<Real> where = line.kinematics_points();
        c->advance(0, dt * s, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom, &cons);
        for (int sp = 0; sp < line.num_spans(); ++sp) {
            const unsigned n0 = line.span_first_node(sp), nn = line.span_num_nodes(sp);
            const double A = g * sf * std::sqrt(mean_k(where, n0, nn, k0, kz));
            for (unsigned n = n0; n < n0 + nn; ++n) {
                const double phase = erf_conductors::gust_event_phase(e, where[3*n], where[3*n+1], tm);
                const double gust = A * erf_conductors::gust_event_shape(phase);
                const double expect = static_cast<double>(m.v_at(where[3*n+1], where[3*n+2])) + gust;
                const auto uvw = line.wind_at_point(n);
                ASSERT_NEAR(uvw[1], expect, roundoff * std::max(1.0, expect)) << "step " << s << " node " << n;
                ASSERT_NEAR(uvw[0], 0.0, roundoff) << "along the +y wind only";
                if (gust > 0.5) { ++during; } else if (phase < 0.0) { ++before; }
                if (sp == 0 && n == n0 + (nn - 1) / 2) { mid_gust.push_back(gust); }
            }
        }
        // the towers stand still: their drag nodes where they were built; the tower's k is the mean over its nodes
        // weighted by the face area each stands for times the wind on it, which differs from the plain mean since k
        // and the wind grow with the height
        for (std::size_t t = 0; t < 2; ++t) {
            const auto& nodes = c->towers()[t].nodes();
            double kk = 0.0, area = 0.0, plain = 0.0;
            for (const auto& nd : nodes) {
                const double w = static_cast<double>(nd.drag_width * nd.length) * std::abs(static_cast<double>(m.v_at(nd.pos[1], nd.pos[2])));
                kk += w * static_cast<double>(k0 + kz * nd.pos[2]);
                plain += static_cast<double>(k0 + kz * nd.pos[2]) / nodes.size();
                area += w;
            }
            ASSERT_GT(std::abs(kk / area - plain), 1.0e-3 * plain) << "the weighting must tell: not a vacuous check";
            const double A = g * sf * std::sqrt(kk / area);
            const auto& w = c->tower_winds()[t];
            ASSERT_EQ(w.size(), 3 * nodes.size());
            for (std::size_t j = 0; j < nodes.size(); ++j) {
                const double phase = erf_conductors::gust_event_phase(e, nodes[j].pos[0], nodes[j].pos[1], tm);
                const double gust = A * erf_conductors::gust_event_shape(phase);
                const double expect = static_cast<double>(m.v_at(nodes[j].pos[1], nodes[j].pos[2])) + gust;
                ASSERT_NEAR(w[3*j+1], expect, roundoff * std::max(1.0, expect)) << "tower " << t << " node " << j << " step " << s;
            }
        }
    }
    EXPECT_GT(during, 100) << "the gust must blow on the line for a while";
    EXPECT_GT(before, 100) << "and not before its front";
    // gust_series.dat: a row per step, the gust at each span's middle node and at each tower's highest drag node
    std::string header;
    const auto rows = read_rows(dir + "/gust_series.dat", header);
    EXPECT_EQ(header, "time Tev_span1 Tev_span2 Tev_span3 Tev_t1 Tev_t2");
    ASSERT_EQ(rows.size(), static_cast<std::size_t>(steps));
    for (int s = 0; s < steps; ++s) {
        const auto ss = static_cast<std::size_t>(s);
        ASSERT_EQ(rows[ss].size(), 6u);
        EXPECT_NEAR(rows[ss][0], dt * s, 1e-12) << "stamped with the step's start";
        EXPECT_NEAR(rows[ss][1], mid_gust[ss], 1e-6 * std::max(1.0, mid_gust[ss]) + roundoff) << "step " << s;
    }
    clear_towered_section("Tev");
}

TEST(Conductors, RandomGustsFollowTheirProcessesAndACircuitSharesThem)
{
    const std::string dir = scratch("gust_random");
    set_circuit(dir);   // C2 owns the towers; C1, C3 and the shield wire CS share them
    amrex::ParmParse pp("conductors");
    pp.remove("prescribed_velocity");
    pp.add("gust_type", std::string("random"));
    pp.add("gust_seed", 11);
    Mesh m(false);
    m.fill_crosswind();
    const Real k0 = 0.5, kz = 0.01;
    amrex::MultiFab cons(m.ba, m.dm, RhoKE_comp + 1, 1);
    fill_k(cons, m, k0, kz);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_closure(true, 0.5562);
    c->set_ground(nullptr, m.geom);
    ASSERT_EQ(c->lines().size(), 4u);
    ASSERT_EQ(c->towers().size(), 2u);
    EXPECT_NE(erf_gtest::abort_message([&] { c->advance(0, 0.0, 0.5, m.u, m.v, m.w, nullptr, nullptr, m.geom); })
                  .find("gust_type = random needs the conserved state"), std::string::npos);
    const double sf = c->gust_sigma_factor(), Ls = 67.056, dt = 0.5;
    // the processes: C2's three spans (0 to 2), whose time scales come from C2's nodes, then the two towers (3, 4)
    std::vector<double> z(5, 0.0);
    const auto& C2 = *c->lines()[1];
    double zmin = 0.0, zmax = 0.0;
    const int steps = 12;
    for (int s = 0; s < steps; ++s) {
        const auto step = static_cast<std::uint64_t>(s + 1);   // the step count the draws take, from 1
        std::vector<std::vector<Real>> where(4);
        for (std::size_t i = 0; i < 4; ++i) { where[i] = c->lines()[i]->kinematics_points(); }
        auto draw = [&] (std::size_t p, double height, double U) {
            const double T = erf_conductors::gust_integral_length(height) / U;
            const double xi = erf_conductors::gust_normal(11, p, step);
            z[p] = (s == 0) ? xi : erf_conductors::gust_ou_step(z[p], dt, T, xi);
        };
        for (int sp = 0; sp < 3; ++sp) {
            const unsigned n0 = C2.span_first_node(sp), nn = C2.span_num_nodes(sp);
            double u2 = 0.0;
            for (unsigned n = n0; n < n0 + nn; ++n) {
                const double v = static_cast<double>(m.v_at(where[1][3*n+1], where[1][3*n+2]));
                u2 += v * v;
            }
            // the span's height: the mean of the heights it hangs from (the strings' lower ends at the towers), the ground flat at 0
            const double height = 0.5 * static_cast<double>(C2.inputs().conductor_point(sp)[2] + C2.inputs().conductor_point(sp + 1)[2]);
            draw(static_cast<std::size_t>(sp), height, std::sqrt(u2 / nn));
        }
        for (std::size_t t = 0; t < 2; ++t) {
            // the tower's wind and height: the means over its drag nodes weighted by the face area each stands for
            // times the wind on it (its share of the fluctuating load)
            const auto& tw = c->towers()[t];
            double u2 = 0.0, zm = 0.0, area = 0.0;
            for (const auto& nd : tw.nodes()) {
                const double v = static_cast<double>(m.v_at(nd.pos[1], nd.pos[2]));
                const double w = static_cast<double>(nd.drag_width * nd.length) * std::abs(v);
                u2 += w * v * v;
                zm += w * static_cast<double>(nd.pos[2] - tw.base()[2]);
                area += w;
            }
            draw(3 + t, zm / area, std::sqrt(u2 / area));
        }
        c->advance(0, dt * s, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom, &cons);
        for (std::size_t i = 0; i < 4; ++i) {
            const auto& L = *c->lines()[i];
            std::vector<double> gust(3);
            for (int sp = 0; sp < 3; ++sp) {
                const unsigned n0 = L.span_first_node(sp), nn = L.span_num_nodes(sp);
                const double sigma = sf * std::sqrt(mean_k(where[i], n0, nn, k0, kz));
                const double B = erf_conductors::gust_background_factor(static_cast<double>(L.inputs().conductor_chord(sp)), Ls);
                // every line's span takes C2's process of that span: the circuit's phases move together
                gust[static_cast<std::size_t>(sp)] = sigma * std::sqrt(B) * z[static_cast<std::size_t>(sp)];
                for (unsigned n = n0; n < n0 + nn; ++n) {
                    const double expect = static_cast<double>(m.v_at(where[i][3*n+1], where[i][3*n+2])) +
                                          gust[static_cast<std::size_t>(sp)];
                    ASSERT_NEAR(L.wind_at_point(n)[1], expect, roundoff * std::max(1.0, expect))
                        << L.name() << " step " << s << " node " << n;
                }
            }
            // a string takes the mean of the spans beside it
            for (int j = 0; j < L.num_insulators(); ++j) {
                const double A = 0.5 * (gust[static_cast<std::size_t>(j)] + gust[static_cast<std::size_t>(j) + 1]);
                for (unsigned n = L.string_first_node(j); n < L.string_first_node(j) + L.string_num_nodes(j); ++n) {
                    const double expect = static_cast<double>(m.v_at(where[i][3*n+1], where[i][3*n+2])) + A;
                    ASSERT_NEAR(L.wind_at_point(n)[1], expect, roundoff * std::max(1.0, expect))
                        << L.name() << " string " << j << " node " << n;
                }
            }
        }
        for (std::size_t t = 0; t < 2; ++t) {
            const auto& tw = c->towers()[t];
            // ASCE 74's h of the background factor is the tower's height, the cross-arm's plus the peak's; the
            // tower's k the load-weighted (area x wind) mean over its drag nodes
            double kk = 0.0, area = 0.0;
            for (const auto& nd : tw.nodes()) {
                const double w = static_cast<double>(nd.drag_width * nd.length) * std::abs(static_cast<double>(m.v_at(nd.pos[1], nd.pos[2])));
                kk += w * static_cast<double>(k0 + kz * nd.pos[2]);
                area += w;
            }
            const double top = static_cast<double>(tw.arm_height() + tw.type().peak);
            const double gust = sf * std::sqrt(kk / area) *
                                std::sqrt(erf_conductors::tower_background_factor(top, Ls)) * z[3 + t];
            for (std::size_t j = 0; j < tw.nodes().size(); ++j) {
                const double expect = static_cast<double>(m.v_at(tw.nodes()[j].pos[1], tw.nodes()[j].pos[2])) + gust;
                ASSERT_NEAR(c->tower_winds()[t][3*j+1], expect, roundoff * std::max(1.0, expect)) << "tower " << t << " node " << j;
            }
        }
        for (const double v : z) { zmin = std::min(zmin, v); zmax = std::max(zmax, v); }
    }
    EXPECT_LT(zmin, -0.3) << "the processes must swing both ways for the test to see the sign";
    EXPECT_GT(zmax, 0.3);
    std::string header;
    const auto rows = read_rows(dir + "/gust_series.dat", header);
    EXPECT_EQ(header, "time C1_span1 C1_span2 C1_span3 C2_span1 C2_span2 C2_span3 C3_span1 C3_span2 C3_span3 CS_span1 CS_span2 "
                      "CS_span3 C2_t1 C2_t2");
    ASSERT_EQ(rows.size(), static_cast<std::size_t>(steps));
    for (const auto& r : rows) { EXPECT_EQ(r.size(), 15u); }
    clear_circuit();
}

TEST(Conductors, RandomGustsContinueAcrossARestart)
{
    Mesh m(false);
    m.fill_crosswind();
    const Real k0 = 0.5, kz = 0.01;
    amrex::MultiFab cons(m.ba, m.dm, RhoKE_comp + 1, 1);
    fill_k(cons, m, k0, kz);
    const double dt = 0.25;
    auto setup = [] (const std::string& dir, const char* type) {
        set_inputs(dir, "Trg", false);
        amrex::ParmParse pp("conductors");
        pp.add("gust_type", std::string(type));
        pp.remove("gust_seed");
        if (std::string(type) == "random") { pp.add("gust_seed", 3); }
    };
    auto winds = [] (const Conductors& c) {
        const auto& line = *c.lines().front();
        std::vector<Real> w;
        for (unsigned n = 0; n < line.num_nodes(); ++n) { const auto u = line.wind_at_point(n); w.insert(w.end(), {u[0], u[1], u[2]}); }
        return w;
    };
    // six steps straight through
    const std::string dir_a = scratch("gust_restart_a");
    setup(dir_a, "random");
    auto a = Conductors::create(0);
    ASSERT_TRUE(a);
    a->set_closure(true, 0.5562);
    a->set_ground(nullptr, m.geom);
    for (int s = 0; s < 6; ++s) { a->advance(0, dt * s, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom, &cons); }
    const std::vector<Real> wind_a = winds(*a);
    // three steps, a checkpoint, and a restart for the other three
    const std::string dir_b = scratch("gust_restart_b");
    setup(dir_b, "random");
    auto b = Conductors::create(0);
    ASSERT_TRUE(b);
    b->set_closure(true, 0.5562);
    b->set_ground(nullptr, m.geom);
    for (int s = 0; s < 3; ++s) { b->advance(0, dt * s, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom, &cons); }
    const std::string chk = dir_b + "/chk00003";
    std::filesystem::create_directories(chk);
    b->write_checkpoint(chk);
    EXPECT_NE(slurp(chk + "/conductors/gust_processes").find("processes 1 set 1"), std::string::npos);
    auto b2 = Conductors::create(0);
    ASSERT_TRUE(b2);
    b2->set_closure(true, 0.5562);
    b2->set_ground(nullptr, m.geom, chk);
    ASSERT_TRUE(b2->restored());
    for (int s = 3; s < 6; ++s) { b2->advance(0, dt * s, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom, &cons); }
    const std::vector<Real> wind_b = winds(*b2);
    ASSERT_EQ(wind_a.size(), wind_b.size());
    for (std::size_t i = 0; i < wind_a.size(); ++i) {
        EXPECT_NEAR(wind_b[i], wind_a[i], roundoff * std::max(Real(1.0), std::abs(wind_a[i]))) << "component " << i;
    }
    std::string ha, hb;
    const auto rows_a = read_rows(dir_a + "/gust_series.dat", ha);
    const auto rows_b = read_rows(dir_b + "/gust_series.dat", hb);
    EXPECT_EQ(ha, hb);
    ASSERT_EQ(rows_b.size(), rows_a.size()) << "the restarted run continues the series from the checkpoint";
    for (std::size_t r = 0; r < rows_a.size(); ++r) {
        ASSERT_EQ(rows_b[r].size(), 2u);
        EXPECT_NEAR(rows_b[r][1], rows_a[r][1], 1e-6) << "row " << r;
    }
    EXPECT_GT(std::abs(rows_a.back()[1]), 1.0e-3) << "a gust to compare";
    // a checkpoint written without gusts: the restart draws the process afresh, from its stationary law, at its first step
    const std::string dir_c = scratch("gust_restart_c");
    setup(dir_c, "none");
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_ground(nullptr, m.geom);
    for (int s = 0; s < 3; ++s) { c->advance(0, dt * s, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom); }
    const std::string chk_c = dir_c + "/chk00003";
    std::filesystem::create_directories(chk_c);
    c->write_checkpoint(chk_c);
    EXPECT_FALSE(std::filesystem::exists(chk_c + "/conductors/gust_processes"));
    {
        amrex::ParmParse pp("conductors");
        pp.add("gust_type", std::string("random"));
        pp.add("gust_seed", 3);
    }
    auto c2 = Conductors::create(0);
    ASSERT_TRUE(c2);
    c2->set_closure(true, 0.5562);
    c2->set_ground(nullptr, m.geom, chk_c);
    ASSERT_TRUE(c2->restored());
    const auto& line = *c2->lines().front();
    const std::vector<Real> where = line.kinematics_points();
    c2->advance(0, dt * 3, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom, &cons);
    const unsigned nn = line.span_num_nodes(0);
    const double gust = c2->gust_sigma_factor() * std::sqrt(mean_k(where, 0, nn, k0, kz)) *
                        std::sqrt(erf_conductors::gust_background_factor(static_cast<double>(line.inputs().chord(0)), 67.056)) *
                        erf_conductors::gust_normal(3, 0, 4);
    for (unsigned n = 0; n < nn; ++n) {
        const double expect = static_cast<double>(m.v_at(where[3*n+1], where[3*n+2])) + gust;
        ASSERT_NEAR(line.wind_at_point(n)[1], expect, roundoff * std::max(1.0, expect)) << "node " << n;
    }
}

TEST(Conductors, TheGustFactorRunsOnLinesWithLatticeTowers)
{
    // gust_type = factor samples k at the line nodes only; the towers' drag nodes are not sampled
    const std::string dir = scratch("gust_factor_towers");
    set_towered_section(dir, "Tft");
    amrex::ParmParse pp("conductors");
    pp.remove("prescribed_velocity");
    pp.add("gust_type", std::string("factor"));
    Mesh m(false);
    m.fill_crosswind();
    amrex::MultiFab cons(m.ba, m.dm, RhoKE_comp + 1, 1);
    fill_k(cons, m, 0.5, 0.01);
    auto c = Conductors::create(0);
    ASSERT_TRUE(c);
    c->set_closure(true, 0.5562);
    c->set_ground(nullptr, m.geom);
    ASSERT_EQ(c->towers().size(), 2u);
    for (int s = 0; s < 3; ++s) { c->advance(0, 0.25 * s, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom, &cons); }
    std::ifstream f(dir + "/gusts.csv");
    std::string row;
    int rows = 0;
    while (std::getline(f, row)) { ++rows; }
    EXPECT_EQ(rows, 4) << "the header and one row per span";
    EXPECT_FALSE(std::filesystem::exists(dir + "/gust_series.dat")) << "the factor adds no gust to the wind";
    clear_towered_section("Tft");
}

TEST(Conductors, ALineTakesAnotherLinesGustsAndAChangedSharingDrawsAfresh)
{
    // three parallel lines 6 m apart without lattice towers: Tgb takes Tga's random gusts (gust_with), Tgc its own
    const std::string dir = scratch("gust_with");
    set_inputs(dir, "Tga", false);
    amrex::ParmParse pp("conductors");
    pp.addarr("lines", std::vector<std::string>{"Tga", "Tgb", "Tgc"});
    for (const char* name : {"Tgb", "Tgc"}) {
        amrex::ParmParse ps(std::string("conductors.") + name);
        const Real y = (std::string(name) == "Tgb") ? Real(506.0) : Real(512.0);
        ps.addarr("end_a", std::vector<Real>{300.0, y, 30.0});
        ps.addarr("end_b", std::vector<Real>{600.0, y, 30.0});
        ps.add("length", 301.5);
        ps.add("diameter", 0.0281);
        ps.add("mass_per_length", 1.628);
        ps.add("axial_stiffness", 3.0e7);
    }
    amrex::ParmParse("conductors.Tgb").add("gust_with", std::string("Tga"));
    pp.add("gust_type", std::string("random"));
    pp.add("gust_seed", 5);
    Mesh m(false);
    m.fill_crosswind();
    // a wind along x that grows with the height as well: the wind veers with the height, so a node's own direction
    // differs from its span's mean direction
    auto u_at = [] (Real z) { return Real(2.0) + Real(0.1) * z; };
    {
        const Real dz = m.H / m.nz;
        for (amrex::MFIter mfi(m.u); mfi.isValid(); ++mfi) {
            auto ua = m.u.array(mfi);
            amrex::LoopOnCpu(mfi.growntilebox(), [&](int i, int j, int k) {
                amrex::ignore_unused(i, j);
                ua(i,j,k) = u_at((k + Real(0.5)) * dz);
            });
        }
    }
    const Real k0 = 0.5, kz = 0.01;
    amrex::MultiFab cons(m.ba, m.dm, RhoKE_comp + 1, 1);
    fill_k(cons, m, k0, kz);
    const double dt = 0.25, Ls = 67.056;
    // the expected wind at a line's nodes, the flow's plus gust along the span's mean horizontal wind
    auto check = [&] (const Conductors& c, std::size_t i, const std::vector<Real>& where, double z, const std::string& what) {
        const auto& L = *c.lines()[i];
        const unsigned nn = L.span_num_nodes(0);
        double mu = 0.0, mv = 0.0;
        for (unsigned n = 0; n < nn; ++n) {
            mu += static_cast<double>(u_at(where[3*n+2]));
            mv += static_cast<double>(m.v_at(where[3*n+1], where[3*n+2]));
        }
        const double h = std::hypot(mu, mv);
        const double gust = c.gust_sigma_factor() * std::sqrt(mean_k(where, 0, nn, k0, kz)) *
                            std::sqrt(erf_conductors::gust_background_factor(static_cast<double>(L.inputs().chord(0)), Ls)) * z;
        for (unsigned n = 0; n < nn; ++n) {
            const double eu = static_cast<double>(u_at(where[3*n+2])) + gust * mu / h;
            const double ev = static_cast<double>(m.v_at(where[3*n+1], where[3*n+2])) + gust * mv / h;
            ASSERT_NEAR(L.wind_at_point(n)[0], eu, roundoff * std::max(1.0, std::abs(eu))) << what << " node " << n;
            ASSERT_NEAR(L.wind_at_point(n)[1], ev, roundoff * std::max(1.0, std::abs(ev))) << what << " node " << n;
        }
    };
    auto a = Conductors::create(0);
    ASSERT_TRUE(a);
    a->set_closure(true, 0.5562);
    a->set_ground(nullptr, m.geom);
    double z = 0.0;
    for (int s = 0; s < 3; ++s) {
        std::vector<std::vector<Real>> where{a->lines()[0]->kinematics_points(), a->lines()[1]->kinematics_points()};
        // one process, Tga's span, its time scale from Tga's nodes
        const unsigned nn = a->lines()[0]->span_num_nodes(0);
        double u2 = 0.0;
        for (unsigned n = 0; n < nn; ++n) {
            const double u = static_cast<double>(u_at(where[0][3*n+2])), v = static_cast<double>(m.v_at(where[0][3*n+1], where[0][3*n+2]));
            u2 += u * u + v * v;
        }
        const double T = erf_conductors::gust_integral_length(30.0) / std::sqrt(u2 / nn);
        const double xi = erf_conductors::gust_normal(5, 0, static_cast<std::uint64_t>(s + 1));
        z = (s == 0) ? xi : erf_conductors::gust_ou_step(z, dt, T, xi);
        a->advance(0, dt * s, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom, &cons);
        check(*a, 0, where[0], z, "Tga step " + std::to_string(s));
        check(*a, 1, where[1], z, "Tgb step " + std::to_string(s));
    }
    const std::string chk = dir + "/chk00003";
    std::filesystem::create_directories(chk);
    a->write_checkpoint(chk);
    EXPECT_NE(slurp(chk + "/conductors/gust_processes").find("Tgc_span1 "), std::string::npos);
    // the restart has Tgc take Tga's gusts instead of Tgb: two processes again, but not the same ones (Tga's and
    // Tgb's where the checkpoint holds Tga's and Tgc's), so both are drawn afresh, from the stationary law, at the
    // restart's first step
    amrex::ParmParse("conductors.Tgb").remove("gust_with");
    amrex::ParmParse("conductors.Tgc").add("gust_with", std::string("Tga"));
    auto b = Conductors::create(0);
    ASSERT_TRUE(b);
    b->set_closure(true, 0.5562);
    b->set_ground(nullptr, m.geom, chk);
    ASSERT_TRUE(b->restored());
    std::vector<std::vector<Real>> where{b->lines()[0]->kinematics_points(), b->lines()[1]->kinematics_points(),
                                         b->lines()[2]->kinematics_points()};
    b->advance(0, dt * 3, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom, &cons);
    check(*b, 0, where[0], erf_conductors::gust_normal(5, 0, 4), "Tga after the restart");
    check(*b, 1, where[1], erf_conductors::gust_normal(5, 1, 4), "Tgb after the restart");
    check(*b, 2, where[2], erf_conductors::gust_normal(5, 0, 4), "Tgc after the restart, with Tga's gusts");
}

namespace {
// the still-air end tension of the 300 m span of set_inputs (301.5 m unstretched) by the elastic catenary (N)
Real still_air_end_tension ()
{
    const Real w = (Real(1.628) - Real(1.2) * Real(0.25) * Real(3.14159265358979323846) * Real(0.0281) * Real(0.0281)) * Real(CONST_GRAV);
    const auto cat = erf_conductors::elastic_catenary(300.0, 301.5, w, 3.0e7);
    EXPECT_TRUE(cat.solved);
    return cat.end_tension;
}
// still air: the line keeps its catenary, whose tension the elastic catenary gives
void set_still_air (const std::string& dir, const std::string& name)
{
    set_inputs(dir, name, false);
    amrex::ParmParse("conductors").addarr("prescribed_velocity", std::vector<Real>{0.0, 0.0, 0.0});
}
void clear_capacities (const std::string& name)
{
    amrex::ParmParse ps("conductors." + name);
    for (const char* k : {"rated_breaking_strength", "insulator_rated_load", "hardware_rated_load"}) { ps.remove(k); }
    amrex::ParmParse pp("conductors");
    pp.remove("failure_stop");
    pp.remove("prescribed_velocity");
}
std::vector<std::string> data_rows (const std::string& file)
{
    std::ifstream f(file);
    std::vector<std::string> rows;
    for (std::string l; std::getline(f, l);) { if (!l.empty() && l.rfind("time", 0) != 0) { rows.push_back(l); } }
    return rows;
}
} // namespace

TEST(Conductors, CapacitiesAreCheckedAfterEveryStepOfTheAnchorLevelOnly)
{
    // Ann's multi-level point (#5009) for the capacities: the lines step on their anchor level, so the
    // check runs once per anchor-level step and is stamped with that step's end, whatever the other
    // levels do. Two levels, the anchor level 1: a level-0 call records nothing.
    const std::string dir = scratch("capacities_anchor");
    set_still_air(dir, "Tcap");
    const Real T = still_air_end_tension();
    amrex::ParmParse ps("conductors.Tcap");
    ps.add("rated_breaking_strength", 0.99 * T);   // just below the still-air end tension: exceeded at once
    ps.add("hardware_rated_load", 0.95 * T);       // the dead ends carry the end tension and half a segment's weight
    Mesh m(false);
    auto c = Conductors::create(1);
    ASSERT_TRUE(c);
    c->set_ground(nullptr, m.geom);
    const auto& cap = c->capacities();
    ASSERT_EQ(cap.size(), 3u) << "the span's conductor and the hardware at its two dead ends";
    EXPECT_EQ(cap.component(0).name, "Tcap");
    EXPECT_EQ(cap.component(1).name, "Tcap_end_a");
    EXPECT_EQ(cap.component(2).name, "Tcap_end_b");
    // the still-air span against the elastic catenary's end tension: the pin of the demand. 0.5 % holds both
    // libraries: the stub's largest tension includes the end node (-0.03 %); MoorDyn-C's is the largest segment
    // tension, -0.16 % below the end's with 20 segments, and its line settles from its initial shape (-0.28 %
    // measured with MoorDyn-C 2.7.1); taking H instead would miss by 1.6 %
    EXPECT_NEAR(c->lines()[0]->max_tension(0) / T, 1.0, 5.0e-3);
    c->advance(0, 0.0, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom);
    for (const auto& comp : cap.components()) { EXPECT_FALSE(comp.seen) << comp.name << ": level 0 is not the anchor level"; }
    const double dt = 0.125;
    for (int s = 0; s < 2; ++s) { c->advance(1, dt * s, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom); }
    for (const auto& comp : cap.components()) {
        EXPECT_TRUE(comp.exceeded) << comp.name;
        EXPECT_EQ(comp.first_time, dt) << comp.name << ": the end of the first anchor-level step";
        EXPECT_EQ(comp.peak_time == dt || comp.peak_time == 2 * dt, true) << comp.name;
        EXPECT_EQ(comp.demand >= comp.first_capacity, true);
    }
    EXPECT_NEAR(cap.component(0).demand, static_cast<double>(c->lines()[0]->max_tension(0)), 1.0e-9 * T);
    const auto ea = c->lines()[0]->end_force(0);
    EXPECT_NEAR(cap.component(1).demand, std::sqrt(double(ea[0]) * ea[0] + double(ea[1]) * ea[1] + double(ea[2]) * ea[2]), 1.0e-9 * T);
    // each first exceedance logged once, in the order the components were checked
    const auto rows = data_rows(dir + "/failures.dat");
    ASSERT_EQ(rows.size(), 3u);
    EXPECT_EQ(rows[0].rfind("0.125 conductor Tcap ", 0), 0u) << rows[0];
    EXPECT_EQ(rows[1].rfind("0.125 hardware Tcap_end_a ", 0), 0u) << rows[1];
    EXPECT_TRUE(std::filesystem::exists(dir + "/capacities.csv"));
    clear_capacities("Tcap");
}

TEST(Conductors, ACapacityAboveTheDemandIsNeverExceeded)
{
    const std::string dir = scratch("capacities_held");
    set_still_air(dir, "Tcaph");
    const Real T = still_air_end_tension();
    amrex::ParmParse("conductors.Tcaph").add("rated_breaking_strength", 1.01 * T);
    Mesh m(false);
    auto c = Conductors::create(0);
    c->set_ground(nullptr, m.geom);
    for (int s = 0; s < 4; ++s) { c->advance(0, 0.25 * s, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom); }
    const auto& comp = c->capacities().component(0);
    EXPECT_TRUE(comp.seen);
    EXPECT_FALSE(comp.exceeded);
    EXPECT_NEAR(comp.peak_ratio, 1.0 / 1.01, 5.0e-3);   // the demand's 0.5 % above (stub -0.03 %, MoorDyn-C -0.28 %)
    EXPECT_TRUE(data_rows(dir + "/failures.dat").empty()) << "the header only";
    clear_capacities("Tcaph");
}

TEST(Conductors, FailureStopEndsTheRunAtTheFirstExceedance)
{
    const std::string dir = scratch("capacities_stop");
    set_still_air(dir, "Tstop");
    const Real T = still_air_end_tension();
    amrex::ParmParse("conductors.Tstop").add("rated_breaking_strength", 0.99 * T);
    amrex::ParmParse("conductors").add("failure_stop", true);
    Mesh m(false);
    auto c = Conductors::create(0);
    c->set_ground(nullptr, m.geom);
    const std::string msg = erf_gtest::abort_message([&] { c->advance(0, 0.0, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom); });
    EXPECT_NE(msg.find("conductors.failure_stop"), std::string::npos) << msg;
    EXPECT_NE(msg.find("conductor Tstop at t = 0.25"), std::string::npos) << msg;
    // the logs are written before the run stops
    EXPECT_TRUE(std::filesystem::exists(dir + "/capacities.csv"));
    const auto rows = data_rows(dir + "/failures.dat");
    ASSERT_EQ(rows.size(), 1u);
    EXPECT_EQ(rows[0].rfind("0.25 conductor Tstop ", 0), 0u) << rows[0];
    clear_capacities("Tstop");
}

TEST(Conductors, FailureStopNeedsACapacity)
{
    const std::string dir = scratch("capacities_none");
    set_still_air(dir, "Tnone");
    amrex::ParmParse("conductors").add("failure_stop", true);
    Mesh m(false);
    auto c = Conductors::create(0);
    const std::string msg = erf_gtest::abort_message([&] { c->set_ground(nullptr, m.geom); });
    EXPECT_NE(msg.find("conductors.failure_stop = true needs a capacity"), std::string::npos) << msg;
    clear_capacities("Tnone");
}

TEST(Conductors, CapacityRecordsContinueAcrossARestart)
{
    const std::string dir = scratch("capacities_restart");
    set_still_air(dir, "Tcr");
    const Real T = still_air_end_tension();
    amrex::ParmParse("conductors.Tcr").add("rated_breaking_strength", 0.99 * T);
    Mesh m(false);
    const double dt = 0.25;
    auto a = Conductors::create(0);
    a->set_ground(nullptr, m.geom);
    for (int s = 0; s < 3; ++s) { a->advance(0, dt * s, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom); }
    const std::string chk = dir + "/chk00003";
    std::filesystem::create_directories(chk);
    a->write_checkpoint(chk);
    const auto saved = a->capacities().component(0);
    // the first run goes on past the checkpoint; the restarted run drops its rows after it
    for (int s = 3; s < 5; ++s) { a->advance(0, dt * s, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom); }
    // the restart adds the hardware's capacity: those components start at the restart
    amrex::ParmParse("conductors.Tcr").add("hardware_rated_load", 0.95 * T);
    auto b = Conductors::create(0);
    b->set_ground(nullptr, m.geom, chk);
    ASSERT_TRUE(b->restored());
    ASSERT_EQ(b->capacities().size(), 3u);
    const auto& r = b->capacities().component(0);
    EXPECT_TRUE(r.exceeded);
    EXPECT_EQ(r.first_time, saved.first_time);
    EXPECT_EQ(r.peak_ratio, saved.peak_ratio);
    EXPECT_FALSE(b->capacities().component(1).seen);
    for (int s = 3; s < 5; ++s) { b->advance(0, dt * s, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom); }
    EXPECT_EQ(b->capacities().component(0).first_time, dt) << "not reported again after the restart";
    const auto rows = data_rows(dir + "/failures.dat");
    ASSERT_EQ(rows.size(), 3u) << "the conductor's row from before the checkpoint, then the two dead ends' at the restart";
    EXPECT_EQ(rows[0].rfind("0.25 conductor Tcr ", 0), 0u) << rows[0];
    EXPECT_EQ(rows[1].rfind("1 hardware Tcr_end_a ", 0), 0u) << rows[1];
    // a restart with failure_stop: the conductor's exceedance restored from the checkpoint (0.25 s) does not stop it,
    // the dead ends' new one at the first step after it (1 s) does, and is the one named
    amrex::ParmParse("conductors").add("failure_stop", true);
    auto s = Conductors::create(0);
    s->set_ground(nullptr, m.geom, chk);
    const std::string msg = erf_gtest::abort_message([&] { s->advance(0, 3 * dt, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom); });
    EXPECT_NE(msg.find("hardware Tcr_end_a at t = 1 s"), std::string::npos) << msg;
    clear_capacities("Tcr");
}

TEST(Conductors, AStandingTowersFootingsAndAttachmentsAreChecked)
{
    // still air on a section over two towers that stand still: each tower's legs share its weight and the
    // lines' downward pull equally, (W - F_z) / 4 in compression, no uplift
    const std::string dir = scratch("capacities_towers");
    set_towered_section(dir, "Tft2");
    amrex::ParmParse("conductors").remove("prescribed_velocity");
    amrex::ParmParse("conductors").addarr("prescribed_velocity", std::vector<Real>{0.0, 0.0, 0.0});
    amrex::ParmParse pt("conductors.lat");
    pt.add("allowable_uplift", 1.0e5);
    // the legs carry about (9e4 + 4.8e3) / 4 = 2.37e4 N: an allowable of 2.3e4 N is exceeded, 1e5 N of uplift is not
    pt.add("allowable_compression", 2.3e4);
    amrex::ParmParse("conductors.Tft2").add("hardware_rated_load", 1.0e6);
    Mesh m(false);
    auto c = Conductors::create(0);
    c->set_ground(nullptr, m.geom);
    const auto& cap = c->capacities();
    // per line: four attachments (end_a, two towers, end_b); per tower: uplift and compression
    ASSERT_EQ(cap.size(), 4u + 2u * 2u);
    EXPECT_EQ(cap.component(2).name, "Tft2_tower2");
    EXPECT_EQ(cap.component(4).name, "Tft2_t1_footing_uplift");
    EXPECT_EQ(cap.component(5).name, "Tft2_t1_footing_compression");
    for (int s = 0; s < 2; ++s) { c->advance(0, 0.25 * s, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom); }
    for (int t = 0; t < 2; ++t) {
        const auto f = c->lines()[0]->tower_force(t);
        const double comp = (9.0e4 - double(f[2])) / 4.0;
        const auto& up = cap.component(4 + 2 * static_cast<std::size_t>(t));
        const auto& co = cap.component(5 + 2 * static_cast<std::size_t>(t));
        EXPECT_NEAR(up.demand, 0.0, 1.0e-6) << up.name;
        // 1e-6: MoorDyn-C's line leaves a moment of a few N m on the tower (2e-8 of the compression), the stub none
        EXPECT_NEAR(co.demand, comp, 1.0e-6 * comp) << co.name;
        EXPECT_TRUE(comp > 2.3e4) << "the case is set up for the compression to exceed its allowable";
        EXPECT_TRUE(co.exceeded) << co.name;
        EXPECT_FALSE(up.exceeded) << up.name;
        const auto& hw = cap.component(1 + static_cast<std::size_t>(t));
        EXPECT_NEAR(hw.demand, std::sqrt(double(f[0]) * f[0] + double(f[1]) * f[1] + double(f[2]) * f[2]), 1.0e-6) << hw.name;
    }
    pt.remove("allowable_uplift");
    pt.remove("allowable_compression");
    clear_capacities("Tft2");
    clear_towered_section("Tft2");
}

TEST(Conductors, AStandingTowerIsCheckedUnderTheLinesPullAtTheStepsEnd)
{
    // a crosswind blows the section out in the first step: a tower that stands still carries, over the next step,
    // the lines' pull where the first step left them (tower_foundation() after the next step starts), and that is
    // what the check after the first step must have taken, not the still-air pull the tower held during it
    const std::string dir = scratch("capacities_tower_end");
    set_towered_section(dir, "Tte");
    amrex::ParmParse pt("conductors.lat");
    pt.add("allowable_compression", 1.0e6);
    Mesh m(false);
    auto c = Conductors::create(0);
    c->set_ground(nullptr, m.geom);
    const auto& cap = c->capacities();
    ASSERT_EQ(cap.size(), 2u);
    const double still = c->tower_foundation(0).max_compression;
    c->advance(0, 0.0, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom);
    const double checked = cap.component(0).demand;
    c->advance(0, 0.25, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom);
    const double next = c->tower_foundation(0).max_compression;
    EXPECT_GT(std::abs(next - still), 1.0e-3 * still) << "the case must move the lines in the first step";
    EXPECT_NEAR(checked, next, 1.0e-9 * next);
    pt.remove("allowable_compression");
    clear_capacities("Tte");
    clear_towered_section("Tte");
}

TEST(Conductors, AMovingTowerIsCheckedAfterEveryCouplingStep)
{
    // towers that bend at 2 Hz: the lines and towers step in N = ceil(dt f 20) = 10 coupling steps of h = dt / 10 per
    // step of 0.25 s, and the check follows each of them: an allowable below the still-air compression is first
    // exceeded at the end of the first coupling step, t = h, not at the end of the ERF step
    const std::string dir = scratch("capacities_moving");
    set_towered_section(dir, "Tmc");
    amrex::ParmParse pt("conductors.lat");
    pt.add("frequency", 2.0);
    pt.add("allowable_compression", 1.0e3);
    Mesh m(false);
    auto c = Conductors::create(0);
    c->set_ground(nullptr, m.geom);
    ASSERT_TRUE(c->lines()[0]->towers_move());
    const auto& cap = c->capacities();
    ASSERT_EQ(cap.size(), 2u);
    const double dt = 0.25, h = dt / 10.0;
    c->advance(0, 0.0, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom);
    for (const auto& comp : cap.components()) {
        EXPECT_TRUE(comp.exceeded) << comp.name;
        EXPECT_NEAR(comp.first_time, h, 1.0e-12) << comp.name << ": the end of the first coupling step";
    }
    pt.remove("allowable_compression");
    clear_capacities("Tmc");
    clear_towered_section("Tmc");
}

namespace {
// a second section, 200 m from set_towered_section's, on its own towers of type lat2: lat's dimensions, bending at 8 Hz
void add_second_towered_section (const std::string& first, const std::string& name)
{
    amrex::ParmParse pp("conductors");
    pp.addarr("lines", std::vector<std::string>{first, name});
    pp.addarr("tower_types", std::vector<std::string>{"lat", "lat2"});
    amrex::ParmParse pt("conductors.lat2");
    pt.add("base_width", 6.0); pt.add("top_width", 1.5); pt.add("solidity", 0.2); pt.add("arm_length", 12.0);
    pt.add("arm_depth", 1.2); pt.add("weight", 9.0e4); pt.add("frequency", 8.0);
    amrex::ParmParse ps("conductors." + name);
    ps.addarr("end_a", std::vector<Real>{100.0, 300.0, 30.0});
    ps.addarr("end_b", std::vector<Real>{1000.0, 300.0, 30.0});
    ps.addarr("towers", std::vector<Real>{400.0, 300.0, 30.0, 700.0, 300.0, 30.0});
    ps.addarr("length", std::vector<Real>{301.5, 301.5, 301.5});
    ps.add("diameter", 0.0281);
    ps.add("mass_per_length", 1.628);
    ps.add("axial_stiffness", 3.0e7);
    ps.add("tower_type", std::string("lat2"));
}
void clear_second_towered_section (const std::string& name)
{
    amrex::ParmParse pt("conductors.lat2");
    for (const char* k : {"base_width", "top_width", "solidity", "arm_length", "arm_depth", "weight", "frequency"}) { pt.remove(k); }
    amrex::ParmParse ps("conductors." + name);
    for (const char* k : {"end_a", "end_b", "towers", "length", "diameter", "mass_per_length", "axial_stiffness", "tower_type",
                          "rated_breaking_strength"}) { ps.remove(k); }
}
} // namespace

TEST(Conductors, AStepsFirstExceedancesAreLoggedInTimeOrderAcrossGroups)
{
    // two sections on their own bending towers: two groups, Tga's (2 Hz: 10 coupling steps of 0.025 s per step of
    // 0.25 s) stepped before Tgb's (8 Hz: 40 of 0.00625 s). Both spans' ratings are far below their tension, so each is
    // first exceeded at its group's first coupling step: Tgb's at 0.00625 s, before Tga's at 0.025 s, though Tga's group
    // stepped first. The step's rows must be in time order (Tgb's first) and failure_stop must name Tgb.
    const std::string dir = scratch("capacities_order");
    set_towered_section(dir, "Tga");
    add_second_towered_section("Tga", "Tgb");
    amrex::ParmParse("conductors.lat").add("frequency", 2.0);
    amrex::ParmParse("conductors.Tga").add("rated_breaking_strength", 1.0);
    amrex::ParmParse("conductors.Tgb").add("rated_breaking_strength", 1.0);
    amrex::ParmParse("conductors").add("failure_stop", true);
    Mesh m(false);
    auto c = Conductors::create(0);
    c->set_ground(nullptr, m.geom);
    const double dt = 0.25;
    const std::string msg = erf_gtest::abort_message([&] { c->advance(0, 0.0, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom); });
    EXPECT_NE(msg.find("conductor Tgb_span1 at t = 0.00625"), std::string::npos) << msg;
    const auto rows = data_rows(dir + "/failures.dat");
    ASSERT_EQ(rows.size(), 6u) << "three spans of each line";
    double tprev = -1.0;
    for (const auto& r : rows) {
        const double t = std::stod(r);
        EXPECT_GE(t, tprev) << "rows in time order: " << r;
        tprev = t;
    }
    EXPECT_NE(rows.front().find(" Tgb_span1 "), std::string::npos) << rows.front();
    EXPECT_NEAR(std::stod(rows.front()), dt / 40.0, 1.0e-12);
    EXPECT_NE(rows.back().find(" Tga_span3 "), std::string::npos) << rows.back();
    EXPECT_NEAR(std::stod(rows.back()), dt / 10.0, 1.0e-12);
    clear_second_towered_section("Tgb");
    clear_capacities("Tga");
    clear_towered_section("Tga");
}

TEST(Conductors, ACircuitsStringsHardwareAndSharedFootingsAreChecked)
{
    // the circuit of set_circuit on towers that bend: C1's insulator strings, the shield wire CS's hardware on the shared
    // towers and each tower's footings, all checked at the end of every coupling step against capacities below their
    // still-air demands, so each is first exceeded at the end of the first coupling step, t = h
    const std::string dir = scratch("capacities_circuit");
    set_circuit(dir, 0.02);
    amrex::ParmParse("conductors.C1").add("insulator_rated_load", 1.0);
    amrex::ParmParse("conductors.CS").add("hardware_rated_load", 1.0);
    amrex::ParmParse("conductors.lat").add("allowable_compression", 1.0);
    Mesh m(true);
    auto c = Conductors::create(0);
    c->set_ground(m.znd.get(), m.geom);
    const auto& cap = c->capacities();
    // C1: two strings; CS: end_a, two towers, end_b; the two shared towers' compression, once each
    ASSERT_EQ(cap.size(), 2u + 4u + 2u);
    EXPECT_EQ(cap.component(0).name, "C1_string1");
    EXPECT_EQ(cap.component(3).name, "CS_tower1");
    EXPECT_EQ(cap.component(6).name, "C2_t1_footing_compression");
    const double dt = 0.25, h = dt / 10.0;   // 2 Hz: ceil(0.25 x 2 x 20) = 10 coupling steps
    c->advance(0, 0.0, dt, m.u, m.v, m.w, m.znd.get(), nullptr, m.geom);
    for (const auto& comp : cap.components()) {
        EXPECT_TRUE(comp.exceeded) << comp.name;
        EXPECT_NEAR(comp.first_time, h, 1.0e-12) << comp.name;
    }
    // the last records are the state the step ended at
    const auto& c1 = *c->lines()[0];
    const auto& cs = *c->lines()[3];
    ASSERT_EQ(c1.name(), "C1");
    ASSERT_EQ(cs.name(), "CS");
    EXPECT_NEAR(cap.component(1).demand, static_cast<double>(c1.insulator_tension(1)), 1.0e-9 * cap.component(1).demand);
    const auto f = cs.tower_force(0);
    EXPECT_NEAR(cap.component(3).demand, std::sqrt(double(f[0]) * f[0] + double(f[1]) * f[1] + double(f[2]) * f[2]),
                1.0e-9 * cap.component(3).demand);
    EXPECT_NEAR(cap.component(6).demand, static_cast<double>(c->tower_foundation(0).max_compression), 1.0e-9 * cap.component(6).demand)
        << "the tower carries all four lines";
    amrex::ParmParse("conductors.C1").remove("insulator_rated_load");
    amrex::ParmParse("conductors.CS").remove("hardware_rated_load");
    amrex::ParmParse("conductors.lat").remove("allowable_compression");
    clear_circuit();
}

TEST(Conductors, AFrameTowersMembersAreCheckedAgainstTheirStrength)
{
    // generated frame towers with angle design data: every member is a component, its demand-to-capacity ratio the
    // member's ASCE 10-15 utilisation on the side that governs, at the state the step ended at
    const std::string dir = scratch("capacities_members");
    set_towered_section(dir, "Tfm");
    amrex::ParmParse pt("conductors.lat");
    pt.remove("weight");
    pt.add("frame_panels", 8);
    pt.addarr("leg_angle", std::vector<Real>{0.15, 0.012});
    pt.addarr("brace_angle", std::vector<Real>{0.10, 0.008});
    pt.add("damping_ratio", 0.02);
    Mesh m(false);
    auto c = Conductors::create(0);
    c->set_ground(nullptr, m.geom);
    const auto checks0 = c->member_checks(0);
    ASSERT_FALSE(checks0.empty());
    const auto& cap = c->capacities();
    ASSERT_EQ(cap.size(), 2 * checks0.size()) << "every member of the two towers";
    EXPECT_EQ(cap.component(0).name, "Tfm_t1_m" + std::to_string(checks0[0].member));
    c->advance(0, 0.0, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom);
    for (std::size_t t = 0; t < 2; ++t) {
        const auto checks = c->member_checks(t);
        for (std::size_t k = 0; k < checks.size(); ++k) {
            const auto& comp = cap.component(t * checks0.size() + k);
            EXPECT_EQ(comp.kind, erf_conductors::ComponentKind::Member);
            EXPECT_NEAR(comp.demand / comp.capacity, checks[k].utilisation, 1.0e-12 + 1.0e-9 * checks[k].utilisation) << comp.name;
        }
    }
    for (const char* k : {"frame_panels", "leg_angle", "brace_angle", "damping_ratio"}) { pt.remove(k); }
    clear_capacities("Tfm");
    clear_towered_section("Tfm");
}

namespace {
// an inclined span in still air over a flat mesh, its end b 10 m above its end a: with MoorDyn-C the tension at b is
// higher by the conductor's weight over those 10 m (160 N), so the choice of the end does not hang on roundoff. The
// stub's parabola carries the same tension at both ends, a tie, which releases end a
int higher_end () { return erf_moordyn::is_stub() ? 0 : 1; }
std::string higher_end_name () { return erf_moordyn::is_stub() ? "span1_end_a" : "span1_end_b"; }
void set_break (const std::string& dir, const std::string& name)
{
    set_inputs(dir, name, false, {{300.0, 500.0, 30.0}}, {{600.0, 500.0, 40.0}});
    amrex::ParmParse pp("conductors");
    pp.addarr("prescribed_velocity", std::vector<Real>{0.0, 0.0, 0.0});
    pp.add("ground_contact", true);
    pp.add("failure_break", true);
}
void clear_break (const std::string& name)
{
    clear_capacities(name);
    amrex::ParmParse pp("conductors");
    pp.remove("ground_contact");
    pp.remove("failure_break");
}
} // namespace

TEST(Conductors, FailureBreakLetsGoAtTheHigherEndAndTheRunGoesOn)
{
    const std::string dir = scratch("break");
    set_break(dir, "Tbr");
    // the span pulls about 13 kN: a strength of 5 kN is exceeded at the first step
    amrex::ParmParse("conductors.Tbr").add("rated_breaking_strength", 5000.0);
    Mesh m(false);
    auto c = Conductors::create(0);
    c->set_ground(nullptr, m.geom);
    EXPECT_TRUE(data_rows(dir + "/breaks.dat").empty());
    const double dt = 0.25;
    c->advance(0, 0.0, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom);
    const auto& line = *c->lines().front();
    ASSERT_EQ(c->breaks().size(), 1u);
    const auto& b = c->breaks().record(0);
    EXPECT_EQ(b.time, dt) << "the end of the step of the exceedance";
    EXPECT_EQ(b.line, "Tbr");
    EXPECT_EQ(b.component, "Tbr");
    EXPECT_EQ(b.release.name(), higher_end_name()) << "the higher end carries the higher tension (the stub: a tie, end a)";
    EXPECT_GT(b.tension, 10000.0);
    if (!erf_moordyn::is_stub()) {
        // the tensions are still those of the step of the exceedance: b's is the higher, by the weight of 10 m of conductor
        EXPECT_NEAR(b.tension - static_cast<double>(line.tension_a()), 15.96 * 10.0, 25.0);
    }
    EXPECT_EQ(b.tension, static_cast<double>(higher_end() == 0 ? line.tension_a() : line.tension_b()));
    EXPECT_TRUE(line.span_end_attached(0, 1 - higher_end()));
    EXPECT_FALSE(line.span_end_attached(0, higher_end()));
    auto rows = data_rows(dir + "/breaks.dat");
    ASSERT_EQ(rows.size(), 1u);
    EXPECT_EQ(rows[0].rfind("0.25 Tbr " + higher_end_name() + " Tbr ", 0), 0u) << rows[0];
    {
        std::ifstream f(dir + "/breaks.dat");
        std::string header;
        std::getline(f, header);
        EXPECT_EQ(header, "time line released component tension");
    }
    // the run goes on: the span falls and comes to hang from end a
    const int nsteps = erf_moordyn::is_stub() ? 4 : 160;
    for (int s = 1; s <= nsteps; ++s) { c->advance(0, dt * s, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom); }
    EXPECT_EQ(c->breaks().size(), 1u) << "no second release";
    EXPECT_EQ(data_rows(dir + "/breaks.dat").size(), 1u);
    EXPECT_EQ(data_rows(dir + "/failures.dat").size(), 1u) << "the conductor's exceedance, once";
    EXPECT_LT(line.max_tension(), Real(1500.0)) << "30 or 40 m of conductor hang from the other end: 480 or 640 N";
    EXPECT_GE(line.nodes_on_ground(0), 15) << "most of the span lies on the ground";
    EXPECT_TRUE(c->ground_contacts().span(0).touched);
    const auto& comp = c->capacities().component(0);
    EXPECT_TRUE(comp.exceeded);
    EXPECT_EQ(comp.first_time, dt);
    EXPECT_GT(comp.peak_ratio, 2.0) << "the peak is the whole span's, before it let go";
    clear_break("Tbr");
}

TEST(Conductors, FailureBreakNeedsGroundContactALineCapacityAndNoFailureStop)
{
    const std::string dir = scratch("break_needs");
    set_break(dir, "Tbn");
    Mesh m(false);
    // no capacity on a line
    auto c = Conductors::create(0);
    std::string msg = erf_gtest::abort_message([&] { c->set_ground(nullptr, m.geom); });
    EXPECT_NE(msg.find("conductors.failure_break = true needs a line's capacity to check"), std::string::npos) << msg;
    amrex::ParmParse("conductors.Tbn").add("rated_breaking_strength", 5000.0);
    amrex::ParmParse pp("conductors");
    pp.remove("ground_contact");
    msg = erf_gtest::abort_message([&] { Conductors::create(0); });
    EXPECT_NE(msg.find("conductors.failure_break = true needs conductors.ground_contact = true"), std::string::npos) << msg;
    pp.add("ground_contact", true);
    pp.add("failure_stop", true);
    msg = erf_gtest::abort_message([&] { Conductors::create(0); });
    EXPECT_NE(msg.find("conductors.failure_break and conductors.failure_stop cannot both be true"), std::string::npos) << msg;
    pp.remove("failure_stop");
    // off by default: with the capacity alone the span stays whole and no release log is kept
    pp.remove("failure_break");
    auto o = Conductors::create(0);
    o->set_ground(nullptr, m.geom);
    o->advance(0, 0.0, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom);
    EXPECT_TRUE(o->breaks().empty());
    EXPECT_TRUE(o->lines().front()->span_end_attached(0, 0) && o->lines().front()->span_end_attached(0, 1));
    EXPECT_TRUE(o->capacities().component(0).exceeded);
    EXPECT_FALSE(std::filesystem::exists(dir + "/breaks.dat"));
    clear_break("Tbn");
}

TEST(Conductors, ASpanIsReleasedOnceInAStepAndADeadEndReleasesItsOwnEnd)
{
    Mesh m(false);
    {
        // the conductor and both dead ends over their capacities in the same step: the conductor's exceedance (the
        // first logged) releases the higher end; the dead end's at the other end would release the same span again,
        // the one at the released end finds it gone
        const std::string dir = scratch("break_once");
        set_break(dir, "Tbo");
        amrex::ParmParse ps("conductors.Tbo");
        ps.add("rated_breaking_strength", 5000.0);
        ps.add("hardware_rated_load", 5000.0);
        auto c = Conductors::create(0);
        c->set_ground(nullptr, m.geom);
        c->advance(0, 0.0, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom);
        EXPECT_EQ(data_rows(dir + "/failures.dat").size(), 3u);
        ASSERT_EQ(c->breaks().size(), 1u);
        EXPECT_EQ(c->breaks().record(0).release.name(), higher_end_name());
        EXPECT_EQ(c->breaks().record(0).component, "Tbo");
        EXPECT_TRUE(c->lines().front()->span_end_attached(0, 1 - higher_end()));
        c->advance(0, 0.25, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom);
        EXPECT_EQ(c->breaks().size(), 1u);
        clear_break("Tbo");
    }
    {
        // the hardware's capacity alone: end a's exceedance is logged first and releases end a, though b pulls harder
        const std::string dir = scratch("break_end");
        set_break(dir, "Tbd");
        amrex::ParmParse("conductors.Tbd").add("hardware_rated_load", 5000.0);
        auto c = Conductors::create(0);
        c->set_ground(nullptr, m.geom);
        c->advance(0, 0.0, 0.25, m.u, m.v, m.w, nullptr, nullptr, m.geom);
        ASSERT_EQ(c->breaks().size(), 1u);
        EXPECT_EQ(c->breaks().record(0).release.name(), "span1_end_a");
        EXPECT_EQ(c->breaks().record(0).component, "Tbd_end_a");
        EXPECT_FALSE(c->lines().front()->span_end_attached(0, 0));
        EXPECT_TRUE(c->lines().front()->span_end_attached(0, 1));
        clear_break("Tbd");
    }
}

TEST(Conductors, ReleasesContinueAcrossARestart)
{
    const std::string dir = scratch("break_restart");
    set_break(dir, "Tbrr");
    amrex::ParmParse("conductors.Tbrr").add("rated_breaking_strength", 5000.0);
    Mesh m(false);
    const double dt = 0.25;
    auto step = [&] (Conductors& c, int s) { c.advance(0, dt * s, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom); };
    auto a = Conductors::create(0);
    a->set_ground(nullptr, m.geom);
    step(*a, 0);
    ASSERT_EQ(a->breaks().size(), 1u);
    // checkpoints at the step of the release, before the released span has moved, and four steps later
    const std::string chk1 = dir + "/chk00001", chk5 = dir + "/chk00005";
    std::filesystem::create_directories(chk1);
    a->write_checkpoint(chk1);
    ASSERT_TRUE(std::filesystem::exists(chk1 + "/conductors/breaks"));
    for (int s = 1; s < 5; ++s) { step(*a, s); }
    std::filesystem::create_directories(chk5);
    a->write_checkpoint(chk5);
    const int last = 12;
    for (int s = 5; s < last; ++s) { step(*a, s); }
    const std::vector<Real> pa = a->lines().front()->node_positions();
    const Real ta = a->lines().front()->tension_a();
    struct Case { std::string chk; int from; bool with_break; };
    for (const Case& k : {Case{chk1, 1, true}, Case{chk5, 5, true}, Case{chk5, 5, false}}) {
        amrex::ParmParse pp("conductors");
        pp.remove("failure_break");
        if (k.with_break) { pp.add("failure_break", true); }
        auto b = Conductors::create(0);
        b->set_ground(nullptr, m.geom, k.chk);
        ASSERT_TRUE(b->restored());
        ASSERT_EQ(b->breaks().size(), 1u) << k.chk;
        const auto& ra = a->breaks().record(0);
        const auto& rb = b->breaks().record(0);
        EXPECT_EQ(rb.time, ra.time);
        EXPECT_EQ(rb.line, ra.line);
        EXPECT_TRUE(rb.release == ra.release);
        EXPECT_EQ(rb.component, ra.component);
        EXPECT_EQ(rb.tension, ra.tension) << "to the last bit";
        const auto& line = *b->lines().front();
        EXPECT_FALSE(line.span_end_attached(0, higher_end())) << "the release is made again";
        EXPECT_EQ(line.num_kinematics_points(), a->lines().front()->num_kinematics_points());
        for (int s = k.from; s < last; ++s) { step(*b, s); }
        const std::vector<Real> pb = line.node_positions();
        ASSERT_EQ(pb.size(), pa.size());
        Real worst = 0.0;
        for (std::size_t i = 0; i < pa.size(); ++i) { worst = std::max(worst, std::abs(pa[i] - pb[i])); }
        EXPECT_LE(worst, roundoff * Real(1000.0)) << "restarted from " << k.chk << (k.with_break ? "" : " without failure_break");
        EXPECT_NEAR(line.tension_a(), ta, roundoff * std::max(ta, Real(1.0)) * Real(100.0));
        EXPECT_EQ(b->breaks().size(), 1u) << "not released again";
    }
    // a released span needs the ground: a restart without ground_contact is refused, naming it
    {
        amrex::ParmParse pp("conductors");
        pp.remove("ground_contact");
        auto g = Conductors::create(0);
        const std::string msg = erf_gtest::abort_message([&] { g->set_ground(nullptr, m.geom, chk5); });
        EXPECT_NE(msg.find("holds 1 release(s) of the lines"), std::string::npos) << msg;
        EXPECT_NE(msg.find("the restart needs conductors.ground_contact = true"), std::string::npos) << msg;
        pp.add("ground_contact", true);
    }
    clear_break("Tbrr");
    // failure_break first set at a restart, after the conductor went over its strength: nothing is released
    const std::string dir2 = scratch("break_late");
    set_break(dir2, "Tbl");
    amrex::ParmParse("conductors.Tbl").add("rated_breaking_strength", 5000.0);
    amrex::ParmParse("conductors").remove("failure_break");
    auto w = Conductors::create(0);
    w->set_ground(nullptr, m.geom);
    step(*w, 0);
    ASSERT_TRUE(w->capacities().component(0).exceeded);
    const std::string chkw = dir2 + "/chk00001";
    std::filesystem::create_directories(chkw);
    w->write_checkpoint(chkw);
    EXPECT_FALSE(std::filesystem::exists(chkw + "/conductors/breaks")) << "no release log without the option";
    amrex::ParmParse("conductors").add("failure_break", true);
    auto l = Conductors::create(0);
    l->set_ground(nullptr, m.geom, chkw);
    for (int s = 1; s < 4; ++s) { step(*l, s); }
    EXPECT_TRUE(l->breaks().empty()) << "a release is made at a first exceedance only";
    EXPECT_TRUE(l->lines().front()->span_end_attached(0, 0) && l->lines().front()->span_end_attached(0, 1));
    clear_break("Tbl");
}

TEST(Conductors, StringsOverTheirLoadLetGoAndTheirTowersTakeNothing)
{
    const std::string dir = scratch("break_strings");
    set_towered_section(dir, "Tbs");
    amrex::ParmParse ps("conductors.Tbs");
    ps.add("insulator_length", 2.5);
    ps.add("insulator_mass", 60.0);
    // each string carries a span's weight, about 5 kN
    ps.add("insulator_rated_load", 1000.0);
    amrex::ParmParse pp("conductors");
    pp.remove("prescribed_velocity");
    pp.addarr("prescribed_velocity", std::vector<Real>{0.0, 0.0, 0.0});
    pp.add("ground_contact", true);
    pp.add("failure_break", true);
    Mesh m(false);
    auto c = Conductors::create(0);
    c->set_ground(nullptr, m.geom);
    const double dt = 0.25;
    c->advance(0, 0.0, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom);
    ASSERT_EQ(c->breaks().size(), 2u);
    EXPECT_EQ(c->breaks().record(0).release.name(), "string1");
    EXPECT_EQ(c->breaks().record(0).component, "Tbs_string1");
    EXPECT_EQ(c->breaks().record(1).release.name(), "string2");
    EXPECT_EQ(c->breaks().record(1).component, "Tbs_string2");
    EXPECT_GT(c->breaks().record(0).tension, 4000.0) << "the string's tension when it let go";
    const auto& line = *c->lines().front();
    EXPECT_FALSE(line.string_attached(0));
    EXPECT_FALSE(line.string_attached(1));
    EXPECT_TRUE(line.span_end_attached(0, 0) && line.span_end_attached(2, 1)) << "the dead ends hold";
    const int nsteps = erf_moordyn::is_stub() ? 2 : 12;
    for (int s = 1; s <= nsteps; ++s) { c->advance(0, dt * s, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom); }
    // the towers stand with no line on them
    for (std::size_t t = 0; t < 2; ++t) {
        const auto& F = c->towers()[t].line_force();
        EXPECT_EQ(F[0], Real(0.0)) << "tower " << t + 1;
        EXPECT_EQ(F[1], Real(0.0)) << "tower " << t + 1;
        EXPECT_EQ(F[2], Real(0.0)) << "tower " << t + 1;
    }
    EXPECT_EQ(c->breaks().size(), 2u);
    ps.remove("insulator_length"); ps.remove("insulator_mass");
    pp.remove("ground_contact"); pp.remove("failure_break");
    clear_capacities("Tbs");
    clear_towered_section("Tbs");
}

TEST(Conductors, LinesOnMovingTowersLetGoAtTheEndOfTheLevelsStep)
{
    const std::string dir = scratch("break_moving");
    set_towered_section(dir, "Tbm");
    make_towers_move(0.3);
    // the clamped section pulls about 13 kN on its dead ends and a span's weight, 4.8 kN, down on each tower
    amrex::ParmParse("conductors.Tbm").add("hardware_rated_load", 4000.0);
    amrex::ParmParse pp("conductors");
    pp.remove("prescribed_velocity");
    pp.addarr("prescribed_velocity", std::vector<Real>{0.0, 0.0, 0.0});
    pp.add("ground_contact", true);
    pp.add("failure_break", true);
    Mesh m(false);
    auto c = Conductors::create(0);
    c->set_ground(nullptr, m.geom);
    ASSERT_TRUE(c->lines().front()->towers_move());
    const double dt = 0.25;
    c->advance(0, 0.0, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom);
    // the four attachments are over their load after the first coupling step, before the level's step ends
    int exceeded = 0;
    for (const auto& comp : c->capacities().components()) {
        if (comp.kind != erf_conductors::ComponentKind::Hardware) { continue; }
        ++exceeded;
        EXPECT_TRUE(comp.exceeded) << comp.name;
        EXPECT_GT(comp.first_time, 0.0) << comp.name;
        EXPECT_LT(comp.first_time, dt) << comp.name << ": found inside the step";
    }
    EXPECT_EQ(exceeded, 4);
    // end_a releases span 1 there; tower 1 would release span 1 again (a tie: the span before it); tower 2 releases
    // span 2 at it; end_b releases span 3
    ASSERT_EQ(c->breaks().size(), 3u);
    const char* expected[3] = {"span1_end_a", "span2_end_b", "span3_end_b"};
    for (std::size_t b = 0; b < 3; ++b) {
        EXPECT_EQ(c->breaks().record(b).time, dt) << "made at the end of the level's step";
        EXPECT_EQ(c->breaks().record(b).release.name(), expected[b]);
    }
    // the run goes on with the lines released and the towers free of them where they let go
    for (int s = 1; s <= 3; ++s) { c->advance(0, dt * s, dt, m.u, m.v, m.w, nullptr, nullptr, m.geom); }
    EXPECT_EQ(c->breaks().size(), 3u);
    pp.remove("ground_contact"); pp.remove("failure_break");
    clear_capacities("Tbm");
    clear_towered_section("Tbm");
}
