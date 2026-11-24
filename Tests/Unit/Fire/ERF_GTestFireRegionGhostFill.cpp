#include <gtest/gtest.h>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_BoxList.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "ERF_FireRegionGhostFill.H"

/**
 * @file ERF_GTestFireRegionGhostFill.cpp
 * @brief fire_fill_region_ghosts on a region of two boxes that covers part of
 *        a 16 x 8 slab: every ghost cell takes the value one box over the
 *        whole region would give it (the nearest cell of the region), whatever
 *        box it belongs to; in particular a corner ghost next to the other box
 *        takes that box's cell, not the diagonal cell of its own box. (The open
 *        fraction and roof height of build_open_fraction kept their defaults
 *        at the edge of such a region until Nov 2026.)
 */

using namespace amrex;

namespace {
Real field (int i, int j) { return Real(i + 100 * j); }

struct Region
{
    Geometry geom; BoxArray ba; DistributionMapping dm; MultiFab mf;
    /// Boxes [ilo, ilo + 3] and [ilo + 4, ilo + 7] by rows [jlo, jlo + 3] of a 16 x 8 slab
    Region (int ilo, int jlo, int per_x)
    {
        const Box domain(IntVect(0, 0, 0), IntVect(15, 7, 0));
        geom = Geometry(domain, RealBox(0.0, 0.0, 0.0, 16.0, 8.0, 1.0), CoordSys::cartesian, {per_x, 0, 0});
        BoxList bl;
        bl.push_back(Box(IntVect(ilo, jlo, 0), IntVect(ilo + 3, jlo + 3, 0)));
        bl.push_back(Box(IntVect(ilo + 4, jlo, 0), IntVect(ilo + 7, jlo + 3, 0)));
        ba = BoxArray(bl);
        dm = DistributionMapping(ba);
        mf.define(ba, dm, 1, 1);
        mf.setVal(1.0_rt);   // build_open_fraction's default before the fill
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            auto a = mf.array(mfi);
            const Box& vb = mfi.validbox();
            for (int j = vb.smallEnd(1); j <= vb.bigEnd(1); ++j)
                for (int i = vb.smallEnd(0); i <= vb.bigEnd(0); ++i) { a(i, j, 0) = field(i, j); }
        }
        fire_fill_region_ghosts(mf, geom, -1.0_rt);
    }
};
} // namespace

TEST(FireRegionGhostFill, TheEdgeOfAPartialRegionTakesItsEdgeColumn)
{
    // region x 2-9, y 2-5, inside the domain on every side
    Region r(2, 2, 0);
    for (MFIter mfi(r.mf); mfi.isValid(); ++mfi) {
        auto a = r.mf.const_array(mfi);
        const Box gb = mfi.fabbox();
        for (int k = gb.smallEnd(2); k <= gb.bigEnd(2); ++k)
        for (int j = gb.smallEnd(1); j <= gb.bigEnd(1); ++j)
            for (int i = gb.smallEnd(0); i <= gb.bigEnd(0); ++i) {
                // as one box over the region x 2-9, y 2-5, whatever the box a ghost belongs to
                const int ic = amrex::min(amrex::max(i, 2), 9);
                const int jc = amrex::min(amrex::max(j, 2), 5);
                EXPECT_EQ(a(i, j, k), field(ic, jc)) << "cell " << i << " " << j << " " << k;
            }
    }
}

TEST(FireRegionGhostFill, AnUncoveredPeriodicImageTakesTheEdgeColumn)
{
    // region x 0-7 against the low face, periodic in x, so the ghost column at
    // i = -1 is covered by nothing (the image of i = 15 is outside the region)
    // and takes the edge column; non-periodic, the same as fire_fill_boundary
    for (int per : {0, 1}) {
        Region r(0, 0, per);
        for (MFIter mfi(r.mf); mfi.isValid(); ++mfi) {
            auto a = r.mf.const_array(mfi);
            const Box& vb = mfi.validbox();
            if (vb.smallEnd(0) != 0) { continue; }
            for (int j = -1; j <= 4; ++j) {
                const int jc = amrex::min(amrex::max(j, 0), 3);
                EXPECT_EQ(a(-1, j, 0), field(0, jc)) << "periodic " << per << " row " << j;
            }
        }
    }
}
