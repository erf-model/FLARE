#ifndef ERF_FIRE_REGION_GHOST_FILL_H_
#define ERF_FIRE_REGION_GHOST_FILL_H_

/**
 * @file ERF_FireRegionGhostFill.H
 * @brief Ghost cells at the edge of a fire region that does not cover its level.
 *
 * On a fine level that covers part of the domain, a MultiFab on the region's
 * boxes has ghost cells inside the domain that no box covers. FillBoundary has
 * nothing to copy into them and fire_fill_boundary (ERF_FireGrid.H) fills only
 * those outside a non-periodic face of the domain, so they keep whatever was
 * set before the exchange.
 *
 * fire_fill_region_ghosts() exchanges the ghost cells as fire_fill_boundary
 * does and then extends the field into each ghost cell that neither exchange
 * reached from the edge of the region, as fire_fill_boundary does across a
 * non-periodic face of the domain: the value is the one a single box over the
 * whole (rectangular) region would get, whatever the box layout. The ghost cells are marked with `unset` before the exchange, so
 * `unset` must be a value the field never takes.
 */

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <ERF_FireGrid.H>

inline void fire_fill_region_ghosts (amrex::MultiFab& mf, const amrex::Geometry& geom, amrex::Real unset)
{
    const int ncomp = mf.nComp();
    // mark the ghost cells, keeping the valid ones
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box& vbx = mfi.validbox();
        auto a = mf.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (!vbx.contains(i, j, k)) { a(i, j, k, n) = unset; }
        });
    }
    fire_fill_boundary(mf, geom);
    // One direction at a time: a ghost still unset takes the cell of its own
    // box's range in x (which may be a ghost FillBoundary filled from a
    // neighbouring box), then y, then z. The answer is that of one box over
    // the whole region, whatever the box layout: a corner ghost next to a
    // neighbouring box takes that box's cell, not the diagonal one of its own.
    // Each pass reads only cells inside its direction's range and writes only
    // cells outside it, so no cell is read and written in the same pass.
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const amrex::Box& vbx = mfi.validbox();
            const int lo = vbx.smallEnd(dir);
            const int hi = vbx.bigEnd(dir);
            auto a = mf.array(mfi);
            amrex::ParallelFor(mfi.fabbox(), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (a(i, j, k, n) != unset) { return; }
                int c[3] = {i, j, k};
                c[dir] = amrex::min(amrex::max(c[dir], lo), hi);
                const amrex::Real v = a(c[0], c[1], c[2], n);
                if (v != unset) { a(i, j, k, n) = v; }
            });
        }
    }
}

#endif
