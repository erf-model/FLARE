#include <ERF_FrozenWind.H>

#include <AMReX_FileSystem.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>
#include <AMReX_VisMF.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <sstream>

using namespace amrex;

namespace {

[[noreturn]] void frozen_abort (const std::string& dir, const std::string& what)
{
    amrex::Abort("erf.fire.frozen_wind.file = " + dir + ": " + what);
    std::abort(); // amrex::Abort is not declared noreturn
}

/// The first word of <dir>/Header, read on the I/O rank and broadcast
std::string header_magic (const std::string& dir)
{
    const std::string hdr = dir + "/Header";
    if (!amrex::FileSystem::Exists(hdr)) {
        frozen_abort(dir, "no Header file: not a GUST plotfile or checkpoint directory");
    }
    Vector<char> buf;
    ParallelDescriptor::ReadAndBcastFile(hdr, buf);
    std::istringstream is(std::string(buf.dataPtr()));
    std::string word;
    is >> word;
    return word;
}

/// At least one cell in each direction and a positive extent
void check_extent (const std::string& dir, const Box& dom, const RealBox& rb)
{
    for (int d = 0; d < 3; ++d) {
        if (!dom.ok() || dom.length(d) < 1 || !(rb.hi(d) > rb.lo(d))) {
            frozen_abort(dir, "the mesh needs at least one cell and prob_hi > prob_lo in every direction");
        }
    }
}

/// A MultiFab's BoxArray covers exactly the cells of dom (after conversion to cells)
bool covers_exactly (const MultiFab& mf, const Box& dom)
{
    const BoxArray ba = amrex::convert(mf.boxArray(), IndexType::TheCellType());
    return ba.minimalBox() == dom && ba.numPts() == dom.numPts();
}

} // namespace

void FrozenWind::read (const std::string& dir)
{
    m_file = dir;
    const std::string magic = header_magic(dir);
    if (magic == "HyperCLaw-V1.1") {
        m_format = "plotfile";
        read_plotfile(dir);
    } else if (magic == "gust-checkpoint" || magic == "ablsteady-checkpoint") {
        m_format = "checkpoint";
        read_checkpoint(dir);
    } else {
        frozen_abort(dir, "the Header starts with '" + magic + "', not a GUST plotfile (HyperCLaw-V1.1)"
                          " or checkpoint (gust-checkpoint)");
    }

    // Every value finite, heights rising with k, the lowest centre above the ground
    for (int c = 0; c < NCOMP; ++c) {
        if (m_cc.contains_nan(c, 1) || m_cc.contains_inf(c, 1)) {
            static const char* names[NCOMP] = {"u", "v", "w", "z", "ground"};
            frozen_abort(dir, std::string("the field ") + names[c] + " holds a value that is not finite");
        }
    }
    const int kz0 = m_geom.Domain().smallEnd(2);
    const int kz1 = m_geom.Domain().bigEnd(2);
    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    for (MFIter mfi(m_cc); mfi.isValid(); ++mfi) {
        const Box bx2 = makeSlab(mfi.validbox(), 2, kz0);
        const auto a = m_cc.const_array(mfi);
        reduce_op.eval(bx2, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int) -> GpuTuple<Real>
            {
                Real dz_min = a(i, j, kz0, Z) - a(i, j, kz0, ZG);
                for (int k = kz0; k < kz1; ++k) {
                    dz_min = amrex::min(dz_min, a(i, j, k + 1, Z) - a(i, j, k, Z));
                }
                return {dz_min};
            });
    }
    Real dz_min = amrex::get<0>(reduce_data.value(reduce_op));
    ParallelDescriptor::ReduceRealMin(dz_min);
    if (!(dz_min > Real(0.0))) {
        frozen_abort(dir, "the cell-centre heights must rise with k from above the ground in every column");
    }

    const Box& dom = m_geom.Domain();
    amrex::Print() << "[FIRE] Frozen wind: read the GUST " << m_format << " " << dir
                   << " (" << dom.length(0) << " x " << dom.length(1) << " x " << dom.length(2)
                   << " cells, x " << m_geom.ProbLo(0) << " to " << m_geom.ProbHi(0)
                   << " m, y " << m_geom.ProbLo(1) << " to " << m_geom.ProbHi(1)
                   << " m); the fire ignores the atmosphere's wind\n";
}

void FrozenWind::read_plotfile (const std::string& dir)
{
    PlotFileData pf(dir);
    if (pf.spaceDim() != 3) {
        frozen_abort(dir, "the plotfile is not 3D");
    }
    if (pf.finestLevel() != 0) {
        frozen_abort(dir, "the plotfile has more than one level; GUST writes one");
    }
    const auto& names = pf.varNames();
    for (const char* want : {"u", "v", "w", "z"}) {
        if (std::find(names.begin(), names.end(), want) == names.end()) {
            frozen_abort(dir, std::string("the plotfile has no variable '") + want + "' (a GUST solution plotfile has u v w p theta tke z)");
        }
    }
    const std::string znd_file = dir + "/Level_0/z_phys_nd";
    if (!amrex::FileSystem::Exists(znd_file + "_H")) {
        frozen_abort(dir, "no Level_0/z_phys_nd (GUST writes the node heights next to every plotfile)");
    }

    const Box dom = pf.probDomain(0);
    const auto plo = pf.probLo();
    const auto phi = pf.probHi();
    const RealBox rb(plo.data(), phi.data());
    check_extent(dir, dom, rb);
    m_geom.define(dom, rb, 0, {AMREX_D_DECL(0, 0, 0)});

    const MultiFab u = pf.get(0, "u");
    const MultiFab v = pf.get(0, "v");
    const MultiFab w = pf.get(0, "w");
    const MultiFab z = pf.get(0, "z");
    MultiFab z_nd;
    VisMF::Read(z_nd, znd_file);
    if (z_nd.ixType() != IndexType::TheNodeType()) {
        frozen_abort(dir, "Level_0/z_phys_nd is not nodal");
    }
    if (!covers_exactly(u, dom) || !covers_exactly(z_nd, dom)) {
        frozen_abort(dir, "the fields do not cover the plotfile's domain exactly");
    }
    build_columns(u, v, w, z, z_nd);
}

void FrozenWind::read_checkpoint (const std::string& dir)
{
    Vector<char> buf;
    ParallelDescriptor::ReadAndBcastFile(dir + "/Header", buf);
    std::istringstream is(std::string(buf.dataPtr()));

    auto expect_key = [&] (const char* want) {
        std::string key;
        is >> key;
        if (!is || key != want) {
            frozen_abort(dir, std::string("corrupt checkpoint Header: expected '") + want + "'");
        }
    };
    auto check_stream = [&] (const char* what) {
        if (!is) { frozen_abort(dir, std::string("corrupt checkpoint Header: bad value for '") + what + "'"); }
    };
    auto read_real = [&] (const char* what) {
        std::string tok;
        is >> tok;
        char* end = nullptr;
        const double x = std::strtod(tok.c_str(), &end);
        if (!is || tok.empty() || *end != '\0') {
            frozen_abort(dir, std::string("corrupt checkpoint Header: bad value for '") + what + "'");
        }
        return static_cast<Real>(x);
    };

    std::string magic;
    int version = 0;
    is >> magic >> version;
    if (!is || version != 1) {
        frozen_abort(dir, "a GUST checkpoint of format version " + std::to_string(version)
                          + "; this reader knows version 1");
    }
    int real_bytes = 0;
    expect_key("real_bytes");
    is >> real_bytes;
    check_stream("real_bytes");
    if (real_bytes != static_cast<int>(sizeof(Real))) {
        frozen_abort(dir, "the checkpoint holds " + std::to_string(real_bytes) + "-byte reals; this build uses "
                          + std::to_string(sizeof(Real)) + " (read GUST's plotfile instead)");
    }
    int step = 0;
    expect_key("step");
    is >> step;
    check_stream("step");
    IntVect n_cell;
    expect_key("n_cell");
    is >> n_cell[0] >> n_cell[1] >> n_cell[2];
    check_stream("n_cell");
    std::array<Real, 3> lo{}, hi{};
    expect_key("prob_lo");
    for (int d = 0; d < 3; ++d) { lo[d] = read_real("prob_lo"); }
    expect_key("prob_hi");
    for (int d = 0; d < 3; ++d) { hi[d] = read_real("prob_hi"); }

    const Box dom(IntVect(0), n_cell - IntVect(1));
    const RealBox rb(lo.data(), hi.data());
    check_extent(dir, dom, rb);
    m_geom.define(dom, rb, 0, {AMREX_D_DECL(0, 0, 0)});

    // The faces and nodes as written (VisMF defines each from its own header)
    MultiFab uf, vf, wf, z_nd_f;
    VisMF::Read(uf, MultiFabFileFullPrefix(0, dir, "Level_", "u"));
    VisMF::Read(vf, MultiFabFileFullPrefix(0, dir, "Level_", "v"));
    VisMF::Read(wf, MultiFabFileFullPrefix(0, dir, "Level_", "w"));
    VisMF::Read(z_nd_f, MultiFabFileFullPrefix(0, dir, "Level_", "z_phys_nd"));
    if (uf.ixType() != IndexType(IntVect(1, 0, 0)) || vf.ixType() != IndexType(IntVect(0, 1, 0))
        || wf.ixType() != IndexType(IntVect(0, 0, 1)) || z_nd_f.ixType() != IndexType::TheNodeType()) {
        frozen_abort(dir, "Level_0/u, v, w must be on x, y, z faces and z_phys_nd on nodes");
    }
    if (uf.nComp() != 1 || vf.nComp() != 1 || wf.nComp() != 1 || z_nd_f.nComp() != 1) {
        frozen_abort(dir, "Level_0/u, v, w and z_phys_nd must each hold one component");
    }
    if (!covers_exactly(uf, dom) || !covers_exactly(vf, dom) || !covers_exactly(wf, dom)
        || !covers_exactly(z_nd_f, dom)) {
        frozen_abort(dir, "the fields do not cover the checkpoint's n_cell exactly");
    }

    // One cell layout; the face and node data copied onto it (valid points only)
    const BoxArray ba = amrex::convert(uf.boxArray(), IndexType::TheCellType());
    const DistributionMapping dm(ba);
    auto on_layout = [&] (const MultiFab& src) {
        MultiFab dst(amrex::convert(ba, src.ixType()), dm, 1, 0);
        dst.ParallelCopy(src, 0, 0, 1);
        return dst;
    };
    const MultiFab xf = on_layout(uf);
    const MultiFab yf = on_layout(vf);
    const MultiFab zf = on_layout(wf);
    const MultiFab z_nd = on_layout(z_nd_f);

    MultiFab u(ba, dm, 1, 0), v(ba, dm, 1, 0), w(ba, dm, 1, 0), z(ba, dm, 1, 0);
    for (MFIter mfi(u, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto xa = xf.const_array(mfi);
        const auto ya = yf.const_array(mfi);
        const auto za = zf.const_array(mfi);
        const auto nd = z_nd.const_array(mfi);
        const auto ua = u.array(mfi);
        const auto va = v.array(mfi);
        const auto wa = w.array(mfi);
        const auto zc = z.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            ua(i, j, k) = Real(0.5) * (xa(i, j, k) + xa(i + 1, j, k));
            va(i, j, k) = Real(0.5) * (ya(i, j, k) + ya(i, j + 1, k));
            wa(i, j, k) = Real(0.5) * (za(i, j, k) + za(i, j, k + 1));
            // the mean of the eight nodes, as GUST's plotfile writes z
            zc(i, j, k) = Real(0.125) * (nd(i, j, k)         + nd(i + 1, j, k)
                                       + nd(i, j + 1, k)     + nd(i + 1, j + 1, k)
                                       + nd(i, j, k + 1)     + nd(i + 1, j, k + 1)
                                       + nd(i, j + 1, k + 1) + nd(i + 1, j + 1, k + 1));
        });
    }
    Gpu::streamSynchronize();   // xf, yf, zf are freed on return
    build_columns(u, v, w, z, z_nd);
}

void FrozenWind::define_from_fields (const Geometry& geom, const std::string& label,
                                     const MultiFab& u, const MultiFab& v, const MultiFab& w,
                                     const MultiFab& z, const MultiFab& z_nd, int max_box)
{
    m_file   = label;
    m_format = "fields";
    m_geom   = geom;
    build_columns(u, v, w, z, z_nd, max_box);
}

void FrozenWind::build_columns (const MultiFab& u, const MultiFab& v, const MultiFab& w,
                                const MultiFab& z, const MultiFab& z_nd, int max_box)
{
    const Box& dom = m_geom.Domain();
    const int nz = dom.length(2);
    if (nz < 2) {
        frozen_abort(m_file, "the mesh needs at least two cells in z");
    }
    BoxArray ba(dom);
    AMREX_ALWAYS_ASSERT(max_box >= 1);
    ba.maxSize(IntVect(AMREX_D_DECL(max_box, max_box, nz)));   // full-height columns, never split in z
    const DistributionMapping dm(ba);

    m_cc.define(ba, dm, NCOMP, 0);
    m_cc.ParallelCopy(u, 0, U, 1);
    m_cc.ParallelCopy(v, 0, V, 1);
    m_cc.ParallelCopy(w, 0, W, 1);
    m_cc.ParallelCopy(z, 0, Z, 1);

    MultiFab nd(amrex::convert(ba, IndexType::TheNodeType()), dm, 1, 0);
    nd.ParallelCopy(z_nd, 0, 0, 1);
    const int kz0 = dom.smallEnd(2);
    for (MFIter mfi(m_cc, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto a = m_cc.array(mfi);
        const auto n = nd.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i, j, k, ZG) = Real(0.25) * (n(i, j, kz0) + n(i + 1, j, kz0) + n(i, j + 1, kz0) + n(i + 1, j + 1, kz0));
        });
    }
    Gpu::streamSynchronize();   // nd is freed on return
}

void FrozenWind::sample_reference_wind (MultiFab& wind_ref,
                                        MultiFab& extract_z,
                                        MultiFab* ground,
                                        const Geometry& fire_geom,
                                        Real z_ref,
                                        const MultiFab* z_ref_cell,
                                        int wind_interp,
                                        Real z_sample,
                                        Real z0_log,
                                        bool below_first_log) const
{
    AMREX_ALWAYS_ASSERT(defined());
    const Box& gdom = m_geom.Domain();
    const int gilo = gdom.smallEnd(0), gihi = gdom.bigEnd(0);
    const int gjlo = gdom.smallEnd(1), gjhi = gdom.bigEnd(1);
    const int kz0  = gdom.smallEnd(2), kz1  = gdom.bigEnd(2);
    const Real gx0 = m_geom.ProbLo(0), gy0 = m_geom.ProbLo(1);
    const Real gdx = m_geom.CellSize(0), gdy = m_geom.CellSize(1);

    const Box& fdom = fire_geom.Domain();
    const int fi0 = fdom.smallEnd(0), fj0 = fdom.smallEnd(1);
    const Real fx0 = fire_geom.ProbLo(0), fy0 = fire_geom.ProbLo(1);
    const Real fdx = fire_geom.CellSize(0), fdy = fire_geom.CellSize(1);

    // The fire region must lie inside GUST's horizontal extent
    for (int d = 0; d < 2; ++d) {
        const Real tol = Real(1.0e-6) * amrex::max(Real(1.0), m_geom.ProbLength(d));
        if (fire_geom.ProbLo(d) < m_geom.ProbLo(d) - tol || fire_geom.ProbHi(d) > m_geom.ProbHi(d) + tol) {
            std::ostringstream os;
            os << "the GUST mesh spans " << (d == 0 ? "x" : "y") << " = " << m_geom.ProbLo(d) << " to "
               << m_geom.ProbHi(d) << " m but the fire region spans " << fire_geom.ProbLo(d) << " to "
               << fire_geom.ProbHi(d) << " m; GUST's domain must cover the fire region (same x, y origin)";
            frozen_abort(m_file, os.str());
        }
    }

    // The GUST columns each fire box reads, with a cell of margin, on that box's rank
    BoxList bl;
    const BoxArray& fba = wind_ref.boxArray();
    for (int n = 0; n < static_cast<int>(fba.size()); ++n) {
        const Box& fb = fba[n];
        const Real xlo = fx0 + Real(fb.smallEnd(0) - fi0) * fdx;
        const Real xhi = fx0 + Real(fb.bigEnd(0) + 1 - fi0) * fdx;
        const Real ylo = fy0 + Real(fb.smallEnd(1) - fj0) * fdy;
        const Real yhi = fy0 + Real(fb.bigEnd(1) + 1 - fj0) * fdy;
        const int ilo = amrex::max(gilo, gilo + static_cast<int>(std::floor((xlo - gx0) / gdx - Real(0.5))) - 1);
        const int ihi = amrex::min(gihi, gilo + static_cast<int>(std::floor((xhi - gx0) / gdx - Real(0.5))) + 2);
        const int jlo = amrex::max(gjlo, gjlo + static_cast<int>(std::floor((ylo - gy0) / gdy - Real(0.5))) - 1);
        const int jhi = amrex::min(gjhi, gjlo + static_cast<int>(std::floor((yhi - gy0) / gdy - Real(0.5))) + 2);
        bl.push_back(Box(IntVect(AMREX_D_DECL(ilo, jlo, kz0)), IntVect(AMREX_D_DECL(ihi, jhi, kz1))));
    }
    const BoxArray lba(std::move(bl));
    MultiFab loc(lba, wind_ref.DistributionMap(), 4, 0);   // u, v, z, ground
    loc.ParallelCopy(m_cc, U, 0, 2);
    loc.ParallelCopy(m_cc, Z, 2, 2);

    const bool bilinear = (wind_interp == 1);
    const bool resample = (z_sample > Real(0.0));
    const bool per_cell = (z_ref_cell != nullptr);
    const bool want_ground = (ground != nullptr);

    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    for (MFIter mfi(wind_ref, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto a   = loc.const_array(mfi);
        const auto fw  = wind_ref.array(mfi);
        const auto ez  = extract_z.array(mfi);
        const Array4<Real> gr = want_ground ? ground->array(mfi) : Array4<Real>{};
        const Array4<const Real> zc = per_cell ? z_ref_cell->const_array(mfi) : Array4<const Real>{};
        FrozenColumnsView cols;
        cols.a = a;
        cols.gilo = gilo; cols.gihi = gihi; cols.gjlo = gjlo; cols.gjhi = gjhi;
        cols.kz0 = kz0; cols.kz1 = kz1;
        cols.gx0 = gx0; cols.gy0 = gy0; cols.gdx = gdx; cols.gdy = gdy;

        reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GpuTuple<Real>
            {
                const Real z_cell = per_cell ? zc(i, j, k) : z_ref;
                constexpr Real z_tiny = Real(1.0e-30);
                const Real z0s   = std::abs(z0_log) + z_tiny;
                const Real z_at  = resample ? z_sample : z_cell;
                const Real scale = resample ? std::log((std::abs(z_cell) + z_tiny) / z0s)
                                              / (std::log((std::abs(z_sample) + z_tiny) / z0s) + z_tiny)
                                            : Real(1.0);

                const Real x = fx0 + (Real(i - fi0) + Real(0.5)) * fdx;
                const Real y = fy0 + (Real(j - fj0) + Real(0.5)) * fdy;

                Real u_sum, v_sum, z_sum, g_sum, excess;
                frozen_point_wind(cols, x, y, z_at, bilinear, below_first_log, z0_log,
                                  u_sum, v_sum, z_sum, g_sum, excess);
                fw(i, j, k, 0) = scale * u_sum;
                fw(i, j, k, 1) = scale * v_sum;
                ez(i, j, k)    = z_sum;
                if (want_ground) { gr(i, j, k) = g_sum; }
                return {excess};
            });
    }
    Real excess = amrex::get<0>(reduce_data.value(reduce_op));
    ParallelDescriptor::ReduceRealMax(excess);
    if (excess > Real(0.0)) {
        std::ostringstream os;
        os << "the wind sampling height (erf.fire.wind_ref_ht, erf.fire.wind_sample_ht or the per-fuel height)"
           << " lies " << excess << " m above the highest cell centre of a GUST column under the fire region;"
           << " lower it by more than that (GUST's top centre is lowest above high ground)";
        frozen_abort(m_file, os.str());
    }
}

void FrozenWind::gather_columns (FrozenWindHost& h, const Geometry& fire_geom, bool with_w) const
{
    AMREX_ALWAYS_ASSERT(defined());
    const Box& gdom = m_geom.Domain();
    const int gilo = gdom.smallEnd(0), gihi = gdom.bigEnd(0);
    const int gjlo = gdom.smallEnd(1), gjhi = gdom.bigEnd(1);
    const Real gx0 = m_geom.ProbLo(0), gy0 = m_geom.ProbLo(1);
    const Real gdx = m_geom.CellSize(0), gdy = m_geom.CellSize(1);

    // The columns the fire region reads, with the margin of sample_reference_wind
    const Real xlo = fire_geom.ProbLo(0), xhi = fire_geom.ProbHi(0);
    const Real ylo = fire_geom.ProbLo(1), yhi = fire_geom.ProbHi(1);
    const int ilo = amrex::max(gilo, gilo + static_cast<int>(std::floor((xlo - gx0) / gdx - Real(0.5))) - 1);
    const int ihi = amrex::min(gihi, gilo + static_cast<int>(std::floor((xhi - gx0) / gdx - Real(0.5))) + 2);
    const int jlo = amrex::max(gjlo, gjlo + static_cast<int>(std::floor((ylo - gy0) / gdy - Real(0.5))) - 1);
    const int jhi = amrex::min(gjhi, gjlo + static_cast<int>(std::floor((yhi - gy0) / gdy - Real(0.5))) + 2);
    const Box region(IntVect(AMREX_D_DECL(ilo, jlo, gdom.smallEnd(2))),
                     IntVect(AMREX_D_DECL(ihi, jhi, gdom.bigEnd(2))));
    AMREX_ALWAYS_ASSERT(region.ok() && gdom.contains(region));

    const int  nc   = with_w ? 5 : 4;   // u, v, z, ground (, w)
    const Long npts = region.numPts();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nc * npts <= Long(std::numeric_limits<int>::max()),
                                     "erf.fire.spotting.drift_wind = frozen or frozen_uvw: the GUST columns under the fire region are too many to copy to every rank");
    h.data.assign(static_cast<std::size_t>(nc * npts), Real(0.0));
    const Dim3 hi = amrex::ubound(region);
    const Dim3 hi1{hi.x + 1, hi.y + 1, hi.z + 1};   // Array4's end is one past the last cell
    Array4<Real> out(h.data.data(), amrex::lbound(region), hi1, nc);
    for (MFIter mfi(m_cc); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox() & region;   // every cell has one owner: m_cc has no ghosts
        if (!bx.ok()) { continue; }
#ifdef AMREX_USE_GPU
        FArrayBox host(bx, NCOMP, The_Pinned_Arena());
        host.copy<RunOn::Device>(m_cc[mfi], bx, 0, bx, 0, NCOMP);
        Gpu::streamSynchronize();
        const auto in = host.const_array();
#else
        const auto in = m_cc.const_array(mfi);
#endif
        LoopOnCpu(bx, [&] (int i, int j, int k) {
            out(i, j, k, 0) = in(i, j, k, U);
            out(i, j, k, 1) = in(i, j, k, V);
            out(i, j, k, 2) = in(i, j, k, Z);
            out(i, j, k, 3) = in(i, j, k, ZG);
            if (with_w) { out(i, j, k, 4) = in(i, j, k, W); }
        });
    }
    ParallelDescriptor::ReduceRealSum(h.data.data(), static_cast<int>(h.data.size()));

    FrozenColumnsView& c = h.view;
    c.a = Array4<const Real>(h.data.data(), amrex::lbound(region), hi1, nc);
    c.gilo = gilo; c.gihi = gihi;   // GUST's bounds: the index origin and the edge clamp
    c.gjlo = gjlo; c.gjhi = gjhi;
    c.kz0  = gdom.smallEnd(2); c.kz1  = gdom.bigEnd(2);
    c.gx0  = m_geom.ProbLo(0); c.gy0  = m_geom.ProbLo(1);
    c.gdx  = m_geom.CellSize(0); c.gdy = m_geom.CellSize(1);
}
