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

/**
 * The wind of one GUST column at a target height: the column sampler of
 * ERF_FireWindExtract.cpp (column_wind_at_height) on cell-centred u, v
 * (components 0, 1 of a) with the cell-centre heights in component 2 and
 * levels kz0..kz1. A target above the highest centre takes that centre's
 * wind; the caller has already refused it.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void frozen_column_wind (const Array4<const Real>& a, int ia, int ja, int kz0, int kz1,
                         Real z_target, Real z_ground, bool below_first_log, Real z0_log,
                         Real& u_out, Real& v_out) noexcept
{
    int k_lo;
    Real below_scale = Real(1.0);
    if (z_target <= a(ia, ja, kz0, 2)) {
        k_lo = kz0;
        if (below_first_log) {
            constexpr Real z_tiny = Real(1.0e-30);
            const Real z0 = std::abs(z0_log) + z_tiny;
            const Real z1 = amrex::max(a(ia, ja, kz0, 2) - z_ground, z0 * Real(1.000001));
            const Real zt = amrex::max(amrex::min(z_target - z_ground, z1), z0);
            below_scale = std::log(zt / z0) / std::log(z1 / z0);
        }
    } else if (z_target >= a(ia, ja, kz1, 2)) {
        k_lo = kz1 - 1;
    } else {
        int lo = kz0;
        int hi = kz1;
        while (hi - lo > 1) {
            const int mid = (lo + hi) / 2;
            if (a(ia, ja, mid, 2) <= z_target) { lo = mid; } else { hi = mid; }
        }
        k_lo = lo;
    }
    k_lo = amrex::max(kz0, amrex::min(k_lo, kz1 - 1));
    const int k_hi = k_lo + 1;

    const Real z_lo = a(ia, ja, k_lo, 2);
    const Real z_hi = a(ia, ja, k_hi, 2);
    Real alpha = Real(0.0);
    if (z_hi > z_lo) {
        alpha = amrex::max(Real(0.0), amrex::min(Real(1.0), (z_target - z_lo) / (z_hi - z_lo)));
    }
    u_out = below_scale * (a(ia, ja, k_lo, 0) + alpha * (a(ia, ja, k_hi, 0) - a(ia, ja, k_lo, 0)));
    v_out = below_scale * (a(ia, ja, k_lo, 1) + alpha * (a(ia, ja, k_hi, 1) - a(ia, ja, k_lo, 1)));
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

void FrozenWind::build_columns (const MultiFab& u, const MultiFab& v, const MultiFab& w,
                                const MultiFab& z, const MultiFab& z_nd)
{
    const Box& dom = m_geom.Domain();
    const int nz = dom.length(2);
    if (nz < 2) {
        frozen_abort(m_file, "the mesh needs at least two cells in z");
    }
    BoxArray ba(dom);
    ba.maxSize(IntVect(AMREX_D_DECL(32, 32, nz)));   // full-height columns, never split in z
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

                int  ia[4] = {0, 0, 0, 0};
                int  ja[4] = {0, 0, 0, 0};
                Real wt[4] = {Real(0.0), Real(0.0), Real(0.0), Real(0.0)};
                int  ncol  = 1;
                if (bilinear) {
                    const Real gx = (x - gx0) / gdx - Real(0.5);
                    const Real gy = (y - gy0) / gdy - Real(0.5);
                    const int i0 = static_cast<int>(std::floor(gx));
                    const int j0 = static_cast<int>(std::floor(gy));
                    const Real wx = gx - Real(i0);
                    const Real wy = gy - Real(j0);
                    for (int c = 0; c < 4; ++c) {
                        const int di = (c & 1);
                        const int dj = ((c >> 1) & 1);
                        ia[c] = amrex::max(gilo, amrex::min(gilo + i0 + di, gihi));
                        ja[c] = amrex::max(gjlo, amrex::min(gjlo + j0 + dj, gjhi));
                        wt[c] = (di ? wx : (Real(1.0) - wx)) * (dj ? wy : (Real(1.0) - wy));
                    }
                    ncol = 4;
                } else {
                    ia[0] = amrex::max(gilo, amrex::min(gilo + static_cast<int>(std::floor((x - gx0) / gdx)), gihi));
                    ja[0] = amrex::max(gjlo, amrex::min(gjlo + static_cast<int>(std::floor((y - gy0) / gdy)), gjhi));
                    wt[0] = Real(1.0);
                    ncol = 1;
                }

                Real u_sum = Real(0.0), v_sum = Real(0.0), z_sum = Real(0.0), g_sum = Real(0.0);
                Real excess = std::numeric_limits<Real>::lowest();
                for (int c = 0; c < ncol; ++c) {
                    const Real zg = a(ia[c], ja[c], kz0, 3);
                    const Real z_target = zg + z_at;
                    if (wt[c] > Real(0.0)) {   // a column with no weight cannot stop the run
                        excess = amrex::max(excess, z_target - a(ia[c], ja[c], kz1, 2));
                    }
                    Real uc, vc;
                    frozen_column_wind(a, ia[c], ja[c], kz0, kz1, z_target, zg, below_first_log, z0_log, uc, vc);
                    u_sum += wt[c] * uc;
                    v_sum += wt[c] * vc;
                    z_sum += wt[c] * z_target;
                    g_sum += wt[c] * zg;
                }
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
