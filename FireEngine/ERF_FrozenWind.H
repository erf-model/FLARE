#ifndef ERF_FROZEN_WIND_H
#define ERF_FROZEN_WIND_H

/**
 * @file ERF_FrozenWind.H
 * @brief A steady 3D wind field read from a GUST-CFD file and held frozen
 *        for the fire.
 *
 * With erf.fire.frozen_wind.file = <dir>, the fire layer reads the wind once
 * at start-up from a GUST-CFD (erf-model/gust-abl) plotfile or checkpoint and
 * never reads the atmosphere's wind again. The fire then spreads in that
 * steady wind while its own state (front, fuel, moisture, spotting) evolves in
 * time. ERF still advances the atmosphere; the fire only ignores its wind.
 *
 * What is read (the format is told apart by the first word of <dir>/Header):
 *  - plotfile ("HyperCLaw-V1.1"): the cell-centred u, v, w and the cell-centre
 *    height z, and the node heights <dir>/Level_0/z_phys_nd (GUST writes them
 *    next to every plotfile). Level 0 only.
 *  - checkpoint ("gust-checkpoint", or "ablsteady-checkpoint" before the
 *    rename; format version 1): the face velocities Level_0/u, v, w, averaged
 *    here to cell centres, and Level_0/z_phys_nd; the cell-centre height is the
 *    mean of the cell's eight nodes, as GUST's plotfile writes it. The
 *    checkpoint must hold reals of this build's size.
 *
 * The whole 3D field (u, v, w, z and the column ground) is kept, laid out in
 * full-height columns, for models that need the wind above the reference
 * height (cc()). The fire's reference wind is sampled from it once
 * (sample_reference_wind()).
 *
 * Coordinates: GUST's x, y are taken as the fire grid's x, y (the same
 * origin). Its horizontal extent must cover the fire region; the ground of
 * each GUST column is the mean of its four k = 0 nodes.
 */

#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

/**
 * @brief The GUST columns as the samplers read them: u, v, the cell-centre
 *        height z and the column ground (components 0-3 of a), and, when
 *        gathered with it, the vertical wind w (component 4), on GUST's index
 *        space, with its horizontal origin and cell sizes.
 *
 * Trivially copyable: captured by value in device kernels, or pointing to a
 * host copy (FrozenWindHost) for the host trajectories of ember spotting.
 */
struct FrozenColumnsView
{
    amrex::Array4<const amrex::Real> a;
    int gilo = 0, gihi = 0, gjlo = 0, gjhi = 0, kz0 = 0, kz1 = 0;
    amrex::Real gx0 = 0.0, gy0 = 0.0, gdx = 1.0, gdy = 1.0;
};

/**
 * The wind of one GUST column at a target height: the column sampler of
 * ERF_FireWindExtract.cpp (column_wind_at_height) on cell-centred u, v
 * (components 0, 1 of a) with the cell-centre heights in component 2 and
 * levels kz0..kz1. A target above the highest centre takes that centre's
 * wind; the caller decides whether that is allowed. If w_out is not null,
 * the vertical wind w (component 4, which a must hold) is sampled the same
 * way: interpolated in height, and below the first centre clamped or
 * log-scaled with u and v, so the wind keeps its direction there.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void frozen_column_wind (const amrex::Array4<const amrex::Real>& a, int ia, int ja, int kz0, int kz1,
                         amrex::Real z_target, amrex::Real z_ground, bool below_first_log,
                         amrex::Real z0_log, amrex::Real& u_out, amrex::Real& v_out,
                         amrex::Real* w_out = nullptr) noexcept
{
    using amrex::Real;
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
    if (w_out) {
        *w_out = below_scale * (a(ia, ja, k_lo, 4) + alpha * (a(ia, ja, k_hi, 4) - a(ia, ja, k_lo, 4)));
    }
}

/**
 * The frozen wind at a point (x, y) and a height h above the ground: each
 * GUST column is sampled at its own ground + h (frozen_column_wind), and the
 * four columns around (x, y) are blended bilinearly (indices clamped to the
 * GUST domain), or, with bilinear = false, the column holding (x, y) is used.
 *
 * @param[out] u, v     Wind [m/s]
 * @param[out] z_sum    Sampling height (column ground + h), weighted as the wind [m]
 * @param[out] g_sum    Column ground, weighted as the wind [m]
 * @param[out] excess   Largest (ground + h - highest centre) over the columns with weight [m];
 *                      > 0 when h lies above the top centre of such a column
 * @param[out] w        If not null: the vertical wind [m/s], blended as u and v
 *                      (the view must be gathered with w)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void frozen_point_wind (const FrozenColumnsView& c, amrex::Real x, amrex::Real y, amrex::Real h,
                        bool bilinear, bool below_first_log, amrex::Real z0_log,
                        amrex::Real& u, amrex::Real& v, amrex::Real& z_sum, amrex::Real& g_sum,
                        amrex::Real& excess, amrex::Real* w = nullptr) noexcept
{
    using amrex::Real;
    int  ia[4] = {0, 0, 0, 0};
    int  ja[4] = {0, 0, 0, 0};
    Real wt[4] = {Real(0.0), Real(0.0), Real(0.0), Real(0.0)};
    int  ncol  = 1;
    if (bilinear) {
        const Real gx = (x - c.gx0) / c.gdx - Real(0.5);
        const Real gy = (y - c.gy0) / c.gdy - Real(0.5);
        const int i0 = static_cast<int>(std::floor(gx));
        const int j0 = static_cast<int>(std::floor(gy));
        const Real wx = gx - Real(i0);
        const Real wy = gy - Real(j0);
        for (int n = 0; n < 4; ++n) {
            const int di = (n & 1);
            const int dj = ((n >> 1) & 1);
            ia[n] = amrex::max(c.gilo, amrex::min(c.gilo + i0 + di, c.gihi));
            ja[n] = amrex::max(c.gjlo, amrex::min(c.gjlo + j0 + dj, c.gjhi));
            wt[n] = (di ? wx : (Real(1.0) - wx)) * (dj ? wy : (Real(1.0) - wy));
        }
        ncol = 4;
    } else {
        ia[0] = amrex::max(c.gilo, amrex::min(c.gilo + static_cast<int>(std::floor((x - c.gx0) / c.gdx)), c.gihi));
        ja[0] = amrex::max(c.gjlo, amrex::min(c.gjlo + static_cast<int>(std::floor((y - c.gy0) / c.gdy)), c.gjhi));
        wt[0] = Real(1.0);
        ncol = 1;
    }

    u = v = z_sum = g_sum = Real(0.0);
    if (w) { *w = Real(0.0); }
    excess = std::numeric_limits<Real>::lowest();
    for (int n = 0; n < ncol; ++n) {
        const Real zg = c.a(ia[n], ja[n], c.kz0, 3);
        const Real z_target = zg + h;
        if (wt[n] > Real(0.0)) {   // a column with no weight cannot stop the run
            excess = amrex::max(excess, z_target - c.a(ia[n], ja[n], c.kz1, 2));
        }
        Real uc, vc, wc = Real(0.0);
        frozen_column_wind(c.a, ia[n], ja[n], c.kz0, c.kz1, z_target, zg, below_first_log, z0_log, uc, vc,
                           w ? &wc : nullptr);
        if (w) { *w += wt[n] * wc; }
        u     += wt[n] * uc;
        v     += wt[n] * vc;
        z_sum += wt[n] * z_target;
        g_sum += wt[n] * zg;
    }
}

/// A host copy of the GUST columns (u, v, z, ground, and w if asked for), the same on every rank
struct FrozenWindHost
{
    std::vector<amrex::Real> data;
    FrozenColumnsView        view;   ///< points into data
    FrozenWindHost () = default;
    FrozenWindHost (const FrozenWindHost&) = delete;             // view points into data
    FrozenWindHost& operator= (const FrozenWindHost&) = delete;
};

/**
 * @brief How ember spotting reads the frozen wind (erf.fire.spotting.drift_wind = frozen or frozen_uvw)
 *
 * The brand drifts on u, v of the 3D field at its own height above the
 * ground, sampled as the reference wind is (the same column interpolation,
 * horizontal blend and rule below the first cell centre). With vertical
 * (drift_wind = frozen_uvw) it also rises or falls with the field's w: its
 * vertical speed is w - v_t instead of -v_t.
 */
struct FrozenDrift
{
    FrozenColumnsView cols;
    bool              bilinear        = true;    ///< erf.fire.wind_interp = 1
    bool              below_first_log = false;   ///< erf.fire.wind_below_first_cell = log
    amrex::Real       z0_log          = 0.1;     ///< erf.fire.wind_sample_z0 [m]
    bool              vertical        = false;   ///< drift_wind = frozen_uvw: cols holds w (component 4)
};

// FrozenColumnsView and FrozenDrift are passed by value into device code
static_assert(std::is_trivially_copyable_v<FrozenColumnsView>, "FrozenColumnsView must be trivially copyable");
static_assert(std::is_trivially_copyable_v<FrozenDrift>, "FrozenDrift must be trivially copyable");

class FrozenWind
{
public:
    /// Components of cc()
    enum Comp : int { U = 0, V = 1, W = 2, Z = 3, ZG = 4, NCOMP = 5 };

    /**
     * @brief Read a GUST plotfile or checkpoint directory
     *
     * Aborts, naming erf.fire.frozen_wind.file, when the directory is not one
     * of the two formats, a field is missing, the mesh is not 3D, or a value
     * read is not finite.
     */
    void read (const std::string& dir);

    [[nodiscard]] bool defined () const noexcept { return m_cc.ok(); }

    /// The file read and its format ("plotfile" or "checkpoint")
    [[nodiscard]] const std::string& file () const noexcept { return m_file; }
    [[nodiscard]] const std::string& format () const noexcept { return m_format; }

    /// GUST's geometry (index space and physical extent of its mesh)
    [[nodiscard]] const amrex::Geometry& geom () const noexcept { return m_geom; }

    /**
     * @brief The 3D field in full-height columns (every box spans the whole
     *        domain in z), no ghost cells.
     *
     * Components: U, V, W [m/s] at cell centres, Z the cell-centre height [m],
     * ZG the ground height of the cell's column [m] (the same at every k).
     */
    [[nodiscard]] const amrex::MultiFab& cc () const noexcept { return m_cc; }

    /**
     * @brief The reference wind on the fire grid, sampled from the frozen field
     *
     * The same algorithm as the atmospheric sampler fill_fire_wind_from_interpolation
     * (ERF_FireWindExtract.H): each GUST column is sampled at the target height
     * above its own ground, linearly in height; a target below the lowest
     * centre takes that centre's wind, as it is or scaled down a neutral log
     * profile (below_first_log); with z_sample > 0 the column is sampled at
     * z_sample and brought down to the target height with the log profile.
     * With wind_interp = 1 the four columns around the fire cell centre are
     * blended bilinearly (a fire cell within half a GUST cell of the GUST
     * domain edge takes the edge columns); with 0 the column that holds the
     * fire cell centre is used.
     *
     * Unlike the atmospheric sampler, a target above the highest cell centre
     * of a column aborts instead of taking the highest centre's wind.
     *
     * @param[out] wind_ref   u, v [m/s] on the fire grid (2 components)
     * @param[out] extract_z  Height [m] the wind was sampled at (GUST ground + height), weighted as the wind
     * @param[out] ground     GUST ground [m] at the fire cell centre, weighted as the wind (1 component), or nullptr
     * @param[in]  fire_geom  Geometry of the fire grid
     * @param[in]  z_ref      Target height above ground [m]
     * @param[in]  z_ref_cell Per-cell target height on the fire grid, or nullptr (z_ref everywhere)
     * @param[in]  wind_interp 1 bilinear, 0 nearest column
     * @param[in]  z_sample   Sampling height [m] of the two-height resampling, 0 = off
     * @param[in]  z0_log     Roughness [m] of the log profiles
     * @param[in]  below_first_log Log profile below the lowest centre (true) or its wind as it is
     */
    void sample_reference_wind (amrex::MultiFab& wind_ref,
                                amrex::MultiFab& extract_z,
                                amrex::MultiFab* ground,
                                const amrex::Geometry& fire_geom,
                                amrex::Real z_ref,
                                const amrex::MultiFab* z_ref_cell,
                                int wind_interp,
                                amrex::Real z_sample,
                                amrex::Real z0_log,
                                bool below_first_log) const;

    /**
     * @brief Copy the columns (u, v, z, ground, and w with with_w) under the
     *        fire region to the host of every rank, for the host trajectories
     *        of ember spotting
     *
     * The columns a point of the fire region reads, with a cell of margin (as
     * sample_reference_wind takes them), clipped to GUST's domain: a brand
     * stops when it leaves the fire region, so it reads no other column. The
     * view keeps GUST's own index bounds, so the edge columns are clamped as
     * in sample_reference_wind. 4 n_x n_y n_z reals of that region on every
     * rank (256 x 256 x 64 columns take 128 MiB in double precision), 5 with w.
     *
     * @param[out] h          The host copy and its view (h.view points into h.data)
     * @param[in]  fire_geom  Geometry of the fire grid
     * @param[in]  with_w     Also copy the vertical wind w, as component 4
     *                        (erf.fire.spotting.drift_wind = frozen_uvw)
     */
    void gather_columns (FrozenWindHost& h, const amrex::Geometry& fire_geom, bool with_w = false) const;

    /**
     * @brief Define the field from cell-centred u, v, w, z and the node
     *        heights on GUST's geometry, as read() does from a file
     *
     * For tests: max_box sets the x, y size of the column boxes (read() uses 32).
     */
    void define_from_fields (const amrex::Geometry& geom, const std::string& label,
                             const amrex::MultiFab& u, const amrex::MultiFab& v,
                             const amrex::MultiFab& w, const amrex::MultiFab& z,
                             const amrex::MultiFab& z_nd, int max_box = 32);

    // Public because they launch device lambdas (CUDA forbids extended
    // lambdas in private or protected member functions); called by read().
    void read_plotfile (const std::string& dir);
    void read_checkpoint (const std::string& dir);
    /// Lay the cell-centred u, v, w, z and the node heights out in full-height
    /// columns (m_cc) and fill the column ground
    void build_columns (const amrex::MultiFab& u, const amrex::MultiFab& v,
                        const amrex::MultiFab& w, const amrex::MultiFab& z,
                        const amrex::MultiFab& z_nd, int max_box = 32);

private:
    std::string     m_file;
    std::string     m_format;
    amrex::Geometry m_geom;
    amrex::MultiFab m_cc;
};

#endif
