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

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <string>

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

    // Public because they launch device lambdas (CUDA forbids extended
    // lambdas in private or protected member functions); called by read().
    void read_plotfile (const std::string& dir);
    void read_checkpoint (const std::string& dir);
    /// Lay the cell-centred u, v, w, z and the node heights out in full-height
    /// columns (m_cc) and fill the column ground
    void build_columns (const amrex::MultiFab& u, const amrex::MultiFab& v,
                        const amrex::MultiFab& w, const amrex::MultiFab& z,
                        const amrex::MultiFab& z_nd);

private:
    std::string     m_file;
    std::string     m_format;
    amrex::Geometry m_geom;
    amrex::MultiFab m_cc;
};

#endif
