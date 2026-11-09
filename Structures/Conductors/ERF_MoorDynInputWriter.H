#ifndef ERF_MOORDYN_INPUT_WRITER_H
#define ERF_MOORDYN_INPUT_WRITER_H

/** The MoorDyn-C input file of one conductor line, generated from the ERF inputs, and the MoorDyn -> ERF frame conversion. */

#include <array>
#include <string>
#include <vector>

#include <AMReX_REAL.H>

#include "ERF_ConductorInputs.H"

namespace erf_conductors {

/**
 * The MoorDyn-C input file of one conductor line, generated from the ERF inputs.
 *
 * With N = s.num_spans(), points 1 .. N+1 are the attachment points in order (end_a, the towers,
 * end_b). end_a and end_b are Fixed; the towers are Fixed, or Coupled when the line's towers move
 * (ConductorInputs::towers_move), and ERF then drives them. With insulator strings, points
 * N+2 .. 2N are the Free points at the lower end of the string at each tower, where the conductor
 * hangs. MoorDyn lines 1 .. N are the spans, from end_a to end_b; with insulator strings, MoorDyn
 * lines N+1 .. 2N-1 are the strings, each from its tower down to its free point.
 *
 * MoorDyn applies its fluid loads to nodes below its free surface at z = 0 and puts a flat
 * bottom at -WtrDpth, so the line is written in MoorDyn's frame, every height lowered by
 * surface_offset (z_md = z_erf - surface_offset); the bottom is placed another surface_offset
 * lower, at ERF z = -surface_offset, where nothing reaches it. The fluid density is the air
 * density, the fluid kinematics are external (WaveKin = 1, the wind ERF samples at the nodes),
 * the initial shape comes from MoorDyn's stationary solver, and MoorDyn's own outputs are off
 * (ERF writes the diagnostics).
 *
 *
 * With conductors.ground_contact the bottom is the terrain instead: the options name the seafloor
 * file (seafloor_file()) that the caller writes before MoorDyn reads the input, with
 * ConductorInputs::ground_stiffness as kbot and the conductor's critical ground_damping() as cbot.
 *
 * @param s the line; its attachment points already at their absolute heights (the terrain added)
 * @param in the shared settings (air density, surface offset, moordyn_dt, and moordyn_cfl, MoorDyn's
 *           Courant (CFL) factor) and the tower types, which decide whether the towers are Coupled points
 * @param gravity ERF's gravity (m/s^2)
 */
std::string moordyn_input_text (const LineInputs& s, const ConductorInputs& in, amrex::Real gravity,
                                const std::string& seafloor = std::string());

/**
 * Write moordyn_input_text() to a file (creating its directory); aborts when the file cannot be written.
 * With conductors.ground_contact the input names seafloor_file(fname).
 */
void write_moordyn_input (const std::string& fname, const LineInputs& s, const ConductorInputs& in, amrex::Real gravity);

/**
 * The seafloor file that goes with the MoorDyn input fname: its ".moordyn.txt" ending replaced by
 * ".seafloor.txt" (".seafloor.txt" appended to any other name). MoorDyn opens it from the run's
 * directory, as it opens the input.
 */
std::string seafloor_file (const std::string& fname);

/**
 * The text of a MoorDyn-C 3-D seafloor file: "nx ny", the nx x ticks, the ny y ticks, then one
 * "x y z" row per grid point, x the slower index. z is the ground's height in MoorDyn's frame,
 * z_erf - surface_offset; MoorDyn interpolates bilinearly between the points and holds the edge
 * values beyond the grid.
 *
 * @param x the x ticks (m), increasing
 * @param y the y ticks (m), increasing
 * @param z_erf the ground's height in ERF's frame (m) at (x[i], y[j]), at z_erf[i * y.size() + j]
 * @param surface_offset conductors.surface_offset (m)
 */
/**
 * The node indices of a seafloor grid from nlo to nhi (inclusive) in x and y: every node, or every
 * stride-th with the smallest stride that keeps the grid at max_points or fewer, the last node always
 * included (so the grid covers the whole range).
 */
void seafloor_ticks (const std::array<int,2>& nlo, const std::array<int,2>& nhi, long max_points,
                     std::vector<int>& ii, std::vector<int>& jj, int& stride);

/**
 * Replace the NaN heights in h (nodes no grid holds) by the lowest of the others, set in lowest.
 * Returns how many were replaced; when all are NaN, h is left as it is and h.size() returned.
 */
std::size_t fill_missing_heights (std::vector<double>& h, double& lowest);

std::string seafloor_text (const std::vector<double>& x, const std::vector<double>& y,
                           const std::vector<double>& z_erf, double surface_offset);

/** A position (m) in ERF's frame from MoorDyn's: z_erf = z_md + surface_offset. */
inline std::array<amrex::Real,3> to_erf_frame (const std::array<double,3>& p, amrex::Real surface_offset)
{
    // the offset added in double: in a single-precision build MoorDyn's z, near -surface_offset, holds only centimetres
    return {{static_cast<amrex::Real>(p[0]), static_cast<amrex::Real>(p[1]),
             static_cast<amrex::Real>(p[2] + static_cast<double>(surface_offset))}};
}

} // namespace erf_conductors

#endif
