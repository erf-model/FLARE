#include "ERF_MoorDynInputWriter.H"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

#include <AMReX.H>
#include <AMReX_Utility.H>

using namespace amrex;

namespace erf_conductors {

std::string moordyn_input_text (const LineInputs& s, const ConductorInputs& in, Real gravity, const std::string& seafloor)
{
    const int N = s.num_spans();
    const bool strings = s.has_insulators();
    const std::string ins_type = s.name + "_insulator";
    // the point a span ends on at attachment k: the fixed point itself, or the free point under the
    // tower's string
    auto hang = [&] (int k) { return (strings && k > 0 && k < N) ? N + 1 + k : k + 1; };
    std::ostringstream out;
    // material constants to twelve digits, or as many as a Real holds: a single-precision build writes
    // 0.0281, not 0.0280999993
    out << std::setprecision(std::min(12, std::numeric_limits<Real>::digits10));
    // positions and lengths in double to twelve digits, so a single-precision build keeps every bit of
    // them (at 30 km six digits would move a point by 5 cm, and 300.0104 m would become 300.01)
    auto exact = [] (double v) {
        std::ostringstream o;
        o << std::setprecision(12) << v;
        return o.str();
    };
    const double offset = static_cast<double>(in.surface_offset);
    out << "MoorDyn-C input written by ERF for conductor line " << s.name << " (conductors." << s.name << ".*)\n"
        << "----------------------- LINE TYPES ------------------------------------------\n"
        << "TypeName   Diam     Mass/m     EA         BA/-zeta    EI         Cd     Ca     CdAx    CaAx\n"
        << "(name)     (m)      (kg/m)     (N)        (N-s/-)     (N-m^2)    (-)    (-)    (-)     (-)\n"
        << s.name << "   " << s.diameter << "   " << s.mass_per_length << "   " << s.axial_stiffness
        << "   " << -s.damping_ratio << "   0   " << s.drag_coefficient << "   0.0   0.0   0.0\n";
    if (strings) {
        out << ins_type << "   " << s.insulator_diameter << "   " << s.insulator_mass / s.insulator_length << "   "
            << LineInputs::insulator_axial_stiffness << "   " << -s.damping_ratio << "   0   "
            << LineInputs::insulator_drag_coefficient << "   0.0   0.0   0.0\n";
    }
    out << "---------------------- POINT PROPERTIES --------------------------------\n"
        << "ID    Type      X       Y       Z       Mass   Volume  CdA    Ca\n"
        << "(#)   (-)       (m)     (m)     (m)     (kg)   (m^3)   (m^2)  (-)\n";
    // the towers' cross-arms are coupled points when the towers move: ERF drives them and takes their pull
    const bool moving = in.towers_move(s);
    for (int k = 0; k <= N; ++k) {
        const auto& p = s.point(k);
        const char* type = (moving && k > 0 && k < N) ? "Coupled" : "Fixed  ";
        out << k + 1 << "     " << type << "   " << exact(p[0]) << "   " << exact(p[1]) << "   "
            << exact(static_cast<double>(p[2]) - offset) << "   0   0   0   0\n";
    }
    if (strings) {
        for (int k = 1; k < N; ++k) {
            const auto& p = s.point(k);
            out << N + 1 + k << "     Free      " << exact(p[0]) << "   " << exact(p[1]) << "   "
                << exact(static_cast<double>(p[2]) - static_cast<double>(s.insulator_length) - offset) << "   0   0   0   0\n";
        }
    }
    out << "---------------------- LINES ----------------------------------------\n"
        << "ID   LineType   AttachA  AttachB  UnstrLen  NumSegs  LineOutputs\n"
        << "(#)   (name)     (#)      (#)       (m)       (-)     (-)\n";
    for (int k = 0; k < N; ++k) {
        out << k + 1 << "     " << s.name << "      " << hang(k) << "        " << hang(k + 1) << "         "
            << exact(s.lengths[static_cast<std::size_t>(k)]) << "   " << s.segments << "   -\n";
    }
    if (strings) {
        for (int k = 1; k < N; ++k) {
            out << N + k << "     " << ins_type << "      " << k + 1 << "        " << N + 1 + k << "         "
                << s.insulator_length << "   " << LineInputs::insulator_segments << "   -\n";
        }
    }
    out << "---------------------- OPTIONS -----------------------------------------\n"
        << "0             writeLog      ERF writes the diagnostics\n";
    if (in.moordyn_dt > 0.0) {
        out << in.moordyn_dt << "   dtM           upper bound on MoorDyn's internal step (s, conductors.moordyn_dt)\n";
    }
    out << in.moordyn_cfl << "   CFL           Courant factor that sets MoorDyn's internal step (conductors.moordyn_cfl)\n"
        << gravity << "   g             gravity (m/s^2)\n"
        << in.air_density << "   WtrDnsty      the fluid is air (kg/m^3)\n"
        << exact(2.0 * offset) << "   WtrDpth       flat bottom, far below the ground (m)\n";
    if (in.ground_contact) {
        // MoorDyn splits an option line at blanks: the file name must be one word
        if (seafloor.empty() || seafloor.find_first_of(" \t") != std::string::npos) {
            Abort("conductors." + s.name + ": with ground_contact the seafloor file name must be given and hold no blanks (it is '" +
                  seafloor + "'; conductors.diagnostics_dir names its directory)");
        }
        out << seafloor << "   SeafloorFile  the terrain under the line (conductors.ground_contact)\n"
            << ConductorInputs::ground_stiffness << "   kbot          the ground's contact stiffness (Pa/m)\n"
            << ConductorInputs::ground_damping(s.mass_per_length, s.diameter)
            << "   cbot          the ground's contact damping, critical for the conductor (Pa-s/m)\n";
    }
    out
        << "1             WaveKin       the fluid kinematics come through the API\n"
        << "0             ICgenDynamic  stationary initial-condition solver\n"
        << "1             disableOutput\n"
        << "1             disableOutTime\n"
        << "------------------------- need this line --------------------------------------\n";
    return out.str();
}

void write_moordyn_input (const std::string& fname, const LineInputs& s, const ConductorInputs& in, Real gravity)
{
    const auto slash = fname.rfind('/');
    if (slash != std::string::npos) { UtilCreateDirectory(fname.substr(0, slash), 0755); }
    std::ofstream out(fname, std::ios::trunc);
    if (!out) { Abort("conductors." + s.name + ": cannot write the MoorDyn input file '" + fname + "'"); }
    out << moordyn_input_text(s, in, gravity, in.ground_contact ? seafloor_file(fname) : std::string());
}

std::string seafloor_file (const std::string& fname)
{
    const std::string ending = ".moordyn.txt";
    if (fname.size() > ending.size() && fname.compare(fname.size() - ending.size(), ending.size(), ending) == 0) {
        return fname.substr(0, fname.size() - ending.size()) + ".seafloor.txt";
    }
    return fname + ".seafloor.txt";
}

void seafloor_ticks (const std::array<int,2>& nlo, const std::array<int,2>& nhi, long max_points,
                     std::vector<int>& ii, std::vector<int>& jj, int& stride)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nhi[0] > nlo[0] && nhi[1] > nlo[1] && max_points >= 4,
                                     "seafloor_ticks: needs at least two nodes each way and four points");
    stride = 1;
    auto ticks = [&] (int d) { return static_cast<long>((nhi[d] - nlo[d] + stride - 1) / stride) + 1; };
    while (ticks(0) * ticks(1) > max_points) { ++stride; }
    ii.clear(); jj.clear();
    for (int i = nlo[0]; i < nhi[0]; i += stride) { ii.push_back(i); }
    ii.push_back(nhi[0]);
    for (int j = nlo[1]; j < nhi[1]; j += stride) { jj.push_back(j); }
    jj.push_back(nhi[1]);
}

std::size_t fill_missing_heights (std::vector<double>& h, double& lowest)
{
    std::size_t missing = 0;
    lowest = std::numeric_limits<double>::infinity();
    for (const double v : h) { if (std::isnan(v)) { ++missing; } else { lowest = std::min(lowest, v); } }
    if (missing == h.size()) { return missing; }
    for (double& v : h) { if (std::isnan(v)) { v = lowest; } }
    return missing;
}

std::string seafloor_text (const std::vector<double>& x, const std::vector<double>& y,
                           const std::vector<double>& z_erf, double surface_offset)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(x.size() >= 2 && y.size() >= 2 && z_erf.size() == x.size() * y.size(),
                                     "seafloor_text: needs at least 2 x 2 ticks and one height per grid point");
    std::ostringstream out;
    // MoorDyn reads the numbers as float: twelve digits keep every bit it can hold
    out << std::setprecision(12) << x.size() << " " << y.size() << "\n";
    for (std::size_t i = 0; i < x.size(); ++i) { out << (i ? " " : "") << x[i]; }
    out << "\n";
    for (std::size_t j = 0; j < y.size(); ++j) { out << (j ? " " : "") << y[j]; }
    out << "\n";
    for (std::size_t i = 0; i < x.size(); ++i) {
        for (std::size_t j = 0; j < y.size(); ++j) {
            out << x[i] << " " << y[j] << " " << z_erf[i * y.size() + j] - surface_offset << "\n";
        }
    }
    return out.str();
}

} // namespace erf_conductors
