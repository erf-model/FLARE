// Stand-in for libmoordyn, implementing the subset of the MoorDyn-C v2 C API declared in
// Stub/moordyn/MoorDyn2.h. It lets ERF's MoorDyn wrapper and the body models on it be built and
// tested where MoorDyn is not installed (CI). It reads the same input file (LINE TYPES, POINTS or POINT PROPERTIES, LINES,
// OPTIONS) and has the geometry, the numbering and the order of the kinematics points of the real
// API, but no line dynamics. "Line" in this file is a MoorDyn line (for overhead conductor lines,
// one span or one insulator string).
//
// What it approximates: every MoorDyn line hangs between its two end points as an elastic
// parabola, as long as its unstretched length stretched by its tension, slack or taut. It swings
// about the chord to the quasi-static blowout angle atan(q / w) set by the fluid velocity it is
// given (q the drag per unit length, w the weight minus buoyancy per unit length, both N/m),
// relaxing towards that angle with a one-second time constant so that the state depends on time,
// and carries the parabola's tension (H along the chord, w_e c^2/(8 sag)). Fixed, coupled and free
// points are accepted; a coupled point moves as MoorDyn moves it, from the position it is given at
// the start of a step at the velocity it is given over the step. A free point hangs from the
// shortest line that ties it to a held point, fixed or coupled (an insulator string from its
// tower), plumb below that point and swung across by the wind across the other line attached to
// it, as that line was last placed; the hanging line runs straight from the held point to it. The
// fluid loads act on nodes below z = 0, as in MoorDyn (the drag on each node is
// 0.5 rho Cd D |u_n| u_n over its share of the chord), and the kinematics points follow MoorDyn-C
// 2.7.1's order: the line nodes, the points, then one entry at the origin (the ground body).
//
// Differences from MoorDyn-C 2.7.1:
//  - no dynamics: each step places every line quasi-statically; dtM is stored but never used for
//    sub-steps, and the CFL option is ignored;
//  - MoorDyn_Init and MoorDyn_Init_NoIC are identical: both place the static shape (no
//    initial-condition solve), and the coupled velocities xd are ignored;
//  - the swing angle relaxes to atan(q / w) with a fixed 1 s time constant; q comes from the normal
//    fluid velocity averaged over the line, and only when the chord's midpoint is below z = 0; w is
//    the weight minus buoyancy in the fluid of density rho (WtrDnsty, default 1025 kg/m^3);
//  - node tensions are the parabola's, H (ec + slope es), with no internal damping, bending
//    stiffness; node velocities are finite differences of successive placements;
//  - seabed contact only with a SeafloorFile (MoorDyn-C's 3-D seafloor, same format, interpolated
//    bilinearly and held at its edge values beyond the grid as MoorDyn-C does): a span's node
//    placed below the seafloor is lifted onto it, its tension unchanged; there is no contact
//    stiffness, damping or friction (kbot and cbot are read and unused), no contact for free points
//    (MoorDyn-C's free points rest through their lines' end nodes), no contact on the flat bottom at
//    -WtrDpth, and the file is read as a stream (MoorDyn-C needs nx entries on line 2, ny on line 3 and
//    three per row);
//  - MoorDyn_GetLineNodeForce is zero at inner nodes and plus or minus the end tension at the end
//    nodes (no end-node weight, drag or inertia);
//  - MoorDyn_GetLineMaxTen includes the last node (MoorDyn-C stops at node N-1);
//  - point forces have no point mass, volume or CdA terms;
//  - MoorDyn_ExternalWaveKinInit returns 0 unless WaveKin = 1 (MoorDyn-C returns the full count
//    whatever WaveKin is); MoorDyn_ExternalWaveKinSet ignores Ud and t;
//  - MoorDyn_SetLineUnstretchedLength changes the length the next placement uses; MoorDyn_Save and
//    MoorDyn_Serialize store the length with each line (MoorDyn-C stores the segment lengths but not the
//    unstretched length, so a caller that changes it sets it again on a restart);
//  - MoorDyn_Save stores the swing angles, the node positions, the last fluid velocity and the free
//    and coupled points in a text format; MoorDyn_Load restarts the node velocities from rest;
//    MoorDyn_Serialize stores the whole stub state as doubles in its own format, the fluid velocity
//    excluded;
//  - MoorDyn_BreakLine detaches a line's end from a point onto a new free point, as MoorDyn-C does, and
//    the kinematics points grow by one at the next MoorDyn_ExternalWaveKinInit. What the detached line
//    does is a quasi-static stand-in for MoorDyn-C's fall: a span with one end detached hangs plumb
//    from its held end down to the seafloor (its whole length without a SeafloorFile) and lies on the
//    seafloor beyond, towards where its other end was, with the tension of the hanging length's
//    weight below each node and none on the seafloor, and no drag; a span with both ends detached lies
//    on the seafloor where it was; a string detached from its held point hangs below its free point
//    with no tension, and that free point moves onto the elastic parabola the two spans it joins make
//    together between their far ends, where the first span's length ends (it does not slide along the
//    line, so two spans of unequal length carry unequal tensions there). A saved or serialized state fits only a system with the same
//    breaks (the caller repeats them before MoorDyn_Load or MoorDyn_Deserialize, as with MoorDyn-C);
//  - input: v2 or v1 LINE TYPES columns; point types fixed/anchor, coupled/vessel/cpld and
//    free/connect; line ends on rods or bodies (R1, B1) are not supported, and a malformed
//    number makes MoorDyn_Create fail.

#include "moordyn/MoorDyn2.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <deque>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr double pi = 3.14159265358979323846;

struct LineType {
    std::string name;
    double diam = 0.0;   // m
    double mass = 0.0;   // kg/m
    double EA = 0.0;     // N
    double Cd = 1.0;     // - (normal drag coefficient)
};

struct Point {
    int type = 1;                         // -1 coupled, 0 free, 1 fixed
    double pos[3] = {0.0, 0.0, 0.0};      // position (m, MoorDyn's frame)
    double force[3] = {0.0, 0.0, 0.0};    // net force of the lines' ends on the point (N)
    int hang_line = -1;   // a free point: the MoorDyn line it hangs from (0-based index into lines)
    int swing_line = -1;  // a free point: the MoorDyn line whose swing direction it follows (0-based)
    bool loose = false;   // made by MoorDyn_BreakLine: it follows the end of the line detached onto it
};

struct Line {
    int type_index = 0;             // 0-based index into types
    int attachA = 0, attachB = 0;   // MoorDyn point ids of the line's ends, from 1
    double length = 0.0;            // unstretched length (m)
    unsigned nseg = 1;              // segments; nodes = nseg + 1
    double phi = 0.0;               // swing angle about the chord (rad), 0 = hanging straight down
    double sag = 0.0;               // sag at mid-span (m)
    double H = 0.0;                 // tension component along the chord (N): the parabola's H in the chord frame
    double w_eff = 0.0;             // effective load per unit length, sqrt(w^2 + q^2) (N/m)
    std::vector<double> pos, vel, ten, drag;   // 3*(nseg+1) each: m, m/s, N, N
    bool hanging = false;            // ties a free point to a held point: placed straight between them
    double es[3] = {0.0, 0.0, -1.0}; // unit direction of the sag, as last placed
    double across[3] = {0.0, 1.0, 0.0};  // unit horizontal normal of its chord, as last placed
    double phi_across = 0.0;             // the angle the across-line wind swings a string it hangs from (rad)
    // a span with a detached end: the horizontal unit direction from its held end to where the other was. It is in
    // the saved file (MoorDyn_Save), which a new system loads; not in the serialized state, which returns to this system
    double drop[2] = {1.0, 0.0};
    double drop_load = 0.0;              // a span with a detached end: the weight of its hanging length (N)
};

struct StubSystem;
struct LineHandle { StubSystem* s; unsigned index; };    // index 0-based
struct PointHandle { StubSystem* s; unsigned index; };   // index 0-based

// MoorDyn-C's 3-D seafloor: the depth z (m, MoorDyn's frame) at the grid points (x[i], y[j]), z[i * ny + j]
struct Seafloor {
    std::vector<double> x, y, z;
    bool active () const { return !z.empty(); }
    // the index i - 1 below x (clamped) and the weight f of x[i], as MoorDyn-C's interp_factor
    static void factor (const std::vector<double>& t, double v, std::size_t& i, double& f)
    {
        if (v <= t.front()) { i = 0; f = 0.0; return; }
        if (v >= t.back()) { i = t.size() - 2; f = 1.0; return; }
        i = static_cast<std::size_t>(std::upper_bound(t.begin(), t.end(), v) - t.begin()) - 1;
        f = (v - t[i]) / (t[i+1] - t[i]);
    }
    double depth_at (double px, double py) const
    {
        std::size_t i = 0, j = 0;
        double fx = 0.0, fy = 0.0;
        factor(x, px, i, fx);
        factor(y, py, j, fy);
        const std::size_t ny = y.size();
        return (1.0 - fx) * ((1.0 - fy) * z[i*ny + j] + fy * z[i*ny + j + 1]) + fx * ((1.0 - fy) * z[(i+1)*ny + j] + fy * z[(i+1)*ny + j + 1]);
    }
};

bool read_seafloor (const std::string& fname, Seafloor& sf)
{
    std::ifstream in(fname);
    if (!in) { std::fprintf(stderr, "MoorDyn stub: cannot open the seafloor file '%s'\n", fname.c_str()); return false; }
    std::size_t nx = 0, ny = 0;
    if (!(in >> nx >> ny) || nx < 2 || ny < 2) { std::fprintf(stderr, "MoorDyn stub: '%s' needs nx ny >= 2 first\n", fname.c_str()); return false; }
    sf.x.resize(nx); sf.y.resize(ny); sf.z.assign(nx * ny, 0.0);
    for (auto* t : {&sf.x, &sf.y}) {
        for (std::size_t k = 0; k < t->size(); ++k) {
            if (!(in >> (*t)[k]) || !std::isfinite((*t)[k]) || (k > 0 && !((*t)[k] > (*t)[k-1]))) {
                std::fprintf(stderr, "MoorDyn stub: the ticks in '%s' must be finite and increasing\n", fname.c_str());
                return false;
            }
        }
    }
    for (std::size_t k = 0; k < nx * ny; ++k) {
        double px = 0.0, py = 0.0, pz = 0.0;
        if (!(in >> px >> py >> pz) || !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz)) {
            std::fprintf(stderr, "MoorDyn stub: '%s' holds fewer than nx ny finite points\n", fname.c_str());
            return false;
        }
        // the grid point the row is at: as MoorDyn-C, the tick at or below it
        std::size_t i = 0, j = 0;
        double f = 0.0;
        Seafloor::factor(sf.x, px, i, f); if (f >= 1.0 - 1.0e-9) { ++i; }
        Seafloor::factor(sf.y, py, j, f); if (f >= 1.0 - 1.0e-9) { ++j; }
        sf.z[i * ny + j] = pz;
    }
    return true;
}

struct StubSystem {
    Seafloor seafloor;
    std::vector<LineType> types;
    std::vector<Point> points;
    std::vector<Line> lines;
    double dtM = 0.001;      // s
    double g = 9.80665;      // m/s^2
    double rho = 1025.0;     // kg/m^3: MoorDyn's water default; ERF writes the air density
    double depth = 0.0;      // m, read and unused
    int wave_kin = 0;
    bool initialised = false;
    unsigned nkin = 0;
    std::vector<double> U;          // 3*nkin, the fluid velocity last set (m/s)
    double t = 0.0;                 // the clock (s)
    // one handle per MoorDyn line and point, made once the input is read: MoorDyn_GetLine and
    // MoorDyn_GetPoint return these, so repeated calls allocate nothing
    std::vector<LineHandle> line_handles;
    std::deque<PointHandle> point_handles;   // a deque: MoorDyn_BreakLine adds a point and the handles given out stay valid
};

std::string lower (std::string s)
{
    for (auto& c : s) { c = static_cast<char>(std::tolower(static_cast<unsigned char>(c))); }
    return s;
}

std::vector<std::string> tokens (const std::string& line)
{
    std::istringstream in(line);
    std::vector<std::string> out;
    std::string tok;
    while (in >> tok) { out.push_back(tok); }
    return out;
}

bool is_header (const std::string& line) { return line.find("---") != std::string::npos; }

bool parse (const std::string& fname, StubSystem& s)
{
    std::ifstream in(fname);
    if (!in) { std::fprintf(stderr, "MoorDyn stub: cannot open '%s'\n", fname.c_str()); return false; }
    std::string section, line;
    int skip = 0;
    while (std::getline(in, line)) {
        if (is_header(line)) {
            const std::string u = lower(line);
            if (u.find("line type") != std::string::npos) { section = "types"; skip = 2; }
            else if (u.find("point") != std::string::npos || u.find("connection") != std::string::npos) { section = "points"; skip = 2; }
            else if (u.find("line") != std::string::npos) { section = "lines"; skip = 2; }
            else if (u.find("option") != std::string::npos) { section = "options"; skip = 0; }
            else { section = "other"; skip = 0; }
            continue;
        }
        if (skip > 0) { --skip; continue; }
        const auto tok = tokens(line);
        if (tok.empty()) { continue; }
        if (section == "types") {
            if (tok.size() < 4) { std::fprintf(stderr, "MoorDyn stub: bad line type row '%s'\n", line.c_str()); return false; }
            LineType t;
            t.name = tok[0]; t.diam = std::stod(tok[1]); t.mass = std::stod(tok[2]); t.EA = std::stod(tok[3]);
            if (tok.size() >= 10) { t.Cd = std::stod(tok[6]); }        // v2 columns: EI Cd Ca CdAx CaAx after BA
            else if (tok.size() >= 8) { t.Cd = std::stod(tok[7]); }   // v1 columns: Can Cat Cdn Cdt after BA
            s.types.push_back(t);
        } else if (section == "points") {
            if (tok.size() < 5) { std::fprintf(stderr, "MoorDyn stub: bad point row '%s'\n", line.c_str()); return false; }
            Point p;
            const std::string ty = lower(tok[1]);
            if (ty.rfind("fix", 0) == 0 || ty.rfind("anch", 0) == 0) { p.type = 1; }
            else if (ty.rfind("coupled", 0) == 0 || ty.rfind("vessel", 0) == 0 || ty.rfind("cpld", 0) == 0) { p.type = -1; }
            else if (ty.rfind("free", 0) == 0 || ty.rfind("conn", 0) == 0) { p.type = 0; }
            else { std::fprintf(stderr, "MoorDyn stub: point type '%s' is not supported (fixed, coupled or free)\n", tok[1].c_str()); return false; }
            for (int d = 0; d < 3; ++d) { p.pos[d] = std::stod(tok[2 + d]); }
            s.points.push_back(p);
        } else if (section == "lines") {
            if (tok.size() < 6) { std::fprintf(stderr, "MoorDyn stub: bad line row '%s'\n", line.c_str()); return false; }
            Line l;
            const std::string ty = tok[1];
            auto it = std::find_if(s.types.begin(), s.types.end(), [&](const LineType& t) { return t.name == ty; });
            if (it == s.types.end()) { std::fprintf(stderr, "MoorDyn stub: unknown line type '%s'\n", ty.c_str()); return false; }
            l.type_index = static_cast<int>(it - s.types.begin());
            l.attachA = std::stoi(tok[2]); l.attachB = std::stoi(tok[3]);
            l.length = std::stod(tok[4]);
            const int nseg = std::stoi(tok[5]);
            if (nseg < 1) { std::fprintf(stderr, "MoorDyn stub: a line needs at least one segment, %d given\n", nseg); return false; }
            l.nseg = static_cast<unsigned>(nseg);
            s.lines.push_back(l);
        } else if (section == "options") {
            if (tok.size() < 2) { continue; }
            const std::string key = lower(tok[1]);
            if (key == "dtm" || key == "dt") { s.dtM = std::stod(tok[0]); }
            else if (key == "g" || key == "gravity") { s.g = std::stod(tok[0]); }
            else if (key == "rho" || key == "wtrdnsty") { s.rho = std::stod(tok[0]); }
            else if (key == "wtrdpth") { s.depth = std::stod(tok[0]); }
            else if (key == "wavekin") { s.wave_kin = std::stoi(tok[0]); }
            else if (key == "seafloorfile") { if (!read_seafloor(tok[0], s.seafloor)) { return false; } }
        }
    }
    if (s.lines.empty()) { std::fprintf(stderr, "MoorDyn stub: no lines in '%s'\n", fname.c_str()); return false; }
    for (const auto& l : s.lines) {
        if (l.attachA < 1 || l.attachA > static_cast<int>(s.points.size()) || l.attachB < 1 || l.attachB > static_cast<int>(s.points.size())) {
            std::fprintf(stderr, "MoorDyn stub: a line attaches to a point that does not exist\n"); return false;
        }
    }
    // each free point hangs from the shortest line that ties it to a held point, and swings with
    // another line attached to it
    for (std::size_t ip = 0; ip < s.points.size(); ++ip) {
        Point& p = s.points[ip];
        if (p.type != 0) { continue; }
        const int id = static_cast<int>(ip) + 1;
        for (std::size_t il = 0; il < s.lines.size(); ++il) {
            const Line& l = s.lines[il];
            const int other = (l.attachA == id) ? l.attachB : ((l.attachB == id) ? l.attachA : 0);
            if (other == 0) { continue; }
            if (s.points[static_cast<std::size_t>(other - 1)].type != 0 &&
                (p.hang_line < 0 || l.length < s.lines[static_cast<std::size_t>(p.hang_line)].length)) {
                p.hang_line = static_cast<int>(il);
            }
        }
        for (std::size_t il = 0; il < s.lines.size(); ++il) {
            const Line& l = s.lines[il];
            if (static_cast<int>(il) != p.hang_line && (l.attachA == id || l.attachB == id)) { p.swing_line = static_cast<int>(il); break; }
        }
        if (p.hang_line < 0) { std::fprintf(stderr, "MoorDyn stub: free point %d hangs from no fixed or coupled point\n", id); return false; }
        s.lines[static_cast<std::size_t>(p.hang_line)].hanging = true;
    }
    return true;
}

// Place MoorDyn line l (shape, tension and drag) from its end points, its swing angles and the
// fluid velocity U_nodes (m/s, 3 per node; null: none) over dt (s; 0: no relaxation, no velocity).
void place (StubSystem& s, Line& l, const double* U_nodes, double dt)
{
    const LineType& ty = s.types[static_cast<std::size_t>(l.type_index)];
    const Point& A = s.points[static_cast<std::size_t>(l.attachA - 1)];
    const Point& B = s.points[static_cast<std::size_t>(l.attachB - 1)];
    const unsigned nn = l.nseg + 1;
    double chord[3], c = 0.0;
    for (int d = 0; d < 3; ++d) { chord[d] = B.pos[d] - A.pos[d]; c += chord[d] * chord[d]; }
    c = std::sqrt(c);
    double ec[3] = {1.0, 0.0, 0.0};
    if (c > 0.0) { for (int d = 0; d < 3; ++d) { ec[d] = chord[d] / c; } }

    // net weight per unit length in the fluid (N/m), and the drag per unit length from the mean normal velocity (N/m)
    const double area = 0.25 * pi * ty.diam * ty.diam;
    const double w = (ty.mass - s.rho * area) * s.g;
    double Um[3] = {0.0, 0.0, 0.0};
    const double zmid = 0.5 * (A.pos[2] + B.pos[2]);
    if (U_nodes != nullptr && zmid <= 0.0) {
        for (unsigned i = 0; i < nn; ++i) { for (int d = 0; d < 3; ++d) { Um[d] += U_nodes[3*i+d] / nn; } }
    }
    const double Udotc = Um[0]*ec[0] + Um[1]*ec[1] + Um[2]*ec[2];
    double Un[3];
    for (int d = 0; d < 3; ++d) { Un[d] = Um[d] - Udotc * ec[d]; }
    const double Un2 = Un[0]*Un[0] + Un[1]*Un[1] + Un[2]*Un[2];
    const double q = 0.5 * s.rho * ty.Cd * ty.diam * Un2;

    // sag: the elastic parabola, slack or taut alike. Under the tension H along the chord the line is
    // stretched to L (1 + H / EA), and the parabola of sag w_e c^2 / (8 H) under the effective weight
    // w_e is c + 8 sag^2 / (3 c) long; the first grows with H and the second shrinks, so the H at which
    // they agree is found by bisection (on log H)
    const double w_e = std::sqrt(w * w + q * q);
    double sag = 1.0e-6 * std::max(c, 1.0);
    if (c > 0.0 && w_e > 0.0) {
        auto excess = [&](double H) {
            const double d = w_e * c * c / (8.0 * H);
            return c + 8.0 * d * d / (3.0 * c) - l.length * (1.0 + H / ty.EA);
        };
        double lo = std::log(1.0e-6), hi = std::log(1.0e3 * ty.EA);
        for (int it = 0; it < 200; ++it) {
            const double mid = 0.5 * (lo + hi);
            if (excess(std::exp(mid)) > 0.0) { lo = mid; } else { hi = mid; }
        }
        sag = std::max(sag, w_e * c * c / (8.0 * std::exp(0.5 * (lo + hi))));
    }

    // the sag direction: down, rotated about the chord towards the normal fluid velocity by the swing angle
    double down[3] = {0.0, 0.0, -1.0};
    const double ddotc = down[0]*ec[0] + down[1]*ec[1] + down[2]*ec[2];
    for (int d = 0; d < 3; ++d) { down[d] -= ddotc * ec[d]; }
    double dn = std::sqrt(down[0]*down[0] + down[1]*down[1] + down[2]*down[2]);
    if (dn < 1.0e-12) { down[0] = 0.0; down[1] = 1.0; down[2] = 0.0; dn = 1.0; }
    for (int d = 0; d < 3; ++d) { down[d] /= dn; }
    double side[3] = {0.0, 0.0, 0.0};
    if (Un2 > 0.0) {
        const double un = std::sqrt(Un2);
        for (int d = 0; d < 3; ++d) { side[d] = Un[d] / un; }
    }
    const double phi_static = (w > 0.0) ? std::atan2(q, w) : 0.0;
    const double tau = 1.0;
    l.phi = (dt > 0.0) ? phi_static + (l.phi - phi_static) * std::exp(-dt / tau) : l.phi;
    double* es = l.es;
    for (int d = 0; d < 3; ++d) { es[d] = std::cos(l.phi) * down[d] + std::sin(l.phi) * side[d]; }
    // a string the line hangs from swings across the line, by the angle of the wind across it on the
    // line's weight, with the same lag; the wind along an inclined chord swings the span in its own
    // plane but not the string, which the spans either side hold along the line
    const double nh = std::hypot(ec[0], ec[1]);
    if (nh > 1.0e-12) { l.across[0] = -ec[1] / nh; l.across[1] = ec[0] / nh; l.across[2] = 0.0; }
    const double ua = Um[0] * l.across[0] + Um[1] * l.across[1];
    const double qa = 0.5 * s.rho * ty.Cd * ty.diam * ua * std::abs(ua);
    const double phia_static = (w > 0.0) ? std::atan2(qa, w) : 0.0;
    l.phi_across = (dt > 0.0) ? phia_static + (l.phi_across - phia_static) * std::exp(-dt / tau) : l.phi_across;
    l.sag = sag;
    l.w_eff = std::sqrt(w * w + q * q);
    l.H = (sag > 0.0) ? l.w_eff * c * c / (8.0 * sag) : 0.0;

    std::vector<double> newpos(3 * nn), newten(3 * nn);
    for (unsigned i = 0; i < nn; ++i) {
        const double xi = static_cast<double>(i) / l.nseg;
        const double y = 4.0 * sag * xi * (1.0 - xi);
        const double slope = 4.0 * sag * (1.0 - 2.0 * xi) / std::max(c, 1.0e-12);
        for (int d = 0; d < 3; ++d) {
            newpos[3*i+d] = A.pos[d] + xi * chord[d] + y * es[d];
            newten[3*i+d] = l.H * (ec[d] + slope * es[d]);   // the tension along the tangent; its component along the chord is H (N)
        }
    }
    // a node below the seafloor rests on it
    if (s.seafloor.active()) {
        for (unsigned i = 0; i < nn; ++i) {
            newpos[3*i+2] = std::max(newpos[3*i+2], s.seafloor.depth_at(newpos[3*i], newpos[3*i+1]));
        }
    }
    // the drag on each node: 0.5 rho Cd D |u_n| u_n (N/m) over its share of the chord (half a segment
    // at the ends), only below z = 0
    l.drag.assign(3 * nn, 0.0);
    if (U_nodes != nullptr) {
        const double share = c / l.nseg;
        for (unsigned i = 0; i < nn; ++i) {
            if (newpos[3*i+2] > 0.0) { continue; }   // above MoorDyn's surface: no fluid
            const double udc = U_nodes[3*i]*ec[0] + U_nodes[3*i+1]*ec[1] + U_nodes[3*i+2]*ec[2];
            double un[3];
            for (int d = 0; d < 3; ++d) { un[d] = U_nodes[3*i+d] - udc * ec[d]; }
            const double mag = std::sqrt(un[0]*un[0] + un[1]*un[1] + un[2]*un[2]);
            const double len = (i == 0 || i == l.nseg) ? 0.5 * share : share;
            for (int d = 0; d < 3; ++d) { l.drag[3*i+d] = 0.5 * s.rho * ty.Cd * ty.diam * mag * un[d] * len; }
        }
    }
    l.vel.assign(3 * nn, 0.0);
    if (dt > 0.0 && l.pos.size() == newpos.size()) {
        for (std::size_t k = 0; k < newpos.size(); ++k) { l.vel[k] = (newpos[k] - l.pos[k]) / dt; }
    }
    l.pos = newpos;
    l.ten = newten;
}

void update_point_forces (StubSystem& s)
{
    // the pull of each line on its end points is the end tension directed into the line; MoorDyn
    // reports a net force for free and coupled points only, so fixed points keep zero
    for (auto& p : s.points) { for (int d = 0; d < 3; ++d) { p.force[d] = 0.0; } }
    for (const auto& l : s.lines) {
        const unsigned last = l.nseg;
        Point& A = s.points[static_cast<std::size_t>(l.attachA - 1)];
        Point& B = s.points[static_cast<std::size_t>(l.attachB - 1)];
        for (int d = 0; d < 3; ++d) {
            if (A.type != 1) { A.force[d] += l.ten[d]; }
            if (B.type != 1) { B.force[d] -= l.ten[3*last+d]; }
        }
    }
}

// the coupled points where a step of dt (s) leaves them: from x (m) at xd (m/s), as MoorDyn moves them
void set_coupled (StubSystem& s, const double* x, const double* xd, double dt)
{
    if (x == nullptr) { return; }
    unsigned ix = 0;
    for (auto& p : s.points) {
        if (p.type == -1) {
            for (int d = 0; d < 3; ++d) { p.pos[d] = x[ix + d] + (xd != nullptr ? xd[ix + d] * dt : 0.0); }
            ix += 3;
        }
    }
}

unsigned coupled_dof (const StubSystem& s)
{
    unsigned n = 0;
    for (const auto& p : s.points) { if (p.type == -1) { n += 3; } }
    return n;
}

unsigned kin_points (const StubSystem& s)
{
    // as MoorDyn-C 2.7.1: the line nodes, the points, then one entry at the origin
    unsigned n = 0;
    for (const auto& l : s.lines) { n += l.nseg + 1; }
    return n + 1 + static_cast<unsigned>(s.points.size());
}

// A hanging line (an insulator string) from its held point straight to its free point. Its
// tension, uniform along it, is its own net weight plus half the effective load (weight and drag,
// w_eff) of every other line at its free point; the normal wind drags its nodes.
void place_hanging (StubSystem& s, Line& l, const double* U_nodes)
{
    const LineType& ty = s.types[static_cast<std::size_t>(l.type_index)];
    const bool a_fixed = s.points[static_cast<std::size_t>(l.attachA - 1)].type != 0;   // held: fixed or coupled
    const Point& top = s.points[static_cast<std::size_t>((a_fixed ? l.attachA : l.attachB) - 1)];
    const Point& bot = s.points[static_cast<std::size_t>((a_fixed ? l.attachB : l.attachA) - 1)];
    const int bot_id = a_fixed ? l.attachB : l.attachA;
    const unsigned nn = l.nseg + 1;
    double e[3], len = 0.0;
    for (int d = 0; d < 3; ++d) { e[d] = bot.pos[d] - top.pos[d]; len += e[d] * e[d]; }
    len = std::sqrt(len);
    for (int d = 0; d < 3; ++d) { e[d] /= std::max(len, 1.0e-12); }
    double load = (ty.mass - s.rho * 0.25 * pi * ty.diam * ty.diam) * s.g * l.length;
    for (const auto& o : s.lines) {
        if (&o == &l || (o.attachA != bot_id && o.attachB != bot_id)) { continue; }
        if (s.points[static_cast<std::size_t>(o.attachA - 1)].loose || s.points[static_cast<std::size_t>(o.attachB - 1)].loose) {
            load += o.drop_load;   // a span with its other end detached hangs from here
            continue;
        }
        double c = 0.0;
        for (int d = 0; d < 3; ++d) {
            const double dd = s.points[static_cast<std::size_t>(o.attachB - 1)].pos[d] - s.points[static_cast<std::size_t>(o.attachA - 1)].pos[d];
            c += dd * dd;
        }
        load += 0.5 * o.w_eff * std::sqrt(c);
    }
    std::vector<double> newpos(3 * nn);
    l.ten.assign(3 * nn, 0.0);
    l.drag.assign(3 * nn, 0.0);
    for (unsigned i = 0; i < nn; ++i) {
        const double xi = static_cast<double>(i) / l.nseg;
        for (int d = 0; d < 3; ++d) {
            const double from = a_fixed ? top.pos[d] : bot.pos[d];
            const double to = a_fixed ? bot.pos[d] : top.pos[d];
            newpos[3*i+d] = from + xi * (to - from);
            // MoorDyn's tension points along the line from A to B: down the string from a fixed top
            l.ten[3*i+d] = load * (a_fixed ? e[d] : -e[d]);
        }
        if (U_nodes == nullptr || newpos[3*i+2] > 0.0) { continue; }
        const double udc = U_nodes[3*i]*e[0] + U_nodes[3*i+1]*e[1] + U_nodes[3*i+2]*e[2];
        double un[3];
        for (int d = 0; d < 3; ++d) { un[d] = U_nodes[3*i+d] - udc * e[d]; }
        const double mag = std::sqrt(un[0]*un[0] + un[1]*un[1] + un[2]*un[2]);
        const double share = (i == 0 || i == l.nseg) ? 0.5 * len / l.nseg : len / l.nseg;
        for (int d = 0; d < 3; ++d) { l.drag[3*i+d] = 0.5 * s.rho * ty.Cd * ty.diam * mag * un[d] * share; }
    }
    l.vel.assign(3 * nn, 0.0);
    l.pos = newpos;
    l.w_eff = 0.0;
}

// The seafloor under (x, y) (m, MoorDyn's frame); without a SeafloorFile there is none (-infinity).
double floor_at (const StubSystem& s, double x, double y)
{
    return s.seafloor.active() ? s.seafloor.depth_at(x, y) : -std::numeric_limits<double>::infinity();
}

bool is_loose (const StubSystem& s, int point_id) { return s.points[static_cast<std::size_t>(point_id - 1)].loose; }

// A span with a detached end (MoorDyn_BreakLine): from its held end it hangs plumb down to the seafloor and
// lies on it beyond, along l.drop; the tension at a hanging node is the weight of the hanging length
// below it, zero on the seafloor. With both ends detached it lies on the seafloor where it was. The
// loose points follow its ends. No drag.
void place_dropped (StubSystem& s, Line& l, double dt)
{
    const LineType& ty = s.types[static_cast<std::size_t>(l.type_index)];
    const unsigned nn = l.nseg + 1;
    const double w = std::max(0.0, (ty.mass - s.rho * 0.25 * pi * ty.diam * ty.diam) * s.g);
    Point& A = s.points[static_cast<std::size_t>(l.attachA - 1)];
    Point& B = s.points[static_cast<std::size_t>(l.attachB - 1)];
    std::vector<double> newpos(3 * nn);
    l.ten.assign(3 * nn, 0.0);
    l.drag.assign(3 * nn, 0.0);
    l.drop_load = 0.0;
    if (A.loose && B.loose) {
        if (l.pos.size() != newpos.size()) {
            for (unsigned i = 0; i < nn; ++i) {
                const double xi = static_cast<double>(i) / l.nseg;
                for (int d = 0; d < 3; ++d) { newpos[3*i+d] = A.pos[d] + xi * (B.pos[d] - A.pos[d]); }
            }
        } else {
            newpos = l.pos;
        }
        for (unsigned i = 0; i < nn; ++i) {
            const double f = floor_at(s, newpos[3*i], newpos[3*i+1]);
            if (std::isfinite(f)) { newpos[3*i+2] = f; }
        }
    } else {
        const bool a_held = !A.loose;
        const Point& held = a_held ? A : B;
        const double f0 = floor_at(s, held.pos[0], held.pos[1]);
        const double h = std::isfinite(f0) ? std::min(l.length, std::max(0.0, held.pos[2] - f0)) : l.length;
        l.drop_load = w * h;
        for (unsigned i = 0; i < nn; ++i) {
            const unsigned m = a_held ? i : l.nseg - i;   // the node at arc length sarc from the held end
            // the last node at the length itself: length * nseg / nseg may round above it, and a span that hangs
            // free (h = length) would then lay its end on a seafloor it does not reach
            const double sarc = (i == l.nseg) ? l.length : l.length * static_cast<double>(i) / l.nseg;
            if (sarc <= h || !std::isfinite(f0)) {
                newpos[3*m] = held.pos[0]; newpos[3*m+1] = held.pos[1]; newpos[3*m+2] = held.pos[2] - sarc;
                // MoorDyn's tension points from A to B: down the line from a held A, up it to a held B
                l.ten[3*m+2] = (a_held ? -1.0 : 1.0) * w * (h - sarc);
            } else {
                newpos[3*m] = held.pos[0] + l.drop[0] * (sarc - h);
                newpos[3*m+1] = held.pos[1] + l.drop[1] * (sarc - h);
                newpos[3*m+2] = floor_at(s, newpos[3*m], newpos[3*m+1]);
            }
        }
    }
    for (int d = 0; d < 3; ++d) {
        if (A.loose) { A.pos[d] = newpos[d]; }
        if (B.loose) { B.pos[d] = newpos[3*l.nseg+d]; }
    }
    l.vel.assign(3 * nn, 0.0);
    if (dt > 0.0 && l.pos.size() == newpos.size()) {
        for (std::size_t k = 0; k < newpos.size(); ++k) { l.vel[k] = (newpos[k] - l.pos[k]) / dt; }
    }
    l.pos = newpos;
    l.sag = 0.0; l.H = 0.0; l.w_eff = 0.0;
}

// A string detached from its held point: it hangs below its free point, down to the seafloor at most,
// with no tension and no drag; its loose point is its lower end.
void place_released (StubSystem& s, Line& l)
{
    const unsigned nn = l.nseg + 1;
    Point& A = s.points[static_cast<std::size_t>(l.attachA - 1)];
    Point& B = s.points[static_cast<std::size_t>(l.attachB - 1)];
    const bool a_loose = A.loose;
    const Point& from = a_loose ? B : A;
    const double f = floor_at(s, from.pos[0], from.pos[1]);
    l.pos.assign(3 * nn, 0.0);
    for (unsigned i = 0; i < nn; ++i) {
        const unsigned m = a_loose ? l.nseg - i : i;   // the node at arc length sarc below the free point
        const double sarc = l.length * static_cast<double>(i) / l.nseg;
        l.pos[3*m] = from.pos[0]; l.pos[3*m+1] = from.pos[1];
        l.pos[3*m+2] = std::max(from.pos[2] - sarc, std::min(f, from.pos[2]));
    }
    Point& end = a_loose ? A : B;
    const unsigned e = a_loose ? 0 : l.nseg;
    for (int d = 0; d < 3; ++d) { end.pos[d] = l.pos[3*e+d]; }
    l.ten.assign(3 * nn, 0.0);
    l.drag.assign(3 * nn, 0.0);
    l.vel.assign(3 * nn, 0.0);
    l.w_eff = 0.0;
}

// Where the free point of a string detached from its held point goes: the two spans it joins make one
// parabola between their far ends, as long as both together, and the point sits on it where the first
// span's length ends, on the seafloor at the lowest. With other than two whole spans at the point it
// drops onto the seafloor where it is.
void place_released_point (StubSystem& s, Point& p, int id, const Line& string)
{
    const Line* sp[2] = {nullptr, nullptr};
    int n = 0;
    for (const auto& o : s.lines) {
        if (&o == &string || (o.attachA != id && o.attachB != id)) { continue; }
        if (n < 2) { sp[n] = &o; }
        ++n;
    }
    bool whole = (n == 2);
    for (int k = 0; whole && k < 2; ++k) { whole = !is_loose(s, sp[k]->attachA) && !is_loose(s, sp[k]->attachB); }
    if (whole) {
        const double* P0 = s.points[static_cast<std::size_t>((sp[0]->attachA == id ? sp[0]->attachB : sp[0]->attachA) - 1)].pos;
        const double* P2 = s.points[static_cast<std::size_t>((sp[1]->attachA == id ? sp[1]->attachB : sp[1]->attachA) - 1)].pos;
        const double Lt = sp[0]->length + sp[1]->length;
        double C = 0.0;
        for (int d = 0; d < 3; ++d) { C += (P2[d] - P0[d]) * (P2[d] - P0[d]); }
        C = std::sqrt(C);
        // the sag of the elastic parabola of both spans together, as place() finds a span's: stretched to
        // Lt (1 + H / EA) under H = w C^2 / (8 S), by bisection on log H (the first span's type for both)
        const LineType& ty = s.types[static_cast<std::size_t>(sp[0]->type_index)];
        const double w = std::max(0.0, (ty.mass - s.rho * 0.25 * pi * ty.diam * ty.diam) * s.g);
        double S = 0.0;
        if (C > 0.0 && w > 0.0) {
            auto excess = [&] (double H) {
                const double sag = w * C * C / (8.0 * H);
                return C + 8.0 * sag * sag / (3.0 * C) - Lt * (1.0 + H / ty.EA);
            };
            double lo = std::log(1.0e-6), hi = std::log(1.0e3 * ty.EA);
            for (int it = 0; it < 200; ++it) {
                const double mid = 0.5 * (lo + hi);
                if (excess(std::exp(mid)) > 0.0) { lo = mid; } else { hi = mid; }
            }
            S = w * C * C / (8.0 * std::exp(0.5 * (lo + hi)));
        }
        const double xi = sp[0]->length / Lt;
        for (int d = 0; d < 3; ++d) { p.pos[d] = P0[d] + xi * (P2[d] - P0[d]); }
        p.pos[2] -= 4.0 * S * xi * (1.0 - xi);
    }
    const double f = floor_at(s, p.pos[0], p.pos[1]);
    if (std::isfinite(f)) { p.pos[2] = whole ? std::max(p.pos[2], f) : f; }
}

// Place every line: the free points first (hanging plumb below their held point, swung across by the
// angle phi_across that the span attached to them had when last placed) unless move_free_points is
// false, then the other lines, then the hanging ones.
void place_all (StubSystem& s, double dt, bool move_free_points)
{
    if (move_free_points) {
        for (std::size_t ip = 0; ip < s.points.size(); ++ip) {
            Point& p = s.points[ip];
            if (p.type != 0 || p.loose) { continue; }   // a loose point follows its line's end
            const Line& h = s.lines[static_cast<std::size_t>(p.hang_line)];
            if (is_loose(s, h.attachA) || is_loose(s, h.attachB)) {
                place_released_point(s, p, static_cast<int>(ip) + 1, h);   // its string was detached from its held point
                continue;
            }
            const Point& top = s.points[static_cast<std::size_t>((s.points[static_cast<std::size_t>(h.attachA - 1)].type != 0 ? h.attachA : h.attachB) - 1)];
            // a string hangs plumb, swung across the line by the wind across the line it carries: on a
            // sloping span the sag direction leans along the chord, and a string hung along it would
            // pull the span's end along the line and stretch a short span taut
            const Line& sw = (p.swing_line >= 0) ? s.lines[static_cast<std::size_t>(p.swing_line)] : h;
            for (int d = 0; d < 3; ++d) {
                const double plumb = (d == 2) ? -1.0 : 0.0;
                p.pos[d] = top.pos[d] + h.length * (std::cos(sw.phi_across) * plumb + std::sin(sw.phi_across) * sw.across[d]);
            }
        }
    }
    std::vector<unsigned> off(s.lines.size(), 0);
    for (std::size_t i = 1; i < s.lines.size(); ++i) { off[i] = off[i-1] + s.lines[i-1].nseg + 1; }
    const bool wind = (s.nkin > 0 && !s.U.empty());
    auto detached = [&s] (const Line& l) { return is_loose(s, l.attachA) || is_loose(s, l.attachB); };
    for (std::size_t i = 0; i < s.lines.size(); ++i) {
        Line& l = s.lines[i];
        if (l.hanging) { continue; }
        if (detached(l)) { place_dropped(s, l, dt); } else { place(s, l, wind ? &s.U[3 * off[i]] : nullptr, dt); }
    }
    for (std::size_t i = 0; i < s.lines.size(); ++i) {
        Line& l = s.lines[i];
        if (!l.hanging) { continue; }
        if (detached(l)) { place_released(s, l); } else { place_hanging(s, l, wind ? &s.U[3 * off[i]] : nullptr); }
    }
    update_point_forces(s);
}

void advance (StubSystem& s, double dt) { place_all(s, dt, true); }

void coupled_forces (const StubSystem& s, double* f)
{
    if (f == nullptr) { return; }
    unsigned ix = 0;
    for (const auto& p : s.points) {
        if (p.type == -1) { for (int d = 0; d < 3; ++d) { f[ix + d] = p.force[d]; } ix += 3; }
    }
}

StubSystem* sys (MoorDyn h) { return reinterpret_cast<StubSystem*>(h); }

} // namespace

// The whole state in memory, as doubles: the clock, then every line's angles, shape and the arrays of
// its nodes (each with its length first), then every point's position and force. The fluid velocity
// is not part of it: like MoorDyn's, it stays as last set.
namespace {
void pack (std::vector<double>& b, const std::vector<double>& v) { b.push_back(static_cast<double>(v.size())); b.insert(b.end(), v.begin(), v.end()); }
bool unpack (const double*& p, const double* end, std::vector<double>& v)
{
    if (p >= end) { return false; }
    const auto n = static_cast<std::size_t>(*p++);
    if (p + n > end) { return false; }
    v.assign(p, p + n);
    p += n;
    return true;
}
std::vector<double> stub_state (const StubSystem& s)
{
    std::vector<double> b{s.t};
    for (const auto& l : s.lines) {
        b.insert(b.end(), {l.phi, l.phi_across, l.sag, l.H, l.w_eff, l.es[0], l.es[1], l.es[2], l.across[0], l.across[1], l.across[2],
                           l.length});
        pack(b, l.pos); pack(b, l.vel); pack(b, l.ten); pack(b, l.drag);
    }
    for (const auto& pt : s.points) { b.insert(b.end(), {pt.pos[0], pt.pos[1], pt.pos[2], pt.force[0], pt.force[1], pt.force[2]}); }
    return b;
}
} // namespace

extern "C" {

MoorDyn MoorDyn_Create (const char* infilename)
{
    if (infilename == nullptr) { return nullptr; }
    auto s = std::make_unique<StubSystem>();
    // a malformed number throws in std::stod or std::stoi: report it as a failed creation, since
    // an exception must not leave this C interface
    try {
        if (!parse(infilename, *s)) { return nullptr; }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "MoorDyn stub: a malformed number in '%s' (%s)\n", infilename, e.what());
        return nullptr;
    }
    for (unsigned i = 0; i < s->lines.size(); ++i) { s->line_handles.push_back(LineHandle{s.get(), i}); }
    for (unsigned i = 0; i < s->points.size(); ++i) { s->point_handles.push_back(PointHandle{s.get(), i}); }
    return reinterpret_cast<MoorDyn>(s.release());
}

int MoorDyn_NCoupledDOF (MoorDyn system, unsigned int* n)
{
    if (system == nullptr || n == nullptr) { return MOORDYN_INVALID_VALUE; }
    *n = coupled_dof(*sys(system));
    return MOORDYN_SUCCESS;
}

int MoorDyn_SetVerbosity (MoorDyn system, int) { return system ? MOORDYN_SUCCESS : MOORDYN_INVALID_VALUE; }
int MoorDyn_SetLogFile (MoorDyn system, const char* path)
{
    if (system == nullptr) { return MOORDYN_INVALID_VALUE; }
    std::ofstream out(path, std::ios::trunc);
    return out ? MOORDYN_SUCCESS : MOORDYN_INVALID_OUTPUT_FILE;
}
int MoorDyn_SetLogLevel (MoorDyn system, int) { return system ? MOORDYN_SUCCESS : MOORDYN_INVALID_VALUE; }

static int stub_init (MoorDyn system, const double* x, bool)
{
    if (system == nullptr) { return MOORDYN_INVALID_VALUE; }
    StubSystem& s = *sys(system);
    if (coupled_dof(s) > 0 && x == nullptr) { return MOORDYN_INVALID_VALUE; }
    set_coupled(s, x, nullptr, 0.0);
    s.t = 0.0;
    for (auto& l : s.lines) { l.phi = 0.0; l.phi_across = 0.0; l.pos.clear(); l.es[0] = 0.0; l.es[1] = 0.0; l.es[2] = -1.0; }
    advance(s, 0.0);
    s.initialised = true;
    return MOORDYN_SUCCESS;
}

int MoorDyn_Init (MoorDyn system, const double* x, const double*) { return stub_init(system, x, true); }
int MoorDyn_Init_NoIC (MoorDyn system, const double* x, const double*) { return stub_init(system, x, false); }

int MoorDyn_Step (MoorDyn system, const double* x, const double* xd, double* f, double* t, double* dt)
{
    if (system == nullptr || t == nullptr || dt == nullptr) { return MOORDYN_INVALID_VALUE; }
    StubSystem& s = *sys(system);
    if (!s.initialised) { return MOORDYN_INVALID_VALUE; }
    if (coupled_dof(s) > 0 && (x == nullptr || f == nullptr)) { return MOORDYN_INVALID_VALUE; }
    if (*dt <= 0.0) { coupled_forces(s, f); return MOORDYN_SUCCESS; }
    set_coupled(s, x, xd, *dt);
    advance(s, *dt);
    s.t += *dt;
    *t = s.t;
    coupled_forces(s, f);
    return MOORDYN_SUCCESS;
}

int MoorDyn_Close (MoorDyn system)
{
    if (system == nullptr) { return MOORDYN_INVALID_VALUE; }
    delete sys(system);
    return MOORDYN_SUCCESS;
}

int MoorDyn_ExternalWaveKinInit (MoorDyn system, unsigned int* n)
{
    if (system == nullptr || n == nullptr) { return MOORDYN_INVALID_VALUE; }
    StubSystem& s = *sys(system);
    // unlike MoorDyn-C 2.7.1, which counts the points whatever WaveKin is: the stub reports 0 so
    // that a missing WaveKin = 1 is caught where only the stub runs
    s.nkin = (s.wave_kin == 1) ? kin_points(s) : 0;
    // (a point added by MoorDyn_BreakLine comes last: the velocities set before it are kept, as MoorDyn-C keeps them)
    s.U.resize(3 * static_cast<std::size_t>(s.nkin), 0.0);
    *n = s.nkin;
    return MOORDYN_SUCCESS;
}

int MoorDyn_ExternalWaveKinGetN (MoorDyn system, unsigned int* n)
{
    if (system == nullptr || n == nullptr) { return MOORDYN_INVALID_VALUE; }
    *n = sys(system)->nkin;
    return MOORDYN_SUCCESS;
}

int MoorDyn_ExternalWaveKinGetCoordinates (MoorDyn system, double* r)
{
    if (system == nullptr || r == nullptr) { return MOORDYN_INVALID_VALUE; }
    const StubSystem& s = *sys(system);
    if (s.nkin == 0) { return MOORDYN_INVALID_VALUE; }
    unsigned k = 0;
    for (const auto& l : s.lines) {
        for (std::size_t j = 0; j < l.pos.size(); ++j) { r[k++] = l.pos[j]; }
    }
    for (const auto& p : s.points) { for (int d = 0; d < 3; ++d) { r[k++] = p.pos[d]; } }
    for (int d = 0; d < 3; ++d) { r[k++] = 0.0; }   // the entry MoorDyn adds at its origin
    return MOORDYN_SUCCESS;
}

int MoorDyn_ExternalWaveKinSet (MoorDyn system, const double* U, const double*, double)
{
    if (system == nullptr || U == nullptr) { return MOORDYN_INVALID_VALUE; }
    StubSystem& s = *sys(system);
    if (s.nkin == 0) { return MOORDYN_INVALID_VALUE; }
    s.U.assign(U, U + 3 * static_cast<std::size_t>(s.nkin));
    return MOORDYN_SUCCESS;
}

int MoorDyn_GetNumberPoints (MoorDyn system, unsigned int* n)
{
    if (system == nullptr || n == nullptr) { return MOORDYN_INVALID_VALUE; }
    *n = static_cast<unsigned>(sys(system)->points.size());
    return MOORDYN_SUCCESS;
}

MoorDynPoint MoorDyn_GetPoint (MoorDyn system, unsigned int c)
{
    if (system == nullptr) { return nullptr; }
    StubSystem& s = *sys(system);
    if (c < 1 || c > s.points.size()) { return nullptr; }
    return reinterpret_cast<MoorDynPoint>(&s.point_handles[c - 1]);
}

int MoorDyn_BreakLine (MoorDyn system, MoorDynPoint point, MoorDynLine line)
{
    if (system == nullptr || point == nullptr || line == nullptr) { return MOORDYN_INVALID_VALUE; }
    StubSystem& s = *sys(system);
    const auto* ph = reinterpret_cast<PointHandle*>(point);
    const auto* lh = reinterpret_cast<LineHandle*>(line);
    if (ph->s != &s || lh->s != &s) { return MOORDYN_INVALID_VALUE; }
    Line& l = s.lines[lh->index];
    const int id = static_cast<int>(ph->index) + 1;
    if (l.attachA != id && l.attachB != id) { return MOORDYN_INVALID_VALUE; }   // the line does not end on the point
    const bool at_a = (l.attachA == id);
    // the new free point, a copy of the one the line leaves
    Point np;
    np.type = 0;
    np.loose = true;
    for (int d = 0; d < 3; ++d) { np.pos[d] = s.points[ph->index].pos[d]; }
    s.points.push_back(np);
    s.point_handles.push_back(PointHandle{&s, static_cast<unsigned>(s.points.size()) - 1});
    const int new_id = static_cast<int>(s.points.size());
    const int other = at_a ? l.attachB : l.attachA;
    (at_a ? l.attachA : l.attachB) = new_id;
    if (!l.hanging) {
        // a span: it will lie towards where this end was
        const double* from = s.points[static_cast<std::size_t>(other - 1)].pos;
        const double dx = np.pos[0] - from[0], dy = np.pos[1] - from[1];
        const double h = std::hypot(dx, dy);
        l.drop[0] = (h > 0.0) ? dx / h : 1.0;
        l.drop[1] = (h > 0.0) ? dy / h : 0.0;
    }
    // a free point that swung with this line and no longer holds it swings with another of its lines
    for (std::size_t ip = 0; ip < s.points.size(); ++ip) {
        Point& p = s.points[ip];
        if (p.type != 0 || p.loose || p.swing_line != static_cast<int>(lh->index)) { continue; }
        const int pid = static_cast<int>(ip) + 1;
        if (l.attachA == pid || l.attachB == pid) { continue; }
        p.swing_line = -1;
        for (std::size_t il = 0; il < s.lines.size(); ++il) {
            const Line& o = s.lines[il];
            if (static_cast<int>(il) != p.hang_line && (o.attachA == pid || o.attachB == pid)) { p.swing_line = static_cast<int>(il); break; }
        }
    }
    return MOORDYN_SUCCESS;
}

int MoorDyn_GetNumberLines (MoorDyn system, unsigned int* n)
{
    if (system == nullptr || n == nullptr) { return MOORDYN_INVALID_VALUE; }
    *n = static_cast<unsigned>(sys(system)->lines.size());
    return MOORDYN_SUCCESS;
}

MoorDynLine MoorDyn_GetLine (MoorDyn system, unsigned int l)
{
    if (system == nullptr) { return nullptr; }
    StubSystem& s = *sys(system);
    if (l < 1 || l > s.lines.size()) { return nullptr; }
    return reinterpret_cast<MoorDynLine>(&s.line_handles[l - 1]);
}

int MoorDyn_GetDt (MoorDyn system, double* dt)
{
    if (system == nullptr || dt == nullptr) { return MOORDYN_INVALID_VALUE; }
    *dt = sys(system)->dtM;
    return MOORDYN_SUCCESS;
}

int MoorDyn_SetDt (MoorDyn system, double dt)
{
    if (system == nullptr || dt <= 0.0) { return MOORDYN_INVALID_VALUE; }
    sys(system)->dtM = dt;
    return MOORDYN_SUCCESS;
}

int MoorDyn_Save (MoorDyn system, const char* filepath)
{
    if (system == nullptr || filepath == nullptr) { return MOORDYN_INVALID_VALUE; }
    const StubSystem& s = *sys(system);
    std::ofstream out(filepath, std::ios::trunc);
    if (!out) { return MOORDYN_INVALID_OUTPUT_FILE; }
    out.precision(17);
    out << "moordyn-stub-state-2 " << s.t << " " << s.lines.size() << " " << s.U.size() << "\n";
    for (const auto& l : s.lines) {
        // the unstretched length too: it may have changed since the input (MoorDyn_SetLineUnstretchedLength)
        out << l.phi << " " << l.phi_across << " " << l.length << " " << l.drop[0] << " " << l.drop[1];
        for (double v : l.pos) { out << " " << v; }
        out << "\n";
    }
    // the fluid velocity of the last step sets the direction the lines swing to: without it a
    // restored line would hang in the plane of its swing angle but towards no wind
    for (double v : s.U) { out << v << " "; }
    out << "\n";
    // where the free points hang (they follow their lines' swing of the step before) and where the
    // coupled points were moved to
    for (const auto& p : s.points) { if (p.type != 1) { out << p.pos[0] << " " << p.pos[1] << " " << p.pos[2] << "\n"; } }
    return MOORDYN_SUCCESS;
}

int MoorDyn_Load (MoorDyn system, const char* filepath)
{
    if (system == nullptr || filepath == nullptr) { return MOORDYN_INVALID_VALUE; }
    StubSystem& s = *sys(system);
    if (!s.initialised) { return MOORDYN_INVALID_VALUE; }
    std::ifstream in(filepath);
    if (!in) { return MOORDYN_INVALID_INPUT_FILE; }
    std::string tag;
    std::size_t nl = 0, nu = 0;
    // "moordyn-stub-state": a file written before MoorDyn_BreakLine was added, without the drop directions (no line was detached)
    if (!(in >> tag >> s.t >> nl >> nu) || (tag != "moordyn-stub-state-2" && tag != "moordyn-stub-state") || nl != s.lines.size()) {
        return MOORDYN_INVALID_INPUT;
    }
    const bool with_drop = (tag == "moordyn-stub-state-2");
    for (auto& l : s.lines) {
        if (!(in >> l.phi >> l.phi_across >> l.length) || !(l.length > 0.0)) { return MOORDYN_INVALID_INPUT; }
        if (with_drop && !(in >> l.drop[0] >> l.drop[1])) { return MOORDYN_INVALID_INPUT; }
        for (double& v : l.pos) { if (!(in >> v)) { return MOORDYN_INVALID_INPUT; } }
    }
    // the fluid velocity of the last step, for the points ExternalWaveKinInit set up
    if (nu != s.U.size()) { return MOORDYN_INVALID_INPUT; }
    for (double& v : s.U) { if (!(in >> v)) { return MOORDYN_INVALID_INPUT; } }
    for (auto& p : s.points) {
        if (p.type != 1) { for (double& v : p.pos) { if (!(in >> v)) { return MOORDYN_INVALID_INPUT; } } }
    }
    // the shape, tensions and drag follow from the restored angles, points and fluid velocity;
    // velocities restart from rest
    place_all(s, 0.0, false);
    return MOORDYN_SUCCESS;
}



int MoorDyn_Serialize (MoorDyn system, size_t* size, uint64_t* data)
{
    if (system == nullptr) { return MOORDYN_INVALID_VALUE; }
    const std::vector<double> b = stub_state(*sys(system));
    // a leading count, then the doubles' bits
    if (size != nullptr) { *size = (b.size() + 1) * sizeof(uint64_t); }
    if (data != nullptr) {
        data[0] = static_cast<uint64_t>(b.size());
        std::memcpy(data + 1, b.data(), b.size() * sizeof(double));
    }
    return MOORDYN_SUCCESS;
}

int MoorDyn_Deserialize (MoorDyn system, const uint64_t* data)
{
    if (system == nullptr || data == nullptr) { return MOORDYN_INVALID_VALUE; }
    StubSystem& s = *sys(system);
    std::vector<double> b(static_cast<std::size_t>(data[0]));
    std::memcpy(b.data(), data + 1, b.size() * sizeof(double));
    const double* p = b.data();
    const double* end = b.data() + b.size();
    if (p >= end) { return MOORDYN_INVALID_INPUT; }
    s.t = *p++;
    for (auto& l : s.lines) {
        if (p + 12 > end) { return MOORDYN_INVALID_INPUT; }
        l.phi = p[0]; l.phi_across = p[1]; l.sag = p[2]; l.H = p[3]; l.w_eff = p[4];
        for (int d = 0; d < 3; ++d) { l.es[d] = p[5+d]; l.across[d] = p[8+d]; }
        l.length = p[11];
        p += 12;
        if (!unpack(p, end, l.pos) || !unpack(p, end, l.vel) || !unpack(p, end, l.ten) || !unpack(p, end, l.drag)) { return MOORDYN_INVALID_INPUT; }
    }
    for (auto& pt : s.points) {
        if (p + 6 > end) { return MOORDYN_INVALID_INPUT; }
        for (int d = 0; d < 3; ++d) { pt.pos[d] = p[d]; pt.force[d] = p[3+d]; }
        p += 6;
    }
    return (p == end) ? MOORDYN_SUCCESS : MOORDYN_INVALID_INPUT;
}


#define STUB_LINE(h) if ((h) == nullptr) { return MOORDYN_INVALID_VALUE; } \
    const Line& L = reinterpret_cast<LineHandle*>(h)->s->lines[reinterpret_cast<LineHandle*>(h)->index]

int MoorDyn_GetLineN (MoorDynLine l, unsigned int* n) { STUB_LINE(l); if (!n) { return MOORDYN_INVALID_VALUE; } *n = L.nseg; return MOORDYN_SUCCESS; }
int MoorDyn_GetLineNumberNodes (MoorDynLine l, unsigned int* n) { STUB_LINE(l); if (!n) { return MOORDYN_INVALID_VALUE; } *n = L.nseg + 1; return MOORDYN_SUCCESS; }
int MoorDyn_GetLineUnstretchedLength (MoorDynLine l, double* ul) { STUB_LINE(l); if (!ul) { return MOORDYN_INVALID_VALUE; } *ul = L.length; return MOORDYN_SUCCESS; }
// the new length shapes the line from its next placement on (the next step), as MoorDyn-C's new segment lengths act from its next step
int MoorDyn_SetLineUnstretchedLength (MoorDynLine l, double ul)
{
    if (l == nullptr || !(ul > 0.0) || !std::isfinite(ul)) { return MOORDYN_INVALID_VALUE; }
    reinterpret_cast<LineHandle*>(l)->s->lines[reinterpret_cast<LineHandle*>(l)->index].length = ul;
    return MOORDYN_SUCCESS;
}

int MoorDyn_GetLineNodePos (MoorDynLine l, unsigned int i, double pos[3])
{
    STUB_LINE(l);
    if (pos == nullptr || i > L.nseg) { return MOORDYN_INVALID_VALUE; }
    for (int d = 0; d < 3; ++d) { pos[d] = L.pos[3*i+d]; }
    return MOORDYN_SUCCESS;
}

int MoorDyn_GetLineNodeVel (MoorDynLine l, unsigned int i, double vel[3])
{
    STUB_LINE(l);
    if (vel == nullptr || i > L.nseg) { return MOORDYN_INVALID_VALUE; }
    for (int d = 0; d < 3; ++d) { vel[d] = L.vel[3*i+d]; }
    return MOORDYN_SUCCESS;
}

int MoorDyn_GetLineNodeTen (MoorDynLine l, unsigned int i, double t[3])
{
    STUB_LINE(l);
    if (t == nullptr || i > L.nseg) { return MOORDYN_INVALID_VALUE; }
    for (int d = 0; d < 3; ++d) { t[d] = L.ten[3*i+d]; }
    return MOORDYN_SUCCESS;
}

int MoorDyn_GetLineNodeDrag (MoorDynLine l, unsigned int i, double f[3])
{
    STUB_LINE(l);
    if (f == nullptr || i > L.nseg) { return MOORDYN_INVALID_VALUE; }
    for (int d = 0; d < 3; ++d) { f[d] = L.drag.empty() ? 0.0 : L.drag[3*i+d]; }
    return MOORDYN_SUCCESS;
}

// the stub's lines are in equilibrium: an inner node carries no net force, and an end node the
// pull of the line on its support, the tension along the line's tangent there (which carries the
// line's weight and drag as the parabola's slope does)
int MoorDyn_GetLineNodeForce (MoorDynLine l, unsigned int i, double f[3])
{
    STUB_LINE(l);
    if (f == nullptr || i > L.nseg) { return MOORDYN_INVALID_VALUE; }
    const double sign = (i == 0) ? 1.0 : (i == L.nseg ? -1.0 : 0.0);
    for (int d = 0; d < 3; ++d) { f[d] = sign * L.ten[3*i+d]; }
    return MOORDYN_SUCCESS;
}

int MoorDyn_GetLineFairTen (MoorDynLine l, double* t)
{
    STUB_LINE(l);
    if (t == nullptr) { return MOORDYN_INVALID_VALUE; }
    const unsigned i = L.nseg;
    *t = std::sqrt(L.ten[3*i]*L.ten[3*i] + L.ten[3*i+1]*L.ten[3*i+1] + L.ten[3*i+2]*L.ten[3*i+2]);
    return MOORDYN_SUCCESS;
}

int MoorDyn_GetLineMaxTen (MoorDynLine l, double* t)
{
    STUB_LINE(l);
    if (t == nullptr) { return MOORDYN_INVALID_VALUE; }
    double m = 0.0;
    for (unsigned i = 0; i <= L.nseg; ++i) {
        m = std::max(m, std::sqrt(L.ten[3*i]*L.ten[3*i] + L.ten[3*i+1]*L.ten[3*i+1] + L.ten[3*i+2]*L.ten[3*i+2]));
    }
    *t = m;
    return MOORDYN_SUCCESS;
}

#define STUB_POINT(h) if ((h) == nullptr) { return MOORDYN_INVALID_VALUE; } \
    const Point& P = reinterpret_cast<PointHandle*>(h)->s->points[reinterpret_cast<PointHandle*>(h)->index]

int MoorDyn_GetPointType (MoorDynPoint p, int* t) { STUB_POINT(p); if (!t) { return MOORDYN_INVALID_VALUE; } *t = P.type; return MOORDYN_SUCCESS; }
int MoorDyn_GetPointPos (MoorDynPoint p, double pos[3]) { STUB_POINT(p); if (!pos) { return MOORDYN_INVALID_VALUE; } for (int d = 0; d < 3; ++d) { pos[d] = P.pos[d]; } return MOORDYN_SUCCESS; }
int MoorDyn_GetPointForce (MoorDynPoint p, double f[3]) { STUB_POINT(p); if (!f) { return MOORDYN_INVALID_VALUE; } for (int d = 0; d < 3; ++d) { f[d] = P.force[d]; } return MOORDYN_SUCCESS; }

} // extern "C"
