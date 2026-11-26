#!/usr/bin/env python3
"""Pick transmission towers along OpenStreetMap power lines and write the FLARE conductor block of an inputs deck.

    fire_lines_towers.py lines_config.json --osm osm_lines.json --terrain terrain.txt [--meta lines_meta.json]
                         [--check inputs]

lines_config.json (one per case folder):
    {"air_density": 1.0, "elevation_m": 1700,
     "lines": [{"name": "L1", "osm_way": 1467202772, "box_m": [x0, x1, y0, y1], "min_span_m": 150, "frames": true}, ...]}
osm_lines.json is the output of osm_power_lines.py (domain metres); terrain.txt is the ERF terrain file of the case.

For each line: keeps the longest run of OSM nodes inside box_m and, from each tower, picks the farthest node whose
span is at most 550 m long, at least min_span_m long and clears the terrain by 8 m in still air (catenary sag of
Drake ACSR at 30 kN). Towers are 30 m above the local ground; both ends of a span that still clears by less than 8 m
are raised (at most to 60 m). Tower points are heights above the local terrain (ConductorInputs.rst).
Prints the conductor block (from the "TRANSMISSION LINES" banner to the end of the deck) and writes lines_meta.json.
With --check, compares the block with the same section of an existing deck, and the towers with the existing
--meta file (which is then not written), and exits 1 if either differs.
Exits with an error if a line's way is missing, fewer than 2 of its nodes are in the box, or a span is still
longer than 550 m or clears the terrain by less than 8 m after the towers are raised.
The terrain file is the ERF format: nx, ny, the x and y coordinates, then z with y varying fastest.
"""
import argparse, json, sys
import numpy as np

TOWER_H, TENSION, W = 30.0, 30000.0, 1.628 * 9.81
SPAN_MAX, CLEAR_MIN, H_MAX = 550.0, 8.0, 60.0


def read_terrain(fn):
    with open(fn) as f:
        nx = int(f.readline()); ny = int(f.readline())
        v = np.array(f.read().split(), dtype=float)
    if v.size != nx + ny + nx * ny:
        sys.exit(f"{fn}: expected {nx + ny + nx * ny} values after nx, ny; found {v.size}")
    x, y, z = v[:nx], v[nx:nx + ny], v[nx + ny:].reshape(nx, ny)      # z[i, j]: y varies fastest
    return x, y, z


def ground(T, px, py):
    x, y, z = T
    i = np.clip(np.searchsorted(x, px) - 1, 0, len(x) - 2); j = np.clip(np.searchsorted(y, py) - 1, 0, len(y) - 2)
    fx = (px - x[i]) / (x[i + 1] - x[i]); fy = (py - y[j]) / (y[j + 1] - y[j])
    return (z[i, j] * (1 - fx) * (1 - fy) + z[i + 1, j] * fx * (1 - fy) + z[i, j + 1] * (1 - fx) * fy + z[i + 1, j + 1] * fx * fy)


def pick(lines, wid, box, min_span, T):
    l = next((l for l in lines if l["id"] == wid), None)
    if l is None:
        sys.exit(f"OSM way {wid} is not in the --osm file")
    p = np.array(l["pts"])[:, :2]
    ins = (p[:, 0] > box[0]) & (p[:, 0] < box[1]) & (p[:, 1] > box[2]) & (p[:, 1] < box[3])
    idx = np.where(ins)[0]
    # the longest consecutive run inside the box
    runs = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)
    q = p[max(runs, key=len)]
    if len(q) < 2:
        sys.exit(f"OSM way {wid}: fewer than 2 nodes inside box_m {box}")
    # greedy: from each tower, the farthest OSM node whose span clears the terrain by 8 m and is at most
    # 550 m long; the next node when none does
    keep, i = [q[0]], 0
    while i < len(q) - 1:
        j = i + 1
        for k in range(i + 1, len(q)):
            Lk = np.hypot(*(q[k] - q[i]))
            if Lk > SPAN_MAX:
                break
            if Lk >= min_span or k == len(q) - 1:
                if clearance(T, np.array([q[i], q[k]]))[0][1] >= CLEAR_MIN:
                    j = k
        if j == i + 1 and np.hypot(*(q[j] - q[i])) < min_span and j < len(q) - 1:
            print(f"way {wid}: no span from node {i} meets min_span and the 8 m clearance; using the next node",
                  file=sys.stderr)
        keep.append(q[j]); i = j
    return l, np.array(keep)


def clearance(T, pts, h=None):
    out = []
    h = [TOWER_H] * len(pts) if h is None else h
    for n, (a, b) in enumerate(zip(pts[:-1], pts[1:])):
        Lh = np.hypot(*(b - a)); s = np.linspace(0, 1, 41)
        xs, ys = a[0] + s * (b[0] - a[0]), a[1] + s * (b[1] - a[1])
        g = ground(T, xs, ys)
        za, zb = g[0] + h[n] - 2.5, g[-1] + h[n + 1] - 2.5
        cable = za + s * (zb - za) - W * (s * Lh) * (Lh - s * Lh) / (2 * TENSION)
        out.append((Lh, float((cable - g).min())))
    return out


def line_block(nm, pts, h, frames):
    f = []
    xyz = [(p[0], p[1], h[k]) for k, p in enumerate(pts)]
    f.append(f"conductors.{nm}.end_a  = {xyz[0][0]:.1f} {xyz[0][1]:.1f} {xyz[0][2]:g}")
    f.append(f"conductors.{nm}.towers = " + "  ".join(f"{x:.1f} {y:.1f} {z:g}" for x, y, z in xyz[1:-1]))
    f.append(f"conductors.{nm}.end_b  = {xyz[-1][0]:.1f} {xyz[-1][1]:.1f} {xyz[-1][2]:g}")
    f += [f"conductors.{nm}.stringing_tension = {TENSION:g}",
          f"conductors.{nm}.diameter         = 0.0281     # Drake ACSR, IEEE Std 738-2023 4.6",
          f"conductors.{nm}.mass_per_length  = 1.628",
          f"conductors.{nm}.axial_stiffness  = 3.0e7",
                    (f"conductors.{nm}.tower_type       = lattice    # frame towers: steel temperature" if frames else
           f"# {nm}: rigid point towers (no frame, no steel temperature), to keep the run short"),
          f"conductors.{nm}.insulator_length = 2.5",
          f"conductors.{nm}.insulator_mass   = 60.",
          f"conductors.{nm}.segments         = 12",
          f"conductors.{nm}.output_root      = {nm}",
          f"conductors.{nm}.absorptivity     = 0.5",
          f"conductors.{nm}.emissivity       = 0.5",
          f"conductors.{nm}.heat_capacity    = 1310.0",
          f"conductors.{nm}.thermal_expansion = 1.93e-5",
          f"conductors.{nm}.current          = 600.0",
          f"conductors.{nm}.resistance       = 7.283e-5 8.688e-5",
          f"conductors.{nm}.resistance_temperatures = 25.0 75.0"]
    return "\n".join(f)


COMMON_STRUCT = """
# ============================================================================
# TRANSMISSION LINES (OpenStreetMap ways, (c) OpenStreetMap contributors, ODbL)
# ============================================================================
amrex.fpe_trap_invalid = 0          # MoorDyn's initial-condition solver overflows otherwise
erf.structures.model   = conductors
conductors.lines       = {names}
conductors.tower_types = lattice
conductors.diagnostics_dir = conductors
conductors.diagnostics_int = 5
conductors.node_output_int = 30
conductors.air_density     = {rho}
conductors.surface_offset  = 20000     # above the domain top
conductors.drag_on_flow    = false
conductors.lattice.base_width    = 6
conductors.lattice.top_width     = 1.5
conductors.lattice.solidity      = 0.2
conductors.lattice.arm_length    = 12
conductors.lattice.arm_depth     = 1.2
conductors.lattice.frame_panels  = 4
conductors.lattice.leg_angle     = 0.15 0.012
conductors.lattice.brace_angle   = 0.10 0.008
conductors.lattice.damping_ratio = 0.02
conductors.thermal_model              = ieee738
conductors.thermal.radiative_fraction = 0.3
conductors.thermal.cutoff_radius      = 300.0
conductors.thermal.source_height      = 1.0
conductors.thermal.elevation          = {elev}
conductors.thermal.towers             = true
{blocks}
"""


def build(cfg, lines, T):
    blocks, meta = [], []
    for c in cfg["lines"]:
        nm, wid, box, ms, frames = c["name"], int(c["osm_way"]), c["box_m"], float(c["min_span_m"]), bool(c["frames"])
        l, pts = pick(lines, wid, box, ms, T)
        h = [TOWER_H] * len(pts)
        for _ in range(20):
            cl = clearance(T, pts, h)
            bad = [n for n, (_, k) in enumerate(cl) if k < CLEAR_MIN]
            if not bad:
                break
            for n in bad:
                d = CLEAR_MIN - cl[n][1]
                h[n] = min(H_MAX, h[n] + d); h[n + 1] = min(H_MAX, h[n + 1] + d)
        cl = clearance(T, pts, h)
        for n, (s_, k_) in enumerate(cl):
            if s_ > SPAN_MAX or k_ < CLEAR_MIN:
                sys.exit(f"{nm}: span {n} is {s_:.0f} m long and clears the terrain by {k_:.1f} m with towers at most "
                         f"{H_MAX:g} m; limits are {SPAN_MAX:g} m and {CLEAR_MIN:g} m")
        print(f"{nm} way {wid} {l['voltage']} {l['name']}: {len(pts) - 1} spans, span {min(s for s, _ in cl):.0f}-"
              f"{max(s for s, _ in cl):.0f} m, min clearance {min(k for _, k in cl):.1f} m, tower heights "
              f"{[round(x) for x in h]}", file=sys.stderr)
        blocks.append(f"# {nm}: OSM way {wid}, {l['voltage']} V, {l['name']}\n" + line_block(nm, pts, h, frames))
        meta.append(dict(name=nm, way=wid, voltage=l["voltage"], osm_name=l["name"], pts=pts.tolist(), heights=h,
                         frames=frames, spans=[s for s, _ in cl], clearance=[k for _, k in cl]))
    block = COMMON_STRUCT.format(names=" ".join(c["name"] for c in cfg["lines"]), rho=cfg["air_density"],
                                 elev=cfg["elevation_m"], blocks="\n".join(blocks))
    return block, meta


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config"); ap.add_argument("--osm", required=True); ap.add_argument("--terrain", required=True)
    ap.add_argument("--meta", default="lines_meta.json"); ap.add_argument("--check", default="")
    a = ap.parse_args()
    cfg = json.load(open(a.config))
    for k in ("air_density", "elevation_m", "lines"):
        if k not in cfg:
            sys.exit(f"{a.config}: missing key '{k}'")
    if not 0.5 <= cfg["air_density"] <= 1.5:
        sys.exit(f"{a.config}: air_density must be in [0.5, 1.5] kg/m3")
    if not -500 <= cfg["elevation_m"] <= 5000:
        sys.exit(f"{a.config}: elevation_m must be in [-500, 5000] m")
    names = [c.get("name") for c in cfg["lines"]]
    if not names or len(set(names)) != len(names):
        sys.exit(f"{a.config}: lines must be a non-empty list with unique names")
    for c in cfg["lines"]:
        for k in ("name", "osm_way", "box_m", "min_span_m", "frames"):
            if k not in c:
                sys.exit(f"{a.config}: line {c.get('name', '?')}: missing key '{k}'")
        if not isinstance(c["frames"], bool):
            sys.exit(f"{a.config}: line {c['name']}: frames must be true or false")
        if len(c["box_m"]) != 4:
            sys.exit(f"{a.config}: line {c['name']}: box_m must have 4 numbers")
        if not (c["box_m"][0] < c["box_m"][1] and c["box_m"][2] < c["box_m"][3]):
            sys.exit(f"{a.config}: line {c['name']}: box_m must be [x0, x1, y0, y1] with x0 < x1 and y0 < y1")
        if not 0 < float(c["min_span_m"]) <= SPAN_MAX:
            sys.exit(f"{a.config}: line {c['name']}: min_span_m must be in (0, {SPAN_MAX:g}]")
    block, meta = build(cfg, json.load(open(a.osm)), read_terrain(a.terrain))
    if not a.check:
        json.dump(meta, open(a.meta, "w"), indent=1)
    else:
        old = json.load(open(a.meta))
        if json.dumps(old, sort_keys=True) != json.dumps(json.loads(json.dumps(meta)), sort_keys=True):
            sys.exit(f"{a.meta}: towers differ from the ones built from {a.config}")
        deck = open(a.check).read()
        banner = "# " + "=" * 76 + "\n# TRANSMISSION LINES"
        i = deck.find(banner)
        if i < 0 or deck[i:].rstrip() != block[block.index(banner):].rstrip():
            sys.exit(f"{a.check}: conductor block differs from the one built from {a.config}")
        print(f"{a.check}: conductor block and {a.meta} match", file=sys.stderr)
    if not a.check:
        sys.stdout.write(block)
