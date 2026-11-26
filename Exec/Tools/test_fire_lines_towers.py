#!/usr/bin/env python3
"""Check fire_lines_towers.py: the terrain read, the run's mesh surface, the still-air clearance and the tower picks.

    python3 test_fire_lines_towers.py

Writes the bilinear surface z = 0.1 x + 0.01 y + 1e-4 x y on a 11 x 21 grid (nx != ny, so a transposed read
gives a wrong height) and checks the ground height against it at points off the grid nodes, on the grid edges
and at the corners, where bilinear interpolation is exact. Then builds the mesh surface for a deck with
prob_hi = 1000 2000 and n_cell = 4 5 and checks its 5 x 6 nodes against the surface.
Over flat ground the clearance of a span of length L is the closed form 30 - 2.5 - W L^2 / (8 T) (tower height,
insulator, parabolic sag of Drake ACSR, W = 1.628 x 9.81 N/m, at T = 30 kN): 10.864 m at 500 m. The picks along a
straight line of nodes 100 m apart are the farthest node at most 550 m away (spans 500, 500, 100 m). With a
second surface 25 m above the ground to clear, no span of 150 m or more clears 8 m, so every node is a tower
(11 spans of 100 m, with a warning), and every tower is raised by 8 - (27.5 - 0.665 - 25) = 6.165 m to 36.165 m,
where every span clears exactly 8 m. A repeated node and a line with fewer than 3 towers stop the tool. Exits 1 on any mismatch.
"""
import importlib.util, os, sys, tempfile
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("towers", os.path.join(here, "fire_lines_towers.py"))
t = importlib.util.module_from_spec(spec); spec.loader.exec_module(t)

x = np.linspace(0.0, 1000.0, 11); y = np.linspace(0.0, 2000.0, 21)
X, Y = np.meshgrid(x, y, indexing="ij")
surf = lambda px, py: 0.1 * px + 0.01 * py + 1e-4 * px * py
Z = surf(X, Y)
with tempfile.TemporaryDirectory() as d:
    fn = os.path.join(d, "terrain.txt")
    with open(fn, "w") as f:
        f.write(f"{len(x)}\n{len(y)}\n")
        f.write("\n".join(map(str, x)) + "\n" + "\n".join(map(str, y)) + "\n")
        f.write("\n".join(map(str, Z.ravel())) + "\n")        # index i * ny + j: y varies fastest
    T = t.read_terrain(fn)
px = np.array([735.0, 12.0, 0.0, 1000.0, 450.0, 1000.0])
py = np.array([1310.0, 1987.0, 0.0, 2000.0, 0.0, 777.0])
g = t.ground(T, px, py); want = surf(px, py)
for a, b, u, w in zip(px, py, g, want):
    print(f"({a:6.0f}, {b:6.0f}): ground {u:.6f} m, surface {w:.6f} m")
ok = np.allclose(g, want, rtol=0, atol=1e-9)
with tempfile.TemporaryDirectory() as d:
    fn = os.path.join(d, "inputs")
    open(fn, "w").write("geometry.prob_lo = 0.0 0.0 0.0\ngeometry.prob_hi = 1000.0 2000.0 900.0   # top\namr.n_cell = 4 5 9\n")
    xn, yn, zn = t.mesh_surface(T, *t.read_deck(fn))
XN, YN = np.meshgrid(xn, yn, indexing="ij")
print(f"mesh surface: {zn.shape[0]} x {zn.shape[1]} nodes, max error {np.abs(zn - surf(XN, YN)).max():.2e} m")
ok = ok and zn.shape == (5, 6) and np.allclose(xn, np.linspace(0, 1000, 5)) and np.allclose(zn, surf(XN, YN), rtol=0, atol=1e-9)

# still-air clearance over flat ground: closed form
flat = (np.array([-100.0, 2000.0]), np.array([-100.0, 100.0]), np.zeros((2, 2)))
L = 500.0
want = t.TOWER_H - 2.5 - t.W * L ** 2 / (8 * t.TENSION)
got = t.clearance((flat, []), np.array([[0.0, 0.0], [L, 0.0]]))[0]
print(f"flat {L:.0f} m span: clearance {got[1]:.6f} m, closed form {want:.6f} m")
ok = ok and abs(got[0] - L) < 1e-9 and abs(got[1] - want) < 1e-9 and abs(want - 10.8639) < 1e-4

# picks along a straight line, then tower raising over a surface 25 m higher
way = dict(id=1, voltage="230000", name="test", pts=[[100.0 * k, 0.0] for k in range(12)])
cfg = dict(air_density=1.0, elevation_m=0.0, lines=[dict(name="T1", osm_way=1, box_m=[-50.0, 1150.0, -50.0, 50.0],
                                                         min_span_m=150.0, frames=False)])
_, pts = t.pick([way], 1, cfg["lines"][0]["box_m"], 150.0, (flat, []))
spans = list(np.hypot(*np.diff(pts, axis=0).T))
print(f"picked spans {spans}")
ok = ok and np.allclose(spans, [500.0, 500.0, 100.0])
_, meta = t.build(cfg, [way], (flat, []))
ok = ok and meta[0]["heights"] == [t.TOWER_H] * 4
high = (flat[0], flat[1], flat[2] + 25.0)
_, meta = t.build(cfg, [way], (flat, [high]))
cl, h = meta[0]["clearance"], meta[0]["heights"]
h_want = t.TOWER_H + t.CLEAR_MIN - (t.TOWER_H - 2.5 - t.W * 100.0 ** 2 / (8 * t.TENSION) - 25.0)
print(f"over a surface 25 m higher: {len(cl)} spans, tower heights {min(h):.4f}-{max(h):.4f} m (closed form "
      f"{h_want:.4f} m), clearance {min(cl):.6f}-{max(cl):.6f} m")
ok = ok and len(cl) == 11 and np.allclose(h, h_want, rtol=0, atol=1e-9) and np.allclose(cl, t.CLEAR_MIN, rtol=0, atol=1e-9)
ok = ok and abs(h_want - 36.1654) < 1e-4

# guards: a repeated node, fewer than 3 towers
for pts_, what in (([[0.0, 0.0], [300.0, 0.0], [300.2, 0.0], [600.0, 0.0]], "repeated node"),
                   ([[0.0, 0.0], [400.0, 0.0]], "two towers")):
    try:
        t.pick([dict(way, pts=pts_)], 1, cfg["lines"][0]["box_m"], 150.0, (flat, []))
        print(f"{what}: not stopped"); ok = False
    except SystemExit as e:
        print(f"{what}: stopped ({e})")
sys.exit(0 if ok else 1)
