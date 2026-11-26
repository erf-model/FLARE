#!/usr/bin/env python3
"""Check that fire_lines_towers.py reads the ERF terrain file with y varying fastest.

    python3 test_fire_lines_towers.py

Writes the bilinear surface z = 0.1 x + 0.01 y + 1e-4 x y on a 11 x 21 grid (nx != ny, so a transposed read
gives a wrong height) and checks the ground height against it at points off the grid nodes, on the grid edges
and at the corners, where bilinear interpolation is exact. Exits 1 on a mismatch.
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
sys.exit(0 if np.allclose(g, want, rtol=0, atol=1e-9) else 1)
