#!/usr/bin/env python3
"""Check that fire_lines_towers.py reads the ERF terrain file with y varying fastest.

    python3 test_fire_lines_towers.py

Writes a plane z = 0.1 x + 0.01 y on a 11 x 21 grid (nx != ny, so a transposed read gives a wrong height or
fails) and checks the bilinear ground height at one point against the plane. Exits 1 on a mismatch.
"""
import importlib.util, os, sys, tempfile
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("towers", os.path.join(here, "fire_lines_towers.py"))
t = importlib.util.module_from_spec(spec); spec.loader.exec_module(t)

x = np.linspace(0.0, 1000.0, 11); y = np.linspace(0.0, 2000.0, 21)
X, Y = np.meshgrid(x, y, indexing="ij")
Z = 0.1 * X + 0.01 * Y
with tempfile.TemporaryDirectory() as d:
    fn = os.path.join(d, "terrain.txt")
    with open(fn, "w") as f:
        f.write(f"{len(x)}\n{len(y)}\n")
        f.write("\n".join(map(str, x)) + "\n" + "\n".join(map(str, y)) + "\n")
        f.write("\n".join(map(str, Z.ravel())) + "\n")        # index i * ny + j: y varies fastest
    g = t.ground(t.read_terrain(fn), np.array([700.0]), np.array([1300.0]))[0]
want = 0.1 * 700.0 + 0.01 * 1300.0
print(f"ground {g:.6f} m, plane {want:.6f} m")
sys.exit(0 if abs(g - want) < 1e-9 else 1)
