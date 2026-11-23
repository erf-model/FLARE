#!/usr/bin/env python3
"""Write the Gaussian hill of the FireFrozenWind case as an ERF terrain file.

    python3 make_hill_terrain.py [terrain_gaussian_hill.txt]

z(x, y) = H exp(-((x - xc)^2 + (y - yc)^2) / (2 s^2)), H = 40 m, s = 150 m,
centred at (640, 640) m on 0 <= x, y <= 1280 m, nodes every 20 m (65 x 65).
The hill is 0.004 m high at the domain edge, so the surface is periodic to
well under a millimetre. Format (ERF's erf.terrain_file_name, GUST's
terrain.type = erf_grid): nx, ny, the nx x values, the ny y values, then the
nx*ny heights with x varying fastest, one number per line. The hill is
symmetric in x and y, so the order of the heights does not matter here.
Standard library only.
"""
import math
import sys

H, S, XC, YC = 40.0, 150.0, 640.0, 640.0
L, N = 1280.0, 65


def height(x, y):
    return H * math.exp(-((x - XC) ** 2 + (y - YC) ** 2) / (2.0 * S * S))


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "terrain_gaussian_hill.txt"
    xs = [L * i / (N - 1) for i in range(N)]
    with open(out, "w") as f:
        f.write(f"{N}\n{N}\n")
        for x in xs:
            f.write(f"{x:.6f}\n")
        for y in xs:
            f.write(f"{y:.6f}\n")
        for y in xs:
            for x in xs:
                f.write(f"{height(x, y):.6f}\n")
    print(f"wrote {out}: {N} x {N} nodes, peak {height(XC, YC):.1f} m")


if __name__ == "__main__":
    main()
