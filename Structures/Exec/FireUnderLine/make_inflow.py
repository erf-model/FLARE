#!/usr/bin/env python3
"""Write the inflow profile and input sounding of the fire-under-a-line benchmark.

A neutral log law up to a capping inversion, held constant above it, with theta rising at a fixed lapse
rate above the inversion's base:

    u(z) = u*/kappa ln((min(z, z_i) + z0) / z0),   u* = kappa U_ref / ln((z_ref + z0) / z0)

Outputs (in the current directory unless --out is given):
  inflow_profile   the x-low Dirichlet profile, columns z u v w theta
  input_sounding   the same profile as ERF's input sounding (surface pressure in hPa)

The rows are the ground (u = 0), every cell centre (dz/2, 3 dz/2, ...) and the top. ERF interpolates both files
linearly in z, so a row at every cell centre gives each cell the log law's value there; rows at the cell faces
would give the first cell the mean of 0 and the 10 m value (3.0 m/s instead of 5.2 m/s).

The inflow turbulent kinetic energy printed for xlo.KE is the neutral wall value of the k-equation closure
(Axell and Liungman 2001, eq 16; Docs/sphinx_doc/theory/RANS.rst), k = u*^2 / (c_mu0)^2 with ERF's
c_mu0 = 0.5562 (Source/DataStructs/ERF_TurbStruct.H).

The defaults give the committed files: python3 make_inflow.py
"""
import argparse
import math
import os

import numpy as np

KAPPA = 0.41
CMU0 = 0.5562   # ERF_TurbStruct.H


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--u_ref", type=float, default=6.0, help="wind speed at z_ref (m/s)")
    ap.add_argument("--z_ref", type=float, default=10.0, help="reference height (m)")
    ap.add_argument("--z0", type=float, default=0.05, help="roughness length (m); must match erf.most.z0")
    ap.add_argument("--inversion", type=float, default=400.0, help="base of the capping inversion (m)")
    ap.add_argument("--lapse", type=float, default=0.01, help="theta gradient above the inversion (K/m)")
    ap.add_argument("--top", type=float, default=600.0, help="domain top (m)")
    ap.add_argument("--dz", type=float, default=10.0, help="the mesh's cell height (m): rows at the cell centres")
    ap.add_argument("--out", default=".", help="output directory")
    a = ap.parse_args()
    for name in ("u_ref", "z_ref", "z0", "inversion", "top", "dz"):
        if getattr(a, name) <= 0.0:
            raise SystemExit(f"--{name} must be > 0")

    ustar = a.u_ref * KAPPA / math.log((a.z_ref + a.z0) / a.z0)
    zs = np.concatenate(([0.0], np.arange(0.5 * a.dz, a.top, a.dz), [a.top]))
    u = ustar / KAPPA * np.log((np.minimum(zs, a.inversion) + a.z0) / a.z0)
    th = 300.0 + a.lapse * np.maximum(zs - a.inversion, 0.0)
    with open(os.path.join(a.out, "inflow_profile"), "w") as f:
        for zz, uu, tt in zip(zs, u, th):
            f.write(f"{zz:.2f} {uu:.6f} 0.0 0.0 {tt:.4f}\n")
    with open(os.path.join(a.out, "input_sounding"), "w") as f:
        f.write("1000.0 300.0 0.0\n")
        for zz, uu, tt in zip(zs, u, th):
            f.write(f"{zz:.2f} {tt:.4f} 0.0 {uu:.6f} 0.0\n")
    print(f"u* = {ustar:.4f} m/s, inflow KE u*^2/c_mu0^2 = {ustar ** 2 / CMU0 ** 2:.4f} m^2/s^2, "
          f"max_geom_lscale 0.1 kappa z_i = {0.1 * KAPPA * a.inversion:.1f} m")


if __name__ == "__main__":
    main()
