#!/usr/bin/env python3
"""Convert FARSITE gridded winds (erf.fire.farsite_winds, or WindNinja) to
viewing formats: VTK, shapefile, KMZ, PDF and GeoPDF.

    python3 farsite_winds_convert.py winds/flare.atm --to vtk shp kmz
    python3 farsite_winds_convert.py winds/flare_000061_vel.asc --units metric --to pdf
    python3 farsite_winds_convert.py --self-test

Input is either a .atm file (every grid pair it lists is converted; the units
come from its ENGLISH/METRIC line) or one _vel.asc file (its _ang.asc partner
is found by name; give --units). The grids are the ESRI ASCII layout WindNinja
writes: speed in mph (ENGLISH) or km/h (METRIC), direction in degrees the wind
blows FROM, clockwise from north, north row first, NODATA_value cells skipped.

Outputs, one per grid pair, named after the _vel.asc without "_vel.asc":

  vtk     legacy ASCII STRUCTURED_POINTS (.vtk): speed (input unit),
          direction (deg) and the vector (u, v, 0) in m/s. Standard library.
  shp     wind arrows as polylines (.shp/.shx/.dbf, plus the .prj copied):
          fields SPEED (input unit) and DIR (deg). Needs pyshp.
  kmz     Google Earth: arrows as coloured line strings in lon/lat. Needs
          pyproj and a projection (the grid's .prj, or --prj / --epsg).
  pdf     a map of the speed with arrows (not georeferenced). Needs matplotlib.
  geopdf  a geospatial PDF of the speed raster through GDAL's PDF driver.
          Needs the GDAL Python bindings (osgeo) and a projection.

--stride N keeps every Nth cell in x and y for the arrow outputs (shp, kmz,
pdf); the VTK and GeoPDF outputs keep every cell.
"""
import argparse
import math
import os
import re
import struct
import sys
import tempfile
import zipfile

MPS_PER = {"english": 1609.344 / 3600.0, "metric": 1.0 / 3.6}
UNIT_NAME = {"english": "mph", "metric": "km/h"}


# --------------------------------------------------------------- reading
def read_asc(path):
    """(header dict, rows north-first as lists of floats)."""
    with open(path) as f:
        text = f.read().split("\n")
    hdr = {}
    i = 0
    while i < len(text) and len(hdr) < 6:
        parts = text[i].split()
        if len(parts) == 2 and re.match(r"[A-Za-z_]+$", parts[0]):
            hdr[parts[0].lower()] = float(parts[1])
            i += 1
        else:
            break
    for k in ("ncols", "nrows", "cellsize"):
        if k not in hdr:
            sys.exit(f"{path}: header has no {k}")
    if "xllcorner" not in hdr:
        if "xllcenter" in hdr:
            hdr["xllcorner"] = hdr["xllcenter"] - 0.5 * hdr["cellsize"]
            hdr["yllcorner"] = hdr["yllcenter"] - 0.5 * hdr["cellsize"]
        else:
            sys.exit(f"{path}: header has no xllcorner")
    hdr.setdefault("nodata_value", -9999.0)
    vals = [float(x) for ln in text[i:] for x in ln.split()]
    nc, nr = int(hdr["ncols"]), int(hdr["nrows"])
    if len(vals) != nc * nr:
        sys.exit(f"{path}: {len(vals)} values for {nc} x {nr} cells")
    rows = [vals[r * nc:(r + 1) * nc] for r in range(nr)]
    return hdr, rows


def read_atm(path):
    """(units, [(vel path, ang path, 'MM DD HHMM'), ...])."""
    with open(path) as f:
        lines = [ln.split() for ln in f if ln.strip()]
    if not lines or lines[0][0].upper() not in ("WINDS", "WINDS_AND_CLOUDS"):
        sys.exit(f"{path}: not a FARSITE .atm file (first line {lines[0] if lines else 'empty'})")
    if len(lines) < 2:
        sys.exit(f"{path}: no ENGLISH/METRIC line")
    unit = lines[1][0].upper()
    if unit not in ("ENGLISH", "METRIC"):
        sys.exit(f"{path}: second line must be ENGLISH or METRIC, got {unit}")
    base = os.path.dirname(os.path.abspath(path))
    pairs = []
    for ln in lines[2:]:
        if len(ln) < 5:
            sys.exit(f"{path}: line '{' '.join(ln)}' needs MM DD HHMM vel ang")
        pairs.append((os.path.join(base, ln[3]), os.path.join(base, ln[4]), " ".join(ln[:3])))
    if not pairs:
        sys.exit(f"{path}: no grids listed")
    return unit.lower(), pairs


def partner_ang(vel):
    if not vel.endswith("_vel.asc"):
        sys.exit(f"{vel}: expected a name ending in _vel.asc")
    return vel[:-len("_vel.asc")] + "_ang.asc"


def find_prj(vel, override):
    if override:
        with open(override) as f:
            return f.read()
    p = vel[:-4] + ".prj"
    if os.path.isfile(p):
        with open(p) as f:
            return f.read()
    return None


class WindGrid:
    """Speed and direction on one grid; cell (i, j) from the south-west corner."""

    def __init__(self, vel, ang, units, label=""):
        hv, sv = read_asc(vel)
        ha, sa = read_asc(ang)
        for k in ("ncols", "nrows", "xllcorner", "yllcorner", "cellsize"):
            if abs(hv[k] - ha[k]) > 1e-6 * max(1.0, abs(hv[k])):
                sys.exit(f"{vel} and {ang} differ in {k}")
        self.nx, self.ny = int(hv["ncols"]), int(hv["nrows"])
        self.x0, self.y0, self.dx = hv["xllcorner"], hv["yllcorner"], hv["cellsize"]
        self.nodata = hv["nodata_value"]
        self.units = units
        self.label = label
        self.root = vel[:-len("_vel.asc")] if vel.endswith("_vel.asc") else vel[:-4]
        # store south row first
        self.spd = sv[::-1]
        self.dir = sa[::-1]

    def valid(self, i, j):
        return self.spd[j][i] != self.nodata and self.dir[j][i] != self.nodata

    def uv(self, i, j):
        """(u east, v north) in m/s; the direction is where the wind comes FROM."""
        s = self.spd[j][i] * MPS_PER[self.units]
        d = math.radians(self.dir[j][i])
        return -s * math.sin(d), -s * math.cos(d)

    def centre(self, i, j):
        return self.x0 + (i + 0.5) * self.dx, self.y0 + (j + 0.5) * self.dx

    def arrows(self, stride):
        """(x0, y0, x1, y1, xl, yl, xr, yr, speed, dir) per kept valid cell: a
        shaft of 0.8 stride cells at full speed scaled by the grid maximum, and
        two head barbs."""
        smax = max((self.spd[j][i] for j in range(self.ny) for i in range(self.nx)
                    if self.valid(i, j)), default=0.0)
        out = []
        for j in range(stride // 2, self.ny, stride):
            for i in range(stride // 2, self.nx, stride):
                if not self.valid(i, j):
                    continue
                u, v = self.uv(i, j)
                m = math.hypot(u, v)
                xc, yc = self.centre(i, j)
                if m == 0.0 or smax == 0.0:
                    continue
                L = 0.8 * stride * self.dx * self.spd[j][i] / smax
                ex, ey = u / m, v / m
                x0, y0 = xc - 0.5 * L * ex, yc - 0.5 * L * ey
                x1, y1 = xc + 0.5 * L * ex, yc + 0.5 * L * ey
                h = 0.3 * L
                c, s = math.cos(math.radians(25)), math.sin(math.radians(25))
                xl = x1 - h * (ex * c - ey * s); yl = y1 - h * (ey * c + ex * s)
                xr = x1 - h * (ex * c + ey * s); yr = y1 - h * (ey * c - ex * s)
                out.append((x0, y0, x1, y1, xl, yl, xr, yr, self.spd[j][i], self.dir[j][i]))
        return out


# --------------------------------------------------------------- writers
def write_vtk(g, out):
    path = out + ".vtk"
    n = g.nx * g.ny
    with open(path, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write(f"FARSITE winds {g.label} ({UNIT_NAME[g.units]})\nASCII\nDATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {g.nx} {g.ny} 1\n")
        f.write(f"ORIGIN {g.x0 + 0.5 * g.dx:.6f} {g.y0 + 0.5 * g.dx:.6f} 0\n")
        f.write(f"SPACING {g.dx:.6f} {g.dx:.6f} 1\n")
        f.write(f"POINT_DATA {n}\n")
        f.write(f"SCALARS speed_{UNIT_NAME[g.units].replace('/', 'p')} double 1\nLOOKUP_TABLE default\n")
        for j in range(g.ny):
            f.write(" ".join(f"{g.spd[j][i]:.6g}" if g.valid(i, j) else "nan" for i in range(g.nx)) + "\n")
        f.write("SCALARS direction_from_deg double 1\nLOOKUP_TABLE default\n")
        for j in range(g.ny):
            f.write(" ".join(f"{g.dir[j][i]:.6g}" if g.valid(i, j) else "nan" for i in range(g.nx)) + "\n")
        f.write("VECTORS wind_mps double\n")
        for j in range(g.ny):
            for i in range(g.nx):
                if g.valid(i, j):
                    u, v = g.uv(i, j)
                    f.write(f"{u:.6g} {v:.6g} 0\n")
                else:
                    f.write("nan nan 0\n")
    return path


def write_shp(g, out, stride, prj):
    try:
        import shapefile
    except ImportError:
        sys.exit("shp output needs pyshp: pip install pyshp")
    w = shapefile.Writer(out, shapeType=shapefile.POLYLINE)
    w.field("SPEED", "N", decimal=2)
    w.field("DIR", "N", decimal=0)
    for (x0, y0, x1, y1, xl, yl, xr, yr, s, d) in g.arrows(stride):
        w.line([[[x0, y0], [x1, y1]], [[xl, yl], [x1, y1], [xr, yr]]])
        w.record(round(s, 2), round(d))
    w.close()
    if prj:
        with open(out + ".prj", "w") as f:
            f.write(prj)
    return out + ".shp"


def _to_lonlat(prj, epsg):
    try:
        from pyproj import CRS, Transformer
    except ImportError:
        sys.exit("kmz output needs pyproj: pip install pyproj")
    if epsg:
        src = CRS.from_epsg(int(epsg))
    elif prj:
        src = CRS.from_wkt(prj)
    else:
        sys.exit("kmz output needs a projection: a .prj beside the grid, --prj FILE or --epsg CODE")
    t = Transformer.from_crs(src, CRS.from_epsg(4326), always_xy=True)
    return lambda x, y: t.transform(x, y)


def _kml_colour(s, smax):
    """aabbggrr from blue (calm) through green to red (the grid maximum)."""
    f = 0.0 if smax <= 0 else min(1.0, s / smax)
    r = int(255 * min(1.0, 2 * f)); b = int(255 * min(1.0, 2 * (1 - f))); gch = int(255 * (1 - abs(2 * f - 1)))
    return f"ff{b:02x}{gch:02x}{r:02x}"


def write_kmz(g, out, stride, prj, epsg):
    ll = _to_lonlat(prj, epsg)
    arr = g.arrows(stride)
    smax = max((a[8] for a in arr), default=0.0)
    unit = UNIT_NAME[g.units]
    kml = ['<?xml version="1.0" encoding="UTF-8"?>',
           '<kml xmlns="http://www.opengis.net/kml/2.2"><Document>',
           f"<name>{os.path.basename(out)}</name>",
           f"<description>FARSITE winds {g.label}, speed in {unit}, direction from (deg)</description>"]
    for (x0, y0, x1, y1, xl, yl, xr, yr, s, d) in arr:
        pts = [ll(x0, y0), ll(x1, y1), ll(xl, yl), ll(x1, y1), ll(xr, yr)]
        coords = " ".join(f"{lon:.7f},{lat:.7f},0" for lon, lat in pts)
        kml.append(f"<Placemark><name>{s:.1f} {unit}</name><description>from {d:.0f} deg</description>"
                   f"<Style><LineStyle><color>{_kml_colour(s, smax)}</color><width>2</width></LineStyle></Style>"
                   f"<LineString><tessellate>1</tessellate><coordinates>{coords}</coordinates></LineString></Placemark>")
    kml.append("</Document></kml>")
    path = out + ".kmz"
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("doc.kml", "\n".join(kml))
    return path


def write_pdf(g, out, stride):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        sys.exit("pdf output needs matplotlib and numpy")
    S = np.array([[g.spd[j][i] if g.valid(i, j) else np.nan for i in range(g.nx)] for j in range(g.ny)])
    ext = [g.x0, g.x0 + g.nx * g.dx, g.y0, g.y0 + g.ny * g.dx]
    fig, ax = plt.subplots(figsize=(8, 8 * g.ny / max(g.nx, 1) + 1))
    im = ax.imshow(S, origin="lower", extent=ext, cmap="viridis")
    fig.colorbar(im, ax=ax, label=f"speed [{UNIT_NAME[g.units]}]")
    js = range(stride // 2, g.ny, stride)
    is_ = range(stride // 2, g.nx, stride)
    X, Y, U, V = [], [], [], []
    for j in js:
        for i in is_:
            if g.valid(i, j):
                x, y = g.centre(i, j); u, v = g.uv(i, j)
                X.append(x); Y.append(y); U.append(u); V.append(v)
    if X:
        ax.quiver(X, Y, U, V, color="w")
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    ax.set_title(f"FARSITE winds {g.label}")
    path = out + ".pdf"
    fig.savefig(path)
    plt.close(fig)
    return path


def write_geopdf(g, out, prj, epsg):
    try:
        from osgeo import gdal, osr
    except ImportError:
        sys.exit("geopdf output needs the GDAL Python bindings (osgeo): conda install gdal")
    srs = osr.SpatialReference()
    if epsg:
        srs.ImportFromEPSG(int(epsg))
    elif prj:
        srs.ImportFromWkt(prj)
    else:
        sys.exit("geopdf output needs a projection: a .prj beside the grid, --prj FILE or --epsg CODE")
    mem = gdal.GetDriverByName("MEM").Create("", g.nx, g.ny, 1, gdal.GDT_Float32)
    mem.SetGeoTransform([g.x0, g.dx, 0.0, g.y0 + g.ny * g.dx, 0.0, -g.dx])
    mem.SetProjection(srs.ExportToWkt())
    band = mem.GetRasterBand(1)
    band.SetNoDataValue(g.nodata)
    band.WriteRaster(0, 0, g.nx, g.ny,
                     struct.pack(f"{g.nx * g.ny}f", *[g.spd[j][i] for j in range(g.ny - 1, -1, -1)
                                                      for i in range(g.nx)]))
    path = out + "_speed.pdf"
    gdal.Translate(path, mem, format="PDF", outputType=gdal.GDT_Byte, scaleParams=[[]])
    return path


# --------------------------------------------------------------- driver
def convert(inp, formats, out_dir=None, units=None, stride=1, prj_file=None, epsg=None):
    if inp.lower().endswith(".atm"):
        u, pairs = read_atm(inp)
    else:
        if units is None:
            sys.exit("a single _vel.asc needs --units english|metric")
        u, pairs = units, [(inp, partner_ang(inp), "")]
    written = []
    for vel, ang, when in pairs:
        g = WindGrid(vel, ang, u, when)
        prj = find_prj(vel, prj_file)
        out = os.path.join(out_dir, os.path.basename(g.root)) if out_dir else g.root
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        for fmt in formats:
            if fmt == "vtk":
                written.append(write_vtk(g, out))
            elif fmt == "shp":
                written.append(write_shp(g, out, stride, prj))
            elif fmt == "kmz":
                written.append(write_kmz(g, out, stride, prj, epsg))
            elif fmt == "pdf":
                written.append(write_pdf(g, out, stride))
            elif fmt == "geopdf":
                written.append(write_geopdf(g, out, prj, epsg))
    return written


def self_test():
    """Round-trip a 3 x 2 grid with known winds through every format whose
    library is installed; vtk always. Exit 1 on a mismatch."""
    fails = []
    with tempfile.TemporaryDirectory() as d:
        # u = 3, v = 4 m/s (from 216.87 deg) everywhere but one NODATA cell
        mph = 5.0 / MPS_PER["english"]
        hdr = "ncols\t3\nnrows\t2\nxllcorner\t500000.0\nyllcorner\t4400000.0\ncellsize\t30.0\nNODATA_value\t-9999\n"
        with open(os.path.join(d, "t_vel.asc"), "w") as f:
            f.write(hdr + f"{mph:.2f}\t{mph:.2f}\t-9999\t\n{mph:.2f}\t{2 * mph:.2f}\t{mph:.2f}\t\n")
        with open(os.path.join(d, "t_ang.asc"), "w") as f:
            f.write(hdr + "217\t217\t-9999\t\n217\t90\t217\t\n")
        with open(os.path.join(d, "t.atm"), "w") as f:
            f.write("WINDS\nENGLISH\n12 30 1100 t_vel.asc t_ang.asc\n")
        with open(os.path.join(d, "t_vel.prj"), "w") as f:
            f.write('PROJCS["WGS_1984_UTM_Zone_13N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",'
                    'SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],'
                    'UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],'
                    'PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],'
                    'PARAMETER["Central_Meridian",-105.0],PARAMETER["Scale_Factor",0.9996],'
                    'PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]')
        g = WindGrid(os.path.join(d, "t_vel.asc"), os.path.join(d, "t_ang.asc"), "english")
        # south row first: (1, 0) is the 2x wind from the east (90 deg) -> u = -10, v = 0
        u, v = g.uv(1, 0)
        if abs(u + 10.0) > 0.01 or abs(v) > 0.01:
            fails.append(f"uv(1,0) = {u}, {v}, expected -10, 0 (row order or direction)")
        u, v = g.uv(0, 1)
        if abs(u - 3.0) > 0.01 or abs(v - 4.0) > 0.01:
            fails.append(f"uv(0,1) = {u}, {v}, expected 3, 4")
        if g.valid(2, 1):
            fails.append("NODATA cell (2, 1) read as valid")
        formats = ["vtk"]
        for mod, fmt in (("shapefile", "shp"), ("pyproj", "kmz"), ("matplotlib", "pdf"), ("osgeo", "geopdf")):
            try:
                __import__(mod)
                formats.append(fmt)
            except ImportError:
                print(f"  skip {fmt}: {mod} not installed")
        out = convert(os.path.join(d, "t.atm"), formats)
        vtk = open(os.path.join(d, "t.vtk")).read().split("VECTORS wind_mps double\n")[1].split("\n")
        if [round(float(x), 2) for x in vtk[1].split()[:2]] != [-10.0, 0.0] or vtk[5].split()[0] != "nan":
            fails.append(f"vtk vectors {vtk[:6]}")
        if "shp" in formats:
            import shapefile
            r = shapefile.Reader(os.path.join(d, "t"))
            recs = [tuple(rec) for rec in r.records()]
            if len(recs) != 5 or (round(2 * mph, 2), 90) not in [(float(a), int(b)) for a, b in recs]:
                fails.append(f"shp records {recs}")
            if not os.path.isfile(os.path.join(d, "t.prj")):
                fails.append("shp .prj not written")
        if "kmz" in formats:
            with zipfile.ZipFile(os.path.join(d, "t.kmz")) as z:
                k = z.read("doc.kml").decode()
            lons = [float(c.split(",")[0]) for c in re.findall(r"(-?\d+\.\d+,-?\d+\.\d+),0", k)]
            if k.count("<Placemark>") != 5 or not all(-105.01 < x < -104.99 for x in lons):
                fails.append(f"kmz: {k.count('<Placemark>')} placemarks, lon range "
                             f"{min(lons) if lons else None}..{max(lons) if lons else None}")
        print("  wrote " + ", ".join(os.path.basename(p) for p in out))
    for f in fails:
        print("  FAIL " + f)
    print("farsite_winds_convert self-test " + ("FAILED" if fails else "passed"))
    return 1 if fails else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", nargs="?", help=".atm file or one _vel.asc")
    ap.add_argument("--to", nargs="+", choices=["vtk", "shp", "kmz", "pdf", "geopdf"], default=["vtk"])
    ap.add_argument("--out", help="output directory (default: beside each grid)")
    ap.add_argument("--units", choices=["english", "metric"], help="units of a single _vel.asc")
    ap.add_argument("--stride", type=int, default=1, help="keep every Nth cell for arrows (>= 1)")
    ap.add_argument("--prj", help="projection (.prj WKT) instead of the grid's own")
    ap.add_argument("--epsg", help="EPSG code instead of a .prj")
    ap.add_argument("--self-test", action="store_true", help="round-trip a synthetic grid and exit")
    a = ap.parse_args()
    if a.self_test:
        sys.exit(self_test())
    if not a.input:
        ap.error("give an input file, or --self-test")
    if a.stride < 1:
        ap.error("--stride must be >= 1")
    for p in convert(a.input, a.to, a.out, a.units, a.stride, a.prj, a.epsg):
        print(p)


if __name__ == "__main__":
    main()
