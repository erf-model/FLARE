#!/usr/bin/env python3
"""Palisades Fire and two transmission lines: prepare the run directory, or post-process a finished run.

    python3 case.py prepare [--run DIR]                 # terrain, LANDFIRE fuel map, decks and small files into DIR
    python3 case.py post --basemap TIF [--run DIR] [--figs DIR]   # maps, movie, line temperatures, fragility, head distance

DIR defaults to ./run. 'prepare' copies this folder's decks and small files into DIR, runs gen_palisades.py in the Real_Terrain folder when its terrain files are missing, fetches the
LANDFIRE LF2024 FBFM40 window and maps it at 30 m (make_landfire_fuel_map.py), and checks that the
conductor block of the deck is the one fire_lines_towers.py builds from lines_config.json and osm_lines.json over
the terrain (exit 1 if not). 'post' runs fire_lines_plot.py, fire_lines_fragility.py and fire_head_distance.py.
The basemap is any RGB GeoTIFF over the domain, for example the one make_map_movie.py --fetch usgs writes in
Exec/CanonicalTests/Fire/Real_Terrain/Palisades_Fire.
"""
import argparse, os, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
FLARE = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
TOOLS = os.path.join(FLARE, "Exec", "Tools")
RT = os.path.join(FLARE, "Exec", "CanonicalTests", "Fire", "Real_Terrain", "Palisades_Fire")
CASE, TITLE = "palisades", "Palisades Fire, 7 Jan 2025"
EDITION, DX = "LF2024", "30"
TERRAIN = ["terrain_palisades.txt", "terrain_palisades_fire.txt"]              # files the decks read, from the Real_Terrain folder
SMALL = ["inputs", "inflow_hrrr.txt", "sounding_hrrr.txt", "ignitions_spinup.csv", "lines_meta.json"]
CROP, TRANGE = "5.0,13.5,4.8,13.3", "40,115"


def sh(*cmd, cwd):
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("step", choices=["prepare", "post"])
ap.add_argument("--run", default=os.path.join(HERE, "run"))
ap.add_argument("--basemap", default="")
ap.add_argument("--figs", default="")
a = ap.parse_args()
run = os.path.abspath(a.run)
dom = os.path.join(RT, f"{CASE}_domain.json")

if a.step == "prepare":
    os.makedirs(run, exist_ok=True)
    if not all(os.path.exists(os.path.join(RT, t)) for t in TERRAIN):
        sh(sys.executable, "gen_palisades.py", cwd=RT)
    for t in TERRAIN:
        if not os.path.exists(os.path.join(RT, t)):
            sys.exit(f"missing {os.path.join(RT, t)}")
        shutil.copy(os.path.join(RT, t), run)
    for f in SMALL:
        shutil.copy(os.path.join(HERE, f), run)
    shutil.copy(dom, run)
    fuel = f"fuel_{CASE}_fbfm40.asc"
    if not os.path.exists(os.path.join(run, fuel)):
        sh(sys.executable, os.path.join(TOOLS, "make_landfire_fuel_map.py"), "--fetch", EDITION, "--dx", DX,
           "--domain", os.path.basename(dom), cwd=run)
    sh(sys.executable, os.path.join(TOOLS, "fire_lines_towers.py"), os.path.join(HERE, "lines_config.json"),
       "--osm", os.path.join(HERE, "osm_lines.json"), "--terrain", TERRAIN[0],
       "--meta", "lines_meta_check.json", "--check", "inputs", cwd=run)
    print(f"ready: cd {run} && mpiexec -n 4 <FLARE build>/flare inputs")
else:
    if not a.basemap:
        sys.exit("post needs --basemap <RGB GeoTIFF over the domain>")
    figs = os.path.abspath(a.figs or os.path.join(run, "figs"))
    os.makedirs(figs, exist_ok=True)
    sh(sys.executable, os.path.join(TOOLS, "fire_lines_plot.py"), CASE, run, os.path.abspath(a.basemap), figs,
       "--domain", dom, "--title", TITLE, "--every", "4", "--crop", CROP, "--trange", TRANGE, cwd=run)
    sh(sys.executable, os.path.join(TOOLS, "fire_lines_fragility.py"), figs, f"{CASE}={run}", cwd=run)
    sh(sys.executable, os.path.join(TOOLS, "fire_head_distance.py"), figs, f"{CASE}={run}", cwd=run)
