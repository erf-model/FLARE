#!/usr/bin/env python3
"""Mean-tension check of the single span of Ma, Bocchini and Christou 2020 (Struct. Saf.
82:101890, Section 3) through the conductor model: Drake on 300 m at a sag of 1/50, in still air
with a uniform wind normal to the span handed to every node (make_case.py --mode prescribed).

    python3 run_check.py --exe path/to/flare --out check_dir [--launcher "mpiexec -n 1"] [--winds 20 60 100] [--without-fft]

For each wind the run's end tensions, averaged over its last 20 s of 120 s, must lie

  * within --tol-moordyn percent (default 0.01) of the lumped chain at rest with MoorDyn's
    conventions (static_reference.chain(moordyn=True): drag on the unstretched length, weight less
    the displaced air's), which the run should reproduce to its convergence; the runs agree to
    0.0005 %. The chain's values are pinned below for 20, 60 and 100 m/s; with numpy present they
    are recomputed and the pins are checked against them to 0.001 %.
  * within --tol-catenary percent (default 1.5) of the elastic catenary under the resultant of
    the weight and the full drag (static_reference.catenary, standard library only), which takes
    the whole wind as normal to the cable; the chain lies 0.03 % (20 m/s) to 1.0 % (100 m/s)
    below it (static_reference.py's docstring) and the run, with MoorDyn's conventions, 0.07 % to
    1.28 % below it.

and the swing of the middle of the span within 0.05 degrees of atan(f / w), f = 0.5 rho C_D D U^2,
w the weight less the displaced air's (MoorDyn's). Prints "SingleSpanMa2020: all checks hold" and
returns 0 when every check holds, 1 otherwise. The runs go in --out/u<wind>; they are removed after a
pass and kept after a failure.
Needs a FLARE built with the real MoorDyn-C (the stub has no dynamics). Standard library only.
"""
import argparse
import importlib.util
import math
import os
import shlex
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("static_reference", os.path.join(HERE, "static_reference.py"))
ref = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ref)

SPAN, SEGMENTS, AVERAGE = 300.0, 80, 20.0
# static_reference.chain(300, L0, 15.95, U, 3e7, 80, moordyn=True)[0] (N), L0 = 300.0204483 m
CHAIN_MOORDYN = {20.0: 31766.49, 60.0: 76506.12, 100.0: 148786.03}


def rows(path):
    with open(path) as f:
        head = f.readline().split()
        return head, [list(map(float, ln.split())) for ln in f if ln.strip()]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, help="the FLARE executable")
    ap.add_argument("--out", required=True, help="directory for the runs: new, empty, or holding only u<wind> run directories")
    ap.add_argument("--launcher", default="", help='prefix of every run, e.g. "mpiexec -n 1"')
    ap.add_argument("--winds", type=float, nargs="+", default=sorted(CHAIN_MOORDYN),
                    help="uniform winds (m/s), each one of %s" % sorted(CHAIN_MOORDYN))
    ap.add_argument("--tol-moordyn", type=float, default=0.01, help="percent, in (0, 1]")
    ap.add_argument("--tol-catenary", type=float, default=1.5, help="percent, in (0, 5]")
    ap.add_argument("--without-fft", action="store_true", help="erf.use_fft = false, for an ERF built without FFT")
    a = ap.parse_args(argv)
    if any(u not in CHAIN_MOORDYN for u in a.winds):
        ap.error("--winds: each must be one of %s (the pinned chain's)" % sorted(CHAIN_MOORDYN))
    if not (0.0 < a.tol_moordyn <= 1.0 and 0.0 < a.tol_catenary <= 5.0):
        ap.error("--tol-moordyn must be in (0, 1] and --tol-catenary in (0, 5]")
    exe = shutil.which(a.exe) or a.exe
    if not os.path.isfile(exe):
        ap.error("--exe: %s is not a file" % a.exe)
    exe = os.path.abspath(exe)      # the runs start in their own directories
    if os.path.isdir(a.out):
        stray = [n for n in os.listdir(a.out) if not (n.startswith("u") and os.path.isdir(os.path.join(a.out, n)))]
        if stray:
            ap.error("--out: %s holds %s besides u<wind> run directories; give a new or empty directory" % (a.out, ", ".join(sorted(stray)[:3])))
    L0 = ref.rest_length(SPAN, ref.SAG_RATIO, ref.DRAKE_WEIGHT, ref.EA_DEFAULT, SEGMENTS)
    try:
        import numpy  # noqa: F401
        have_numpy = True
    except ImportError:
        have_numpy = False
        print("SingleSpanMa2020: numpy not found; the chain's pinned values are used as they are")
    for n in (os.listdir(a.out) if os.path.isdir(a.out) else []):
        shutil.rmtree(os.path.join(a.out, n))
    ok = True
    w_wet = ref.DRAKE_WEIGHT - ref.RHO * ref.G * math.pi * ref.DIAMETER ** 2 / 4.0   # MoorDyn: weight less the displaced air's
    print("%6s %12s %12s %9s %12s %9s %9s %9s" % ("U", "T_run(N)", "T_chain(N)", "diff%", "T_cat(N)", "diff%", "swing", "atan(f/w)"))
    for U in a.winds:
        pin = CHAIN_MOORDYN[U]
        if have_numpy:
            live = ref.chain(SPAN, L0, ref.DRAKE_WEIGHT, U, ref.EA_DEFAULT, SEGMENTS, moordyn=True)[0]
            if abs(live - pin) > 1.0e-5 * pin:
                print("SingleSpanMa2020: FAIL the pinned chain tension at %g m/s, %.2f N, is not the chain's %.2f N" % (U, pin, live))
                ok = False
        d = os.path.join(a.out, "u%g" % U)
        subprocess.run([sys.executable, "-I", os.path.join(HERE, "make_case.py"), "--out", d, "--mode", "prescribed",
                        "--span", str(SPAN), "--wind", str(U), "--segments", str(SEGMENTS)] + (["--without-fft"] if a.without_fft else []),
                       check=True, stdout=subprocess.DEVNULL)
        with open(os.path.join(d, "run.log"), "w") as log:
            rc = subprocess.run(shlex.split(a.launcher) + [exe, "inputs", "amrex.call_addr2line=0"], cwd=d,
                                stdout=log, stderr=subprocess.STDOUT).returncode
        if rc != 0:
            print("SingleSpanMa2020: FAIL the run at %g m/s ended with return code %d (see %s)" % (U, rc, os.path.join(d, "run.log")))
            ok = False
            continue
        try:
            head, data = rows(os.path.join(d, "conductors", "S.dat"))
            it, ia, ib, isw = (head.index(k) for k in ("time", "tension_a", "tension_b", "swing_deg"))
            if not data:
                raise ValueError("no rows")
        except (OSError, ValueError) as e:
            print("SingleSpanMa2020: FAIL at %g m/s the run wrote no readable conductors/S.dat (%s)" % (U, e))
            ok = False
            continue
        end = data[-1][it]
        last = [r for r in data if r[it] >= end - AVERAGE]
        T = sum(0.5 * (r[ia] + r[ib]) for r in last) / len(last)
        swing = abs(sum(r[isw] for r in last) / len(last))
        f = 0.5 * ref.RHO * ref.CD * ref.DIAMETER * U * U
        Tc = ref.catenary(SPAN, L0, ref.DRAKE_WEIGHT, f, ref.EA_DEFAULT, SEGMENTS)[0]
        expect_swing = math.degrees(math.atan2(f, w_wet))
        dm, dc = 100.0 * (T - pin) / pin, 100.0 * (T - Tc) / Tc
        print("%6g %12.2f %12.2f %9.4f %12.2f %9.4f %9.3f %9.3f" % (U, T, pin, dm, Tc, dc, swing, expect_swing))
        if abs(dm) > a.tol_moordyn:
            print("SingleSpanMa2020: FAIL at %g m/s the run's tension is %.4f %% from the MoorDyn-convention chain (tolerance %g %%)" % (U, dm, a.tol_moordyn))
            ok = False
        if abs(dc) > a.tol_catenary:
            print("SingleSpanMa2020: FAIL at %g m/s the run's tension is %.4f %% from the catenary (tolerance %g %%)" % (U, dc, a.tol_catenary))
            ok = False
        if abs(swing - expect_swing) > 0.05:
            print("SingleSpanMa2020: FAIL at %g m/s the swing is %.3f deg, atan(f/w) %.3f deg" % (U, swing, expect_swing))
            ok = False
    if ok:
        print("SingleSpanMa2020: all checks hold")
        for U in a.winds:
            shutil.rmtree(os.path.join(a.out, "u%g" % U), ignore_errors=True)
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
