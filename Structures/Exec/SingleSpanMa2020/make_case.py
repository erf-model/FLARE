#!/usr/bin/env python3
"""Write the deck of the single-span case of Ma, Bocchini and Christou 2020 (Struct. Saf.
82:101890, Section 3): one ACSR Drake span on hinged supports at the same level over flat open
terrain, the wind normal to the span.

    python3 make_case.py --out run_dir --mode prescribed --wind 40      # still air, uniform wind
    python3 make_case.py --out run_dir --mode flow --wind 40 --seed 1   # ERF flow, random gusts

The span runs along +y and the wind blows along +x. Its unstretched length gives a sag of
span / 50 under its weight alone (static_reference.rest_length). Ma et al. take the wind at
40 m; we read that as the reference height of the span, two thirds of the sag below the supports,
so the supports stand at 40 m + (2/3) span / 50.

--mode prescribed  The air is at rest and the line takes --wind (m/s) as a uniform wind at every
                   node (conductors.prescribed_velocity): the steady case that static_reference.py
                   solves. Writes inputs.
--mode flow        A neutral boundary layer over flat ground with roughness z0 (Ma's 0.05 m), k-equation
                   RANS, entering through x-low as the log law with --wind at 10 m and leaving
                   through x-high; the line takes the flow's wind with random gusts
                   (conductors.gust_type = random, conductors.gust_seed = --seed). The wind at 10 m
                   averaged over --average seconds is recorded beside the conductor's largest
                   tension (conductors.intensity_measure = wind). Ma et al.'s intensity measure is
                   V_ms, the highest 1-min mean wind at 10 m; a 10-min mean is lower than the
                   highest 1-min mean in the same 10 min by a gust factor, which is not applied here.
                   Writes inputs, inflow_profile and input_sounding: the base deck of the fragility
                   ensemble driver.

The conductor's tensions are recorded against Drake's rated strength (140.18 kN, Ma 2020 Table 3);
nothing lets go (failure_break and failure_stop are off).
"""
import argparse
import importlib.util
import math
import os
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("static_reference", os.path.join(HERE, "static_reference.py"))
ref = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ref)

KAPPA = 0.4
DRAKE_RATED = 140180.0        # N, Ma 2020 Table 3
REFERENCE_HEIGHT = 40.0       # m, Ma 2020 Section 3 ("conductor height z = 40 m")
THETA = 300.0                 # K, neutral
KE_OVER_USTAR2 = 3.3          # inflow k = 3.3 u*^2, as Structures/Exec/PowerLines/flow.inputs


def parse(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--mode", choices=("prescribed", "flow"), required=True)
    ap.add_argument("--span", type=float, default=300.0, help="span (m), in [50, 1000]")
    ap.add_argument("--wind", type=float, required=True,
                    help="prescribed: the uniform wind at the line (m/s); flow: the inflow's wind at 10 m (m/s); in [0, 120]")
    ap.add_argument("--EA", type=float, default=ref.EA_DEFAULT, help="axial stiffness (N); default the conductor model's tests' value")
    ap.add_argument("--segments", type=int, default=80, help="segments of the span (Ma 2020: at least 80 elements); >= 2")
    ap.add_argument("--z0", type=float, default=None, help="flow: roughness length (m), default Ma 2020's 0.05 (Section 3); in (0, 1]")
    ap.add_argument("--seed", type=int, default=None, help="flow: conductors.gust_seed, default 1; >= 0")
    ap.add_argument("--dt", type=float, default=0.0, help="time step (s); 0: prescribed 0.05, flow Courant 0.5 on the top wind")
    ap.add_argument("--time", type=float, default=0.0, help="length of the run (s); 0: prescribed 120, flow spin-up + 60 + average")
    ap.add_argument("--spinup", type=float, default=None, help="flow: the flow's spin-up before the statistics (s), > 0; default three flow-throughs at 10 m")
    ap.add_argument("--average", type=float, default=None, help="flow: averaging time of the wind at 10 m and of the statistics (s), default 600; in (0, 600]")
    ap.add_argument("--without-fft", action="store_true", help="erf.use_fft = false, for an ERF built without FFT")
    a = ap.parse_args(argv)
    flow_only = [k for k in ("z0", "seed", "spinup", "average") if getattr(a, k) is not None]
    if a.mode == "prescribed" and flow_only:
        ap.error("--%s: only with --mode flow" % ", --".join(flow_only))
    a.z0 = 0.05 if a.z0 is None else a.z0
    a.seed = 1 if a.seed is None else a.seed
    a.average = 600.0 if a.average is None else a.average
    if a.spinup is not None and not a.spinup > 0.0:
        ap.error("--spinup must be > 0 s")
    if not 50.0 <= a.span <= 1000.0:
        ap.error("--span must be in [50, 1000] m")
    if not 0.0 <= a.wind <= 120.0:
        ap.error("--wind must be in [0, 120] m/s")
    if a.mode == "flow" and not a.wind > 0.0:
        ap.error("--wind must be > 0 with --mode flow: the inflow needs a wind")
    if not (a.EA > 0.0 and a.segments >= 2):
        ap.error("--EA must be > 0 and --segments >= 2")
    if not 0.0 < a.z0 <= 1.0:
        ap.error("--z0 must be in (0, 1] m")
    if a.seed < 0:
        ap.error("--seed must be >= 0")
    if a.dt < 0.0 or a.time < 0.0:
        ap.error("--dt and --time must be >= 0")
    if not 0.0 < a.average <= 600.0:
        ap.error("--average must be in (0, 600] s")
    return a


def log_law(v10, z, z0):
    return v10 * math.log((z + z0) / z0) / math.log((10.0 + z0) / z0)


def main(argv=None):
    a = parse(argv)
    sag = a.span * ref.SAG_RATIO
    support = REFERENCE_HEIGHT + 2.0 * sag / 3.0
    L0 = ref.rest_length(a.span, ref.SAG_RATIO, ref.DRAKE_WEIGHT, a.EA, a.segments)
    mass = ref.DRAKE_WEIGHT / 9.81

    if a.mode == "prescribed":
        cell = 50.0
        lx, ly, lz = 8.0 * cell, a.span + 8.0 * cell, 4.0 * cell
        nx, ny, nz = 8, int(math.ceil(ly / cell)), 8
        ly = ny * cell
        dt = a.dt or 0.05
        time = a.time or 120.0
    else:
        cell, dz = 20.0, 10.0
        lx, lz = 600.0, 240.0
        ly = 40.0 * math.ceil((a.span + 400.0) / 40.0)   # even cell counts (amr.blocking_factor = 2)
        nx, ny, nz = int(lx / cell), int(ly / cell), int(lz / dz)
        top = log_law(a.wind, lz, a.z0)
        dt = a.dt or float("%.3g" % (0.5 * cell / top))
        spin = a.spinup or 3.0 * lx / a.wind
        time = a.time or spin + 60.0 + a.average
        if time < spin + 60.0 + a.average:
            raise SystemExit("make_case: --time %g s ends before the statistics do (spin-up %g s + 60 s + average %g s)"
                             % (time, spin, a.average))
    x0, y0 = 0.5 * lx, 0.5 * (ly - a.span)
    steps = int(math.ceil(time / dt - 1.0e-9))    # the run reaches time

    d = []
    d.append("# Single span after Ma, Bocchini and Christou 2020 (Struct. Saf. 82:101890, Section 3): Drake on %g m," % a.span)
    d.append("# hinged supports %.3f m up (the span's reference height 40 m, two thirds of the %g m sag below them)," % (support, sag))
    if a.mode == "prescribed":
        d.append("# flat ground, still air, a uniform %g m/s wind normal to the span handed to every node." % a.wind)
    else:
        d.append("# flat open terrain (z0 = %g m), neutral log-law inflow with %g m/s at 10 m, k-equation RANS, random gusts." % (a.z0, a.wind))
    d.append("# Written by Structures/Exec/SingleSpanMa2020/make_case.py %s" % " ".join(
        "--%s %s" % (k.replace("_", "-"), v) for k, v in sorted(vars(a).items())
        if k != "out" and v is not None and v is not False and not (a.mode == "prescribed" and k in ("z0", "seed", "spinup", "average"))))
    d.append("max_step  = %d" % steps)
    d.append("stop_time = %g" % (time + 1.0))
    d.append("erf.fixed_dt = %g" % dt)
    d.append("")
    d.append("geometry.prob_lo     = 0. 0. 0.")
    d.append("geometry.prob_hi     = %g %g %g" % (lx, ly, lz))
    d.append("amr.n_cell           = %d %d %d" % (nx, ny, nz))
    d.append("amr.max_level        = 0")
    if a.mode == "prescribed":
        d.append("geometry.is_periodic = 1 1 0")
        d.append("zlo.type = \"SlipWall\"")
        d.append("zhi.type = \"SlipWall\"")
        d.append("")
        d.append("# the air is at rest; the line takes conductors.prescribed_velocity")
        d.append("erf.compressibility_type = Anelastic")
        d.append("erf.timestepping_type    = MidPoint")
        d.append("erf.vert_implicit        = true")
        d.append("erf.use_fft        = %s" % ("false" if a.without_fft else "true"))
        d.append("erf.molec_diff_type = \"None\"")
        d.append("erf.les_type        = \"None\"")
        d.append("erf.init_type = \"uniform\"")
        d.append("erf.prob_name = \"ABL\"")
        d.append("prob.rho_0 = 1.0")
        d.append("prob.T_0   = %g" % THETA)
        d.append("prob.U_0   = 0.0")
    else:
        ustar = a.wind * KAPPA / math.log((10.0 + a.z0) / a.z0)
        d.append("amr.max_grid_size    = 64")
        d.append("amr.blocking_factor  = 2")
        d.append("geometry.is_periodic = 0 0 0")
        d.append("")
        d.append("xlo.type           = \"Inflow\"")
        d.append("xlo.dirichlet_file = \"inflow_profile\"      # z u v w theta")
        d.append("xlo.KE             = %.6g                # %g u*^2, u* = %.4f m/s (%g m/s at 10 m, z0 = %g m)"
                 % (KE_OVER_USTAR2 * ustar * ustar, KE_OVER_USTAR2, ustar, a.wind, a.z0))
        d.append("xhi.type           = \"Outflow\"")
        d.append("ylo.type           = \"SlipWall\"")
        d.append("yhi.type           = \"SlipWall\"")
        d.append("zhi.type           = \"SlipWall\"")
        d.append("zlo.type           = \"surface_layer\"")
        d.append("erf.most.z0        = %g" % a.z0)
        d.append("")
        d.append("erf.compressibility_type = Anelastic")
        d.append("erf.timestepping_type    = MidPoint        # RK2 ignores the implicit vertical solve; MidPoint honours it")
        d.append("erf.vert_implicit        = true")
        d.append("erf.use_fft        = %s" % ("false" if a.without_fft else "true"))
        d.append("")
        d.append("erf.init_type           = \"input_sounding\"")
        d.append("erf.input_sounding_file = \"input_sounding\"    # the inflow profile; surface pressure in hPa")
        d.append("erf.prob_name           = \"ABL\"")
        d.append("erf.use_gravity         = true")
        d.append("")
        d.append("erf.dycore_horiz_adv_type  = \"Upwind_5th\"")
        d.append("erf.dycore_vert_adv_type   = \"Upwind_3rd\"")
        d.append("erf.dryscal_horiz_adv_type = \"Upwind_5th\"")
        d.append("erf.dryscal_vert_adv_type  = \"Upwind_3rd\"")
        d.append("erf.molec_diff_type        = \"None\"")
        d.append("erf.les_type               = \"None\"")
        d.append("erf.rans_type              = \"kEqn\"")
        d.append("erf.dirichlet_k            = true")
        d.append("erf.init_tke_from_ustar    = true")
        d.append("erf.rans_consistent_diffusivities = true")
        d.append("erf.rans_lscale_min        = 1.0")
        d.append("erf.theta_ref              = %g" % THETA)
    d.append("amrex.fpe_trap_invalid = 0      # AMReX's default; the traps are masked while MoorDyn starts and steps the lines")
    d.append("erf.v            = 1")
    d.append("erf.sum_interval = -1")
    d.append("erf.plot_int_1   = -1")
    d.append("erf.check_file   = chk")
    d.append("erf.check_int    = -1")
    d.append("")
    d.append("erf.structures.model = conductors")
    d.append("conductors.lines = S")
    d.append("conductors.S.end_a  = %.3f %.3f %.6f" % (x0, y0, support))
    d.append("conductors.S.end_b  = %.3f %.3f %.6f" % (x0, y0 + a.span, support))
    d.append("conductors.S.length = %.6f      # sag span/50 = %g m under its weight alone (Ma 2020 Section 3)" % (L0, sag))
    d.append("conductors.S.diameter         = %g" % ref.DIAMETER)
    d.append("conductors.S.mass_per_length  = %.6f   # %g N/m, Drake (Ma 2020 Table 3)" % (mass, ref.DRAKE_WEIGHT))
    d.append("conductors.S.axial_stiffness  = %g" % a.EA)
    d.append("conductors.S.drag_coefficient = %g      # Ma 2020 Section 3" % ref.CD)
    d.append("conductors.S.segments         = %d" % a.segments)
    d.append("conductors.S.rated_breaking_strength = %g   # Drake (Ma 2020 Table 3): recorded, nothing lets go" % DRAKE_RATED)
    d.append("conductors.air_density   = %g      # Ma 2020 Section 3" % ref.RHO)
    d.append("conductors.diagnostics_dir = conductors")
    if a.mode == "prescribed":
        d.append("conductors.prescribed_velocity = %g 0. 0." % a.wind)
        d.append("conductors.diagnostics_int     = 20")
    else:
        d.append("conductors.diagnostics_int     = %d" % max(1, int(round(1.0 / dt))))
        d.append("conductors.gust_type = random")
        d.append("conductors.gust_seed = %d" % a.seed)
        d.append("conductors.stats_start          = %g   # spin-up %g s + 60 s" % (spin + 60.0, spin))
        d.append("conductors.capacity_check_start = %g" % (spin + 60.0))
        d.append("conductors.intensity_measure       = wind")
        d.append("conductors.intensity_height        = 10.0")
        d.append("conductors.intensity_averaging_time = %g    # Ma 2020's V_ms is the highest 1-min mean at 10 m" % a.average)
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "inputs"), "w") as f:
        f.write("\n".join(d) + "\n")
    if a.mode == "flow":
        # 1 m rows: ERF interpolates both files linearly (Initialization/ERF_InitBCs.cpp, interpolate_1d), so rows
        # 10 m apart from u = 0 at the ground would hand the 5 m cell 20 m/s for a log-law 35 m/s at 40 m/s
        zs = [float(k) for k in range(int(lz) + 1)]
        with open(os.path.join(a.out, "inflow_profile"), "w") as f:
            for z in zs:
                f.write("%.1f %.6f 0.0 0.0 %.4f\n" % (z, log_law(a.wind, z, a.z0), THETA))
        with open(os.path.join(a.out, "input_sounding"), "w") as f:
            f.write("1000.0 %g 0.0\n" % THETA)
            for z in zs:
                f.write("%.1f %.4f 0.0 %.6f 0.0\n" % (z, THETA, log_law(a.wind, z, a.z0)))
    print("%s: span %g m, L0 %.6f m, supports %.3f m, dt %g s, %d steps (%g s)" % (a.out, a.span, L0, support, dt, steps, time))
    return 0


if __name__ == "__main__":
    sys.exit(main())
