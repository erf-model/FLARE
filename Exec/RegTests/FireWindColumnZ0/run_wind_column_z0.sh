#!/bin/bash
# Run the FireWindColumnZ0 variants and check them.
#
#   [MPIRUN="mpirun -np 2"] ./run_wind_column_z0.sh /path/to/erf_exec
#   SKIP_RUN=1 ./run_wind_column_z0.sh x     # only rerun the checks on existing output
#
# Each pair runs one deck with erf.fire.wind_sample_z0_source = scalar (_sc) and
# surface_layer (_sl):
# uni     inputs_base: z0 = 0.1 m everywhere, the scalar's value
# map     inputs_map: z0 = 0.03 / 0.5 m over the western / eastern half
# map10   inputs_map with erf.fire.wind_ref_ht = 10 (the 20 ft conversion)
# cg      inputs_map with erf.fire.ros_model = cheney_gould (the 10 m wind)
# fmap    map10 on a two-fuel map (the spatial-fuel branch of the WAF)
# bil     inputs_map, bilinear, the 2 m reference wind below the 2.5 m first
#         centre brought down with wind_below_first_cell = log
# Restart parity with bulk_coeff fluxes, whose z0 the surface layer recovers
# every step: rs_a runs 4 steps with a checkpoint at step 2, rs_b restarts
# from it.

set -u
EXE=${1:?usage: run_wind_column_z0.sh /path/to/erf_exec}
PY=${PYTHON:-python3}
SL="erf.fire.wind_sample_z0_source=surface_layer"
BIL="erf.fire.wind_interp=bilinear erf.fire.wind_sample_ht=0 erf.fire.wind_below_first_cell=log erf.fire.wind_ref_ht=2"
FMAP="erf.fire.wind_ref_ht=10 erf.fire.fuel_map.file=fuel_stripes.asc"
BULK="erf.surface_layer.flux_type=bulk_coeff erf.most.Cd=0.004 erf.most.Ch=0.0 erf.most.Cq=0.0 erf.most.surf_moist=0.0"

run () {   # variant deck [extra erf args...]
    v=$1; deck=$2; shift 2
    ${MPIRUN:-} "$EXE" "$deck" erf.fire.plot_file="plt_fire_${v}_" "$@" > "run_${v}.log" 2>&1
}

must_run () {
    if ! run "$@"; then
        echo "run $1 failed; the end of run_$1.log:"
        tail -n 30 "run_$1.log"
        exit 1
    fi
}

if [ "${SKIP_RUN:-0}" != "1" ]; then
    rm -rf plt_fire_* run_*.log Backtrace.* roughness_halves.txt fuel_stripes.asc chk*
    $PY check_firewindcolumnz0.py write_map || exit 1
    must_run uni_sc   inputs_base
    must_run uni_sl   inputs_base $SL
    must_run map_sc   inputs_map
    must_run map_sl   inputs_map  $SL
    must_run map10_sc inputs_map  erf.fire.wind_ref_ht=10
    must_run map10_sl inputs_map  erf.fire.wind_ref_ht=10 $SL
    must_run cg_sc    inputs_map  erf.fire.ros_model=cheney_gould
    must_run cg_sl    inputs_map  erf.fire.ros_model=cheney_gould $SL
    must_run fmap_sc  inputs_map  $FMAP
    must_run fmap_sl  inputs_map  $FMAP $SL
    must_run bil_sc   inputs_map  $BIL
    must_run bil_sl   inputs_map  $BIL $SL
    must_run rs_a     inputs_base $BULK $SL erf.check_int=2
    must_run rs_b     inputs_base $BULK $SL erf.restart=chk00002
fi

status=0
$PY check_firewindcolumnz0.py parity sample ft20 m10 fmap bil restart || status=1
[ $status -eq 0 ] && echo "ALL PASS" || echo "SOME CHECKS FAILED"
exit $status
