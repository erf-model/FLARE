#!/bin/bash
# FireFrozenWind restart: a restart re-reads the GUST file at start-up, so the
# fire after a restart equals the fire of a straight run.
#
#   [MPIRUN="mpirun -np 2"] ./run_restart.sh /path/to/erf_exec
#
# 1. inputs_plt to step 20 with a checkpoint at step 10; plt_fire_00020 is
#    kept as straight_00020.
# 2. The restart from chk00010 to step 20.
# 3. fire_wind_ref_u/v, fire_extract_z and fire_phi of the two step-20 fire
#    plotfiles must be equal bit for bit (compare_restart.py).

set -u
EXE=${1:?usage: run_restart.sh /path/to/erf_exec}
rm -rf chk000* plt_fire_* straight_00020

${MPIRUN:-} "$EXE" inputs_plt max_step=20 erf.check_int=10 > run_straight.log 2>&1 \
    || { echo "straight run failed"; tail -n 30 run_straight.log; exit 1; }
[ -d plt_fire_00020 ] || { echo "no plt_fire_00020"; exit 1; }
[ -d chk00010 ] || { echo "no checkpoint chk00010"; exit 1; }
mv plt_fire_00020 straight_00020

${MPIRUN:-} "$EXE" inputs_plt max_step=20 erf.restart=chk00010 > run_restart.log 2>&1 \
    || { echo "restart run failed"; tail -n 30 run_restart.log; exit 1; }
if [ "$(grep -c 'Frozen wind: read the GUST plotfile' run_restart.log)" = "1" ]; then
    echo "reread   PASS  the restart read the GUST plotfile again"
    status=0
else
    echo "reread   FAIL  the restart did not read the GUST plotfile"
    status=1
fi
${PYTHON:-python3} compare_restart.py straight_00020 plt_fire_00020 || status=1
rm -rf chk000*
exit $status
