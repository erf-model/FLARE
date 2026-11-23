#!/bin/bash
# FireFrozenWind ember spotting: each deck runs to step 20 and check_spot3d.py
# replays the flight of every brand launched at step 20.
#
#   [MPIRUN="mpirun -np 2"] ./run_spot3d.sh /path/to/erf_exec
#
# inputs_spot3d        brands drift in the frozen 3-D GUST wind (drift_wind = frozen)
# inputs_spot_profile  the same run on the log-scaled reference wind (the default)
# inputs_spot_nearest_log  frozen, in the nearest GUST column, log profile below
#                      its lowest cell centre (wind_interp, wind_below_first_cell)

set -u
EXE=${1:?usage: run_spot3d.sh /path/to/erf_exec}
status=0

run_case () {
    deck=$1
    drift=$2
    shift 2
    rm -rf plt_fire_* fire_stats*
    ${MPIRUN:-} "$EXE" "$deck" max_step=20 > "run_$deck.log" 2>&1 \
        || { echo "$deck: run failed"; tail -n 30 "run_$deck.log"; status=1; return; }
    echo "== $deck"
    ${PYTHON:-python3} check_spot3d.py --drift "$drift" "$@" || status=1
}

run_case inputs_spot3d frozen
run_case inputs_spot_profile profile
run_case inputs_spot_nearest_log frozen --interp nearest --below log
rm -rf plt_fire_*
exit $status
