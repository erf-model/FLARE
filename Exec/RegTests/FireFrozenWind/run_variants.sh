#!/bin/bash
# FireFrozenWind sampling variants: each deck runs to step 20 and is checked
# against the independent sampling of check_frozen_wind.py with its options.
#
#   [MPIRUN="mpirun -np 2"] ./run_variants.sh /path/to/erf_exec
#
# inputs_nearest_log  nearest column, 2 m (below the first centre), log below it
# inputs_sample_ht    bilinear, sampled at 20 m and brought to 6.1 m (checkpoint)

set -u
EXE=${1:?usage: run_variants.sh /path/to/erf_exec}
status=0

run_case () {
    deck=$1
    shift
    rm -rf plt_fire_* fire_stats*
    ${MPIRUN:-} "$EXE" "$deck" max_step=20 > "run_$deck.log" 2>&1 \
        || { echo "$deck: run failed"; tail -n 30 "run_$deck.log"; status=1; return; }
    echo "== $deck"
    ${PYTHON:-python3} check_frozen_wind.py "$@" || status=1
}

run_case inputs_nearest_log --interp nearest --zref 2.0 --below log
run_case inputs_sample_ht --zref 6.1 --zsample 20.0 --z0 0.1
rm -rf plt_fire_*
exit $status
