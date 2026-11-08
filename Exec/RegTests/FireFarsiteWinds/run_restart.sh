#!/bin/bash
# FireFarsiteWinds restart: the .atm file and the grids after a restart from
# an older checkpoint equal those of a straight run.
#
#   [MPIRUN="mpirun -np 2"] ./run_restart.sh /path/to/erf_exec
#
# 1. inputs_prescribed to step 125, with a checkpoint at step 61 (t = 61 s).
#    The straight .atm is kept as straight.atm (lines for steps 1, 61, 121).
# 2. The restart from chk00061 to step 125 in the same directory. Its first
#    write (step 121, t = 120 s) must drop the straight run's step-121 line
#    before appending its own, so the file ends as the straight one did.
# The prescribed wind is uniform, so the step-121 grids must also be equal.

set -u
EXE=${1:?usage: run_restart.sh /path/to/erf_exec}
D=winds_prescribed
rm -rf chk000* "$D" straight.atm straight_121_vel.asc

${MPIRUN:-} "$EXE" inputs_prescribed max_step=125 erf.check_int=61 > run_straight.log 2>&1 \
    || { echo "straight run failed"; tail -n 30 run_straight.log; exit 1; }
cp "$D/prescribed.atm" straight.atm
cp "$D/prescribed_000121_vel.asc" straight_121_vel.asc
[ -d chk00061 ] || { echo "no checkpoint chk00061"; exit 1; }

${MPIRUN:-} "$EXE" inputs_prescribed max_step=125 erf.restart=chk00061 > run_restart.log 2>&1 \
    || { echo "restart run failed"; tail -n 30 run_restart.log; exit 1; }

status=0
echo "straight .atm:"; cat straight.atm
echo "after restart:"; cat "$D/prescribed.atm"
if cmp -s straight.atm "$D/prescribed.atm"; then
    echo "atm      PASS  the restart's .atm equals the straight run's"
else
    echo "atm      FAIL  the restart's .atm differs from the straight run's"; status=1
fi
if [ "$(grep -c 'FARSITE winds t=' run_restart.log)" = "1" ]; then
    echo "writes   PASS  the restart wrote one grid pair (step 121)"
else
    echo "writes   FAIL  the restart wrote $(grep -c 'FARSITE winds t=' run_restart.log) grid pairs, expected 1"; status=1
fi
if cmp -s straight_121_vel.asc "$D/prescribed_000121_vel.asc"; then
    echo "grid     PASS  step-121 speed grid equal"
else
    echo "grid     FAIL  step-121 speed grid differs"; status=1
fi
rm -rf chk000*
exit $status
