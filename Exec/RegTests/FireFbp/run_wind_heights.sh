#!/bin/bash
# The empirical spread models take the wind at a fixed height: FBP the 10 m open
# wind, the Rothermel wind adjustment factor the 20 ft (6.1 m) wind. The fire
# layer brings the wind sampled at erf.fire.wind_ref_ht to that height with the
# neutral log profile of erf.fire.wind_sample_z0 (0.1 m). Until November 2026 it
# did not: FBP read the 6.1 m wind as its 10 m wind (11 % low), and the WAF took
# a 10 m wind as the 20 ft wind (12 % high).
#
#   [MPIRUN="mpirun -np 1"] ./run_wind_heights.sh /path/to/erf_exec
#
# The wind is a uniform 8 m/s at every height, so only the conversion separates
# the runs below; each is one step and the head rate at step 1 is read.
#   1. FBP C-2 sampled at 6.1 m: the head rate equals the independent
#      fbp_reference.py at the 10 m wind 8 ln(10/0.1) / ln(6.1/0.1) = 8.962 m/s
#      to 1e-8 (without the conversion it is the 8 m/s rate, 15 % lower).
#   2. Rothermel sampled at 10 m spreads slower than at 6.1 m (the 10 m wind is
#      brought down to 20 ft before the WAF), by more than 5 %; without the
#      conversion the two are equal.
set -u
EXE=${1:?usage: run_wind_heights.sh /path/to/erf_exec}
shift || true
run() { ${MPIRUN:-} "$EXE" "$2" max_step=1 erf.fire.plot_int=-1 erf.plot_int_1=-1 erf.check_int=-1 $3 > "wh_$1.log" 2>&1 \
          || { echo "run $1 failed (see wh_$1.log)"; tail -20 "wh_$1.log"; exit 1; }; }
run fbp6   inputs_fbp_c2    "erf.fire.wind_ref_ht=6.1"
run roth10 inputs_rothermel "erf.fire.wind_ref_ht=10.0"
run roth6  inputs_rothermel "erf.fire.wind_ref_ht=6.1"

ros1()  { grep 'max_ROS=' "wh_$1.log" | head -1 | sed 's/.*max_ROS=\([^ ]*\) .*/\1/'; }
wind1() { grep 'Max reference wind' "wh_$1.log" | head -1 | sed 's/.*Max reference wind: \([^ ]*\) m.*/\1/'; }

# every value read must be a positive number, or the check fails (an empty value
# would otherwise compare as a string and pass)
num() { awk -v x="$1" 'BEGIN { print (x ~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/ && x + 0 > 0) ? 1 : 0 }'; }
status=0
check() { if [ "$2" = "1" ]; then echo "  $1: PASS"; else echo "  $1: FAIL"; status=1; fi; }
for v in "$(wind1 fbp6)" "$(ros1 fbp6)" "$(ros1 roth10)" "$(ros1 roth6)"; do
    [ "$(num "$v")" = "1" ] || { echo "  a value read from the logs is missing or not positive: '$v'"; echo "SOME CHECKS FAILED"; exit 1; }
done
w10=$(awk -v w="$(wind1 fbp6)" 'BEGIN { printf "%.12g", w * log(10.0/0.1) / log(6.1/0.1) }')
ref=$(${PYTHON:-python3} fbp_reference.py C2 90 60 "$w10")
[ "$(num "$ref")" = "1" ] || { echo "  fbp_reference.py gave no rate: '$ref'"; echo "SOME CHECKS FAILED"; exit 1; }
d=$(awk -v a="$(ros1 fbp6)" -v b="$ref" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; printf "%.2e", d }')
check "FBP C-2 at a 6.1 m wind of $(wind1 fbp6) m/s reads the 10 m wind $w10 m/s ($(ros1 fbp6) vs $ref m/s, rel. diff $d)" \
      "$(awk -v d="$d" 'BEGIN { print (d + 0 <= 1e-8) ? 1 : 0 }')"
r=$(awk -v a="$(ros1 roth10)" -v b="$(ros1 roth6)" 'BEGIN { printf "%.4f", a / b }')
check "Rothermel with the wind sampled at 10 m spreads slower than at 6.1 m (ratio $r, want < 0.95)" \
      "$(awk -v r="$r" 'BEGIN { print (r + 0 < 0.95) ? 1 : 0 }')"
[ $status -eq 0 ] && echo "ALL PASS" || echo "SOME CHECKS FAILED"
exit $status
