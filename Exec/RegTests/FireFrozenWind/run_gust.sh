#!/bin/bash
# Remake gust_hill_plt and gust_hill_chk, the frozen wind of this case, with
# GUST-CFD (erf-model/gust-abl). Not part of CTest: the committed directories
# are what the tests read.
#
#   [MPIRUN="mpirun -np 4"] ./run_gust.sh /path/to/gust
#
# 1. make_hill_terrain.py writes terrain_gaussian_hill.txt.
# 2. GUST runs inputs_gust_hill to steady (about 6500 iterations, 1.5 min on
#    4 ranks) and writes gust_hill_plt<step> and gust_hill_chk<step>.
# 3. Those are moved to gust_hill_plt and gust_hill_chk (the names the decks
#    read; GUST's files do not record their directory name).

set -u
GUST=${1:?usage: run_gust.sh /path/to/gust}
python3 make_hill_terrain.py terrain_gaussian_hill.txt || exit 1
${MPIRUN:-} "$GUST" inputs_gust_hill > gust.log 2>&1 \
    || { echo "GUST failed"; tail -n 30 gust.log; exit 1; }
tail -n 6 gust.log
plt=$(ls -d gust_hill_plt[0-9]* | tail -n 1)
chk=$(ls -d gust_hill_chk[0-9]* | tail -n 1)
[ -n "$plt" ] && [ -n "$chk" ] || { echo "no GUST plotfile or checkpoint"; exit 1; }
rm -rf gust_hill_plt gust_hill_chk
mv "$plt" gust_hill_plt
mv "$chk" gust_hill_chk
echo "wrote gust_hill_plt (from $plt) and gust_hill_chk (from $chk)"
