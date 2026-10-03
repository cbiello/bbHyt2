#!/bin/bash
# Stage 1-4 run of the M1M0 comparison, 5 seeds in parallel per stage.
#   ./run_M1M0.sh <dir>
# virtinborn puts the virtual in place of the Born, so the events are
# M1M0 events; bornonly + LOevents, no radiation.
set -u
DIR=${1:?usage: run_M1M0.sh <rundir>}
EXE=$(cd "$(dirname "$0")" && pwd)/pwhg_main
NSEED=$(wc -l < "$DIR/pwgseeds.dat")
cd "$DIR" || exit 1

setkey () {   # setkey <key> <value>: drop every occurrence, append one
  sed -i -E "/^[[:space:]]*#?$1[[:space:]]/d" powheg.input
  echo "$1 $2" >> powheg.input
}

runstage () {  # runstage <label>
  echo "=== $DIR : $1  ($(date +%H:%M:%S)) ==="
  for i in $(seq 1 "$NSEED"); do
     echo "$i" | $EXE > run-$1-$i.log 2>&1 &
  done
  wait
  grep -l -iE "exiting|forrtl|Error in" run-$1-*.log 2>/dev/null && \
     echo "  !! check the logs above"
}

# iteration 2 builds on iteration 1, so use-old-grid must be on;
# stages 2-4 look for the grid of the LAST iteration, hence xg 2.
setkey use-old-grid 1
setkey parallelstage 1; setkey xgriditeration 1; runstage st1-xg1
setkey parallelstage 1; setkey xgriditeration 2; runstage st1-xg2
setkey parallelstage 2; setkey xgriditeration 2; runstage st2
setkey parallelstage 3; setkey xgriditeration 2; runstage st3
setkey parallelstage 4; setkey xgriditeration 2; runstage st4
echo "=== $DIR : done ($(date +%H:%M:%S)) ==="
ls -1 pwgPOWHEG*.top pwgNLO*.top 2>/dev/null | head
