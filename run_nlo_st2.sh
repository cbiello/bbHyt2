#!/bin/bash
# Stage-2 NLO plots (testplots 1) on existing stage-1 grids, then merge.
#   ./run_nlo_st2.sh <dir>     (dir holds powheg.input, pwgseeds.dat, xg grids)
set -u
DIR=${1:?usage: run_nlo_st2.sh <dir>}
EXE=$(cd "$(dirname "$0")" && pwd)/pwhg_main
cd "$DIR" || exit 1
NSEED=$(wc -l < pwgseeds.dat)
rm -f pwgfullgrid-* pwggrid-* pwgbtlupb-* pwg-*-NLO.top pwg-*-st2-stat.dat \
      pwg-st2-xgrid-* sigborn_equiv-* pwgcounters-st2-* run-st2-*.log
for i in $(seq 1 "$NSEED"); do echo "$i" | $EXE > run-st2-$i.log 2>&1 & done
wait
grep -l -iE "exiting|forrtl|Error in" run-st2-*.log && echo "  !! check the logs above"
rm -f fort.12
ls pwg-*-NLO.top | /scratch/cbiello/bin/mergedata 1 $(ls pwg-*-NLO.top) > merge.log 2>&1
mv fort.12 pwg-NLO-merged.top
echo "$DIR: merged $(ls pwg-*-NLO.top | wc -l) files -> pwg-NLO-merged.top"
