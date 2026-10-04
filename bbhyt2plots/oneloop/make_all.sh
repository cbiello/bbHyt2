#!/bin/bash
# M1M0 comparison, exact vs massified (MA), LHE level.
#   ./make_all.sh <TAG> <exact run dir> <massified run dir>
#   e.g. ./make_all.sh HT4 run-M1M0-exact run-M1M0-massified-flow
#        ./make_all.sh mH2 run-M1M0-exact-mH2 run-M1M0-massified-mH2
# (run dirs relative to RUNS, the bbHyt2 POWHEG directory)
#   0. lhef_analysis on every pwgevents-*.lhe whose -W1.top is missing
#   1. merge the per-seed tops (W1 central, W3 eta=2, W4 eta=1/2,
#      and, if present, W5 lambda=2, W6 lambda=1/2)
#   2. rebin, 3. convert to datfiles, 4. stat error bars, 5. plot
set -e
TAG=${1:?usage: make_all.sh <TAG> <exact dir> <mass dir>}
EXA=${2:?}; MAS=${3:?}
RUNS=${RUNS:-/scratch/cbiello/POWHEG-BOX-RES/bbHyt2}
MERGE=${MERGE:-/scratch/cbiello/bin/mergedata}
REBIN="pt_H=4 ptzoom_H=5 y_H=2 eta_H=2"
cd "$(dirname "$0")"; HERE=$PWD

for d in $EXA $MAS; do
   ( cd $RUNS/$d
     for ev in pwgevents-*.lhe; do
        s=${ev#pwgevents-}; s=${s%.lhe}
        [ -f pwgLHEF_analysis-$s-W1.top ] || \
           (echo $ev | ../lhef_analysis > lhefana-$s.log 2>&1 &)
     done; wait
     while pgrep -f "^../lhef_analysis" > /dev/null; do sleep 2; done )
done

T=topfiles/$TAG-lhe; mkdir -p $T
merge () { rm -f fort.12; $MERGE 1 "${@:2}" > /dev/null; mv fort.12 $T/$1; }
merge exact-central.top $RUNS/$EXA/pwgLHEF_analysis-[0-9]*-W1.top
merge mass-central.top  $RUNS/$MAS/pwgLHEF_analysis-[0-9]*-W1.top
merge mass-eta2.top     $RUNS/$MAS/pwgLHEF_analysis-[0-9]*-W3.top
merge mass-etaH.top     $RUNS/$MAS/pwgLHEF_analysis-[0-9]*-W4.top
VARS="mass-eta2 mass-etaH"
if ls $RUNS/$MAS/pwgLHEF_analysis-[0-9]*-W6.top > /dev/null 2>&1; then
   merge mass-lam2.top $RUNS/$MAS/pwgLHEF_analysis-[0-9]*-W5.top
   merge mass-lamH.top $RUNS/$MAS/pwgLHEF_analysis-[0-9]*-W6.top
   VARS="$VARS mass-lam2 mass-lamH"
fi
echo "$TAG: merged $(ls $RUNS/$EXA/pwgLHEF_analysis-[0-9]*-W1.top | wc -l) exact," \
     "$(ls $RUNS/$MAS/pwgLHEF_analysis-[0-9]*-W1.top | wc -l) massified seeds"

R=topfiles/$TAG-lhe-rebin; mkdir -p $R
for f in exact-central mass-central $VARS; do
   python3 topfiles/rebin_top.py $T/$f.top $R/$f.top $REBIN
done
# MA band: bin-by-bin envelope of all the variations
python3 topfiles/envelope_top.py $R/mass-min.top $R/mass-max.top \
   $(for v in $VARS; do echo $R/$v.top; done)

mkdir -p datfiles; cd datfiles
python3 ../topfiles/convert_top_files_to_plots.py ../$R/exact-central.top \
   ../$R/exact-central.top ../$R/exact-central.top $TAG-exact-scales > /dev/null
python3 ../topfiles/convert_top_files_to_plots.py ../$R/mass-central.top \
   ../$R/mass-min.top ../$R/mass-max.top $TAG-mass-scales > /dev/null
cd ..

python3 make_ebars.py $R $TAG
python3 myscprit-oneloop.py $TAG 2>/dev/null | grep -i "plot" || true

# one file with all plots of this TAG
G=gnuplot_minnlo
pdfunite $G/ptzoom_H__$TAG-.pdf $G/pt_H__$TAG-.pdf $G/y_H__$TAG-.pdf \
         $G/eta_H__$TAG-.pdf $G/allplots-$TAG-lhe.pdf
