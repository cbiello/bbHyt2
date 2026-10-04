# Running bbHyt2 on a cluster (SLURM)

How to run the bbH (y_t^2) generator in POWHEG-BOX-RES parallel mode on a
batch system: full NLO+PS at LHE level, or the M1M0 (virtual-in-Born) runs.
The local driver `run_M1M0.sh` runs the same stages one after the other on a
single machine; this file shows the same stages as SLURM array jobs.

## 1. Build (once, on the cluster)

```bash
cd POWHEG-BOX-RES/bbHyt2
make ANALYSIS=MA pwhg_main       # MA = Higgs/b-jet analysis used for the plots
make ANALYSIS=MA lhef_analysis   # LHE-level analysis of the event files
```
Needs Recola2 (HEFT model), COLLIER, LHAPDF with `NNPDF31_nlo_as_0118_nf_4`
(LHAID 320500) and FastJet; paths are set at the top of the `Makefile`.
Check that `pwhg_main` links `pwhg_analysis-MA.o` (the default `ANALYSIS`
is MiNNLO, which is the wrong analysis for this study).

## 2. Run directory

```bash
mkdir -p run-NLOPS-exact-HT4/output_cll        # Recola writes its log here
cp <card> run-NLOPS-exact-HT4/powheg.input
seq 8001 8200 > run-NLOPS-exact-HT4/pwgseeds.dat   # one line per job
```
The number of lines of `pwgseeds.dat` must be at least the array size.

Cards in this repository:

| run | card | what |
|---|---|---|
| full NLO+PS, exact virtual, H_T/4 | `run-NLOPS-exact-HT4/powheg.input` | Born+virtual+real+POWHEG radiation, LHE level |
| M1M0 exact, H_T/4 | `run-M1M0-exact/powheg.input` | `virtinborn 1 bornonly 1 LOevents 1` |
| M1M0 massified, H_T/4, eta+lambda | `run-M1M0-massified-lam/powheg.input` | `MAapprox 1`, `MAapprox_mapping 2`, weights eta, lambda = 1/2, 2 |
| M1M0 exact / massified, mH/2 | `run-M1M0-*-mH2/powheg.input` | `runningscales 5` |

The full NLO+PS card differs from the benchmark card
(`powhegrun-benchmark2307.09992/powheg.input`, validated against
2307.09992) only in statistics. Essential physics keys:
`runningscales 4` (H_T/4), `btlscalereal 1`, `wilsonlog 1`, and for NLO+PS
a real `nubound` (the M1M0 cards have `nubound 100`, fine without
radiation, far too small for NLO+PS).

Statistics knobs (per job): `ncall1`, `ncall2`, `nubound`, `numevts`.
The grids of all jobs are combined, so the total grid statistics is
(number of jobs) x ncall.

## 3. Stages

POWHEG parallel mode needs, for every job, `manyseeds 1` and the keys
below; each job reads its seed index from stdin. All jobs of a stage must
finish before the next stage starts (they read each other's grid files),
hence one array job per stage, chained with `--dependency=afterok`.

| step | keys set in powheg.input | output |
|---|---|---|
| 1a | `parallelstage 1  xgriditeration 1  use-old-grid 1` | `pwg-xg1-xgrid-*.dat` |
| 1b | `parallelstage 1  xgriditeration 2  use-old-grid 1` | `pwg-xg2-xgrid-*.dat` |
| 2  | `parallelstage 2  xgriditeration 2` | `pwggrid-*.dat`, `pwg-*-st2-stat.dat` |
| 3  | `parallelstage 3  xgriditeration 2` | upper bounds `pwgubound-*`, `pwg*upb-*` |
| 4  | `parallelstage 4  xgriditeration 2` | `pwgevents-NNNN.lhe` |

`use-old-grid 1` is needed so that xgrid iteration 2 builds on iteration 1.
Stages 2-4 look for the grids of the LAST xgrid iteration (hence
`xgriditeration 2`). Stage 2 must not find old `*fullgrid*.dat` files.

Since the keys change between stages, give each stage its own copy of the
card and point POWHEG at it through the run directory, or (simplest)
edit the card in a short serial step between the array jobs, as below.

### Job script `pwg_stage.sh`

```bash
#!/bin/bash
#SBATCH --job-name=bbHyt2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=slurm-%x-%A_%a.out
# usage: sbatch --array=1-N pwg_stage.sh <rundir> <label>
cd "$1" || exit 1
echo "$SLURM_ARRAY_TASK_ID" | ../pwhg_main > "run-$2-$SLURM_ARRAY_TASK_ID.log" 2>&1
```

### Card editing `setstage.sh`

```bash
#!/bin/bash
# usage: setstage.sh <rundir> <parallelstage> <xgriditeration>
cd "$1" || exit 1
setkey () { sed -i -E "/^[[:space:]]*#?$1[[:space:]]/d" powheg.input; echo "$1 $2" >> powheg.input; }
setkey use-old-grid 1
setkey parallelstage "$2"
setkey xgriditeration "$3"
```

### Submission `submit.sh`

```bash
#!/bin/bash
# usage: ./submit.sh <rundir> <njobs>
D=$1; N=$2
./setstage.sh $D 1 1
j1=$(sbatch --parsable --array=1-$N --time=01:00:00 pwg_stage.sh $D st1-xg1)
e1=$(sbatch --parsable --dependency=afterok:$j1 --wrap "./setstage.sh $D 1 2")
j2=$(sbatch --parsable --dependency=afterok:$e1 --array=1-$N --time=01:00:00 pwg_stage.sh $D st1-xg2)
e2=$(sbatch --parsable --dependency=afterok:$j2 --wrap "./setstage.sh $D 2 2")
j3=$(sbatch --parsable --dependency=afterok:$e2 --array=1-$N --time=02:00:00 pwg_stage.sh $D st2)
e3=$(sbatch --parsable --dependency=afterok:$j3 --wrap "./setstage.sh $D 3 2")
j4=$(sbatch --parsable --dependency=afterok:$e3 --array=1-$N --time=01:00:00 pwg_stage.sh $D st3)
e4=$(sbatch --parsable --dependency=afterok:$j4 --wrap "./setstage.sh $D 4 2")
j5=$(sbatch --parsable --dependency=afterok:$e4 --array=1-$N --time=24:00:00 pwg_stage.sh $D st4)
echo "stage-4 array job: $j5"
```

Stages 1-3 only need to be run once; for more events, add seeds to
`pwgseeds.dat` and submit further stage-4 arrays (`--array=N+1-M`): new
event files only, same grids. Keep `numevts` per job moderate (a few
thousand) so that jobs fit in the queue limits; statistics come from the
number of jobs.

## 4. Timings measured on yueh (full NLO+PS, exact virtual)

Per job, 28 jobs in parallel, card `run-NLOPS-exact-HT4`:

| step | settings | wall time per job |
|---|---|---|
| 1a / 1b | `ncall1 20000` | ~10 min each |
| 2 | `ncall2 50000` | ~25 min |
| 3 | `nubound 20000` | ~8 min |
| 4 | `numevts 30000` | **see below** |

Results of stages 1-3 with 28 jobs: sigma_NLO = 0.7115 +- 0.0024 pb
(chi2/ndf over jobs 1.5), stage-3 total 0.720 pb, |negative| fraction 7%;
the benchmark gives 0.681 +- 0.010 pb (2307.09992: 0.689 pb).

The M1M0 runs (no radiation) generate ~4000 events/hour/job (exact
virtual) and ~1000-1300 events/hour/job (massified, 4 Recola one-loop calls
per point plus reweighting).

## 5. WARNING: stage 4 of the full NLO+PS run

On yueh, stage 4 of `run-NLOPS-exact-HT4` ran 8.5 h on 28 jobs at full CPU
without writing a single event (event file = header only; events are
flushed in ~1 MB blocks, i.e. fewer than ~800 events per job in 8.5 h), and
a 3-event test did not finish in minutes either. Stages 1-3 are fine. This
is not a question of statistics: before submitting stage 4 at scale, run a
smoke test

```bash
sed -i -E 's/^numevts .*/numevts 10/' powheg.input
echo 1 | ../pwhg_main > st4-test.log 2>&1      # with the stage-3 output present
```
and check that it ends with `pwgevents-0001.lhe` holding 10 events in a
reasonable time (minutes). If it does not, the cause has to be found first
(candidates: radiation generation / upper-bound vetoes with `bornktmin 0`
and massive b emitters, or the btilde unweighting with the mint upper
bounds); `pwgcounters-st4-*.dat` at the end of a finished test lists the
upper-bound violations and vetoed calls.

## 6. After stage 4

```bash
cd <rundir>
for ev in pwgevents-*.lhe; do echo $ev | ../lhef_analysis > lhefana-${ev#pwgevents-}.log 2>&1; done
```
(one job each, or an array job) gives `pwgLHEF_analysis-NNNN-W*.top`
(`W1` = central, then one file per weight of the card). Plots:
`git@github.com:cbiello/bbhytplots.git`, `oneloop/make_all.sh` (merging,
rebinning, error bars; see its README).
