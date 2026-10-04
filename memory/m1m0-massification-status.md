---
name: m1m0-massification-status
description: "Current state (2026-10-04) of the bbHyt2 M1M0 exact-vs-massification study and the full NLO+PS run - code, runs, results, open problem"
metadata:
  node_type: memory
  type: project
  originSessionId: 660413b7-bdbc-4081-adcf-1fd98bdd9ba7
  modified: 2026-10-04T08:32:36.060Z
---

Study: M1M0 (one-loop virtual in place of the Born: `virtinborn 1 bornonly 1 LOevents 1`), exact (`MAapprox 0`) vs massified (`MAapprox 1`), on Higgs observables at LHE level. All code committed in git@github.com:cbiello/bbHyt2.git (87fface code, f4f0e90 plots + submitjobs.md + NLO+PS card).

**Massification code** (virtual.f / recola.f, follows MiNNLOPS_res/ttH_NLO):
- one Recola session: massive processes + massless images g g/q q~ -> H s s~ (d d~ for s s~ initial state), MS=0, massive b in loops; F1 = Mitov-Moch without loop-content terms.
- at Q = eta*sqrt(s_bbH) of the massive point: VQ = (2F1 + V0/B0)*B_massive, then exact RGE flow to muR with massive poles (MA_runQtomuR; c2 analytic). See [[massification-at-Q-then-rge]].
- massless map: `MAapprox_mapping 2` = lambda family, `lambdascvar` (see [[lambdascvar-mapping]]); reweighting keys `etascfact`, `lambdascvar` (rwl_setup_param_weights_user.f, pwhg_MAeta.h).
- `runningscales 5` = fixed muR=muF=mH/2 (Born_phsp.f).
- checks: MAflowcheck 1 (exact flow, <1e-3; residual ~m_b^2, unexplained), MAdebug 2/3 (fort.78/79 pointwise vs exact, lambda scan).
- build `make ANALYSIS=MA pwhg_main lhef_analysis` (default MiNNLO analysis is wrong here). analysis-MA had two bugs, fixed (multi_plot_setup args; NLO mode must read HEPEVT).
- user wants only central + eta/lambda weights in these runs, no muR/muF variation weights.

**Runs (yueh, /scratch/cbiello/POWHEG-BOX-RES/bbHyt2)**, sigma in pb:
- HT4: run-M1M0-exact (5x20k) 0.4582; run-M1M0-massified-flow (eta weights) 0.4467; run-M1M0-massified-lam (eta+lambda weights) 0.4467, band [0.4364, 0.4600].
- mH2: run-M1M0-{exact,massified}-mH2 (10x15k): 0.4184 vs 0.4089, eta band +-0.16%.
- MA/exact: ~0.975 inclusive, 0.92-0.95 at pT_H < 20 GeV; lambda band reaches exact there within ~1-1.5 sigma stat.
- run-M1M0-massified (first attempt) is obsolete and was not committed.
- Plots: see [[bbhytplots-repo]].

**Open problem:** full NLO+PS exact HT4 (run-NLOPS-exact-HT4, 28 seeds): stages 1-3 fine (sigma_NLO 0.7115+-0.0024 pb; benchmark 0.681, 2307.09992 0.689), but stage 4 wrote no event in 8.5 h and a 3-event test also stalled. Cause not found (no ptrace/perf on yueh). User stopped it and will run on their cluster with submitjobs.md (SLURM chain + stage-4 smoke test).

**How to apply:** start from here; for NLO+PS the stage-4 stall must be understood before producing events.
