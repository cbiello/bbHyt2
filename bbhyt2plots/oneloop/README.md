# oneloop: M1M0 for bbH (y_t^2), exact vs massification

Same organisation as `tthplots/oneloop`.

M1M0 = one-loop virtual put in place of the Born (`virtinborn 1`, POWHEG-BOX-RES/bbHyt2),
LHE level, 13 TeV. Two central scales (TAG):
* `HT4`: mu_R = mu_F = H_T/4, 5 x 20k events per run;
* `mH2`: mu_R = mu_F = m_H/2 (fixed), 10 x 15k events per run;
* `HT4lam`: as HT4, with the lambda family of massless maps (`MAapprox_mapping 2`,
  `lambdascvar`): pT'^2 = pT^2 + lambda m_b^2, p_z kept, E'^2 = E^2 + (lambda-1) m_b^2
  (gg in the lab frame, q qbar in the partonic frame), central lambda = 1, lambda > 0 for IR
  safety. Band = bin-by-bin envelope of eta = 1/2, 2 and lambda = 1/2, 2.

* **exact**: massive one-loop amplitude (Recola, HEFT).
* **MA**: massification evaluated at Q = eta * sqrt(s_bbH) on the massless (b -> s) image,
  multiplied by exact Born / massless Born, then moved to mu_R with the exact one-loop RGE
  (poles of the massive amplitude), as in MiNNLOPS_res/ttH_NLO.
  Band: bin-by-bin envelope of the stage-4 weights (eta = 1/2, 2; for HT4lam also
  lambda = 1/2, 2, `topfiles/envelope_top.py`); the exact result is independent of both.
* Statistical errors: vertical bars at the bin centres (`ebars/`, `make_ebars.py`),
  as in the bbH GMVFNS plots. Ratio panel normalised to exact.

Chain (`./make_all.sh <TAG> <exact run dir> <massified run dir>`, dirs relative to bbHyt2):
```
./make_all.sh HT4 run-M1M0-exact     run-M1M0-massified-flow
./make_all.sh mH2 run-M1M0-exact-mH2 run-M1M0-massified-mH2
```
0. `lhef_analysis` (ANALYSIS=MA) on every event file not yet analysed
1. merge per-seed `pwgLHEF_analysis-*-W*.top` (`W1` central, `W3` eta=2, `W4` eta=1/2) -> `topfiles/<TAG>-lhe/`
2. rebin (`topfiles/rebin_top.py`) -> `topfiles/<TAG>-lhe-rebin/`
3. `topfiles/convert_top_files_to_plots.py` -> `datfiles/<TAG>-{exact,mass}-scales-run/`
4. `make_ebars.py` -> `ebars/<obs>__<TAG>-{exact,mass}.{ebar,rebar}`
5. `myscprit-oneloop.py <TAG>` -> `gnuplot_minnlo/*__<TAG>-.pdf`, all in `gnuplot_minnlo/allplots-<TAG>-lhe.pdf`
