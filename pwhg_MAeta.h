c     MAapprox: massification scale factor eta, Q = eta*sqrt(s_bbH).
c     Set from etascfact in init_processes, changed per weight by the
c     reweighting (rwl_setup_param_weights_user.f, key etascfact).
c     ma_lam: lambdascvar, share of m_b^2 moved into the b transverse
c     momentum by the massless mapping (map_lambda_bbH), MAapprox_mapping
c     2. Must be > 0: lambda = 0 lets pT' -> 0, collinear (IR) unsafe
c     for the massless amplitude.
      real * 8 ma_eta,ma_lam
      common/cma_eta/ma_eta,ma_lam
