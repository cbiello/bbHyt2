c     Process-specific reweighting keys (see rwl_setup_param_weights.f
c     of POWHEG-BOX-RES), as in MiNNLOPS_res/ttH_NLO:
c        etascfact   massification scale factor eta (MAapprox 1)
c        lambdascvar share of m_b^2 moved into pT by the massless
c                    mapping (MAapprox_mapping 2), must be > 0
c     count = 0: store, count = -1: restore, count > 0: set weight count
      subroutine rwl_setup_params_weights_user(count)
      implicit none
      integer count
      include 'pwhg_rwl.h'
      include 'pwhg_MAeta.h'
      real * 8, save :: old_ma_eta,old_ma_lam
      logical rwl_keypresent
      real * 8 val
      if(count==0) then
         old_ma_eta = ma_eta
         old_ma_lam = ma_lam
      elseif(count == -1) then
         ma_eta = old_ma_eta
         ma_lam = old_ma_lam
      else
         if(rwl_keypresent(count,'etascfact',val)) then
            ma_eta = val
         endif
         if(rwl_keypresent(count,'lambdascvar',val)) then
            if(val.le.0d0) then
               write(*,*) ' lambdascvar must be > 0 (IR safety)'
               call exit(-1)
            endif
            ma_lam = val
         endif
      endif
      end
