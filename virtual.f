      subroutine setvirtual(p,vflav,virtual)
c     Recola one-loop amplitude. The coupling powers selected in
c     recola_generate_process (QCD = res_powst+2, QED = res_powew) pick
c     out the NLO QCD correction to the yt^2 Born.
      use recola_powheg, only: recola_virtual,recola_check_virtual_poles
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      integer, parameter :: nlegs=nlegbornexternal
      real * 8, intent(in)  :: p(0:3,nlegs)
      integer,  intent(in)  :: vflav(nlegs)
      real * 8, intent(out) :: virtual
      real * 8 powheginput
      external powheginput
      logical, save :: ini = .true., inigg = .true.

c     checkvirtpoles 1 : one-off validation of the virtual (IR poles vs
c     the Casimir prediction, UV finiteness, nf dependence). Run it on
c     the first point and on the first gg point.
      if(powheginput("#checkvirtpoles").eq.1) then
         if(ini) then
            ini = .false.
            call recola_check_virtual_poles(p,vflav)
         endif
         if(inigg .and. vflav(1).eq.0 .and. vflav(2).eq.0) then
            inigg = .false.
            call recola_check_virtual_poles(p,vflav)
         endif
      endif

      call recola_virtual(p,vflav,virtual)

      end
