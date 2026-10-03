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

c     checkvirtpoles 1 : validate the virtual (IR poles, UV finiteness,
c     nf dependence) on the first point and on the first gg point.
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

c     MAapprox 1 : massification at Q = sqrt(s_bbH) on the massless
c     (b -> s) image, normalised to the massive Born and run to muR
c     with the poles of the massive amplitude (as in ttH_NLO).
      if(powheginput("#MAapprox").eq.1) then
         call MAapprox_virtual(p,vflav,virtual)
         return
      endif

c     MAflowcheck 1 : exact amplitude evaluated at Q = sqrt(s_bbH),
c     run to muR with its own poles, against the direct muR result.
      if(powheginput("#MAflowcheck").eq.1) call MA_flowcheck(p,vflav)

      call recola_virtual(p,vflav,virtual)

      end


c=======================================================================
c     Running of a one-loop virtual (BLHA finite part, POWHEG units)
c     from the scale Q to st_muren2, as in ttH_NLO:
c        V(muR) = V(Q) + c1(Q) L + c2 L^2/2 + n b0 L B ,
c        L = log(muR^2/Q^2), b0 = (11 CA - 2 nf)/6, n = res_powst,
c     with c1, c2 the IR poles of the MASSIVE amplitude at Q. The HEFT
c     Wilson log (wilsonlog, recola_virtual) is not part of the Recola
c     amplitude: it is removed at Q and put back at muR.
c=======================================================================
      subroutine MA_runQtomuR(Vq,B,c1,c2,Q2,virtual)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_res.h'
      include 'PhysPars.h'
      real * 8, intent(in)  :: Vq,B,c1,c2,Q2
      real * 8, intent(out) :: virtual
      real * 8 L,b0
      real * 8 powheginput
      external powheginput

      L  = dlog(st_muren2/Q2)
      b0 = (11d0*CA - 2d0*st_nlight)/6d0
      virtual = Vq + c1*L + c2*L**2/2d0 + res_powst*b0*L*B
      if(powheginput("#wilsonlog").ne.0d0) virtual = virtual
     $     + 2d0/3d0*(dlog(st_muren2/ph_tmass**2)
     $               -dlog(Q2/ph_tmass**2))*B

      end


      subroutine MA_flowcheck(p,vflav)
      use recola_powheg, only: recola_born,recola_virtual,
     $     recola_virtual_poles
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      integer, parameter :: nlegs=nlegbornexternal
      real * 8, intent(in) :: p(0:3,nlegs)
      integer,  intent(in) :: vflav(nlegs)
      real * 8 B,Vmu,Vq,Vrun,c1,c2,Q2,pbbH(0:3),save2
      real * 8 bjk(nlegs,nlegs),bmn(0:3,0:3,nlegs)
      integer, save :: npt = 0
      real * 8 sq
      external sq

      pbbH = p(:,3) + p(:,4) + p(:,5)
      Q2 = sq(pbbH)
      call recola_born(p,vflav,B,bjk,bmn)
      call recola_virtual(p,vflav,Vmu)
      save2 = st_muren2
      st_muren2 = Q2
      call recola_virtual(p,vflav,Vq)
      call recola_virtual_poles(p,vflav,c1,c2)
      st_muren2 = save2
      call MA_runQtomuR(Vq,B,c1,c2,Q2,Vrun)
      npt = npt + 1
      write(*,'(a,2i3,f9.2,f9.2,4es15.7,es11.3)') ' FLOWCHECK ',
     $     vflav(1),vflav(2),dsqrt(st_muren2),dsqrt(Q2),
     $     Vmu/B,Vq/B,Vrun/B,c1/B,(Vrun-Vmu)/B/dlog(st_muren2/Q2)
      if(npt.ge.40) stop

      end


c=======================================================================
c     MAapprox : massification of the one-loop amplitude,
c
c         M1^(m)  =  F1 * M0^(0)  +  M1^(0)                       (*)
c
c     with the massless amplitudes evaluated on the massless image of
c     the phase-space point and F1 the one-loop Mitov-Moch constant.
c     As in ttH_NLO:
c       - everything is evaluated at Q = eta*sqrt(s_bbH) of the massive
c         point, with POWHEG's alpha_s;
c       - the massified M1M0/M0 is multiplied by the MASSIVE Born;
c       - the result is run from Q to muR with the exact one-loop RGE,
c         using the IR poles of the MASSIVE amplitude (MA_runQtomuR).
c     The massless image is g g / q q~ -> H s s~ (d d~ if the initial
c     state is s s~), registered in the same Recola session next to the
c     massive processes, with MS = 0 and the massive b in the loops.
c     Checked against the exact amplitude: MAflowcheck 1 (flow alone)
c     and MAdebug 2 (massified vs exact, same point, fort.78).
c
c     powheg.input keys:
c        MAapprox         1 to activate
c        MAapprox_mapping 1 MATRIX maps (default), 0 MiNNLOPS maps,
c                         2 lambda family (map_lambda_bbH)
c        lambdascvar      lambda of the family, > 0 (default 1)
c        etascfact        eta in the massification scale Q = eta*sqrt(Q^2)
c        MAapprox_scale   0 massify at muR instead of Q (no running)
c        MAdebug          2 dump massified and exact virtual to fort.78
c=======================================================================
      subroutine MAapprox_virtual(p,vflav,virtual)
      use recola_powheg, only: recola_born,recola_virtual,
     $     recola_virtual_poles
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_res.h'
      include 'PhysPars.h'
      include 'pwhg_MAeta.h'
      integer, parameter :: nlegs=nlegbornexternal
      real * 8, intent(in)  :: p(0:3,nlegs)
      integer,  intent(in)  :: vflav(nlegs)
      real * 8, intent(out) :: virtual
      real * 8 p0(0:3,nlegs),pc(0:3,nlegs),pbbH(0:3)
      real * 8 Bm,c1,c2,F0,F1,VQ,Q2,Qscale,st_muren2_save
      real * 8 Vex,lam,dmax
      real * 8 lamscan(4),Vlam(4)
      data lamscan/0.25d0,0.5d0,1d0,2d0/
      integer vflav0(nlegs),imap,iqim,il,ipart
      logical gg
      real * 8 powheginput
      external powheginput
      real * 8 sq
      external sq

      imap = nint(powheginput("#MAapprox_mapping"))
      if(imap.lt.0) imap = 1
      gg = vflav(1).eq.0 .and. vflav(2).eq.0

c     massive -> massless mapping (quarks sit at legs 4,5)
      call MA_map(p,vflav,imap,p0)

c     the massless image is the s s~ process (MS = 0, b massive in the
c     loops, so the loop content is the one of the exact amplitude)
c     (d d~ for an s s~ initial state, see recola_generate_process)
      iqim = 3
      if(abs(vflav(1)).eq.3) iqim = 1
      vflav0 = vflav
      vflav0(4) = sign(iqim,vflav(4))
      vflav0(5) = sign(iqim,vflav(5))

c     massification scale Q = eta*sqrt(s_bbH) of the MASSIVE point, as
c     in ttH_NLO (scale_M1M0)
      pbbH = p(0:3,3) + p(0:3,4) + p(0:3,5)
      Qscale = ma_eta*dsqrt(sq(pbbH))
c     MAapprox_scale 0 : massify directly at muR (no running)
      if(nint(powheginput("#MAapprox_scale")).eq.0)
     $     Qscale = dsqrt(st_muren2)
      Q2 = Qscale**2

c     at Q: massive Born and massive poles
      st_muren2_save = st_muren2
      st_muren2 = Q2
      call recola_born(p,vflav,Bm)
      call recola_virtual_poles(p,vflav,c1,c2)
      st_muren2 = st_muren2_save

      call massif_factors(ph_bmass,Qscale,vflav(1),F0,F1)

c     (*) massified M1M0 at Q, normalised to the massive Born, then the
c     exact RGE flow Q -> muR with the massive poles
      call MA_image(p0,vflav0,Q2,Bm,F1,VQ)
      call MA_runQtomuR(VQ,Bm,c1,c2,Q2,virtual)

c     MAdebug 2 : massified vs exact at muR on the same point -> fort.78
c     MAdebug 3 : as 2, plus the lambda family of mappings -> fort.79
      if(powheginput("#MAdebug").ge.2) then
         call recola_virtual(p,vflav,Vex)
         write(78,'(2i4,3e16.8,3e18.10)') vflav(1),vflav(2),
     $        dsqrt(p(1,4)**2+p(2,4)**2),dsqrt(p(1,5)**2+p(2,5)**2),
     $        Qscale,Bm,Vex,virtual
      endif
      if(powheginput("#MAdebug").eq.3) then
c        lambda scan in the frame of the channel (gg lab, q qbar
c        partonic); dmax checks the lambda map against p0 at ma_lam
         ipart = 1
         if(gg) ipart = 0
         call map_lambda_bbH(p,pc,ma_lam,ipart)
         dmax = maxval(abs(pc-p0))
         do il=1,4
            call map_lambda_bbH(p,pc,lamscan(il),ipart)
            call MA_image(pc,vflav0,Q2,Bm,F1,VQ)
            call MA_runQtomuR(VQ,Bm,c1,c2,Q2,Vlam(il))
         enddo
         write(79,'(2i4,2e16.8,e11.3,6e18.10)') vflav(1),vflav(2),
     $        dsqrt(p(1,4)**2+p(2,4)**2),dsqrt(p(1,5)**2+p(2,5)**2),
     $        dmax,Bm,Vex,Vlam
      endif

      end


c=======================================================================
c     Massless image at Q2: Born and virtual of the b -> s process on
c     p0, put in the Neubert/MSbar scheme of F1, and the massified
c     M1M0 normalised to the massive Born Bm,
c        VQ = ( 2*F1 + V0/B0 ) * Bm .
c=======================================================================
      subroutine MA_image(p0,vflav0,Q2,Bm,F1,VQ)
      use recola_powheg, only: recola_born,recola_virtual
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      integer, parameter :: nlegs=nlegbornexternal
      real * 8, intent(in)  :: p0(0:3,nlegs),Q2,Bm,F1
      integer,  intent(in)  :: vflav0(nlegs)
      real * 8, intent(out) :: VQ
      real * 8 B0,V0,st_muren2_save

      st_muren2_save = st_muren2
      st_muren2 = Q2
      call recola_born(p0,vflav0,B0)
      call recola_virtual(p0,vflav0,V0)
      st_muren2 = st_muren2_save

c     Move the Recola massless virtual to the Neubert/MSbar scheme that
c     F1 is written in. Recola returns the BLHA finite part (delta_IR=0,
c     delta_IR2 = zeta2); the two differ by c2*zeta2/2 on each side and
c     the massive and massless double poles differ by exactly 2*C_F,
c     so the net shift of the massless input is + C_F*pi^2/6 * B0.
      V0 = V0 + CF*Pi**2/6d0*B0
      VQ = (2d0*F1 + V0/B0)*Bm

      end


c=======================================================================
c     Massive -> massless mapping for MAapprox.
c        imap 0 : MiNNLOPS maps, 1 : MATRIX maps (gg map 2, q qbar map 1)
c        imap 2 : lambda family (map_lambda_bbH), lambda = lambdascvar
c                 (ma_lam, central 1) for both channels: gg in the lab
c                 frame (lambda = 1 is MATRIX map 2), q qbar in the
c                 partonic frame
c=======================================================================
      subroutine MA_map(p,vflav,imap,p0)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_MAeta.h'
      integer, parameter :: nlegs=nlegbornexternal
      real * 8, intent(in)  :: p(0:3,nlegs)
      integer,  intent(in)  :: vflav(nlegs),imap
      real * 8, intent(out) :: p0(0:3,nlegs)
      logical gg

      gg = vflav(1).eq.0 .and. vflav(2).eq.0
      if(imap.eq.0) then
         if(gg) then
            call map_massive_to_massless_momenta_3(p,p0)
         else
            call map_massive_to_massless_momenta_1(p,p0)
         endif
      elseif(imap.eq.1) then
         if(gg) then
            call mapping2_momenta_massiveTOmassless_bbH(p,p0)
         else
            call mapping1new_momenta_massiveTOmassless_bbH(p,p0)
         endif
      else
         if(gg) then
            call map_lambda_bbH(p,p0,ma_lam,0)
         else
            call map_lambda_bbH(p,p0,ma_lam,1)
         endif
      endif

      end


c=======================================================================
c     One-parameter family of massive -> massless maps. In the chosen
c     frame (ipart 0: lab, 1: partonic c.m.) the Higgs is kept and each
c     b quark keeps p_z and its azimuth, while the share lambda of m_b^2
c     goes into the transverse momentum:
c        pT'^2 = pT^2 + lambda m^2 ,  E'^2 = E^2 + (lambda-1) m^2 .
c     The initial states are rebuilt back-to-back in the rest frame of
c     the final state, along the original beam direction (as map 2).
c        gg, lab, lambda = 1      : MATRIX map 2 (E and p_z kept)
c        q qbar, partonic, lambda = 0 : MATRIX map 1 (3-momentum kept)
c     lambda must be > 0: then pT' >= sqrt(lambda) m_b and the massless
c     b quarks stay away from the beam-collinear singularity of the
c     massless amplitude; lambda = 0 (map 1) is IR unsafe.
c=======================================================================
      subroutine map_lambda_bbH(pin,pout,lam,ipart)
      implicit none
      include 'nlegborn.h'
      integer, parameter :: nl=nlegbornexternal
      real * 8, intent(in)  :: pin(0:3,nl),lam
      integer,  intent(in)  :: ipart
      real * 8, intent(out) :: pout(0:3,nl)
      real * 8 w(0:3,nl),wo(0:3,nl),u(0:3),pf(0:3),pp(0:3)
      real * 8 Qp,Qf,m2,pT,pTn,un
      integer i,k
      real * 8 sq
      external sq

      if(lam.le.0d0) then
         write(*,*) ' map_lambda_bbH: lambda must be > 0 (IR safety)'
         call exit(-1)
      endif
      pp = pin(:,1)+pin(:,2)
      Qp = dsqrt(sq(pp))
      if(ipart.eq.1) then
         do i=1,nl
            call boost_recola_bbH(.true.,Qp,pp,pin(:,i),w(:,i))
         enddo
      else
         w = pin
      endif

      wo = 0d0
      wo(:,3) = w(:,3)
      do k=4,5
         m2  = sq(w(:,k))
         pT  = dsqrt(w(1,k)**2+w(2,k)**2)
         pTn = dsqrt(max(pT**2+lam*m2,0d0))
         wo(1,k) = w(1,k)/pT*pTn
         wo(2,k) = w(2,k)/pT*pTn
         wo(3,k) = w(3,k)
         wo(0,k) = dsqrt(pTn**2+w(3,k)**2)
      enddo

      pf = wo(:,3)+wo(:,4)+wo(:,5)
      Qf = dsqrt(sq(pf))
      call boost_recola_bbH(.true.,Qf,pf,w(:,1),u)
      un = dsqrt(u(1)**2+u(2)**2+u(3)**2)
      u(0)   = Qf/2d0
      u(1:3) = u(1:3)/un*Qf/2d0
      call boost_recola_bbH(.false.,Qf,pf,u,wo(:,1))
      u(1:3) = -u(1:3)
      call boost_recola_bbH(.false.,Qf,pf,u,wo(:,2))

      if(ipart.eq.1) then
         do i=1,nl
            call boost_recola_bbH(.false.,Qp,pp,wo(:,i),pout(:,i))
         enddo
      else
         pout = wo
      endif

      end


c=======================================================================
c     One-loop Mitov-Moch massification constant for the TWO massive
c     bottom legs, in units of as/(2pi). Same form as the ttH_NLO study
c     (MiNNLOPS_res/ttH_NLO/virtual.f, massif_factors),
c
c        2*F1 = C_F ( L^2 + L + 4 + pi^2/6 ) ,   L = log(mu^2/m^2)
c
c     which is the Neubert/MSbar-scheme constant: it is what the
c     massification gives once BOTH amplitudes are in that scheme.
c     MAapprox_virtual therefore shifts the Recola massless virtual by
c     C_F pi^2/6 * B0 before using this F1 (see there).
c
c     The lmt terms below account for the bottom being a massless
c     flavour in the loops of the massless amplitude (Recola is
c     generated with m_b = 0, so n_f = 5) while the exact calculation
c     has it decoupled (dZgs_QCD2 = Nf4). Measured in ma-test/:
c     the leftover is exactly (2/3) L for gg and (4/3) L for q qbar,
c     with a vanishing constant -- the channel difference being the
c     "-2*lmt/3" that ttH_NLO already carries for fullflav != 0.
c=======================================================================
      subroutine massif_factors(mass,scale,iflav,F0,F1)
      implicit none
      include 'pwhg_math.h'
      double precision mass, scale, F0, F1, lmt
      integer iflav

      lmt = dlog(mass/scale)
      F0 = 1d0
      F1 = 2d0*CF + (CF*Pi**2)/12d0 - CF*lmt + 2d0*CF*lmt**2
c     No loop-content terms: the massless image (b -> s) keeps the
c     massive b in the loops, as the exact amplitude does. With m_b = 0
c     in the loops (old setup) one needed -2*lmt/3, and another -2*lmt/3
c     for q qbar (ma-test/RESULT.txt).

      end


c=======================================================================
c     ColombaPS : Born and virtual on a fixed massless phase-space point
c     (sqrt(s) = 1000 GeV, m_b = 0, m_H = 125 GeV), for comparison with
c     an external calculation. Runs once and stops.
c
c     Called from init_processes, before the production processes are
c     generated: Recola's reset_recola_rcl does not fully restore the
c     model state, and re-using the session shifts the amplitude by a
c     few per mille.
c     Run in the on-shell scheme with vanishing W/Z/H/t widths.
c     The EW input is traded for a fixed vev: the amplitude carries one
c     power of 1/vev (QED power 1) and Recola builds vev = 2 MW sw/ee,
c     so aEW is chosen to reproduce the requested value. The W/Z/H/t
c     widths are zeroed so that sw, and hence vev, come out real.
c
c     powheg.input keys:
c        ColombaPS        1 to activate
c        ColombaPS_mtop   top mass (default 173.2)
c        ColombaPS_vev    Higgs vev (default 246.220569073)
c        ColombaPS_muren  renormalisation scale (default m_t)
c        ColombaPS_alphas a_s at that scale (default 0.118)
c        ColombaPS_chan   0 = gg (default), 1 = d d~
c=======================================================================
      subroutine ColombaPS_virtual
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_res.h'
      include 'PhysPars.h'
      integer, parameter :: nlegs=nlegbornexternal
      real * 8 p0(0:3,nlegs),pbbH(0:3)
      integer vflav(nlegs)
      real * 8 born0,virt0,Q2,Qscale,as0,vev,sw2
      complex * 16 cvev
      character * 100 proc
      integer i,inf,ichan
      real * 8 powheginput,pwhg_alphas
      external powheginput,pwhg_alphas
      real * 8 sq
      external sq

c     Reference point, reordered into POWHEG's i1 i2 -> H b b~ layout
c     (the quoted list is g g -> b~ b H).
      p0(:,1) = (/ 500.00000000000000000d0,   0.00000000000000000d0,
     $               0.00000000000000000d0, 500.00000000000000000d0 /)
      p0(:,2) = (/ 500.00000000000000000d0,   0.00000000000000000d0,
     $               0.00000000000000000d0,-500.00000000000000000d0 /)
c     H
      p0(:,3) = (/ 211.07738134101575156d0,-144.92833555604772755d0,
     $             -30.63965724497451859d0,  83.58020023094401552d0 /)
c     b
      p0(:,4) = (/ 349.14239937112841972d0, -17.57847102882445967d0,
     $            -333.45082223855939674d0, 101.99000707591839898d0 /)
c     b~
      p0(:,5) = (/ 439.78021928785585715d0, 162.50680658487215169d0,
     $             364.09047948353401125d0,-185.57020730686238608d0 /)

      ichan = nint(powheginput("#ColombaPS_chan"))
      if(ichan.eq.1) then
         vflav = (/ 1, -1, 25, 5, -5 /)
      else
         vflav = (/ 0,  0, 25, 5, -5 /)
      endif

      pbbH = p0(0:3,3) + p0(0:3,4) + p0(0:3,5)
      Q2 = sq(pbbH)

c     top mass, and the renormalisation scale which defaults to it
      ph_tmass = powheginput("#ColombaPS_mtop")
      if(ph_tmass.le.0d0) ph_tmass = 173.2d0
      Qscale = powheginput("#ColombaPS_muren")
      if(Qscale.le.0d0) Qscale = ph_tmass
      as0 = powheginput("#ColombaPS_alphas")
      if(as0.le.0d0) as0 = 0.118d0

c     zero widths and the on-shell (not complex-mass) scheme, so that
c     sw and hence vev are real. Then invert
c        vev = 2 MW sw/ee ,  ee = 2 sqrt(pi aEW)
      ph_Wwidth = 0d0
      ph_Zwidth = 0d0
      ph_Hwidth = 0d0
      ph_twidth = 0d0
      vev = powheginput("#ColombaPS_vev")
      if(vev.le.0d0) vev = 246.220569073d0
      sw2 = 1d0 - ph_Wmass**2/ph_Zmass**2
      ph_alphaem = ph_Wmass**2*sw2/(pi*vev**2)

c     always n_f = 5 for this check against Bayu
      inf = 5
      call massless_recola_amp(vflav,p0,Qscale,as0,inf,
     $     'output_cll/InfOut_ColombaPS.rcl',.true.,.false.,proc,
     $     born0,virt0,cvev)

      write(*,*) ''
      write(*,*) '=================== ColombaPS ==================='
      write(*,*) ' process           : ',trim(proc)
      write(*,'(a,5i5)') '  flavour channel   : ',vflav
      write(*,*) ' m_H               : ',ph_Hmass
      write(*,*) ' m_t               : ',ph_tmass
      write(*,*) ' vev (requested)   : ',vev
      write(*,*) ' vev (from Recola) : ',cvev
      write(*,*) ' alpha_em          : ',ph_alphaem
      write(*,*) ' sqrt(Q^2) bbH     : ',dsqrt(Q2)
      write(*,*) ' renorm. scale     : ',Qscale
      write(*,*) ' alpha_s(scale)    : ',as0
      write(*,*) ' nf (dZgs_QCD2)    : ',inf
      write(*,*) ' -- phase-space point (POWHEG order: i1 i2 H b b~) --'
      do i=1,nlegs
         write(*,'(i3,4(1x,e22.15),2x,a,e12.5)') i,p0(0:3,i),
     $        'p^2=',sq(p0(0:3,i))
      enddo
      write(*,*) ' Born, O(as^',res_powst,')   : ',born0
      write(*,*) ' Virt, O(as^',res_powst+1,')   : ',virt0
      write(*,*) ' V/(as/2pi) [POWHEG] : ',virt0/(as0/2d0/pi)
      write(*,*) ' V/B (relative)      : ',virt0/born0
      write(*,*) '================================================='
      write(*,*) ''
      call flush(6)
      stop

      end


c=======================================================================
c     Born and virtual for a MASSLESS-b point. Recola fixes the masses
c     at process generation, so this replaces the massive processes and
c     can only be used once, at the end of a run.
c     born0/virt0 come back with a_s(Qscale) restored, in POWHEG's
c     normalisation for the flavour structure vflav.
c=======================================================================
      subroutine massless_recola_amp(vflav,p0,Qscale,as0,inf,logfile,
     $     onshell,doreset,proc,born0,virt0,cvev)
      use recola
      use recola_powheg, only: flav_to_string
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_res.h'
      include 'PhysPars.h'
      integer, parameter :: nlegs=nlegbornexternal
      integer,  intent(in)    :: vflav(nlegs)
      real * 8, intent(in)    :: p0(0:3,nlegs),Qscale,as0
      integer,  intent(inout) :: inf
      character(*), intent(in):: logfile
c     doreset: drop any processes already generated. Recola's reset does
c     not fully restore the model state, so a check that needs exact
c     numbers must run BEFORE the production processes are generated and
c     pass doreset = .false.
      logical,  intent(in)    :: onshell,doreset
      character * 100, intent(out) :: proc
      real * 8, intent(out)   :: born0,virt0
      complex * 16, intent(out) :: cvev
      real * 8 cgghfin
      character * 4 nfstr
      real * 8 powheginput
      external powheginput

      if(doreset) call reset_recola_rcl

      call set_output_file_rcl(trim(logfile))
      call set_print_level_squared_amplitude_rcl(0)
      call set_print_level_parameters_rcl(2)

      call set_parameter_rcl('MZ' ,dcmplx(ph_Zmass))
      call set_parameter_rcl('WZ' ,dcmplx(ph_Zwidth))
      call set_parameter_rcl('MW' ,dcmplx(ph_Wmass))
      call set_parameter_rcl('WW' ,dcmplx(ph_Wwidth))
      call set_parameter_rcl('MH' ,dcmplx(ph_Hmass))
      call set_parameter_rcl('WH' ,dcmplx(ph_Hwidth))
      call set_parameter_rcl('MT' ,dcmplx(ph_tmass))
      call set_parameter_rcl('WT' ,dcmplx(ph_twidth))
c     massless bottom, vanishing bottom Yukawa
      call set_parameter_rcl('MB' ,dcmplx(0d0))
      call set_parameter_rcl('WB' ,dcmplx(0d0))
      call set_parameter_rcl('ymb',dcmplx(0d0))
      call set_parameter_rcl('aEW',dcmplx(ph_alphaem))

      if(onshell .or. powheginput('#complexscheme').eq.0)then
         call set_on_shell_scheme_rcl
      else
         call set_complex_mass_scheme_rcl
      endif
      call use_dim_reg_soft_rcl
      call set_delta_ir_rcl(0d0,pi**2/6d0)
      call set_dynamic_settings_rcl(1)

c     HEFT Wilson coefficient, as in init_heft_parameters
      cgghfin = powheginput("#cgghfin")
      if(cgghfin.lt.0d0) cgghfin = 0d0
      call set_parameter_rcl('cgghfin',dcmplx(cgghfin))

c     Recola normalises the coupling-power coefficients to aS = 1, so
c     the aS parameter must be left at 1 (as recola_alphas does in the
c     production run) and a_s(Q) restored by hand below.
      call set_parameter_rcl('aS',dcmplx(1d0))
c     muMS enters the a_s counterterm as beta0*log(muUV^2/muMS^2)
      call set_mu_ms_rcl(Qscale)
      if(inf.lt.3 .or. inf.gt.6) inf = 5
      write(nfstr,'(a2,i1)') 'Nf',inf
      call set_renoscheme_rcl('dZgs_QCD2',trim(nfstr))

      call flav_to_string(vflav,proc)
      call define_process_rcl(1,proc(:),'NLO')
      call unselect_all_powers_BornAmpl_rcl(1)
      call select_power_BornAmpl_rcl(1,'QCD',res_powst)
      call select_power_BornAmpl_rcl(1,'QED',res_powew)
      call unselect_all_powers_LoopAmpl_rcl(1)
      call select_power_LoopAmpl_rcl(1,'QCD',res_powst+2)
      call select_power_LoopAmpl_rcl(1,'QED',res_powew)
      call generate_processes_rcl

      call get_parameter_rcl('vev',cvev)

      call set_mu_ir_rcl(Qscale)
      call set_mu_uv_rcl(Qscale)
      call set_mu_ms_rcl(Qscale)

      call compute_process_rcl(1,p0,'NLO')
      call get_squared_amplitude_rcl(1,
     $     [2*res_powst,2*res_powew],'LO',born0)
      call get_squared_amplitude_rcl(1,
     $     [2*(res_powst+1),2*res_powew],'NLO',virt0)

c     restore a_s(Q): the virtual carries one power more than the Born
      born0 = born0*as0**res_powst
      virt0 = virt0*as0**(res_powst+1)

      end


c=======================================================================
c     Massive -> massless mappings from MiNNLOPS_res/ttH_NLO/virtual.f.
c     The massive quarks are expected at positions 4 and 5.
c=======================================================================

      double precision function sp(p1,p2)
      implicit none
      double precision p1(0:3), p2(0:3)
      sp = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)
      end

      double precision function sp3(p1,p2)
      implicit none
      double precision p1(3), p2(3)
      sp3 = p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)
      end

      double precision function sq(p1)
      implicit none
      double precision p1(0:3)
      double precision, external :: sp
      sq = sp(p1,p1)
      end

      subroutine boost_momentum(p_direction, p_in, p_out)
      implicit none
      double precision p_direction(0:3),p_in(0:3),p_out(0:3)
      double precision beta(3),r(3)
      double precision bn,gam,gamm1
      double precision, external :: sp,sq,sp3

      beta = -p_direction(1:3)/p_direction(0)
      r = p_in(1:3)
      bn = sp3(beta,beta)
      gam = 1/sqrt(1-bn)
      p_out(0) = gam * ( p_in(0) - sp3(r,beta) )
      if (abs(bn) > 1.d-6) then
         gamm1 = (gam-1)/bn
      else
         gamm1 = 1d0/2d0 + 3d0/8d0*bn
      endif
      p_out(1:3) = r + ( gamm1 * sp3(r,beta) - gam * p_in(0) ) * beta
      end

c     map p_in onto the light cone using the massive reference qref_in
      subroutine map_to_lightcone(p_in, qref_in, p_out)
      implicit none
      double precision p_in(0:3), qref_in(0:3), p_out(0:3)
      double precision msq, qsq, sppq
      double precision, external :: sp,sq

      qsq = sq(qref_in)
      msq = sq(p_in)
      sppq = sp(p_in, qref_in)
      if (qsq >= 0) then
         print*, "ERROR: map_to_lightcone not implemented, qsq = ", qsq
         stop
      endif
      p_out = p_in + ((sqrt(1 - (msq*qsq)/sppq**2) - 1)*sppq/qsq)*qref_in
      if (abs(sq(p_out)) > 1d-5) then
         print*, "ERROR: projected momentum is not massless: ",sq(p_out)
         stop
      endif
      end

c     simultaneous lightcone decomposition of two massive momenta
      subroutine simultaneous_lightcone_decompose(p1_in,p2_in,
     $     p1_out,p2_out)
      implicit none
      double precision p1_in(0:3), p2_in(0:3), p1_out(0:3), p2_out(0:3)
      double precision mQ, mQQ, bp, bm, beta, m3, m4
      double precision, external :: sp,sq

      mQ = sqrt(sq(p1_in))
      mQQ = sqrt(sq(p1_in+p2_in))
      beta = sqrt(1d0 - 4d0 * mQ**2/mQQ**2)
      bp = (1d0 + beta)/(2d0 * beta)
      bm = (1d0 - beta)/(2d0 * beta)
      p1_out = bp * p1_in - bm * p2_in
      p2_out = bp * p2_in - bm * p1_in
      m3 = abs(sq(p1_out))
      m4 = abs(sq(p2_out))
      if(m3.gt.1.d-5 .or. m4.gt.1.d-5) then
         print*, "ERROR: projected momenta are not massless"
         print*, "m3, m4 = ",m3, m4
         stop
      endif
      end

c     MiNNLOPS map 1: preserves p_b + p_bbar; the quarks may end up
c     collinear to the beams.
      subroutine map_massive_to_massless_momenta_1(p_in,p_out)
      implicit none
      include 'nlegborn.h'
      double precision p_in(0:3,nlegbornexternal)
      double precision p_out(0:3,nlegbornexternal), qq(0:3)
      integer ii

      p_out = p_in
      call simultaneous_lightcone_decompose(p_in(0:3,4), p_in(0:3,5),
     $     p_out(0:3,4), p_out(0:3,5))
      qq = -p_out(0:3,1)-p_out(0:3,2)
      do ii = 3, nlegbornexternal
         qq = qq + p_out(0:3,ii)
      end do
      if (maxval(abs(qq)) > 1d-6) then
         print*, "ERROR: momentum conservation violated: ", qq
         print*, '--- continue but be careful ---'
      endif
      end

c     MiNNLOPS map 3: the quarks stay off the beam axis, the residual
c     momentum is reabsorbed in the initial state. Assumes the partonic
c     CM frame, where POWHEG builds the Born phase space.
      subroutine map_massive_to_massless_momenta_3(p_in,p_out)
      implicit none
      include 'nlegborn.h'
      double precision p_in(0:3,nlegbornexternal)
      double precision p_out(0:3,nlegbornexternal)
      double precision n1(0:3), n2(0:3), p1(0:3), p2(0:3), x
      double precision, external :: sp,sq
      integer b3i,b4i

      p_out = p_in
      b3i = 4
      b4i = 5
      p1 = p_in(0:3,1)
      p2 = p_in(0:3,2)

c     seed momenta transverse to p1,p2
      n1 = p_in(0:3,b3i)
      n2 = p_in(0:3,b4i)
      n1 = n1 - (sp(p1,n1)/sp(p1,p2))*p2 - (sp(p2,n1)/sp(p1,p2))*p1
      n2 = n2 - (sp(p1,n2)/sp(p1,p2))*p2 - (sp(p2,n2)/sp(p1,p2))*p1

      call map_to_lightcone(p_in(0:3,b3i), n1, p_out(0:3,b3i))
      call map_to_lightcone(p_in(0:3,b4i), n2, p_out(0:3,b4i))

c     residual momentum reabsorbed into the initial state
      n1 = -p1 - p2 - p_out(0:3,b3i) - p_out(0:3,b4i)
     $     + p_in(0:3,b3i) + p_in(0:3,b4i)
      n2 = p1 + p2
      if(abs(n2(1)) + abs(n2(2)) + abs(n2(3)) > 1d-8) then
         print*, "ERROR: not in the rest frame of initial momenta!"
         stop
      endif
      call boost_momentum(n1, p1, p_out(0:3,1))
      call boost_momentum(n1, p2, p_out(0:3,2))
      x = sqrt(sq(n1)/sq(n2))
      p_out(:,1) = x * p_out(:,1)
      p_out(:,2) = x * p_out(:,2)
      if (abs(sq(p_out(0:3,1))) + abs(sq(p_out(0:3,2))) > 1d-5) then
         print*, "ERROR: projected momentum is not massless"
         stop
      endif
      n1 = -p_out(0:3,1)-p_out(0:3,2)
      do b3i = 3, nlegbornexternal
         n1 = n1 + p_out(0:3,b3i)
      end do
      if (maxval(abs(n1)) > 1d-6) then
         print*, "ERROR: momentum conservation violated: ", n1
         print*, '--- continue but be careful ---'
      endif
      end

c     MATRIX map 1 (q qbar in ttH), re-indexed for i1 i2 -> H b b~: the
c     Higgs is untouched, the quarks are put on the light cone in the
c     partonic CM.
      subroutine mapping1new_momenta_massiveTOmassless_bbH(pin,pout)
      implicit none
      include 'nlegborn.h'
      integer i, nl
      parameter (nl=nlegbornexternal)
      double precision pin(0:3,nl), pout(0:3,nl)
      double precision pin_part(0:3,nl), pout_part(0:3,nl)
      double precision Q
      double precision, external :: sq

      pout = 0d0
      pout_part = 0d0
c     the Higgs momentum is kept unchanged
      pout(:,3) = pin(:,3)

c     boost to the partonic CM frame
      Q = sqrt(sq(pin(:,1)+pin(:,2)))
      do i=1,nl
         call boost_recola_bbH(.true.,Q,pin(:,1)+pin(:,2),
     $        pin(:,i),pin_part(:,i))
      enddo
      pout_part(:,3) = pin_part(:,3)

      pout_part(:,4) = pin_part(:,4)
      pout_part(0,4) = dsqrt(pin_part(1,4)**2+pin_part(2,4)**2
     $     +pin_part(3,4)**2)
      pout_part(:,5) = pin_part(:,5)
      pout_part(0,5) = dsqrt(pin_part(1,5)**2+pin_part(2,5)**2
     $     +pin_part(3,5)**2)

      pout_part(0:3,1) = (/ 1d0, 0d0, 0d0, 1d0 /)
      pout_part(0:3,1) = (pout_part(0,3)+pout_part(0,4)+pout_part(0,5)
     $     +pout_part(3,3)+pout_part(3,4)+pout_part(3,5))/2d0
     $     *pout_part(0:3,1)
      pout_part(0:3,2) = (/ 1d0, 0d0, 0d0, -1d0 /)
      pout_part(0:3,2) = (pout_part(0,3)+pout_part(0,4)+pout_part(0,5)
     $     -pout_part(3,3)-pout_part(3,4)-pout_part(3,5))/2d0
     $     *pout_part(0:3,2)

      do i=1,5
         call boost_recola_bbH(.false.,Q,pin(:,1)+pin(:,2),
     $        pout_part(:,i),pout(:,i))
      enddo
      end

c     MATRIX map 2 (gg in ttH), same re-indexing: the quark energy and
c     longitudinal momentum are preserved and p_T is rescaled to put
c     them on shell, avoiding beam-collinear points.
      subroutine mapping2_momenta_massiveTOmassless_bbH(pin,pout)
      implicit none
      include 'nlegborn.h'
      integer i, nu, nl
      parameter (nl=nlegbornexternal)
      double precision pin(0:3,nl), pout(0:3,nl)
      double precision pin_part(0:3,nl), pout_part(0:3,nl)
      double precision Q, E1, p4T, p5T, mQsq
      double precision, external :: sq

      pout = 0d0
      pout_part = 0d0
      mQsq = sq(pin(:,4))
c     the Higgs momentum is kept unchanged
      pout(:,3) = pin(:,3)

      p4T = dsqrt(pin(1,4)**2 + pin(2,4)**2)
      pout(0,4) = pin(0,4)
      pout(3,4) = pin(3,4)
      pout(1,4) = pin(1,4)/p4T * dsqrt(p4T**2+mQsq)
      pout(2,4) = pin(2,4)/p4T * dsqrt(p4T**2+mQsq)

      p5T = dsqrt(pin(1,5)**2 + pin(2,5)**2)
      pout(0,5) = pin(0,5)
      pout(3,5) = pin(3,5)
      pout(1,5) = pin(1,5)/p5T * dsqrt(p5T**2+mQsq)
      pout(2,5) = pin(2,5)/p5T * dsqrt(p5T**2+mQsq)

      Q = sqrt(sq(pout(:,3)+pout(:,4)+pout(:,5)))
      do i=1,2
         call boost_recola_bbH(.true.,Q,pout(:,3)+pout(:,4)+pout(:,5),
     $        pin(:,i),pin_part(:,i))
      enddo
      E1 = dsqrt(pin_part(1,1)**2+pin_part(2,1)**2+pin_part(3,1)**2)
      do nu=1,3
         pout_part(nu,1) = pin_part(nu,1)/E1 * Q/2d0
         pout_part(nu,2) = -pout_part(nu,1)
      enddo
      pout_part(0,1) = Q/2d0
      pout_part(0,2) = Q/2d0
      do i=1,2
         call boost_recola_bbH(.false.,Q,pout(:,3)+pout(:,4)+pout(:,5),
     $        pout_part(:,i),pout(:,i))
      enddo
      end

      subroutine boost_recola_bbH(bool,mass,p1,p_in,p_out)
      implicit none
      double precision mass,p1(0:3),p_in(0:3),p_out(0:3)
      double precision gam,beta(1:3),bdotp,one
      parameter(one=1d0)
      integer j,k,sign
      logical bool
      if(bool) then
         sign = 1
      else
         sign = -1
      endif
      gam=p1(0)/mass
      bdotp=0d0
      do j=1,3
         beta(j) = sign*p1(j)/p1(0)
         bdotp=bdotp+p_in(j)*beta(j)
      enddo
      p_out(0)=gam*(p_in(0)-bdotp)
      do k=1,3
         p_out(k)=p_in(k)+gam*beta(k)*(gam/(gam+one)*bdotp-p_in(0))
      enddo
      end
