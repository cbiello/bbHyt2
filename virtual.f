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

c     MAapprox 1 : massification test, see MAapprox_virtual below.
      if(powheginput("#MAapprox").eq.1) then
c        MAapprox_chan: 0 = gg, 1 = q qbar, absent = first point
         if(powheginput("#MAapprox_chan").lt.0d0 .or.
     $      (powheginput("#MAapprox_chan").eq.0d0 .and.
     $       vflav(1).eq.0 .and. vflav(2).eq.0) .or.
     $      (powheginput("#MAapprox_chan").eq.1d0 .and.
     $       vflav(1).ne.0)) then
            call MAapprox_virtual(p,vflav)
         endif
      endif

      call recola_virtual(p,vflav,virtual)

      end


c=======================================================================
c     MAapprox : one-loop amplitude with a MASSLESS bottom quark, on the
c     massless image of the current point built with the ttH mappings.
c     Recola fixes the masses at process generation, so the massless
c     processes replace the massive ones: this runs once and stops.
c
c     powheg.input keys:
c        MAapprox         1 to activate
c        MAapprox_mapping 1 MATRIX maps (default), 0 MiNNLOPS maps
c        MAapprox_nf      n_f in the a_s renormalisation (default 5)
c        MAapprox_scale   1 evaluate at Q (default), 0 at POWHEG's muR
c        MAapprox_chan    0 = gg, 1 = q qbar, absent = first point
c        etascfact        eta in Q = eta*sqrt(Q^2), as in ttH_NLO
c=======================================================================
      subroutine MAapprox_virtual(p,vflav)
      use recola
      use recola_powheg, only: flav_to_string
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_res.h'
      include 'PhysPars.h'
      integer, parameter :: nlegs=nlegbornexternal
      real * 8, intent(in) :: p(0:3,nlegs)
      integer,  intent(in) :: vflav(nlegs)
      real * 8 p0(0:3,nlegs)
      real * 8 born0,virt0,cgghfin
      real * 8 pbbH(0:3),Q2,Qscale,as0,eta
      character * 100 proc
      integer i,imap,inf,iscale
      character * 4 nfstr
      real * 8 powheginput,pwhg_alphas
      external powheginput,pwhg_alphas
      real * 8 sq
      external sq

      imap = nint(powheginput("#MAapprox_mapping"))
      if(imap.lt.0) imap = 1

c     1) massive -> massless mapping. The Born ordering i1 i2 -> H b b~
c     puts the quarks at 4,5, as the ttH maps expect. As in ttH, gg uses
c     the map keeping the quarks off the beam axis, q qbar the one
c     preserving p_b + p_bbar.
      if(vflav(1).eq.0 .and. vflav(2).eq.0) then
         if(imap.eq.0) then
            call map_massive_to_massless_momenta_3(p,p0)
         else
            call mapping2_momenta_massiveTOmassless_bbH(p,p0)
         endif
      else
         if(imap.eq.0) then
            call map_massive_to_massless_momenta_1(p,p0)
         else
            call mapping1new_momenta_massiveTOmassless_bbH(p,p0)
         endif
      endif

c     2) scale of the massless amplitude, as in ttH_NLO:
c     Q = eta*sqrt((p_H+p_b+p_bbar)^2), the invariant mass of the
c     massless bbH system. The amplitude is renormalised there and
c     evaluated with a_s(Q).
      pbbH = p0(0:3,3) + p0(0:3,4) + p0(0:3,5)
      Q2 = sq(pbbH)
      eta = powheginput("#etascfact")
      if(eta.le.0d0) eta = 1d0
      Qscale = eta*dsqrt(Q2)

      iscale = nint(powheginput("#MAapprox_scale"))
      if(iscale.lt.0) iscale = 1
      if(iscale.eq.0) then
         Qscale = dsqrt(st_muren2)
      endif
      as0 = pwhg_alphas(Qscale**2,st_lambda5MSB,st_nlight)

c     3) restart Recola with m_b = 0 and register this channel only.
c     Settings as in recola_init (HEFT needs set_parameter_rcl).
      call reset_recola_rcl

      call set_output_file_rcl('output_cll/InfOut_MAapprox.rcl')
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

      if(powheginput('#complexscheme').eq.0)then
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
      inf = nint(powheginput("#MAapprox_nf"))
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

      call set_mu_ir_rcl(Qscale)
      call set_mu_uv_rcl(Qscale)
      call set_mu_ms_rcl(Qscale)

c     4) evaluate and report
      call compute_process_rcl(1,p0,'NLO')
      call get_squared_amplitude_rcl(1,
     $     [2*res_powst,2*res_powew],'LO',born0)
      call get_squared_amplitude_rcl(1,
     $     [2*(res_powst+1),2*res_powew],'NLO',virt0)

c     restore a_s(Q): the virtual carries one power more than the Born
      born0 = born0*as0**res_powst
      virt0 = virt0*as0**(res_powst+1)

      write(*,*) ''
      write(*,*) '==================== MAapprox ===================='
      write(*,*) ' mapping           : ',imap,' (1=MATRIX, 0=MiNNLOPS)'
      write(*,*) ' process           : ',trim(proc)
      write(*,'(a,5i5)') '  flavour channel   : ',vflav
      write(*,*) ' POWHEG muR        : ',dsqrt(st_muren2)
      write(*,*) ' POWHEG alpha_s    : ',st_alpha
      write(*,*) ' sqrt(Q^2) bbH     : ',dsqrt(Q2)
      write(*,*) ' eta (etascfact)   : ',eta
      write(*,*) ' massless scale    : ',Qscale
      write(*,*) ' alpha_s(scale)    : ',as0
      write(*,*) ' nf (dZgs_QCD2)    : ',inf
      write(*,*) ' -- massive phase-space point (m_b = ',ph_bmass,') --'
      do i=1,nlegs
         write(*,'(i3,4(1x,e22.15),2x,a,e12.5)') i,p(0:3,i),
     $        'p^2=',sq(p(0:3,i))
      enddo
      write(*,*) ' -- massless image --'
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
