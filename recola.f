
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!                    ===POWHEG+Recola Interface=== 
!                 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc           
      module recola_powheg
      use recola
      implicit none
      
      character*100,dimension(:),allocatable,private::map

      type process_rcl
      integer :: i_mast
      integer, allocatable :: flav(:)
!      character*100 :: name
      integer,allocatable::cols(:,:)
      integer,allocatable::hels(:,:)
      end type
      type(process_rcl), allocatable, save :: processes(:)
      character(len=100) :: model_rcl
      common/physicsmodel/model_rcl

c     Whether the active model file supports Recola's fast a_s
c     rescaling (set_alphas_rcl after generate_processes_rcl).
c     The HEFT model file sets has_feature_mdl('qcd_rescaling')=.false.
c     -- the effective Hgg Wilson coefficient cggh is itself
c     proportional to a_s, so the amplitude is not a plain power of gs
c     and cannot be rescaled a posteriori. With such a model file a_s
c     must be fixed BEFORE the processes are generated, i.e. the run
c     is restricted to a fixed renormalisation scale.
      logical, save :: als_rescaling = .true.
c     Value alpha_s was frozen to (used only when als_rescaling false)
      double precision, save :: als_frozen = -1d0
      
      contains
      
      subroutine recola_init()
      use modelfile, only: has_feature_mdl
      implicit none
      real * 8 powheginput
      external powheginput
      include 'PhysPars.h'
      include 'pwhg_st.h'
      integer, parameter :: dp = kind (23d0) ! double precision
      real(dp), parameter :: pi = 3.141592653589793238462643d0
      character(len=50) :: out_path
c     Recola's convenience setters (set_pole_mass_*_rcl,
c     use_gfermi_scheme_rcl, ...) are implemented only for the SM, THDM,
c     HS and TGC model files. With any other model file they abort with
c        "Model file <name> not supported.
c         Supported models are: SM,THDM,HS,TGC"
c     and parameters must go through the generic set_parameter_rcl
c     interface instead. That applies to SM_ATGC and equally to HEFT,
c     which is the model file bbH_yt2 uses ("SM (QCD) + HEFT").
c     NOTE: this only concerns how the INPUT PARAMETERS are set. The
c     ATGC-specific 'LAM' coupling-order blocks further down stay keyed
c     on the ATGC model name, since HEFT has only QCD and QED orders.
      logical :: generic_par

      ! Retrieve active model                                       
      call get_modelname_rcl(model_rcl)
      generic_par = (index(model_rcl,'ATGC').ne.0) .or.
     $              (index(model_rcl,'HEFT').ne.0)
      als_rescaling = has_feature_mdl('qcd_rescaling')
      if(.not.als_rescaling) then
         write(*,*) ' ############################################'
         write(*,*) ' # Model file "',trim(model_rcl),'"'
         write(*,*) ' # does NOT support qcd_rescaling, so alpha_s'
         write(*,*) ' # cannot be changed after process generation.'
         write(*,*) ' # alpha_s is FROZEN to the value set below and'
         write(*,*) ' # the run is only meaningful with a FIXED'
         write(*,*) ' # renormalisation scale (runningscales 0).'
         write(*,*) ' ############################################'
      endif
      if(index(model_rcl,'ATGC').ne.0)then
         ! See arXiv:1804.01477
         ! atgc_powlamm2 = 0 SM^2
         ! atgc_powlamm2 = 1 SMxEFT6 + EFT6^2
         ! atgc_powlamm2 = 2 SMxEFT8   
         atgc_powlamm2 = powheginput('#atgc_powlamm2')
         if(atgc_powlamm2<0) atgc_powlamm2 = 0
         call recola_atgc
      endif

!      call set_cache_mode_rcl(0) ! need to turn off cache in RCL to avoid segfault

      out_path='output_cll'
!      call set_output_file_rcl('*') ! Standard output
      call set_output_file_rcl(trim(out_path)//'/InfOut.rcl')
      call set_print_level_squared_amplitude_rcl (0) ! By default, amplitude is not printed
      call set_print_level_parameters_rcl (2) !  By default, parameters are not printed

      !call set_crossing_symmetry_rcl(.true./.false.) Switch on(default)/off crossing symmetry

      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Physical parameters           
      if(generic_par)then
         if (ph_Zmass/= 0)call set_parameter_rcl('MZ', dcmplx(ph_Zmass))
         if (ph_Zwidth/= 0)call set_parameter_rcl('WZ', dcmplx(ph_Zwidth))
         if (ph_Wmass/= 0)call set_parameter_rcl('MW', dcmplx(ph_Wmass))
         if (ph_Wwidth/= 0)call set_parameter_rcl('WW', dcmplx(ph_Wwidth))
         if (ph_Zmass/= 0)call set_parameter_rcl('MH', dcmplx(ph_Hmass))
         if (ph_Zwidth/= 0)call set_parameter_rcl('WH', dcmplx(ph_Hwidth))
         if (ph_tmass/= 0)call set_parameter_rcl('MT', dcmplx(ph_tmass))
         if (ph_twidth/= 0)call set_parameter_rcl('WT', dcmplx(ph_twidth))
         if (ph_bmass/= 0)call set_parameter_rcl('MB', dcmplx(ph_bmass))
         if (ph_bwidth/= 0)call set_parameter_rcl('WB', dcmplx(ph_bwidth))
         if (ph_cmass/= 0)call set_parameter_rcl('MC', dcmplx(ph_cmass))
         if (ph_cwidth/= 0)call set_parameter_rcl('WC', dcmplx(ph_cwidth))
         if (ph_taumass/= 0)call set_parameter_rcl('MTA', dcmplx(ph_taumass))
         if (ph_tauwidth/= 0)call set_parameter_rcl('WTA', dcmplx(ph_tauwidth))
         if (ph_mumass/= 0)call set_parameter_rcl('MM', dcmplx(ph_mumass))
         if (ph_muwidth/= 0)call set_parameter_rcl('WM', dcmplx(ph_muwidth))
      else
         if ((ph_Zmass/= 0).or.(ph_Zwidth/= 0))call set_pole_mass_Z_rcl(ph_Zmass,ph_Zwidth)
         if ((ph_Wmass/= 0).or.(ph_Wwidth/= 0))call set_pole_mass_W_rcl(ph_Wmass,ph_Wwidth)
         if ((ph_Hmass/= 0).or.(ph_Hwidth/= 0))call set_pole_mass_H_rcl(ph_Hmass,ph_Hwidth)
         if ((ph_tmass/= 0).or.(ph_twidth/= 0))call set_pole_mass_top_rcl(ph_tmass,ph_twidth)
         if ((ph_bmass/= 0).or.(ph_bwidth/= 0))call set_pole_mass_bottom_rcl(ph_bmass,ph_bwidth)
         if ((ph_cmass/= 0).or.(ph_cwidth/= 0))call set_pole_mass_charm_rcl(ph_cmass,ph_cwidth)
         if ((ph_taumass/= 0).or.(ph_tauwidth/= 0))call set_pole_mass_tau_rcl(ph_taumass,ph_tauwidth)
         if ((ph_mumass/= 0).or.(ph_muwidth/= 0))call set_pole_mass_muon_rcl(ph_mumass,ph_muwidth)
      endif

      if(index(model_rcl,'ATGC').ne.0)then
         if(atgc_cwwwl2/=0) call set_parameter_rcl('CWWWL2',dcmplx(atgc_cwwwl2))
         if(atgc_cwl2/=0) call set_parameter_rcl('CWL2',dcmplx(atgc_cwl2))
         if(atgc_cbl2/=0) call set_parameter_rcl('CBL2',dcmplx(atgc_cbl2))
         if(atgc_cpwwwl2/=0) call set_parameter_rcl('CPWWWL2',dcmplx(atgc_cpwwwl2))
         if(atgc_cpwl2/=0) call set_parameter_rcl('CPWL2',dcmplx(atgc_cpwl2))
         if(atgc_cbwl4/=0) call set_parameter_rcl('CBWL4',dcmplx(atgc_cbwl4))
         if(atgc_cwwl4/=0) call set_parameter_rcl('CWWL4',dcmplx(atgc_cwwl4))
         if(atgc_cbbl4/=0) call set_parameter_rcl('CBBL4',dcmplx(atgc_cbbl4))
         if(atgc_cbtwl4/=0) call set_parameter_rcl('CBtWL4',dcmplx(atgc_cbtwl4))
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!   NEEDED ???
      ! Sets a masscut mcut below which massive fermions are treated as light particles
      ! for dimensional and mass regularization
c$$$      call set_light_fermions_rcl(1d-3)

      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! EW scheme
      if(powheginput('#ewscheme').eq.0)then
        !  alpha(0)  scheme (by default alpha set to value of Thomson scattering)
        !  alpha_qed_0     alpha =alpha_qed_0
         if(generic_par)then
            if(ph_alphaem/=0)then
               call set_parameter_rcl('aEW', dcmplx(ph_alphaem))
            else
               call set_parameter_rcl('aEW', dcmplx(1d0/137.035999679d0))
            endif
         else
            if(ph_alphaem/=0)then
               call use_alpha0_scheme_rcl(a=ph_alphaem)
            else
               call use_alpha0_scheme_rcl()
            endif     
         endif
      else if(powheginput('#ewscheme').eq.1.or.powheginput('#ewscheme').lt.0)then
        ! G_F scheme, obtained providing Gmu or alpha directly
        ! Gmu     alpha = sqrt2/pi*Gmu*|MW2*sw2|
         if(generic_par)then
            if(ph_alphaem/=0)then
               call set_parameter_rcl('aEW', dcmplx(ph_alphaem))
            else if(ph_gfermi/=0) then
               ph_alphaem=ph_gfermi/pi*ph_Wmass**2d0*(1d0-ph_Wmass**2d0/ph_Zmass**2d0)*sqrt(2d0)
               call set_parameter_rcl('aEW', dcmplx(ph_alphaem))
            else
               ph_alphaem=(1.16637d-5)/pi*ph_Wmass**2d0*(1d0-ph_Wmass**2d0/ph_Zmass**2d0)*sqrt(2d0)   
               call set_parameter_rcl('aEW', dcmplx(ph_alphaem))
            endif
         else
            if(ph_alphaem/=0)then
               call use_gfermi_scheme_rcl(a=ph_alphaem)
            else if(ph_gfermi/=0) then
               call use_gfermi_scheme_rcl(g=ph_gfermi)
            else
               call use_gfermi_scheme_rcl(g=1.16637d-5)
            endif
         endif
      else if(powheginput('#ewscheme').eq.2)then
         ! alpha(MZ) scheme
         ! alpha_qed_mz    alpha =alpha_qed_mz 
         if(generic_par)then
            if(ph_alphaem/=0)then
               call set_parameter_rcl('aEW', dcmplx(ph_alphaem))
            else
               call set_parameter_rcl('aEW', dcmplx(1d0/128.936d0))
            endif
         else
            if(ph_alphaem/=0)then
               call use_alphaz_scheme_rcl(a=ph_alphaem)
            else
               call use_alphaz_scheme_rcl()
            endif     
         endif
      endif
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Unstable particle renormalization scheme
      if(powheginput('#complexscheme').eq.0)then
         call set_on_shell_scheme_rcl
      else
         call set_complex_mass_scheme_rcl
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! For NLO
      call use_dim_reg_soft_rcl
      ! The double IR-pole DeltaIR2 is set to Zeta(2) = pi^2/6.
      ! In this way the result of the amplitude corresponds to the finite 
      ! part of the amplitude in the conventions of the Binoth-Les Houches
      ! Accord (arXiv:1001.1307 [hep-ph])
      !  call set_delta_ir_rcl (d,d2)
      call set_delta_ir_rcl(0d0,pi**2/6d0)
      ! In addition to alphas also lambda, DeltaUV, muUV, DeltaIR2, 
      ! DeltaIR and muIR can be reset after generate_processes_rcl
      call set_dynamic_settings_rcl(1)




      !!!!!!!debug
      print*, 'I am in recola_init inside recola.f'
      end subroutine recola_init


      subroutine recola_atgc
      implicit none
      include 'PhysPars.h'
      real * 8 powheginput
      external powheginput

      ! Default values from arXiv:1804.01477
      atgc_cwwwl2  = powheginput('#cwwwl2')
      if(atgc_cwwwl2<0d0) atgc_cwwwl2 = 3e-6 
      atgc_cwl2    = powheginput('#cwl2')
      if(atgc_cwl2<0d0) atgc_cwl2 = 3e-6
      atgc_cbl2    = powheginput('#cbl2')
      if(atgc_cbl2<0d0) atgc_cbl2 = 1.5e-5
      atgc_cpwwwl2 = powheginput('#cpwwwl2')
      if(atgc_cpwwwl2<0d0) atgc_cpwwwl2 = 3e-6
      atgc_cpwl2   = powheginput('#cpwl2')
      if(atgc_cpwl2<0d0) atgc_cpwl2 = 1e-6
      atgc_cbwl4   = powheginput('#cbwl4')
      if(atgc_cbwl4<0d0) atgc_cbwl4 = 2e-12
      atgc_cwwl4   = powheginput('#cwwl4')
      if(atgc_cwwl4<0d0) atgc_cwwl4 = 3.5e-12
      atgc_cbbl4   = powheginput('#cbbl4')
      if(atgc_cbbl4<0d0) atgc_cbbl4 = 2e-12
      atgc_cbtwl4  = powheginput('#cbtwl4')
      if(atgc_cbtwl4<0d0) atgc_cbtwl4 = 2e-12

      end subroutine recola_atgc

      
      subroutine recola_generate_process()
      implicit none
      include 'nlegborn.h'
      include 'PhysPars.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_res.h'
      integer i,k,j
      character*100 proc
c     is_fs is filled by getisfsparticles both for the Born and for the
c     real flavour structures, so it must be dimensioned to nlegreal,
c     not nlegborn (the original overflowed by one element on the reals;
c     harmless without bounds checking, caught here by -fbounds-check).
      integer is_fs(nlegreal),isfslength
      integer bflav(nlegborn,maxprocborn),rflav(nlegreal,maxprocreal)
      logical skip
      integer bflav_ordered1(nlegbornexternal),bflav_ordered2(nlegbornexternal)
      integer rflav_ordered1(nlegrealexternal),rflav_ordered2(nlegrealexternal)
      integer bflav_clean(nlegbornexternal,maxprocborn),rflav_clean(nlegrealexternal,maxprocreal)
      integer ubflav(nlegbornexternal-1,maxprocborn),ubflav_clean(nlegbornexternal-1,maxprocborn)
      integer nprborn_loop
      integer ubpr, bpr, rpr, shift, npruborn, nprborn, nprreal
      common/nprocesses_rcl/npruborn, nprborn, nprreal
      integer qed_qcd_rcl
      real * 8 powheginput
      external powheginput
      common/pertord/qed_qcd_rcl
      integer uub_st, uub_ew
     
 
      qed_qcd_rcl=int(powheginput("#qed_qcd"))
      ! CB+CS: qed_qcd is flag for the choice of the correction
      ! If you choose qed_qcd=0, you are performing a NLO QCD calculation:
      ! the real will have an additional alphas.
      ! If you choose qed_qcd=0, you are performing a NLO QED calculation:
      ! the real will have an additional alphaew.
      

      
      ! This subroutine should be called after build_resonance_histories
      ! in init_processes, since the subroutine building the resonance 
      ! histories might also discard some input flavour configuration
      ! depending on the model used (i.e a diagonal ckm matrix requirement
      ! will discard all inconsistent flavour string inputs). Therefore,
      ! one has to make sure to:
      ! 1) REMOVE RESONANCES from flavour string
      do k=1,flst_nborn
         call getisfsparticles(flst_bornlength(k),flst_born(:,k),flst_bornres(:,k),
     1     isfslength,is_fs)
         do j=1,isfslength
            bflav(j,k)=flst_born(is_fs(j),k)
         enddo
      enddo

      if(.not.flg_bornonly)then
         do k=1,flst_nreal
            call getisfsparticles(flst_reallength(k),flst_real(:,k),flst_realres(:,k),
     1           isfslength,is_fs)
            do j=1,isfslength
               rflav(j,k)=flst_real(is_fs(j),k)
            enddo
         enddo
      endif

      ! 2) REMOVE EQUAL FLAVOUR CONFIGURATION, which differed only for the 
      !    resonance structure (this is not fully needed but will avoid 
      !    collier to store processes which are anyway releted via crossing
      !    and never used). Moreover, the flavours of final state partons should
      !    be ordered with a given convention (to account for final state 
      !    particle swapping). The subroutine flav_order orders them 
      !    in decreasing order of flav
      nprborn=flst_nborn
      bpr=0
      do k=1,flst_nborn
         skip=.false.
         call flav_order(bflav(:nlegbornexternal,k),bflav_ordered1)
         do j=k+1,flst_nborn
            call flav_order(bflav(:nlegbornexternal,j),bflav_ordered2)
            if(all(bflav_ordered1==bflav_ordered2))then
               nprborn=nprborn-1
               skip=.true.
               goto 9
            endif
         enddo
 9       continue
         if(.not.skip)then
            bpr=bpr+1
            bflav_clean(:nlegbornexternal,bpr)=bflav_ordered1(:)
            if(qed_qcd_rcl.eq.2)then
               bflav_clean(:nlegbornexternal,bpr+flst_nborn)=bflav_ordered1(:)
            endif
         endif
      enddo

      if(.not.flg_bornonly)then
         nprreal=flst_nreal
         rpr=0 
         do k=1,flst_nreal
            skip=.false.
            call flav_order(rflav(:nlegrealexternal,k),rflav_ordered1)
            do j=k+1,flst_nreal
               call flav_order(rflav(:nlegrealexternal,j),rflav_ordered2)
               if(all(rflav_ordered1==rflav_ordered2))then
                  nprreal=nprreal-1
                  skip=.true.
                  goto 99
               endif
            enddo
 99         continue
            if(.not.skip)then
               rpr=rpr+1
               rflav_clean(:nlegrealexternal,rpr)=rflav_ordered1(:)
            endif
         enddo
      else
         nprreal = 0
      endif

      ! Recola requires the processes to be used
      ! to be defined at the very beginning. So,
      ! if MiNLO/MiNNLO is needed, the uuborn
      ! processes need to be allocated now
      npruborn=0
      if(flg_minlo.or.flg_minnlo)then
         do k=1,nprborn
            call find_uub_for_minlo(bflav_clean(:nlegbornexternal,k),ubflav(:,k))
         enddo
         ! Remove same uub configuration
         npruborn=nprborn 
         ubpr=0
         do k=1,nprborn
            skip=.false.
            do j=k+1,flst_nborn               
               if(all(ubflav(:,k)==ubflav(:,j)))then
                  npruborn=npruborn-1
                  skip=.true.
                  goto 999
               endif
            enddo
 999        continue
            if(.not.skip)then
               ubpr=ubpr+1
               ubflav_clean(:,ubpr)=ubflav(:,k)
            endif
         enddo         
      endif

      if(qed_qcd_rcl.eq.2) then
         allocate(processes(npruborn+nprborn*2+nprreal))
      else
         allocate(processes(npruborn+nprborn+nprreal))
      endif

      if(flg_minlo.or.flg_minnlo)then
         if(res_powst==0)then
c     purely EW Born: the underlying Born removes a photon
            uub_st = 0
            uub_ew = 1
         else
c     QCD Born (res_powst>=1): the underlying Born removes one parton.
c     The original interface only handled res_powst<=1 and left these
c     uninitialised otherwise; bbH_yt2 has res_powst=4.
            uub_st = 1
            uub_ew = 0
         endif
         do k=1,npruborn
            call flav_to_string(ubflav_clean(:nlegbornexternal-1,k),proc)
            processes(k)%i_mast=k
            allocate(processes(k)%flav(nlegbornexternal-1))
            processes(k)%flav=ubflav_clean(:nlegbornexternal-1,k)
            call define_process_rcl(k,proc(:),'NLO')
            call unselect_all_powers_BornAmpl_rcl(k)
            call select_power_BornAmpl_rcl(k, 'QCD', res_powst-uub_st)
            call select_power_BornAmpl_rcl(k, 'QED', res_powew-uub_ew)
            call unselect_all_powers_LoopAmpl_rcl(k)
            call select_power_LoopAmpl_rcl(k, 'QCD', res_powst+2-uub_st)
            call select_power_LoopAmpl_rcl(k, 'QED', res_powew-uub_ew)
            if(model_rcl .eq. 'SM (QCD) + ATGC')then
               ! See arXiv:1804.01477
               ! i = 0 SM^2
               ! i = 1 SMxEFT6 + EFT6^2
               ! i = 2 SMxEFT8
               do i=0,atgc_powlamm2
                  call select_power_BornAmpl_rcl(k, 'LAM', i)
                  call select_power_LoopAmpl_rcl(k, 'LAM', i)
               enddo
            endif
!     call split_collier_cache_rcl(k,10) ! This is added to reduce the memory consuption.       
         enddo         
      endif

      shift = npruborn
      if(qed_qcd_rcl.eq.2)then
         nprborn_loop = nprborn*2
      else
         nprborn_loop = nprborn
      endif

      do k=1,nprborn_loop


            call flav_to_string(bflav_clean(:nlegbornexternal,k),proc)
            print*, '**index=', k+shift, ':', proc
         
         !call flav_to_string(bflav_clean(:nlegbornexternal,k),proc)
         processes(k+shift)%i_mast=k+shift
         allocate(processes(k+shift)%flav(nlegbornexternal))
         processes(k+shift)%flav=bflav_clean(:nlegbornexternal,k)
         call define_process_rcl(k+shift,proc(:),'NLO')
         call unselect_all_powers_BornAmpl_rcl(k+shift)
         call unselect_all_powers_LoopAmpl_rcl(k+shift)
         call select_power_BornAmpl_rcl(k+shift, 'QCD', res_powst)
         call select_power_BornAmpl_rcl(k+shift, 'QED', res_powew)
         if(qed_qcd_rcl.eq.0) then
            call select_power_LoopAmpl_rcl(k+shift, 'QCD', res_powst+2)
            call select_power_LoopAmpl_rcl(k+shift, 'QED', res_powew)
         elseif(qed_qcd_rcl.eq.1) then
            call select_power_LoopAmpl_rcl(k+shift, 'QCD', res_powst)
            call select_power_LoopAmpl_rcl(k+shift, 'QED', res_powew+2)
         elseif(qed_qcd_rcl.eq.2) then
            if(k.le.nprborn)then
               call select_power_LoopAmpl_rcl(k+shift, 'QCD', res_powst+2)
               call select_power_LoopAmpl_rcl(k+shift, 'QED', res_powew)
            else
               call select_power_LoopAmpl_rcl(k+shift, 'QCD', res_powst)
               call select_power_LoopAmpl_rcl(k+shift, 'QED', res_powew+2)
            endif
         endif
         if(model_rcl .eq. 'SM (QCD) + ATGC')then
            ! See arXiv:1804.01477
            ! i = 0 SM^2
            ! i = 1 SMxEFT6 + EFT6^2
            ! i = 2 SMxEFT8
            do i=0,atgc_powlamm2
               call select_power_BornAmpl_rcl(k+shift, 'LAM', i)
               call select_power_LoopAmpl_rcl(k+shift, 'LAM', i)
            enddo
         endif
!         call split_collier_cache_rcl(k,10) ! This is added to reduce the memory consuption.       
      enddo

      if(qed_qcd_rcl.eq.2)then
         shift=shift+nprborn*2
      else
         shift=shift+nprborn
      endif

      do k=1,nprreal

         call flav_to_string(rflav_clean(:nlegrealexternal,k),proc)
         print*, '**index=', k+shift, ':', proc
          
         
 !        call flav_to_string(rflav_clean(:nlegrealexternal,k),proc)
         processes(k+shift)%i_mast=k+shift
         allocate(processes(k+shift)%flav(nlegrealexternal))
         processes(k+shift)%flav=rflav_clean(:nlegrealexternal,k)
         call define_process_rcl(k+shift,proc,'LO')
         call unselect_all_powers_BornAmpl_rcl(k+shift)
         if(qed_qcd_rcl.eq.0) then         
            call select_power_BornAmpl_rcl(k+shift, 'QCD', res_powst+1)
            call select_power_BornAmpl_rcl(k+shift, 'QED', res_powew)
         elseif(qed_qcd_rcl.eq.1) then
            call select_power_BornAmpl_rcl(k+shift, 'QCD', res_powst)
            call select_power_BornAmpl_rcl(k+shift, 'QED', res_powew+1)
         elseif(qed_qcd_rcl.eq.2) then
!            if(processes(k+shift)%flav(nlegrealexternal).ne.22) then
            if(any(processes(k+shift)%flav.eq.22)) then               
               call select_power_BornAmpl_rcl(k+shift, 'QCD', res_powst)
               call select_power_BornAmpl_rcl(k+shift, 'QED', res_powew+1)
            else
               call select_power_BornAmpl_rcl(k+shift, 'QCD', res_powst+1)
               call select_power_BornAmpl_rcl(k+shift, 'QED', res_powew)
            endif

         endif
         if(model_rcl .eq. 'SM (QCD) + ATGC')then
            ! See arXiv:1804.01477
            ! i = 0 SM^2
            ! i = 1 SMxEFT6 + EFT6^2
            ! i = 2 SMxEFT8
            do i=0,atgc_powlamm2
               call select_power_BornAmpl_rcl(k+shift, 'LAM', i)              
            enddo
         endif
      enddo


      if(.not.als_rescaling) call recola_freeze_alphas

      call generate_processes_rcl

      ! Initialize helicity and colorflow configurations
      call complete_process_registration

      contains

      subroutine find_uub_for_minlo(bflav,ubflav) 
      ! MiNLO is supposed to be applied to processes
      ! having a Born with 1jet only
      implicit none
      integer bflav(:),ubflav(:)
      integer i

      do i=1,size(bflav)
         if(.not.is_parton(bflav(i)))then
            ubflav(i)=bflav(i)
         endif
      enddo
      
      if(bflav(1).ne.0 .and. bflav(2).ne.0)then
         ubflav(1:2)=bflav(1:2)
      else if(bflav(1).eq.0 .and. bflav(2).ne.0)then
         ubflav(1)=-bflav(size(bflav))
         ubflav(2)=bflav(2)
      else if(bflav(1).ne.0 .and. bflav(2).eq.0)then
         ubflav(1)=bflav(1)
         ubflav(2)=-bflav(size(bflav))
      else if(bflav(1).eq.0 .and. bflav(2).eq.0)then
         ubflav(1:2)=bflav(1:2)
      endif

      end subroutine find_uub_for_minlo

      subroutine complete_process_registration
      use globals_rcl, only: prs, c0EffMax, cEffMax!,get_pr
      implicit none
      integer i_mast,i_cross
      integer,allocatable::cols(:,:),cols_cross(:,:)
      integer,allocatable::hels(:,:),hels_cross(:,:)
      integer,allocatable::cols_Nc(:,:)
      integer ncols,nhels
      integer nlength
      integer :: i,j,k,l
 
      do i=1,size(processes)
        i_mast=processes(i)%i_mast
        i_cross=i_mast
!        call get_pr(i_mast,"recola_colour",i_mast_in)
        if(prs(i_mast)%crosspr.gt.0)then
           ! Use prs(i_mast)%crosspr as i_mast to extract colour 
           ! and helicity configurations, which are empty for
           ! processes related by symmetry. Then, put i_mast to
           ! the original value for the amplitude computation
           i_mast=prs(i_mast)%crosspr
        endif

        ! Check if process exists in Recola
        if(prs(i_mast)%loop.and.c0EffMax(i_mast).eq.0.and.cEffMax(i_mast).eq.0)then ! NLO process
           write(*,*) " The process ",processes(i)%flav
           write(*,*) " is not available in Recola neither"
           write(*,*) " at LO nor at NLO !!!"
           stop
        else if(c0EffMax(i_mast).eq.0)then ! LO process 
           write(*,*) " The process ",processes(i)%flav
           write(*,*) " is not available in Recola at LO"
           stop
        endif

        nlength=size(processes(i)%flav)

        if(allocated(cols))deallocate(cols)
        call get_colour_configurations_rcl(i_mast,cols)
        ncols=size(cols,2)

c        print*, 'before large NC'
c        print*, 'processes(i)%flav= ', processes(i)%flav
c        print*, 'processes(i)%cols= ', cols
        

        ! Processes related by crossing turn out to have
        ! empty cols(:,:) ... fill them!
        if(prs(i_cross)%crosspr.gt.0)then
           if(allocated(cols_cross))deallocate(cols_cross)
           allocate(cols_cross(size(cols,1),size(cols,2)))
           do l=1,ncols
              ! Ex. 060001->600001 if 1<->2
              do j=1,nlength!size(flasv)
                 cols_cross(j,l)=cols(prs(i_cross)%relperm(j),l)
              enddo
              ! Ex. 600001->600002 to get correct result
              do j=1,nlength
                 if(prs(i_cross)%relperm(j).ne.j)then
                    do k=1,nlength
                       if(cols(prs(i_cross)%relperm(k),l)==prs(i_cross)%relperm(j))then
                          cols_cross(k,l)=j
                       endif
                    enddo
                 endif
              enddo
           enddo
           cols=cols_cross
        endif
        
        ! Recola works in color flow basis, but cols also contains the
        ! colour suppressed configurations. Parton showers work in the 
        ! large N_c limit, so these configurations must be removed

c        if(processes(i)%flav(1).eq.0) then       
           call large_Nc_limit(processes(i)%flav,cols,cols_Nc)
           ncols=size(cols_Nc,2)        
           allocate(processes(i)%cols(size(cols_Nc,1),size(cols_Nc,2)))
           processes(i)%cols=cols_Nc
c        else
c           cols_Nc=cols
c           processes(i)%cols=cols_Nc
c        endif

c        print*, '****CB+CS debug****'
c        print*, 'processes(i)%flav= ', processes(i)%flav
c        print*, 'processes(i)%cols= ', processes(i)%cols

        if(allocated(hels))deallocate(hels)
        call get_helicity_configurations_rcl(i_mast,hels)
        nhels=size(hels,2)
        ! Processes related by crossing turn out to have
        ! empty hels(:,:) ... fill them!
        if(prs(i_cross)%crosspr.gt.0)then
           if(allocated(hels_cross))deallocate(hels_cross)
           allocate(hels_cross(size(hels,1),size(hels,2)))
           do j=1,nhels
              do k=1,nlength
                 hels_cross(k,j)=hels(prs(i_cross)%relperm(k),j)
              enddo
           enddo
           hels=hels_cross
        endif
        
        allocate(processes(i)%hels(size(hels,1),size(hels,2)))
        processes(i)%hels=hels
        
      enddo
      
      end subroutine complete_process_registration

      subroutine large_Nc_limit(flav,cols,cols_Nc)
      implicit none
      integer ::cols(:,:),i,j,l,m,count,qcount
      integer ::flav(:),qflav(100)
      integer,allocatable,intent(out)::cols_Nc(:,:)
      integer,allocatable:: mark(:),index_sh(:)
      logical qqbpair
      logical chain,chain_ini
      integer chain_length,colour_length,curr,chain_head

      logical myflag
      
      ! For colorflow formalism see hep-ph/0209271
      if(allocated(mark))deallocate(mark)
      allocate(mark(size(cols,2)))
      mark=1

      !!!!!!!!!!!!!!!!!!!!!
      ! (1) FIRST CHECK   !
      !!!!!!!!!!!!!!!!!!!!!

c      print*, 'first check'
      
      outer: do i=1,size(cols,2)
         inner: do j=1,size(cols,1)
            if(cols(j,i)==j)then
               ! This just suffices for diagrams having only external
               ! gluons or like like qqb + m gluons
               ! with one qqb pair, where one just spots subleading
               ! diagrams by looking for external U(1) gluons
               mark(i)=0
               exit inner
            endif
         end do inner
      end do outer

c      print*, 'mark= ', mark(:)

      ! Now check how many pairs of qqb we have,
      ! starting by counting quarks
      qcount=0
      qflav=0
      do i=1,size(flav)
         if(abs(flav(i)).le.6.and.flav(i).ne.0)then
            qcount=qcount+1
            qflav(qcount)=flav(i)
         endif
      enddo

c      print*, 'qcount= ', qcount
c      print*, 'qflav(qcount)= ', qflav(qcount) 

      ! If we have only external gluons or just one qqb pair,
      ! nothing more must be done
      if(qcount.le.2)goto 111

      !!!!!!!!!!!!!!!!!!!!!!
      ! (2) SECOND CHECK   !
      !!!!!!!!!!!!!!!!!!!!!!

      ! If we have at least two qqb pairs of indistinguishable
      ! quarks, then we keep all colorflows, which may correspond
      ! both to leanding and subleading configurations (qqb pairs also
      ! include s-channel diagrams -> this is not true for qq pairs)
      qqbpair=.false.
      do i=1,qcount
         do j=i+1,qcount
            if(qflav(i)==-qflav(j))then
               ! One qqb pair found ... now look for 
               ! a second pair
               do l=1,qcount
                  if(l.eq.i.or.l.eq.j)cycle
                  do m=1,qcount
                     if(m.eq.i.or.m.eq.j.or.m.eq.l)cycle
                     if(qflav(l)==-qflav(m))then
                        ! Found second pair. Is it distinguishable
                        ! from the first one?
                        if(abs(qflav(i)).eq.abs(qflav(l)))qqbpair=.true.
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo

c      print*, 'qqbpair= ', qqbpair
c      print*, 'size(cols,2)= ', size(cols,2)
c      print*, 'size(cols,1)= ', size(cols,1)
      
      ! If all qqb pairs (if any) are distinguishable, some extra subleading
      ! configurations must be removed
      if(.not.qqbpair)then

      !CB+CS
         do i=1,size(cols,2)
            myflag=.false.
            do j=1,size(cols,1)
               if(cols(j,i).ne.0) then
                  curr=cols(j,i)
                  if(abs(flav(curr)) .ne. abs(flav(j))) then
                     myflag=.true.
                  endif
               endif
            enddo
            if(.not. myflag) then
               mark(i)=0
            endif
         enddo

      endif
c$$$         
c$$$         ! We need to check for broken colour chains, representing 
c$$$         ! subleading colour configurations. For a chain to be broken,
c$$$         ! we need to look for separate color flows, each beginning 
c$$$         ! with an incoming antiquark (outgoing quark) and ending with an
c$$$         ! outgoing antiqu ark or incoming quark of the same flavour.
c$$$         do i=1,size(cols,2)
c$$$            colour_length=0
c$$$            chain_ini=.true.
c$$$            ! Count number of non-zero cols entry
c$$$            do j=1,size(cols,1)
c$$$               print*, 'cols(', j,',',i,')= ', cols(j,i)
c$$$               if(cols(j,i).ne.0)then
c$$$                  colour_length=colour_length+1
c$$$                  print*, 'colour_length= ', colour_length
c$$$                  if(chain_ini.and.flav(j).ne.0)then
c$$$                     chain_ini=.false.
c$$$                     chain_head=j
c$$$                     print*, 'chain_ini now is', chain_ini
c$$$                     print*, 'chain_head= ', chain_head
c$$$                  endif
c$$$               endif
c$$$            end do
c$$$            ! Now let's follow a colour chain. If its number of
c$$$            ! rings is less than colour_length, we have at least 
c$$$            ! two colour chains and the colour configuration is 
c$$$            ! at least 1/N suppressed
c$$$            curr=cols(chain_head,i)
c$$$            print*, 'curr= ', curr
c$$$            chain_length=1
c$$$            chain=.true.
c$$$            do while(chain)
c$$$               if(flav(curr).eq.0)then
c$$$                  curr=cols(curr,i)
c$$$                  chain_length=chain_length+1
c$$$               else
c$$$                  chain=.false.
c$$$               endif
c$$$            enddo
c$$$            print*, 'j= ', j
c$$$            print*, 'chainlenght= ', chain_length
c$$$            print*, 'colourlenght= ', colour_length
c$$$            if(chain_length.ne.colour_length) mark(i)=0
c$$$         end do 
c$$$      endif

c      print*, 'mark after check2= ', mark(:)
      
 111  continue

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! (3) NOW REMOVE SUBLEADING CONFIGURATIONS  !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !print*, 'I am here'
      !print*, 'size(mark)=', size(mark)
      !print*, 'size(cols,1)= ', size(cols,1)

      if(allocated(index_sh))deallocate(index_sh)
      !allocate(index_sh(size(cols,1)))

      allocate(index_sh(size(cols,1)+1))

      index_sh=0
      count=0
      do i=1,size(mark)
         
         !print*, 'mark(', i, ')', mark(i)

         if(mark(i)==1)then
            count=count+1
            index_sh(i)=count

            !print*, 'index_sh(', i, ')', index_sh(i)

         endif
      enddo
      if(allocated(cols_Nc))deallocate(cols_Nc)
      allocate(cols_Nc(size(cols,1),count))
      do i=1,size(cols,2)
         if(mark(i)==1)then
            do j=1,size(cols,1)
               cols_Nc(j,index_sh(i))=cols(j,i)
            end do
         endif
      end do        
      end subroutine large_Nc_limit


      end subroutine recola_generate_process

      logical function is_parton(i)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_res.h'
      integer i
      is_parton=.false.

      if(res_powst==0)then
! Case for purely EW born: underlying Born 
! obtained removing a photon
         if(i.eq.22)is_parton=.true.
      endif

      if(i.eq.0.or.abs(i).lt.6)is_parton=.true.
      end

      subroutine flav_to_string(flav_ordered,proc)
      implicit none
      integer, intent(in):: flav_ordered(:)
      character,intent(inout):: proc(100)
      integer i,fl,l,ip

      ! The flavours of final state partons should
      ! be ordered with a given convention (to account
      ! for final state particle swapping). They are
      ! passed already ordered in decreasing order of flav

      proc=''
      l=0
      do i=1,sizeof(flav_ordered)/4
         l=l+1
         fl=abs(flav_ordered(i))
         SELECT CASE (fl)
         CASE (1)
            proc(l)='d'
         CASE (2)
            proc(l)='u'
         CASE (3)
            proc(l)='s'
         CASE (4)
            proc(l)='c'
         CASE (5)
            proc(l)='b'
         CASE (6)
            proc(l)='t'
         CASE (22)
            proc(l)='A'
         CASE (23)
            proc(l)='Z'
         CASE (24)
            proc(l)='W'
         CASE (25)
c     Higgs boson: needed by bbH_yt2, absent from the bb4l original
            proc(l)='H'
         CASE (0)
            proc(l)='g'            
         CASE (11)
            proc(l)='e'            
         CASE (12)
            proc(l:l+3)=(/ 'n','u','_','e'/)
            l=l+3
         CASE (13)
            proc(l:l+1)=(/ 'm','u'/)
            l=l+1
         CASE (14)
            proc(l:l+4)=(/ 'n','u','_','m','u'/)
            l=l+4
         CASE (15)
            proc(l:l+2)=(/ 't','a','u'/)
            l=l+2
         CASE (16)
            proc(l:l+5)=(/ 'n','u','_','t','a','u'/)
            l=l+5
         CASE DEFAULT
            WRITE(*,*)  "Hmmmm, I don't know ",fl
            stop
         END SELECT
         if( (abs(flav_ordered(i))<=6 .or. abs(flav_ordered(i))==12 
     &        .or. abs(flav_ordered(i))== 14 .or. abs(flav_ordered(i))==16 )
     &        .and. flav_ordered(i)<0) then
            l=l+1
            proc(l)='~'
         else if(abs(flav_ordered(i))==11.or. abs(flav_ordered(i))==13 
     &           .or.abs(flav_ordered(i))==15 ) then
            l=l+1
            if(flav_ordered(i).gt.0) then
               proc(l)='-'
            else
               proc(l)='+'
            endif
         endif
         if(abs(flav_ordered(i))==24) then
            l=l+1
            if(flav_ordered(i).lt.0) then
               proc(l)='-'
            else
               proc(l)='+'
            endif
         endif
         l=l+1
         proc(l)=' '
         if(i==2)then 
            proc(l+1)='-'
            proc(l+2)='>'
            proc(l+3)=' '
            l=l+3
         endif
      enddo
      
      end subroutine flav_to_string


      
      subroutine flav_order(flav,flav_ordered)
      implicit none
      integer, intent(in):: flav(:)
      integer  flav_ordered(:)
      integer i, ip

      flav_ordered(1:sizeof(flav)/4)=flav(1:sizeof(flav)/4)
      ! get first final state parton position
      do i=3,sizeof(flav)/4
         if(is_parton(flav(i)))then
            ip=i
            exit               
         endif
      enddo

      
c      ip=5
            
c      call sort(flav_ordered(ip:sizeof(flav)/4))

      contains

      subroutine sort(array) 
      implicit none
      integer,intent(inout)::array(:)
      integer temp,i,j
      temp=0
      do i = 1,sizeof(array)/4
         do j=1,sizeof(array)/4
            if (array(j).lt.array(i)) then
               temp = array(j)
               array(j) = array(i)
               array(i) = temp
            endif
         enddo
      enddo
      end

      end subroutine flav_order

      integer function get_rcl_index(flav) result(index)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer, intent(in):: flav(:)
      integer, allocatable :: flav_ordered(:)
      integer k, k0, kn
      integer npruborn, nprborn, nprreal
      common/nprocesses_rcl/npruborn, nprborn, nprreal
      logical found
      integer qed_qcd_rcl
      common/pertord/qed_qcd_rcl


      ! If the processes is not register, it might be a problem 
      ! of ordering of the final state partons ... order the input
      ! string according to the convention used in the registration 
      ! of the process ... 
      if(allocated(flav_ordered))deallocate(flav_ordered)
      allocate(flav_ordered(size(flav)))
      call flav_order(flav,flav_ordered)
      
      if(size(flav)==nlegbornexternal-1)then
         k0=1
         kn=npruborn
      else if(size(flav)==nlegbornexternal)then
         ! For qed_qcd_rcl = 2, only the first set of i_mast (corresponding to
         ! the ones including QCD corrections) is searched. That does not
         ! make a difference for the Born. For the virtual, that is accounted
         ! for in the proper subroutine.
         k0=npruborn+1
         kn=npruborn+nprborn         
      else if (size(flav)==nlegrealexternal)then
         if(qed_qcd_rcl.eq.2)then
            k0=npruborn+2*nprborn+1
            kn=npruborn+2*nprborn+nprreal
         else
            k0=npruborn+nprborn+1
            kn=npruborn+nprborn+nprreal
         endif
      endif

      found=.false.
      do k=k0, kn
         if(all(processes(k)%flav == flav_ordered))then
            index=processes(k)%i_mast
            found=.true.
            exit
         endif
      enddo

      if(.not.found)then
         write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
         write(*,*) " Error in get_rcl_index: the process  "
         write(*,*) " ",flav(:)
         write(*,*) " was not found!!!"
         write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
         stop
      endif

      end function get_rcl_index
      
      
      subroutine recola_born(p,bflav,born,bornjk,bmunu)
      include 'pwhg_st.h'
      include 'PhysPars.h'
      include 'pwhg_math.h'      
      include 'nlegborn.h'
      include 'pwhg_res.h'
      integer, parameter :: nlegs=nlegbornexternal      
      double precision,intent(in) :: p(:,:)
      double precision  :: p_rcl(0:3,nlegs)
      integer,intent(in) :: bflav(:)
      integer :: bflav_rcl(nlegs)
      double precision,intent(out) :: born
      double precision born_tmp
      double precision,intent(out),optional :: bornjk(nlegs,nlegs)
      double precision bornjk_tmp(nlegs,nlegs)
      double precision,intent(out),optional :: bmunu(0:3,0:3,nlegs)
      double precision bmunu_tmp(0:3,0:3,nlegs)
      integer :: i_mast,mu,nu,i,j,k,l
      integer uub_st, uub_ew !order
      double precision asf
      complex *16 v(0:3)


      
      i_mast=get_rcl_index(bflav)
      
      if(size(bflav)==nlegs)then
         uub_st = 0
         uub_ew = 0
      else if(size(bflav)==nlegs-1)then
         ! Get amplitude for underlying born
         if(res_powst==0)then
            ! purely EW Born: underlying Born removes a photon
            uub_st = 0
            uub_ew = 1
         else
            ! QCD Born (res_powst>=1): underlying Born removes a parton
            uub_st = 1
            uub_ew = 0
         endif
      else
         write(*,*) "unknow order in recola_born ... "
         stop
      endif
      if(als_rescaling) then
         call set_alphas_rcl(st_alpha,dsqrt(st_muren2),st_nlight)
      endif
      
      call remap_momenta(i_mast, bflav, p, p_rcl)

      
      call compute_process_rcl(i_mast,p_rcl(0:3, 1 : (nlegs-1 + 1 - uub_st - uub_ew) ),'LO')
      
      born=0d0
      if(model_rcl .eq. 'SM (QCD) + ATGC')then
         do i=0,atgc_powlamm2
            call get_squared_amplitude_rcl(i_mast,[2*res_powst,2*res_powew,2*i],'LO',born_tmp)
            born=born+born_tmp
         enddo
      else
         call get_squared_amplitude_rcl(i_mast,[2*(res_powst-uub_st),2*(res_powew-uub_ew)],'LO',born)
      endif

      if (present(bornjk) ) then ! color-correlated
!       call compute_all_colour_correlations_rcl(i_mast,p,'LO')
         bornjk=0d0
         do j=1,size(bflav)
            if (abs(bflav(j)).gt.6)  cycle
            do k=j,size(bflav)
               if (abs(bflav(k)).gt.6)  cycle
               if(als_rescaling) then
                  call rescale_colour_correlation_rcl(i_mast,j,k,'LO')
               else
c     Without qcd_rescaling the fast rescaling variants are refused by
c     Recola; the correlation has to be recomputed from the momenta.
                  call compute_colour_correlation_rcl(i_mast,
     $                 p_rcl(0:3,1:(nlegs-1+1-uub_st-uub_ew)),j,k,'LO')
               endif
               if(model_rcl .eq. 'SM (QCD) + ATGC')then
                  do i=0,atgc_powlamm2
                     call get_colour_correlation_rcl(i_mast,[2*res_powst,2*res_powew,2*i],j,k,'LO',bornjk_tmp(j,k))
                     bornjk(j,k)=bornjk(j,k)+bornjk_tmp(j,k)
                  enddo
               else
                  if(res_powst==0)then
                     bornjk(j,k)=0d0
                  else
                     call compute_colour_correlation_rcl(i_mast,p_rcl,j,k)
                     call get_colour_correlation_rcl(i_mast,[2*res_powst,2*res_powew],j,k,'LO',bornjk(j,k))
                  endif
               endif
               if (bflav(j) == 0) then
                  ! extra minus and color factor
                  ! related to recola conventions
                  bornjk(j,k) = -bornjk(j,k)*ca
               else
                  bornjk(j,k) = -bornjk(j,k)*cf
               endif
               bornjk(k,j) = bornjk(j,k)
            enddo
         enddo
      endif

      if (present(bmunu)) then  ! spin-correlated
       bmunu=0d0
       v = (/ 0d0,0d0,0d0,0d0 / )
       do j=1,size(bflav)
          if (bflav(j) == 0) then
             do i = 0,3
                do k = 0,3
                   if(k.eq.i)then
                      v(i) = (1d0,0d0)
                   else
                      v(i) = (1d0,0d0)
                      v(k) = (1d0,0d0)
                   endif
                   if(als_rescaling) then
                      call rescale_spin_correlation_rcl(i_mast,j,v,'LO')
                   else
                      call compute_spin_correlation_rcl(i_mast,
     $                     p_rcl(0:3,1:(nlegs-1+1-uub_st-uub_ew)),j,v,'LO')
                   endif
                   if(model_rcl .eq. 'SM (QCD) + ATGC')then
                      do l=0,atgc_powlamm2
                         call get_spin_correlation_rcl(i_mast,[2*res_powst,2*res_powew,2*l],'LO',bmunu_tmp(i,k,j))
                         bmunu(i,k,j)=bmunu(i,k,j)+bmunu_tmp(i,k,j)
                      enddo
                   else
                      call get_spin_correlation_rcl(i_mast,[2*res_powst,2*res_powew],'LO',bmunu(i,k,j))
                   endif
                   v(:)=0d0
                enddo
             enddo
             do i = 0,3
                do k = 0,3
                   if(k.ne.i)then
                      bmunu(i,k,j)=(bmunu(i,k,j)-bmunu(i,i,j)-bmunu(k,k,j))/2d0
                   endif
                enddo
             enddo
          endif
       enddo
      endif

c     exact running-alpha_s restoration (see recola_asfact)
      asf = recola_asfact(res_powst-uub_st)
      born = born*asf
      if(present(bornjk)) bornjk = bornjk*asf
      if(present(bmunu))  bmunu  = bmunu*asf

      end subroutine recola_born


      
!      subroutine recola_virtual(p,vflav,virtual,Vqcd,Vqed)
      subroutine recola_virtual(p,vflav,virtual)
      include 'pwhg_st.h'
      include 'PhysPars.h'
      include 'nlegborn.h'
      include 'pwhg_res.h'
      include 'pwhg_em.h'
      integer, parameter :: nlegs=nlegbornexternal      
      double precision,intent(in) :: p(:,:)
      double precision p_rcl(0:3,1:nlegs)
      integer,intent(in) :: vflav(:)
      double precision,intent(out) :: virtual
!      double precision,intent(out) :: Vqcd,Vqed
      double precision :: virtual_qcd,virtual_qed
      double precision virtual_tmp
      double precision b0wc
      real * 8 powheginput
      external powheginput
      integer :: i_mast,i,j,k
      integer uub_st, uub_ew
      integer, parameter :: dp = kind (23d0) ! double precision
      real(dp), parameter :: pi = 3.141592653589793238462643d0
      integer qed_qcd_rcl
      common/pertord/qed_qcd_rcl
      integer npruborn, nprborn, nprreal
      common/nprocesses_rcl/npruborn, nprborn, nprreal

      i_mast=get_rcl_index(vflav)
      if(size(vflav)==nlegs)then
         uub_st = 0
         uub_ew = 0
      else if(size(vflav)==nlegs-1)then
         ! Get amplitude for underlying born
         if(res_powst==0)then
            ! purely EW Born: underlying Born removes a photon
            uub_st = 0
            uub_ew = 1
         else if (res_powst.ge.1)then
            ! QCD Born (res_powst>=1): underlying Born removes a parton
            uub_st = 1
            uub_ew = 0
         else
            write(*,*) " Unknow coupling combination for underlying "
            write(*,*) " born in recola_born ... "
            stop
         endif
      else
         write(*,*) "unknow order in recola_virtual ... "
         stop
      endif

      if(als_rescaling) then
         call set_alphas_rcl(st_alpha,dsqrt(st_muren2),st_nlight)
      endif  
      call set_mu_ir_rcl (dsqrt(st_muren2))
      call set_mu_uv_rcl (dsqrt(st_muren2))
c     muMS is Recola's Qren (the scale alpha_s is defined at). It
c     enters the alpha_s counterterm as beta0*log(muUV^2/muMS^2) and
c     must track muR together with muUV/muIR.
      call set_mu_ms_rcl (dsqrt(st_muren2))




      call remap_momenta(i_mast, vflav, p, p_rcl)

      virtual=0d0
      

      if(model_rcl .eq. 'SM (QCD) + ATGC')then
         call compute_process_rcl(i_mast,p_rcl(:,:nlegs-1 + 1 - uub_st - uub_ew),'NLO')
         do i=0,atgc_powlamm2
            call get_squared_amplitude_rcl(i_mast,[2*(res_powst+1),2*res_powew,2*i],'NLO',virtual_tmp)
            virtual=virtual+virtual_tmp
         enddo
      else
!call get_squared_amplitude_rcl(i_mast,order,'NLO',virtual)
         if(qed_qcd_rcl.eq.0) then
            call compute_process_rcl(i_mast,p_rcl(:,:nlegs-1 + 1 - uub_st - uub_ew),'NLO')
            call get_squared_amplitude_rcl(i_mast,[2*(res_powst+1-uub_st),2*(res_powew-uub_ew)],'NLO',virtual)
         elseif(qed_qcd_rcl.eq.1) then
            call compute_process_rcl(i_mast,p_rcl(:,:nlegs-1 + 1 - uub_st - uub_ew),'NLO')
            call get_squared_amplitude_rcl(i_mast,[2*(res_powst-uub_st),2*(res_powew+1-uub_ew)],'NLO',virtual)

         elseif(qed_qcd_rcl.eq.2) then
            virtual_qcd=0d0
            virtual_qed=0d0
!            if(i_mast.le.nprborn)then
               ! QCD corrections are placed before the EW ones and 
               ! i_mast can only be less or equal to nprborn for the virtual
               call compute_process_rcl(i_mast,p_rcl(:,:nlegs-1+ 1 - uub_st - uub_ew),'NLO')
               call get_squared_amplitude_rcl(i_mast,[2*(res_powst+1-uub_st),2*(res_powew-uub_ew)],'NLO',virtual_qcd )
               call compute_process_rcl(i_mast+nprborn,p_rcl(:,:nlegs-1 + 1 - uub_st - uub_ew),'NLO')
               call get_squared_amplitude_rcl(i_mast+nprborn,[2*(res_powst-uub_st),2*(res_powew+1-uub_ew)],'NLO',virtual_qed )
!            else
!               write(*,*) "Wrong i_mast=",i_mast,"in recola_virtual"
!               write(*,*) "subroutine for qed_qcd_rcl=2!!!"
!               stop
!            endif
!            Vqcd = virtual_qcd/(st_alpha/2d0/pi)
!            Vqed = virtual_qed/(st_alpha/2d0/pi)
            virtual = virtual_qcd + virtual_qed
         endif
         virtual=virtual/(st_alpha/2d0/pi)            
      endif


c     exact running-alpha_s restoration (see recola_asfact)
      virtual = virtual*recola_asfact(res_powst+1-uub_st)

c     Missing piece of the HEFT Wilson coefficient. The model file
c     supplies dcgghfin_QCD2 = 11*aS/(4pi), i.e. C1 = 1+(aS/pi)*11/4,
c     but the matched coefficient (1808.01660 eq. B.12) also has
c     (aS/pi)*(1/6)*log(muR^2/mt^2). A shift dC gives d(sigma) =
c     2*dC*sigma_LO, and POWHEG's V is Born-dimensionful, so
c        dV = 2*dC/(aS/2pi) = (2/3)*log(muR^2/mt^2) * Born,
c     alpha_s independent. It scales like the Born, not the loop.
c     wilsonlog 0 switches it off.
      if(powheginput("#wilsonlog").ne.0d0) then
         call get_squared_amplitude_rcl(i_mast,
     $        [2*(res_powst-uub_st),2*(res_powew-uub_ew)],'LO',b0wc)
         virtual = virtual
     $        + 2d0/3d0*log(st_muren2/ph_tmass**2)*b0wc
     $        *recola_asfact(res_powst-uub_st)
      endif

      end subroutine recola_virtual

      subroutine recola_real(p,rflav,amp2)
      include 'pwhg_st.h'
      include 'PhysPars.h'
      include 'nlegborn.h'
      include 'pwhg_res.h'
      include 'pwhg_em.h'
      integer, parameter :: nlegs=nlegrealexternal      
      double precision,intent(in) :: p(0:3,nlegs)
      double precision p_rcl(0:3,nlegs)
      integer,intent(in) :: rflav(nlegs)
      double precision,intent(out) :: amp2
      double precision amp2_tmp
      integer :: i_mast,i,j,k
      integer, parameter :: dp = kind (23d0) ! double precision
      real(dp), parameter :: pi = 3.141592653589793238462643d0
      integer qed_qcd_rcl
      common/pertord/qed_qcd_rcl      

      i_mast=get_rcl_index(rflav)
      if(als_rescaling) then
         call set_alphas_rcl(st_alpha,dsqrt(st_muren2),st_nlight)
      endif  

      call remap_momenta(i_mast, rflav, p, p_rcl)
      
      call compute_process_rcl(i_mast,p_rcl,'LO')
      amp2=0d0
      if(model_rcl .eq. 'SM (QCD) + ATGC')then
         do i=0,atgc_powlamm2
            call get_squared_amplitude_rcl(i_mast,[4,2*res_powew,2*i],'LO',amp2_tmp)
            amp2=amp2+amp2_tmp
         enddo
      else
         !call get_squared_amplitude_rcl(i_mast,2,'LO',amp2)
         if(qed_qcd_rcl.eq.0) then
            call get_squared_amplitude_rcl(i_mast,[2*res_powst+2,2*res_powew],'LO',amp2)
         elseif(qed_qcd_rcl.eq.1) then
               call get_squared_amplitude_rcl(i_mast,[2*res_powst,2*res_powew+2],'LO',amp2)
         elseif(qed_qcd_rcl.eq.2) then
            if(any(rflav.eq.22)) then
               call get_squared_amplitude_rcl(i_mast,[2,2*res_powew+2],'LO',amp2)
            else
               call get_squared_amplitude_rcl(i_mast,[4,2*res_powew],'LO',amp2)
            endif
         endif
         amp2=amp2/(st_alpha/2d0/pi)            
      endif


c     exact running-alpha_s restoration (see recola_asfact)
      amp2 = amp2*recola_asfact(res_powst+1)

      end subroutine recola_real


! Color flow
      subroutine recola_colour(p,flav,color,ifl)
      use globals_rcl, only: prs!,get_pr
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_res.h'
      integer :: nlegs   
      double precision,intent(in) :: p(:,:)
      double precision p_rcl(0:3,nlegrealexternal)
      integer :: flav(:)
      integer i_mast
      integer ic,ih
      complex*16 amp
      real*8,allocatable,save::ccamp(:)
      double precision, allocatable,save :: m2arr(:)
      integer, intent(out) ::  color(:,:)
      integer, intent(out) :: ifl
      integer :: randomflow, i,j,k
c     ampshift = 0 for the Born (nlegbornexternal legs), 1 for the real
c     emission (nlegrealexternal legs), matching the coupling powers
c     selected when the processes were registered.
      integer :: ampshift

      double precision random, toss, cumm2arr
      external random
      

      randomflow=1

      i_mast=get_rcl_index(flav)
      if(als_rescaling) then
         call set_alphas_rcl(st_alpha,dsqrt(st_muren2),st_nlight)
      endif  
      call remap_momenta(i_mast, flav, p, p_rcl)
      nlegs=size(flav)
      if(nlegs.eq.nlegrealexternal)then
         ampshift = 1
      else
         ampshift = 0
      endif
      call compute_process_rcl(i_mast,p_rcl(:,1:nlegs),'LO')

      if(allocated(m2arr))deallocate(m2arr)
      allocate(m2arr(size(processes(i_mast)%cols,2)))
      if(allocated(ccamp))deallocate(ccamp)
      allocate(ccamp(size(processes(i_mast)%cols,2)))

      m2arr=0d0
      ccamp=0d0

      do ic=1,size(processes(i_mast)%cols,2)
         do ih=1,size(processes(i_mast)%hels,2)
c     get_amplitude_rcl resolves to the Recola1-compatibility wrapper
c     get_amplitude_r1_rcl, which derives the QED power from the leg
c     count and is restricted to the SM, THDM and HS model files
c     ("Model file SM (QCD) + HEFT not supported"). Call the general
c     routine with the explicit (QCD,QED) powers instead -- the same
c     ones selected in recola_generate_process. Note these are the
c     AMPLITUDE powers, not doubled as for the squared amplitude.
            call get_amplitude_general_rcl(i_mast,
     $           [res_powst+ampshift,res_powew],'LO',
     $           processes(i_mast)%cols(:,ic),
     $           processes(i_mast)%hels(:,ih),amp)
            ccamp(ic)=ccamp(ic)+conjg(amp)*amp
         enddo
         m2arr(ic)=ccamp(ic)
      enddo

c$$$        m2arrCOM=0d0
c$$$        ! Factor 2 account for Recola normalization of 
c$$$        ! Gell-Mann Matrices (a factor of 2 for any external
c$$$        ! gluon) compared to OL. Not relevant for normalized ratio
c$$$        norm=1d0
c$$$        do i=1,size(flav)
c$$$           if(flav(i)==0)norm=norm*2
c$$$        enddo
c$$$        m2arrCOM(:size(m2arr))=m2arr*norm 
      m2arr  = m2arr/sum(m2arr) ! normalize

      toss = random()
      cumm2arr=0
      do i=1,size(m2arr)
         cumm2arr = cumm2arr + m2arr(i)
         if (toss < cumm2arr) then
            randomflow = i
            exit
         end if
      end do

c      print*, 'randflow= ', randomflow
c      print*, 'call color!'

      color(1:2,:)=0
      do i=1,size(flav)
         if(flav(i)==0)then
            ! gluon
            if(i.le.2)then
               color(:,i)=(/i,processes(i_mast)%cols(i,randomflow)/) 
            else
               color(:,i)=(/processes(i_mast)%cols(i,randomflow),i/) 
            endif
         endif
         if(abs(flav(i)).le.6)then
            ! quarks
            if(flav(i).lt.0.and.i.le.2)then
               ! incoming antiquark
               ! cols_Nc contains the gluon position
c               print*, 'i_mast= ', i_mast
c               print*, 'randomflow= ', randomflow
c               print*, 'i= ', i
               color(:,i)=(/0,processes(i_mast)%cols(i,randomflow)/)
c               color(:,i)=(/ 0, 0 /)
            else if(flav(i).gt.0.and.i.gt.2)then
               ! outgoing quark
               ! cols_Nc contains the gluon position
                color(:,i)=(/processes(i_mast)%cols(i,randomflow),0/)
c               color(:,i)=(/ 0, 0 /)
            elseif(flav(i).gt.0.and.i.le.2)then
               ! incoming quark
               ! They have cols_Nc=0 in Recola
               color(:,i)=(/i,0/)
            else if(flav(i).lt.0.and.i.gt.2)then
               ! outgoing antiquark
               ! They have cols_Nc=0 in Recola
               color(:,i)=(/0,i/)
            endif
         endif
      enddo
      
      call map_flow_basis(color,size(flav))
      ifl =  randomflow
      deallocate(m2arr)
        
      contains

      subroutine map_flow_basis(color,n)
      implicit none
      integer, intent(inout) ::  color(2,n)
      integer, intent(in) :: n
      integer :: i, j
      color = color
      do i=1,2
        do j=1,n
          if (color(i,j)/=0) color(i,j)=color(i,j)+500
        end do
      end do
      end subroutine map_flow_basis

      end subroutine recola_colour


      subroutine recola_borncolour(p,bflav,color)
        implicit none
        include "nlegborn.h"
        double precision, intent(in) ::  p(:,:)
        integer, intent(in) :: bflav(:)
        integer, intent(out) ::  color(:,:)
        integer ifl,j
        integer flav(1:nlegbornexternal)
        do j=1,nlegbornexternal
           flav(j)=bflav(j)
        enddo
        call recola_colour(p,flav,color,ifl)
      end subroutine recola_borncolour

      subroutine recola_realcolour(p,rflav,color)
        implicit none
        include "nlegborn.h"
        double precision, intent(in) ::  p(:,:)
        integer, intent(in) :: rflav(:)
        integer, intent(out) ::  color(:,:)
        integer ifl,j
        integer flav(1:nlegrealexternal)
        do j=1,nlegrealexternal
           flav(j)=rflav(j)
        enddo
        call recola_colour(p,flav,color,ifl)
      end subroutine recola_realcolour

      subroutine recola_borncolour_ifl(p,bflav,color,ifl)
        implicit none
        include "nlegborn.h"
        double precision, intent(in) ::  p(:,:)
        integer, intent(in) :: bflav(:)
        integer, intent(out) ::  color(:,:), ifl
        integer j
        integer flav(1:nlegbornexternal)
        do j=1,nlegbornexternal
           flav(j)=bflav(j)
        enddo
        call recola_colour(p,flav,color,ifl)
      end subroutine recola_borncolour_ifl

      subroutine recola_realcolour_ifl(p,rflav,color,ifl)
        implicit none
        include "nlegborn.h"
        double precision, intent(in) ::  p(:,:)
        integer, intent(in) :: rflav(:)
        integer, intent(out) ::  color(:,:), ifl
        integer j
        integer flav(1:nlegrealexternal)
        do j=1,nlegrealexternal
           flav(j)=rflav(j)
        enddo
        call recola_colour(p,flav,color,ifl)
      end subroutine recola_realcolour_ifl


      subroutine remap_momenta(i_mast, flav, p, p_rcl)
      implicit none
      integer, parameter ::maxlength=40
      integer i, j, ip, i_mast
      integer, intent(in) :: flav(:)
      integer flav_tmp(maxlength)
      double precision, intent(in) :: p(:,:)
      double precision, intent(inout) :: p_rcl(:,:)

      if(size(flav).gt.maxlength)then
         write(*,*) " In recola interface, remap_momenta: "
         write(*,*) " flavour string length is too long!!!"
         stop
      endif

      p_rcl(:,1:size(flav))=p(:,1:size(flav))
      ! get first final state parton position
      do i=3,size(flav)
         if(is_parton(flav(i)))then
            ip=i
            exit
         endif
      enddo

      flav_tmp(1:size(flav))=flav(:)
      do i=ip,size(flav)
         do j=ip,size(flav)
            if(processes(i_mast)%flav(i)==flav_tmp(j))then
               flav_tmp(j)=-999
               p_rcl(:,i)=p(:,j)
               exit
            endif
         enddo
      enddo

      end subroutine remap_momenta

      

      subroutine recola_freeze_alphas
c     Fix a_s and the renormalisation scale once, before the processes
c     are generated. Only used for model files without qcd_rescaling
c     (e.g. HEFT). The values are taken from powheg.input:
c        recola_muren   renormalisation scale in GeV
c                       (default (mH + 2 mb)/2, the fixed scale that
c                        Born_phsp.f uses when runningscales 0)
c        recola_alphas  a_s at that scale
c                       (default 0.118)
c     They MUST be consistent with the scale POWHEG itself uses, i.e.
c     run with runningscales 0.
      implicit none
      include 'PhysPars.h'
      include 'pwhg_st.h'
      real * 8 powheginput
      external powheginput
      real * 8 mu0, as0
      character * 4 nfstr
      integer inf

      mu0 = powheginput("#recola_muren")
      if(mu0.le.0d0) mu0 = 0.5d0*(ph_Hmass + 2d0*ph_bmass)
      as0 = powheginput("#recola_alphas")
      if(as0.le.0d0) as0 = 0.118d0

      call set_parameter_rcl('aS',dcmplx(as0))
      call set_mu_ms_rcl(mu0)
c     Number of active flavours in the alpha_s renormalisation
c     counterterm dZgs_QCD2. Defaults to st_nlight (=4 here, the
c     consistent 4FS choice: both b and t decoupled). It can be
c     overridden with recola_nfalphas, because the heavy-quark
c     decoupling logs sitting in dZgs are the same objects that appear
c     in the HEFT Wilson coefficient when converting its alpha_s
c     between flavour schemes -- so this is exactly the bookkeeping
c     ambiguity between our setup and a published one.
      inf = nint(powheginput("#recola_nfalphas"))
      if(inf.lt.3 .or. inf.gt.6) inf = st_nlight
      if(inf.eq.3) then
         nfstr = 'Nf3'
      elseif(inf.eq.4) then
         nfstr = 'Nf4'
      elseif(inf.eq.5) then
         nfstr = 'Nf5'
      else
         nfstr = 'Nf6'
      endif
      call set_renoscheme_rcl('dZgs_QCD2',trim(nfstr))

      als_frozen = as0
      write(*,*) ' Recola: alpha_s frozen to ',as0,' at mu = ',mu0
      write(*,*) ' Recola: the frozen value is arbitrary -- amplitudes'
      write(*,*) '         are rescaled by (alphaS(muR)/als_frozen)**n'
      write(*,*) '         and muMS/muUV/muIR track muR per point.'
      write(*,*) ' Recola: a_s renormalisation scheme ',trim(nfstr)

      end subroutine recola_freeze_alphas


      subroutine recola_check_virtual_poles(p,vflav)
c     Diagnostic: extract the IR pole coefficients of the one-loop
c     amplitude from Recola and compare them with the known QCD values.
c
c     Recola returns  V = c2*DeltaIR2 + c1*DeltaIR + finite, so varying
c     (DeltaIR,DeltaIR2) linearly exposes c1 and c2. In the POWHEG/BLHA
c     normalisation, V = (as/2pi)[c2/eps^2 + c1/eps + c0], the double
c     pole must be
c         c2 = - sum_i C_i   over the MASSLESS coloured legs
c     (massive quarks have no collinear pole), i.e.
c         g g -> H b bbar :  -(CA+CA)   = -6
c         q qbar -> H b bbar: -(CF+CF)  = -8/3
c     A mismatch means the virtual is not correctly normalised or the
c     IR convention does not match POWHEG's.
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_res.h'
      integer, parameter :: nlegs=nlegbornexternal
      double precision, intent(in) :: p(0:3,nlegs)
      integer, intent(in) :: vflav(nlegs)
      double precision p_rcl(0:3,nlegs)
      double precision v00,v10,v01,b0,c1,c2,expect,vuv,vnf
      double precision mms0,vms
      integer inf
      character*4 nfstr
      integer i_mast,i,ncol

      i_mast=get_rcl_index(vflav)
      call remap_momenta(i_mast, vflav, p, p_rcl)
      call set_mu_ir_rcl (dsqrt(st_muren2))
      call set_mu_uv_rcl (dsqrt(st_muren2))
c     muMS is Recola's Qren (the scale alpha_s is defined at). It
c     enters the alpha_s counterterm as beta0*log(muUV^2/muMS^2) and
c     must track muR together with muUV/muIR.
      call set_mu_ms_rcl (dsqrt(st_muren2))

      call set_delta_ir_rcl(0d0,0d0)
      call compute_process_rcl(i_mast,p_rcl(:,1:nlegs),'NLO')
      call get_squared_amplitude_rcl(i_mast,
     $     [2*res_powst,2*res_powew],'LO',b0)
      call get_squared_amplitude_rcl(i_mast,
     $     [2*(res_powst+1),2*res_powew],'NLO',v00)

      call set_delta_ir_rcl(1d0,0d0)
      call compute_process_rcl(i_mast,p_rcl(:,1:nlegs),'NLO')
      call get_squared_amplitude_rcl(i_mast,
     $     [2*(res_powst+1),2*res_powew],'NLO',v10)

      call set_delta_ir_rcl(0d0,1d0)
      call compute_process_rcl(i_mast,p_rcl(:,1:nlegs),'NLO')
      call get_squared_amplitude_rcl(i_mast,
     $     [2*(res_powst+1),2*res_powew],'NLO',v01)

c     restore the BLHA convention used for the physics run
      call set_delta_ir_rcl(0d0,pi**2/6d0)

      c1 = (v10-v00)/b0/(st_alpha/2d0/pi)
      c2 = (v01-v00)/b0/(st_alpha/2d0/pi)

c     expected double pole: -sum of Casimirs of the massless coloured legs
      expect = 0d0
      do i=1,nlegs
         if(vflav(i).eq.0) then
            expect = expect - ca
         elseif(abs(vflav(i)).le.6 .and. abs(vflav(i)).ne.5) then
            expect = expect - cf
         endif
      enddo

c     Can muMS (Recola's Qren, the scale alpha_s is defined at) be
c     changed AFTER generate_processes_rcl? If yes, a fully dynamic
c     renormalisation scale is possible: set muUV=muIR=muMS=muR per
c     point and rescale by (alphas(muR)/als_frozen)**n, which is exact
c     because the amplitudes are homogeneous in alpha_s at fixed muMS.
      call get_mu_ms_rcl(mms0)
      call set_mu_ms_rcl(2d0*mms0)
      call compute_process_rcl(i_mast,p_rcl(:,1:nlegs),'NLO')
      call get_squared_amplitude_rcl(i_mast,
     $     [2*(res_powst+1),2*res_powew],'NLO',vms)
      call set_mu_ms_rcl(mms0)
      call compute_process_rcl(i_mast,p_rcl(:,1:nlegs),'NLO')

      write(*,*) '================ VIRTUAL POLE CHECK ================'
      write(*,*) '  muMS            : ',mms0
      write(*,*) '  V(2*muMS)-V(muMS): ',vms-v00
      write(*,*) '  relative        : ',(vms-v00)/v00,
     $           '   (nonzero => muMS is dynamic)'
      write(*,'(a,10i4)') '  flavours      : ',vflav
      write(*,*) '  Born          : ',b0
      write(*,*) '  V(0,0)        : ',v00
      write(*,*) '  single pole c1: ',c1
      write(*,*) '  double pole c2: ',c2
      write(*,*) '  expected c2   : ',expect
      write(*,*) '  c2/expected   : ',c2/expect
      write(*,*) '  finite V/B    : ',v00/b0/(st_alpha/2d0/pi)
      write(*,*) '  V as POWHEG sees it: ',v00/(st_alpha/2d0/pi)

c     UV-finiteness scan over the number of active flavours in the
c     alpha_s renormalisation. A complete renormalisation -- which for
c     the HEFT includes the counterterm of the effective Hgg Wilson
c     coefficient, dcggh_QCD4 = aS*dcgghfin_QCD2/(12pi)
c     + dZgs_QCD2*gs^2/(24pi^2) -- must give a result INDEPENDENT of
c     DeltaUV. Any residual dependence is a leftover UV pole.
      do inf=4,6
         if(inf.eq.4) nfstr='Nf4'
         if(inf.eq.5) nfstr='Nf5'
         if(inf.eq.6) nfstr='Nf6'
         call set_renoscheme_rcl('dZgs_QCD2',nfstr)
         call set_delta_uv_rcl(0d0)
         call compute_process_rcl(i_mast,p_rcl(:,1:nlegs),'NLO')
         call get_squared_amplitude_rcl(i_mast,
     $        [2*(res_powst+1),2*res_powew],'NLO',vnf)
         call set_delta_uv_rcl(1d0)
         call compute_process_rcl(i_mast,p_rcl(:,1:nlegs),'NLO')
         call get_squared_amplitude_rcl(i_mast,
     $        [2*(res_powst+1),2*res_powew],'NLO',vuv)
         write(*,'(a,a4,a,es14.6,a,es14.6)') '   ',nfstr,
     $        ':  V/B = ',vnf/b0/(st_alpha/2d0/pi),
     $        '   dV(DeltaUV)/B = ',(vuv-vnf)/b0
      enddo
      call set_delta_uv_rcl(0d0)
      call set_renoscheme_rcl('dZgs_QCD2','Nf4')
      write(*,*) '  st_alpha, als_frozen: ',st_alpha,als_frozen
      write(*,*) '===================================================='

      end subroutine recola_check_virtual_poles

      double precision function recola_asfact(npow)
c     Restores a running alpha_s. The HEFT model file has no
c     qcd_rescaling, so set_alphas_rcl is refused and Recola evaluates
c     everything at als_frozen. At fixed muMS the squared amplitudes
c     are exactly homogeneous in alpha_s, so this rescaling is exact.
c     npow: res_powst for the Born, res_powst+1 for real and virtual.
      implicit none
      include 'pwhg_st.h'
      integer npow
      if(als_rescaling) then
         recola_asfact = 1d0
      else
         recola_asfact = (st_alpha/als_frozen)**npow
      endif
      end function recola_asfact

      end module recola_powheg




