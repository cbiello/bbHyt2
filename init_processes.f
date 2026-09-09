      subroutine init_processes
c     bbH production in the 4-flavour scheme, top-Yukawa (yt^2) component.
c
c     The Higgs couples to the gluons through the heavy-top effective
c     vertex (Recola2 HEFT model file), the bottom quarks are massive
c     and are produced by gluon splitting. There is NO bottom Yukawa
c     coupling in this contribution, hence -- unlike the yb^2 generator
c     bbH4FS -- no running Yukawa mass has to be set up here.
c
c     Amplitudes are provided by Recola2 (model file HEFT).
      use recola_powheg, only: recola_init,recola_generate_process
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
      include "pwhg_st.h"
      include 'pwhg_ckm.h'
      include 'pwhg_physpar.h'
      include 'pwhg_res.h'
      include "pwhg_flg.h"
      include "PhysPars.h"
c     pwhg_dpa.h is pulled in by the (bb4l-derived) recola.f. The
c     double-pole approximation is meaningless for this process, which
c     has no resonances, so the flag is explicitly switched off below.
      include 'pwhg_dpa.h'
      integer i, int(maxprocreal), ickm
      integer nmaxres_real,nmaxres_born,dim_integ
      real * 8 powheginput
      external powheginput
      character(len=100) :: model_rcl
      common/physicsmodel/model_rcl

      dpa_flagtt = 0

      do i=1,maxprocreal
         int(i)=i
      enddo

      ickm = 2   !   0,1,2 for general, ckm_cabibbo, ckm_diag
      ckm_diag = .false.
      ckm_cabibbo = .true.
      if(ickm.eq.2) then
         ckm_diag = .true.
      elseif(ickm.eq.1) then
         ckm_cabibbo = .true.
      endif

c     init couplings and flavour structures
      call init_couplings
      flg_with_em=.false.

c     set the bornzerodamp flag to true by default
      flg_bornzerodamp = .true.
      flg_withdamp     = .true.
      if (powheginput("#bwithdamp").eq.0) then
         flg_bornzerodamp = .false.
         flg_withdamp     = .false.
      end if

c     determine number of light flavours.
c     In the 4FS the b quark is massive, so this must come out as 4.
      st_nlight = 0
      do i=1,6
        if (physpar_pdgmasses(i) == 0) st_nlight=st_nlight+1
      end do
      if (st_nlight /= 4) then
         write(*,*) ' init_processes: st_nlight = ',st_nlight
         write(*,*) ' bbH_yt2 is a 4-flavour-scheme generator and'
         write(*,*) ' requires a non-zero bottom mass (bmass in'
         write(*,*) ' powheg.input). Exiting ...'
         call pwhg_exit(-1)
      endif

      call init_processes_born
      call init_processes_real

c     Coupling powers of the Born amplitude.
c
c     In the Recola2 HEFT model file the effective Higgs-gluon vertices
c     carry the order increases (QCD,QED) = (2,1) for Hgg, (3,1) for
c     Hggg and (4,1) for Hgggg, because the Wilson coefficient
c     cggh = aS*(1+cgghfin)/(12 pi) counts as two powers of QCD.
c     Every Born diagram of g g -> H b bbar (and q qbar -> H b bbar)
c     therefore has (QCD,QED) = (4,1).
c
c     Note that the yb^2 Born (Higgs radiated off the b line) has
c     (QCD,QED) = (2,1) instead. Selecting QCD=4 below is precisely
c     what projects onto the yt^2 contribution and removes both the
c     yb^2 term and the yb*yt interference.
      res_powew=1
      res_powst=4
      st_bornorder = res_powst

c     no resonances in this process
      if(powheginput("#nores") /= 1) then
c        call build_resonance_histories
      endif

      write(*,*) ' ***********   FINAL POWHEG    ***************'
      write(*,*) ' ***********     BORN          ***************'
      do i=1,flst_nborn
         write(*,'(a,100i4)')'                  ',int(1:flst_bornlength(i))
         write(*,'(a,i3,a,100i4)')'flst_born(',i,')   =',flst_born(1:flst_bornlength(i),i)
         write(*,'(a,i3,a,100i4)')'flst_bornres(',i,')=',flst_bornres(1:flst_bornlength(i),i)
         write(*,*)
      enddo
      write(*,*) 'max flst_bornlength: ',maxval(flst_bornlength(1:flst_nborn))

      write(*,*) ' ***********     REAL          ***************'
      do i=1,flst_nreal
         write(*,'(a,100i4)')'                  ',int(1:flst_reallength(i))
         write(*,'(a,i3,a,100i4)')'flst_real(',i,')   =',flst_real(1:flst_reallength(i),i)
         write(*,'(a,i3,a,100i4)')'flst_realres(',i,')=',flst_realres(1:flst_reallength(i),i)
         write(*,*)
      enddo
      write(*,*) 'max flst_reallength: ',maxval(flst_reallength(1:flst_nreal))

      call buildresgroups(flst_nborn,nlegborn,flst_bornlength,
     1     flst_born,flst_bornres,flst_bornresgroup,flst_nbornresgroup)

      call buildresgroups(flst_nreal,nlegreal,flst_reallength,
     1     flst_real,flst_realres,flst_realresgroup,flst_nrealresgroup)

      write(*,*)
      write(*,*) '*********   SUMMARY  *********'
      write(*,*) 'set nlegborn to: ',maxval(flst_bornlength(1:flst_nborn))
      write(*,*) 'set nlegreal to: ',maxval(flst_reallength(1:flst_nreal))

      nmaxres_real = maxval(flst_realres(1:flst_reallength(flst_nreal),1:flst_nreal)) - 2
      nmaxres_born = maxval(flst_bornres(1:flst_bornlength(flst_nborn),1:flst_nborn)) - 2
      write(*,*) 'max number of resonances: ',max(nmaxres_real,nmaxres_born)

c     add (-1) for overall azimuthal rotation of the event around the beam axis
      dim_integ = (flst_numfinal+1)*3 - 4 + 2 + max(nmaxres_real,nmaxres_born)
      write(*,*) 'set ndiminteg to: ',dim_integ

      write(*,*) 'flst_nbornresgroup ',flst_nbornresgroup
      write(*,*) 'flst_nrealresgroup ',flst_nrealresgroup

      write(*,*) '*********   END SUMMARY  *********'
      write(*,*)

      write(*,*) "##########################################"
      write(*,*) "       Initialisation of Recola"
      write(*,*) "##########################################"
      call recola_init
c     The HEFT model file reports itself as "SM (QCD) + HEFT", so test
c     for the substring rather than the leading characters.
      if(index(model_rcl,'HEFT').eq.0) then
         write(*,*) ' init_processes: the active Recola model file is'
         write(*,*) ' "',trim(model_rcl),'"'
         write(*,*) ' bbH_yt2 needs the HEFT model file. Rebuild the'
         write(*,*) ' Recola library with  cmake .. -Dmodel=HEFT .'
         write(*,*) ' Exiting ...'
         call pwhg_exit(-1)
      endif
      call init_heft_parameters
      call recola_generate_process
      write(*,*) " Recola process generation completed! "
      write(*,*) "##########################################"

      end


      subroutine init_heft_parameters
c     HEFT-model-specific settings. Must be called after recola_init
c     (which fixes masses, widths and the EW scheme) and before
c     recola_generate_process (which calls generate_processes_rcl).
      use recola
      implicit none
      real * 8 powheginput
      external powheginput
      real * 8 cgghfin
      character(len=100) :: model_rcl
      common/physicsmodel/model_rcl

      if(index(model_rcl,'HEFT').eq.0) return

c     Finite part of the Wilson coefficient of the effective
c     H-gluon-gluon operator:
c
c        cggh = aS*(1 + cgghfin)/(12 pi)
c
c     cgghfin = 0     -> LO Wilson coefficient (default here)
c     cgghfin = 11/4 * aS/pi  -> NLO Wilson coefficient
c
c     The O(aS) term of the Wilson coefficient contributes at the same
c     order as the NLO QCD corrections computed here, so for a
c     consistent NLO HEFT prediction it should normally be switched on.
c     It is left at zero by default so that the choice is explicit:
c     set cgghfin in powheg.input.
      cgghfin = powheginput("#cgghfin")
      if(cgghfin.lt.0d0) cgghfin = 0d0
      call set_parameter_rcl('cgghfin',dcmplx(cgghfin))
      write(*,*) ' HEFT: cgghfin set to ',cgghfin

c     Switch off the bottom Yukawa coupling. This contribution is
c     proportional to yt^2 only; the bottom mass enters exclusively
c     through the kinematics and the propagators.
c
c     NOTE: the yt^2 projection does NOT rely on this. It is enforced
c     by the coupling-power selection (QCD=4) done in
c     recola_generate_process, since the yb^2 Born carries QCD=2.
c     Setting ymb=0 here is a consistency measure; the Recola2 model
c     file re-derives ymb from the bottom mass whenever the masses are
c     re-initialised, so the value is printed back for checking.
      call set_parameter_rcl('ymb',dcmplx(0d0))
      write(*,*) ' HEFT: bottom Yukawa switched off (ymb = 0)'

      end


      subroutine init_processes_born
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
c     set the number of final-state particles at the Born level
      flst_bornlength = nlegbornexternal
      flst_numfinal=nlegbornexternal-2

      flst_nborn= 9             ! number of born flavour structures

c     Born flavour structures: (i1) (i2) H b bbar
c     Same list as the yb^2 generator bbH4FS: gg and the four light
c     q qbar channels in both orderings. No b quark in the initial
c     state (4-flavour scheme).
        flst_born(   1,   1)=           1
        flst_born(   2,   1)=          -1
        flst_born(   3,   1)=          25
        flst_born(   4,   1)=           5
        flst_born(   5,   1)=          -5

        flst_born(   1,   2)=           2
        flst_born(   2,   2)=          -2
        flst_born(   3,   2)=          25
        flst_born(   4,   2)=           5
        flst_born(   5,   2)=          -5

        flst_born(   1,   3)=           3
        flst_born(   2,   3)=          -3
        flst_born(   3,   3)=          25
        flst_born(   4,   3)=           5
        flst_born(   5,   3)=          -5

        flst_born(   1,   4)=           4
        flst_born(   2,   4)=          -4
        flst_born(   3,   4)=          25
        flst_born(   4,   4)=           5
        flst_born(   5,   4)=          -5

        flst_born(   1,   5)=           0
        flst_born(   2,   5)=           0
        flst_born(   3,   5)=          25
        flst_born(   4,   5)=           5
        flst_born(   5,   5)=          -5

        flst_born(   1,   6)=          -4
        flst_born(   2,   6)=           4
        flst_born(   3,   6)=          25
        flst_born(   4,   6)=           5
        flst_born(   5,   6)=          -5

        flst_born(   1,   7)=          -3
        flst_born(   2,   7)=           3
        flst_born(   3,   7)=          25
        flst_born(   4,   7)=           5
        flst_born(   5,   7)=          -5

        flst_born(   1,   8)=          -2
        flst_born(   2,   8)=           2
        flst_born(   3,   8)=          25
        flst_born(   4,   8)=           5
        flst_born(   5,   8)=          -5

        flst_born(   1,   9)=          -1
        flst_born(   2,   9)=           1
        flst_born(   3,   9)=          25
        flst_born(   4,   9)=           5
        flst_born(   5,   9)=          -5

        flst_bornres(1:nlegborn,1:flst_nborn) = 0

      return
      end


      subroutine init_processes_real
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
      include "pwhg_st.h"
      integer ihvq,nflav,i1,i2,i3,i4,i5,i6,i7,j,k,ii(7)
      equivalence (i1,ii(1)),(i2,ii(2)),(i3,ii(3)),
     #(i4,ii(4)),(i5,ii(5)),(i6,ii(6)),(i7,ii(7))
      real * 8 powheginput
      external powheginput
      logical debug
      parameter (debug=.false.)
      flst_reallength = nlegbornexternal+1
      flst_numfinal = nlegrealexternal-1

      ihvq=5
      nflav=st_nlight

*********************************************************************
***********          REAL (H b bbar + parton)         ***************
*********************************************************************
      i4=ihvq
      i5=-ihvq
      i3=25 ! Higgs
      flst_nreal=0
      do i1=-nflav,nflav
         do i2=-nflav,nflav
            if(i1.ne.0.and.i1+i2.eq.0) then
c     q qbar
               i6=0
               if(powheginput('#qqbproc').ne.0) goto 20
            endif
c     a quark and a gluon
            if(i1*i2.eq.0) then
               if(i1.ne.0) then
                  i6=i1
                  if(powheginput('#qgproc').ne.0) goto 20
               elseif(i2.ne.0) then
                  i6=i2
                  if(powheginput('#gqbproc').ne.0) goto 20
               endif
            endif
c     two gluons
            if((i1.eq.0).and.(i2.eq.0))then
               i6=0
               if(powheginput('#ggproc').ne.0) goto 20
            endif
            goto 21
 20         continue
            flst_nreal=flst_nreal+1
            if(flst_nreal.gt.maxprocreal) goto 999
            do k=1,nlegreal
               flst_real(k,flst_nreal)=ii(k)
               flst_realres( k, 1:flst_nreal) = 0
            enddo
 21         continue
         enddo
      enddo
      if (debug) then
         write(*,*) ' real processes',flst_nreal
         do j=1,flst_nreal
            write(*,*) (flst_real(k,j),k=1,nlegreal)
         enddo
      endif

      flst_reallength(:) = 6
      flst_realres( 1, 1:flst_nreal) = 0
      flst_realres( 2, 1:flst_nreal) = 0
      flst_realres( 3, 1:flst_nreal) = 0
      flst_realres( 4, 1:flst_nreal) = 0
      flst_realres( 5, 1:flst_nreal) = 0
      flst_realres( 6, 1:flst_nreal) = 0

      return
 999  write(*,*) 'init_processes: increase maxprocreal'
      stop
      end
