      subroutine born_phsp(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'pwhg_flst.h'
      include 'pwhg_physpar.h'
      include 'PhysPars.h'
      real * 8  :: p(0:3,nlegborn)
      real * 8  :: cmp(0:3,nlegborn)
      real * 8  :: masses(nlegborn)
      real * 8  :: xborn(*)
      real * 8  :: powheginput
      integer   :: iborn,j
      logical, save :: ini=.true., nores
      real * 8 :: ptH, ptHmin

      if(ini) then
         if(powheginput("#nores") == 1) then
            nores = .true.
         else
            nores = .false.
         endif
         ini = .false.
      endif

      if(nores) then
c provide Born phase space for resonance unaware integration

      else
c generate Born phase space from resoance information
         do iborn=1,flst_nborn
            if(flst_bornresgroup(iborn).eq.flst_ibornresgroup) exit
         enddo

         flst_ibornlength = flst_bornlength(iborn)
         flst_ireallength = flst_ibornlength + 1

         do j=1,flst_ibornlength
            kn_masses(j) = physpar_phspmasses(flst_born(j,iborn))
         enddo

         call genphasespace(xborn,flst_ibornlength,flst_born(:,iborn),
     1        flst_bornres(:,iborn),kn_beams,kn_jacborn,
     1        kn_xb1,kn_xb2,kn_sborn,kn_cmpborn,kn_pborn)

         kn_minmass = 2*ph_bmass
      endif

      ptHmin=powheginput("#BornPTHmin")
      ptH = dsqrt(kn_pborn(1,3)**2+kn_pborn(2,3)**2)
      if(ptH.lt.ptHmin) then
           kn_jacborn = 0d0
      endif

      end

      subroutine born_suppression(fact)
      implicit none
      real * 8 fact
      fact=1d0
      end

      subroutine rmn_suppression(fact)
      implicit none
      real * 8 fact
      fact=1d0
      end


      subroutine regular_suppression(fact)
      implicit none
      real * 8 fact
      call rmn_suppression(fact)
      end


      subroutine global_suppression(c,fact)
      implicit none
      character * 1 c
      real * 8 fact
      fact=1d0
      end



      subroutine set_fac_ren_scales(muf,mur)
      implicit none
      include 'PhysPars.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'pwhg_physpar.h'
      
      real*8 phiggs(0:3),pb(0:3),pbbar(0:3), p(0:3)
      real*8 mt_higgs,mt_b,mt_bbar,pt_i
      real * 8, intent(out) ::  muf,mur
      logical, save        :: ini=.true., nores
      real *8              :: mu0, MbbH
      real *8              :: dotp,powheginput
      real * 8 renscfact,facscfact
      integer  :: i
      external dotp,powheginput
      logical, save        :: fixedscale,runningscales


      logical gscale,hscale,Mscale,HT4scale
      save gscale, hscale,Mscale,HT4scale
      integer nfs,ileg
      real*8 ptsq,mm,htsum

      if (ini) then
         runningscales = .false.
         gscale = .false.
         hscale = .false.
         Mscale = .false.
         HT4scale = .false.

         renscfact=powheginput("#renscfact")
         facscfact=powheginput("#facscfact")

          if(powheginput('#runningscales').eq.0) then
            runningscales=.false.
            hscale = .false.
            gscale = .false.
          endif
         
         if(powheginput('#runningscales').eq.1) then
            runningscales=.true.   
            gscale = .true.
         elseif(powheginput('#runningscales').eq.2) then
            runningscales=.true.   
            hscale = .true.
         elseif(powheginput('#runningscales').eq.3) then
            runningscales=.true.
            Mscale = .true.
         elseif(powheginput('#runningscales').eq.4) then
c     HT/4 over all final-state particles, including the emitted
c     parton. Central scale of 2307.09992 eq.(3).
c     NOTE: needs "btlscalereal 1" in powheg.input, otherwise POWHEG
c     never sets flg_btildepart='r' and the real uses the Born HT.
            runningscales=.true.
            HT4scale = .true.
         endif   

         write(*,*) '*************************************'
         write(*,*) 'Factorization and renormalization '
         write(*,*) 'scales for the Born alphs (mur, muf) set to '
         if (gscale) then 
            write(*,*) '[ mT(Higgs) + mT(b) + mT(bbar) ]/2'
         elseif (hscale) then 
            write(*,*) '[ mT(Higgs) + mT(b) + mT(bbar) + pT(i) ]/2'
         elseif (Mscale) then 
            write(*,*) 'MbbH'
         elseif (HT4scale) then
            write(*,*) 'HT/4 over all final-state particles'
         else
            write(*,*) '(mH+2mb)/2'
         endif   
         if (renscfact .gt. 0d0) 
     &        write(*,*) 'Renormalization scale rescaled by', renscfact
         if (facscfact .gt. 0d0) 
     &        write(*,*) 'Factorization scale rescaled by  ', facscfact
         write(*,*) '***********************************************'

         if(powheginput("#nores") == 1) then
            nores = .true.
         else
            nores = .false.
         endif
         ini=.false.

      endif

c     Fixed-scale mode. Must be assigned on EVERY call: muf/mur are
c     intent(out) dummies, and setting them only inside "if (ini)"
c     (as in bbH4FS) leaves them undefined from the second call on.
      if (.not.runningscales) then
         muf = 0.5d0*(ph_Hmass + 2d0*ph_bmass)
         mur = muf
         return
      endif

      if (runningscales) then 

         MbbH = 0d0
         p(:) = 0d0
         if(flg_btildepart.eq.'r') then
            phiggs(:) = kn_preal(:,3)   !Higgs
            pb(:)     = kn_preal(:,4)   !b
            pbbar(:)  = kn_preal(:,5)   !bbar

          do i=3,5
            p(:) = p(:) + kn_preal(:,i)
          enddo
         else   
c take Born momenta:
            phiggs(:) = kn_pborn(:,3)   !Higgs
            pb(:)     = kn_pborn(:,4)   !b
            pbbar(:)  = kn_pborn(:,5)   !bbar
          do i=3,5
            p(:) = p(:) + kn_pborn(:,i)
          enddo
         endif

        
         if (Mscale) then
         mu0=sqrt(p(0)**2-p(1)**2-p(2)**2-p(3)**2) !The invariant mass of the bbH system
         muf = mu0 
         mur = mu0
         endif

c        endif
     
         if (HT4scale) then
            htsum = 0d0
            if(flg_btildepart.eq.'r') then
               nfs = flst_ireallength
            else
               nfs = flst_ibornlength
            endif
            do ileg=3,nfs
               if(flg_btildepart.eq.'r') then
                  ptsq = kn_preal(1,ileg)**2+kn_preal(2,ileg)**2
                  mm   = physpar_phspmasses(flst_real(ileg,1))
               else
                  ptsq = kn_pborn(1,ileg)**2+kn_pborn(2,ileg)**2
                  mm   = physpar_phspmasses(flst_born(ileg,1))
               endif
               htsum = htsum + dsqrt(ptsq + mm**2)
            enddo
            muf = htsum/4d0
            mur = muf
            return
         endif

         call gettransmass(ph_Hmass,phiggs,mt_higgs)
         call gettransmass(ph_bmass,pb,mt_b)
         call gettransmass(ph_bmass,pbbar,mt_bbar)

         if (gscale) then
            muf = 0.5d0*(mt_higgs+mt_b+mt_bbar)
         elseif (hscale) then 
            muf = 0.5d0*(mt_higgs+mt_b+mt_bbar + pt_i)
         endif

         mur = muf

      endif


      contains

      subroutine no_scales
      write(*,*) "Error in scale setting, exiting ..."
      call pwhg_exit(-1)
      end subroutine no_scales

      end

      subroutine gettransmass(m,p,mt)
      implicit none
      real * 8 m,p(0:3),mt
      mt=dsqrt(abs(m**2+p(1)**2+p(2)**2))
      end
