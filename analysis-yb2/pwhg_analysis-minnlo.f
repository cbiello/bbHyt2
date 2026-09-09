c  The next subroutines open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'
      character * 9 pr
      character * 6 suffix 
      common/pwhgprocess/pr
      character * 6 whcprg      
      common/cwhcprg/whcprg
      real * 8 dy,dphi,dpt,dr,dptzoom,dM,dmm,dpt3
      real*8 ptmin,ptmax,mmin,mmax
      real*8 ymin,ymax,phimin,phimax,rmin,rmax
      integer j,i

      call inihists
         
      pr='inc-bbH'

      dy=0.2d0
      dphi=0.2d0
      dpt=10d0
      dM=20d0
      dptzoom=1d0
      dmm = 5d0
      dr = 0.2d0

      ptmin = 0d0
      ptmax = 800d0

      ymin = -4d0
      ymax = 4d0

      phimin = 0d0
      phimax = 3.2d0

      mmin = 0d0
      mmax = 1000d0

      rmin = 0d0
      rmax = 5d0

cccccccccccccccccccccccccc
c
c inclusive setup:
c
c incl. xsec
      call bookupeqbins('xsec',1d0,0d0,1d0)
      call bookupeqbins('xsec-ptHge10',1d0,0d0,1d0)
      call bookupeqbins('xsec-ptHge15',1d0,0d0,1d0)
      call bookupeqbins('xsec-ptHge20',1d0,0d0,1d0)

c Higgs
      call bookupeqbins('pt_Higgs',dpt,ptmin,ptmax)
      call bookupeqbins('ptzoom_Higgs',dpt*0.1d0,ptmin,ptmax*0.1d0)
      call bookupeqbins('eta_Higgs',dy,ymin,ymax)
      call bookupeqbins('y_Higgs',dy,ymin,ymax)
      call bookupeqbins('mass_Higgs',dmm,0d0,200d0)

c Higgs-bottom pair:
      call bookupeqbins('deta_Higgs_b',dy,ymin,ymax)
      call bookupeqbins('dy_Higgs_b',dy,ymin,ymax)
      call bookupeqbins('dphi_Higgs_b',dphi,phimin,phimax)
      call bookupeqbins('Reta_Higgs_b',dr,rmin,rmax)
      call bookupeqbins('Ry_Higgs_b',dr,rmin,rmax)

c Higgs-bbar pair:
      call bookupeqbins('deta_Higgs_bbar',dy,ymin,ymax)
      call bookupeqbins('dy_Higgs_bbar',dy,ymin,ymax)
      call bookupeqbins('dphi_Higgs_bbar',dphi,phimin,phimax)
      call bookupeqbins('Reta_Higgs_bbar',dr,rmin,rmax)
      call bookupeqbins('Ry_Higgs_bbar',dr,rmin,rmax)

c Higgs-b-bbar system:
      call bookupeqbins('pt_Higgsbbbar',dpt,ptmin,1d3)
      call bookupeqbins('ptzoom_Higgsbbbar',dpt*0.1d0,ptmin,1d2)
      call bookupeqbins('mass_Higgsbbbar',dmm,mmin,mmax)

c bottom and anti-bottom:
c pt:
      call bookupeqbins('pt_b',dpt,ptmin,ptmax)
      call bookupeqbins('ptzoom_b',dpt*0.1d0,ptmin,ptmax*0.1d0)
      call bookupeqbins('pt_bbar',dpt,ptmin,ptmax)
      call bookupeqbins('ptzoom_bbar',dpt*0.1d0,ptmin,ptmax*0.1d0)

c eta:
      call bookupeqbins('eta_b',dy,ymin,ymax)
      call bookupeqbins('eta_bbar',dy,ymin,ymax)

c y: 
      call bookupeqbins('y_b',dy,ymin,ymax)
      call bookupeqbins('y_bbar',dy,ymin,ymax)

c bbbar pair:
      call bookupeqbins('mass_bbbar',dM,mmin,mmax)
      call bookupeqbins('pt_bbbar',dpt,ptmin,ptmax)
      call bookupeqbins('ptzoom_bbbar',dpt*0.1d0,ptmin,ptmax*0.1d0)
c
c b,bbar separation:
      call bookupeqbins('deta_b_bbar',dy,ymin,ymax)
      call bookupeqbins('dy_b_bbar',dy,ymin,ymax)
      call bookupeqbins('dphi_b_bbar',dphi,phimin,phimax)
      call bookupeqbins('Reta_b_bbar',dr,rmin,rmax)     
      call bookupeqbins('Ry_b_bbar',dr,rmin,rmax)     

c bbbar, H separation:
      call bookupeqbins('Reta_bbbar_Higgs',dr,rmin,rmax)
      call bookupeqbins('Ry_bbbar_Higgs',dr,rmin,rmax)
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
      end

      subroutine analysis(dsig0)
      implicit none
      real * 8 dsig0
      real * 8, allocatable, save :: dsig(:)
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include  'LesHouches.h'
      include 'pwhg_weights.h'
      include 'pwhg_rwl.h'
      logical ini
      data ini/.true./
      save ini
      character * 9 pr
      common/pwhgprocess/pr
      real*8 ppairbbbar(4),ptb,yb,etab,ptbbar,ybbar,etabbar,dphi
      real*8 rbbbar,phibbbar,ybbbar,etabbbar,ptbbbar,mbbbar
      real*8 pth,etah,yh,etabh,ybh,phibh,rbh,mhiggs
      real*8 etabbarh,ybbarh,phibbarh,rbbarh
      real*8 ybbh,etabbh,phibbh,rbbh
      real*8 pbbh(4),ptbbh,mbbh
      integer i

      real*8 retabbbar,retabh,retabbarh,retabbh
      real * 8 binsize(700)
      common/pwhghistcommon/binsize
      integer ihep,ib,ibbar,ihiggs,mu
      character * 6 whcprg      
      common/cwhcprg/whcprg
      data whcprg/'NLO   '/

      integer iib,iibbar,iihiggs
      logical write_cuts
      
c================================================

      write_cuts=.false.     
      
      if(ini) then

         write(*,*) ''
         write(*,*) '*****************************'
         write(*,*) '** weights_num     = ',weights_num
         write(*,*) '** rwl_num_weights = ',rwl_num_weights
         write(*,*) '** rwl_num_groups = ',rwl_num_groups
         write(*,*) '*****************************'
         write(*,*) ''

         write_cuts=.true.

         if(weights_num.eq.0.and.rwl_num_weights.eq.0) then
            call setupmulti(1)
         else if(weights_num.ne.0.and.rwl_num_weights.eq.0) then
            call setupmulti(weights_num)
         else if(weights_num.eq.0.and.rwl_num_weights.ne.0) then
            call setupmulti(rwl_num_weights)
         else
            call setupmulti(rwl_num_weights)
         endif

         if(.not. allocated(dsig)) then
            allocate(dsig(max(10,rwl_num_weights+1)))
         endif
          
         ini=.false.
      endif

      dsig(:)=0d0

      if(weights_num.eq.0.and.rwl_num_weights.eq.0) then
         dsig(1)=dsig0
      else if(weights_num.ne.0.and.rwl_num_weights.eq.0) then
         dsig(1:weights_num)=weights_val(1:weights_num)
      else if(weights_num.eq.0.and.rwl_num_weights.ne.0) then
         dsig(1:rwl_num_weights)=rwl_weights(1:rwl_num_weights)
      else
         dsig(1:rwl_num_weights)=rwl_weights(1:rwl_num_weights)
      endif

      if(sum(abs(dsig)).eq.0) return
   
C     ---------------------------------------------------------------
C     ---------------------------------------------------------------

      if(write_cuts) then
         write(*,*) '********************************************'
         write(*,*) '********************************************'
         write(*,*) '                ANALYSIS CUTS               '
         write(*,*) '********************************************'
         write(*,*) '********************************************'
         write(*,*) ''
         write(*,*) 'INCLUSIVE ANALYSIS:'
         write(*,*) 'no cuts' 
         write(*,*) ''
         write(*,*) '*******************************************'
         write(*,*) '*******************************************'
      endif

      ib=0
      iib=0
      ibbar=0
      iibbar=0
      ihiggs=0
      iihiggs=0

c Parton level analysis
      if(whcprg.eq.'NLO'.or.whcprg.eq.'LHE') then
         do ihep=3,nhep
            if(idhep(ihep).eq.5) then
               ib=ihep
               iib=iib+1
            elseif(idhep(ihep).eq.-5) then
               ibbar=ihep
               iibbar=iibbar+1
            elseif(idhep(ihep).eq.25) then
               ihiggs=ihep
               iihiggs=iihiggs+1
            endif
         enddo !ihep

c check that exactly 1 bottom, 1 anti-bottom, and 1 Higgs are selected
         if (iib.ne.1) then
            write(*,*) "Error in pwhg_analysis: ",iib," bottoms"
c            call printleshouches
            call exit(1)
         endif
         if (iibbar.ne.1) then
            write(*,*) "Error in pwhg_analysis: ",iibbar," anti-bottoms"
c            call printleshouches
            call exit(1)
         endif
         if (iihiggs.ne.1) then
            write(*,*) "Error in pwhg_analysis: ",iihiggs,"Higgses"
c            call printleshouches
            call exit(1)
         endif

c Hadron level analysis
      else

         if (WHCPRG.eq.'HERWIG') then 
            
         do ihep=1,nhep

c searching stategies according to the shower Monte Carlo
c     program used      
c simplistic analysis:
!            print*, 'Check isthep if you want to user HERWIG for bbH analysis' 
            if (isthep(ihep).eq.155) then
               if (idhep(ihep).eq.25) then 
                  ihiggs=ihep
                  iihiggs=iihiggs+1
               endif   
               if (idhep(ihep).eq.5) then !bottom
                  ib = ihep
                  iib = iib+1
               endif   
               if (idhep(ihep).eq.-5) then !bbar
                  ibbar = ihep
                  iibbar = iibbar+1
               endif 
            endif !isthep
            
         enddo !ihep

c         elseif (WHCPRG.eq.'PYTHIA') then
          elseif (WHCPRG.eq.'PY8') then
            
         do ihep=1,nhep
            if(idhep(ihep).eq.25) then                
               if(isthep(ihep).eq.1) then !'stable' Higgs
                  ihiggs=ihep
                  iihiggs=iihiggs+1
               elseif(isthep(ihep).eq.2) then ! particle has decayed
                  ihiggs=ihep
                  iihiggs=iihiggs+1
               endif
            endif !idhep
            if (idhep(ihep).eq.5) then !bottom
               ib = ihep
               iib = iib+1
            endif   
            if (idhep(ihep).eq.-5) then !bbar
               ibbar = ihep
               iibbar = iibbar+1
            endif   
         enddo !nhep


         endif !pythia or herwig

c check that at least one top and one anti-top are selected
c during showering there may be more than one entry with PDG=5,
c with different IST. 
         if (iib.lt.1) then
            write(*,*) "Error in pwhg_analysis: ",iib," bottoms"
c            call printleshouches
c            call exit(1)
         endif
         if (iibbar.lt.1) then
            write(*,*) "Error in pwhg_analysis: ",iibbar," anti-bottoms"
c            call printleshouches
            call exit(1)
         endif

      endif                     !parton or hadron level analysis

      if(rwl_num_weights.eq.0) then
         if(dsig0+1 .eq. dsig0) then
             write(*,*) "LARGE weight. DISCARDING EVENT, weight = ",dsig0
             return
         endif
      else
         do i=1,rwl_num_weights
           if(abs(dsig(i))>1d3 .or. dsig(i)+1 .eq. dsig(i)) then
             write(*,*) "LARGE weight. DISCARDING EVENT, i, weight = ",i, dsig(i)
             return
           endif
         enddo
      endif 


c bottoms:
      call ptyeta(phep(1,ib),ptb,yb,etab)
      call ptyeta(phep(1,ibbar),ptbbar,ybbar,etabbar)

      call getdydetadphidr(phep(1,ib),phep(1,ibbar),
     %     ybbbar,etabbbar,phibbbar,rbbbar)

      call getdydetadphidrr(phep(1,ib),phep(1,ibbar),
     %     retabbbar)

c Higgs:
      call ptyeta(phep(1,ihiggs),pth,yh,etah)

      mhiggs = dsqrt(abs(phep(4,ihiggs)**2-phep(1,ihiggs)**2-
     &                   phep(2,ihiggs)**2-phep(3,ihiggs)**2))

      call getdydetadphidr(phep(1,ib),phep(1,ihiggs),
     %     ybh,etabh,phibh,rbh)
      call getdydetadphidr(phep(1,ibbar),phep(1,ihiggs),
     %     ybbarh,etabbarh,phibbarh,rbbarh)

      call getdydetadphidrr(phep(1,ib),phep(1,ihiggs),
     %     retabh)
      call getdydetadphidrr(phep(1,ibbar),phep(1,ihiggs),
     %     retabbarh)

c mass of the pair
      do mu=1,4
         ppairbbbar(mu)=phep(mu,ib)+phep(mu,ibbar)
      enddo
      mbbbar=dsqrt(abs(ppairbbbar(4)**2
     1            -ppairbbbar(1)**2-ppairbbbar(2)**2-ppairbbbar(3)**2))
      ptbbbar=dsqrt(abs(ppairbbbar(1)**2+ppairbbbar(2)**2))

c bbbar+Higgs system:
      pbbh(1:4) = ppairbbbar(1:4)+phep(1:4,ihiggs)
      ptbbh = dsqrt(abs(pbbh(1)**2+pbbh(2)**2))
      mbbh = dsqrt(abs(pbbh(4)**2-pbbh(1)**2-pbbh(2)**2-pbbh(3)**2))

c bbbar vs Higgs:
      call getdydetadphidr(ppairbbbar(1),phep(1,ihiggs),
     %     ybbh,etabbh,phibbh,rbbh)

      call getdydetadphidrr(ppairbbbar(1),phep(1,ihiggs),
     %     retabbh)


cccccccccccccccccccccccccccccc
c
c inclusive analysis:
c
c xsec:
      call filld('xsec',0.5d0,dsig)


      if(pth.ge.10) call filld('xsec-ptHge10',0.5d0,dsig)
      if(pth.ge.15) call filld('xsec-ptHge15',0.5d0,dsig)
      if(pth.ge.20) call filld('xsec-ptHge20',0.5d0,dsig)
      
c Higgs:
      call filld('pt_Higgs',pth,dsig)
      call filld('ptzoom_Higgs',pth,dsig)
      call filld('eta_Higgs',etah,dsig)
      call filld('y_Higgs',yh,dsig)
      call filld('mass_Higgs',mhiggs,dsig)

c H-b pair:
      call filld('deta_Higgs_b',etabh,dsig)
      call filld('dy_Higgs_b',ybh,dsig)
      call filld('dphi_Higgs_b',phibh,dsig)
      call filld('Reta_Higgs_b',retabh,dsig)
      call filld('Ry_Higgs_b',rbh,dsig)

c H-bbar pair:
      call filld('deta_Higgs_bbar',etabbarh,dsig)
      call filld('dy_Higgs_bbar',ybbarh,dsig)
      call filld('dphi_Higgs_bbar',phibbarh,dsig)
      call filld('Reta_Higgs_bbar',retabbarh,dsig)
      call filld('Ry_Higgs_bbar',rbbarh,dsig)

c Higgs-b-bbar system:
      call filld('pt_Higgsbbbar',ptbbh,dsig)
      call filld('ptzoom_Higgsbbbar',ptbbh,dsig)
      call filld('mass_Higgsbbbar',mbbh,dsig)

c bottom and anti-bottom:
c pt 
      call filld('pt_b',ptb,dsig)
      call filld('ptzoom_b',ptb,dsig)
      call filld('pt_bbar',ptbbar,dsig)
      call filld('ptzoom_bbar',ptbbar,dsig)

c eta:
      call filld('eta_b',etab,dsig)
      call filld('eta_bbar',etabbar,dsig)

c y:
      call filld('y_b',yb,dsig)
      call filld('y_bbar',ybbar,dsig)

c bbbar pair:
      call filld('mass_bbbar',mbbbar,dsig)
      call filld('pt_bbbar',ptbbbar,dsig)
      call filld('ptzoom_bbbar',ptbbbar,dsig)
c
c b,bbar separation:
      call filld('deta_b_bbar',etabbbar,dsig)
      call filld('dy_b_bbar',ybbbar,dsig)
      call filld('dphi_b_bbar',phibbbar,dsig)
      call filld('Reta_b_bbar',retabbbar,dsig)
      call filld('Ry_b_bbar',rbbbar,dsig)


c bbbar, Higgs separation:
      call filld('Reta_bbbar_Higgs',retabbh,dsig)
      call filld('Ry_bbbar_Higgs',rbbh,dsig)
   
   
      end

cccccccccccccccccccccccccccc

      subroutine ptyeta(p,pt,y,eta)
      implicit none
      real * 8 p(1:4),pt,y,eta
      real * 8 pp,tiny
      parameter (tiny=1d-12)
      pt=sqrt(p(1)**2+p(2)**2)
      y=log((p(4)+p(3))/(p(4)-p(3)))/2
      pp=sqrt(pt**2+p(3)**2)*(1+tiny)
      eta=log((pp+p(3))/(pp-p(3)))/2
      end
 
      subroutine getdydetadphidr(p1,p2,dy,deta,dphi,dr)
      implicit none
      include 'pwhg_math.h' 
      real * 8 p1(1:4),p2(1:4),dy,deta,dphi,dr
      real * 8 y1,eta1,pt1,mass1,phi1
      real * 8 y2,eta2,pt2,mass2,phi2
      call ptyeta(p1,pt1,y1,eta1)
      call ptyeta(p2,pt2,y2,eta2)
      dy=y1-y2
      deta=eta1-eta2
      phi1=atan2(p1(1),p1(2))
      phi2=atan2(p2(1),p2(2))
      dphi=abs(phi1-phi2)
      dphi=min(dphi,2d0*pi-dphi)
      dr=sqrt(dy**2+dphi**2)
      end

      subroutine getdydetadphidrr(p1,p2,dr)
      implicit none
      include 'pwhg_math.h' 
      real * 8 p1(1:4),p2(1:4),dy,deta,dphi,dr
      real * 8 y1,eta1,pt1,mass1,phi1
      real * 8 y2,eta2,pt2,mass2,phi2

      call ptyeta(p1,pt1,y1,eta1)
      call ptyeta(p2,pt2,y2,eta2)

      deta=eta1-eta2
      phi1=atan2(p1(1),p1(2))
      phi2=atan2(p2(1),p2(2))
      dphi=abs(phi1-phi2)
      dphi=min(dphi,2d0*pi-dphi)
      dr=sqrt(deta**2+dphi**2)
      end
