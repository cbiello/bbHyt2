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
      real*8 ptzoommin,ptzoommax
      integer j,i

      call inihists
      
      dy=0.2d0
      dphi=0.4d0
      dmm = 5d0

      dpt=20d0
      ptmin = 0d0
      ptmax = 340d0
      
      dptzoom=2d0
      ptzoommin= 0d0
      ptzoommax= 68d0
      
      dr = 0.4d0
      ymin = -2.6d0
      ymax = 2.6d0

      phimin = 0d0
      phimax = 3.2d0

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
c
c incl. Higgs
      call bookupeqbins('pt_Higgs',dpt,ptmin,ptmax)

cccccccccccccccc
      call bookupeqbins('ptzoom_Higgs_bin1',4.0d0,0.0d0,64d0)
      call bookupeqbins('ptzoom_Higgs_bin2',8.0d0,0.0d0,64d0)
      call bookupeqbins('absy_Higgs_bin1',0.2d0,0.0d0,2.4d0)
      call bookupeqbins('absy_Higgs_bin2',0.4d0,0.0d0,2.4d0)
cccccccccccccccc

      call bookupeqbins('eta_Higgs',dy,ymin,ymax)
      call bookupeqbins('y_Higgs-ptHge10',dy,ymin,ymax)
      call bookupeqbins('y_Higgs-ptHge15',dy,ymin,ymax)
      call bookupeqbins('y_Higgs-ptHge20',dy,ymin,ymax)
c
c Higgs - analysis A1: at least an anti-kt any-flavour jet
      call bookupeqbins('pt_Higgs-A1',dpt,ptmin,ptmax)
      call bookupeqbins('ptzoom_Higgs-A1',dptzoom,ptzoommin,ptzoommax)
      call bookupeqbins('y_Higgs-A1',dy,ymin,ymax)
c jet - analysis A1
      call bookupeqbins('pt_jet-A1',dpt,ptmin,ptmax)
      call bookupeqbins('y_jet-A1',dy,ymin,ymax)
      call bookupeqbins('dy_Higgs_jet-A1',dy,ymin,ymax)
      call bookupeqbins('dphi_Higgs_jet-A1',dphi,phimin,phimax)
      call bookupeqbins('Ry_Higgs_jet-A1',dr,rmin,rmax)
c
c Higgs - analysis B1: at least an anti-kt naive-b jet
      call bookupeqbins('pt_Higgs-B1',dpt,ptmin,ptmax)
      call bookupeqbins('ptzoom_Higgs-B1',dptzoom,ptzoommin,ptzoommax)
      call bookupeqbins('y_Higgs-B1',dy,ymin,ymax)
c jet - analysis B1
      call bookupeqbins('pt_jet-B1',dpt,ptmin,ptmax)
      call bookupeqbins('y_jet-B1',dy,ymin,ymax)
      call bookupeqbins('dy_Higgs_jet-B1',dy,ymin,ymax)
      call bookupeqbins('dphi_Higgs_jet-B1',dphi,phimin,phimax)
      call bookupeqbins('Ry_Higgs_jet-B1',dr,rmin,rmax)
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

      integer maxjet
      parameter (maxjet=2048)
      integer idbottom,idabottom
      parameter (idbottom=5,idabottom=-5)
      logical is_bjet_array(maxjet)
      real *8 p_b_jets(4,maxjet),p_all_jets(4,maxjet)
      real *8  kt_b_jet(maxjet),eta_b_jet(maxjet),rap_b_jet(maxjet), phi_b_jet(maxjet)
      real *8  kt_all_jet(maxjet),eta_all_jet(maxjet),rap_all_jet(maxjet), phi_all_jet(maxjet)
      logical condition
      integer counter_all_jets, counter_b_jets, ihardest, ijet, j
      integer mjets, total_b_jets
      real *8  ktj(maxjet),etaj(maxjet),rapj(maxjet), phij(maxjet),pj(4
     $     ,maxjet),rr,ptrel(4)
      real *8 etamax, ptsave, tmp, ptmin
      real *8 etajh,phijh,rjh,yjh


      
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

         write(*,*) 'ERROR: This analysis is for FO studies'
         write(*,*) 'Why we are not at stage 2??'
         stop  

      endif                     !parton or hadron level analysis

ccc   CB(20Dec2025): I have removed the check of large weight, since
ccc   it was based on the xsection at physics mb, i.e. we were removing
ccc   large weights of two order of magnitude bigger than 500fb, but
ccc   when I take the small mb limit, I can get large weights due to larger
ccc   xsection, therefore let me keep any large wright and remove only
ccc   the INFINITY ones.
      if(rwl_num_weights.eq.0) then
         if(dsig0+1 .eq. dsig0) then
             write(*,*) "INF weight. DISCARDING EVENT, weight = ",dsig0
             return
         endif
      else
         do i=1,rwl_num_weights
           if(abs(dsig(i))>1d3 .or. dsig(i)+1 .eq. dsig(i)) then
             write(*,*) "INF weight. DISCARDING EVENT, i, weight = ",i, dsig(i)
             return
           endif
         enddo
      endif 

cccccccc Higgs:
      call ptyeta(phep(1,ihiggs),pth,yh,etah)
      mhiggs = dsqrt(abs(phep(4,ihiggs)**2-phep(1,ihiggs)**2-
     &                   phep(2,ihiggs)**2-phep(3,ihiggs)**2))

c xsec:
      call filld('xsec',0.5d0,dsig)
      if(pth.ge.10) call filld('xsec-ptHge10',0.5d0,dsig)
      if(pth.ge.15) call filld('xsec-ptHge15',0.5d0,dsig)
      if(pth.ge.20) call filld('xsec-ptHge20',0.5d0,dsig)      
c Higgs fully-inclusive:
      call filld('pt_Higgs',pth,dsig)
      call filld('ptzoom_Higgs_bin1',pth,dsig)
      call filld('ptzoom_Higgs_bin2',pth,dsig)
      call filld('eta_Higgs',etah,dsig)
      call filld('absy_Higgs_bin1',abs(yh),dsig)
      call filld('absy_Higgs_bin2',abs(yh),dsig)
      if(pth.ge.10) call filld('y_Higgs-ptHge10',0.5d0,dsig)
      if(pth.ge.15) call filld('y_Higgs-ptHge15',0.5d0,dsig)
      if(pth.ge.20) call filld('y_Higgs-ptHge20',0.5d0,dsig)

ccccccc Build jet:
      
c     find jets
      rr=0.4d0
      ptmin=20d0
      etamax=2.4d0

c This buildjet is for NAIVE b-jets
      call buildjets(1,rr,ptmin,mjets,ktj,etaj,rapj,phij,ptrel,pj,is_bjet_array,total_b_jets)

      if(mjets.eq.0) return 

      counter_b_jets=0
      counter_all_jets=0

      do i=1,mjets

         if(dabs(etaj(i)).lt.etamax) then
            if(is_bjet_array(i)) then
               counter_b_jets=counter_b_jets+1
               do mu = 1, 4
                  p_b_jets(mu,counter_b_jets) = pj(mu,i)
               enddo
            endif
            counter_all_jets=counter_all_jets+1
            do mu = 1, 4
               p_all_jets(mu,counter_all_jets) = pj(mu,i)
            enddo
         endif

      enddo

      do j=1,counter_all_jets
         call getyetaptmass2(p_all_jets(:,j),rap_all_jet(j),eta_all_jet(j),kt_all_jet(j),tmp)
         phi_all_jet(j)=atan2(p_all_jets(2,j),p_all_jets(1,j))
      enddo

      do j=1,counter_b_jets
         call getyetaptmass2(p_b_jets(:,j),rap_b_jet(j),eta_b_jet(j),kt_b_jet(j),tmp)
         phi_b_jet(j)=atan2(p_b_jets(2,j),p_b_jets(1,j))
      enddo
      
c Analysis A1: flavour-blind, anyjet: lets choose the hardest
      ihardest=-1
      if(counter_all_jets.ge.1) then
         ptsave=-1d0
         do ijet=1,counter_all_jets
            condition=(dabs(eta_all_jet(ijet)).le.etamax).and.(kt_all_jet(ijet).ge.ptsave)
            if(condition) then
               ihardest=ijet
               ptsave=kt_all_jet(ihardest)
            endif
            if(kt_all_jet(ijet).lt.ptmin) then
               write(*,*) 'ERROR1: this cannot happen'
            endif
         enddo
         if(ihardest.ne.-1) then
cccccccccccccccccc

         call filld('pt_Higgs-A1',pth,dsig)
         call filld('ptzoom_Higgs-A1',pth,dsig)
         call filld('y_Higgs-A1',yh,dsig)
         
         call filld('pt_jet-A1',kt_all_jet(ihardest),dsig)
         call filld('y_jet-A1',rap_all_jet(ihardest),dsig)

         !getdydetadphidr returns the dR as separation in the y-phi plane
         !we should use getdydetadphidrr for Reta eventually
         call getdydetadphidr(p_all_jets(1,ihardest),phep(1,ihiggs),yjh,etajh,phijh,rjh)
         call filld('dy_Higgs_jet-A1',yjh,dsig)
         call filld('dphi_Higgs_jet-A1',phijh,dsig)
         call filld('Ry_Higgs_jet-A1',rjh,dsig)
              
ccccccccccccccccc
         endif
      endif

c Analysis B1: naive (bottoms modulo 2): lets choose the hardest bjet
      ihardest=-1
      if(counter_b_jets.ge.1) then
         ptsave=-1d0
         do ijet=1,counter_b_jets
            condition=(dabs(eta_b_jet(ijet)).le.etamax).and.(kt_b_jet(ijet).ge.ptsave)
            if(condition) then
               ihardest=ijet
               ptsave=kt_b_jet(ihardest)
            endif
            if(kt_b_jet(ijet).lt.ptmin) then
               write(*,*) 'ERROR1: this cannot happen'
            endif
         enddo
         if(ihardest.ne.-1) then
cccccccccccccccccc

         call filld('pt_Higgs-B1',pth,dsig)
         call filld('ptzoom_Higgs-B1',pth,dsig)
         call filld('y_Higgs-B1',yh,dsig)
         
         call filld('pt_jet-B1',kt_b_jet(ihardest),dsig)
         call filld('y_jet-B1',rap_b_jet(ihardest),dsig)

         call getdydetadphidr(p_b_jets(1,ihardest),phep(1,ihiggs),yjh,etajh,phijh,rjh)
         call filld('dy_Higgs_jet-B1',yjh,dsig)
         call filld('dphi_Higgs_jet-B1',phijh,dsig)
         call filld('Ry_Higgs_jet-B1',rjh,dsig)
         
ccccccccccccccccc
         endif
      endif

   
      end

cccccccccccccccccccccccccccc
ccc   General routines imported for other analyses
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

      subroutine getyetaptmass2(p,y,eta,pt,mass)
      implicit none
      real * 8 p(4),y,eta,pt,mass,pv
      real *8 tiny
      parameter (tiny=1.d-5)
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      pt=sqrt(p(1)**2+p(2)**2)
      pv=sqrt(pt**2+p(3)**2)
      if(pt.lt.tiny)then
         eta=sign(1.d0,p(3))*1.d8
      else
         eta=0.5d0*log((pv+p(3))/(pv-p(3)))
      endif
      mass=sqrt(abs(p(4)**2-pv**2))
      end
      
cccccccccccccccccccc
ccc   BUILD ANTI-KT JETS
ccc   WITH NAIVE B-JET TAGGING (NAIVE=odd number of bottoms)
cccccccccccccccccccc
      
      subroutine buildjets(iflag,rr,ptmin,mjets,kt,eta,rap,phi,
     $     ptrel,pjet,is_bjet_array,total_b_jets)
c     arrays to reconstruct jets, radius parameter rr
      implicit none
c     tell to the analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      integer iflag,mjets
      real * 8  rr,ptmin,kt(*),eta(*),rap(*),
     1     phi(*),ptrel(3),pjet(4,*)
      logical is_bjet_array(*)
      include   'hepevt.h'
      include  'LesHouches.h'
      integer   maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real * 8  ptrack(4,maxtrack),pj(4,maxjet)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      integer   ntracks,njets
      integer   j,k,mu,i
      real * 8 r,palg,tmp
      real * 8 vec(3),pjetin(0:3),pjetout(0:3),beta,
     $     ptrackin(0:3),ptrackout(0:3)
      real * 8 get_ptrel
      external get_ptrel
      logical is_B_hadron,is_BBAR_hadron,is_neutrino
      external is_B_hadron,is_BBAR_hadron,is_neutrino
      logical is_b_track(maxtrack)
      integer const_indices(maxtrack),nconst,ijet
      integer nbjet_array(maxjet),
     $     nbbarjet_array(maxjet),jetinfo(maxjet),id,nb,
     $     nbbar,nbjet,nbbarjet,total_b_jets
C - Initialize arrays and counters for output jets
      do j=1,maxtrack
         do mu=1,4
            ptrack(mu,j)=0d0
         enddo
         jetvec(j)=0
      enddo      
      ntracks=0
      do j=1,maxjet
         do mu=1,4
            pjet(mu,j)=0d0
            pj(mu,j)=0d0
         enddo
         is_bjet_array(j) = .false. 
      enddo
      if(iflag.eq.1) then
C     - Extract final state particles to feed to jet finder
         if(WHCPRG.eq.'PY8   ') then
            do j=1,nhep
c all but top/anti-top
c     exclude leptons, gauge and higgs bosons, but include gluons             
               if(isthep(j).eq.1.and.((abs(idhep(j)).le.5.or.abs(idhep(j)).ge.40
               ! changed isthep == 1 by is_final_state, check! (see also further changes below)
               !if(is_final_state(j).and.((abs(idhep(j)).le.4.or.abs(idhep(j)).ge.40
     &              .or.abs(idhep(j)).eq.21))) then
                  if(ntracks.eq.maxtrack) then
                     write(*,*) 'analyze: need to increase maxtrack!'
                     write(*,*) 'ntracks: ',ntracks
                     stop
                  endif
                  ntracks=ntracks+1
                  do mu=1,4
                     ptrack(mu,ntracks)=phep(mu,j)
                  enddo
                  itrackhep(ntracks)=j
               endif
            enddo
         else
           do j=1,nhep
c all but top/anti-top
               if(isthep(j).eq.1.and.((abs(idhep(j)).le.5.or.abs(idhep(j)).ge.40
     &              .or.abs(idhep(j)).eq.21))) then
                  if(ntracks.eq.maxtrack) then
                     write(*,*) 'analyze: need to increase maxtrack!'
                     write(*,*) 'ntracks: ',ntracks
                     stop
                  endif
                  ntracks=ntracks+1
                  do mu=1,4
                     ptrack(mu,ntracks)=phep(mu,j)
                  enddo
                  itrackhep(ntracks)=j
               endif
            enddo
         endif
      else
         do j=1,nup
            if(istup(j).eq.1.and.((abs(idhep(j)).le.5.or.abs(idhep(j)).ge.40
     &           .or.abs(idhep(j)).eq.21))) then
               if(ntracks.eq.maxtrack) then
                  write(*,*) 'analyze: need to increase maxtrack!'
                  write(*,*) 'ntracks: ',ntracks
                  stop
               endif
               ntracks=ntracks+1
               do mu=1,4
                  ptrack(mu,ntracks)=pup(mu,j)
               enddo
               itrackhep(ntracks)=j
            endif
         enddo
      endif
      if (ntracks.eq.0) then
         mjets=0
         return
      endif
C --------------------------------------------------------------------- C
C     R = 0.7   radius parameter
c palg=1 is standard kt, -1 is antikt
      palg=-1
      r=rr
c      ptmin=20d0 
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,
     $                        jetvec)
      mjets=njets

c----------------------------------------------------------------------
c     find in which ptrack the B hadrons ended up
      nbjet_array = 0
      nbbarjet_array = 0
      nbjet=0
      nbbarjet=0
c     loop over tracks
      do i=1,ntracks
         id=idhep(itrackhep(i))
         if (is_B_hadron(id)) then
            nbjet=nbjet+1
            nbjet_array(nbjet)=jetvec(i)            
         elseif (is_BBAR_hadron(id)) then   
            nbbarjet=nbbarjet+1
            nbbarjet_array(nbbarjet)=jetvec(i)                        
         endif
      enddo

c CB(30Jan26): I modify it in order to compute if there is an odd number of bs. In that case let's call it
      total_b_jets = 0
      do i=1,njets
         jetinfo(i)=0
         is_bjet_array(i)=.false.
         do j=1,nbjet
            if (i.eq.nbjet_array(j)) then
               if(is_bjet_array(i)) then
                  is_bjet_array(i)=.false.
               else
                  is_bjet_array(i)=.true.
               endif
            endif
         enddo
         do j=1,nbbarjet
            if (i.eq.nbbarjet_array(j)) then
               if(is_bjet_array(i)) then
                  is_bjet_array(i)=.false.
               else
                  is_bjet_array(i)=.true.
               endif
             endif
          enddo
          if(is_bjet_array(i)) total_b_jets = total_b_jets + 1 
       enddo

c----------------------------------------------------------------------

      if(njets.eq.0) return
c check consistency
      do k=1,ntracks
         if(jetvec(k).gt.0) then
            do mu=1,4
               pj(mu,jetvec(k))=pj(mu,jetvec(k))+ptrack(mu,k)
            enddo
         endif
      enddo
      tmp=0
      do j=1,mjets
         do mu=1,4
            tmp=tmp+abs(pj(mu,j)-pjet(mu,j))
         enddo
      enddo
      if(tmp.gt.1d-4) then
         write(*,*) ' bug!'
      endif
C --------------------------------------------------------------------- C
C - Computing arrays of useful kinematics quantities for hardest jets - C
C --------------------------------------------------------------------- C
      do j=1,mjets
         call getyetaptmass2(pjet(:,j),rap(j),eta(j),kt(j),tmp)
         phi(j)=atan2(pjet(2,j),pjet(1,j))
      enddo

c     loop over the hardest 3 jets
      do j=1,min(njets,3)
         do mu=1,3
            pjetin(mu) = pjet(mu,j)
         enddo
         pjetin(0) = pjet(4,j)         
         vec(1)=0d0
         vec(2)=0d0
         vec(3)=1d0
         beta = -pjet(3,j)/pjet(4,j)
         call mboost(1,vec,beta,pjetin,pjetout)         
c     write(*,*) pjetout
         ptrel(j) = 0
         do i=1,ntracks
            if (jetvec(i).eq.j) then
               do mu=1,3
                  ptrackin(mu) = ptrack(mu,i)
               enddo
               ptrackin(0) = ptrack(4,i)
               call mboost(1,vec,beta,ptrackin,ptrackout) 
               ptrel(j) = ptrel(j) + get_ptrel(ptrackout,pjetout)
            endif
         enddo
      enddo
      end

      function get_ptrel(pin,pjet)
      implicit none
      real * 8 get_ptrel,pin(0:3),pjet(0:3)
      real * 8 pin2,pjet2,cth2,scalprod
      pin2  = pin(1)**2 + pin(2)**2 + pin(3)**2
      pjet2 = pjet(1)**2 + pjet(2)**2 + pjet(3)**2
      scalprod = pin(1)*pjet(1) + pin(2)*pjet(2) + pin(3)*pjet(3)
      cth2 = scalprod**2/pin2/pjet2
      get_ptrel = sqrt(pin2*abs(1d0 - cth2))
      end

      function is_B_hadron(id)
      implicit none
      logical is_B_hadron
      integer id
      is_B_hadron=((id.gt.-600).and.(id.lt.-500)).or.
     $     ((id.gt.5000).and.(id.lt.6000)).or.(id.eq.5)
      end

      function is_BBAR_hadron(id)
      implicit none
      logical is_BBAR_hadron
      integer id
      is_BBAR_hadron=((id.gt.500).and.(id.lt.600)).or.
     $     ((id.gt.-6000).and.(id.lt.-5000)).or.(id.eq.-5)
      end
