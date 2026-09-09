ccccccccccc bbH analysis with b-tagging ccccccccccc

c  The next subroutines open some histograms and prepare them 
c  to receive data 
c  You can substitute these  with your favourite ones

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'
      character * 9 pr
      integer, parameter :: ncuts=10, ncuts_CMS=2
      real *8 jetcut(ncuts),mbbhcut(ncuts)
      real *8 pt_bins(13)
      character * 10 suffix(ncuts),suffix_exp(2)
      character * 10 suffix_CMS(ncuts_CMS)
      common/jcut/jetcut,mbbhcut,suffix,suffix_CMS,suffix_exp
      common/pwhgprocess/pr
      character * 6 whcprg      
      common/cwhcprg/whcprg
      real * 8 dy,dphi,dpt,dr,dptzoom,dM,dmm,dpt3
      real*8 ptmin,ptmax,mmin,mmax
      real*8 ymin,ymax,phimin,phimax,rmin,rmax
      integer j,i
      integer,parameter :: nbins_pT_H=14, nbins_pT_bjet1_2=6,nbins_pT_bjet1=7, nbins_eta_bjet1=6,
     f     nbins_dphi_H_bjet1=7, nbins_dy_H_bjet1=6,nbins_dr_H_bjet1=14, nbins_pT_H_2=8, nbins_eta_bjet1_2=6,
     f     nbins_pT_bjet2_2=6, nbins_dr_bb_2=8, nbins_dr_Hbb_2=9,nbins_A_Hbb_2=9, nbins_m_bb_2=10, nbins_m_Hbb_2=8
      real *8  xx_pT_H(nbins_pT_H+1), xx_pt_bjet1_2(nbins_pT_bjet1_2+1),xx_pt_bjet1(nbins_pT_bjet1+1),
     f     xx_eta_bjet1(nbins_eta_bjet1+1),xx_dphi_H_bjet1(nbins_dphi_H_bjet1+1),xx_dy_H_bjet1(nbins_dy_H_bjet1+1),
     f     xx_dr_H_bjet1(nbins_dr_H_bjet1+1),xx_pt_H_2(nbins_pT_H_2+1),xx_eta_bjet1_2(nbins_eta_bjet1_2+1),
     f     xx_pt_bjet2_2(nbins_pT_bjet2_2+1),xx_dr_bb_2(nbins_dr_bb_2+1), xx_dr_Hbb_2(nbins_dr_Hbb_2+1),
     f     xx_A_Hbb_2(nbins_A_Hbb_2+1), xx_m_bb_2(nbins_m_bb_2+1),xx_m_Hbb_2(nbins_m_Hbb_2+1)
      data xx_pt_H /0d0, 20d0, 30d0, 40d0, 50d0, 65d0 , 80d0, 95d0,110d0, 130d0, 165d0, 200d0, 240d0, 300d0, 600d0/
      data xx_pt_bjet1_2 /30d0, 50d0, 70d0, 90d0, 110d0, 150d0, 200d0/
      data xx_pt_bjet1 /30d0, 50d0, 75d0, 110d0, 155d0, 220d0, 360d0,600d0/
      data xx_eta_bjet1 /0d0, 0.4d0, 0.8d0, 1.2d0, 1.6d0, 2d0, 2.4d0/
      data xx_dphi_H_bjet1 /0d0, 0.6d0, 1.2d0, 1.8d0, 2.2d0, 2.6d0,2.9d0, 3.142d0/
      data xx_dy_H_bjet1 /0d0, 0.4d0, 1d0, 1.6d0, 2.3d0, 3.4d0, 4.6d0/
      data xx_dr_H_bjet1 /0d0, 1d0, 1.6d0, 2d0, 2.4d0, 2.6d0, 2.8d0,3d0, 3.2d0, 3.4d0, 3.6d0, 4d0, 4.4d0, 5d0, 5.6d0/
      data xx_pt_H_2 /0d0, 10d0, 20d0, 40d0, 60d0, 80d0, 100d0, 150d0,200d0/
      data xx_eta_bjet1_2 /0d0, 0.4d0, 0.8d0, 1.2d0, 1.6d0, 2d0, 2.4d0/
      data xx_pt_bjet2_2 /30d0, 40d0, 60d0, 80d0, 110d0, 150d0, 200d0/
      data xx_dr_bb_2 /0.4d0, 1.2d0, 2d0, 2.4d0, 3d0, 3.4d0, 4d0, 4.5d0,5d0/
      data xx_dr_Hbb_2 /0d0, 0.6d0, 1.2d0, 1.8d0, 2.4d0, 2.7d0, 3d0,3.3d0, 4d0, 5d0/
      data xx_A_Hbb_2 /0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0,0.7d0, 0.8d0, 1d0/
      data xx_m_bb_2 /20d0, 50d0, 75d0, 100d0, 125d0, 150d0, 200d0,250d0, 300d0, 400d0, 500d0/
      data xx_m_Hbb_2 /150d0, 200d0, 250d0, 300d0, 350d0, 400d0, 450d0,500d0, 600d0/
      integer icut
      
      call inihists
         
c      pr='inc-bbH'

      dy=0.1d0
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

c Photon
      suffix_exp(1) = '-1photon'
      suffix_exp(2) = '-2photon'
c CMS-like cuts      
      suffix_CMS(1) = '-1bjet'
      suffix_CMS(2) = '-2bjet'
c pt and mbbH cuts
      suffix(1) = '-ptj20'
      jetcut(1) = 20d0
      suffix(2) = '-ptj30'
      jetcut(2) = 30d0
      suffix(3) = '-ptj60'
      jetcut(3) = 60d0
      suffix(4) = '-ptj120'
      jetcut(4) = 120d0
      suffix(5) = '-ptj160'
      jetcut(5) = 160d0
      suffix(6) = '-ptj200'
      jetcut(6) = 200d0
      suffix(7) = '-mbbH150'
      mbbhcut(7) = 150d0
      suffix(8) = '-mbbH250'
      mbbhcut(8) = 250d0
      suffix(9) = '-mbbH350'
      mbbhcut(9) = 350d0
      suffix(10) = '-mbbH450'
      mbbhcut(10) = 450d0

cccccccccccccccccccccccccc
c
c inclusive setup:
c
c incl. xsec
      call bookupeqbins('xsec',1d0,0d0,1d0)

c Higgs
      call bookupeqbins('pt_Higgs',dpt,ptmin,ptmax)
      call bookupeqbins('ptzoom_Higgs',dpt*0.1d0,ptmin,ptmax*0.1d0)
      call bookupeqbins('eta_Higgs',dy,ymin,ymax)
      call bookupeqbins('y_Higgs',dy,ymin,ymax)
      call bookupeqbins('mass_Higgs',dmm,0d0,200d0)

c Photon

      call bookupeqbins('etaphoton1',dy,-8d0,8d0)
      call bookupeqbins('etaphoton2',dy,-8d0,8d0)
      call bookupeqbins('ptphoton1',dpt,0d0,500d0)
      call bookupeqbins('ptphoton1zoom',dptzoom,0d0,50d0)
      call bookupeqbins('ptphoton2',dpt,0d0,500d0)
      call bookupeqbins('ptphoton2zoom',dptzoom,0d0,50d0)


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
c b-jets
      call bookupeqbins('n_b_jets',1d0,-0.5d0,6.5d0)
      call bookupeqbins('n_light_jets',1d0,-0.5d0,6.5d0)

c jet1
      call bookupeqbins('pt_j1',dpt,0d0,500d0)
      call bookupeqbins('pt_j1-zoom',dptzoom,0d0,50d0)

c jet2
      call bookupeqbins('pt_j2',dpt,0d0,500d0)
      call bookupeqbins('pt_j2-zoom',dptzoom,0d0,50d0)
     
c Higgs-jet1
      call bookupeqbins('dyHj1',dy,ymin,ymax)
      call bookupeqbins('detaHj1',dy,ymin,ymax)
      call bookupeqbins('dphiHj1',dphi,phimin,phimax)
      call bookupeqbins('drHj1',dr,rmin,rmax)

cccccc Cut analysis cccccc
      do icut=1,ncuts
         call bookupeqbins('total'//trim(suffix(icut)),1d0,0d0,1d0)
         call bookupeqbins('y_bbH'//trim(suffix(icut)),dy,-5d0,5d0)
         call bookupeqbins('pt_bbH'//trim(suffix(icut)),dpt,0d0,1000d0)
         call bookupeqbins('pt_bbH-zoom'//trim(suffix(icut)),dptzoom,0d0,50d0)
         call bookupeqbins('m_bbH'//trim(suffix(icut)),1d0,0d0,1000d0)
         call bookupeqbins('pt_H'//trim(suffix(icut)),dpt,0d0,1000d0)
         call bookupeqbins('pt_bb'//trim(suffix(icut)),dpt,0d0,1000d0)
      enddo
      do icut=1,6
         call bookupeqbins('y_j1'//trim(suffix(icut)),dy,-5d0,5d0)
         call bookupeqbins('dy_bbH-j1'//trim(suffix(icut)),dy,-8d0,8d0)
         call bookupeqbins('dphi_bbH-j1'//trim(suffix(icut)),dphi,0d0,pi)
      enddo

cccccc CMS analysis cccccc
      call bookupeqbins('n_b_jets-CMS',1d0,-0.5d0,6.5d0)
      call bookupeqbins('n_light_jets-CMS',1d0,-0.5d0,6.5d0)
      call bookupeqbins('total-CMS',1d0,0d0,1d0)
      do icut=1,ncuts_CMS
         call bookupeqbins('total-CMS'//trim(suffix_CMS(icut)),1d0,0d0,1d0)
         call bookup('pt_bjet1_2-CMS'//trim(suffix_CMS(icut)),nbins_pT_bjet1_2,xx_pt_bjet1_2)
      enddo
      ! distributions 1-bjet
      call bookup('pt_H-CMS'//trim(suffix_CMS(1)),nbins_pT_H,xx_pt_H)
      call bookup('pt_bjet1-CMS'//trim(suffix_CMS(1)),nbins_pT_bjet1,xx_pt_bjet1)
      call bookup('eta_bjet1-CMS'//trim(suffix_CMS(1)),nbins_eta_bjet1,xx_eta_bjet1)
      call bookup('dphi_H_bjet1-CMS'//trim(suffix_CMS(1)),nbins_dphi_H_bjet1,xx_dphi_H_bjet1)
      call bookup('dy_H_bjet1-CMS'//trim(suffix_CMS(1)),nbins_dy_H_bjet1,xx_dy_H_bjet1)
      call bookup('dR_H_bjet1-CMS'//trim(suffix_CMS(1)),nbins_dr_H_bjet1,xx_dr_H_bjet1)
      ! distributions 2-bjet
      call bookup('pt_H-CMS'//trim(suffix_CMS(2)),nbins_pT_H_2,xx_pt_H_2)
      call bookup('eta_bjet1-CMS'//trim(suffix_CMS(2)),nbins_eta_bjet1_2,xx_eta_bjet1_2)
      call bookup('pt_bjet2-CMS'//trim(suffix_CMS(2)),nbins_pT_bjet2_2,xx_pt_bjet2_2)
      call bookup('dR_bb-CMS'//trim(suffix_CMS(2)),nbins_dr_bb_2,xx_dr_bb_2)
      call bookup('dR_Hbb-CMS'//trim(suffix_CMS(2)),nbins_dr_Hbb_2,xx_dr_Hbb_2)
      call bookup('A_Hbb-CMS'//trim(suffix_CMS(2)),nbins_A_Hbb_2,xx_A_Hbb_2)
      call bookup('m_bb-CMS'//trim(suffix_CMS(2)),nbins_m_bb_2,xx_m_bb_2)
      call bookup('m_Hbb-CMS'//trim(suffix_CMS(2)),nbins_m_Hbb_2,xx_m_Hbb_2)                  

c  Photon exp cuts
      do icut=1,2
         call bookupeqbins('total'//trim(suffix_exp(icut)),1d0,0d0,1d0)

         call bookupeqbins('y_Higgs'//trim(suffix_exp(icut)),dy,-8d0,8d0)
         call bookupeqbins('pt_Higgs'//trim(suffix_exp(icut)),dpt,0d0,500d0)
         call bookupeqbins('ptHzoom'//trim(suffix_exp(icut)),dptzoom,0d0,50d0)

         call bookupeqbins('etaphoton1'//trim(suffix_exp(icut)),dy,-8d0,8d0)
         call bookupeqbins('etaphoton2'//trim(suffix_exp(icut)),dy,-8d0,8d0)
         call bookupeqbins('ptphoton1'//trim(suffix_exp(icut)),dpt,0d0,500d0)
         call bookupeqbins('ptphoton1zoom'//trim(suffix_exp(icut)),dptzoom,0d0,50d0)
         call bookupeqbins('ptphoton2'//trim(suffix_exp(icut)),dpt,0d0,500d0)
         call bookupeqbins('ptphoton2zoom'//trim(suffix_exp(icut)),dptzoom,0d0,50d0)

         call bookupeqbins('ptj1'//trim(suffix_exp(icut)),dpt,0d0,500d0)
      enddo


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
c      include 'pwhg_rad.h' ! AS

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
      real*8 retabbbar,retabh,retabbarh,retabbh
      real * 8 binsize(700)
      common/pwhghistcommon/binsize
      integer ihep,ib,ibbar,ihiggs,mu
      character * 6 whcprg      
      common/cwhcprg/whcprg
      data whcprg/'NLO   '/
      integer, parameter :: ncuts=10, ncuts_CMS=2
      real *8 jetcut(ncuts),mbbHcut(ncuts)
      real *8 pt_bins(13)
      character * 10 suffix(ncuts),suffix_exp(2)
      character * 10 suffix_CMS(ncuts_CMS)
      common/jcut/jetcut,mbbHcut,suffix,suffix_CMS,suffix_exp
C     bottom-antibottom variables
      real*8 p_bbx(4),m_bbx,pt_bbx,y_bbx,eta_bbx
c     arrays to reconstruct jets
      integer maxjet
      parameter (maxjet=2048)
      real *8 ptmin,etamax,ptsave,tmp
      real *8  ktj(maxjet),etaj(maxjet),rapj(maxjet), phij(maxjet),pj(4
     $     ,maxjet),rr,ptrel(4), yj1,etaj1,ptj1,mj1,dybbHj1,detabbHj1
     $     ,dphibbHj1,drbbHj1,dybbx,detabbx,dphibbx,drbbx
      integer j1,found,mjets,icut,i,j,ihardest,isecond,ijet
      integer maxnumbottom
      parameter (maxnumbottom=10)
      integer bottomvec(maxnumbottom),abottomvec(maxnumbottom),iabottom,ibottom,nabottom,nbottom 

      integer iib,iibbar,iihiggs
      logical write_cuts

      logical isabottom,isbottom,condition
      integer idbottom,idabottom
      parameter (idbottom=5,idabottom=-5)
      logical is_bjet_array(maxjet)
      real*8 phiggs(4),pthiggs,pH(4)
      real*8 ptb1, yb1, etab1, ptb2, yb2, etab2 
      real* 8 pt_bjet_cut_CMS,eta_bjet_cut_CMS,dy_bjet_em,deta_bjet_em,dphi_bjet_em,dr_bjet_em
     f     ,dy_bjet_ep,deta_bjet_ep,dphi_bjet_ep,dr_bjet_ep,pT_bjet1,eta_bjet1
     f     ,dy_H_bjet1,deta_H_bjet1,dphi_H_bjet1,dr_H_bjet1,pT_bjet2,dy_H_bjet2,deta_H_bjet2,dphi_H_bjet2
     f     ,dr_H_bjet2,dr_Hb1b2_min,dr_Hb1b2_max,A_Hbb,dy_b1_b2,deta_b1_b2,dphi_b1_b2,dr_b1_b2,y_b1b2
     f     ,eta_b1b2,pt_b1b2,m_b1b2,y_b1b2H,eta_b1b2H,pt_b1b2H,m_b1b2H,detahj1
     f     ,dphihj1,drhj1,dyhj1
      integer total_b_jets,counter_b_jets,counter_light_jets,counter_all_jets
      real *8 p_b1b2(4), p_b1b2h(4)
c      integer  ihardest,isecond,ijet
      real *8  kt_light_jet(maxjet),eta_light_jet(maxjet),rap_light_jet(maxjet), phi_light_jet(maxjet)
      real *8 p_ep(4),y_ep,eta_ep,pt_ep,m_ep
      real *8 p_em(4),y_em,eta_em,pt_em,m_em
c      real *8  ktj(maxjet),etaj(maxjet),rapj(maxjet), phij(maxjet),pj(4
c     $     ,maxjet),rr,ptrel(4), yj1,etaj1,ptj1,mj1,dybbx,detabbx,dphibbx,drbbx
 
      integer n_bjets,i_bjet1,i_bjet2
 
      real *8 p_b_jets(4,maxjet),p_light_jets(4,maxjet),p_all_jets(4,maxjet)
      real *8  kt_b_jet(maxjet),eta_b_jet(maxjet),rap_b_jet(maxjet), phi_b_jet(maxjet)
      real *8 y_photon_a,eta_photon_a,pt_photon_a,m_photon_a,
     1     y_photon_b,eta_photon_b,pt_photon_b,m_photon_b,
     2     y_photon1,eta_photon1,pt_photon1,m_photon1,
     3     y_photon2,eta_photon2,pt_photon2,m_photon2
      integer n_photon,maxnumphoton
      parameter (maxnumphoton=10)
      integer photonvec(maxnumphoton)
      real *8 BRHaa

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

      BRHaa=0.0027
      
      if(write_cuts) then
         write(*,*) '********************************************'
         write(*,*) '********************************************'
         write(*,*) '                ANALYSIS CUTS               '
         write(*,*) '********************************************'
         write(*,*) '********************************************'
         write(*,*) ''
         write(*,*) 'INCLUSIVE ANALYSIS:'
         write(*,*) 'Higgs decays into two photons with'
         write(*,*) 'BRHaa= ', BRHaa
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
      n_photon=0
      
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

         elseif (WHCPRG.eq.'PY8') then
c       Loop over final state particles to find Higgs, bottom quark, photon
         do ihep=1,nhep
            if((idhep(ihep).eq.25).and.(isthep(ihep).eq.1)) then !'stable' Higgs
                  ihiggs=ihep
                  pH(1:4) = phep(1:4,ihiggs)
                  iihiggs=iihiggs+1
            endif  
           if(idhep(ihep).eq.22) then
               n_photon=n_photon+1
               photonvec(n_photon)=ihep
            endif
c             print*,'n_photon--------',n_photon
c             print*,'iihiggs--------',iihiggs
c                elseif(isthep(ihep).eq.2) then ! particle has decayed
c                  ihiggs=ihep
c                  iihiggs=iihiggs+1
c               endif
c            endif !idhep
            if (idhep(ihep).eq.5) then !bottom
               ib = ihep
               iib = iib+1
            endif   
            if (idhep(ihep).eq.-5) then !bbar
               ibbar = ihep
               iibbar = iibbar+1
            endif   
         enddo !nhep
      else
         do ihep=1,nhep
            if (((isthep(ihep).eq.1).or.(isthep(ihep).eq.2)
     #.or.(isthep(ihep).eq.155).or.(isthep(ihep).eq.195))
     #.and.(idhep(ihep).eq.25)) then
               ihiggs=ihep
               iihiggs=iihiggs+1
               pH(1:4) = phep(1:4,ihiggs)
            endif
         enddo


         endif !pythia or herwig

c check that at least one bottom and one anti-bottom are selected
c during showering there may be more than one entry with PDG=5,
c with different IST. 
         if (iib.lt.1) then
            write(*,*) "Error in pwhg_analysis: ",iib," bottoms"
c            call printleshouches
            call exit(1)
         endif
         if (iibbar.lt.1) then
            write(*,*) "Error in pwhg_analysis: ",iibbar," anti-bottoms"
c            call printleshouches
            call exit(1)
         endif

      endif                     !parton or hadron level analysis

      if(rwl_num_weights.eq.0) then
         if(abs(dsig0)>1d2 .or. dsig0+1 .eq. dsig0) then
             write(*,*) "LARGE weight. DISCARDING EVENT, weight = ",dsig0
             return
         endif
      else
         do i=1,rwl_num_weights
           if(abs(dsig(i))>1d2 .or. dsig(i)+1 .eq. dsig(i)) then
             write(*,*) "LARGE weight. DISCARDING EVENT, i, weight = ",i, dsig(i)
             return
           endif
         enddo
      endif

c     Find the number of bottoms
      nbottom=0
      nabottom=0
      bottomvec=0
      abottomvec=0
      do ihep=1,nhep
         isbottom=idhep(ihep).eq.idbottom
         isabottom=idhep(ihep).eq.idabottom
         if (isthep(ihep).eq.1 .or. isthep(ihep).eq.62) then
            if(isbottom) then
               nbottom=nbottom+1
               bottomvec(nbottom)=ihep
            elseif(isabottom) then
               nabottom=nabottom+1
               abottomvec(nabottom)=ihep
            endif
         endif
      enddo
c     Check number of bottoms
      if (WHCPRG.eq.'PY8') then                                                                                                                 
         if (nbottom.lt.1.or.nabottom.lt.1) then                                                                                                   
            write(*,*) "Not enough bottom quarks found. PROGRAM ABORT ",nabottom,nbottom                                                           
            call exit(1)                                                                                                                           
         else                                                                                                                                      
            ibottom=bottomvec(1)                                                                                                                   
            iabottom=abottomvec(1)                                                                                                                 
            do i=1,nbottom                                                                                                                         
               if(phep(1,i)**2+phep(2,i)**2 .gt. phep(1,ibottom)**2+phep(2,ibottom)**2) ibottom=bottomvec(i)                                       
            enddo                                                                                                                                  
            do i=1,nabottom                                                                                                                        
               if(phep(1,i)**2+phep(2,i)**2 .gt. phep(1,iabottom)**2+phep(2,iabottom)**2) iabottom=bottomvec(i)                                    
            enddo                                                                                                                                  
         endif                                                                                                                                     
      else                                                                                                                                         
         if (nbottom.ne.1.or.nabottom.ne.1) then                                                                                                   
            write(*,*) "Too many bottom quarks found. PROGRAM ABORT ",nabottom,nbottom                                                             
            call exit(1)                                                                                                                           
         else                                                                                                                                      
            ibottom=bottomvec(1)                                                                                                                   
            iabottom=abottomvec(1)                                                                                                                 
         endif                                                                                                                                     
      endif


c bottoms:
      call ptyeta(phep(1,ib),ptb,yb,etab)
      call ptyeta(phep(1,ibbar),ptbbar,ybbar,etabbar)

      call getdydetadphidr(phep(1,ib),phep(1,ibbar),
     %     ybbbar,etabbbar,phibbbar,rbbbar)

      call getdydetadphidrr(phep(1,ib),phep(1,ibbar),
     %     retabbbar)

      if(ptb.gt.ptbbar)then                                                                                                                        
        ptb1  = ptb                                                                                                                              
        yb1   = yb                                                                                                                               
        etab1 = etab                                                                                                                             
        ptb2  = ptbbar                                                                                                                             
        yb2   = ybbar                                                                                                                              
        etab2 = etabbar                                                                                                                            
      else                                                                                                                                         
        ptb1  = ptbbar                                                                                                                             
        yb1   = ybbar                                                                                                                              
        etab1 = etabbar                                                                                                                            
        ptb2  = ptb                                                                                                                              
        yb2   = yb                                                                                                                               
        etab2 = etab                                                                                                                             
      endif

c AS
      if(IIhiggs.gt.0) then
         write(*,*) 'CB: Stable Higgs in the event.'
         write(*,*) 'Why the Higgs did not decay into two photons?'
         print*, 'IIhiggs= ', IIhiggs
         call exit(-1)
      else
         dsig=dsig*BRHaa
      endif
      
      if(IIhiggs.eq.0) then     ! no Higgs found, so reconstruct from photons (if found)
         if(n_photon.lt.2) then
            write(*,*) 'ERROR: Higgs not found and less than two photons found to reconstuct Higgs'
            call exit(1)
         else if(n_photon.gt.2) then
            write(*,*) 'ERROR: Higgs not found and more than two photons found to reconstuct Higgs (needs to be adapted when QED shower or MPI on)'
            call exit(1)
         else
            call getyetaptmass2(phep(:,photonvec(1)),y_photon_a,eta_photon_a,pt_photon_a,m_photon_a)
            call getyetaptmass2(phep(:,photonvec(2)),y_photon_b,eta_photon_b,pt_photon_b,m_photon_b)

            if(pt_photon_a .gt. pt_photon_b) then
               y_photon1   = y_photon_a
               eta_photon1 = eta_photon_a
               pt_photon1  = pt_photon_a
               m_photon1   = m_photon_a

               y_photon2   = y_photon_b
               eta_photon2 = eta_photon_b
               pt_photon2  = pt_photon_b
               m_photon2   = m_photon_b
            else
               y_photon1   = y_photon_b
               eta_photon1 = eta_photon_b
               pt_photon1  = pt_photon_b
               m_photon1   = m_photon_b

               y_photon2   = y_photon_a
               eta_photon2 = eta_photon_a
               pt_photon2  = pt_photon_a
               m_photon2   = m_photon_a
            endif

c  reconstruct Higgs momentum from hardest and second hardest photon (only two photons should be there)
            do mu=1,4
               pH(mu)=phep(mu,photonvec(1))+phep(mu,photonvec(2))
            enddo
            iihiggs = iihiggs + 1
         endif
      endif
c              print*,'iihiggs',iihiggs
c              print*,'pH',pH

      if(iihiggs.lt.1) then
         write(*,*) 'ERROR: Higgs not found'
         call exit(1)
      elseif(iihiggs.gt.1) then
         write(*,*) 'ERROR: more Higgs-like particles found'
         call exit(1)
      endif


      
c Higgs:

cc      phiggs(:)=phep(:,ihiggs)
c      call ptyeta(phep(1,ihiggs),pth,yh,etah)

c      mhiggs = dsqrt(abs(phep(4,ihiggs)**2-phep(1,ihiggs)**2-
c     &                   phep(2,ihiggs)**2-phep(3,ihiggs)**2))

      call getyetaptmass2(pH,yh,etah,pth,mHiggs)
 

      call getdydetadphidr(phep(1,ib),pH,
     %     ybh,etabh,phibh,rbh)
c      call getdydetadphidr(phep(1,ib),phep(1,ihiggs),
c     %     ybh,etabh,phibh,rbh)

       call getdydetadphidr(phep(1,ibbar),pH,
     %     ybbarh,etabbarh,phibbarh,rbbarh)
c       call getdydetadphidr(phep(1,ibbar),phep(1,ihiggs),
c     %     ybbarh,etabbarh,phibbarh,rbbarh)

       call getdydetadphidrr(phep(1,ib),pH,
     %     retabh)
 
c      call getdydetadphidrr(phep(1,ib),phep(1,ihiggs),
c    %     retabh)

       call getdydetadphidrr(phep(1,ibbar),pH,
     %     retabbarh)

c      call getdydetadphidrr(phep(1,ibbar),phep(1,ihiggs),
c     %     retabbarh)

c mass of the pair
      do mu=1,4
         ppairbbbar(mu)=phep(mu,ib)+phep(mu,ibbar)
      enddo
      mbbbar=dsqrt(abs(ppairbbbar(4)**2
     1            -ppairbbbar(1)**2-ppairbbbar(2)**2-ppairbbbar(3)**2))
      ptbbbar=dsqrt(abs(ppairbbbar(1)**2+ppairbbbar(2)**2))

c bbbar+Higgs system:
c      pbbh(1:4) = ppairbbbar(1:4)+phep(1:4,ihiggs)
      do mu=1,4
      pbbh(mu) = ppairbbbar(mu)+pH(mu)
      enddo
      ptbbh = dsqrt(abs(pbbh(1)**2+pbbh(2)**2))
      mbbh = dsqrt(abs(pbbh(4)**2-pbbh(1)**2-pbbh(2)**2-pbbh(3)**2))

c bbbar vs Higgs:

c      call getdydetadphidr(ppairbbbar(1),phep(1,ihiggs),
c     %     ybbh,etabbh,phibbh,rbbh)

c      call getdydetadphidrr(ppairbbbar(1),phep(1,ihiggs),
c     %     retabbh)


      call getdydetadphidr(ppairbbbar(1),pH(1),
     %     ybbh,etabbh,phibbh,rbbh)

      call getdydetadphidrr(ppairbbbar(1),pH(1),
     %     retabbh)
cccccccccccccccccccccccccccccc
c
c inclusive analysis:
c
c xsec:
      call filld('xsec',0.5d0,dsig)

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

c transverse invariant mass cuts
      do icut= 7, ncuts
         if(mbbh .lt. mbbHcut(icut)) then
            call filld('total'//trim(suffix(icut)),0.5d0,dsig)
            call filld('y_bbH'//trim(suffix(icut)),ybbh,dsig)
            call filld('pt_bbH'//trim(suffix(icut)),ptbbh,dsig)
            call filld('pt_bbH-zoom'//trim(suffix(icut)),ptbbh,dsig)
            call filld('m_bbH'//trim(suffix(icut)),mbbh,dsig)
            call filld('pt_H'//trim(suffix(icut)),pth,dsig)
            call filld('pt_bb'//trim(suffix(icut)),ptbbbar,dsig)
         endif
      enddo

c     photon histograms
      call filld('etaphoton1',eta_photon1,dsig)
      call filld('etaphoton2',eta_photon2,dsig)

      call filld('ptphoton1',pt_photon1,dsig)
      call filld('ptphoton1zoom',pt_photon1,dsig)
      call filld('ptphoton2',pt_photon2,dsig)
      call filld('ptphoton2zoom',pt_photon2,dsig)



cccccccccccccccccccc b-JET ANALYSIS ccccccccccccccccccccc
c     find jets
      rr=0.4d0
      ptmin=0d0
      etamax=4.5d0

      call buildjets(1,rr,ptmin,mjets,ktj,etaj,rapj,phij,ptrel,pj,is_bjet_array,total_b_jets)

      if(mjets.eq.0) return

c     now I loop over the jets and check                                                                                                           
c     that they fulfill the eta requirement                                                                                                        
c     I group them into light jets and b jets                                                                                                      
      counter_b_jets=0                                                                                                                             
      counter_light_jets=0                                                                                                                         
      counter_all_jets=0                                                                                                                           
      do i=1,mjets                                                                                                                                 
        if(dabs(etaj(i)).lt.etamax) then                                                                                                           
c     increase the counter of the corresponding jet                                                                                                
          if(is_bjet_array(i)) then                                                                                                                
            counter_b_jets=counter_b_jets+1                                                                                                        
            p_b_jets(:,counter_b_jets) = pj(:,i)                                                                                                   
          else                                                                                                                                     
            counter_light_jets=counter_light_jets+1                                                                                                
            p_light_jets(:,counter_light_jets) = pj(:,i)                                                                                           
          endif                                                                                                                                    
          counter_all_jets=counter_all_jets+1                                                                                                      
          p_all_jets(:,counter_all_jets) = pj(:,i)                                                                                                 
        endif                                                                                                                                      
      enddo

      call filld('n_b_jets',1d0*counter_b_jets,dsig)
      call filld('n_light_jets',1d0*counter_light_jets,dsig)

      !     compute useful kinematical quantities for each of these categories                                                                           
      do j=1,counter_light_jets                                                                                                                    
         call getyetaptmass2(p_light_jets(:,j),rap_light_jet(j),eta_light_jet(j),kt_light_jet(j),tmp)                                              
         phi_light_jet(j)=atan2(p_light_jets(2,j),p_light_jets(1,j))                                                                               
      enddo                                                                                                                                        
                                                                                                                                                   
      do j=1,counter_b_jets                                                                                                                        
         call getyetaptmass2(p_b_jets(:,j),rap_b_jet(j),eta_b_jet(j),kt_b_jet(j),tmp)                                                              
         phi_b_jet(j)=atan2(p_b_jets(2,j),p_b_jets(1,j))                                                                           
      enddo                                                                                                                                   

ccc JET1
      ihardest=-1
      if(counter_light_jets.ge.1) then
c     if at least a jet was found, loop over jets and find the hardest within rapidity cut
         ptsave=-1d0
         do ijet=1,counter_light_jets
            condition=(dabs(eta_light_jet(ijet)).le.etamax).and.(kt_light_jet(ijet).ge.ptsave)
            if(condition) then
               ihardest=ijet
               ptsave=kt_light_jet(ihardest)
            endif
            if(kt_light_jet(ijet).lt.ptmin) then
               write(*,*) 'ERROR1: this cannot happen'
            endif
         enddo
c     if a good jet is found, fill histograms
         if(ihardest.ne.-1) then
c    OBS for j1
            call filld('pt_j1',kt_light_jet(ihardest),dsig)
            call filld('pt_j1-zoom',kt_light_jet(ihardest),dsig)

c    OBS for Higgs-j1
c            call getdydetadphidr2(phep(1,ihiggs),p_light_jets(:,ihardest),dyHj1,detaHj1,dphiHj1,drHj1)
            call getdydetadphidr2(pH,p_light_jets(:,ihardest),dyHj1,detaHj1,dphiHj1,drHj1)
            call filld('dyHj1',dyHj1,dsig)
            call filld('detaHj1',detaHj1,dsig)
            call filld('dphiHj1',dphiHj1,dsig)
            call filld('drHj1',drHj1,dsig)
            
            call getdydetadphidr2(pbbh,p_light_jets(:,ihardest),dybbHj1,detabbHj1,dphibbHj1,drbbHj1) 

            do icut = 1,6
               if(kt_light_jet(ihardest) .gt. jetcut(icut)) then
                  call filld('total'//trim(suffix(icut)),0.5d0,dsig)
                  call filld('y_bbH'//trim(suffix(icut)),ybbh,dsig)
                  call filld('pt_bbH'//trim(suffix(icut)),ptbbh,dsig)
                  call filld('pt_bbH-zoom'//trim(suffix(icut)),ptbbh,dsig)
                  call filld('m_bbH'//trim(suffix(icut)),mbbh,dsig)
                  call filld('pt_H'//trim(suffix(icut)),pth,dsig)
                  call filld('pt_bb'//trim(suffix(icut)),ptbbbar,dsig)
                  call filld('y_j1'//trim(suffix(icut)),rapj(ihardest),dsig)
                  call filld('dy_bbH-j1'//trim(suffix(icut)),dybbHj1,dsig)
                  call filld('dphi_bbH-j1'//trim(suffix(icut)),dphibbHj1,dsig) 
               endif
            enddo
            
         endif
      endif

ccc JET2
      isecond=-1      
      if(counter_light_jets.ge.2) then
c     if at least 2 jets found, loop over jets, find the second-hardest within rapidity cut
         ptsave=-1d0
         do ijet=1,counter_light_jets
            condition=(dabs(eta_light_jet(ijet)).le.etamax).and.(kt_light_jet(ijet).ge.ptsave).and.(ijet.ne.ihardest)
            if(condition) then
               isecond=ijet
               ptsave=kt_light_jet(isecond)
            endif
            if(kt_light_jet(ijet).lt.ptmin) then
               write(*,*) 'ERROR2: this cannot happen'
            endif
         enddo
c     if a good jet is found, fill histograms
         if(isecond.ne.-1) then
c         print*, '2 jets found'
            call filld('pt_j2',kt_light_jet(isecond),dsig)
            call filld('pt_j2-zoom',kt_light_jet(isecond),dsig)
         endif
      endif

      if((isecond.eq.ihardest).and.(isecond.ne.-1)) then
         write(*,*) 'ERROR3: this cannot happen'
      endif

cccc     CMS-like analysis    cccc

      pt_bjet_cut_CMS = 30d0
      eta_bjet_cut_CMS = 2.4d0

c Search for bjet satisfying the cuts
      pT_bjet1 = 0d0
      pT_bjet2 = 0d0
      i_bjet1 = -10000
      i_bjet2 = -10000
      n_bjets = 0
      do j=1,counter_b_jets
         call getdydetadphidr2(p_b_jets(:,j),p_em,dy_bjet_em,deta_bjet_em,dphi_bjet_em,dr_bjet_em)
         call getdydetadphidr2(p_b_jets(:,j),p_ep,dy_bjet_ep,deta_bjet_ep,dphi_bjet_ep,dr_bjet_ep)
         if(kt_b_jet(j) > pt_bjet_cut_CMS .and. abs(eta_b_jet(j)) < eta_bjet_cut_CMS ) then
            n_bjets=n_bjets+1
            if(pT_bjet1 < kt_b_jet(j)) then
               i_bjet2 = i_bjet1
               i_bjet1 = j
               pT_bjet1 = kt_b_jet(i_bjet1)
             if(i_bjet2.gt.0) pT_bjet2 = kt_b_jet(i_bjet2)
            else
               if(i_bjet2.lt.0) i_bjet2 = j
            endif

         endif
      enddo
      if(n_bjets.ge.1) then
         pT_bjet1 = kt_b_jet(i_bjet1)
         eta_bjet1 = eta_b_jet(i_bjet1)
c         call getdydetadphidr2(phiggs,p_b_jets(:,i_bjet1),dy_H_bjet1,deta_H_bjet1,dphi_H_bjet1,dr_H_bjet1)
         call getdydetadphidr2(pH,p_b_jets(:,i_bjet1),dy_H_bjet1,deta_H_bjet1,dphi_H_bjet1,dr_H_bjet1)
      endif
      if(n_bjets.ge.2) then
         pT_bjet2 = kt_b_jet(i_bjet2)
c         call getdydetadphidr2(phiggs,p_b_jets(:,i_bjet2),dy_H_bjet2,deta_H_bjet2,dphi_H_bjet2,dr_H_bjet2)
         call getdydetadphidr2(pH,p_b_jets(:,i_bjet2),dy_H_bjet2,deta_H_bjet2,dphi_H_bjet2,dr_H_bjet2)
         dr_Hb1b2_min = min(dr_H_bjet1,dr_H_bjet2)
         dr_Hb1b2_max = max(dr_H_bjet1,dr_H_bjet2)
         A_Hbb = (dr_Hb1b2_max - dr_Hb1b2_min) / (dr_Hb1b2_max + dr_Hb1b2_min)
         
         call getdydetadphidr2(p_b_jets(:,i_bjet1),p_b_jets(:,i_bjet2),dy_b1_b2,deta_b1_b2,dphi_b1_b2,dr_b1_b2)
c        bjet1bjet2 momentum
         do mu=1,4
            p_b1b2(mu)=p_b_jets(mu,i_bjet1)+p_b_jets(mu,i_bjet2)
         enddo
         call getyetaptmass2(p_b1b2,y_b1b2,eta_b1b2,pt_b1b2,m_b1b2)

c        bjet1bjet2H momentum
         do mu=1,4
c            p_b1b2H(mu)=p_b1b2(mu)+phiggs(mu)
            p_b1b2H(mu)=p_b1b2(mu)+pH(mu)
         enddo
         call getyetaptmass2(p_b1b2H,y_b1b2H,eta_b1b2H,pt_b1b2H,m_b1b2H)         
      endif

      call filld('total-CMS',0.5d0,dsig)
      call filld('n_b_jets-CMS',1d0*n_bjets,dsig)
      call filld('n_light_jets-CMS',1d0*counter_light_jets,dsig)

      do icut=1,ncuts_CMS
             if(n_bjets .ge. icut) then
                call filld('total-CMS'//trim(suffix_CMS(icut)),0.5d0,dsig)
                call filld('pt_bjet1_2-CMS'//trim(suffix_CMS(icut)),pT_bjet1,dsig)
             endif
      enddo

      if(n_bjets .ge. 1) then
            ! distributions with at least 1-bjet
            call filld('pt_H-CMS'//trim(suffix_CMS(1)), pth, dsig)
            call filld('pt_bjet1-CMS'//trim(suffix_CMS(1)), pt_bjet1, dsig)
            call filld('eta_bjet1-CMS'//trim(suffix_CMS(1)), abs(eta_bjet1), dsig)
            call filld('dphi_H_bjet1-CMS'//trim(suffix_CMS(1)), dphi_H_bjet1, dsig)
            call filld('dy_H_bjet1-CMS'//trim(suffix_CMS(1)), dy_H_bjet1, dsig)
            call filld('dR_H_bjet1-CMS'//trim(suffix_CMS(1)), dr_H_bjet1, dsig)
         endif
         if(n_bjets .ge. 2) then
            ! distributions with at least 2-bjet
            call filld('pt_H-CMS'//trim(suffix_CMS(2)), pth , dsig)
            call filld('eta_bjet1-CMS'//trim(suffix_CMS(2)), abs(eta_bjet1), dsig)
            call filld('pt_bjet2-CMS'//trim(suffix_CMS(2)), pt_bjet2, dsig)
            call filld('dR_bb-CMS'//trim(suffix_CMS(2)), dr_b1_b2, dsig)
            call filld('dR_Hbb-CMS'//trim(suffix_CMS(2)), dr_Hb1b2_min, dsig)
            call filld('A_Hbb-CMS'//trim(suffix_CMS(2)), A_Hbb, dsig)
            call filld('m_bb-CMS'//trim(suffix_CMS(2)), m_b1b2, dsig)
            call filld('m_Hbb-CMS'//trim(suffix_CMS(2)), m_b1b2H, dsig)
         endif

c     exp cuts
      do icut=1,2
         if( (abs(eta_photon1) .le. 1.37 .or. (abs(eta_photon1) .le. 2.37 .and. abs(eta_photon1) .ge. 1.52) )
     $ .and. (abs(eta_photon2) .le. 1.37 .or. (abs(eta_photon2) .le. 2.37 .and. abs(eta_photon2) .ge. 1.52) ) ) then
            if( (icut.eq.1 .and. pt_photon1 .ge. 0.35*mHiggs .and. pt_photon2.ge. 0.25*mHiggs)
     $     .or. (icut.eq.2 .and. dsqrt(pt_photon1*pt_photon2) .ge.0.35*mHiggs.and. pt_photon2 .ge. 0.25*mHiggs)) then 
               call filld('total'//trim(suffix_exp(icut)),0.5d0,dsig)
               call filld('y_Higgs'//trim(suffix_exp(icut)),yH,dsig)
               call filld('etaphoton1'//trim(suffix_exp(icut)),eta_photon1,dsig)
               call filld('etaphoton2'//trim(suffix_exp(icut)),eta_photon1,dsig)

               call filld('pt_Higgs'//trim(suffix_exp(icut)),ptH,dsig)
               call filld('ptHzoom'//trim(suffix_exp(icut)),ptH,dsig)
               call filld('ptphoton1'//trim(suffix_exp(icut)),pt_photon1,dsig)
               call filld('ptphoton1zoom'//trim(suffix_exp(icut)),pt_photon1,dsig)
               call filld('ptphoton2'//trim(suffix_exp(icut)),pt_photon2,dsig)
               call filld('ptphoton2zoom'//trim(suffix_exp(icut)),pt_photon2,dsig)

               call filld('ptj1'//trim(suffix_exp(icut)),ktj(1),dsig)
            endif
         endif
      enddo
      
cccccccccccccccccccc end of b_JET ANALYSIS ccccccccccccccccccccc
      
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

c CB: due to multiple definitions, I prefer to change the name of the routine (see below)      
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

c If you want only dr use this get...drr
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
c     If you want also dy, deta, dphi use this get..dr2
      subroutine getdydetadphidr2(p1,p2,dy,deta,dphi,dr)
      implicit none
      include 'pwhg_math.h' 
      real * 8 p1(*),p2(*),dy,deta,dphi,dr
      real * 8 y1,eta1,pt1,mass1,phi1
      real * 8 y2,eta2,pt2,mass2,phi2
      call getyetaptmass2(p1,y1,eta1,pt1,mass1)
      call getyetaptmass2(p2,y2,eta2,pt2,mass2)
      dy=y1-y2
      deta=eta1-eta2
      phi1=atan2(p1(1),p1(2))
      phi2=atan2(p2(1),p2(2))
      dphi=abs(phi1-phi2)
      dphi=min(dphi,2d0*pi-dphi)
      dr=sqrt(deta**2+dphi**2)
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
               ! JM: changed isthep == 1 by is_final_state, check! (see also further changes below)
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

c JM: added      
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

c      write(*,*)"nbjet = ",nbjet
c     jets are ordered in decreasing pt. Set up array of info on jets
c     if jetinfo=0 then non-b jet
c     if jetinfo=5 then b jet
c     if jetinfo=-5 then bbar jet
c      do i=1,njets
c         jetinfo(i)=0
c      enddo

      total_b_jets = 0
      do i=1,njets
         jetinfo(i)=0
         is_bjet_array(i)=.false.
         do j=1,nbjet
            if (i.eq.nbjet_array(j)) then
               !jetinfo(i)=5
               is_bjet_array(i)=.true.
            endif
         enddo
         do j=1,nbbarjet
            if (i.eq.nbbarjet_array(j)) then
               !jetinfo(i)=-5
               is_bjet_array(i)=.true.
             endif
          enddo
          if(is_bjet_array(i)) total_b_jets = total_b_jets + 1 
       enddo

c       write(*,*)"total_b_jets = ",total_b_jets 

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

