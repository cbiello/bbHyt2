c     Analysis for the M1M0 study (virtinborn runs, exact vs massified).
c
c     Higgs:   pt_H, ptzoom_H, y_H, eta_H  and  pt_HbbX
c     b jets:  three flavour definitions, booked with the suffixes
c                -exp   a jet is a b jet if it holds at least one b OR bbar
c                -nai   naive/net flavour: n_b - n_bbar  /=  0
c                -ifn   interleaved flavour neutralisation (fastjet plugin)
c              and for each of them
c                pt_bj1, y_bj1, dy_H_bj1, deta_H_bj1, dphi_H_bj1,
c                dRy_H_bj1, dR_H_bj1, m_bj1bj2, m_Hbj1, pt_bj1bj2
c
c     IFN needs the IFNPlugin; if it is not linked the -ifn histograms
c     stay empty and a warning is printed once (see ifn_dummy.cc).
      subroutine init_hist
      implicit none
      include 'pwhg_bookhist-multi.h'
      integer i
      character * 4 cdef(3)
      logical dobjets
      common/cbjetdefs/cdef,dobjets
      data cdef/'-exp','-nai','-ifn'/
      real * 8 powheginput
      external powheginput
      real * 8 dpt,dy,dphi,dr,dm
      parameter (dpt=5d0, dy=0.25d0, dphi=0.1d0, dr=0.2d0, dm=10d0)

      call inihists

      call bookupeqbins('xsec',1d0,0d0,1d0)

c     Higgs
      call bookupeqbins('pt_H'     ,dpt      ,0d0  ,500d0)
      call bookupeqbins('ptzoom_H' ,dpt*0.2d0,0d0  ,100d0)
      call bookupeqbins('y_H'      ,dy       ,-5d0 ,5d0  )
      call bookupeqbins('eta_H'    ,dy       ,-5d0 ,5d0  )
c     transverse momentum of the H b bbar system
      call bookupeqbins('pt_bbH'   ,dpt      ,0d0  ,500d0)

c     b jets, one block per flavour definition.
c     bjetobs 1 switches them on; off by default, the M1M0 comparison
c     is being looked at with Higgs observables first.
      dobjets = (powheginput("#bjetobs").eq.1)
      if(.not.dobjets) return
      do i=1,3
         call bookupeqbins('pt_bj1'   //cdef(i),dpt ,0d0  ,500d0)
         call bookupeqbins('y_bj1'    //cdef(i),dy  ,-5d0 ,5d0  )
         call bookupeqbins('dy_H_bj1' //cdef(i),dy  ,-5d0 ,5d0  )
         call bookupeqbins('deta_H_bj1'//cdef(i),dy ,-5d0 ,5d0  )
         call bookupeqbins('dphi_H_bj1'//cdef(i),dphi,0d0 ,3.15d0)
         call bookupeqbins('dRy_H_bj1'//cdef(i),dr  ,0d0  ,8d0  )
         call bookupeqbins('dR_H_bj1' //cdef(i),dr  ,0d0  ,8d0  )
         call bookupeqbins('m_bj1bj2' //cdef(i),dm  ,0d0  ,1000d0)
         call bookupeqbins('m_Hbj1'   //cdef(i),dm  ,0d0  ,1000d0)
         call bookupeqbins('pt_bj1bj2'//cdef(i),dpt ,0d0  ,500d0)
      enddo

      end


      subroutine analysis(dsig0)
      implicit none
      include 'hepevt.h'
      include 'pwhg_math.h'
      include 'LesHouches.h'
      include 'pwhg_weights.h'
      include 'pwhg_bookhist-multi.h'
      real * 8 dsig0,dsig(weights_max)
      logical ini
      data ini/.true./
      save ini
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      character * 4 cdef(3)
      logical dobjets
      common/cbjetdefs/cdef,dobjets

      integer maxjet
      parameter (maxjet=128)
      real * 8 pH(0:3),pbb(0:3),psum(0:3)
      real * 8 pjet(0:3,maxjet,3)
      integer  njets(3),nbj(3)
      integer  ibj1(3),ibj2(3)
      real * 8 ptj,yj,etaj,mj
      real * 8 ptH,yH,etaH,mH
      real * 8 dy,deta,dphi,dry,drr
      real * 8 rr,ptminj
      parameter (rr=0.4d0, ptminj=25d0)
      integer i,k,mu
      logical foundH

c     multi_plot_setup (POWHEG-BOX-RES/multi_plot.f) takes three
c     arguments and zeroes dsig beyond the active weights
      call multi_plot_setup(dsig0,dsig,weights_max)
      if(all(dsig.eq.0d0)) return

      if(ini) then
         write(*,*) '*****************************************'
         if(WHCPRG.eq.'NLO   ') then
            write(*,*) '    NLO analysis (M1M0 study)         '
         elseif(WHCPRG.eq.'LHE   ') then
            write(*,*) '    LHE analysis (M1M0 study)         '
         else
            write(*,*) '    ',WHCPRG,' analysis (M1M0 study)  '
         endif
         write(*,*) '    anti-kt R = ',rr,'  pt_j > ',ptminj
         write(*,*) '*****************************************'
         ini=.false.
      endif

      call filld('xsec',0.5d0,dsig)

c     the Higgs
      call find_higgs(pH,foundH)
      if(.not.foundH) return
      call getyetaptmass(pH,yH,etaH,ptH,mH)
      call filld('pt_H'    ,ptH ,dsig)
      call filld('ptzoom_H',ptH ,dsig)
      call filld('y_H'     ,yH  ,dsig)
      call filld('eta_H'   ,etaH,dsig)

c     pt of the H b bbar system = pt of everything that is not extra
c     radiation; built from H plus all b-flavoured final-state partons
      call sum_bflavour(pbb)
      do mu=0,3
         psum(mu) = pH(mu) + pbb(mu)
      enddo
      call filld('pt_bbH',dsqrt(psum(1)**2+psum(2)**2),dsig)

c     jets, three flavour definitions
      if(.not.dobjets) return
      call build_bjets(rr,ptminj,maxjet,pjet,njets,nbj,ibj1,ibj2)

      do i=1,3
         if(nbj(i).lt.1) cycle
         call getyetaptmass(pjet(0,ibj1(i),i),yj,etaj,ptj,mj)
         call filld('pt_bj1'//cdef(i),ptj,dsig)
         call filld('y_bj1' //cdef(i),yj ,dsig)
         call getdydetadphidr(pH,pjet(0,ibj1(i),i),dy,deta,dphi,dry,drr)
         call filld('dy_H_bj1'  //cdef(i),dy  ,dsig)
         call filld('deta_H_bj1'//cdef(i),deta,dsig)
         call filld('dphi_H_bj1'//cdef(i),dphi,dsig)
         call filld('dRy_H_bj1' //cdef(i),dry ,dsig)
         call filld('dR_H_bj1'  //cdef(i),drr ,dsig)
         call invmass2(pH,pjet(0,ibj1(i),i),mj)
         call filld('m_Hbj1'//cdef(i),mj,dsig)
         if(nbj(i).ge.2) then
            call invmass2(pjet(0,ibj1(i),i),pjet(0,ibj2(i),i),mj)
            call filld('m_bj1bj2'//cdef(i),mj,dsig)
            do mu=0,3
               psum(mu)=pjet(mu,ibj1(i),i)+pjet(mu,ibj2(i),i)
            enddo
            call filld('pt_bj1bj2'//cdef(i),
     $           dsqrt(psum(1)**2+psum(2)**2),dsig)
         endif
      enddo

      end


c=======================================================================
c     Final-state bookkeeping
c=======================================================================
      subroutine get_final_state(np,pp,idp)
c     final-state partons (|id| <= 5 or gluon) and their flavour
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      integer np,idp(*)
      real * 8 pp(0:3,*)
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      integer j,mu,id
      np = 0
c     NLO: the analysis driver fills HEPEVT, not the Les Houches block
      if(WHCPRG.eq.'LHE   ') then
         do j=3,nup
            if(istup(j).ne.1) cycle
            id = idup(j)
            if(id.eq.21) id = 0
            if(abs(id).gt.5) cycle
            np = np + 1
            do mu=0,3
               pp(mu,np) = pup(mu+1,j)
            enddo
            if(mu.eq.0) continue
            pp(0,np) = pup(4,j)
            pp(1,np) = pup(1,j)
            pp(2,np) = pup(2,j)
            pp(3,np) = pup(3,j)
            idp(np) = id
         enddo
      else
         do j=1,nhep
            if(isthep(j).ne.1) cycle
            id = idhep(j)
            if(id.eq.21) id = 0
            if(abs(id).gt.5) cycle
            np = np + 1
            pp(0,np) = phep(4,j)
            pp(1,np) = phep(1,j)
            pp(2,np) = phep(2,j)
            pp(3,np) = phep(3,j)
            idp(np) = id
         enddo
      endif
      end


      subroutine find_higgs(pH,found)
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      real * 8 pH(0:3)
      logical found
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      integer j
      found = .false.
c     NLO: the analysis driver fills HEPEVT, not the Les Houches block
      if(WHCPRG.eq.'LHE   ') then
         do j=3,nup
            if(idup(j).eq.25) then
               pH(0)=pup(4,j); pH(1)=pup(1,j)
               pH(2)=pup(2,j); pH(3)=pup(3,j)
               found = .true.
               return
            endif
         enddo
      else
         do j=1,nhep
            if(idhep(j).eq.25) then
               pH(0)=phep(4,j); pH(1)=phep(1,j)
               pH(2)=phep(2,j); pH(3)=phep(3,j)
               found = .true.
               return
            endif
         enddo
      endif
      end


      subroutine sum_bflavour(pbb)
c     sum of all final-state b-flavoured partons
      implicit none
      real * 8 pbb(0:3)
      integer maxp
      parameter (maxp=512)
      real * 8 pp(0:3,maxp)
      integer idp(maxp),np,j,mu
      pbb = 0d0
      call get_final_state(np,pp,idp)
      do j=1,np
         if(abs(idp(j)).eq.5) then
            do mu=0,3
               pbb(mu) = pbb(mu) + pp(mu,j)
            enddo
         endif
      enddo
      end


c=======================================================================
c     Jets with the three b-flavour definitions
c
c     EXP (1): at least one b or bbar among the constituents
c     NAI (2): net flavour  n_b - n_bbar  /= 0
c     IFN (3): interleaved flavour neutralisation, from the plugin
c=======================================================================
      subroutine build_bjets(rr,ptmin,mxj,pjet,njets,nbj,ibj1,ibj2)
      implicit none
      real * 8 rr,ptmin
      integer mxj
      real * 8 pjet(0:3,mxj,3)
      integer njets(3),nbj(3),ibj1(3),ibj2(3)
      integer maxp
      parameter (maxp=512)
      real * 8 pp(0:3,maxp),ptrack(4,maxp),pj(4,mxj)
      integer idp(maxp),np,jetvec(maxp),isb(mxj)
      integer nb(mxj),nbb(mxj)
      integer i,j,mu,nj,idef
      real * 8 palg,pt1,pt2,ptj
      logical isbjet

      njets = 0
      nbj   = 0
      ibj1  = 0
      ibj2  = 0
      call get_final_state(np,pp,idp)
      if(np.lt.1) return
      do j=1,np
         ptrack(1,j)=pp(1,j)
         ptrack(2,j)=pp(2,j)
         ptrack(3,j)=pp(3,j)
         ptrack(4,j)=pp(0,j)
      enddo

c     ---- flavour-blind anti-kt, used by EXP and NAI ----
      palg = -1d0
      nj = 0
      call fastjetppgenkt(ptrack,np,rr,palg,ptmin,pj,nj,jetvec)
      nb  = 0
      nbb = 0
      do j=1,np
         if(jetvec(j).lt.1 .or. jetvec(j).gt.nj) cycle
         if(idp(j).eq. 5) nb (jetvec(j)) = nb (jetvec(j)) + 1
         if(idp(j).eq.-5) nbb(jetvec(j)) = nbb(jetvec(j)) + 1
      enddo
      do idef=1,2
         njets(idef) = nj
         do i=1,nj
            pjet(0,i,idef)=pj(4,i)
            pjet(1,i,idef)=pj(1,i)
            pjet(2,i,idef)=pj(2,i)
            pjet(3,i,idef)=pj(3,i)
         enddo
         do i=1,nj
            if(idef.eq.1) then
               isbjet = (nb(i)+nbb(i)).ge.1
            else
               isbjet = (nb(i)-nbb(i)).ne.0
            endif
            if(isbjet) call push_bjet(i,nbj(idef),ibj1(idef),ibj2(idef),
     $           pjet(0,1,idef),mxj)
         enddo
      enddo

c     ---- IFN ----
      nj = 0
      call fastjetifn(ptrack,idp,np,rr,ptmin,pj,isb,nj,jetvec)
      njets(3) = nj
      do i=1,nj
         pjet(0,i,3)=pj(4,i)
         pjet(1,i,3)=pj(1,i)
         pjet(2,i,3)=pj(2,i)
         pjet(3,i,3)=pj(3,i)
      enddo
      do i=1,nj
         if(isb(i).eq.1) call push_bjet(i,nbj(3),ibj1(3),ibj2(3),
     $        pjet(0,1,3),mxj)
      enddo

      end


      subroutine push_bjet(i,nbjet,i1,i2,pjet,mxj)
c     keep the indices of the two hardest b jets (fastjet orders by pt,
c     so the first two b jets found are already the hardest ones)
      implicit none
      integer i,nbjet,i1,i2,mxj
      real * 8 pjet(0:3,mxj)
      nbjet = nbjet + 1
      if(nbjet.eq.1) then
         i1 = i
      elseif(nbjet.eq.2) then
         i2 = i
      endif
      end


c=======================================================================
c     small kinematic helpers
c=======================================================================
      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(0:3),y,eta,pt,mass,pv,tiny
      parameter (tiny=1d-8)
      pt=dsqrt(p(1)**2+p(2)**2)
      pv=dsqrt(pt**2+p(3)**2)
      if(pt.lt.tiny) then
         y   = sign(1d8,p(3))
         eta = y
      else
         y   = 0.5d0*dlog((p(0)+p(3))/(p(0)-p(3)))
         eta = 0.5d0*dlog((pv+p(3))/(pv-p(3)))
      endif
      mass=dsqrt(max(0d0,p(0)**2-pv**2))
      end

      subroutine getdydetadphidr(p1,p2,dy,deta,dphi,dry,drr)
      implicit none
      include 'pwhg_math.h'
      real * 8 p1(0:3),p2(0:3),dy,deta,dphi,dry,drr
      real * 8 y1,eta1,pt1,m1,y2,eta2,pt2,m2
      call getyetaptmass(p1,y1,eta1,pt1,m1)
      call getyetaptmass(p2,y2,eta2,pt2,m2)
      dy   = y1-y2
      deta = eta1-eta2
      dphi = (p1(1)*p2(1)+p1(2)*p2(2))/dsqrt(
     $     (p1(1)**2+p1(2)**2)*(p2(1)**2+p2(2)**2))
      dphi = dacos(max(-1d0,min(1d0,dphi)))
      dry  = dsqrt(dy**2  + dphi**2)
      drr  = dsqrt(deta**2+ dphi**2)
      end

      subroutine invmass2(p1,p2,m)
      implicit none
      real * 8 p1(0:3),p2(0:3),m,q(0:3)
      integer mu
      do mu=0,3
         q(mu)=p1(mu)+p2(mu)
      enddo
      m = dsqrt(max(0d0,q(0)**2-q(1)**2-q(2)**2-q(3)**2))
      end
