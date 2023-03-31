      subroutine dshigwid1
c_______________________________________________________________________
c  needs chasct, neusct, sfesct, higsct, vertx.
c  merging of dshwidths and dshigferqcd
c  author: Piero Ullio (ullio@sissa.it) 020917
c          partly based on dshwidth by P. Gondolo and J. Edsjo
c       formulas from higgs hunters guide, 
c       Djouadi, Spira and Zerwas, hep-ph/9511344
c       and Spira, hep-ph/9705337
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsio.h'
      include 'dsidtag.h'
      real*8 mn(4),mi,mj,vvv(12),sqlam,mh
      integer f,i,j,g,k(12),k1,k2
      real*8 hwidthqcd,gratio,gsquared,mscale,fact
      real*8 dsralph3,dsrmq
      real*8 DeltaQCD,Deltat,DeltaPhi,DeltaA,kf
      integer kfer,khiggs,ksfer,ksfer1,ksfer2,khico
      integer kferu,kferd,kckmu,kckmd
      real*8 gAusquared,gAdsquared,muu,mud,lambda,Deltap,Deltam,
     & factr,factl,dsabsq,partial,parsum,mhiggso2
      real*8 gcou,gcod,tau
      complex*16 ampl,AQhH,AsQhH,AQA
      integer nf
      real*8 BigEnf,hwidthqcdnl
      complex*16 sqampl,parampl
      real*8 pi
      parameter (pi=3.141592653589793238d0)


10    format (1x,a2,' -> dsff(',i2,')  partial=',e8.2)
11    format (1x,a2,' -> ww      partial=',e8.2)
12    format (1x,a2,' -> zz      partial=',e8.2)
13    format (1x,a2,' -> n',i1,'n',i1,'    partial=',e8.2)

      do i=1,4
        mn(i)=dabs(mass(kn(i)))
      enddo
      do i=1,4
        do j=1,32
          hdwidth(j,i)=0.0d0
        enddo
      enddo
      i=1
      do g=1,3
        k(i)=knu(g)
        k(i+1)=kl(g)
        k(i+2)=kqu(g)
        k(i+3)=kqd(g)
        i=i+4
      enddo

c======================================================== h2 decay width
      khiggs=kh2
      khico=2
      mhiggso2 = 0.5d0*mass(khiggs)
      width(khiggs) = 0.0d0
c----------------------------------------------------------- into c c-bar
      kfer=kc
      gratio=1.d0
      gsquared=(dcos(alpha))**2*(1.d0+1.d0/tanbe**2)
      mscale=mass(khiggs)
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        fact=DeltaQCD(mscale)+gratio*Deltat(mscale,khiggs,kfer)
        hwidthqcd=3.d0/32.d0/pi*g2weak**2/(mass(kw))**2*mass(khiggs)
     &    *(dsrmq(mscale,kfer))**2*gsquared*fact
      else  
        hwidthqcd=0.d0
      endif
      mscale=2.d0*mass(lsp)
      fact=0.d0
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        fact=DeltaQCD(mscale)+gratio*Deltat(mscale,khiggs,kfer)
      endif  
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(22,khico)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'qdc width = ',
     &    hwidthqcd
      endif  
      if(fact.gt.0.d0) then
c... redefine the couplings
        gl(khiggs,kfer,kfer)=gl(khiggs,kfer,kfer)*dsqrt(fact)
        gl(kfer,khiggs,kfer)=gl(khiggs,kfer,kfer)
        gl(kfer,kfer,khiggs)=gl(khiggs,kfer,kfer)
        gr(khiggs,kfer,kfer)=gr(khiggs,kfer,kfer)*dsqrt(fact)
        gr(kfer,khiggs,kfer)=gr(khiggs,kfer,kfer)
        gr(kfer,kfer,khiggs)=gr(khiggs,kfer,kfer)
      endif
c
c
c----------------------------------------------------------- into b b-bar
      kfer=kb
      gratio=-dcos(alpha)/dsin(alpha)/tanbe
      gsquared=(dsin(alpha))**2*(1.d0+tanbe**2)
      mscale=mass(khiggs)
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        fact=DeltaQCD(mscale)+gratio*Deltat(mscale,khiggs,kfer)
        hwidthqcd=3.d0/32.d0/pi*g2weak**2/(mass(kw))**2*mass(khiggs)
     &    *(dsrmq(mscale,kfer))**2*gsquared*fact
      else  
        hwidthqcd=0.d0
      endif
      mscale=2.d0*mass(lsp)
      fact=0.d0
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        fact=DeltaQCD(mscale)+gratio*Deltat(mscale,khiggs,kfer)
      endif  
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(25,khico)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'qdc width = ',
     &    hwidthqcd
      endif  
      if(fact.gt.0.d0) then
c... redefine the couplings
        gl(khiggs,kfer,kfer)=gl(khiggs,kfer,kfer)*dsqrt(fact)
        gl(kfer,khiggs,kfer)=gl(khiggs,kfer,kfer)
        gl(kfer,kfer,khiggs)=gl(khiggs,kfer,kfer)
        gr(khiggs,kfer,kfer)=gr(khiggs,kfer,kfer)*dsqrt(fact)
        gr(kfer,khiggs,kfer)=gr(khiggs,kfer,kfer)
        gr(kfer,kfer,khiggs)=gr(khiggs,kfer,kfer)
      endif
c
c
c----------------------------------------------------------- into t t-bar
      kfer=kt
      gratio=1.d0
      gsquared=(dcos(alpha))**2*(1.d0+1.d0/tanbe**2)
      mscale=mass(khiggs)
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        kf=sqrt(1.0d0-(2.0d0*mass(kfer)/mscale)**2)
c... I am not sure you have to use dsralph3(mscale) in fact below (PU)
        fact=(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*DeltaPhi(kf))
        hwidthqcd=3.d0/32.d0/pi*g2weak**2/(mass(kw))**2*mass(khiggs)
     &    *gsquared*(mass(kfer))**2*kf**3*fact
      else  
        hwidthqcd=0.d0
      endif
      fact=0.d0
      mscale=2.d0*mass(lsp)
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        kf=sqrt(1.0d0-(2.0d0*mass(kfer)/mscale)**2)
c... I am not sure you have to use dsralph3(mscale) in fact below (PU)
        fact=(mass(kfer))**2/(dsrmq(mscale,kfer))**2
     &    *(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*DeltaPhi(kf))
      endif
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(24,khico)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'qdc width = ',
     &    hwidthqcd
      endif  
      if(fact.gt.0.d0) then
c... redefine the couplings
        gl(khiggs,kfer,kfer)=gl(khiggs,kfer,kfer)*dsqrt(fact)
        gl(kfer,khiggs,kfer)=gl(khiggs,kfer,kfer)
        gl(kfer,kfer,khiggs)=gl(khiggs,kfer,kfer)
        gr(khiggs,kfer,kfer)=gr(khiggs,kfer,kfer)*dsqrt(fact)
        gr(kfer,khiggs,kfer)=gr(khiggs,kfer,kfer)
        gr(kfer,kfer,khiggs)=gr(khiggs,kfer,kfer)
      endif
c
c
c---------------------------------------------- leptons and light quarks
      do f=1,12
        if(f.ne.7.and.f.ne.12.and.f.ne.11) then
        vvv(f)=0.0d0
        if (mass(k(f)).lt.mhiggso2.and.
     &      mass(k(f)).ne.0.0d0) then
          kf = sqrt(1.0d0-(2.0d0*mass(k(f))/mass(khiggs))**2)
          partial = mass(khiggs) *  ncolor(k(f)) * kf**3 *
     &     dsabsq(gr(khiggs,k(f),k(f))) / (8.0d0*pi)
          width(khiggs) = width(khiggs) + partial
          vvv(f)=partial
c          write(*,*) pname(k(f)),pname(k(f)),partial
        endif
        endif
      enddo
c
      hdwidth(19,khico)=vvv(10)   ! tau+ tau-
      hdwidth(17,khico)=vvv(6)   ! mu+ mu-
      hdwidth(23,khico)=vvv(8)   ! s s-bar
c
c
c-------------------------------------------------------------- into w w
      if (mass(kw).lt.mhiggso2) then
        kf = sqrt(0.25d0-(mass(kw)/mass(khiggs))**2)
        partial = mass(khiggs) * kf *
     &   (mass(kw)/mass(khiggs))**2*dsabsq(gr(khiggs,kw,kw))/8.0d0/pi*
     &    ( 2.0d0 + (1.0d0-0.5d0*(mass(kh2)/mass(kw))**2)**2 )
        width(khiggs) = width(khiggs) + partial
        hdwidth(13,khico)=partial
c        write(*,*) pname(kw),pname(kw),partial
      endif
c-------------------------------------------------------------- into z z
      if (mass(kz).lt.mhiggso2) then
        kf = sqrt(0.25d0-(mass(kz)/mass(khiggs))**2)
        partial = mass(khiggs) * kf *
     &    (mass(kz)/mass(khiggs))**2 * dsabsq(gr(khiggs,kz,kz))*
     &    (mass(kw)/mass(kz))**2  /16.0d0/pi *
     &    ( 2.0d0 + (1.0d0-0.5d0*(mass(khiggs)/mass(kz))**2)**2 )
        width(khiggs) = width(khiggs) + partial
        hdwidth(12,khico)=partial
c        write(*,*) pname(kz),pname(kz),partial
      endif
c------------------------------------------------------------ into ni nj
      parsum=0.d0
      do i=1,4
        do j=1,4
          if (mn(i)+mn(j).lt.mass(khiggs)) then
            kf = 0.5d0 *
     &        sqrt(1.0d0-( (mn(i)+mn(j))/mass(khiggs) )**2) *
     &        sqrt(1.0d0-( (mn(i)-mn(j))/mass(khiggs) )**2)
            partial=kf/(4*pi*mass(khiggs))
     &        *dsabsq(gr(khiggs,kn(i),kn(j))) *
     &        (mass(khiggs)**2-(mass(kn(i))+mass(kn(j)))**2)
            if (i.eq.j) partial = partial/2.0d0
            parsum=parsum+partial
c            write(*,*) pname(kn(i)),pname(kn(j)),partial
          endif
        enddo
      enddo
      width(khiggs) = width(khiggs) + parsum
      hdwidth(31,khico)=parsum

c------------------------------------------------------------ into ci cj
      parsum=0.d0
      do i=1,2
        do j=1,2
          mi = mass(kcha(i))
          mj = mass(kcha(j))
          if (abs(mi)+abs(mj).lt.mass(khiggs)) then
            kf = 0.5d0 *
     &        sqrt(1.0d0-( (mi+mj)/mass(khiggs) )**2) *
     &        sqrt(1.0d0-( (mi-mj)/mass(khiggs) )**2)
            partial = kf/(8*pi*mass(khiggs)) *
     &        ( (dsabsq(gl(khiggs,kcha(i),kcha(j)))+
     &           dsabsq(gr(khiggs,kcha(i),kcha(j)))) *
     &             (mass(khiggs)**2-mi**2-mj**2) -
     &          4.0d0*dreal(conjg(gl(khiggs,kcha(i),kcha(j)))*
     &                 gr(khiggs,kcha(i),kcha(j)))*mi*mj )
            parsum=parsum+partial
c            write(*,*) pname(kcha(i)),pname(kcha(j)),partial
          endif
        enddo
      enddo
      width(khiggs) = width(khiggs) + parsum
      hdwidth(32,khico)=parsum
c--------------------------------------------------------------into z h3
      if (mass(kz)+mass(kh3).lt.mass(khiggs)) then
        kf=0.5d0*sqrt(1.0d0-((mass(kz)+mass(kh3))/mass(khiggs))**2) *
     &          sqrt(1.0d0-((mass(kz)-mass(kh3))/mass(khiggs))**2)
        partial = mass(khiggs) * (mass(khiggs)/mass(kz))**2 * kf**3 *
     &   dsabsq(gr(kz,kh3,khiggs)) / (2.0d0*pi)
        hdwidth(10,khico)=partial
        width(khiggs) = width(khiggs) + partial
c        write(*,*) pname(kz),pname(kh3),partial
      endif
c-------------------------------------------------------------- into w h+
      if (mass(kw)+mass(khc).lt.mass(khiggs)) then
        kf=0.5d0 * sqrt(1.0d0-((mass(kw)+mass(khc))/mass(khiggs))**2) *
     &          sqrt(1.0d0-((mass(kw)-mass(khc))/mass(khiggs))**2)
        partial = mass(khiggs) * (mass(khiggs)/mass(kw))**2 * kf**3 *
     &   dsabsq(gr(kw,khc,khiggs)) / (2.0d0*pi)
        partial=partial*2.0d0 ! count both W+ H- and W- H+, corr 020918
        hdwidth(11,khico)=partial
        width(khiggs) = width(khiggs) + partial
c        write(*,*) pname(kw),pname(khc),partial
      endif
c------------------------------------------------------------ into h1 h1
      if (mass(kh1).lt.mhiggso2) then
         kf = sqrt(0.25d0-(mass(kh1)/mass(khiggs))**2)
         partial = (mass(kw)**2/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,kh1,kh1)) / (16.d0*pi)
         width(khiggs) = width(khiggs) + partial
         hdwidth(1,khico)=partial
c         write(*,*) pname(kh1),pname(kh1),partial
      endif
c------------------------------------------------------------ into h3 h3
      if (mass(kh3).lt.mhiggso2) then
         kf = sqrt(0.25d0-(mass(kh3)/mass(khiggs))**2)
         partial = (mass(kw)**2/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,kh3,kh3)) / (16.d0*pi)
         width(khiggs) = width(khiggs) + partial
         hdwidth(4,khico)=partial
c         write(*,*) pname(kh3),pname(kh3),partial
      endif
c------------------------------------------------------------ into h+ h-
      if (mass(khc).lt.mhiggso2) then
         kf = sqrt(0.25d0-(mass(khc)/mass(khiggs))**2)
         partial = (mass(kw)**2/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,khc,khc)) / (16.d0*pi) *2.d0 
         width(khiggs) = width(khiggs) + partial
         hdwidth(7,khico)=partial
c         write(*,*) pname(khc),pname(khc),partial
      endif
c
c
c------------------------------------------------------- into gluon gluon
      gcou=dcos(alpha)*dsqrt(1.d0+1.d0/tanbe**2)
      gcod=-dsin(alpha)*dsqrt(1.d0+tanbe**2)
      mscale=mass(khiggs)
      ampl=dcmplx(0.d0,0.d0)
c top
      kfer=kt  
      tau=(2.d0*mass(kfer)/mass(khiggs))**2
      ampl=ampl+gcou*AQhH(tau)
c bottom
      kfer=kb  
      tau=(2.d0*mass(kfer)/mass(khiggs))**2
      ampl=ampl+gcod*AQhH(tau)
      parampl=ampl
c up-squarks
      do i=1,6
        ksfer= ksqu(i)
        tau=(2.d0*mass(ksfer)/mass(khiggs))**2
        ampl=ampl-mass(kw)/g2weak*gl(khiggs,ksfer,ksfer)*AsQhH(tau)
     &       /(mass(ksfer))**2
      enddo
c down-squarks
      do i=1,6
        ksfer= ksqd(i)
        tau=(2.d0*mass(ksfer)/mass(khiggs))**2
        ampl=ampl-mass(kw)/g2weak*gl(khiggs,ksfer,ksfer)*AsQhH(tau)
     &       /(mass(ksfer))**2
      enddo
      sqampl=ampl-parampl
      hwidthqcd=g2weak**2/(mass(kw))**2/8.d0*(dsralph3(mscale))**2
     &   *mass(khiggs)**3/pi**3/36.d0*dsabsq(ampl)
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(26,khico)=hwidthqcd
c        write(*,*) 'Gluons in final state'
c        write(*,*) 'qdc width = ',hwidthqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Gluons in final state'
        write(*,*) 'qdc width = ',hwidthqcd
      endif  
      if(mscale.lt.mtmt) then
        nf=5
      else
        nf=6
      endif
      BigEnf=95.d0/4.d0-7.d0/6.d0*nf+17.d0/6.d0*dreal(sqampl/parampl)
      hwidthqcdnl=hwidthqcd*dsralph3(mscale)
     &   *(1.d0+BigEnf*dsralph3(mscale)/pi)
      if(hwidthqcdnl.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcdnl
        hdwidth(26,khico)=hdwidth(26,khico)+hwidthqcdnl
c        write(*,*) 'Gluons in final state'
c        write(*,*) 'next to leading qdc term = ',hwidthqcdnl
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Gluons in final state'
        write(*,*) 'next to leading qdc term = ',hwidthqcdnl
      endif
c
c
c------------------------------------------ into sfermion - sfermion star
      parsum=0.d0
c sneutrinos:
c
      do i=1,3
      ksfer=ksnu(i)  
      if (mass(khiggs).ge.2.d0*mass(ksfer)) then
         kf = sqrt(0.25d0-(mass(ksfer)/mass(khiggs))**2)
         partial = (1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer,ksfer)) / (8.d0*pi)
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer),pname(ksfer),partial
      endif
      enddo
c sleptons, i-i final state:
c
      do i=1,6
      ksfer=ksl(i)
      if (mass(khiggs).ge.2.d0*mass(ksfer)) then
         kf = sqrt(0.25d0-(mass(ksfer)/mass(khiggs))**2)
         partial = (1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer,ksfer)) / (8.d0*pi)
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer),pname(ksfer),partial
      endif
      enddo
c sleptons, i-j final states:
c
      do i=1,3
      ksfer1=ksl(i)
      ksfer2=ksl(i+3)
      if (mass(khiggs).ge.mass(ksfer1)+mass(ksfer2)) then
         kf = dsqrt(0.25d0*(1.d0-(mass(ksfer1)/mass(khiggs))**2
     &      -(mass(ksfer2)/mass(khiggs))**2)**2
     &      -(mass(ksfer1)*mass(ksfer2))**2/(mass(khiggs))**4)
         partial = (1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer1,ksfer2)) / (8.d0*pi)*2.d0
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer1),pname(ksfer2),partial
      endif
      enddo
c up-squarks, i-i final states:
c
      do i=1,6
      ksfer=ksqu(i)  
      if (mass(khiggs).ge.2.d0*mass(ksfer)) then
         kf = sqrt(0.25d0-(mass(ksfer)/mass(khiggs))**2)
         partial = 3.d0*(1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer,ksfer)) / (8.d0*pi)
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer),pname(ksfer),partial
      endif
      enddo
c up-squarks, i-j final states:
c
      do i=1,3
      ksfer1=ksqu(i)
      ksfer2=ksqu(i+3)
      if (mass(khiggs).ge.mass(ksfer1)+mass(ksfer2)) then
         kf = dsqrt(0.25d0*(1.d0-(mass(ksfer1)/mass(khiggs))**2
     &      -(mass(ksfer2)/mass(khiggs))**2)**2
     &      -(mass(ksfer1)*mass(ksfer2))**2/(mass(khiggs))**4)
         partial = 3.d0*(1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer1,ksfer2)) / (8.d0*pi)*2.d0
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer1),pname(ksfer2),partial
      endif
      enddo
c down-squarks, i-i final states:
c
      do i=1,6
      ksfer=ksqd(i)  
      if (mass(khiggs).ge.2.d0*mass(ksfer)) then
         kf = sqrt(0.25d0-(mass(ksfer)/mass(khiggs))**2)
         partial = 3.d0*(1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer,ksfer)) / (8.d0*pi)
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer),pname(ksfer),partial
      endif
      enddo
c down-squarks, i-j final states:
c
      do i=1,3
      ksfer1=ksqd(i)
      ksfer2=ksqd(i+3)
      if (mass(khiggs).ge.mass(ksfer1)+mass(ksfer2)) then
         kf = dsqrt(0.25d0*(1.d0-(mass(ksfer1)/mass(khiggs))**2
     &      -(mass(ksfer2)/mass(khiggs))**2)**2
     &      -(mass(ksfer1)*mass(ksfer2))**2/(mass(khiggs))**4)
         partial = 3.d0*(1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer1,ksfer2)) / (8.d0*pi)*2.d0
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer1),pname(ksfer2),partial
      endif
      enddo
c
      width(khiggs)=width(khiggs)+parsum
      hdwidth(30,khico)=parsum
c
c
c======================================================== h1 decay width
      khiggs=kh1
      khico=1
      mhiggso2 = 0.5d0*mass(khiggs)
      width(khiggs) = 0.0d0
c----------------------------------------------------------- into c c-bar
      kfer=kc
      gratio=1.d0
      gsquared=(dsin(alpha))**2*(1.d0+1.d0/tanbe**2)
      mscale=mass(khiggs)
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        fact=DeltaQCD(mscale)+gratio*Deltat(mscale,khiggs,kfer)
        hwidthqcd=3.d0/32.d0/pi*g2weak**2/(mass(kw))**2*mass(khiggs)
     &    *(dsrmq(mscale,kfer))**2*gsquared*fact
      else  
        hwidthqcd=0.d0
      endif
      mscale=2.d0*mass(lsp)
      fact=0.d0
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        fact=DeltaQCD(mscale)+gratio*Deltat(mscale,khiggs,kfer)
      endif  
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(22,khico)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'qdc width = ',
     &    hwidthqcd
      endif  
      if(fact.gt.0.d0) then
c... redefine the couplings
        gl(khiggs,kfer,kfer)=gl(khiggs,kfer,kfer)*dsqrt(fact)
        gl(kfer,khiggs,kfer)=gl(khiggs,kfer,kfer)
        gl(kfer,kfer,khiggs)=gl(khiggs,kfer,kfer)
        gr(khiggs,kfer,kfer)=gr(khiggs,kfer,kfer)*dsqrt(fact)
        gr(kfer,khiggs,kfer)=gr(khiggs,kfer,kfer)
        gr(kfer,kfer,khiggs)=gr(khiggs,kfer,kfer)
      endif
c
c
c----------------------------------------------------------- into b b-bar
      kfer=kb
      gratio=dsin(alpha)/dcos(alpha)/tanbe
      gsquared=(dcos(alpha))**2*(1.d0+tanbe**2)
      mscale=mass(khiggs)
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        fact=DeltaQCD(mscale)+gratio*Deltat(mscale,khiggs,kfer)
        hwidthqcd=3.d0/32.d0/pi*g2weak**2/(mass(kw))**2*mass(khiggs)
     &    *(dsrmq(mscale,kfer))**2*gsquared*fact
      else  
        hwidthqcd=0.d0
      endif
      mscale=2.d0*mass(lsp)
      fact=0.d0
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        fact=DeltaQCD(mscale)+gratio*Deltat(mscale,khiggs,kfer)
      endif  
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(25,khico)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'qdc width = ',
     &    hwidthqcd
      endif  
      if(fact.gt.0.d0) then
c... redefine the couplings
        gl(khiggs,kfer,kfer)=gl(khiggs,kfer,kfer)*dsqrt(fact)
        gl(kfer,khiggs,kfer)=gl(khiggs,kfer,kfer)
        gl(kfer,kfer,khiggs)=gl(khiggs,kfer,kfer)
        gr(khiggs,kfer,kfer)=gr(khiggs,kfer,kfer)*dsqrt(fact)
        gr(kfer,khiggs,kfer)=gr(khiggs,kfer,kfer)
        gr(kfer,kfer,khiggs)=gr(khiggs,kfer,kfer)
      endif
c
c
c----------------------------------------------------------- into t t-bar
      kfer=kt
      gratio=1.d0
      gsquared=(dsin(alpha))**2*(1.d0+1.d0/tanbe**2)
      mscale=mass(khiggs)
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        kf=sqrt(1.0d0-(2.0d0*mass(kfer)/mscale)**2)
c... I am not sure you have to use dsralph3(mscale) in fact below (PU)
        fact=(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*DeltaPhi(kf))
        hwidthqcd=3.d0/32.d0/pi*g2weak**2/(mass(kw))**2*mass(khiggs)
     &    *gsquared*(mass(kfer))**2*kf**3*fact
      else  
        hwidthqcd=0.d0
      endif
      fact=0.d0
      mscale=2.d0*mass(lsp)
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        kf=sqrt(1.0d0-(2.0d0*mass(kfer)/mscale)**2)
c... I am not sure you have to use dsralph3(mscale) in fact below (PU)
        fact=(mass(kfer))**2/(dsrmq(mscale,kfer))**2
     &    *(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*DeltaPhi(kf))
      endif
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(24,khico)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'qdc width = ',
     &    hwidthqcd
      endif  
      if(fact.gt.0.d0) then
c... redefine the couplings
        gl(khiggs,kfer,kfer)=gl(khiggs,kfer,kfer)*dsqrt(fact)
        gl(kfer,khiggs,kfer)=gl(khiggs,kfer,kfer)
        gl(kfer,kfer,khiggs)=gl(khiggs,kfer,kfer)
        gr(khiggs,kfer,kfer)=gr(khiggs,kfer,kfer)*dsqrt(fact)
        gr(kfer,khiggs,kfer)=gr(khiggs,kfer,kfer)
        gr(kfer,kfer,khiggs)=gr(khiggs,kfer,kfer)
      endif
c
c
c---------------------------------------------- leptons and light quarks
      do f=1,12
        if(f.ne.7.and.f.ne.12.and.f.ne.11) then
        vvv(f)=0.0d0
        if (mass(k(f)).lt.mhiggso2.and.
     &      mass(k(f)).ne.0.0d0) then
          kf = sqrt(1.0d0-(2.0d0*mass(k(f))/mass(khiggs))**2)
          partial = mass(khiggs) *  ncolor(k(f)) * kf**3 *
     &     dsabsq(gr(khiggs,k(f),k(f))) / (8.0d0*pi)
          width(khiggs) = width(khiggs) + partial
          vvv(f)=partial
c          write(*,*) pname(k(f)),pname(k(f)),partial
        endif
        endif
      enddo
c
      hdwidth(19,khico)=vvv(10)   ! tau+ tau-
      hdwidth(17,khico)=vvv(6)   ! mu+ mu-
      hdwidth(23,khico)=vvv(8)   ! s s-bar
c
c
c-------------------------------------------------------------- into w w
      if (mass(kw).lt.mhiggso2) then
        kf = sqrt(0.25d0-(mass(kw)/mass(khiggs))**2)
        partial = mass(khiggs) * kf *
     &    (mass(kw)/mass(khiggs))**2*dsabsq(gr(khiggs,kw,kw))/8.0d0/pi*
     &    ( 2.0d0 + (1.0d0-0.5d0*(mass(khiggs)/mass(kw))**2)**2 )
        width(khiggs) = width(khiggs) + partial
        hdwidth(13,khico)=partial
c        write(*,*) pname(kw),pname(kw),partial
      endif
c-------------------------------------------------------------- into z z
      if (mass(kz).lt.mhiggso2) then
        kf = sqrt(0.25d0-(mass(kz)/mass(khiggs))**2)
        partial = mass(khiggs) * kf *
     &    (mass(kz)/mass(khiggs))**2 * dsabsq(gr(khiggs,kz,kz))*
     &    (mass(kw)/mass(kz))**2 /16.0d0/pi *
     &    ( 2.0d0 + (1.0d0-0.5d0*(mass(khiggs)/mass(kz))**2)**2 )
        width(khiggs) = width(khiggs) + partial
        hdwidth(12,khico)=partial
c        write(*,*) pname(kz),pname(kz),partial
      endif
c------------------------------------------------------------ into ni nj
      parsum=0.d0
      do i=1,4
        do j=1,4
          if (mn(i)+mn(j).lt.mass(khiggs)) then
            kf = 0.5d0 *
     &        sqrt(1.0d0-( (mn(i)+mn(j))/mass(khiggs) )**2) *
     &        sqrt(1.0d0-( (mn(i)-mn(j))/mass(khiggs) )**2)
            partial=kf/(4*pi*mass(khiggs))
     &        *dsabsq(gr(khiggs,kn(i),kn(j))) *
     &        (mass(khiggs)**2-(mass(kn(i))+mass(kn(j)))**2)
            if (i.eq.j) partial = partial/2.0d0
            parsum=parsum+partial
c            write(*,*) pname(kn(i)),pname(kn(j)),partial
          endif
        enddo
      enddo
      width(khiggs) = width(khiggs) + parsum
      hdwidth(31,khico)=parsum
c------------------------------------------------------------ into ci cj
      parsum=0.d0
      do i=1,2
        do j=1,2
          mi = mass(kcha(i))
          mj = mass(kcha(j))
          if (abs(mi)+abs(mj).lt.mass(khiggs)) then
            kf = 0.5d0 *
     &        sqrt(1.0d0-( (mi+mj)/mass(khiggs) )**2) *
     &        sqrt(1.0d0-( (mi-mj)/mass(khiggs) )**2)
            partial = kf/(8*pi*mass(khiggs)) *
     &        ( (dsabsq(gl(khiggs,kcha(i),kcha(j)))+
     &           dsabsq(gr(khiggs,kcha(i),kcha(j)))) *
     &             (mass(khiggs)**2-mi**2-mj**2) -
     &          4.0d0*dreal(conjg(gl(khiggs,kcha(i),kcha(j)))*
     &                 gr(khiggs,kcha(i),kcha(j)))*mi*mj )
            parsum=parsum+partial
c            write(*,*) pname(kcha(i)),pname(kcha(j)),partial
          endif
        enddo
      enddo
      width(khiggs) = width(khiggs) + parsum
      hdwidth(32,khico)=parsum
c------------------------------------------------------------- into z h3
      if (mass(kz)+mass(kh3).lt.mass(khiggs)) then
        kf=0.5d0 * sqrt(1.0d0-((mass(kz)+mass(kh3))/mass(khiggs))**2) *
     &          sqrt(1.0d0-((mass(kz)-mass(kh3))/mass(khiggs))**2)
        partial = mass(khiggs) * (mass(khiggs)/mass(kz))**2 * kf**3 *
     &   dsabsq(gr(kz,kh3,khiggs)) / (2.0d0*pi)
        width(khiggs) = width(khiggs) + partial
        hdwidth(10,khico)=partial
c        write(*,*) pname(kz),pname(kh3),partial
      endif
c-------------------------------------------------------------- into w h+
      if (mass(kw)+mass(khc).lt.mass(khiggs)) then
        kf=0.5d0 * sqrt(1.0d0-((mass(kw)+mass(khc))/mass(khiggs))**2) *
     &          sqrt(1.0d0-((mass(kw)-mass(khc))/mass(khiggs))**2)
        partial = mass(khiggs) * (mass(khiggs)/mass(kw))**2 * kf**3 *
     &   dsabsq(gr(kw,khc,khiggs)) / (2.0d0*pi)
        partial=partial*2.0d0 ! count both W+ H- and W- H+, corr 020918
        width(khiggs) = width(khiggs) + partial
        hdwidth(11,khico)=partial
c        write(*,*) pname(kw),pname(khc),partial
      endif
c------------------------------------------------------------ into h2 h2
      if (mass(kh2).lt.mhiggso2) then
         kf = sqrt(0.25d0-(mass(kh2)/mass(khiggs))**2)
         partial = (mass(kw)**2/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,kh2,kh2)) / (16.d0*pi)
         width(khiggs) = width(khiggs) + partial
         hdwidth(3,khico)=partial
c         write(*,*) pname(kh2),pname(kh2),partial
      endif
c------------------------------------------------------------ into h3 h3
      if (mass(kh3).lt.mhiggso2) then
         kf = sqrt(0.25d0-(mass(kh3)/mass(khiggs))**2)
         partial = (mass(kw)**2/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,kh3,kh3)) / (16.d0*pi)
         width(khiggs) = width(khiggs) + partial
         hdwidth(4,khico)=partial
c         write(*,*) pname(kh3),pname(kh3),partial
      endif
c------------------------------------------------------------- into h+ h-
      if (mass(khc).lt.mhiggso2) then
         kf = sqrt(0.25d0-(mass(khc)/mass(khiggs))**2)
         partial = (mass(kw)**2/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,khc,khc)) / (16.d0*pi) *2.d0
         width(khiggs) = width(khiggs) + partial
         hdwidth(7,khico)=partial
c         write(*,*) pname(khc),pname(khc),partial
      endif
c
c------------------------------------------------------- into gluon gluon
      gcou=dsin(alpha)*dsqrt(1.d0+1.d0/tanbe**2)
      gcod=dcos(alpha)*dsqrt(1.d0+tanbe**2)
      mscale=mass(khiggs)
      ampl=dcmplx(0.d0,0.d0)
c top
      kfer=kt  
      tau=(2.d0*mass(kfer)/mass(khiggs))**2
      ampl=ampl+gcou*AQhH(tau)
c bottom
      kfer=kb  
      tau=(2.d0*mass(kfer)/mass(khiggs))**2
      ampl=ampl+gcod*AQhH(tau)
      parampl=ampl
c up-squarks
      do i=1,6
        ksfer= ksqu(i)
        tau=(2.d0*mass(ksfer)/mass(khiggs))**2
        ampl=ampl-mass(kw)/g2weak*gl(khiggs,ksfer,ksfer)*AsQhH(tau)
     &       /(mass(ksfer))**2
      enddo
c down-squarks
      do i=1,6
        ksfer= ksqd(i)
        tau=(2.d0*mass(ksfer)/mass(khiggs))**2
        ampl=ampl-mass(kw)/g2weak*gl(khiggs,ksfer,ksfer)*AsQhH(tau)
     &       /(mass(ksfer))**2
      enddo
      sqampl=ampl-parampl
      hwidthqcd=g2weak**2/(mass(kw))**2/8.d0*(dsralph3(mscale))**2
     &   *mass(khiggs)**3/pi**3/36.d0*dsabsq(ampl)
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(26,khico)=hwidthqcd
c        write(*,*) 'Gluons in final state'
c        write(*,*) 'qdc width = ',hwidthqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Gluons in final state'
        write(*,*) 'qdc width = ',hwidthqcd
      endif  
      if(mscale.lt.mtmt) then
        nf=5
      else
        nf=6
      endif
      BigEnf=95.d0/4.d0-7.d0/6.d0*nf+17.d0/6.d0*dreal(sqampl/parampl)
      hwidthqcdnl=hwidthqcd*dsralph3(mscale)
     &   *(1.d0+BigEnf*dsralph3(mscale)/pi)
      if(hwidthqcdnl.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcdnl
        hdwidth(26,khico)=hdwidth(26,khico)+hwidthqcdnl
c        write(*,*) 'Gluons in final state'
c        write(*,*) 'next to leading qdc term = ',hwidthqcdnl
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Gluons in final state'
        write(*,*) 'next to leading qdc term = ',hwidthqcdnl
      endif
c
c
c------------------------------------------ into sfermion - sfermion star
      parsum=0.d0
c sneutrinos:
c
      do i=1,3
      ksfer=ksnu(i)  
      if (mass(khiggs).ge.2.d0*mass(ksfer)) then
         kf = sqrt(0.25d0-(mass(ksfer)/mass(khiggs))**2)
         partial = (1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer,ksfer)) / (8.d0*pi)
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer),pname(ksfer),partial
      endif
      enddo
c sleptons, i-i final state:
c
      do i=1,6
      ksfer=ksl(i)
      if (mass(khiggs).ge.2.d0*mass(ksfer)) then
         kf = sqrt(0.25d0-(mass(ksfer)/mass(khiggs))**2)
         partial = (1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer,ksfer)) / (8.d0*pi)
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer),pname(ksfer),partial
      endif
      enddo
c sleptons, i-j final states:
c
      do i=1,3
      ksfer1=ksl(i)
      ksfer2=ksl(i+3)
      if (mass(khiggs).ge.mass(ksfer1)+mass(ksfer2)) then
         kf = dsqrt(0.25d0*(1.d0-(mass(ksfer1)/mass(khiggs))**2
     &      -(mass(ksfer2)/mass(khiggs))**2)**2
     &      -(mass(ksfer1)*mass(ksfer2))**2/(mass(khiggs))**4)
         partial = (1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer1,ksfer2)) / (8.d0*pi)*2.d0
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer1),pname(ksfer2),partial
      endif
      enddo
c up-squarks, i-i final states:
c
      do i=1,6
      ksfer=ksqu(i)  
      if (mass(khiggs).ge.2.d0*mass(ksfer)) then
         kf = sqrt(0.25d0-(mass(ksfer)/mass(khiggs))**2)
         partial = 3.d0*(1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer,ksfer)) / (8.d0*pi)
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer),pname(ksfer),partial
      endif
      enddo
c up-squarks, i-j final states:
c
      do i=1,3
      ksfer1=ksqu(i)
      ksfer2=ksqu(i+3)
      if (mass(khiggs).ge.mass(ksfer1)+mass(ksfer2)) then
         kf = dsqrt(0.25d0*(1.d0-(mass(ksfer1)/mass(khiggs))**2
     &      -(mass(ksfer2)/mass(khiggs))**2)**2
     &      -(mass(ksfer1)*mass(ksfer2))**2/(mass(khiggs))**4)
         partial = 3.d0*(1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer1,ksfer2)) / (8.d0*pi)*2.d0
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer1),pname(ksfer2),partial
      endif
      enddo
c down-squarks, i-i final states:
c
      do i=1,6
      ksfer=ksqd(i)  
      if (mass(khiggs).ge.2.d0*mass(ksfer)) then
         kf = sqrt(0.25d0-(mass(ksfer)/mass(khiggs))**2)
         partial = 3.d0*(1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer,ksfer)) / (8.d0*pi)
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer),pname(ksfer),partial
      endif
      enddo
c down-squarks, i-j final states:
c
      do i=1,3
      ksfer1=ksqd(i)
      ksfer2=ksqd(i+3)
      if (mass(khiggs).ge.mass(ksfer1)+mass(ksfer2)) then
         kf = dsqrt(0.25d0*(1.d0-(mass(ksfer1)/mass(khiggs))**2
     &      -(mass(ksfer2)/mass(khiggs))**2)**2
     &      -(mass(ksfer1)*mass(ksfer2))**2/(mass(khiggs))**4)
         partial = 3.d0*(1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer1,ksfer2)) / (8.d0*pi)*2.d0
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer1),pname(ksfer2),partial
      endif
      enddo
c
      width(khiggs)=width(khiggs)+parsum
      hdwidth(30,khico)=parsum
c
c
c======================================================== h3 decay width
      khiggs=kh3
      khico=3
      mhiggso2 = 0.5d0*mass(khiggs)
      width(khiggs) = 0.0d0
c----------------------------------------------------------- into c c-bar
      kfer=kc
      gratio=1.d0
      gsquared=1.d0/tanbe**2
      mscale=mass(khiggs)
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        fact=DeltaQCD(mscale)+gratio*Deltat(mscale,khiggs,kfer)
        hwidthqcd=3.d0/32.d0/pi*g2weak**2/(mass(kw))**2*mass(khiggs)
     &    *(dsrmq(mscale,kfer))**2*gsquared*fact
      else  
        hwidthqcd=0.d0
      endif
      mscale=2.d0*mass(lsp)
      fact=0.d0
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        fact=DeltaQCD(mscale)+gratio*Deltat(mscale,khiggs,kfer)
      endif  
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(22,khico)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'qdc width = ',
     &    hwidthqcd
      endif  
      if(fact.gt.0.d0) then
c... redefine the couplings
        gl(khiggs,kfer,kfer)=gl(khiggs,kfer,kfer)*dsqrt(fact)
        gl(kfer,khiggs,kfer)=gl(khiggs,kfer,kfer)
        gl(kfer,kfer,khiggs)=gl(khiggs,kfer,kfer)
        gr(khiggs,kfer,kfer)=gr(khiggs,kfer,kfer)*dsqrt(fact)
        gr(kfer,khiggs,kfer)=gr(khiggs,kfer,kfer)
        gr(kfer,kfer,khiggs)=gr(khiggs,kfer,kfer)
      endif
c
c
c----------------------------------------------------------- into b b-bar
      kfer=kb
      gratio=1.d0/tanbe**2
      gsquared=tanbe**2
      mscale=mass(khiggs)
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        fact=DeltaQCD(mscale)+gratio*Deltat(mscale,khiggs,kfer)
        hwidthqcd=3.d0/32.d0/pi*g2weak**2/(mass(kw))**2*mass(khiggs)
     &    *(dsrmq(mscale,kfer))**2*gsquared*fact
      else  
        hwidthqcd=0.d0
      endif
      mscale=2.d0*mass(lsp)
      fact=0.d0
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        fact=DeltaQCD(mscale)+gratio*Deltat(mscale,khiggs,kfer)
      endif  
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(25,khico)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'qdc width = ',
     &    hwidthqcd
      endif  
      if(fact.gt.0.d0) then
c... redefine the couplings
        gl(khiggs,kfer,kfer)=gl(khiggs,kfer,kfer)*dsqrt(fact)
        gl(kfer,khiggs,kfer)=gl(khiggs,kfer,kfer)
        gl(kfer,kfer,khiggs)=gl(khiggs,kfer,kfer)
        gr(khiggs,kfer,kfer)=gr(khiggs,kfer,kfer)*dsqrt(fact)
        gr(kfer,khiggs,kfer)=gr(khiggs,kfer,kfer)
        gr(kfer,kfer,khiggs)=gr(khiggs,kfer,kfer)
      endif
c
c
c----------------------------------------------------------- into t t-bar
      kfer=kt
      gratio=1.d0
      gsquared=1.d0/tanbe**2
      mscale=mass(khiggs)
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        kf=sqrt(1.0d0-(2.0d0*mass(kfer)/mscale)**2)
c... I am not sure you have to use dsralph3(mscale) in fact below (PU)
        fact=(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*DeltaA(kf))
        hwidthqcd=3.d0/32.d0/pi*g2weak**2/(mass(kw))**2*mass(khiggs)
     &    *gsquared*(mass(kfer))**2*kf*fact
      else  
        hwidthqcd=0.d0
      endif
      fact=0.d0
      mscale=2.d0*mass(lsp)
      if(mscale.gt.2.d0*mass(kfer)+1.d0) then
        kf=sqrt(1.0d0-(2.0d0*mass(kfer)/mscale)**2)
c... I am not sure you have to use dsralph3(mscale) in fact below (PU)
        fact=(mass(kfer))**2/(dsrmq(mscale,kfer))**2*
     &    (1.d0+4.d0/3.d0*dsralph3(mscale)/pi*DeltaA(kf))
      endif
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(24,khico)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'qdc width = ',
     &    hwidthqcd
      endif  
      if(fact.gt.0.d0) then
c... redefine the couplings
        gl(khiggs,kfer,kfer)=gl(khiggs,kfer,kfer)*dsqrt(fact)
        gl(kfer,khiggs,kfer)=gl(khiggs,kfer,kfer)
        gl(kfer,kfer,khiggs)=gl(khiggs,kfer,kfer)
        gr(khiggs,kfer,kfer)=gr(khiggs,kfer,kfer)*dsqrt(fact)
        gr(kfer,khiggs,kfer)=gr(khiggs,kfer,kfer)
        gr(kfer,kfer,khiggs)=gr(khiggs,kfer,kfer)
      endif
c
c
c---------------------------------------------- leptons and light quarks
      do f=1,12
        if(f.ne.7.and.f.ne.12.and.f.ne.11) then
        vvv(f)=0.0d0
        if (mass(k(f)).lt.mhiggso2.and.
     &      mass(k(f)).ne.0.0d0) then
          kf = sqrt(1.0d0-(2.0d0*mass(k(f))/mass(khiggs))**2)
          partial = mass(khiggs) *  ncolor(k(f)) * kf *
     &     dsabsq(gr(khiggs,k(f),k(f))) / (8.0d0*pi)
          width(khiggs) = width(khiggs) + partial
          vvv(f)=partial
c          write(*,*) pname(k(f)),pname(k(f)),partial
        endif
        endif
      enddo
c
      hdwidth(19,khico)=vvv(10)   ! tau+ tau-
      hdwidth(17,khico)=vvv(6)   ! mu+ mu-
      hdwidth(23,khico)=vvv(8)   ! s s-bar
c
c
c------------------------------------------------------------ into ni nj
      parsum=0.d0
      do i=1,4
        do j=1,4
          if (mn(i)+mn(j).lt.mass(khiggs)) then
            kf = 0.5d0 *
     &        sqrt(1.0d0-( (mn(i)+mn(j))/mass(khiggs) )**2) *
     &        sqrt(1.0d0-( (mn(i)-mn(j))/mass(khiggs) )**2)
            partial=kf/(4*pi*mass(khiggs))
     &        *dsabsq(gr(khiggs,kn(i),kn(j))) *
     &        (mass(khiggs)**2-(mass(kn(i))-mass(kn(j)))**2)
            if (i.eq.j) partial = partial/2.0d0
            parsum=parsum+partial
c            write(*,*) pname(kn(i)),pname(kn(j)),partial
          endif
        enddo
      enddo
      width(khiggs) = width(khiggs) + parsum
      hdwidth(31,khico)=parsum
c------------------------------------------------------------ into ci cj
      parsum=0.d0
      do i=1,2
        do j=1,2
          mi = mass(kcha(i))
          mj = mass(kcha(j))
          if (abs(mi)+abs(mj).lt.mass(khiggs)) then
            kf = 0.5d0 *
     &        sqrt(1.0d0-( (mi+mj)/mass(khiggs) )**2) *
     &        sqrt(1.0d0-( (mi-mj)/mass(khiggs) )**2)
            partial = kf/(8*pi*mass(khiggs)) *
     &        ( (dsabsq(gl(khiggs,kcha(i),kcha(j)))+
     &           dsabsq(gr(khiggs,kcha(i),kcha(j)))) *
     &             (mass(khiggs)**2-mi**2-mj**2) +
     &          4.0d0*dreal(conjg(gl(khiggs,kcha(i),kcha(j)))*
     &                 gr(khiggs,kcha(i),kcha(j)))*mi*mj )
            parsum=parsum+partial
c            write(*,*) pname(kcha(i)),pname(kcha(j)),partial
          endif
        enddo
      enddo
      width(khiggs) = width(khiggs) + parsum
      hdwidth(32,khico)=parsum
c------------------------------------------------------------- into z h1
      if (mass(kz)+mass(kh1).lt.mass(khiggs)) then
        kf=0.5d0 * sqrt(1.0d0-((mass(kz)+mass(kh1))/mass(khiggs))**2) *
     &          sqrt(1.0d0-((mass(kz)-mass(kh1))/mass(khiggs))**2)
        partial = mass(khiggs) * (mass(khiggs)/mass(kz))**2 * kf**3 *
     &   dsabsq(gr(kz,khiggs,kh1)) / (2.0d0*pi)
        width(khiggs) = width(khiggs) + partial
        hdwidth(8,khico)=partial
c        write(*,*) pname(kz),pname(kh1),partial
      endif
c------------------------------------------------------------- into z h2
      if (mass(kz)+mass(kh2).lt.mass(khiggs)) then
        kf=0.5d0 * sqrt(1.0d0-((mass(kz)+mass(kh2))/mass(khiggs))**2) *
     &          sqrt(1.0d0-((mass(kz)-mass(kh2))/mass(khiggs))**2)
        partial = mass(khiggs) * (mass(khiggs)/mass(kz))**2 * kf**3 *
     &   dsabsq(gr(kz,khiggs,kh2)) / (2.0d0*pi)
        width(khiggs) = width(khiggs) + partial
        hdwidth(9,khico)=partial
c        write(*,*) pname(kz),pname(kh2),partial
      endif


c-------------------------------------------------------------- into w h+
      if (mass(kw)+mass(khc).lt.mass(khiggs)) then
        kf=0.5d0 * sqrt(1.0d0-((mass(kw)+mass(khc))/mass(khiggs))**2) *
     &          sqrt(1.0d0-((mass(kw)-mass(khc))/mass(khiggs))**2)
        partial = mass(khiggs) * (mass(khiggs)/mass(kw))**2 * kf**3 *
     &   dsabsq(gr(kw,khc,khiggs)) / (2.0d0*pi)
        partial=partial*2.0d0 ! count both W+ H- and W- H+, corr 020918
        width(khiggs) = width(khiggs) + partial
        hdwidth(11,khico)=partial
c        write(*,*) pname(kw),pname(khc),partial
      endif
c
c
c------------------------------------------------------- into gluon gluon
      gcou=1.d0/tanbe
      gcod=tanbe
      mscale=mass(khiggs)
      ampl=dcmplx(0.d0,0.d0)
c top
      kfer=kt  
      tau=(2.d0*mass(kfer)/mass(khiggs))**2
      ampl=ampl+gcou*AQA(tau)
c bottom
      kfer=kb  
      tau=(2.d0*mass(kfer)/mass(khiggs))**2
      ampl=ampl+gcod*AQA(tau)
      hwidthqcd=g2weak**2/(mass(kw))**2/8.d0*(dsralph3(mscale))**2
     &   *mass(khiggs)**3/pi**3/16.d0*dsabsq(ampl)
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(26,khico)=hwidthqcd
c        write(*,*) 'Gluons in final state'
c        write(*,*) 'qdc width = ',hwidthqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Gluons in final state'
        write(*,*) 'qdc width = ',hwidthqcd
      endif  
      if(mscale.lt.mtmt) then
        nf=5
      else
        nf=6
      endif
      BigEnf=95.d0/4.d0-7.d0/6.d0*nf
      hwidthqcdnl=hwidthqcd*dsralph3(mscale)
     &   *(1.d0+BigEnf*dsralph3(mscale)/pi)
      if(hwidthqcdnl.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcdnl
        hdwidth(26,khico)=hdwidth(26,khico)+hwidthqcdnl
c        write(*,*) 'Gluons in final state'
c        write(*,*) 'next to leading qdc term = ',hwidthqcdnl
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Gluons in final state'
        write(*,*) 'next to leading qdc term = ',hwidthqcdnl
      endif
c
c
c------------------------------------------ into sfermion - sfermion star
      parsum=0.d0
c sneutrinos:
c
      do i=1,3
      ksfer=ksnu(i)  
      if (mass(khiggs).ge.2.d0*mass(ksfer)) then
         kf = sqrt(0.25d0-(mass(ksfer)/mass(khiggs))**2)
         partial = (1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer,ksfer)) / (8.d0*pi)
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer),pname(ksfer),partial
      endif
      enddo
c sleptons, i-i final state:
c
      do i=1,6
      ksfer=ksl(i)
      if (mass(khiggs).ge.2.d0*mass(ksfer)) then
         kf = sqrt(0.25d0-(mass(ksfer)/mass(khiggs))**2)
         partial = (1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer,ksfer)) / (8.d0*pi)
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer),pname(ksfer),partial
      endif
      enddo
c sleptons, i-j final states:
c
      do i=1,3
      ksfer1=ksl(i)
      ksfer2=ksl(i+3)
      if (mass(khiggs).ge.mass(ksfer1)+mass(ksfer2)) then
         kf = dsqrt(0.25d0*(1.d0-(mass(ksfer1)/mass(khiggs))**2
     &      -(mass(ksfer2)/mass(khiggs))**2)**2
     &      -(mass(ksfer1)*mass(ksfer2))**2/(mass(khiggs))**4)
         partial = (1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer1,ksfer2)) / (8.d0*pi)*2.d0
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer1),pname(ksfer2),partial
      endif
      enddo
c up-squarks, i-i final states:
c
      do i=1,6
      ksfer=ksqu(i)  
      if (mass(khiggs).ge.2.d0*mass(ksfer)) then
         kf = sqrt(0.25d0-(mass(ksfer)/mass(khiggs))**2)
         partial = 3.d0*(1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer,ksfer)) / (8.d0*pi)
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer),pname(ksfer),partial
      endif
      enddo
c up-squarks, i-j final states:
c
      do i=1,3
      ksfer1=ksqu(i)
      ksfer2=ksqu(i+3)
      if (mass(khiggs).ge.mass(ksfer1)+mass(ksfer2)) then
         kf = dsqrt(0.25d0*(1.d0-(mass(ksfer1)/mass(khiggs))**2
     &      -(mass(ksfer2)/mass(khiggs))**2)**2
     &      -(mass(ksfer1)*mass(ksfer2))**2/(mass(khiggs))**4)
         partial = 3.d0*(1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer1,ksfer2)) / (8.d0*pi)*2.d0
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer1),pname(ksfer2),partial
      endif
      enddo
c down-squarks, i-i final states:
c
      do i=1,6
      ksfer=ksqd(i)  
      if (mass(khiggs).ge.2.d0*mass(ksfer)) then
         kf = sqrt(0.25d0-(mass(ksfer)/mass(khiggs))**2)
         partial = 3.d0*(1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer,ksfer)) / (8.d0*pi)
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer),pname(ksfer),partial
      endif
      enddo
c down-squarks, i-j final states:
c
      do i=1,3
      ksfer1=ksqd(i)
      ksfer2=ksqd(i+3)
      if (mass(khiggs).ge.mass(ksfer1)+mass(ksfer2)) then
         kf = dsqrt(0.25d0*(1.d0-(mass(ksfer1)/mass(khiggs))**2
     &      -(mass(ksfer2)/mass(khiggs))**2)**2
     &      -(mass(ksfer1)*mass(ksfer2))**2/(mass(khiggs))**4)
         partial = 3.d0*(1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer1,ksfer2)) / (8.d0*pi)*2.d0
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer1),pname(ksfer2),partial
      endif
      enddo
c
      width(khiggs)=width(khiggs)+parsum
      hdwidth(30,khico)=parsum
c
c
c======================================================== h+ decay width
      khiggs=khc
      khico=4
      mhiggso2 = 0.5d0*mass(khiggs)
      width(khiggs) = 0.0d0
      mh=mass(khiggs)
c----------------------------------------------------------- into quarks
      i=1
      do j=1,2
        vvv((i-1)*3+j)=0.0d0
        k1=kqu(i)
        k2=kqd(j)
        mi=mass(k1)
        mj=mass(k2)
        if ((mi+mj).lt.mh) then
          sqlam=sqrt((mi**2+mj**2-mh**2)**2-4*mi**2*mj**2)
          partial=3.0d0*sqlam/(16.0d0*pi*mh**3)*
     &        ((mh**2-mi**2-mj**2)*
     &        (dsabsq(gr(khiggs,k2,k1))+dsabsq(gl(khiggs,k2,k1)))-
     &        4*mi*mj*gl(khiggs,k2,k1)*gr(khiggs,k2,k1))
          width(khiggs)=width(khiggs)+partial
          vvv((i-1)*3+j)=partial
c            write(*,*) pname(k1),pname(k2),partial
        endif
      enddo
      do i=1,2
        hdwidth(i,khico)=vvv(i)
      enddo
c----------------------------------------------------------- into u b-bar
c i=1 (u index), j=3 (d index), (i-1)*3+j= 3
      kferu=ku
      kckmu=1
      kferd=kb
      kckmd=3
      gAusquared=1.d0/tanbe**2
      gAdsquared=tanbe**2
      mscale=mass(khiggs)
      if(mscale.gt.mass(kferu)+mass(kferd)+1.d0) then
        fact=DeltaQCD(mscale)
        hwidthqcd=3.d0/32.d0/pi*g2weak**2/(mass(kw))**2*mass(khiggs)
     &    *(ckm(kckmu,kckmd))**2
     &    *((mass(kferu))**2*gAusquared
     &      +(dsrmq(mscale,kferd))**2*gAdsquared)*fact
      else  
        hwidthqcd=0.d0
      endif
      mscale=2.d0*mass(lsp)
      fact=0.d0
      if(mscale.gt.mass(kferu)+mass(kferd)+1.d0) then
        fact=DeltaQCD(mscale)
      endif  
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(3,khico)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kferu)
c     &      ,pname(kferd)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kferu)
     &      ,pname(kferd)
        write(*,*) 'qdc width = ',
     &    hwidthqcd
      endif  
      if(fact.gt.0.d0) then
c... redefine the couplings
c        write(*,*) 'vertex corr',dsqrt(fact)
        gl(khiggs,kferd,kferu)=gl(khiggs,kferd,kferu)*dsqrt(fact)
        gl(kferd,khiggs,kferu)=gl(khiggs,kferd,kferu)
        gl(kferd,kferu,khiggs)=gl(khiggs,kferd,kferu)
        gl(khiggs,kferu,kferd)=gl(khiggs,kferu,kferd)*dsqrt(fact)
        gl(kferu,khiggs,kferd)=gl(khiggs,kferu,kferd)
        gl(kferu,kferd,khiggs)=gl(khiggs,kferu,kferd)
        gr(khiggs,kferd,kferu)=gr(khiggs,kferd,kferu)*dsqrt(fact)
        gr(kferd,khiggs,kferu)=gr(khiggs,kferd,kferu)
        gr(kferd,kferu,khiggs)=gr(khiggs,kferd,kferu)
        gr(khiggs,kferu,kferd)=gr(khiggs,kferu,kferd)*dsqrt(fact)
        gr(kferu,khiggs,kferd)=gr(khiggs,kferu,kferd)
        gr(kferu,kferd,khiggs)=gr(khiggs,kferu,kferd)
      endif
c
c
c----------------------------------------------------------- into c d-bar
c i=2 (u index), j=1 (d index), (i-1)*3+j= 4
      kferu=kc
      kckmu=2
      kferd=kd
      kckmd=1
      gAusquared=1.d0/tanbe**2
      gAdsquared=tanbe**2
      mscale=mass(khiggs)
      if(mscale.gt.mass(kferu)+mass(kferd)+1.d0) then
        fact=DeltaQCD(mscale)
        hwidthqcd=3.d0/32.d0/pi*g2weak**2/(mass(kw))**2*mass(khiggs)
     &    *(ckm(kckmu,kckmd))**2
     &    *((dsrmq(mscale,kferu))**2*gAusquared
     &      +(mass(kferd))**2*gAdsquared)*fact
      else  
        hwidthqcd=0.d0
      endif
      mscale=2.d0*mass(lsp)
      fact=0.d0
      if(mscale.gt.mass(kferu)+mass(kferd)+1.d0) then
        fact=DeltaQCD(mscale)
      endif  
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(4,khico)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kferu)
c     &      ,pname(kferd)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kferu)
     &      ,pname(kferd)
        write(*,*) 'qdc width = ',
     &    hwidthqcd
      endif  
      if(fact.gt.0.d0) then
c... redefine the couplings
c        write(*,*) 'vertex corr',dsqrt(fact)
        gl(khiggs,kferd,kferu)=gl(khiggs,kferd,kferu)*dsqrt(fact)
        gl(kferd,khiggs,kferu)=gl(khiggs,kferd,kferu)
        gl(kferd,kferu,khiggs)=gl(khiggs,kferd,kferu)
        gl(khiggs,kferu,kferd)=gl(khiggs,kferu,kferd)*dsqrt(fact)
        gl(kferu,khiggs,kferd)=gl(khiggs,kferu,kferd)
        gl(kferu,kferd,khiggs)=gl(khiggs,kferu,kferd)
        gr(khiggs,kferd,kferu)=gr(khiggs,kferd,kferu)*dsqrt(fact)
        gr(kferd,khiggs,kferu)=gr(khiggs,kferd,kferu)
        gr(kferd,kferu,khiggs)=gr(khiggs,kferd,kferu)
        gr(khiggs,kferu,kferd)=gr(khiggs,kferu,kferd)*dsqrt(fact)
        gr(kferu,khiggs,kferd)=gr(khiggs,kferu,kferd)
        gr(kferu,kferd,khiggs)=gr(khiggs,kferu,kferd)
      endif
c
c
c----------------------------------------------------------- into c s-bar
c i=2 (u index), j=2 (d index), (i-1)*3+j= 5
      kferu=kc
      kckmu=2
      kferd=ks
      kckmd=2
      gAusquared=1.d0/tanbe**2
      gAdsquared=tanbe**2
      mscale=mass(khiggs)
      if(mscale.gt.mass(kferu)+mass(kferd)+1.d0) then
        fact=DeltaQCD(mscale)
        hwidthqcd=3.d0/32.d0/pi*g2weak**2/(mass(kw))**2*mass(khiggs)
     &    *(ckm(kckmu,kckmd))**2
     &    *((dsrmq(mscale,kferu))**2*gAusquared
     &      +(mass(kferd))**2*gAdsquared)*fact
      else  
        hwidthqcd=0.d0
      endif
      mscale=2.d0*mass(lsp)
      fact=0.d0
      if(mscale.gt.mass(kferu)+mass(kferd)+1.d0) then
        fact=DeltaQCD(mscale)
      endif  
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(5,khico)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c     &      ,pname(kferd)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kferu)
     &      ,pname(kferd)
        write(*,*) 'qdc width = ',
     &    hwidthqcd
      endif  
      if(fact.gt.0.d0) then
c... redefine the couplings
c        write(*,*) 'vertex corr',dsqrt(fact)
        gl(khiggs,kferd,kferu)=gl(khiggs,kferd,kferu)*dsqrt(fact)
        gl(kferd,khiggs,kferu)=gl(khiggs,kferd,kferu)
        gl(kferd,kferu,khiggs)=gl(khiggs,kferd,kferu)
        gl(khiggs,kferu,kferd)=gl(khiggs,kferu,kferd)*dsqrt(fact)
        gl(kferu,khiggs,kferd)=gl(khiggs,kferu,kferd)
        gl(kferu,kferd,khiggs)=gl(khiggs,kferu,kferd)
        gr(khiggs,kferd,kferu)=gr(khiggs,kferd,kferu)*dsqrt(fact)
        gr(kferd,khiggs,kferu)=gr(khiggs,kferd,kferu)
        gr(kferd,kferu,khiggs)=gr(khiggs,kferd,kferu)
        gr(khiggs,kferu,kferd)=gr(khiggs,kferu,kferd)*dsqrt(fact)
        gr(kferu,khiggs,kferd)=gr(khiggs,kferu,kferd)
        gr(kferu,kferd,khiggs)=gr(khiggs,kferu,kferd)
      endif
c
c
c----------------------------------------------------------- into c b-bar
c i=2 (u index), j=3 (d index), (i-1)*3+j= 6
      kferu=kc
      kckmu=2
      kferd=kb
      kckmd=3
      gAusquared=1.d0/tanbe**2
      gAdsquared=tanbe**2
      mscale=mass(khiggs)
      if(mscale.gt.mass(kferu)+mass(kferd)+1.d0) then
        fact=DeltaQCD(mscale)
        hwidthqcd=3.d0/32.d0/pi*g2weak**2/(mass(kw))**2*mass(khiggs)
     &    *(ckm(kckmu,kckmd))**2
     &    *((dsrmq(mscale,kferu))**2*gAusquared
     &      +(dsrmq(mscale,kferd))**2*gAdsquared)*fact
      else  
        hwidthqcd=0.d0
      endif
      mscale=2.d0*mass(lsp)
      fact=0.d0
      if(mscale.gt.mass(kferu)+mass(kferd)+1.d0) then
        fact=DeltaQCD(mscale)
      endif  
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(6,khico)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kferu)
c     &      ,pname(kferd)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kferu)
     &      ,pname(kferd)
        write(*,*) 'qdc width = ',
     &    hwidthqcd
      endif  
      if(fact.gt.0.d0) then
c... redefine the couplings
c        write(*,*) 'vertex corr',dsqrt(fact)
        gl(khiggs,kferd,kferu)=gl(khiggs,kferd,kferu)*dsqrt(fact)
        gl(kferd,khiggs,kferu)=gl(khiggs,kferd,kferu)
        gl(kferd,kferu,khiggs)=gl(khiggs,kferd,kferu)
        gl(khiggs,kferu,kferd)=gl(khiggs,kferu,kferd)*dsqrt(fact)
        gl(kferu,khiggs,kferd)=gl(khiggs,kferu,kferd)
        gl(kferu,kferd,khiggs)=gl(khiggs,kferu,kferd)
        gr(khiggs,kferd,kferu)=gr(khiggs,kferd,kferu)*dsqrt(fact)
        gr(kferd,khiggs,kferu)=gr(khiggs,kferd,kferu)
        gr(kferd,kferu,khiggs)=gr(khiggs,kferd,kferu)
        gr(khiggs,kferu,kferd)=gr(khiggs,kferu,kferd)*dsqrt(fact)
        gr(kferu,khiggs,kferd)=gr(khiggs,kferu,kferd)
        gr(kferu,kferd,khiggs)=gr(khiggs,kferu,kferd)
      endif
c
c
c----------------------------------------------------------- into t d-bar
c i=3 (u index), j=1 (d index), (i-1)*3+j= 7
      kferu=kt
      kckmu=3
      kferd=kd
      kckmd=1
      gAusquared=1.d0/tanbe**2
      gAdsquared=tanbe**2
      mscale=mass(khiggs)
      if(mscale.gt.mass(kferu)+mass(kferd)+1.d0) then
        muu=(mass(kferu)/mscale)**2
        mud=(mass(kferd)/mscale)**2
        lambda=(1.d0-muu-mud)**2-4.d0*muu*mud
c... I am not sure you have to use dsralph3(mscale) in fact below (PU)
        fact=(1.d0-muu-mud)*((mass(kferu))**2*gAusquared
     &     *(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*Deltap(lambda,muu,mud))
     &     +(mass(kferd))**2*gAdsquared
     &     *(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*Deltap(lambda,mud,muu)))
     &    -4.d0*(mass(kferu))*(mass(kferd))*dsqrt(muu*mud)
     &     *(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*Deltam(lambda,muu,mud))
        hwidthqcd=3.d0/32.d0/pi*g2weak**2/(mass(kw))**2*mass(khiggs)
     &    *(ckm(kckmu,kckmd))**2*dsqrt(lambda)*fact
      else  
        hwidthqcd=0.d0
      endif
      mscale=2.d0*mass(lsp)
      factr=0.d0
      factl=0.d0
      if(mscale.gt.mass(kferu)+mass(kferd)+1.d0) then
        muu=(mass(kferu)/mscale)**2
        mud=(mass(kferd)/mscale)**2
        lambda=(1.d0-muu-mud)**2-4.d0*muu*mud
c... I am not sure you have to use dsralph3(mscale) in fact below (PU)
        factr=(mass(kferu))**2/(dsrmq(mscale,kferu))**2
     &     *(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*Deltap(lambda,muu,mud))
        factl=
     &     (1.d0+4.d0/3.d0*dsralph3(mscale)/pi*Deltap(lambda,mud,muu))
      endif  
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(7,khico)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kferu)
c     &      ,pname(kferd)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kferu)
     &      ,pname(kferd)
        write(*,*) 'qdc width = ',
     &    hwidthqcd
      endif  
      if(factr.gt.0.d0.and.factl.gt.0.d0) then
c... redefine the couplings
c        write(*,*) 'vertex corr',dsqrt(factr),dsqrt(factl) 
        gl(khiggs,kferd,kferu)=gl(khiggs,kferd,kferu)*dsqrt(factl)
        gl(kferd,khiggs,kferu)=gl(khiggs,kferd,kferu)
        gl(kferd,kferu,khiggs)=gl(khiggs,kferd,kferu)
        gl(khiggs,kferu,kferd)=gl(khiggs,kferu,kferd)*dsqrt(factl)
        gl(kferu,khiggs,kferd)=gl(khiggs,kferu,kferd)
        gl(kferu,kferd,khiggs)=gl(khiggs,kferu,kferd)
        gr(khiggs,kferd,kferu)=gr(khiggs,kferd,kferu)*dsqrt(factr)
        gr(kferd,khiggs,kferu)=gr(khiggs,kferd,kferu)
        gr(kferd,kferu,khiggs)=gr(khiggs,kferd,kferu)
        gr(khiggs,kferu,kferd)=gr(khiggs,kferu,kferd)*dsqrt(factr)
        gr(kferu,khiggs,kferd)=gr(khiggs,kferu,kferd)
        gr(kferu,kferd,khiggs)=gr(khiggs,kferu,kferd)
      endif
c
c
c----------------------------------------------------------- into t s-bar
c i=3 (u index), j=2 (d index), (i-1)*3+j= 8
      kferu=kt
      kckmu=3
      kferd=ks
      kckmd=2
      gAusquared=1.d0/tanbe**2
      gAdsquared=tanbe**2
      mscale=mass(khiggs)
      if(mscale.gt.mass(kferu)+mass(kferd)+1.d0) then
        muu=(mass(kferu)/mscale)**2
        mud=(mass(kferd)/mscale)**2
        lambda=(1.d0-muu-mud)**2-4.d0*muu*mud
c... I am not sure you have to use dsralph3(mscale) in fact below (PU)
        fact=(1.d0-muu-mud)*((mass(kferu))**2*gAusquared
     &     *(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*Deltap(lambda,muu,mud))
     &     +(mass(kferd))**2*gAdsquared
     &     *(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*Deltap(lambda,mud,muu)))
     &    -4.d0*(mass(kferu))*(mass(kferd))*dsqrt(muu*mud)
     &     *(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*Deltam(lambda,muu,mud))
        hwidthqcd=3.d0/32.d0/pi*g2weak**2/(mass(kw))**2*mass(khiggs)
     &    *(ckm(kckmu,kckmd))**2*dsqrt(lambda)*fact
      else  
        hwidthqcd=0.d0
      endif
      mscale=2.d0*mass(lsp)
      factr=0.d0
      factl=0.d0
      if(mscale.gt.mass(kferu)+mass(kferd)+1.d0) then
        muu=(mass(kferu)/mscale)**2
        mud=(mass(kferd)/mscale)**2
        lambda=(1.d0-muu-mud)**2-4.d0*muu*mud
c... I am not sure you have to use dsralph3(mscale) in fact below (PU)
        factr=(mass(kferu))**2/(dsrmq(mscale,kferu))**2
     &     *(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*Deltap(lambda,muu,mud))
        factl=
     &     (1.d0+4.d0/3.d0*dsralph3(mscale)/pi*Deltap(lambda,mud,muu))
      endif  
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(8,khico)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kferu)
c     &      ,pname(kferd)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kferu)
     &      ,pname(kferd)
        write(*,*) 'qdc width = ',
     &    hwidthqcd
      endif  
      if(factr.gt.0.d0.and.factl.gt.0.d0) then
c... redefine the couplings
c        write(*,*) 'vertex corr',dsqrt(factr),dsqrt(factl) 
        gl(khiggs,kferd,kferu)=gl(khiggs,kferd,kferu)*dsqrt(factl)
        gl(kferd,khiggs,kferu)=gl(khiggs,kferd,kferu)
        gl(kferd,kferu,khiggs)=gl(khiggs,kferd,kferu)
        gl(khiggs,kferu,kferd)=gl(khiggs,kferu,kferd)*dsqrt(factl)
        gl(kferu,khiggs,kferd)=gl(khiggs,kferu,kferd)
        gl(kferu,kferd,khiggs)=gl(khiggs,kferu,kferd)
        gr(khiggs,kferd,kferu)=gr(khiggs,kferd,kferu)*dsqrt(factr)
        gr(kferd,khiggs,kferu)=gr(khiggs,kferd,kferu)
        gr(kferd,kferu,khiggs)=gr(khiggs,kferd,kferu)
        gr(khiggs,kferu,kferd)=gr(khiggs,kferu,kferd)*dsqrt(factr)
        gr(kferu,khiggs,kferd)=gr(khiggs,kferu,kferd)
        gr(kferu,kferd,khiggs)=gr(khiggs,kferu,kferd)
      endif
c
c
c----------------------------------------------------------- into t b-bar
c i=3 (u index), j=3 (d index), (i-1)*3+j= 9
      kferu=kt
      kckmu=3
      kferd=kb
      kckmd=3
      gAusquared=1.d0/tanbe**2
      gAdsquared=tanbe**2
      mscale=mass(khiggs)
      if(mscale.gt.mass(kferu)+mass(kferd)+1.d0) then
        muu=(mass(kferu)/mscale)**2
        mud=(mass(kferd)/mscale)**2
        lambda=(1.d0-muu-mud)**2-4.d0*muu*mud
c... I am not sure you have to use dsralph3(mscale) in fact below (PU)
        fact=(1.d0-muu-mud)*((mass(kferu))**2*gAusquared
     &     *(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*Deltap(lambda,muu,mud))
     &     +(mass(kferd))**2*gAdsquared
     &     *(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*Deltap(lambda,mud,muu)))
     &    -4.d0*(mass(kferu))*(mass(kferd))*dsqrt(muu*mud)
     &     *(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*Deltam(lambda,muu,mud))
        hwidthqcd=3.d0/32.d0/pi*g2weak**2/(mass(kw))**2*mass(khiggs)
     &    *(ckm(kckmu,kckmd))**2*dsqrt(lambda)*fact
      else  
        hwidthqcd=0.d0
      endif
      mscale=2.d0*mass(lsp)
      factr=0.d0
      factl=0.d0
      if(mscale.gt.mass(kferu)+mass(kferd)+1.d0) then
        muu=(mass(kferu)/mscale)**2
        mud=(mass(kferd)/mscale)**2
        lambda=(1.d0-muu-mud)**2-4.d0*muu*mud
c... I am not sure you have to use dsralph3(mscale) in fact below (PU)
        factr=(mass(kferu))**2/(dsrmq(mscale,kferu))**2
     &     *(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*Deltap(lambda,muu,mud))
        factl=(mass(kferd))**2/(dsrmq(mscale,kferd))**2
     &     *(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*Deltap(lambda,mud,muu))
      endif  
      if(hwidthqcd.ge.0.d0) then
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(9,khico)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kferu)
c     &      ,pname(kferd)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigwid for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kferu)
     &      ,pname(kferd)
        write(*,*) 'qdc width = ',
     &    hwidthqcd
      endif  
      if(factr.gt.0.d0.and.factl.gt.0.d0) then
c        write(*,*) 'vertex corr',dsqrt(factr),dsqrt(factl) 
c... redefine the couplings
        gl(khiggs,kferd,kferu)=gl(khiggs,kferd,kferu)*dsqrt(factl)
        gl(kferd,khiggs,kferu)=gl(khiggs,kferd,kferu)
        gl(kferd,kferu,khiggs)=gl(khiggs,kferd,kferu)
        gl(khiggs,kferu,kferd)=gl(khiggs,kferu,kferd)*dsqrt(factl)
        gl(kferu,khiggs,kferd)=gl(khiggs,kferu,kferd)
        gl(kferu,kferd,khiggs)=gl(khiggs,kferu,kferd)
        gr(khiggs,kferd,kferu)=gr(khiggs,kferd,kferu)*dsqrt(factr)
        gr(kferd,khiggs,kferu)=gr(khiggs,kferd,kferu)
        gr(kferd,kferu,khiggs)=gr(khiggs,kferd,kferu)
        gr(khiggs,kferu,kferd)=gr(khiggs,kferu,kferd)*dsqrt(factr)
        gr(kferu,khiggs,kferd)=gr(khiggs,kferu,kferd)
        gr(kferu,kferd,khiggs)=gr(khiggs,kferu,kferd)
      endif
c
c
c---------------------------------------------------------- into leptons
      do i=1,3
        vvv(i)=0.0d0
        k1=knu(i)
        k2=kl(i)
        mi=mass(k1)
        mj=mass(k2)
        if ((mi+mj).lt.mh) then
          sqlam=sqrt((mi**2+mj**2-mh**2)**2-4*mi**2*mj**2)
          partial=sqlam/(16.0d0*pi*mh**3)*
     &      ((mh**2-mi**2-mj**2)*
     &      (dsabsq(gr(khiggs,k2,k1))+dsabsq(gl(khiggs,k2,k1)))-
     &      4.0d0*mi*mj*gl(khiggs,k2,k1)*gr(khiggs,k2,k1))
          width(khiggs)=width(khiggs)+partial
          vvv(i)=partial
c          write(*,*) pname(k1),pname(k2),partial
        endif
      enddo
      do i=1,3
        hdwidth(i+9,khico)=vvv(i)
      enddo

c------------------------------------------------------------ into w+ h1
      k1=kw
      k2=kh1
      mi=mass(k1)
      mj=mass(k2)
      if ((mi+mj).lt.mh) then
        sqlam=sqrt((mi**2+mj**2-mh**2)**2-4*mi**2*mj**2)
        partial=sqlam**3*dsabsq(gl(kw,khiggs,kh1))/
     &    (16.0d0*pi*mass(kw)**2*mh**3)
        width(khiggs)=width(khiggs)+partial
c        write(*,*) pname(k1),pname(k2),partial
        hdwidth(13,khico)=partial
      endif

c------------------------------------------------------------ into w+ h2
      k1=kw
      k2=kh2
      mi=mass(k1)
      mj=mass(k2)
      if ((mi+mj).lt.mh) then
        sqlam=sqrt((mi**2+mj**2-mh**2)**2-4*mi**2*mj**2)
        partial=sqlam**3*dsabsq(gl(kw,khiggs,kh2))/
     &    (16.0d0*pi*mass(kw)**2*mh**3)
        width(khiggs)=width(khiggs)+partial
c        write(*,*) pname(k1),pname(k2),partial
        hdwidth(14,khico)=partial
      endif

c------------------------------------------------------------ into w+ h3
      k1=kw
      k2=kh3
      mi=mass(k1)
      mj=mass(k2)
      if ((mi+mj).lt.mh) then
        sqlam=sqrt((mi**2+mj**2-mh**2)**2-4*mi**2*mj**2)
        partial=sqlam**3*dsabsq(gl(kw,khiggs,kh3))/
     &    (16.0d0*pi*mass(kw)**2*mh**3)
        width(khiggs)=width(khiggs)+partial
c        write(*,*) pname(k1),pname(k2),partial
        hdwidth(15,khico)=partial
      endif

c---------------------------------------- into charginos and neutralinos
      parsum=0.d0
      do i=1,2
        do j=1,4
          k1=kcha(i)
          k2=kn(j)
          mi=mass(k1)
          mj=mass(k2)
          if ((mi+mj).lt.mh) then
            sqlam=sqrt((mi**2+mj**2-mh**2)**2-4*mi**2*mj**2)
            partial=sqlam/(16.0d0*pi*mh**3)*
     &        ((mh**2-mi**2-mj**2)*
     &        (dsabsq(gl(khiggs,k2,k1))+dsabsq(gr(khiggs,k2,k1)))-
     &        4.0d0*mi*mj*gl(khiggs,k2,k1)*gr(khiggs,k2,k1))
            parsum=parsum+partial
c            write(*,*) pname(k1),pname(k2),partial
          endif
        enddo
      enddo
      width(khiggs)=width(khiggs)+parsum
      hdwidth(21,khico)=parsum
c
c
c---------------------------------- into up-sfermion - down-sfermion star
      parsum=0.d0
c sneutrino, slepton-star:
c
      do i=1,3
      ksfer1=ksnu(i)
      ksfer2=ksl(i)
      if (mass(khiggs).ge.mass(ksfer1)+mass(ksfer2)) then
         kf = dsqrt(0.25d0*(1.d0-(mass(ksfer1)/mass(khiggs))**2
     &      -(mass(ksfer2)/mass(khiggs))**2)**2
     &      -(mass(ksfer1)*mass(ksfer2))**2/(mass(khiggs))**4)
         partial = 3.d0*(1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer1,ksfer2)) / (8.d0*pi)
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer1),pname(ksfer2),partial
      endif
      ksfer2=ksl(i+3)
      if (mass(khiggs).ge.mass(ksfer1)+mass(ksfer2)) then
         kf = dsqrt(0.25d0*(1.d0-(mass(ksfer1)/mass(khiggs))**2
     &      -(mass(ksfer2)/mass(khiggs))**2)**2
     &      -(mass(ksfer1)*mass(ksfer2))**2/(mass(khiggs))**4)
         partial = 3.d0*(1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer1,ksfer2)) / (8.d0*pi)
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer1),pname(ksfer2),partial
      endif
      enddo
c up-squark, down-squark-star:
c
      do i=1,6
      ksfer1=ksqu(i)
      ksfer2=ksqd(i)
      if (mass(khiggs).ge.mass(ksfer1)+mass(ksfer2)) then
         kf = dsqrt(0.25d0*(1.d0-(mass(ksfer1)/mass(khiggs))**2
     &      -(mass(ksfer2)/mass(khiggs))**2)**2
     &      -(mass(ksfer1)*mass(ksfer2))**2/(mass(khiggs))**4)
         partial = 3.d0*(1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer1,ksfer2)) / (8.d0*pi)
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer1),pname(ksfer2),partial
      endif
      if(i.le.3) then
        ksfer2=ksqd(i+3)
      else
        ksfer2=ksqd(i-3)
      endif  
      if (mass(khiggs).ge.mass(ksfer1)+mass(ksfer2)) then
         kf = dsqrt(0.25d0*(1.d0-(mass(ksfer1)/mass(khiggs))**2
     &      -(mass(ksfer2)/mass(khiggs))**2)**2
     &      -(mass(ksfer1)*mass(ksfer2))**2/(mass(khiggs))**4)
         partial = 3.d0*(1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer1,ksfer2)) / (8.d0*pi)
c... sum the new term
         parsum=parsum+partial
c        write(*,*) pname(ksfer1),pname(ksfer2),partial
      endif
      enddo
      width(khiggs)=width(khiggs)+parsum
      hdwidth(20,khico)=parsum

c=============================================== decay widths calculated

      if (prtlevel.eq.-1) then
        write(*,*) 'suhwidths:'
        write(*,*) '  mh(1) = ',mass(kh1),'  width(kh1) = ',width(kh1)
        write(*,*) '  mh(2) = ',mass(kh2),'  width(kh2) = ',width(kh2)
        write(*,*) '  mh(3) = ',mass(kh3),'  width(kh3) = ',width(kh3)
        write(*,*) '    mhc = ',mass(khc),'  width(khc) = ',width(khc)

        do i=1,4
          do j=1,32
            write(*,*) 'hw: i=',i,' j=',j,': ',hdwidth(j,i)
          enddo
        enddo

      endif

      end


cccccccccccccccccccccc

      real*8 function DeltaQCD(mscale)
      implicit none
      include 'dsmssm.h'
      real*8 mscale,dsralph3
      integer nf
      real*8 pi
      parameter (pi=3.141592653589793238d0)
c
      if(mscale.lt.mtmt) then
        nf=5
      else
        nf=6
      endif  
      DeltaQCD=1.d0+5.67d0*dsralph3(mscale)/pi
     &  +(35.94d0-1.36d0*nf)*(dsralph3(mscale)/pi)**2
      return
      end

cccccccccccccccccccccc

      real*8 function Deltat(mscale,khiggs,kfer)
      implicit none
      include 'dsmssm.h'
      real*8 mscale,dsralph3,dsrmq
      integer khiggs,kfer
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      if(khiggs.eq.kh1.or.khiggs.eq.kh2) then
        Deltat=(dsralph3(mscale)/pi)**2*
     &    (1.57d0-2.d0/3.d0*dlog((mass(khiggs)/mass(kt))**2)+
     &     1.d0/9.d0*(dlog((dsrmq(mscale,kfer)/mass(khiggs))**2))**2)
      elseif(khiggs.eq.kh3) then  
        Deltat=(dsralph3(mscale)/pi)**2*
     &    (3.83d0-dlog((mass(khiggs)/mass(kt))**2)+
     &     1.d0/6.d0*(dlog((dsrmq(mscale,kfer)/mass(khiggs))**2))**2)
      else  
        write(*,*) 'Deltat called for wrong higgs = ',pname(khiggs)
        write(*,*) 'program stopped'
        stop
      endif  
      return
      end

cccccccccccccccccccccc 

      real*8 function DeltaPhi(kf)
      implicit none
      real*8 kf,kffrac,A,dsdilog2
      kffrac=(1.d0-kf)/(1.d0+kf)
      A=(1.d0+kf**2)*(4.d0*dsdilog2(kffrac)+2.d0*dsdilog2(-kffrac)
     &   -3.d0*dlog(1.d0/kffrac)*dlog(2.d0/(1.d0+kf))
     &   -2.d0*dlog(1.d0/kffrac)*dlog(kf))
     &   -3.d0*kf*dlog(4.d0/(1.d0-kf**2))-4.d0*kf*dlog(kf)
c
      DeltaPhi=1.d0/kf*A
     &   +1.d0/16.d0/kf**3*(3.d0+34.d0*kf**2-13.d0*kf**4)
     &   *dlog(1.d0/kffrac)+3.d0/8.d0/kf**2*(7.d0*kf**2-1.d0)
      return
      end

cccccccccccccccccccccc

      real*8 function DeltaA(kf)
      implicit none
      real*8 kf,kffrac,A,dsdilog2
      kffrac=(1.d0-kf)/(1.d0+kf)
      A=(1.d0+kf**2)*(4.d0*dsdilog2(kffrac)+2.d0*dsdilog2(-kffrac)
     &   -3.d0*dlog(1.d0/kffrac)*dlog(2.d0/(1.d0+kf))
     &   -2.d0*dlog(1.d0/kffrac)*dlog(kf))
     &   -3.d0*kf*dlog(4.d0/(1.d0-kf**2))-4.d0*kf*dlog(kf)
c
      DeltaA=1.d0/kf*A
     &   +1.d0/16.d0/kf*(19.d0+2.d0*kf**2+3.d0*kf**4)
     &   *dlog(1.d0/kffrac)+3.d0/8.d0*(7.d0-kf**2)
      return
      end

cccccccccccccccccccccc

      real*8 function Deltap(lambda,mui,muj)
      implicit none
      real*8 lambda,mui,muj,xi,xj,BigBij
c
      xi=2.d0*mui/(1.d0-mui-muj+dsqrt(lambda))
      xj=2.d0*muj/(1.d0-mui-muj+dsqrt(lambda))
      Deltap=9.d0/4.d0+(3.d0-2.d0*mui+2.d0*muj)/4.d0*dlog(mui/muj)
     &  +((1.5d0-mui-muj)*lambda+5.d0*mui*muj)/(2.d0*dsqrt(lambda)
     &    *(1.d0-mui-muj))*dlog(xi*xj)+BigBij(lambda,mui,muj,xi,xj)
      return
      end

cccccccccccccccccccccc

      real*8 function Deltam(lambda,mui,muj)
      implicit none
      real*8 lambda,mui,muj,xi,xj,BigBij
c
      xi=2.d0*mui/(1.d0-mui-muj+dsqrt(lambda))
      xj=2.d0*muj/(1.d0-mui-muj+dsqrt(lambda))
      Deltam=3.d0+(muj-mui)/2.d0*dlog(mui/muj)
     &  +(lambda+2.d0*(1.d0-mui-muj))/(2.d0*dsqrt(lambda))*dlog(xi*xj)
     &  +BigBij(lambda,mui,muj,xi,xj)
      return
      end

cccccccccccccccccccccc

      real*8 function BigBij(lambda,mui,muj,xi,xj)
      implicit none
      real*8 lambda,mui,muj,xi,xj,dsdilog2
c
c      Bij=(1.d0-mui-muj)/dsqrt(lambda)*(4.d0*dsdilog2(xi*xj)
c     &   -2.d0*dsdilog2(xi)-2.d0*dsdilog2(xj)
c     &   +2.d0*dlog(xi*xj)*dlog(1.d0-xi*xj)-dlog(xi)*dlog(1.d0-xi)
c     &   -dlog(xj)*dlog(1.d0-xj))
c     &   -4.d0*(dlog(1.d0-xi*xj)+xi*xj/(1.d0-xi*xj)*dlog(xi*xj))
c     &   +(dsqrt(lambda)+mui-muj)/dsqrt(lambda)
c     &     *(dlog(1.d0-xi)+xi/(1.d0-xi)*dlog(xi))
c     &   +(dsqrt(lambda)-mui+muj)/dsqrt(lambda)
c     &     *(dlog(1.d0-xj)+xj/(1.d0-xj)*dlog(xj))
c
c Spira review:
      BigBij=(1.d0-mui-muj)/dsqrt(lambda)*(4.d0*dsdilog2(xi*xj)
     &   -2.d0*dsdilog2(-xi)-2.d0*dsdilog2(-xj)
     &   +2.d0*dlog(xi*xj)*dlog(1.d0-xi*xj)-dlog(xi)*dlog(1.d0+xi)
     &   -dlog(xj)*dlog(1.d0+xj))
     &   -4.d0*(dlog(1.d0-xi*xj)+xi*xj/(1.d0-xi*xj)*dlog(xi*xj))
     &   +(dsqrt(lambda)+mui-muj)/dsqrt(lambda)
     &     *(dlog(1.d0+xi)-xi/(1.d0+xi)*dlog(xi))
     &   +(dsqrt(lambda)-mui+muj)/dsqrt(lambda)
     &     *(dlog(1.d0+xj)-xj/(1.d0+xj)*dlog(xj))
      return
      end

 
cccccccccccccccccccccc

      complex*16 function AQhH(tau)
      implicit none
      real*8 tau
      complex*16 ftau 
c
      AQhH=1.5d0*tau*(1.d0+(1.d0-tau)*ftau(tau))
      return
      end

cccccccccccccccccccccc

      complex*16 function AsQhH(tau)
      implicit none
      real*8 tau
      complex*16 ftau 
c
      AsQhH=-0.75d0*tau*(1.d0-tau*ftau(tau))
      return
      end

cccccccccccccccccccccc

      complex*16 function AQA(tau)
      implicit none
      real*8 tau
      complex*16 ftau 
c
      AQA=tau*ftau(tau)
      return
      end

cccccccccccccccccccccc

      complex*16 function ftau(tau)
      implicit none
      real*8 tau,frac
      real*8 pi
      parameter (pi=3.141592653589793238d0)
c
      if(tau.ge.1.d0) then
        ftau=(dasin(1.d0/dsqrt(tau)))**2
      else
        frac=(1.d0+dsqrt(1.d0-tau))/(1.d0-dsqrt(1.d0-tau)) 
        ftau=-0.25d0*((dlog(frac))**2-pi**2
     &                 -2.d0*dlog(frac)*dcmplx(0.d0,1.d0))
      endif  
      return
      end


cccccccccccccccccccccc

      real*8 function dsdilog2(x)
      implicit none
      real*8 x,dsdilog,t
      real*8 pi
      parameter (pi=3.141592653589793238d0)
c
      if(x.le.1.d0.and.x.ge.-1.d0) then
        dsdilog2=dsdilog(x)
      else
        t=1.d0/x
        dsdilog2=-dsdilog(t)-1.d0/6.d0*pi**2
        if(t.lt.0.d0) then
          dsdilog2=dsdilog2-0.5d0*(dlog(-t))**2
        else
          dsdilog2=dsdilog2-0.5d0*((dlog(t))**2-pi**2)
        endif  
      endif  
      return
      end

