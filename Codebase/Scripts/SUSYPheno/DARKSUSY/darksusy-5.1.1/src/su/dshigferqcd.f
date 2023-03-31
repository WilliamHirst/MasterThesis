      subroutine dshigferqcd
c_______________________________________________________________________
c  THIS ROUTINE IS OBSOLETE. USE DSHIGWID INSTEAD.
c  qcd corrections to the widths of the decays H^0_i --> c cbar, b bbar,
c  t tbar and to the corresponding vertices
c
c  typed in from formulas in Djouadi, Spira and Zerwas, astro-ph/9511344
c  and Spira, astro-ph/9705337
c
c  author: piero ullio, ullio@sissa.it, 02-09-13
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsidtag.h'
      real*8 hwidthnoqcd,hwidthqcd,gratio,gsquared,mscale,fact
      real*8 dsralph3,dsrmq
      real*8 DeltaQCD,Deltat,DeltaPhi,DeltaA,kf
      integer kfer,khiggs,ksfer,ksfer1,ksfer2,i
      integer kferu,kferd,kckmu,kckmd
      real*8 gAusquared,gAdsquared,muu,mud,lambda,Deltap,Deltam,
     & factr,factl,dsabsq,partial
      real*8 gcou,gcod,tau
      complex*16 ampl,AQhH,AsQhH,AQA
      integer nf
      real*8 BigEnf,hwidthqcdnl
      complex*16 sqampl,parampl
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      logical first
      data first/.true./
      save first

      if (first) then
        write(*,*) ' '
        write(*,*)
     &    'WARNING in dshigferqcd: This routine is obsolete. ',
     &    'Use dshigwid instead.'
        write(*,*) 'This warning message will only be printed once.'
        first=.false.
      endif

c

c======================================================== h2 decay width
      khiggs=kh2
c----------------------------------------------------------- into c c-bar
      hwidthnoqcd=hdwidth(22,2)
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
c... subtract out the contribution in dshwidth
        width(khiggs)=width(khiggs)-hwidthnoqcd
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(22,2)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'no qcd and with qdc width = ',hwidthnoqcd,
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
      hwidthnoqcd=hdwidth(25,2)
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
c... subtract out the contribution in dshwidth
        width(khiggs)=width(khiggs)-hwidthnoqcd
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(25,2)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'no qcd and with qdc width = ',hwidthnoqcd,
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
      hwidthnoqcd=hdwidth(24,2)
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
c... subtract out the contribution in dshwidth
        width(khiggs)=width(khiggs)-hwidthnoqcd
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(24,2)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'no qcd and with qdc width = ',hwidthnoqcd,
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
        hdwidth(26,2)=hwidthqcd
c        write(*,*) 'Gluons in final state'
c        write(*,*) 'qdc width = ',hwidthqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
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
        hdwidth(26,2)=hdwidth(26,2)+hwidthqcdnl
c        write(*,*) 'Gluons in final state'
c        write(*,*) 'next to leading qdc term = ',hwidthqcdnl
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Gluons in final state'
        write(*,*) 'next to leading qdc term = ',hwidthqcdnl
      endif
c
c
c------------------------------------------ into sfermion - sfermion star
      hdwidth(30,2)=0.d0
c sneutrinos:
c
      do i=1,3
      ksfer=ksnu(i)  
      if (mass(khiggs).ge.2.d0*mass(ksfer)) then
         kf = sqrt(0.25d0-(mass(ksfer)/mass(khiggs))**2)
         partial = (1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer,ksfer)) / (8.d0*pi)
c... sum the new term
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,2)=hdwidth(30,2)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,2)=hdwidth(30,2)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,2)=hdwidth(30,2)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,2)=hdwidth(30,2)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,2)=hdwidth(30,2)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,2)=hdwidth(30,2)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,2)=hdwidth(30,2)+partial
c        write(*,*) pname(ksfer1),pname(ksfer2),partial
      endif
      enddo
c
c
c======================================================== h1 decay width
      khiggs=kh1
c----------------------------------------------------------- into c c-bar
      hwidthnoqcd=hdwidth(22,1)
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
c... subtract out the contribution in dshwidth
        width(khiggs)=width(khiggs)-hwidthnoqcd
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(22,1)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'no qcd and with qdc width = ',hwidthnoqcd,
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
      hwidthnoqcd=hdwidth(25,1)
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
c... subtract out the contribution in dshwidth
        width(khiggs)=width(khiggs)-hwidthnoqcd
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(25,1)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'no qcd and with qdc width = ',hwidthnoqcd,
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
      hwidthnoqcd=hdwidth(24,1)
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
     &     *(1.d0+4.d0/3.d0*dsralph3(mscale)/pi*DeltaPhi(kf))
      endif
      if(hwidthqcd.ge.0.d0) then
c... subtract out the contribution in dshwidth
        width(khiggs)=width(khiggs)-hwidthnoqcd
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(24,1)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'no qcd and with qdc width = ',hwidthnoqcd,
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
        hdwidth(26,1)=hwidthqcd
c        write(*,*) 'Gluons in final state'
c        write(*,*) 'qdc width = ',hwidthqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
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
        hdwidth(26,1)=hdwidth(26,1)+hwidthqcdnl
c        write(*,*) 'Gluons in final state'
c        write(*,*) 'next to leading qdc term = ',hwidthqcdnl
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Gluons in final state'
        write(*,*) 'next to leading qdc term = ',hwidthqcdnl
      endif
c
c------------------------------------------ into sfermion - sfermion star
      hdwidth(30,1)=0.d0
c sneutrinos:
c
      do i=1,3
      ksfer=ksnu(i)  
      if (mass(khiggs).ge.2.d0*mass(ksfer)) then
         kf = sqrt(0.25d0-(mass(ksfer)/mass(khiggs))**2)
         partial = (1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer,ksfer)) / (8.d0*pi)
c... sum the new term
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,1)=hdwidth(30,1)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,1)=hdwidth(30,1)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,1)=hdwidth(30,1)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,1)=hdwidth(30,1)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,1)=hdwidth(30,1)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,1)=hdwidth(30,1)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,1)=hdwidth(30,1)+partial
c        write(*,*) pname(ksfer1),pname(ksfer2),partial
      endif
      enddo
c
c
c======================================================== h3 decay width
      khiggs=kh3
c----------------------------------------------------------- into c c-bar
      hwidthnoqcd=hdwidth(22,3)
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
c... subtract out the contribution in dshwidth
        width(khiggs)=width(khiggs)-hwidthnoqcd
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(22,3)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'no qcd and with qdc width = ',hwidthnoqcd,
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
      hwidthnoqcd=hdwidth(25,3)
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
c... subtract out the contribution in dshwidth
        width(khiggs)=width(khiggs)-hwidthnoqcd
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(25,3)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'no qcd and with qdc width = ',hwidthnoqcd,
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
      hwidthnoqcd=hdwidth(24,3)
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
c... subtract out the contribution in dshwidth
        width(khiggs)=width(khiggs)-hwidthnoqcd
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(24,3)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kfer)
        write(*,*) 'no qcd and with qdc width = ',hwidthnoqcd,
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
        hdwidth(26,3)=hwidthqcd
c        write(*,*) 'Gluons in final state'
c        write(*,*) 'qdc width = ',hwidthqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
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
        hdwidth(26,3)=hdwidth(26,3)+hwidthqcdnl
c        write(*,*) 'Gluons in final state'
c        write(*,*) 'next to leading qdc term = ',hwidthqcdnl
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Gluons in final state'
        write(*,*) 'next to leading qdc term = ',hwidthqcdnl
      endif
c
c
c------------------------------------------ into sfermion - sfermion star
      hdwidth(30,3)=0.d0
c sleptons, i-i final state:
c
      do i=1,6
      ksfer=ksl(i)
      if (mass(khiggs).ge.2.d0*mass(ksfer)) then
         kf = sqrt(0.25d0-(mass(ksfer)/mass(khiggs))**2)
         partial = (1.d0/mass(khiggs)) * kf *
     &        dsabsq(gr(khiggs,ksfer,ksfer)) / (8.d0*pi)
c... sum the new term
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,3)=hdwidth(30,3)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,3)=hdwidth(30,3)+partial
c        write(*,*) pname(ksfer),pname(ksfer),partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,3)=hdwidth(30,3)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,3)=hdwidth(30,3)+partial
c        write(*,*) pname(ksfer),pname(ksfer),partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,3)=hdwidth(30,3)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(30,3)=hdwidth(30,3)+partial
c        write(*,*) pname(ksfer),pname(ksfer),partial
      endif
      enddo
c
c
c======================================================== h+ decay width
      khiggs=khc
c----------------------------------------------------------- into u b-bar
c i=1 (u index), j=3 (d index), (i-1)*3+j= 3
      hwidthnoqcd=hdwidth(3,4)
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
c... subtract out the contribution in dshwidth
        width(khiggs)=width(khiggs)-hwidthnoqcd
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(3,4)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kferu)
c     &      ,pname(kferd)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kferu)
     &      ,pname(kferd)
        write(*,*) 'no qcd and with qdc width = ',hwidthnoqcd,
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
      hwidthnoqcd=hdwidth(4,4)
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
c... subtract out the contribution in dshwidth
        width(khiggs)=width(khiggs)-hwidthnoqcd
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(4,4)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kferu)
c     &      ,pname(kferd)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kferu)
     &      ,pname(kferd)
        write(*,*) 'no qcd and with qdc width = ',hwidthnoqcd,
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
      hwidthnoqcd=hdwidth(5,4)
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
c... subtract out the contribution in dshwidth
        width(khiggs)=width(khiggs)-hwidthnoqcd
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(5,4)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kfer)
c     &      ,pname(kferd)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kferu)
     &      ,pname(kferd)
        write(*,*) 'no qcd and with qdc width = ',hwidthnoqcd,
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
      hwidthnoqcd=hdwidth(6,4)
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
c... subtract out the contribution in dshwidth
        width(khiggs)=width(khiggs)-hwidthnoqcd
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(6,4)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kferu)
c     &      ,pname(kferd)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kferu)
     &      ,pname(kferd)
        write(*,*) 'no qcd and with qdc width = ',hwidthnoqcd,
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
      hwidthnoqcd=hdwidth(7,4)
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
c... subtract out the contribution in dshwidth
        width(khiggs)=width(khiggs)-hwidthnoqcd
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(7,4)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kferu)
c     &      ,pname(kferd)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kferu)
     &      ,pname(kferd)
        write(*,*) 'no qcd and with qdc width = ',hwidthnoqcd,
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
      hwidthnoqcd=hdwidth(8,4)
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
c... subtract out the contribution in dshwidth
        width(khiggs)=width(khiggs)-hwidthnoqcd
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(8,4)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kferu)
c     &      ,pname(kferd)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kferu)
     &      ,pname(kferd)
        write(*,*) 'no qcd and with qdc width = ',hwidthnoqcd,
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
      hwidthnoqcd=hdwidth(9,4)
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
c... subtract out the contribution in dshwidth
        width(khiggs)=width(khiggs)-hwidthnoqcd
c... sum the new term
        width(khiggs)=width(khiggs)+hwidthqcd
        hdwidth(9,4)=hwidthqcd
c        write(*,*) 'higgs, fer: ',pname(khiggs),pname(kferu)
c     &      ,pname(kferd)
c        write(*,*) 'qcd, noqcd = ',hwidthqcd,hwidthnoqcd
      else
        write(*,*) 'troubles in dshigferqcd for model ',idtag
        write(*,*) 'Higgs particle = ',pname(khiggs)
        write(*,*) 'Fermion in final state = ',pname(kferu)
     &      ,pname(kferd)
        write(*,*) 'no qcd and with qdc width = ',hwidthnoqcd,
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
c---------------------------------- into up-sfermion - down-sfermion star
      hdwidth(20,4)=0.d0
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(20,4)=hdwidth(20,4)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(20,4)=hdwidth(20,4)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(20,4)=hdwidth(20,4)+partial
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
         width(khiggs)=width(khiggs)+partial
         hdwidth(20,4)=hdwidth(20,4)+partial
c        write(*,*) pname(ksfer1),pname(ksfer2),partial
      endif
      enddo
c
c
      return
      end


