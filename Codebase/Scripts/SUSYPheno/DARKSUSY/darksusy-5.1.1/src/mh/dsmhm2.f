***********************************************************************
*** dsmhm2 returns the full scattering amplitude at zero 
*** momentum transfer squared, averaged over the initial spin states of 
*** the DM particle and summed over all other spin states.
***
***  input: m0     - DM mass (in GeV)
***         SMtype - SM scattering partners
***             "    = 7,8,9 - u,d,s quarks
***                    4,5,6 - e,m,t leptons
***                    1,2,3 - e,m,t neutrinos
***
*** author: torsten bringmann (troms@physto.se), 2010-01-23
***********************************************************************

      real*8 function dsmhm2(omega,m0,SMtype)
      implicit none

      include 'dsmssm.h'
      include 'dsmhcom.h'
      include 'dsidtag.h'

      real*8  m0,omega,omegax,xf,tmpres
      integer SMtype,kf,ksf(2),nsf,i,j

c... propagators and various coupling combinations
      complex*8 ds(2),du(2),dtZ,dtH(2)
      complex*8 omgllrr(2,2),mfglrrl(2,2)
      real*8 gZZtot,gHH(2),gimlr(2),grelr(2),glrsq(2)
      real*8 gZsttot(2),gZtutot(2)
      real*8 mss,muu,mtt,msu,mst,mtu

      character*12 memory    ! to suppress multiple error messages
      save memory
      data memory /'____________'/


      dsmhm2=0d0
      tmpres=0d0
      

c... setup (s)particles
      if (SMtype.le.3) kf=knu(SMtype)
      if (SMtype.ge.4.and.SMtype.le.6) kf=kl(SMtype-3)
      if (SMtype.eq.7) kf=ku
      if (SMtype.eq.8) kf=kd
      if (SMtype.eq.9) kf=ks
      if (SMtype.ge.10) return
      nsf=2
      if (SMtype.le.3) nsf=1
        do 10 i=1,nsf
          if (SMtype.le.3) ksf(i)=ksnu(SMtype)
          if (SMtype.ge.4.and.SMtype.le.6)
     -       ksf(i)=ksl(SMtype+3*i-6)
          if (SMtype.eq.7) ksf(i)=ksu(i)
          if (SMtype.eq.8) ksf(i)=ksd(i)
          if (SMtype.eq.9) ksf(i)=kss(i)
  10    continue

c... from here on, keep everything dimensionless
      omegax=omega/m0
      xf=mass(kf)/m0


      if (mhtype.eq.1) then ! neutralino scattering

c... set up propagators and couplings
        do 20 i=1,nsf
          ds(i)=dcmplx(1.+2*omegax+xf**2-(mass(ksf(i))/m0)**2,
     -           width(ksf(i))/m0)
          du(i)=dcmplx(1.-2*omegax+xf**2-(mass(ksf(i))/m0)**2,
     -           width(ksf(i))/m0)
  20    continue
        dtZ=dcmplx(-(mass(kz)/m0)**2,width(kz)/m0)
        dtH(1)=dcmplx(-(mass(kh1)/m0)**2,width(kh1)/m0)
        dtH(2)=dcmplx(-(mass(kh2)/m0)**2,width(kh2)/m0)

        do 40 i=1,nsf
          do 30 j=1,nsf
            omgllrr(i,j)=dcmplx(omegax,0d0)*(conjg(gl(ksf(i),kf,kn(1)))
     -                    *gl(ksf(j),kf,kn(1)) + gr(ksf(j),kf,kn(1))
     -                    *conjg(gr(ksf(i),kf,kn(1))))
            mfglrrl(i,j)=dcmplx(xf,0d0)*(conjg(gl(ksf(i),kf,kn(1)))
     -                    *gr(ksf(j),kf,kn(1)) + gl(ksf(j),kf,kn(1))
     -                    *conjg(gr(ksf(i),kf,kn(1))))
  30      continue
          gimlr(i)=Imag(conjg(gl(ksf(i),kf,kn(1)))*gr(ksf(i),kf,kn(1)))
          grelr(i)=dble(conjg(gl(ksf(i),kf,kn(1)))*gr(ksf(i),kf,kn(1)))
          glrsq(i)=abs(gl(ksf(i),kf,kn(1)))**2
     -             +abs(gr(ksf(i),kf,kn(1)))**2
          gZsttot(i)=dble(gl(kz,kn(1),kn(1)))*(3*dble(gr(kz,kf,kf))
     -           *abs(dcmplx(omegax,0d0)*conjg(gl(ksf(i),kf,kn(1)))
     -                 -dcmplx(xf,0d0)*gr(ksf(i),kf,kn(1)))**2
     -            -3*dble(gl(kz,kf,kf))
     -              *abs(dcmplx(xf,0d0)*conjg(gl(ksf(i),kf,kn(1)))
     -                 -dcmplx(omegax,0d0)*gr(ksf(i),kf,kn(1)))**2
     -            +(omegax**2-xf**2)*(
     -              dble(gl(kz,kf,kf))*abs(gr(ksf(i),kf,kn(1)))**2
     -             -dble(gr(kz,kf,kf))*abs(gl(ksf(i),kf,kn(1)))**2))
          gZtutot(i)=-dble(gl(kz,kn(1),kn(1)))*(3*dble(gr(kz,kf,kf))
     -           *abs(dcmplx(omegax,0d0)*conjg(gl(ksf(i),kf,kn(1)))
     -                 +dcmplx(xf,0d0)*gr(ksf(i),kf,kn(1)))**2
     -            -3*dble(gl(kz,kf,kf))
     -              *abs(dcmplx(xf,0d0)*conjg(gl(ksf(i),kf,kn(1)))
     -                 +dcmplx(omegax,0d0)*gr(ksf(i),kf,kn(1)))**2
     -            +(omegax**2-xf**2)*(
     -              dble(gl(kz,kf,kf))*abs(gr(ksf(i),kf,kn(1)))**2
     -             -dble(gr(kz,kf,kf))*abs(gl(ksf(i),kf,kn(1)))**2))
  40    continue
        gZZtot=dble(gl(kz,kn(1),kn(1)))**2*(
     -       (2*omegax**2+xf**2)*dble(gl(kz,kf,kf)**2+
     -        gr(kz,kf,kf)**2)-6*xf**2*dble(gr(kz,kf,kf)*gl(kz,kf,kf)))
        gHH(1)=xf*dble(gl(kh1,kf,kf))*dble(gl(kh1,kn(1),kn(1)))
        gHH(2)=xf*dble(gl(kh2,kf,kf))*dble(gl(kh2,kn(1),kn(1)))

        mss=0d0
        muu=0d0
        msu=0d0
        do 60 i=1,nsf
          do 50 j=1,nsf
            mss=mss+dble(dcmplx(1d0,0d0)/ds(i)/conjg(ds(j)))
     -          *abs(omgllrr(i,j)-mfglrrl(i,j))**2
            muu=muu+dble(dcmplx(1d0,0d0)/du(i)/conjg(du(j)))
     -          *abs(omgllrr(i,j)+mfglrrl(i,j))**2
            msu=msu+dble(dcmplx(1d0,0d0)/ds(i)/conjg(du(j)))
     -          *(2*gimlr(i)*gimlr(j)*(omegax**2-xf**2)
     -           -(xf*glrsq(i)-2*omegax*grelr(i))*
     -            (xf*glrsq(j)+2*omegax*grelr(j))/2d0)
  50      continue
  60    continue


        mtt=4*gZZtot/abs(dtZ)**2
        do 80 i=1,2
         do 70 j=1,2
           mtt=mtt+16*gHH(i)*gHH(j)
     -               *dble(dcmplx(1d0,0d0)/dtH(i)/conjg(dtH(j)))
  70     continue
  80    continue
   
        mst=0d0
        mtu=0d0
        do 90 i=1,nsf
          mst=mst+gZsttot(i)*dble(dcmplx(1d0,0d0)/ds(i)/conjg(dtZ))
          mtu=mtu+gZtutot(i)*dble(dcmplx(1d0,0d0)/du(i)/conjg(dtZ))
  90    continue
        do 110 i=1,nsf
         do 100 j=1,2
           mst=mst+dble(dcmplx(1d0,0d0)/ds(i)/conjg(dtH(j)))
     -         *2*gHH(j)*xf*(omegax*glrsq(i)-2*xf*grelr(i))
           mtu=mtu+dble(dcmplx(1d0,0d0)/du(i)/conjg(dtH(j)))
     -         *2*gHH(j)*xf*(omegax*glrsq(i)+2*xf*grelr(i))
 100     continue
 110    continue

c... sum all diagrams and correct for dof
        tmpres=mss+muu+mtt+2*(mst-msu-mtu)        
        tmpres=2.*tmpres*nsf
        if (SMtype.ge.7) tmpres=tmpres*3.


c... mUED
c..  (mf dependence neglected so far)
      elseif (mhtype.eq.2) then  
        if (SMtype.ge.7) return 
        do 120 i=1,nsf
          tmpres=tmpres+
     -      16/3.*gyweak**4*(i/2.)**4*omegax**2*
     -      (((mass(ksf(i))/m0)**2-1.)**2+2*omegax**2*
     -       ((mass(ksf(i))/m0)**4+2.))/
     -      ( (((mass(ksf(i))/m0)**2-1.)**2-2*omegax**2)**2
     -       + 16*(mass(ksf(i))*omegax*width(ksf(i))/m0**2)**2 )
 120    continue
c... finally, take into account dof
        tmpres=2.*tmpres*nsf
        if (SMtype.ge.7) tmpres=tmpres*3.

      
c... user defined
      elseif (mhtype.ge.3) then

      write (*,*) 'Choice ''user'' for dsmhset
     -             not yet implemented in dsmhm2.f!'

      endif


      if (tmpres.lt.0d0) then
        if (memory.ne.idtag) then
          write(*,*) 'WARNING: negative |M|^2 in dsmhm2.f for model',
     &               idtag, ' !'
c          write(*,*) 'SMtype, tmpres = ',SMtype,tmpres
c          write(*,*) mss,muu,mtt,2*msu,2*mst,2*mtu
          memory=idtag
        endif
      else  
        dsmhm2=tmpres
      endif
    
      return

      end

