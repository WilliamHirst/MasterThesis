      subroutine dsanalbe(alph,bet)
c_______________________________________________________________________
c  determine alph and bet for the integration.
c  modified: joakim edsjo (edsjo@physto.se) 97-09-09
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      include 'dsandwcom.h'
      real*8 dsandwdcos,p1,p2,r1,r2,sig11,sig12,sig21,sig22,
     &  cth1, cth2,alph,bet,mx

c-----------------------------------------------------------------------

      mx=mco(1)
      cth1=0.0d0
      cth2=1.0d0
      p1=5.0d0*mx
      sig11=2.0d0*dsandwdcos(p1,cth1)
      sig12=dsandwdcos(p1,cth2)+dsandwdcos(p1,-cth2)
      r1=sig11/sig12

      p2=10.0d0*mx
      sig21=2.0d0*dsandwdcos(p2,cth1)
      sig22=dsandwdcos(p2,cth2)+dsandwdcos(p2,-cth2)
      r2=sig21/sig22
      alph=1+log(r2/r1)/(2.0d0*log(p2/p1))
c...make sure alpha <1 (with some margin)
      alph=min(alph,0.75d0)
      alph=max(alph,-1.0d0)
      if (alph.gt.-0.1d0.and.alph.le.0.0d0)
     &  alph=-0.1d0
      if (alph.lt.0.1d0.and.alph.ge.0.0d0)
     &  alph=0.1d0

      bet=mx**2/p1**2*(r1**(1.0d0/(alph-1.0d0))-1)
c...make sure beta is ok
      bet=min(bet,100.0d0)
      bet=max(bet,1.0d-2)

c      write(*,*) 'sig11 = ',sig11
c      write(*,*) 'sig12 = ',sig12
c      write(*,*) 'sig21 = ',sig21
c      write(*,*) 'sig22 = ',sig22
c      write(*,*) 'r1 = ',r1
c      write(*,*) 'r2 = ',r2
c      write(*,*) 'alph = ',alph,'   bet = ',bet

      return
      end
