      complex*16 function dsb0loop(q,m1,m2)
c_______________________________________________________________________
c  the two-point function b0(q,m1,m1).
c     uses two-point function b_0 from m. drees, k. hagiwara and a.
c     yamada, phys. rev. d45 (1992) 1725.
c  author: joakim edsjo (edsjo@physto.se) 97-02-11
c=======================================================================
      implicit none
      include 'dsge.h'
      real*8 q,m1,m2,q2,m12,m22,bigdelta,delta,sigma,beta
      complex*16 f0,l
      q2=q**2
      m12=m1**2
      m22=m2**2
      bigdelta = log(q2)
      delta = (m12-m22)/q2
      sigma = (m12+m22)/q2
      beta = sqrt(abs(1.0d0-2.0d0*sigma+delta**2))
      if (q2.gt.(m1+m2)**2) then
        l=0.5d0*log((1.0d0+beta-sigma)/(1.0d0-beta-sigma))-
     &    pi*cmplx(.0d0,1.0d0)
      elseif (q2.lt.(m1-m2)**2) then
        l=0.5d0*log((1.0d0+beta-sigma)/(1.0d0-beta-sigma))
      else
        l=atan((1.0d0-delta)/beta) + atan((1.0d0+delta)/beta)
      endif
      f0=log(m1*m2)-delta*log(m2/m1)-2.0d0+beta*l
      dsb0loop=bigdelta-f0
      return
      end
