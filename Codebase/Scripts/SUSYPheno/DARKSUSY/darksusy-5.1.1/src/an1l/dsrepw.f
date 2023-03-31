c====================================================================
c
c   auxiliary function used in:  
c   dsanggrepar.f 
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dsrepw(a,b,sq,dq,signm)
      implicit none
      real*8 a,b,sq,dq,signm,sqr,dspi1,dspiw2,dspiw3
      sqr=signm*dsqrt(a)
      dsrepw=2.d0*sq*(a-b)/(1.d0+a-b)*dspi1(a,b)
     &     +(sq-2.d0*dq*sqr)/(1.d0-a-b)*dspi1(a,1.d0)
     &     +(2.d0*(sq-2.d0*dq*sqr)/(1.d0-a-b)-
     &     (3.d0*sq-4.d0*dq*sqr)/(1.d0-b))*dspiw2(a,b)
     &     +(((b+2.d0)*sq-4.d0*sqr*dq)/(1.d0-b)
     &     -2.d0*sq*(1.d0-a+b)/(1.d0+a-b))*dspiw3(a,b)
      return
      end
