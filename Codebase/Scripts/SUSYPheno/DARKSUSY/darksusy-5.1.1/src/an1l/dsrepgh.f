c====================================================================
c
c   auxiliary function used in:  
c   dsanggrepar.f 
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dsrepgh(a,b,sq,dq,signm)
      implicit none
      real*8 a,b,sq,dq,signm,sqr,dspi1,dspiw2,dspiw3
      sqr=signm*dsqrt(a)
      dsrepgh=-.5d0/(1.d0-a-b)*(sq+sqr*dq)*dspi1(a,1.d0)
     &     +(.5d0*sq/(1.d0-b)-(sq+sqr*dq)/(1.d0-a-b))*dspiw2(a,b)
     &     -.5d0*b*sq/(1.d0-b)*dspiw3(a,b)
      return
      end
