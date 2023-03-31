c====================================================================
c
c   auxiliary function used in:  
c   dsanggrepar.f 
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dsrepfbox(a,b,sq,dq,signm)
      implicit none
      real*8 a,b,sq,dq,sla,signm,sqr,dspi1,dspiw2,dspiw3
      sqr=signm*dsqrt(a*b)
      sla=((b*sq+sqr*dq)/(1.d0+a-b))*dspi1(a,b)/2.d0+
     &     sq/(1.d0-b)*dspiw2(a,b)/2.d0+
     &     ((b*sq+sqr*dq)/(1.d0+a-b)-
     &     b*sq/2.d0/(1.d0-b))*dspiw3(a,b)
      dsrepfbox=sla
      return
      end
