c====================================================================
c
c   auxiliary function used in:  
c   dsanzgpar.f
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dsj_3(a,b,c)
      implicit none
      real*8 a,b,c,root
      root=(a-c/4.d0)*dsqrt(1.d0-4.d0*b/c)
      dsj_3=dlog(dabs((1.d0-b+c/4.d0+root)/(1.d0-b+c/4.d0-root)))
      return
      end
