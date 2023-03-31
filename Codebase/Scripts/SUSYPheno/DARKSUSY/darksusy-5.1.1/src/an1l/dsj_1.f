c====================================================================
c
c   auxiliary function used in:  
c   dsanzgpar.f
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dsj_1(a,b)
      implicit none
      real*8 a,b,root
      root=dsqrt(1.d0-b/a)
      dsj_1=dlog((1.d0+root)/(1.d0-root))
      return
      end
