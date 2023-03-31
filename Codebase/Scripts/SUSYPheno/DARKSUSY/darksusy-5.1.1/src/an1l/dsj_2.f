c====================================================================
c
c   auxiliary function used in:  
c   dsanzgpar.f
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dsj_2(b,c)
      implicit none
      real*8 b,c,root
      root=dsqrt(1.d0-4.d0*b/c)
      dsj_2=dlog((1.d0+root)/(1.d0-root))
      return
      end
