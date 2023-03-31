c====================================================================
c
c   auxiliary function used in:  
c   dspiw3.f 
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dspiw3i(x)
      implicit none
      real*8 a,b,x,arg
      common/dspiw3c/a,b
      arg=dabs((-a*x**2+x*(a-b+1.d0)+b)/(a*x**2+x*(-a-b+1.d0)+b))
      dspiw3i=1.d0/x*dlog(arg)
      return
      end

