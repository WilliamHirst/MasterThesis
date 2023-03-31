c====================================================================
c
c   auxiliary function used in:  
c   dspiw2.f 
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dspiw2i(x)
      implicit none
      real*8 a,b,x,arg
      common/dspiw2c/a,b
      arg=dabs((-a*x**2+x*(a+b-1.d0)+1.d0)/(a*x**2+x*(-a+b-1.d0)+1.d0))
      dspiw2i=1.d0/x*dlog(arg)
      return
      end

