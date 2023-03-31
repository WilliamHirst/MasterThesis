c====================================================================
c
c   auxiliary function used in:  
c   dsanggrepar.f  dsanglglre.f  dsi_12.f  dsrepfbox.f  dsrepgh.f  dsrepw.f 
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dspi1(a,b)
      implicit none
      real*8 a,b
      if (a.lt.b) then
      dspi1=-2.d0*(datan(1.d0/dsqrt(b/a-1.d0)))**2
      else
      dspi1=(dlog((1.d0+dsqrt(1.d0-b/a))/(1.d0-dsqrt(1.d0-b/a))))**2-
     &    (4.d0*datan(1.d0))**2
      dspi1=dspi1/2.d0
      endif
      return
      end
