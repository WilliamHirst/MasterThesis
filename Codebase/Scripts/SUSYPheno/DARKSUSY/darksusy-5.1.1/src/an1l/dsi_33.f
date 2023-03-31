c====================================================================
c
c   auxiliary function used in:  
c   dsanzgpar.f
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dsi_33(r1,r2,r3)
      implicit none
      real*8 dsi_34,r1,r2,r3
      dsi_33=dsi_34(r1,r2,1.d0,r3)
      return
      end
