c====================================================================
c
c   auxiliary function used in:  
c   dsanzgpar.f  dsi_13.f
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dsi_13(r1,r2,r3)
      implicit none
      real*8 dsi_14,r1,r2,r3
      dsi_13=dsi_14(r1,r2,r2,r3)
      return
      end
