c====================================================================
c
c   auxiliary function used in:  
c   not used, this is equivalent to dspiw3.f
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dsi_32(r1,r2)
      implicit none
      real*8 dsi_34,r1,r2
      dsi_32=dsi_34(r1,r2,1.d0,0.d0)
      return
      end
