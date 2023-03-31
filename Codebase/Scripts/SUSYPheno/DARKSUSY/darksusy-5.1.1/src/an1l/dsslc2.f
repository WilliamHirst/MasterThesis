c====================================================================
c
c   auxiliary function used in:  
c   dsi_14.f  dsi_24.f  dsi_34.f 
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dsslc2(r1,r2,r3)
c this function gives the coefficient c2 of slog
      implicit none
      real*8 r1,r2,r3
      dsslc2=-r1/r3
      return
      end
