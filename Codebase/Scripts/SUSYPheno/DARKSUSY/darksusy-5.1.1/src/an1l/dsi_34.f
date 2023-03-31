c====================================================================
c
c   auxiliary function used in:  
c   dsanzgpar.f  dsi_32.f  dsi_33.f
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dsi_34(r1,r2,r3,r4)
      implicit none
      real*8 r1,r2,r3,r4,c1,c2,dsslc1,dsslc2,dslp2
      dsi_34=0.d0
c sum the first term
      c1=dsslc1(r1-2.d0*r4,r3,r2)
      c2=dsslc2(r1-2.d0*r4,r3,r2)
      dsi_34=dsi_34+dslp2(c1,c2)
c subtract the second term
      c1=dsslc1(-r1,r3,r2)
      c2=dsslc2(-r1,r3,r2)
      dsi_34=dsi_34-dslp2(c1,c2)
      return
      end
