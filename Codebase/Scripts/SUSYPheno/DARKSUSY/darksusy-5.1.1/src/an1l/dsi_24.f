      real*8 function dsi_24(r1,r2,r3,r4)
      implicit none
      real*8 r1,r2,r3,r4,c1,c2,dsslc1,dsslc2,dslp2
      dsi_24=0.d0
c sum the first term
      c1=dsslc1(r1-2.d0*r4,r2,r3)
      c2=dsslc2(r1-2.d0*r4,r2,r3)
      dsi_24=dsi_24+dslp2(c1,c2)
c subtract the second term
      c1=dsslc1(-r1,r2,r3)
      c2=dsslc2(-r1,r2,r3)
      dsi_24=dsi_24-dslp2(c1,c2)
      return
      end
