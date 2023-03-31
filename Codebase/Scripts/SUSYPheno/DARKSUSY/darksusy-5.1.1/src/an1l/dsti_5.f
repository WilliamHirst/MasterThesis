c====================================================================
c
c   auxiliary function used in:  
c   dsti_214.f  dsti_224.f  dsti_23.f dsti_33.f 
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dsti_5(r1,r2,r3,r4,r5)
      implicit none
      real*8 r1,r2,r3,r4,r5,c1,c2,dslp2
      real*8 dsfl1c1,dsfl1c2,dsfl2c1,dsfl2c2,dsfl3c1,dsfl3c2,
     &  dsfl4c1,dsfl4c2
      dsti_5=0.d0
c subtrac the first term
      c1=dsfl1c1(r1-2.d0*r5,r2,r3,r4,r5)
      c2=dsfl1c2(r1-2.d0*r5,r2,r3,r4,r5)
      dsti_5=dsti_5-dslp2(c1,c2)
c add the second term
      c1=dsfl2c1(r1-2.d0*r5,r2,r3,r4,r5)
      c2=dsfl2c2(r1-2.d0*r5,r2,r3,r4,r5)
      dsti_5=dsti_5+dslp2(c1,c2)
c add the third term
      c1=dsfl3c1(r1-2.d0*r5,r2,r3,r4,r5)
      c2=dsfl3c2(r1-2.d0*r5,r2,r3,r4,r5)
      dsti_5=dsti_5+dslp2(c1,c2)
c subtract the fourth term
      c1=dsfl4c1(r1-2.d0*r5,r2,r3,r4,r5)
      c2=dsfl4c2(r1-2.d0*r5,r2,r3,r4,r5)
      dsti_5=dsti_5-dslp2(c1,c2)
c add the first term
      c1=dsfl1c1(-r1,r2,r4,r3,r5)
      c2=dsfl1c2(-r1,r2,r4,r3,r5)
      dsti_5=dsti_5+dslp2(c1,c2)
c subtract the second term
      c1=dsfl2c1(-r1,r2,r4,r3,r5)
      c2=dsfl2c2(-r1,r2,r4,r3,r5)
      dsti_5=dsti_5-dslp2(c1,c2)
c subtract the third term
      c1=dsfl3c1(-r1,r2,r4,r3,r5)
      c2=dsfl3c2(-r1,r2,r4,r3,r5)
      dsti_5=dsti_5-dslp2(c1,c2)
c add the fourth term
      c1=dsfl4c1(-r1,r2,r4,r3,r5)
      c2=dsfl4c2(-r1,r2,r4,r3,r5)
      dsti_5=dsti_5+dslp2(c1,c2)
      return
      end












