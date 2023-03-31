c====================================================================
c
c   auxiliary function used in:  
c   dsti_5.f
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dsfl4c1(r1,r2,r3,r4,r5)
      implicit none
      real*8 r1,r2,r3,r4,r5,c0,c1,c2,c3
      c3=((r3+r4)*0.5d0-r2-r5)/(r1+r5)
      c0=0.5d0*(r3+r4)-r5
      c1=0.5d0*(r3-r4)
      c2=r5
      dsfl4c1=(c1*(1.d0-c3)+2*c2*c3*(1.d0-c3))/(c0+c1*c3+c2*c3*c3)
      return
      end
