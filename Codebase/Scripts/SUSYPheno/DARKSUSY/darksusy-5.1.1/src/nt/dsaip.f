

      real*8 function dsaip(x)
c   dsairy function derivative
c   lb 990224
      implicit none
      include 'dsge.h'
      real*8 x,z,k23,dbskr3
      if (x.lt.0.d0) then
         write(*,*) 'dsairy fcn of negative argument - not implemented'
         stop
         else
           z=2.d0/3.d0*x**(1.5d0)
           k23=dbskr3(z,2)  ! cernlib function for k{2/3}
           dsaip=-x/pi/sqrt(3.d0)*k23  ! num recipes (6.7.45)
         endif
      return
      end
