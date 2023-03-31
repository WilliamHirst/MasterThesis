

      real*8 function dsbip(x)
c   dsairy function derivative
c   lb 990224
      implicit none
      include 'dsge.h'
      real*8 x,z,i23,k23,dbsir3,dbskr3
      if (x.lt.0.d0) then
         write(*,*) 'dsairy fcn of negative argument - not implemented'
         stop
         else
           z=2.d0/3.d0*x**(1.5d0)
           i23=dbsir3(z,2)  ! cernlib function for i{2/3}
           k23=dbskr3(z,2)  ! cernlib function for k{2/3}
           dsbip=x*(2.d0*i23/sqrt(3.d0)+k23/pi)  ! num recipes (6.7.45)
         endif
      return
      end
