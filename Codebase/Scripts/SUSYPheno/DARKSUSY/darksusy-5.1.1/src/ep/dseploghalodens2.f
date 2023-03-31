************************************************************************
*** Halo density squared. This function is a function of log(z) and calls
*** the standard halo density routine dshmrho.
*** The jacobian z for integration is also included
*** u = log(z)
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2000-09-03
************************************************************************

************************************************************************
      real*8 function dseploghalodens2(u)
************************************************************************
      implicit none
      include 'dsepcom.h'
      real*8 z,dshmrho,rr,u
      z=exp(u)
      rr=sqrt(r_gc**2+z**2) ! spherical distance to the galactic center
      dseploghalodens2=dshmrho(rr)**2
     &  *z   ! for  the integration in log(z)

c...Rescale for integration to work better
      dseploghalodens2=dseploghalodens2*1d10
      return
      end
