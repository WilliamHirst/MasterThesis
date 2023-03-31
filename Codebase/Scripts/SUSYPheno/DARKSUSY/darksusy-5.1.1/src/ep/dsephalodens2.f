************************************************************************
*** Halo density squared. This function is a function of z and calls
*** the standard halo density routine dshmrho.
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2000-09-03
************************************************************************

************************************************************************
      real*8 function dsephalodens2(z)
************************************************************************
      implicit none
      include 'dsepcom.h'
      real*8 z,dshmrho,rr
      rr=sqrt(r_gc**2+z**2) ! spherical distance to the galactic center
      dsephalodens2=dshmrho(rr)**2

c...Rescale for integration to work better
      dsephalodens2=dsephalodens2*1d10
      return
      end
