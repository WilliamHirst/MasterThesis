***********************************************************************
*** dsntICbglikeprecomp calls the calculation of the p-value for the
*** background in IceCube calculations based on Poissonian statistics.
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: May, July 2011
***********************************************************************

      subroutine dsntICbglikeprecomp

      implicit none
      include 'dsntIC.h'

      real*8 dsntICpval

      BGpvalPoissonian = dsntICpval(sum(nEvents_inEAErrBins),sum(theta_BG),0.d0)
      pvalBGPoisComputed = .true.

      end subroutine dsntICbglikeprecomp
