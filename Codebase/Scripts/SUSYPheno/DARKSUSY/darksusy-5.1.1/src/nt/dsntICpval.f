***********************************************************************
*** dsntICpval performs the calculation of the p-value 
*** in IceCube calculations based on Poissonian statistics.
*** 
*** input:  ntot      observed number of events
***         theta_tot total predicted number of events
***         theta_sig   predicted number of signal events
*** output:           p value
***
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: May, June, July, Dec 2011
***********************************************************************

      double precision function dsntICpval(ntot,theta_tot,theta_sig)

      implicit none
      include 'dsntIC.h'

      integer ntot, i
      real*8 theta_tot, theta_sig, sigma, lnpval, lnpin, lngesum

      if (theta_tot.lt.0.d0.or.ntot.lt.0.or.theta_sig.lt.0.d0) then 
        write(*,*) theta_tot, theta_sig, ntot
        stop 'Error: something has gone negative in dsntICpval!'
      endif

      sigma = dsqrt(EAErr_max*EAErr_max+theoryErr*theoryErr)
      if (sysErrDist_logNorm) then
        !Treat percentage error as log-normal distributed
        call dslnpilnsum(ntot,theta_tot-theta_sig,theta_sig,sigma,
     &                   lnpin,lnpval,lngesum)
      else
        !Treat percentage error as Gaussianly-distributed
        call dslnpinsum(ntot,theta_tot-theta_sig,theta_sig,sigma,
     &                  lnpin,lnpval,lngesum)
      endif
      dsntICpval = dexp(lnpval)

      end function dsntICpval

