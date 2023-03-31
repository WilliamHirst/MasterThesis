***********************************************************************
*** dsntICnlike calculates the binned number likelihood for a given
*** observed and predicted number of events, marginalising over
*** nuisance parameters arising from uncertainties on the exposure/
*** luminosity and in the theoretical prediction of the signal rate.
***
*** input:  n_tot        total number of observed events
***         theta_tot    total number of predicted events
***         theta_sig      predicted number of signal events
***         sigma_eps    fractional systematic error on exposure/effective
***                       area/luminosity 
***         tau          fractional theoretical error
***		      
*** output: ln(likelihood)  (dimensionless)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 22, 2011, Jun 30, 2011, Jul 21, 2011
***********************************************************************

      real*8 function dsntICnlike(n_tot,theta_tot,theta_sig,sigma_eps,
     & tau)

      implicit none
      include 'dsntIC.h'

      integer n_tot
      real*8 theta_tot,theta_sig,sigma_eps
      real*8 lnlike, tau, sigsq, dslnpoisint
      
      !compute Poissonian likelihood

      if (theta_tot.lt.0.d0.or.n_tot.lt.0.or.theta_sig.lt.0.d0) then 
        write(*,*) theta_tot, theta_sig, n_tot
        stop 'Error: something has gone negative in dsntICnlike!'
      endif

      sigsq = sigma_eps*sigma_eps + tau*tau
      lnlike = dslnpoisint(theta_tot,theta_sig,n_tot,sigsq,
     & sysErrDist_logNorm)

      if (lnlike.gt.0.d0) then 
        write(*,*) lnlike
        stop 'Error: negative lnlike in dsntICnlike!'
      endif

      dsntICnlike = lnlike

      end function dsntICnlike
