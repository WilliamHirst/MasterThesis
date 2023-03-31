***********************************************************************
*** dsntICsignal computes the predicted number of neutrino events due
*** to neutralino annihilation in each superbin, saving global variables
*** for later access by likelihood codes.
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 22, 2011
***********************************************************************


      subroutine dsntICsignal

      implicit none
      include 'dsntIC.h'

      real*8 integral, eps, dsf_int, dsntICsigintegrand, upperLimit
      integer i
      parameter (eps = 1.d-3)
      external dsntICsigintegrand
 
      theta_S = 0.d0
      
      do i = 1,nBinsEAError

        if (log10mwimp .lt. EAlogE_inEAErrBins(1,i)) then

          theta_Snu(i) = 0.d0
          theta_Snubar(i) = 0.d0
 
        else

          if (log10mwimp .lt. EAlogE_inEAErrBins(2,i)) then
            upperLimit = log10mwimp
          else
            upperLimit = EAlogE_inEAErrBins(2,i)
          endif

          ptypeshare = 1
          integral = dsf_int(dsntICsigintegrand,EAlogE_inEAErrBins(1,i),
     &     upperLimit,eps)
          theta_Snu(i) = integral * dlog(10.d0) * exp_time * annrateIC

          ptypeshare = 2
          integral = dsf_int(dsntICsigintegrand,EAlogE_inEAErrBins(1,i),
     &     upperLimit,eps)
          theta_Snubar(i) = integral * dlog(10.d0) * exp_time *annrateIC

        endif

        theta_S(i) = theta_Snu(i) + theta_Snubar(i) 

      enddo

      theta_S_total = sum(theta_S)

      end subroutine dsntICsignal

