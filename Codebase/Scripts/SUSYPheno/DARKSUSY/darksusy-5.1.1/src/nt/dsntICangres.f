***********************************************************************
*** dsntICangres calculates the 1 sigma angular resolution of the IceCube
*** detector for neutrinos and anti-neutrinos with energy E, in degrees.
*** Here 1 sigma refers to a single dimension of a 2D Gaussian PSF;
*** in terms of absolute containment angle, this therefore corresponds to
*** 39.3% containment, not 68%.
***
*** Input:	log10E	        log(neutrino energy/GeV)
*** Output:			39.3% containment ang resolution (degrees)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 24, 2011
***********************************************************************

      real*8 function dsntICangres(log10E)

      implicit none
      include 'dsntIC.h'

      real*8 log10E
      integer IER
      
      call TSVAL1(nBinsEA,effArea_logEcentres,effArea_AngRes,
     & effArea_AngResderivs,effArea_AngRessigma,0,1,
     & log10E,dsntICangres,IER)

      if (IER .lt. 0) then
        write(*,*) 'TSVAL1 error from angular'
        write(*,*) 'resolution in dsntICangres, code:', IER
        stop
      endif

      end function dsntICangres
