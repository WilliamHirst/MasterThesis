***********************************************************************
*** dsntICbgangpdf provides the interpolated probability distribution
*** function for the observed angle between the arrival direction
*** of background events and the Solar direction on the sky, normalised
*** to unit integrated probability over the range [0,phi_cut].  
***
*** Input:	cosphi		cos(arrival angle of event)
*** Output:                     pdf (degrees^-1)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 24, 2011
***********************************************************************

      real*8 function dsntICbgangpdf(cosphi)

      implicit none
      include 'dsntIC.h'
      include 'dsge.h'

      real*8 cosphi, CosToInvDeg
      integer IER

      call TSVAL1(nBinsBGAng,cos(BGangdist_phi),BGangdist_prob,
     & BGangdist_derivs,BGangdist_sigma,0,1,cosphi,dsntICbgangpdf,IER)

      CosToInvDeg = sqrt(1.d0 - cosphi*cosphi) * pi / 180.d0
      dsntICbgangpdf = dsntICbgangpdf * CosToInvDeg / BGangdist_conenorm

      if (IER .lt. 0) then
        write(*,*) 'TSVAL1 error from background spectral'
        write(*,*) 'distribution in dsntICbgspec, code:', IER
        stop
      endif

      end function dsntICbgangpdf
