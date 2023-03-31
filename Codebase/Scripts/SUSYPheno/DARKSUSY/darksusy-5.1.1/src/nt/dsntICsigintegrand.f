***********************************************************************
*** dsntICsigintegrad provides the integrand for computing theta_S,
*** the number of expected signal events at IceCube.
*** Input:		log10E		log10(neutrino energy/GeV)
*** Hidden Input:	ptypeshare   1 => return integrand for neutrinos
***                                  2 => for anti-neutrinos 
*** Output:             integrand       dimensionless
***
*** Note that the factor of ln(10) in the logarithmic integral has
*** been left for post-multiplication in order to increase efficiency.      
***
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 24, 2011
***********************************************************************

      real*8 function dsntICsigintegrand(log10E)

      implicit none
      include 'dsntIC.h'

      real*8 log10E, E, dsntICdPhi_SdE

      E = 10.d0**log10E

      !Return either predicted differential neutrino or
      !anti-neutrino signal in IceCube, depending on 
      !ptypeshare
      dsntICsigintegrand = dsntICdPhi_SdE(log10E,E,ptypeshare)*E

      end function dsntICsigintegrand
