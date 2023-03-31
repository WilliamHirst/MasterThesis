***********************************************************************
*** dsntICdP_SdE provides the differential probability for any neutrino
*** or anti-neutrino signal event in IceCube to have the energy E.
***
*** Input:	log10E	  log10(neutrino energy/GeV)
***		tS_tot	  the total number of signal events expected,
***			  of any energy or type (nu or nubar)
*** Output:               differential probability (GeV^-1)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 24, 2011
***********************************************************************

      real*8 function dsntICdP_SdE(log10E,E,tS_tot)

      implicit none
      include 'dsntIC.h'
      include 'dsntcom.h'

      real*8 log10E, tS_tot, dsntICdPhi_SdE, E

      !Obtain the neutrino contribution in yr^-1 GeV^-1
      dsntICdP_SdE = dsntICdPhi_SdE(log10E,E,1)
      !Obtain the anti-neutrino contribution in yr^-1 GeV^-1
      dsntICdP_SdE = dsntICdP_SdE + dsntICdPhi_SdE(log10E,E,2)
      !Normalise the distribution to differential probability per GeV
      !dsntICdP_SdE in GeV^-1, exp_time in s, annrateIC in s^-1
      dsntICdP_SdE = dsntICdP_SdE * exp_time * annrateIC / tS_tot

      end function dsntICdP_SdE

