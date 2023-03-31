***********************************************************************
*** dsntICdPhi_SdE provides the differential spectrum of neutrinos or
*** anti-neutrinos to be expected at IceCube, after weighting 
*** by the effective area and the angular loss factor (accounting for
*** the loss of neutrinos to regions beyond the analysis cone due to the
*** PSF).  Note that the annihilation rate in the Sun is not included.
*** The energy dispersion is also not taken into
*** account here, as this function is intended for use as either
*** a) input into an integral over all energies, which will give the
***    total number of epected events (essentially identical regardless
***     of whether energy dispersion is included or not)
*** b) a prefactor, to be weighted by the energy dispersion 
*** Input:	log10E		log10(neutrino energy/GeV)
***		E		neutrino energy/GeV
*** 		ptype	= 1	neutrinos
***			= 2	anti-neutrinos 
*** Output:                     differential flux (GeV^-1)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 24, 2011
***********************************************************************


      real*8 function dsntICdPhi_SdE(log10E,E,ptype)

      implicit none
      include 'dsntIC.h'

      real*8 E, log10E, effArea, angLossFac, dsf_int
      real*8 sigma, dsntICangintegrand, spec
      real*8 dsntICeffarea, dsntICangres, dsntmuonyield 
      integer ptype, istat, IER, i
      external dsntICangintegrand 

      !Obtain differential neutrino or anti-neutrino 
      !flux spectrum as it arrives at the detector; spec in m^-2 GeV^-1
      spec = 1.d-30 * dsntmuonyield(E,10.d0,'su',3,1,ptype,istat)
      
      !Obtain effective area for relevant species and energy; effArea in m^2
      effArea = dsntICeffarea(log10E, ptype)

      !Obtain angular resolution for nu and nubar with energy E; sigma in degrees
      sigma = dsntICangRes(log10E)

      !Obtain angular loss factor for neutrinos and anti-neutrinos
      !with energy E, to account for events that will leak out of
      !the analysis cone.
      angLossFac = 1.d0-exp(phi_max_deg*phi_max_deg/(-2.d0*
     & sigma*sigma))

      !Put everything together to obtain the spectrum observed
      !by IceCube; dsntICdPhi_SdE in GeV^-1
      dsntICdPhi_SdE = spec * effArea * angLossFac

      end function dsntICdPhi_SdE


