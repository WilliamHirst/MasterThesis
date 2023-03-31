***********************************************************************
*** dsntICspecintegrand provides the integrand for computing the
*** signal component of the spectral likelihood in dsntICspeclike.
***
*** Input:		log10E       log10(neutrino energy/GeV)
*** Hidden Input:	nchanshare   number of hit DOMs for this event
***                     thetashare   total number of predicted signal 
***                                   events (nu + nubar) 
*** Output:             integrand    (chan^-1)
***
*** Note that the factor of ln(10) in the logarithmic integral has
*** been left for post-multiplication in order to increase efficiency.      
***
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 2011
***********************************************************************

      real*8 function dsntICspecintegrand(log10E)

      implicit none
      include 'dsntIC.h'

      real*8 log10E, E, dsntICdP_SdE, edisp, specpdf
      real*8 dsntICedisp

      E = 10.d0**log10E

      !Return energy dispersion
      edisp = dsntICedisp(log10E,nchanshare)
      !Return spectral probability distribution function
      specpdf = dsntICdP_SdE(log10E,E,thetashare)
      !Put them together, weight by E to give full integrand
      dsntICspecintegrand = specpdf * E * edisp

      end function dsntICspecintegrand
