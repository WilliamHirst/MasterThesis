**********************************************************************
*** function dsembg gives the differential flux of electrons from the
*** halo coming from primary and secondary (background) sources.
*** these background fluxes are a parameterization by j. edsjo to
*** the results of moskalenko & strong for a model without
*** reacceleration (08-005). ref.: apj 493 (1998) 694.
*** input: positron energy in gev.
*** output: flux in cm^-2 sec^-1 gev^-1 sr^-1
*** author: joakim edsjo, edsjo@physto.se
*** date: 98-07-21
**********************************************************************

      real*8 function dsembg(egev)
      implicit none

c----------------------------------------------------------------------

      real*8 egev

c----------------------------------------------------------------------

c...kam & turner values
c      dsembg=0.07 *(egev)**(-3.3d0)

c...heat 94 fit
c      dsembg=0.0227*egev**(-3.09d0)

c...parameterization of moskalenko & strong
c...first the primary component
      dsembg=0.16d0*egev**(-1.1d0)/
     &  (1.0d0+11.0d0*egev**0.9d0+3.2d0*egev**2.15d0)
c...then the secondary component
      dsembg=dsembg+0.700d0*egev**(0.7d0)/
     &  (1.0d0+110.0d0*egev**1.5d0+600.0d0*egev**2.9d0+
     &  580.0d0*egev**4.2d0)

      return
      end
