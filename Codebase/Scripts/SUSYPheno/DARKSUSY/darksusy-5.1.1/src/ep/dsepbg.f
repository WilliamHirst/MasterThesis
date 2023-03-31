**********************************************************************
*** function dsepbg gives the differential flux of positrons from the
*** halo coming from secondary sources.
*** these background fluxes are a parameterization by j. edsjo to
*** the results of moskalenko & strong for a model without
*** reacceleration (08-005). ref.: apj 493 (1998) 694.
*** input: positron energy in gev.
*** output: flux in cm^-2 sec^-1 gev^-1 sr^-1
*** author: joakim edsjo, edsjo@physto.se
*** date: 98-07-21
**********************************************************************

      real*8 function dsepbg(egev)
      implicit none

c----------------------------------------------------------------------

      real*8 egev

c----------------------------------------------------------------------

c...standard kam & turner
c      dsepbg=dsembg(egev)*(0.02d0+0.10d0/sqrt(egev))/
c     &  (0.98d0-0.10d0/sqrt(egev))

c...improved to fit protheroe better
c      dsepbg=dsembg(egev)*(0.02d0+0.07d0*egev**(-0.45d0))/
c     &  (0.98d0-0.10d0/sqrt(egev))

c...parameterization of moskalenko & strong
      dsepbg=4.5d0*egev**0.7d0/
     &  (1.0d0+650.0d0*egev**2.3d0+1500.0d0*egev**4.2d0)

      return
      end
