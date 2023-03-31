***********************************************************************
*** dsntcapearthnum calculates the capture rate at present.
*** Intead of using the assumptions of Gould (i.e. capture as in
*** free space), a tabulted velocity distribution based on detailed
*** numerical simulations of Johan Lundberg is used.
*** A numerical intregration has to be performed instead of the
*** convenient expressions in jkg.
*** Input: mx = neutralino mass in GeV
***        sigsi = spin-independent scattering cross section in cm^2
***        type = type of velocity distribution
***          1 = best estimate of distribution at Earth form numerical sims
***          2 = conservative estimate, only including free orbits and
***              jupiter-crossing orbits
***          3 = ultraconservative estimate, only including free orbits
***          4 = as if Earth was in free space, i.e. with a Gaussian
***              (Gaussian provided by dsntsdfoveru)
***          5 = as if Earth was in free space, i.e. with a Gaussian
***              (Gaussian provided by dsntdkfoverugauss)
***          6 = Damour-Krauss population (this is per gtot10), i.e.
***              multiply with gtot10 to get the full capture rate
***          Note: 1 is the best estimate of the distribution at Earth
***          and should be used as a default
*** author: joakim edsjo (edsjo@physto.se)
*** date: July 10, 2003
***********************************************************************
      real*8 function dsntcapearthnum(mx,sigsi)
      implicit none
      include 'dshmcom.h'
      include 'dsntcom.h'
      include 'dsntdkcom.h'

      integer type
      real*8 mx,sigsi,sigsd
      real*8 dsntcapearthnumi
      real*8 dsntfoveruearth,dsntdkfoveru

      real*8 sdmx
      integer sdtype
      common /ntsd/sdmx,sdtype
      save /ntsd/

      external dsntfoveruearth,dsntdkfoveru

c...integrate over the earth according to gould apj 321 (1987) 571.
c...set up things for radial and velocity integration

      ntmx=mx   ! so that internal routines knows the mass

c...perform integration
      if (veldfearth.eq.'dk') then  ! use special integration, option 2
        dsntcapearthnum=dsntcapearthnumi(mx,sigsi,dsntfoveruearth,2)
      else ! general case
        dsntcapearthnum=dsntcapearthnumi(mx,sigsi,dsntfoveruearth,1)
      endif

      return
      end
