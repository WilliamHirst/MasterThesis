***********************************************************************
*** dsntcapsunnum calculates the capture rate at present.
*** Instead of using the approximations in jkg, i.e. a gaussian
*** velocity distribution and approximating all elements as being 
*** at their typical radius, we here integrate numerically over
*** the actual velocity distribution and over the Sun's radius.
*** The velocity distribution used is the one set up by the
*** option veldf in dshmcom.h (see src/hm/dshmudf.f for details)
*** Input: mx = neutralino mass in GeV
***        sigsi = spin-independent scattering cross section (cm^2)
***        sgisd = spin-dependent scattering cross section on protons (cm^2)
*** author: joakim edsjo (edsjo@physto.se)
*** date: 2003-11-26
***********************************************************************
      real*8 function dsntcapsunnum(mx,sigsi,sigsd)
      implicit none
      include 'dshmcom.h'
      include 'dsntcom.h'
      include 'dsntdkcom.h'

      integer type
      real*8 mx,sigsi,sigsd
      real*8 dsntcapsunnumi,dsntfoveru
      real*8 dsntdkgtot10

      real*8 sdmx
      integer sdtype
      common /ntsd/sdmx,sdtype
      save /ntsd/

      external dsntfoveru

c...integrate over the sun according to gould apj 321 (1987) 571.
c...set up things for radial and velocity integration

      ntmx=mx   ! so that internal nt routines know the mass

c...perform integration
      dsntcapsunnum=dsntcapsunnumi(mx,sigsi,sigsd,dsntfoveru)

      return
      end
