***********************************************************************
*** dsntdkfoveru is the velocity distribution for the damour-krauss
*** distribution of wimps divided by u, i.e. f(u)/u.
*** the normalization is such that int_0^infinity f(u) du = 1/gtot10
*** where gtot10 is a normalization factor that depends on how effectively
*** the WIMPs scatter off the Sun into these orbits.
*** Multiplying by (rhox/mx) and gtot10 gives the full f(u)/u.
*** The first factor is added in dsntfoveru.f and the second in
*** dsntdkcapearth.f
*** author: j. edsjo (edsjo@physto.se)
*** input: velocity relative to earth [ km s^-1 ]
*** output: f(u) / u [ (km/s)^2 ]
*** date: april 6, 1999
*** modified: 2004-01-30, gtot10 normalization taken out
***********************************************************************

      real*8 function dsntdkfoveru(u)
      implicit none
      include 'dsntdkcom.h'
      include 'dshmcom.h'
      include 'dsntcom.h'

      real*8 u,n_eps,bigu,dsntdkfbigu,delta_e
      parameter(n_eps=0.456579d0)

      bigu=(u/v_earth)**2

c...first calculate f(u)/u
      dsntdkfoveru=dsntdkfbigu(bigu)*2/v_earth**2

c      write(*,*) 'bigu=',bigu,'  dsntdkfbigu=',dsntdkfbigu(bigu)

c...  
      delta_e=0.212d0/(v_sun/220.0d0)

c...now normalize
      dsntdkfoveru=dsntdkfoveru*
     &  n_eps*delta_e

c...f(u)/u is now in units of (km/s)^(-2)

c      if (u.gt.29.0d0) then
c        write(*,*) 'u=',u,'  foveru=',dsntdkfoveru
c      endif

      return
      end
