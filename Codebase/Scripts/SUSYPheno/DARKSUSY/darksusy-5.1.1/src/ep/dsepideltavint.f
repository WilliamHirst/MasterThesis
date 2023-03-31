************************************************************************
*** positron propagation routines.
*** the integrand in the integration for i(delta v).
*** To speed up the integration, u=ln(r) is used as an integration variable
*** Author: J. Edsjo (edsjo@physto.se), based on routines by
*** e.a. baltz (eabaltz@astron.berkeley.edu)
*** Date: 2004-01-26
************************************************************************

************************************************************************
      real*8 function dsepideltavint(u)
************************************************************************
      implicit none
      include 'dsepcom.h'
      real*8 dsepf,r,dsbessei0,u

      real*8 dv
      common/epint/dv

      r=exp(u)
      dsepideltavint=r*dsepf(r)*
     &  dsbessei0(2.0d0*r_e*r/dv)*
     &  exp(-(r_e-r)**2/dv)
     &  /(0.5d0*dv)   ! 1/dv is the first line in eq. (25), 1/0.5 from sum_+-
     &  *1.d10  ! rescale by 1.d10 for integration to work fine
     &  *r   ! dr/du

      return
      end
