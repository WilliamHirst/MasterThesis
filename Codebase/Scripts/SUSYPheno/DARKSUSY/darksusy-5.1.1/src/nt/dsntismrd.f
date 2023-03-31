***********************************************************************
*** function dsntismrd gives the column density of interstellar matter
*** along the line of sight. the model is from ingelman & thunman with
*** rho=rho_0 exp(z/z_0) with rho_0=1.0 nucleon/cm^3 and z_0=0.26 kpc
*** input: b   = galactic latitude (degrees)
***        psi = angle (in the plane) from the galactic centre (degrees)
*** output: column density in units of nucleons / cm^2 kpc/cm.
*** author: joakim edsjo (edsjo@physto.se)
*** date: 1998-09-20
**********************************************************************

      real*8 function dsntismrd(b,psi)
      implicit none

      real*8 b,psi,amax,lmax,z0,rho_0,r_gc,r_0,deg2rad,cb,sb
      parameter(z0=0.26d0,     ! scale hight (kpc)
     &          rho_0=1.0d0,    ! density in the disk nucleons/cm^3
     &          r_gc=8.5d0,    ! our galactocentric distance (kpc)
     &          r_0=12.0d0)    ! radial size of disk
      parameter(deg2rad=1.74532925199d-2) ! degrees -> radians

      sb=max(1.0d-10,sin(b*deg2rad))
      cb=max(1.0d-10,cos(b*deg2rad))
      amax=r_gc*cos(psi*deg2rad)+
     &  sqrt(r_0**2-r_gc**2*(1.0d0-cos(psi*deg2rad)**2))
c      write(*,*) 'amax = ',amax
      lmax=amax/cb
c      write(*,*) 'lmax = ',lmax
      dsntismrd=rho_0*z0/sb*
     &  (1.0d0-exp(-lmax*sb/z0))
c      write(*,*) 'dsntismrd = ',dsntismrd

      return
      end
