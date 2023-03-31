***********************************************************************
*** The halo velocity profile in the Maxwell-Boltzmann (Gaussian)
*** approximation with the addition of a dark disk.
*** input: velocity relative to earth [ km s^-1 ]
*** output: f(u) / u [ (km/ s)^(-2) ]
*** date: april 6, 1999
*** Modified: 2004-01-29, 2009-04-30
*** Modified for best case dark disc by Lars Rosenstrom 
***********************************************************************

      real*8 function dshmuDFgaussdd(u)
      implicit none
      include 'dshmcom.h'
      real*8 u

      dshmuDFgaussdd=sqrt(3./2.)/(vd_3d*v_sun)/sqrt(3.1415)*( ! halo
     &  exp(-1.5*(u-v_sun)**2/vd_3d**2)                     ! halo
     &  -exp(-1.5*(u+v_sun)**2/vd_3d**2))                   ! halo
     &  +delta_dd*(sqrt(3./2.)/(vd_3d_dd*v_lag_dd)/sqrt(3.1415)*( ! dark disk
     &  exp(-1.5*(u-v_lag_dd)**2/vd_3d_dd**2)    ! dark disk
     &  -exp(-1.5*(u+v_lag_dd)**2/vd_3d_dd**2))) ! dark disk

      return
      end

      

