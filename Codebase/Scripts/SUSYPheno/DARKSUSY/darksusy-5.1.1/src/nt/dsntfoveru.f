***********************************************************************
*** input: velocity relative to Sun [ km s^-1 ]
*** output: f(u) / u [ cm^-3 (cm/s)^(-2) ]
*** Date: 2004-01-28
***********************************************************************

      real*8 function dsntfoveru(u)
      implicit none
      include 'dshmcom.h'
      include 'dsntcom.h'
      real*8 u,dshmuDF

      dsntfoveru=dshmuDF(u) ! use user-defined halo profile
     &  *1.0d-10            ! (km/s)^(-2) -> (cm/s)^(-2)
     &  *(rhox/ntmx)        ! normalize to local number density

      return
      end
