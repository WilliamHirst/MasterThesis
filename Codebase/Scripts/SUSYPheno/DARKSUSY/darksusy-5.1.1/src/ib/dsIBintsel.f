*****************************************************************************
*** auxiliary function that selects integrand for integration routines
*** author: Torsten Bringmann 2008-03-10
*****************************************************************************

      real*8 function dsIBintsel(yint) 

      implicit none
      include 'dsibcom.h'

      real*8 yint, xz
      real*8 dsIBwwdxdy, dsIBwhdxdy, dsIBhhdxdy, dsIBffdxdy,
     &       dsIBfsrdxdy

c... intch = 1..12 for IB photons
c... intch = 101..112 for FSR photons
c... intch = -1...-12 for IB positrons

      if (intch.eq.1) then
         dsIBintsel = dsIBwwdxdy(1,ibcom_x,yint)
      elseif (intch.eq.2) then
         dsIBintsel = dsIBwhdxdy(2,ibcom_x,yint)
      elseif (intch.eq.3) then
         dsIBintsel = dsIBhhdxdy(3,ibcom_x,yint)
      elseif (intch.ge.4.and.intch.le.12) then
         dsIBintsel = dsIBffdxdy(intch,ibcom_x,yint)
      elseif (intch.ge.101.and.intch.le.112) then
         dsIBintsel = dsIBfsrdxdy(intch-100,ibcom_x,yint)
      elseif (intch.eq.-1) then      
         xz=1d0-ibcom_z+yint-
     &           ibcom_mp2**2/4d0/ibcom_mx**2
         dsIBintsel = dsIBwwdxdy(1,xz,yint)
      elseif (intch.le.-4.and.intch.ge.-12) then
         xz=1d0-ibcom_z+yint-
     &           ibcom_mp2**2/4d0/ibcom_mx**2
         dsIBintsel = dsIBffdxdy(-intch,xz,yint)
      else
         dsIBintsel=0.d0
      endif
 
      end



