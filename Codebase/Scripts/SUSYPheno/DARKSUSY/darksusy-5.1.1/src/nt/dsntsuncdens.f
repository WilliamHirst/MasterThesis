***********************************************************************
*** This routine uses a derived column density from the BP2000 model
*** The data in sdcens() is calculated by dsntsunread.f. 
***********************************************************************

***********************************************************************
*** dsntsuncdens gives the column density in the Sun from the
*** centre out the tha given radius r (in meters).
*** The radius should be given in m and the column density is returned in
*** g/cm^2
*** if type = 'N', the total column density (up to that r) is calculated
***         = 'p', the column density on protons is calculated
***         = 'n', the column density on neutrons is calculated
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2005-11-25
***********************************************************************

      real*8 function dsntsuncdens(r,type)
      implicit none

      include 'dssun.h'
      real*8 r,rpl
      character*1 type
      integer i,j,ti

c...Check if data file is loaded
      call dsntsunread

      if (type.eq.'N') then 
        ti=0
      elseif (type.eq.'p') then
        ti=1
      elseif (type.eq.'n') then
        ti=2
      else
        write(*,*) 'ERROR in dsntsundens: wrong type: ',type
        stop
      endif

      if (r.ge.r_sun) then
        dsntsuncdens=cd_sun(ti) ! total column density
        return
       endif

      if (r.le.0.0d0) then
        dsntsuncdens=0.0d0 ! no column density
        return
      endif

c...Interpolate in table

      call dshunt(sdr,sdn,r/r_sun,j)
      if (j.lt.sdn) goto 20

      dsntsuncdens=sdcdens(sdn,ti)

      return

 20   rpl=(r-sdr(j)*r_sun)/(sdr(j+1)*r_sun-sdr(j)*r_sun)

      dsntsuncdens=sdcdens(j,ti)*(1.0d0-rpl)+sdcdens(j+1,ti)*rpl

      return

      end
