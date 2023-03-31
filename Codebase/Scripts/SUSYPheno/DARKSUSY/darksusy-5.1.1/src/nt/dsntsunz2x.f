***********************************************************************
*** The Sun routines uses different variables to describe position
*** in the Sun:
***    r:  radius (in meters)   [0,r_sun]
***    x:  radius in units of r_sun  [0,1]
***    z:  fraction of total column density traversed  [0,1]
***        the column density is either of p or n or the total
***        and the totals are stored in cd_sun
***
*** This routine converts from z (in p, n or total) to x 
***
*** Inputs
***       z =  = fraction of total column density (for chosen type) that
***              has been traversed from the centre of the Sun
***    type = 'N', the total column density (up to that r) is calculated
***         = 'p', the column density on protons is calculated
***         = 'n', the column density on neutrons is calculated
***
*** Outputs
***       x = radius in units of r_sun [0,1] that corresponds to the
***           supplied z value
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2005-11-25
***********************************************************************

      real*8 function dsntsunz2x(z,type)
      implicit none

      include 'dssun.h'
      real*8 z,zpl,dsntsuncdens
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
        write(*,*) 'ERROR in dsntsunx2z: wrong type: ',type
        stop
      endif


c...Interpolate in table

      call dshunt(sdcdens(1,ti),sdn,z*cd_sun(ti),j)
      if (j.lt.sdn) goto 20

      dsntsunz2x=1.0d0

      return

 20   zpl=(z*cd_sun(ti)-sdcdens(j,ti))
     &  /(sdcdens(j+1,ti)-sdcdens(j,ti))

      
      dsntsunz2x=sdr(j)*(1.0d0-zpl)+sdr(j+1)*zpl
      return

      end
