***********************************************************************
*** The Sun routines uses different variables to describe position
*** in the Sun:
***    r:  radius (in meters)   [0,r_sun]
***    x:  radius in units of r_sun  [0,1]
***    z:  fraction of total column density traversed  [0,1]
***        the column density is either of p or n or the total
***        and the totals are stored in cd_sun
***
*** This routine converts from x to z (in p, n or total)
***
*** Inputs
***       x = radius in units of r_sun [0,1]
***    type = 'N', the total column density (up to that r) is calculated
***         = 'p', the column density on protons is calculated
***         = 'n', the column density on neutrons is calculated
***
*** Outputs
***       z = fraction of total column density (for chosen type) that
***           has been traversed from the centre of the Sun out to x.
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2005-11-25
***********************************************************************

      real*8 function dsntsunx2z(x,type)
      implicit none

      include 'dssun.h'
      real*8 x,rpl,dsntsuncdens
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

      dsntsunx2z=dsntsuncdens(x*r_sun,type)/cd_sun(ti)

      return

      end
