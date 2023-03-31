***********************************************************************
*** dsntsundens gives the density in the Sun as a function of radius
*** the radius should be given in m and the density is returned in
*** g/cm^3
*** Density and element mass fractions up to O16 are from the standard
*** solar model BP2000 of Bahcall, Pinsonneault and Basu,
*** ApJ 555 (2001) 990.
*** The mass fractions for heavier elements are from N. Grevesse and
*** A.J. Sauval, Space Science Reviews 85 (1998) 161 normalized such that
*** their total mass fractions matches that of the heavier elements in 
*** the BP2000 model.
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2003-11-25
***********************************************************************

      real*8 function dsntsundens(r)
      implicit none

      include 'dssun.h'
      real*8 r,rpl
      integer i,j

c...Check if data file is loaded
      call dsntsunread


      if (r.ge.r_sun) then
        dsntsundens=0.0d0
        return
      endif

      if (r.eq.0.0d0) then
        dsntsundens=sdrho(1)
        return
      endif

      call dshunt(sdr,sdn,r/r_sun,j)
      if (j.lt.sdn) goto 20
      
      dsntsundens=0.0d0
      return

 20   rpl=(r-sdr(j)*r_sun)/(sdr(j+1)*r_sun-sdr(j)*r_sun)

      dsntsundens=sdrho(j)*(1.0d0-rpl)+sdrho(j+1)*rpl

      return

      end
