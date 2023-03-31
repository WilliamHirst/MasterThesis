***********************************************************************
*** This routine uses a derived potential from the BP2000 model
*** The data in sdphi() is calculated by dsntsunread.f. 
***********************************************************************

***********************************************************************
*** dsntsunpot gives the potential in the Sun as a function of radius
*** the radius should be given in m and the potential is returned in
*** m^2 s^-2
*** Density and element mass fractions up to O16 are from the standard
*** solar model BP2000 of Bahcall, Pinsonneault and Basu,
*** ApJ 555 (2001) 990.
*** The mass fractions for heavier elements are from N. Grevesse and
*** A.J. Sauval, Space Science Reviews 85 (1998) 161 normalized such that
*** their total mass fractions matches that of the heavier elements in 
*** the BP2000 model.
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2003-11-26
***********************************************************************

      real*8 function dsntsunpot(r)
      implicit none

      include 'dssun.h'
      real*8 r,rpl,gn
      integer i,j
      parameter(gn=6.67259d-11) ! m^3 kg^-1 s^-1

c...Check if data file is loaded
      call dsntsunread

      if (r.ge.r_sun) then
        dsntsunpot=-m_sun*gn/r
        return
      endif

      if (r.le.0.0d0) then
        dsntsunpot=sdphi(1)
        return
      endif

c...Interpolate in table

      call dshunt(sdr,sdn,r/r_sun,j)
      if (j.lt.sdn) goto 20

      dsntsunpot=sdphi(sdn)
      return

 20   rpl=(r-sdr(j)*r_sun)/(sdr(j+1)*r_sun-sdr(j)*r_sun)

      dsntsunpot=sdphi(j)*(1.0d0-rpl)+sdphi(j+1)*rpl

      return

      end
