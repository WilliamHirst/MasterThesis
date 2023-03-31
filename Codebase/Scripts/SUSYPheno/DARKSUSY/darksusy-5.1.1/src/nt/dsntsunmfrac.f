***********************************************************************
*** dsntsunmfrac gives the mass fraction of element i (see dsntsunread.f
*** for definition of i) as a function of the solar radius r.
*** the radius should be given in m and returned is the mass fraction.
***
*** Element mass fractions up to O16 are from the standard
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

      real*8 function dsntsunmfrac(r,itype)
      implicit none

      include 'dssun.h'
      real*8 r,rpl
      integer i,j,itype

c...Check if data file is loaded
      call dsntsunread

      if (r.gt.r_sun) then
        dsntsunmfrac=0.0d0
        return
      endif

      if (r.eq.0.0d0) then
        dsntsunmfrac=sdmfr(itype,1)
        return
      endif

      call dshunt(sdr,sdn,r/r_sun,j)
      if (j.lt.sdn) goto 20

      j=sdn-1

 20   rpl=(r-sdr(j)*r_sun)/(sdr(j+1)*r_sun-sdr(j)*r_sun)

      dsntsunmfrac=sdmfr(itype,j)*(1.0d0-rpl)
     &  +sdmfr(itype,j+1)*rpl

      return

      end
