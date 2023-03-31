***********************************************************************
*** dsntsunne gives the number density of electrons as a function
*** of the Sun's radius.
*** Input: solar radius [m]
*** Output: n_e [cm^-3]
*** See dsntsunread for information about which solar model is used.
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2006-03-27
***********************************************************************

      real*8 function dsntsunne(r)
      implicit none

      include 'dssun.h'
      include 'dsmpconst.h'

      real*8 r,rpl
      integer j

c...Check if data file is loaded
      call dsntsunread

      if (r.ge.r_sun*sdrne(sdnne)) then
        dsntsunne=0.0d0
        return
      endif

      if (r.eq.0.0d0) then
        dsntsunne=10**sdne(1)*n_avogadro
        return
      endif

      call dshunt(sdrne,sdnne,r/r_sun,j)
      if (j.lt.sdnne) goto 20
      
      dsntsunne=0.0d0
      return

 20   rpl=(r-sdrne(j)*r_sun)/(sdrne(j+1)*r_sun-sdrne(j)*r_sun)

      dsntsunne=sdne(j)*(1.0d0-rpl)+sdne(j+1)*rpl
c...Tabulated is log(n_e/n_avogadro). Take away log and n_avogadro
      dsntsunne=10**(dsntsunne)*n_avogadro

      return

      end
