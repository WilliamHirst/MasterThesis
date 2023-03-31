      subroutine dswexcl(unit,excl)
c_______________________________________________________________________
c  write reasons for exclusion to specified unit.
c  input:
c    unit - logical unit to write on (integer)
c    excl - code of the reason for exclusion (integer); 0 if allowed
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c=======================================================================
      implicit none
      integer unit,excl
      character*12 dsidtag
      if (unit.le.0) return
      if (excl.eq.0) return
      write (unit,1001) dsidtag(),' excluded by: '
      if (btest(excl,0)) write (unit,1000) 'chargino mass; '
      if (btest(excl,1)) write (unit,1000) 'gluino mass; '
      if (btest(excl,2)) write (unit,1000) 'squark mass; '
      if (btest(excl,3)) write (unit,1000) 'slepton mass; '
      if (btest(excl,4)) write (unit,1000) 'gamma_z(inv); '
      if (btest(excl,5)) write (unit,1000) 'h2 mass; '
      if (btest(excl,6)) write (unit,1000) 'neutralino mass; '
      if (btest(excl,7)) write (unit,1000) 'b->s+gamma; '
      if (btest(excl,8)) write (unit,1000) 'delta rho; '
      if (btest(excl,9)) write (unit,1000) '(g-2)_mu; '
      write (unit,*)
 1000 format (1x,a,$)
 1001 format (1x,'(',a,')',1x,a,$)
      end
