      subroutine dsrdwintprint(unit)
      implicit none
***********************************************************************
*** Print out a the table of invariant rate with the points
*** that are tabulated. The output is printed to unit.
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2005-10-31
***********************************************************************


      include 'dsrdcom.h'
      real*8 dsrdwintp
      external dsrdwintp
      integer unit,i
c------------------------------------------------------------------------
      write(unit,*) '#......p..... ......w.....'
      do i=1,nr
        write(unit,101) pp(indx(i)),dsrdwintp(pp(indx(i)))
 101    format(2(1x,E12.6))
      enddo

      end








