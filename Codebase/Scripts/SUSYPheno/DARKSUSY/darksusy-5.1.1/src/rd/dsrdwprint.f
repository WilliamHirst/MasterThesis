      subroutine dsrdwprint(unit,np,wrate,p_min,p_max)
      implicit none
***********************************************************************
*** Print out a the table of invariant rate starting at p_min, ending
*** at p_max and with np+1 number of points. The rate routine called is
*** wrate. The output is printed to unit.
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2005-10-31
***********************************************************************


      include 'dsrdcom.h'
      real*8 wrate,p_min,p_max,p,ww
      external wrate
      integer unit,np,i
c------------------------------------------------------------------------
      write(unit,*) '#......p..... ......w.....'
      do i=0,np
        p=p_min+(p_max-p_min)*dble(i)/dble(np)
        ww=wrate(p)
        write(unit,101) p,ww
 101    format(2(1x,E12.6))
      enddo

      end








