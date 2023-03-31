      subroutine dsanclearaa
c_______________________________________________________________________
c  clear the amplitude matrix
c  author: joakim edsjo (edsjo@physto.se) 95-10-25
c          paolo gondolo 99-1-15 factor of 3.7 faster
c  called by: dwdcos
c=======================================================================
      implicit none
      include 'dsandiacom.h'

      integer i
      real*8 aaaa(108)
      equivalence (aa,aaaa)

      do i=1,108
         aaaa(i)=0.0d0
      enddo

      end
