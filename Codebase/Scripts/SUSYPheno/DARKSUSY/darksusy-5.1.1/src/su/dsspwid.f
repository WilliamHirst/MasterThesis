      subroutine dsspwid
c_______________________________________________________________________
c  Set or calcualte the widths of sparticles
c  Author: Joakim Edsjo, edsjo@physto.se
c  Date: 2008-07-02
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsprep.h'

      integer i

cc...Widths changed 020613 by Joakim Edsjoc
c      width(kcha(1)) = 1.0d0
c      width(kcha(2)) = 1.0d0
c      width(kn(1)) = 1.0d0
c      width(kn(2)) = 1.0d0
c      width(kn(3)) = 1.0d0
c      width(kn(4)) = 1.0d0
c      do i=1,6
c        width(ksl(i))=0.2d0
c        width(ksnu(i))=0.2d0
c        width(ksqu(i))=0.2d0
c        width(ksqd(i))=0.2d0
c      enddo

c...Widths changed 050907 by Joakim Edsjo
      width(kcha(1)) = mass(kcha(1))*0.005d0
      width(kcha(2)) = mass(kcha(2))*0.005d0
      width(kn(1)) = mass(kn(1))*0.005d0
      width(kn(2)) = mass(kn(2))*0.005d0
      width(kn(3)) = mass(kn(3))*0.005d0
      width(kn(4)) = mass(kn(4))*0.005d0

      do i=1,6
        width(ksl(i))=mass(ksl(i))*0.005d0
        width(ksnu(i))=mass(ksnu(i))*0.005d0
        width(ksqu(i))=mass(ksqu(i))*0.005d0
        width(ksqd(i))=mass(ksqd(i))*0.005d0
      enddo

      end


