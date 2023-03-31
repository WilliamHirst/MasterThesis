      logical function dsrdthtest(i)
c_______________________________________________________________________
c  routine to check if the momentum p with index i is below a
c     threshold and p with index i+1 is above it. if that is the case,
c     .true. is returned, otherwise .false.
c  author: joakim edsjo, edsjo@physto.se
c  date: april 30, 1998
c  modified: april 30, 1998.
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 pp1,pp2
      integer i,j
      logical tmp
c-----------------------------------------------------------------------
      pp1=pp(indx(i))
      pp2=pp(indx(i+1))

      tmp=.false.

      do j=1,nth
        if (pp1.le.pth(j).and.pp2.gt.pth(j)) then
          tmp=.true.
c          write(*,*) pp1,pp2,pth(j)
        endif
      enddo

      dsrdthtest=tmp

      return
      end










