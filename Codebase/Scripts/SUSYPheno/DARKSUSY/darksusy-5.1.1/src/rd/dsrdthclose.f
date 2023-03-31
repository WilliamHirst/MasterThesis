      real*8 function dsrdthclose(p)
c_______________________________________________________________________
c  returns
c  author: joakim edsjo, edsjo@physto.se
c  date: april 30, 1998
c  modified: april 30, 1998.
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 pp1,pp2,p
      integer i
c-----------------------------------------------------------------------

      if (nth.eq.0) then
        dsrdthclose=0.0d0
      endif

      pp1=0.0d0
      do i=1,nth
        if (p.ge.pth(i)) pp1=pth(i)
      enddo

      pp2=pp(indx(nr))
      do i=nth,1,-1
        if (p.le.pth(i)) pp2=pth(i)
      enddo

      if (abs(p-pp1).le.abs(p-pp2)) then
        dsrdthclose=pp1
      else
        dsrdthclose=pp2
      endif

      return
      end










