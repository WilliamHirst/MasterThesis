      subroutine dsrdwdwdcos(p,n)

****************************************************************
** write out a table of dsandwdcos as a function of costheta for  **
** the given p and with n number of steps (n+1 points)        **
****************************************************************

      implicit none
      real*8 p,dsandwdcos,del,cth
      integer n,i
      include 'dsidtag.h'

c---------------------------------------------------------------
      if (n.gt.8192) n=8192
      open (unit=13,file='dsandwdcos.dat',status='unknown',
     &  form='formatted')
      write(13,*) '# p = ',p
      write(13,*) '# n = ',n
      write(13,*) '# model: ',idtag
      del = 2.0d0/dble(n)
      cth=-1.0d0
      write(13,*) cth,dsandwdcos(p,cth)
      do i=1,n
        cth=min(cth+del,1.0d0)
        write(13,*) cth,dsandwdcos(p,cth)
      enddo

      close(13)

      return
      end

