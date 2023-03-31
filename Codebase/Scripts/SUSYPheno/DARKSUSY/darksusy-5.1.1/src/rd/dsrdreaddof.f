      subroutine dsrdreaddof(filename)
c_______________________________________________________________________
c     read table of effective degrees of freedom in the early universe
c     input:
c       filename - character*(*) - name of file containing table of dof
c     common:
c       'dsrdcom.h' - included common blocks
c  author: paolo gondolo (paolo@physics.utah.edu) 2005
c=======================================================================
      implicit none
      character*(*) filename
      character*300 msg
      character*200 scr
      include 'dsrdcom.h'
      integer i
c-----------------------------------------------------------------------
      write (msg,'(a,a)') 'dsrdreaddof: reading ',filename
      call dswrite(0,0,msg)
      open (unit=99,file=filename)
      read (99,'(A)',err=2000) scr
      read (99,'(A)',err=2000) scr
      read (99,'(A)',err=2000) scr
      i=1
 100  continue
c short version without fe fp
c      read (99,*,end=1000,err=2000) tgev(i),fg(i),fh(i),fe(i),fp(i)
      read (99,*,end=1000,err=2000) tgev(i),fg(i),fh(i)
      i=i+1
      goto 100
 1000 continue
      close(99)
      nf=i-1
      return
 2000 continue
      close(99)
      write (*,*) 'dsrdreaddof: error while reading ',filename
      end
