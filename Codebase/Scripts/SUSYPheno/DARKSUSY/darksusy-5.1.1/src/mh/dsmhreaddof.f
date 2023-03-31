************************************************************************
*** read table of effective degrees of freedom in the early universe
***
***     file: name of file containing table of dof
***
*** author: Torsten Bringmann (troms@physto.se), 2010-01-23
***         (modified version of paolo gondolos dsrdreaddof) 
************************************************************************

      subroutine dsmhreaddof(filename)

      implicit none
      character*(*) filename
      character*300 msg,scr
      include 'dsmhcom.h'
      integer i

c      write (msg,'(a,a)') 'dsmhreaddof: reading ',filename
c      call dswrite(0,0,msg)
      open (unit=99,file=filename)
c... skipping first 3 lines
      read (99,*,end=2000,err=2000) scr
      read (99,*,end=2000,err=2000) scr
      read (99,*,end=2000,err=2000) scr
      i=1
 100  continue
      read (99,*,end=1000,err=2000) tmev(i),sqrtgtab(i),sqrtgttab(i)
      i=i+1
      goto 100
 1000 continue
      close(99)
      nf=i-1
      return
 2000 continue
      close(99)
      write (*,*) 'dsmhreaddof: error while reading ',filename
      end
