      subroutine dssetdsinstall
      implicit none
      integer i
      character*300 name
      include 'dsdirver.h'
      call dsdirname(name)
      i=1
 10   if (name(i:i).eq.'\000') goto 20
      i=i+1
      goto 10
 20   i=i-1
c      write (*,*) 'i=',i
cccc      dsinstall=name(:i)//'/'
      write (*,*) 'dssetdsinstall: cannot reset dsinstall'
      stop
c      write (*,*) dsinstall
      end
