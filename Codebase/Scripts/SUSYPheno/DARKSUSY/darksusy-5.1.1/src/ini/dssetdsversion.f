      subroutine dssetdsversion
      implicit none
      integer i
      character*300 name
      include 'dsdirver.h'
      call dsvername(name)
      i=1
 10   if (name(i:i).eq.'\000') goto 20
      i=i+1
      goto 10
 20   i=i-1
c      write (*,*) 'i=',i
      dsversion=name(:i)
c      write (*,*) dsversion
      end
