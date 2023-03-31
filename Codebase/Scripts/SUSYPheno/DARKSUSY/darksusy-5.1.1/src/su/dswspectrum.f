************************************************************************
      subroutine dswspectrum(unit)
c_______________________________________________________________________
c  write out a table of the mass spectrum.
c  input:
c    unit - logical unit to write on (integer)
c  common:
c    'dsmssm.h' - file with susy common blocks
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c=======================================================================
      implicit none
      include 'dsmssm.h'
      integer unit,i,j
      if (unit.le.0) return
      do i=1,16
         write (unit,1000) pname(i),mass(i),width(i)
      enddo
      write (unit,1000) pname(kh1),mass(kh1),width(kh1)
      write (unit,1000) pname(kh2),mass(kh2),width(kh2)
      write (unit,1001) cos(alpha),sin(alpha)
      write (unit,1001) -sin(alpha),cos(alpha)
      write (unit,1000) pname(kh3),mass(kh3),width(kh3)
      write (unit,1000) pname(khc),mass(khc),width(khc)
      write (unit,1000) pname(kgluin),mass(kgluin),width(kgluin)
      write (unit,1000)
     &     (pname(kn(i)),mass(kn(i)),width(kn(i)),i=1,4)
      do i=1,4
         write (unit,1002)
     &        (dimag(neunmx(i,j)),dreal(neunmx(i,j)),j=1,4)
      enddo
      write (unit,1000)
     &     (pname(kcha(i)),mass(kcha(i)),width(kcha(i)),i=1,2)
      do i=1,2
         write (unit,1003)
     &        (dimag(chaumx(i,j)),dreal(chaumx(i,j)),j=1,2)
      enddo
      do i=1,2
         write (unit,1004)
     &        (dimag(chavmx(i,j)),dreal(chavmx(i,j)),j=1,2)
      enddo
      write (unit,1000)
     &     (pname(ksnu(i)),mass(ksnu(i)),width(ksnu(i)),i=1,3)
      do i=1,3
         write (unit,1005)
     &        (dimag(slulmx(i,j)),dreal(slulmx(i,j)),j=1,3)
      enddo
      write (unit,1000)
     &     (pname(ksl(i)),mass(ksl(i)),width(ksl(i)),i=1,6)
      do i=1,6
         write (unit,1006)
     &        (dimag(sldlmx(i,j)),dreal(sldlmx(i,j)),j=1,3)
      enddo
      do i=1,6
         write (unit,1007)
     &        (dimag(sldrmx(i,j)),dreal(sldrmx(i,j)),j=1,3)
      enddo
      write (unit,1000)
     &     (pname(ksqu(i)),mass(ksqu(i)),width(ksqu(i)),i=1,6)
      do i=1,6
         write (unit,1006)
     &        (dimag(squlmx(i,j)),dreal(squlmx(i,j)),j=1,3)
      enddo
      do i=1,6
         write (unit,1007)
     &        (dimag(squrmx(i,j)),dreal(squrmx(i,j)),j=1,3)
      enddo
      write (unit,1000)
     &     (pname(ksqd(i)),mass(ksqd(i)),width(ksqd(i)),i=1,6)
      do i=1,6
         write (unit,1006)
     &        (dimag(sqdlmx(i,j)),dreal(sqdlmx(i,j)),j=1,3)
      enddo
      do i=1,6
         write (unit,1007)
     &        (dimag(sqdrmx(i,j)),dreal(sqdrmx(i,j)),j=1,3)
      enddo
      write (unit,*)
      return
 1000 format (2x,a8,1x,f12.6,1x,f11.6)
 1001 format (sp,7x,'o',2(1x,f8.5))
 1002 format (sp,7x,'n',4(1x,f8.5,'i',f8.5))
 1003 format (sp,7x,'u',2(1x,f8.5,'i',f8.5))
 1004 format (sp,7x,'v',2(1x,f8.5,'i',f8.5))
 1005 format (sp,7x,'l',3(1x,f8.5,'i',f8.5))
 1006 format (sp,7x,'l',6(1x,f8.5,'i',f8.5))
 1007 format (sp,7x,'r',6(1x,f8.5,'i',f8.5))
      end
