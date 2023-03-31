************************************************************************
      subroutine dswvertx(unit)
c_______________________________________________________________________
c  write out a table of the vertices constants.
c  input:
c    unit - logical unit to write on (integer)
c  common:
c    'dsmssm.h' - file with susy common blocks
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c=======================================================================
      implicit none
      include 'dsmssm.h'
      integer unit,i,j,k,l
      complex*16 g,g4p
      if (unit.le.0) return
      do i=1,50
         do j=i,50
            do k=j,50
               if (dreal(gl(i,j,k)).ne.0.d0.or.
     &              dimag(gl(i,j,k)).ne.0.d0.or.
     &              dreal(gr(i,j,k)).ne.0.d0.or.
     &              dimag(gr(i,j,k)).ne.0.d0) then
                  if ((13.le.i.and.i.le.20.or.i.eq.49.or.i.eq.50).and.
     &                 (13.le.j.and.j.le.20.or.j.eq.49.or.j.eq.50).and.
     &                 (13.le.k.and.k.le.20.or.k.eq.49.or.k.eq.50)) then
                                ! gauge or higgs
                     write (unit,1000) pname(i),pname(j),pname(k),
     &                    dimag(gl(i,j,k)),dreal(gl(i,j,k))
                  else
                     write (unit,2000) pname(i),pname(j),pname(k),
     &                    dimag(gl(i,j,k)),dreal(gl(i,j,k)),
     &                    dimag(gr(i,j,k)),dreal(gr(i,j,k))
                  endif
               endif
            enddo
         enddo
      enddo
      do i=1,50
         do j=i,50
            do k=j,50
               do l=k,50
                  g=g4p(i,j,k,l)
                  if (dreal(g).ne.0.d0.or.
     &                 dimag(g).ne.0.d0) then
                     write (unit,3000) pname(i),pname(j),pname(k),
     &                    pname(l),
     &                    dimag(g),dreal(g)
                  endif
               enddo
            enddo
         enddo
      enddo
 1000 format (1x,a8,1x,a8,1x,a8,3x,sp,f11.5,'i',f11.5)
 2000 format (1x,a8,1x,a8,1x,a8,2x,
     &     sp,'(',f11.5,'i',f11.5,')pl + (',f11.5,'i',f11.5,')pr')
 3000 format (1x,a8,1x,a8,1x,a8,1x,a8,3x,sp,f11.5,'i',f11.5)
      end
