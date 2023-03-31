      subroutine dsanwriteaa
c_______________________________________________________________________
c  write out the amplitude matrix
c  author: joakim edsjo (edsjo@physto.se) 95-10-25
c  called by: different routines during debugging
c=======================================================================
      implicit none
      include 'dsandiacom.h'
      integer i,j,k,l

      write(*,*)
      write(*,*) 'the non-zero components of the amplitude',
     &  ' matrix are:'
      do i=0,1
        do j=-1,1
          do k=-1,1
            do l=-1,1
              if (aa(i,j,k,l).ne.cmplx(0.d0,0.d0)) then
                write(*,*) 'aa(',i,',',j,',',k,',',l,') = ',
     &            aa(i,j,k,l)
              endif
            enddo
          enddo
        enddo
      enddo

 1000 format(' ',4(a,i2),a)
      end
