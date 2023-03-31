      subroutine dshealpixave(n,f,ave)
      implicit none
      integer n
      real*8 f,ave
      integer i,ppmax
      real*8 sum,l,b
      real*8 pp(12*n*n,2)
      ppmax=12*n*n
      call dshealpixpoints(n,ppmax,pp)
      sum=0.d0
      do i=1,ppmax
         l=pp(i,1)
         b=pp(i,2)         
         sum=sum+f(l,b)
      enddo
      ave=sum/real(ppmax)
      end
