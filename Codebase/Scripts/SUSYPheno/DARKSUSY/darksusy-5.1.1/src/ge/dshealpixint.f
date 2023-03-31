      subroutine dshealpixint(f,nmin,nmax,eps,eps1,ans,niter,ierr)
      implicit none
      real*8 f,ans
      external f
      integer ierr,niter
      real*8 aveo,ave,eps,eps1,FourPI
      integer n,nmin,nmax
c      parameter (nmin=4,nmax=400)
c      parameter (eps=1.d-4,eps1=1.d-4)
      FourPI = 16.0 * datan(1.d0)
      call dshealpixave(nmin,f,aveo)
      do n=nmin+1,nmax
         call dshealpixave(n,f,ave)
         if (abs(ave-aveo).lt.eps*abs(ave)+eps1) then
            ans=ave*FourPI
            niter=n
            return
         endif
         aveo=ave
      enddo
      niter=n
      ans=ave*FourPI
      write (*,*) 'dshealpixint: max number of sides reached ',nmax
      ierr=1
      return
      end
