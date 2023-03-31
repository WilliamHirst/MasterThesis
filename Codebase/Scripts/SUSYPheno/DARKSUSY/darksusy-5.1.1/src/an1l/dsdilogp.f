c====================================================================
c
c   auxiliary function used in:
c   dsdilog.f
c
c   author: lars bergstrom (lbe@physto.se)
c
c____________________________________________________________________

      real*8 function dsdilogp(x)
      implicit none
      real*8 x, term, sum, psum
      integer i
      sum=0.d0
      psum=0.d0
      term=1.d0
      do 10 i=1,200
      term=term*x
      sum=sum+term/(i*i*(i+1)*(i+1))
      if (abs(sum-psum).le.1.d-20) goto 20
      psum=sum
 10    continue
 20    continue
      dsdilogp=x*(3.d0+sum)+2.d0*(1.d0-x)*dlog(1.d0-x)
      dsdilogp=dsdilogp/(1.d0+x)
      return
      end
