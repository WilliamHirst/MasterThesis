      function dshlf3(p2,y1,y2)
c paolo gondolo
      implicit none
      real*8 dshlf3,p2,y1,y2,delta
      complex*16 r,l,ans
      delta = (y1-y2)/p2
      r = sqrt(dcmplx((1.d0+delta)**2-4.d0*y1/p2,0.d0))
      l = log(((1.d0+r)**2-delta**2)/((1.d0-r)**2-delta**2))
      if (p2.gt.y1+y2) l = l + dcmplx(0.d0,-2.d0*3.141592653589793238d0)
      if (y2-y1.eq.0.d0) then
         ans = -2.d0+0.5d0*r*l
      else
         ans = -1.d0+0.5d0*((y1+y2)/(y1-y2)-delta)*dlog(y2/y1)+0.5d0*r*l
      endif
      dshlf3=dreal(ans)
      return
      end
