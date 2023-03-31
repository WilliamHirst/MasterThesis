      function dshlf2(x,y)
c paolo gondolo
      implicit none
      real*8 dshlf2,x,y
      if (x-y.eq.0.d0) then
         dshlf2 = 0.d0
      else
         dshlf2 = 2.d0 - (x+y)/(x-y)*log(x/y)
      endif
      return
      end
