      function dsg0loop(qsq,m1sq,m2sq)
c     a loop function
      implicit none
      real*8 dsg0loop,qsq,m1sq,m2sq
      if (m1sq-m2sq.eq.0.0d0) then
         dsg0loop = dlog(m1sq*m2sq/qsq**2)
      else
         dsg0loop = dlog(m1sq*m2sq/qsq**2) - 2.d0 +
     &        (m1sq+m2sq)/(m1sq-m2sq)*dlog(m1sq/m2sq)
      endif
      return
      end
