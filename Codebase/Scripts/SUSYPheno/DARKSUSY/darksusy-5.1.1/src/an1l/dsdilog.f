c====================================================================
c
c   dsdilogarithm function
c   argument should be between -1 and 1
c
c   author: lars bergstrom (lbe@physto.se)
c
c____________________________________________________________________

      real*8 function dsdilog(x)
      implicit none
      include 'dsge.h'
      real*8 x,dsdilogp
      if (abs(x).ge.1.d0) then
      write(*,*) 'dsdilog called with wrong argument : x= ',x
      stop
      endif
      if(x.le.0.d0) then
      dsdilog=-dsdilogp(-x/(1.d0-x))-((dlog(1.d0-x))**2)/2.d0
      return
      endif
      if(x.ge.0.5d0) then
      dsdilog=-dsdilogp(1.d0-x)-dlog(x)*dlog(1.d0-x)+pi**2/6.d0
      return
      endif
      dsdilog=dsdilogp(x)
      return
      end
