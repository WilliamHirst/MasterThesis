      function dsbsgg7h(y)

***********************************************************************
* Function G_7^H(y) in (60) of Ciuchini et al.,                       *
* hep-ph/9710335                                                      *
* y must be a positive number                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-03-31                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 y,au,ad
      real*8 ddilog
      real*8 dsbsgg7h

      au=1.d0/tanbe
      ad=-tanbe

      dsbsgg7h=ad*au*(4.d0/3.d0)*y*(
     &         4.d0*(-3.d0+7.d0*y-2.d0*y**2)*
     &         ddilog(1.d0-1.d0/y)/(3.d0*(y-1.d0)**3)
     &        +(8.d0-14.d0*y-3.d0*y**2)*
     &         (log(y))**2/(3.d0*(y-1.d0)**4)
     &        +2.d0*(-3.d0-y+12.d0*y**2-2.d0*y**3)*
     &         log(y)/(3.d0*(y-1.d0)**4)
     &        +(7.d0-13.d0*y+2.d0*y**2)/(y-1.d0)**3)
     &       +au**2*(2.d0/9.d0)*y*(
     &         y*(18.d0-37.d0*y+8.d0*y**2)
     &         *ddilog(1.d0-1.d0/y)/(y-1.d0)**4
     &         +y*(-14.d0+23.d0*y+3.d0*y**2)*
     &          (log(y))**2/(y-1.d0)**5
     &        +(-50.d0+251.d0*y-174.d0*y**2-192.d0*y**3+21*y**4)
     &         *log(y)/(9.d0*(y-1.d0)**5)
     &        +(797.d0-5436.d0*y+7569.d0*y**2-1202.d0*y**3)/
     &         (108.d0*(y-1.d0)**4))

      return
      end


