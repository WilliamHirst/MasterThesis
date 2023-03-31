      function dsbsgg8h(y)

***********************************************************************
* Function G_8^H(y) in (62) of Ciuchini et al.,                       *
* hep-ph/9710335                                                      *
* y must be a positive number                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-03-31                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 y,au,ad
      real*8 ddilog
      real*8 dsbsgg8h

      au=1.d0/tanbe
      ad=-tanbe

      dsbsgg8h=ad*au*(1.d0/3.d0)*y*(
     &         (-36.d0+25.d0*y-17.d0*y**2)*
     &         ddilog(1.d0-1.d0/y)/(2.d0*(y-1.d0)**3)
     &        +(19.d0+17.d0*y)*
     &         (log(y))**2/(y-1.d0)**4
     &        +(-3.d0-187.d0*y+12.d0*y**2-14.d0*y**3)*
     &         log(y)/(4.d0*(y-1.d0)**4)
     &        +3.d0*(143.d0-44.d0*y+29.d0*y**2)/(8.d0*(y-1.d0)**3))
     &       +au**2*(1.d0/6.d0)*y*(
     &         y*(30.d0-17.d0*y+13.d0*y**2)
     &         *ddilog(1.d0-1.d0/y)/(y-1.d0)**4
     &         -y*(31.d0+17.d0*y)*
     &          (log(y))**2/(y-1.d0)**5
     &        +(-226.d0+817.d0*y+1353.d0*y**2+318.d0*y**3+42*y**4)
     &         *log(y)/(36.d0*(y-1.d0)**5)
     &        +(1130.d0-18153.d0*y+7650.d0*y**2-4451*y**3)/
     &         (216.d0*(y-1.d0)**4))

      return
      end


