      function dsbsgd8h(y)

***********************************************************************
* Function \Delta_8^H(y) in (63) of Ciuchini et al.,                  *
* hep-ph/9710335                                                      *
* y must be a positive number                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-03-31                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 y,au,ad
      real*8 dsbsgd8h

      au=1.d0/tanbe
      ad=-tanbe

      dsbsgd8h=ad*au*(1.d0/3.d0)*y*(
     &         (81.d0-16.d0*y+7.d0*y**2)/(2.d0*(y-1.d0)**3)
     &         -(19.d0+17.d0*y)*
     &         log(y)/(y-1.d0)**4)
     &       +au**2*(1.d0/6.d0)*y*(
     &         (-38.d0-261.d0*y+18.d0*y**2-7.d0*y**3)
     &         /(6.d0*(y-1.d0)**4)
     &         +y*(31.d0+17.d0*y)*
     &          log(y)/(y-1.d0)**5)

      return
      end


