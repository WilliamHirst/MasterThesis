      function dsbsgeh(y)

***********************************************************************
* Function E^H(y) in (64) of Ciuchini et al.,                         *
* hep-ph/9710335                                                      *
* y must be a positive number                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-03-31                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 y,au,ad
      real*8 dsbsgeh

      au=1.d0/tanbe
      ad=-tanbe

      dsbsgeh=au**2*(
     &         y*(16.d0-29.d0*y+7.d0*y**2)
     &         /(36.d0*(y-1.d0)**3)
     &        +y*(3.d0*y-2.d0)*
     &          log(y)/(6.d0*(y-1.d0)**4))

      return
      end


