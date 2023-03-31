      function dsbsgd7h(y)

***********************************************************************
* Function \Delta_7^H(y) in (61) of Ciuchini et al.,                  *
* hep-ph/9710335                                                      *
* y must be a positive number                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-03-31                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 y,au,ad
      real*8 dsbsgd7h

      au=1.d0/tanbe
      ad=-tanbe

      dsbsgd7h=ad*au*(2.d0/9.d0)*y*(
     &         (21.d0+47.d0*y+8.d0*y**2)/(y-1.d0)**3
     &         +2.d0*(-8.d0+14.d0*y+3.d0*y**2)*
     &         log(y)/(y-1.d0)**4)
     &       +au**2*(2.d0/9.d0)*y*(
     &         (-31.d0-18.d0*y+135.d0*y**2-14.d0*y**3)
     &         /(6.d0*(y-1.d0)**4)
     &         +y*(14.d0-23.d0*y-3.d0*y**2)*
     &          log(y)/(y-1.d0)**5)

      return
      end


