      function dsbsgf71(y)

***********************************************************************
* Function F_7^(1)(y) in (29) of Ciuchini et al.,                     *
* hep-ph/9710335                                                      *
* y must be a positive number                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-03-31                   *
***********************************************************************

      implicit none
      real*8 y
      real*8 dsbsgf71

      dsbsgf71=y*(7.d0-5.d0*y-8.d0*y**2)/(24.d0*(y-1.d0)**3)
     &         +log(y)*y**2*(3.d0*y-2.d0)/(4.d0*(y-1.d0)**4)

      return
      end


