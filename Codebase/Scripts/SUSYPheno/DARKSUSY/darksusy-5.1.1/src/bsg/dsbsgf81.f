      function dsbsgf81(y)

***********************************************************************
* Function F_8^(1)(y) in (30) of Ciuchini et al.,                     *
* hep-ph/9710335                                                      *
* y must be a positive number                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-03-31                   *
***********************************************************************

      implicit none
      real*8 y
      real*8 dsbsgf81

      dsbsgf81=y*(2.d0+5.d0*y-y**2)/(8.d0*(y-1.d0)**3)
     &         -log(y)*3.d0*y**2/(4.d0*(y-1.d0)**4)

      return
      end


