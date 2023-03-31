      function dsbsgf82(y)

***********************************************************************
* Function F_8^(2)(y) in (55) of Ciuchini et al.,                     *
* hep-ph/9710335                                                      *
* y must be a positive number                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-03-31                   *
***********************************************************************

      implicit none
      real*8 y
      real*8 dsbsgf82

      dsbsgf82=y*(3.d0-y)/(4.d0*(y-1.d0)**2)
     &         -log(y)*y/(2.d0*(y-1.d0)**3)

      return
      end


