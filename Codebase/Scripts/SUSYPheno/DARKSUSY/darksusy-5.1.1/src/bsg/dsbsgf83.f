      function dsbsgf83(y)

***********************************************************************
* Function F_8^(3)(y) in (22) of Degrassi et al.,                     *
* hep-ph/0009337                                                      *
* y must be a positive number                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-03                   *
***********************************************************************

      implicit none
      real*8 y
      real*8 dsbsgf83

      dsbsgf83=(1.d0+y)/(2.d0*(y-1.d0)**2)
     &         -log(y)*y/(y-1.d0)**3

      return
      end


