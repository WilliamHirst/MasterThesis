      function dsbsgf72(y)

***********************************************************************
* Function F_7^(2)(y) in (54) of Ciuchini et al.,                     *
* hep-ph/9710335                                                      *
* y must be a positive number                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-03-31                   *
***********************************************************************

      implicit none
      real*8 y
      real*8 dsbsgf72

      dsbsgf72=y*(3.d0-5.d0*y)/(12.d0*(y-1.d0)**2)
     &         +log(y)*y*(3.d0*y-2.d0)/(6.d0*(y-1.d0)**3)

      return
      end


