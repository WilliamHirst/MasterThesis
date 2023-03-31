      function dsbsgechi(y)

***********************************************************************
* Function E_\chi(y) in (11) of Ciuchini et al.,                      *
* hep-ph/9806308                                                      *
* y must be a positive number                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-04                   *
***********************************************************************

      implicit none
      real*8 y
      real*8 dsbsgechi

      dsbsgechi=y*(11.d0-7.d0*y+2.d0*y**2)/(18.d0*(y-1.d0)**3)
     &         -log(y)*y/(3.d0*(y-1.d0)**4)

      return
      end


