      function dsbsgh3x(x)

***********************************************************************
* Function H_3(x) in app A p. 17 of Ciuchini et al., hep-ph/9806308   *
* x must be a positive number                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-07                   *
***********************************************************************

      implicit none
      real*8 x
      real*8 dsbsgh3x

      dsbsgh3x=2.d0*x*log(x)/(1.d0-x)

      return
      end


