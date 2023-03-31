      function dsbsgh1x(x)

***********************************************************************
* Function H_1(x) in app A p. 16 of Ciuchini et al., hep-ph/9806308   *
* x must be a positive number                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-07                   *
***********************************************************************

      implicit none
      real*8 x
      real*8 dsbsgh1x


      dsbsgh1x=1.d0/(1.d0-x)
     &         +(2.d0*x-x**2)*log(x)/(1.d0-x)**2

      return
      end


