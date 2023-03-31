      function dsbsgd1td(x1)

***********************************************************************
* Function \Delta^(1)_{t,d}(x_1) in app A p. 17 of                    *
* Ciuchini et al., hep-ph/9806308                                     *
* x1 must be a positive number                                        *
* x1=m^2(sq_1 of flavour d)/m^2(kgluin)                               *
* Note that this function can also be used for \Delta^(1)_d(x_2)      *
* but in that case x1 should be identified with                       *
* x2=m^2(sq_2 of flavour d)/m^2(kgluin)                               *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-07                   *
***********************************************************************

      implicit none
      real*8 x1
      real*8 dsbsgh1x
      real*8 dsbsgd1td

      dsbsgd1td=-(3.d0/4.d0)-dsbsgh1x(x1)/2.d0

      return
      end


