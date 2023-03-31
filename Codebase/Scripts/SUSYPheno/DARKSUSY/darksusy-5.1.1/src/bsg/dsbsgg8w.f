      function dsbsgg8w(x)

***********************************************************************
* Function G^W_8(x) in eq. (27) of                                    *
* Ciuchini et al., hep-ph/9806308                                     *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-08                   *
***********************************************************************

      implicit none
      real*8 x
      real*8 dsbsgg8w

      dsbsgg8w=-(4.d0-5.d0*x+5.d0*x**2)/(12.d0*(x-1.d0)**3)
     &   +x*(1.d0-2.d0*x)*log(x)/(2.d0*(x-1.d0)**4)


      return
      end


