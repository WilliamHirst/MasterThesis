      function dsbsgg7w(x)

***********************************************************************
* Function G^W_7(x) in eq. (27) of                                    *
* Ciuchini et al., hep-ph/9806308                                     *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-08                   *
***********************************************************************

      implicit none
      real*8 x
      real*8 dsbsgg7w

      dsbsgg7w=-(23.d0-67.d0*x+50.d0*x**2)/(36.d0*(x-1.d0)**3)
     &   +x*(2.d0-7.d0*x+6.d0*x**2)*log(x)
     &    /(6.d0*(x-1.d0)**4)


      return
      end


