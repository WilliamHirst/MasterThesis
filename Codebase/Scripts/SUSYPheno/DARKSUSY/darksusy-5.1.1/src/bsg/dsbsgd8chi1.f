      function dsbsgd8chi1(x)

***********************************************************************
* Function \Delta_8^{\chi,1} eq (21) of                               *
* Ciuchini et al., hep-ph/9806308                                     *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-08                   *
***********************************************************************

      implicit none
      real*8 x
      real*8 dsbsgf81
      real*8 dsbsgd8chi1

      dsbsgd8chi1=-(28.d0/9.d0)*dsbsgf81(x)


      return
      end


