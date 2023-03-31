      function dsbsgd8chi2(x)

***********************************************************************
* Function \Delta_8^{\chi,2} eq (21) of                               *
* Ciuchini et al., hep-ph/9806308                                     *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-08                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 x
      real*8 dsbsgf83
      real*8 dsbsgd8chi2

      dsbsgd8chi2=-(14.d0/3.d0)*dsbsgf83(x)


      return
      end


