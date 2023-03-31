      function dsbsgd7chi2(x)

***********************************************************************
* Function \Delta_7^{\chi,2} eq (20) of                               *
* Ciuchini et al., hep-ph/9806308                                     *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-08                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 x
      real*8 dsbsgf73,dsbsgf83
      real*8 dsbsgd7chi2

      dsbsgd7chi2=(16.d0/9.d0)*(dsbsgf83(x)
     &        -3.d0*dsbsgf73(x))


      return
      end


