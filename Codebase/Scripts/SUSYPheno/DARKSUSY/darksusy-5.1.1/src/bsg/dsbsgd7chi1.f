      function dsbsgd7chi1(x)

***********************************************************************
* Function \Delta_7^{\chi,1} eq (20) of                               *
* Ciuchini et al., hep-ph/9806308                                     *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-08                   *
***********************************************************************

      implicit none
      real*8 x
      real*8 dsbsgf71,dsbsgf81
      real*8 dsbsgd7chi1

      dsbsgd7chi1=(32.d0/27.d0)*(dsbsgf81(x)
     &        -3.d0*dsbsgf71(x))


      return
      end


