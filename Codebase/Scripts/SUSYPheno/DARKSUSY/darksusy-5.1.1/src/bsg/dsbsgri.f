      function dsbsgri(famd)

***********************************************************************
* Function R_i in eq. (19) of                                         *
* Ciuchini et al., hep-ph/9806308                                     *
* The expression has been extended to large tanbe, by dropping        *
* ln((mu_w)^2)/m^2(kgluin))                                           *
* as explained in Degrassi et al., hep-ph/0009337 p.11                *
* The input parameter famd is the family of the down sector           *
* it should be 1 for the first family, 2 for the 2nd and 3 for the 3rd* 
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-07                   *
***********************************************************************

      implicit none
      integer famd
      real*8 dsbsgd2td
      real*8 dsbsgri

      dsbsgri=dsbsgd2td(famd)-1.d0

      return
      end


