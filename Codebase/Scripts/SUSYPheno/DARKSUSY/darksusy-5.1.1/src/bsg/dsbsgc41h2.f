      function dsbsgc41h2()

***********************************************************************
* The next to leading order contribution to the Wilson coefficient C_4*
* from the two-Higgs doublet model                                    *
* Eq. (58) of Ciuchini et al.,                                        *
* hep-ph/9710335                                                      *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-03-31                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 mt,y
      real*8 dsbsgeh
      real*8 dsbsgc41h2


c     We here set \bar(m)_t(\mu_W)=\bar(m)_t(\bar(m)_t)
c     as these enter in the Kt calculcation where the scale is mt

      mt=mtmt

      y=mt**2/mass(khc)**2

      dsbsgc41h2=dsbsgeh(y)

      return
      end


