      function dsbsgc70h2()

***********************************************************************
* The leading order contribution to the Wilson coefficient C_7        *
* from the two-Higgs doublet model                                    *
* Eq. (53) of Ciuchini et al.,                                        *
* hep-ph/9710335                                                      *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-03-31                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 mt,y,au,ad
      real*8 dsbsgf71,dsbsgf72
      real*8 dsbsgc70h2

      au=1.d0/tanbe
      ad=-tanbe

c     We here set \bar(m)_t(\mu_W)=\bar(m)_t(\bar(m)_t) as these enter
c     in Kt where the scale is mt

      mt=mtmt

      y=mt**2/mass(khc)**2

      dsbsgc70h2=au**2*dsbsgf71(y)/3.d0
     &          -au*ad*dsbsgf72(y)

      return
      end


