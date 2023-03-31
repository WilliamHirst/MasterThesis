      function dsbsgc81h2()  ! (muw)

***********************************************************************
* The next to leading order contribution to the Wilson coefficient C_8*
* from the two-Higgs doublet model                                    *
* Eq. (59) of Ciuchini et al.,                                        *
* hep-ph/9710335                                                      *
* for the calculation of b --> s gamma                                *
* The input parameter muw=\mu_W is the matching scale                 *
* author:Mia Schelke, schelke@physto.se, 2003-03-31                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 mt,y,muw
      real*8 dsbsgg8h,dsbsgd8h,dsbsgeh
      real*8 dsbsgc81h2


c     We here set \bar(m)_t(\mu_W)=\bar(m)_t(\bar(m)_t)
c     as these enter in Kt where the scale is mt

      mt=mtmt

      y=mt**2/mass(khc)**2

      muw=mt


      dsbsgc81h2=dsbsgg8h(y)
     &          +dsbsgd8h(y)*log(muw**2/mass(khc)**2)
     &          -(1.d0/6.d0)*dsbsgeh(y)

      return
      end


