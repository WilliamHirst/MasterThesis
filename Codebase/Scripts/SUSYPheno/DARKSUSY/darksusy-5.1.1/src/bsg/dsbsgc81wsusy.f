      function dsbsgc81wsusy()

***********************************************************************
* The next to leading order contribution to                           *
* the Wilson coefficient C_8 from the susy renormalization effect in  *
* the W coupling                                                      *
* Eq (24) of Ciuchini et al.,                                         *
* hep-ph/9806308                                                      *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-08                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 muw,mtw
      real*8 dsbsgmtmuw
      real*8 dsbsgg8w,dsbsgwud
      real*8 dsbsgc81wsusy,mt

c     Set the running top mass
c     uses the function dsbsgmtmuw with input m

      mt=mtmt  ! mt at scale mt, value from DarkSUSY
      muw=mt

      mtw=dsbsgmtmuw(muw)


      dsbsgc81wsusy=(4.d0/3.d0)*(dsbsgwud(3,3)+dsbsgwud(3,2))
     &              *dsbsgg8w(mtw**2/mass(kw)**2)
     &             -(4.d0/9.d0)*(dsbsgwud(2,3)+dsbsgwud(2,2))

      return
      end



