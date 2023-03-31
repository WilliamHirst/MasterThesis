      function dsbsgc71wsusy()

***********************************************************************
* The next to leading order contribution to                           *
* the Wilson coefficient C_7 from the susy renormalization effect in  *
* the W coupling                                                      *
* Eq (23) of Ciuchini et al.,                                         *
* hep-ph/9806308                                                      *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-08                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 muw,mtw
      real*8 dsbsgmtmuw
      real*8 dsbsgg7w,dsbsgwud
      real*8 dsbsgc71wsusy,mt

c     Set the running top mass
c     uses the function dsbsgmtmuw with input m

      mt=mtmt ! mt mass at scale mt, value from DarkSUSY
      muw=mt

      mtw=dsbsgmtmuw(muw)


      dsbsgc71wsusy=(4.d0/3.d0)*(dsbsgwud(3,3)+dsbsgwud(3,2))
     &              *dsbsgg7w(mtw**2/mass(kw)**2)
     &             -(23.d0/27.d0)*(dsbsgwud(2,3)+dsbsgwud(2,2))

      return
      end



