      function dsbsgc81phi2susy()

***********************************************************************
* The next to leading order contribution to                           *
* the Wilson coefficient C_8 from the susy renormalization effect in  *
* the unphysical charged scalar,\phi_2, coupling                      *
* Eq (26) of Ciuchini et al.,                                         *
* hep-ph/9806308                                                      *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-08                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 muw,mtw
      real*8 dsbsgmtmuw
      real*8 dsbsgf81,dsbsgf82,dsbsgutd,dsbsgud
      real*8 dsbsgc81phi2susy,mt

c     Set the running top mass
c     uses the function dsbsgmtmuw with input m

      mt=mtmt  ! mt at scale mt, value from DarkSUSY
      muw=mt

      mtw=dsbsgmtmuw(muw)


      dsbsgc81phi2susy=(4.d0/9.d0)
     &        *(dsbsgutd(2)+dsbsgutd(3))
     &        *dsbsgf81(mtw**2/mass(kw)**2)
     &      -(4.d0/3.d0)*(dsbsgutd(2)+dsbsgud())
     &        *dsbsgf82(mtw**2/mass(kw)**2)

      return
      end



