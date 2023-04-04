      function dsbsgc71phi1susy()

***********************************************************************
* The next to leading order contribution to                           *
* the Wilson coefficient C_7 from the susy renormalization effect in  *
* the charged scalar,\phi_1, coupling                                 *
* Eq (25) of Ciuchini et al.,                                         *
* hep-ph/9806308                                                      *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-08                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 muw,mtw
      real*8 dsbsgmtmuw
      real*8 dsbsgf71,dsbsgf72,dsbsghtd,dsbsghd
      real*8 dsbsgc71phi1susy,mt

c     Set the running top mass
c     uses the function dsbsgmtmuw with input m

      mt=mtmt  ! mt at scale mt, value from DarkSUSY
      muw=mt

      mtw=dsbsgmtmuw(muw)


      dsbsgc71phi1susy=(4.d0/(9.d0*tanbe**2))
     &        *(dsbsghtd(2)+dsbsghtd(3))
     &        *dsbsgf71(mtw**2/mass(khc)**2)
     &      +(4.d0/3.d0)*(dsbsghtd(2)+dsbsghd())
     &        *dsbsgf72(mtw**2/mass(khc)**2)

      return
      end


