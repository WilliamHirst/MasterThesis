      function dsbsgyt(m)

***********************************************************************
* Function that calculates y_t(m=mu_susy), the top Yukawa coupling    *
* at the scale m=mu_susy                                              *
* Note that m=mu_susy should be of the order of 1TeV, e.g. m_gluino   *
* Uses eq (23) of Degrassi et al.,  hep-ph/0009337                    *
* Note that y_t(susy)=\tilde{y}_t(susy) (see text in the ref.)        *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-03                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 m,mtpole
      real*8 dsbsgalpha3,dsbsgmtmuw
      real*8 dsbsgyt
      real*8 pi
      parameter (pi=3.141592653589793238d0)

c     set m_t to the top pole mass
c     use the value from DarkSUSY

      mtpole=mass(kt)

c     Note that in the following we use m_W
c     for the weak mass scale mu_w

c     Also note that the input scale m is used as
c     input for the functions of alpha3 and m_t

      dsbsgyt=(dsbsgmtmuw(mass(kw))/mass(khc))**2
     & *(dsbsgalpha3(m)/dsbsgalpha3(mtpole))**(4.d0/7.d0)
     & *(dsbsgalpha3(mtpole)
     &  /dsbsgalpha3(mass(kw)))**(12.d0/23.d0)
     & /sqrt(1.d0+
     &  (9.d0*(dsbsgmtmuw(mtpole)/mass(khc))**4
     &   /(8.d0*pi*dsbsgalpha3(mtpole)))
     &  *((dsbsgalpha3(m)/dsbsgalpha3(mtpole))**(1.d0/7.d0)
     &    -1.d0))

      return
      end


