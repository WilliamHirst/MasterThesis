      function dsbsgmtmuwint(mtstart,mstart,m,nf)

***********************************************************************
* The running top mass from value mtstart at scale mstart.            *
* The running is done with nf effective active quark flavours.        *
* Uses eq. (32) of Ciuchini et al. hep-ph/9710335                     *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-03                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      integer nf
      real*8 m
      real*8 mt,mtpole
      real*8 mtstart,mstart
      real*8 b0,b1,g0m,g1m,r 
      real*8 dsbsgalpha3
      real*8 dsbsgmtmuwint
      real*8 pi
      parameter (pi=3.141592653589793238d0)

c     set the value of the running mass m_t(m_t) 
c     Use the DarkSUSY value

      mt=mtmt

c     set the top pole mass
c     Use DarkSUSY value

      mtpole=mass(kt)

c     define the numbers beta_0 and beta_1 from eq. (33)
c     depend on the input number n_f of eff. quark flavours

      g0m=8.d0

      if (nf.le.6) then

        b0=11.d0-2.d0*dble(nf)/3.d0   
        b1=102.d0-38.d0*dble(nf)/3.d0
c     define the numbers gamma_0^m and gamma_1^m from eq.(34)
        g1m=404.d0/3.d0-40.d0*dble(nf)/9.d0

      else

        b0=41.0d0/6.0d0                ! UPDATE what should this be?
        b1=102.d0-38.d0*dble(6)/3.d0   ! UPDATE what should this be?
c     define the numbers gamma_0^m and gamma_1^m from eq.(34)
        g1m=404.d0/3.d0-40.d0*dble(nf)/9.d0  ! UPDATE what should this be?

      endif

c     define a ratio needed below

      r=g0m/(2.d0*b0)

      dsbsgmtmuwint=mtstart*(dsbsgalpha3(m)/dsbsgalpha3(mstart))**r
     &  *(1.d0+(dsbsgalpha3(mstart)/(4.d0*pi))*r*
     &    (g1m/g0m-b1/b0)*
     &    (dsbsgalpha3(m)/dsbsgalpha3(mstart)-1.d0))

      return
      end


