
      function dsbsgmtmuw(m)

***********************************************************************
* The running top mass evaluated at a weak scale m=mu_w               *
* using nf effective quark flavours (taken to be nf=5 here)           *
* Uses eq. (32) of Ciuchini et al. hep-ph/9710335                     *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-03                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      integer nf
      real*8 m
      real*8 mt,mtpole    ! ,mttmp
      real*8 b0,b1,g0m,g1m,r 
      real*8 dsbsgalpha3
      real*8 dsbsgmtmuw,dsbsgmtmuwint


c     set the value of the running mass m_t(m_t) 
c     Use the DarkSUSY value

      mt=mtmt

c     set the top pole mass
c     Use DarkSUSY value

      mtpole=mass(kt)

      if (m.le.mtmt) then   ! 5 active quarks

        dsbsgmtmuw=dsbsgmtmuwint(mt,mtpole,m,5)

      else       ! 6 active quarks  
c      elseif (m.gt.mtmt.and.m.le.1.5d0*mtmt) then ! 6 active quarks
c      the line above should be used if we want to use 6 quarks and 
c      one light squark at high energies
c      We don't need to include a step with 5 active quarks as we 
c      already know the value mtmt 
        
        dsbsgmtmuw=dsbsgmtmuwint(mt,mtpole,m,6)

c      elseif (m.gt.1.5d0*mtmt) then   ! 6 active quarks and one squark

c        mttmp=dsbsgmtmuwint(mt,mtpole,1.5d0*mtmt,6)
c        dsbsgmtmuw=dsbsgmtmuwint(mttmp,1.5d0*mtmt,m,7)

c      else

c        write(*,*) 'ERROR in dsbsgmtmuw: scale out of range: ',m
c        stop

      endif

      return
      end


