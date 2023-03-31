      function dsbsgc41susy()

***********************************************************************
* The next to leading order contribution to                           *
* the Wilson coefficient C_4 from susy                                *
* Eq (10) of Ciuchini et al.,                                         *
* hep-ph/9806308                                                      *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-07                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      integer i
      real*8 thetast,st2(2)
      real*8 beta,muw,mtw,mt
      real*8 dsbsgmtmuw,dsbsgechi
      real*8 dsbsgc41susy,mst2
      real*8 pi
      parameter (pi=3.141592653589793238d0)

      beta=atan(tanbe)                        

      mt=mtmt  ! mt at scale mt, value from DarkSUSY

c     Introduce the mixing angle in the stop sector
c     as it's defined p.7 in the ref stated above

      if(dimag(squlmx(6,3)).ne.0.d0) then
       write(*,*) 'squlmx(6,3) has a complex part
     &            so thetast in dsbsgc41susy might be wrong'
      endif

      if(dimag(squrmx(6,3)).ne.0.d0) then
       write(*,*) 'squrmx(6,3) has a complex part
     &            so thetast in dsbsgc41susy might be wrong'
      endif

c...Their convention is q1= cos(theta)q_L + sin(theta)q_R
c                       q2=-sin(tehta)q_L + cos(theta)q_R
c   with m_q1>m_q2. Hence, to be fully consistent with their
c   notation, we pick the corresponding eigenvalues and angles
c   from DarkSUSY


      if (mass(kst1).ge.mass(kst2)) then   ! DarkSUSY convention matches
        do i=1,2
         st2(i)=mass(kst2)**2/mass(kcha(i))**2
        enddo
        mst2=mass(kst2)
        thetast=atan(-dble(conjg(squlmx(6,3)))
     &             /dble(conjg(squrmx(6,3))))

      else        ! DarkSUSY convention differs. Switch 1 <-> 2
        do i=1,2
         st2(i)=mass(kst1)**2/mass(kcha(i))**2
        enddo
        mst2=mass(kst1)
        thetast=atan(-dble(conjg(squlmx(6,3)))
     &             /dble(conjg(squrmx(6,3))))+pi/2.0d0
      endif

c     Set the running top mass
c     uses the function dsbsgmtmuw with input m

      muw=mt

      mtw=dsbsgmtmuw(muw)

      dsbsgc41susy=0.d0
 
      do i=1,2
       dsbsgc41susy=dsbsgc41susy
     &        +(mass(kw)/mst2)**2
     &         *(chavmx(i,1)*sin(thetast)
     &           +chavmx(i,2)*cos(thetast)*mtw
     &            /(sqrt(2.d0)*sin(beta)*mass(kw)))**2
     &         *dsbsgechi(st2(i))
      enddo

      return
      end



