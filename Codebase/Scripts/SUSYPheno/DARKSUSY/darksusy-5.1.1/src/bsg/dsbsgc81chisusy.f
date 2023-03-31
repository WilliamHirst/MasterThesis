      function dsbsgc81chisusy()

***********************************************************************
* The next to leading order contribution to                           *
* the Wilson coefficient C_8 from chargino (susy)                     *
* Eq (14) of Ciuchini et al.,                                         *
* hep-ph/9806308                                                      *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-08                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      integer i
      real*8 thetast,st2(2),mt
      real*8 beta,muw,mtw
      real*8 yj(2),t2j(2)
      real*8 dsbsgmtmuw,dsbsgechi
      real*8 dsbsgg8chij1,dsbsgg8chi2,dsbsgri
      real*8 dsbsgd8chi1,dsbsgd8chi2
      real*8 dsbsgf81,dsbsgf83 
      real*8 dsbsgc81chisusy,mst2
      real*8 pi
      parameter (pi=3.141592653589793238d0)


      beta=atan(tanbe)                        


c     Introduce the mixing angle in the stop sector
c     as it's defined p.7 in the ref stated above

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

      mt=mtmt  ! mt at scale mt, use value from DarkSUSY
      muw=mt

      mtw=dsbsgmtmuw(muw)

c     Define some expressions used below
c     from eq (6) and (7)

      do i=1,2
       yj(i)=chavmx(i,2)*cos(thetast)*mtw
     &            /(sqrt(2.d0)*sin(beta)*mass(kw))
      enddo

      do i=1,2
       t2j(i)=chavmx(i,1)*sin(thetast)+yj(i)
      enddo



      dsbsgc81chisusy=0.d0
 
      do i=1,2
       dsbsgc81chisusy=dsbsgc81chisusy
     &      +t2j(i)**2*mass(kw)**2/mst2**2
     &       *(dsbsgg8chij1(st2(i),i)
     &         +dsbsgd8chi1(st2(i))*log(muw**2/mass(kcha(i))**2)
     &         -(1.d0/6.d0)*dsbsgechi(st2(i)))
     &      +chaumx(i,2)*sin(thetast)*mass(kw)
     &       /(sqrt(2.d0)*cos(beta)*mass(kcha(i)))
     &       *(t2j(i)*dsbsgg8chi2(st2(i))
     &         +t2j(i)*dsbsgd8chi2(st2(i))
     &          *log(muw**2/mass(kcha(i))**2)
     &         -(4.d0/3.d0)*yj(i)*dsbsgri(2)
     &          *dsbsgf83(st2(i)))
     &      -(8.d0/9.d0)*yj(i)*(mass(kw)**2/mst2**2)
     &       *t2j(i)*(dsbsgri(2)+dsbsgri(3))*dsbsgf81(st2(i))
      enddo

      return
      end



