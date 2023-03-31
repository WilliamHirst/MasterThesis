      function dsbsgc71chisusy()

***********************************************************************
* The next to leading order contribution to                           *
* the Wilson coefficient C_7 from chargino (susy)                     *
* Eq (13) of Ciuchini et al.,                                         *
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
      real*8 dsbsgg7chij1,dsbsgg7chi2,dsbsgri
      real*8 dsbsgd7chi1,dsbsgd7chi2
      real*8 dsbsgf71,dsbsgf73 
      real*8 dsbsgc71chisusy,mst2
      real*8 pi
      parameter (pi=3.141592653589793238d0)


      beta=atan(tanbe)                        


c     Introduce the mixing angle in the stop sector
c     as it's defined p.7 in the ref stated above

      if(dimag(squlmx(6,3)).ne.0.d0) then
       write(*,*) 'squlmx(6,3) has a complex part
     &            so thetast in c71chisusy might be wrong'
      endif

      if(dimag(squrmx(6,3)).ne.0.d0) then
       write(*,*) 'squrmx(6,3) has a complex part
     &            so thetast in c71chisusy might be wrong'
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

      mt=mtmt  ! mt at scale mt, value from DarkSUSY
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



      dsbsgc71chisusy=0.d0
 
      do i=1,2
       dsbsgc71chisusy=dsbsgc71chisusy
     &      +t2j(i)**2*mass(kw)**2/mst2**2
     &       *(dsbsgg7chij1(st2(i),i)
     &         +dsbsgd7chi1(st2(i))*log(muw**2/mass(kcha(i))**2)
     &         -(4.d0/9.d0)*dsbsgechi(st2(i)))
     &      +chaumx(i,2)*sin(thetast)*mass(kw)
     &       /(sqrt(2.d0)*cos(beta)*mass(kcha(i)))
     &       *(t2j(i)*dsbsgg7chi2(st2(i))
     &         +t2j(i)*dsbsgd7chi2(st2(i))
     &          *log(muw**2/mass(kcha(i))**2)
     &         -(4.d0/3.d0)*yj(i)*dsbsgri(2)
     &          *dsbsgf73(st2(i)))
     &      -(8.d0/9.d0)*yj(i)*(mass(kw)**2/mst2**2)
     &       *t2j(i)*(dsbsgri(2)+dsbsgri(3))*dsbsgf71(st2(i))
      enddo

      return
      end



