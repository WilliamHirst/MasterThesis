      function dsbsgud()

***********************************************************************
* Function U_d (with d=b) in app A p. 15-6 of                         *
* Ciuchini et al., hep-ph/9806308                                     *
* Note that we have inserted d=b                                      *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-08                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 thetast,x1,x2,u1,u2,ad
      real*8 dsbsgh1x,dsbsgh2xy
      real*8 dsbsgud
      real*8 pi
      parameter (pi=3.141592653589793238d0)

c     Introduce the mixing angle in the stop sector
c     as it's defined p.7 in the ref stated above

c...Their convention is q1= cos(theta)q_L + sin(theta)q_R
c                       q2=-sin(tehta)q_L + cos(theta)q_R
c   with m_q1>m_q2. Hence, to be fully consistent with their
c   notation, we pick the corresponding eigenvalues and angles
c   from DarkSUSY

      if (mass(kst1).ge.mass(kst2)) then   ! DarkSUSY convention matches
        thetast=atan(-dble(conjg(squlmx(6,3)))
     &             /dble(conjg(squrmx(6,3))))
        u1=mass(kst1)**2/mass(kgluin)**2
        u2=mass(kst2)**2/mass(kgluin)**2        

      else        ! DarkSUSY convention differs. Switch 1 <-> 2
        thetast=atan(-dble(conjg(squlmx(6,3)))
     &             /dble(conjg(squrmx(6,3))))+pi/2.0d0
        u1=mass(kst2)**2/mass(kgluin)**2
        u2=mass(kst1)**2/mass(kgluin)**2        
      endif

c     Mass-squared ratios for the functions H_i
c     as well as a trilinear parameter
c     Note that they will all have to be changed if 
c     d should be identified with another family 
c     than the third one

      if (mass(ksb1).ge.mass(ksb2)) then ! DS conveniton agrees
        x1=mass(ksb1)**2/mass(kgluin)**2
        x2=mass(ksb2)**2/mass(kgluin)**2    
      else  ! conventions differ, switch 1 <-> 2
        x1=mass(ksb2)**2/mass(kgluin)**2
        x2=mass(ksb1)**2/mass(kgluin)**2    
      endif

      ad=asoftd(3)


      dsbsgud=0.5d0*(dsbsgh1x(x1)-(cos(thetast))**2*dsbsgh1x(u1)
     &         -(sin(thetast))**2*dsbsgh1x(u2))
     &    -2.d0*(ad-mu*tanbe)*dsbsgh2xy(x1,x2)/mass(kgluin)
     &    +2.d0*(ad-mu*tanbe)/mass(kgluin)
     &     *((cos(thetast))**2*dsbsgh2xy(u1,x2)
     &         +(sin(thetast))**2*dsbsgh2xy(u2,x2))

      return
      end


