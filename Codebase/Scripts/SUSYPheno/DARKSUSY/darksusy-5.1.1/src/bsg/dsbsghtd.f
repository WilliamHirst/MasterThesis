      function dsbsghtd(famd)

***********************************************************************
* Function H_t^d in app A p. 16 of                                    *
* Ciuchini et al., hep-ph/9806308                                     *
* The input parameter famd is the family of the down sector           *
* it should be 1 for the first family, 2 for the 2nd and 3 for the 3rd*
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-08                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      integer famd
      real*8 x1,au,u1,u2,thetast
      real*8 dsbsgh1x,dsbsgh2xy
      real*8 dsbsghtd
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


c     Mass-squared ratio for the functions H_i
c     as well as a trilinear parameter

      if(famd.eq.1) then

        if (mass(ksd1).ge.mass(ksd2)) then ! DarkSUSY convention matches
          x1=mass(ksd1)**2/mass(kgluin)**2
        else  ! DarkSUSY convention differs
          x1=mass(ksd2)**2/mass(kgluin)**2
        endif

      elseif(famd.eq.2) then

        if (mass(kss1).ge.mass(kss2)) then  ! conventions agree
          x1=mass(kss1)**2/mass(kgluin)**2
        else   ! conventions differ, switch 1 <-> 2
          x1=mass(kss2)**2/mass(kgluin)**2
        endif

      elseif(famd.eq.3) then

        if(mass(ksb1).ge.mass(ksb2)) then  ! conventions agree
          x1=mass(ksb1)**2/mass(kgluin)**2
        else  ! conventions differ, switch 1 <-> 2
          x1=mass(ksb2)**2/mass(kgluin)**2
        endif

      else
       write(*,*) 'dsbsghtd called with wrong family nb.'  
       stop 
      endif 

c     I assumed that the tril. par. should be 
c     from the top sector
       au=asoftu(3)             

      dsbsghtd=-0.5d0*(dsbsgh1x(x1)-(cos(thetast))**2*dsbsgh1x(u1)
     &         -(sin(thetast))**2*dsbsgh1x(u2))
     &    -2.d0*(au-mu/tanbe)*dsbsgh2xy(u1,u2)/mass(kgluin)
     &    +2.d0*(au+mu*tanbe)/mass(kgluin)
     &     *((cos(thetast))**2*dsbsgh2xy(u2,x1)
     &         +(sin(thetast))**2*dsbsgh2xy(u1,x1))   


      return
      end


