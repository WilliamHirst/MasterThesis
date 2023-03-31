      function dsbsgd2d()

***********************************************************************
* Function \Delta^(2)_d (with d=b) in app A p. 17 of                  *
* Ciuchini et al., hep-ph/9806308                                     *
* Note that we have inserted d=b                                      *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-07                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 thetast,x1,x2,m2sd2,ad
      real*8 dsbsgh1x,dsbsgh2xy,dsbsgh3x
      real*8 dsbsgd2d
      real*8 pi
      parameter (pi=3.141592653589793238d0)

c     Introduce the mixing angle in the stop sector
c     as it's defined p.7 in the ref stated above

      if(dimag(squlmx(6,3)).ne.0.d0) then
       write(*,*) 'squlmx(6,3) has a complex part
     &            so thetast in dsbsgd2d might be wrong'
      endif

      if(dimag(squrmx(6,3)).ne.0.d0) then
       write(*,*) 'squrmx(6,3) has a complex part
     &            so thetast in dsbsgd2d might be wrong'
      endif
   

c...Their convention is q1= cos(theta)q_L + sin(theta)q_R
c                       q2=-sin(tehta)q_L + cos(theta)q_R
c   with m_q1>m_q2. Hence, to be fully consistent with their
c   notation, we pick the corresponding eigenvalues and angles
c   from DarkSUSY

      if (mass(kst1).ge.mass(kst2)) then   ! DarkSUSY convention matches
        thetast=atan(-dble(conjg(squlmx(6,3)))
     &             /dble(conjg(squrmx(6,3))))

      else        ! DarkSUSY convention differs. Switch 1 <-> 2
        thetast=atan(-dble(conjg(squlmx(6,3)))
     &             /dble(conjg(squrmx(6,3))))+pi/2.0d0
      endif

c     Mass-squared ratios for the functions H_i
c     as well as a mass sq. and a trilinear parameter
c     Note that they will all have to be changed if 
c     d should be identified with another family 
c     than the third one

      if (mass(ksb1).ge.mass(ksb2)) then  ! DS convention agrees
        x1=mass(ksb1)**2/mass(kgluin)**2
        x2=mass(ksb2)**2/mass(kgluin)**2    
        m2sd2=mass(ksb2)**2
      else ! conventions differ, switch 1 <-> 2
        x1=mass(ksb2)**2/mass(kgluin)**2
        x2=mass(ksb1)**2/mass(kgluin)**2    
        m2sd2=mass(ksb1)**2
      endif

      ad=asoftd(3)

      dsbsgd2d=5.d0/2.d0
     &    +(1.d0-mass(kt)*mass(kgluin)/(tan(thetast)*m2sd2))
     &     *dsbsgh3x(x2)
     &    +dsbsgh1x(x1)/2.d0+dsbsgh1x(x2)/2.d0
     &    -2.d0*(ad-mu*tanbe)*dsbsgh2xy(x1,x2)/mass(kgluin)


      return
      end


