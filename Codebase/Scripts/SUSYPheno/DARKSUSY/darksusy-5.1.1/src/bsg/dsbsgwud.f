      function dsbsgwud(famu,famd)

***********************************************************************
* Function W_u^d in app A p. 15 of                                    *
* Ciuchini et al., hep-ph/9806308                                     *
* The input parameters famu and famd are the families of the up- and  *
* down sector respectively, and they should be 1 for the first family,*
* 2 for the 2nd and 3 for the 3rd                                     *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-08                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      integer famu,famd
      real*8 x1,w(2),thetasu
c      real*8 angle1(3),angle2(3)
      real*8 dsbsgwxy
      real*8 dsbsgwud
      real*8 pi
      parameter (pi=3.141592653589793238d0)

c     Mass-squared ratios for the function W[x,y]
c     as well as the mixing angle for the up squark sector
c     (as defined at p.7 in the reference)
c     We also include a test of the mixing angles
c      angle1(1)=atan(-dble(squlmx(4,1))/dble(squrmx(4,1))) 
c      angle2(1)=atan(dble(squrmx(1,1))/dble(squlmx(1,1))) 
c      angle1(2)=atan(-dble(squlmx(5,2))/dble(squrmx(5,2))) 
c      angle2(2)=atan(dble(squrmx(2,2))/dble(squlmx(2,2)))
c      angle1(3)=atan(-dble(squlmx(6,3))/dble(squrmx(6,3))) 
c      angle2(3)=atan(dble(squrmx(3,3))/dble(squlmx(3,3))) 


      if(dimag(squrmx(4,1)).ne.0.d0) then
       write(*,*) 'squrmx(4,1) has a complex part
     &            so thetast in dsbsgwud might be wrong'
      endif
      if(dimag(squlmx(4,1)).ne.0.d0) then
       write(*,*) 'squlmx(4,1) has a complex part
     &            so thetast in dsbsgwud might be wrong'
      endif
      if(dimag(squrmx(5,2)).ne.0.d0) then
       write(*,*) 'squrmx(5,2) has a complex part
     &            so thetast in dsbsgwud might be wrong'
      endif
      if(dimag(squlmx(5,2)).ne.0.d0) then
       write(*,*) 'squlmx(6,3) has a complex part
     &            so thetast in dsbsgwud might be wrong'
      endif
      if(dimag(squrmx(6,3)).ne.0.d0) then
       write(*,*) 'squrmx(6,3) has a complex part
     &            so thetast in dsbsgwud might be wrong'
      endif
      if(dimag(squlmx(6,3)).ne.0.d0) then
       write(*,*) 'squlmx(6,3) has a complex part
     &            so thetast in dsbsgwud might be wrong'
      endif


c...Their convention is q1= cos(theta)q_L + sin(theta)q_R
c                       q2=-sin(tehta)q_L + cos(theta)q_R
c   with m_q1>m_q2. Hence, to be fully consistent with their
c   notation, we pick the corresponding eigenvalues and angles
c   from DarkSUSY

      if(famu.eq.1) then

        if (mass(ksu1).ge.mass(ksu2)) then  ! DarkSUSY convention matches
          w(1)=mass(ksu1)**2/mass(kgluin)**2
          w(2)=mass(ksu2)**2/mass(kgluin)**2
          thetasu=atan(-dble(conjg(squlmx(4,1)))
     &              /dble(conjg(squrmx(4,1)))) 
        else    ! DarkSUSY convention differs. Switch 1 <-> 2
          w(1)=mass(ksu2)**2/mass(kgluin)**2
          w(2)=mass(ksu1)**2/mass(kgluin)**2
          thetasu=atan(-dble(conjg(squlmx(4,1)))
     &              /dble(conjg(squrmx(4,1)))) +pi/2.0d0
        endif

      elseif(famu.eq.2) then
        if (mass(ksc1).ge.mass(ksc2)) then  ! DarkSUSY convention matches
          w(1)=mass(ksc1)**2/mass(kgluin)**2
          w(2)=mass(ksc2)**2/mass(kgluin)**2
         thetasu=atan(-dble(conjg(squlmx(5,2)))
     &              /dble(conjg(squrmx(5,2)))) 
        else! DarkSUSY convention differs. Switch 1 <-> 2
          w(1)=mass(ksc2)**2/mass(kgluin)**2
          w(2)=mass(ksc1)**2/mass(kgluin)**2
         thetasu=atan(-dble(conjg(squlmx(5,2)))
     &              /dble(conjg(squrmx(5,2)))) +pi/2.0d0
        endif

      elseif(famu.eq.3) then

        if (mass(kst1).ge.mass(kst2)) then   ! DarkSUSY convention matches
          thetasu=atan(-dble(conjg(squlmx(6,3)))
     &             /dble(conjg(squrmx(6,3))))
          w(1)=mass(kst1)**2/mass(kgluin)**2
          w(2)=mass(kst2)**2/mass(kgluin)**2

        else        ! DarkSUSY convention differs. Switch 1 <-> 2
          thetasu=atan(-dble(conjg(squlmx(6,3)))
     &             /dble(conjg(squrmx(6,3))))+pi/2.0d0
          w(1)=mass(kst2)**2/mass(kgluin)**2
          w(2)=mass(kst1)**2/mass(kgluin)**2
        endif

      else
       write(*,*) 'dsbsgwud called with wrong family nb.'  
       stop 
      endif 

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
       write(*,*) 'dsbsgwud called with wrong family nb.'  
       stop 
      endif 



      dsbsgwud=0.5d0*((cos(thetasu))**2*dsbsgwxy(x1,w(1))
     &                +(sin(thetasu))**2*dsbsgwxy(x1,w(2)))


      return
      end


