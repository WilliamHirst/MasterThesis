      function dsbsgd2td(famd)

***********************************************************************
* Function \Delta^(2)_{t,d} in app A p. 17 of                         *
* Ciuchini et al., hep-ph/9806308                                     *
* The input parameter famd is the family of the down sector           *
* it should be 1 for the first family, 2 for the 2nd and 3 for the 3rd*
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-07                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      integer famd
      real*8 x1,au,u1,u2
      real*8 dsbsgh1x,dsbsgh2xy,dsbsgh3x
      real*8 dsbsgd2td

c     Mass-squared ratios for the function H_3
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
       write(*,*) 'dsbsgd2td called with wrong family nb.'  
       stop 
      endif 

c     I assumed that the tril. par. should be 
c     from the top sector
       au=asoftu(3)             

c     Mass-squared ratios for the function H_1 and H_2

      if (mass(kst2).le.mass(kst1)) then
        u1=mass(kst1)**2/mass(kgluin)**2
        u2=mass(kst2)**2/mass(kgluin)**2
      else ! switch 1 <-> 2 to match convention
        u1=mass(kst2)**2/mass(kgluin)**2
        u2=mass(kst1)**2/mass(kgluin)**2
      endif 

      dsbsgd2td=5.d0/2.d0
     &    +dsbsgh3x(x1)
     &    +dsbsgh1x(u1)/2.d0+dsbsgh1x(u2)/2.d0
     &    -2.d0*(au-mu/tanbe)*dsbsgh2xy(u1,u2)/mass(kgluin)


      return
      end


