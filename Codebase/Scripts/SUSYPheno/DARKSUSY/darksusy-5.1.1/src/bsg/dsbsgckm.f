      function dsbsgckm()

***********************************************************************
* The ratio,|V_ts^* V_tb/V_cb|^2, of ckm elements                     *
* with susy corrections                                               *
* Eq (34) of Ciuchini et al., hep-ph/9806308                          * 
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-22                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      integer i,j
      real*8 beta,thetast,muw,mtw,st2(2),t2(2),x(2,2)
      real*8 mtpole,tw,th,hw,delsm,delh,deltil,r
      real*8 smratio,lambda
      real*8 dsbsgmtmuw,dsbsgfxy,dsbsggxy
      real*8 dsbsgckm
      real*8 pi
      parameter (pi=3.141592653589793238d0)

      beta=atan(tanbe)                        

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
c        st2(j) is called y_kj (with k=2) in the ref.
         st2(i)=mass(kst2)**2/mass(kcha(i))**2
        enddo
        thetast=atan(-dble(conjg(squlmx(6,3)))
     &             /dble(conjg(squrmx(6,3))))

      else        ! DarkSUSY convention differs. Switch 1 <-> 2
        do i=1,2
c        st2(j) is called y_kj (with k=2) in the ref.
         st2(i)=mass(kst1)**2/mass(kcha(i))**2
        enddo
        thetast=atan(-dble(conjg(squlmx(6,3)))
     &             /dble(conjg(squrmx(6,3))))+pi/2.0d0
      endif

c     Set the running top mass
c     uses the function dsbsgmtmuw with input m
c     Use scale mass(kw) here?

      muw=mass(kw)
      mtw=dsbsgmtmuw(muw)

c     set the top pole mass, if this is what should be 
c     used in \tilde\Delta in (B5)
c     We use the mass given by DarkSUSY

      mtpole=mass(kt)


c     mass-squared ratios for \Delta_SM and \Delta_H 
c     in (B1) and (B3)

      tw=mtw**2/mass(kw)**2
      th=mtw**2/mass(khc)**2
      hw=mass(khc)**2/mass(kw)**2

c     mass-squared ratios for \tilde\Delta in (B5)
c     they are defined in eq 6+7 and text below (B5)

      do i=1,2
       t2(i)=chavmx(i,1)*sin(thetast)
     &       +chavmx(i,2)*cos(thetast)*mtw
     &        /(sqrt(2.d0)*sin(beta)*mass(kw))
      enddo
    
      do i=1,2
       do j=1,2
        x(i,j)=mass(kcha(i))**2/mass(kcha(j))**2
       enddo
      enddo

c    calculate \Delta_SM (B1), \Delta_H (B3)
c    and \tilde\Delta (B5)

      delsm=(tw**3-12.d0*tw**2+15.d0*tw-4.d0
     &       +6.d0*tw**2*log(tw))/(4.d0*(tw-1.d0)**3)

      delh=(1.d0/tanbe)**4*(th**2-1.d0-2.d0*th*log(th))
     &     *th/(4*(th-1.d0)**3)
     &     +2.d0*tw*(dsbsgfxy(tw,hw)+dsbsggxy(tw,hw)/4.d0)
     &      /tanbe**2
      


      deltil=0.d0
 
      do i=1,2
       do j=1,2
        deltil=deltil
     &        +mass(kw)**4/(mtpole**2*mass(kcha(j))**2)
     &         *t2(i)**2*t2(j)**2*dsbsggxy(st2(j),x(i,j))
c         write(*,*) i,j,dsbsggxy(st2(j),x(i,j))
       enddo
      enddo

c     calculate the effect of supersymmetry as parametrized 
c     in eq (28). \Delta_H + \tilde\Delta = \Delta_susy (B2)

      r=1.d0+(delh+deltil)/delsm

c     set the Standard Model value of the relevant ratio
c     of ckm elements and set the value of \lambda=|V_us|
c     We use the values stated in hep-ph/0104034 sec 4
c     The values in DarkSUSY can be used instead if they are updated
      
      smratio=0.971d0 

      lambda=0.2237d0


c     We calculate the new ratio of the relevant ckm elements
c     as given in eq (34)

c      write(*,*) 'dsbsgckm: ',deltil,t2(1),t2(2),st2(1),st2(2),
c     &  x(1,1),x(1,2),x(2,1),x(2,2)

      dsbsgckm=smratio+(lambda**2+1.d0-smratio)
     &         *(1.d0-1.d0/sqrt(r))


      return
      end



