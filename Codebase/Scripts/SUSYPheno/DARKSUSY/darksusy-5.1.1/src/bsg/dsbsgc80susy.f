      function dsbsgc80susy()

***********************************************************************
* The leading order contribution to the Wilson coefficient C_8        *
* from susy                                                           *
* Eq (31) of Degrassi et al., hep-ph/0009337                          *
* Differs from dsbsgc70susy only by a few changes described p.11      *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-04                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      integer i
      real*8 beta
      real*8 mcomsq,st1(2),st2(2),sq(2)
      real*8 thetast,mt
      real*8 mus,muw,eta,mtw,mts,keb,beta0
      real*8 dsbsgalpha3,dsbsgeb
      real*8 dsbsgf81,dsbsgf83
      real*8 dsbsgmtmuw   ! ,dsbsgyt
      real*8 dsbsgc80susy,mst1,mst2
      real*8 pi
      parameter (pi=3.141592653589793238d0)


      mt=mtmt  ! mt at scale mt, value from DarkSUSY
      beta=atan(tanbe)

c     set a common mass for the two first generations of squarks
c     here set to \tilde{u}_1

      mcomsq=mass(ksu1)


c     Introduce the mixing angle in the stop sector
c     as it's defined p.4 in the ref stated above

c...Their convention is q1= cos(theta)q_L + sin(theta)q_R
c                       q2=-sin(tehta)q_L + cos(theta)q_R
c   with m_q1>m_q2. Hence, to be fully consistent with their
c   notation, we pick the corresponding eigenvalues and angles
c   from DarkSUSY


      if (mass(kst1).ge.mass(kst2)) then   ! DarkSUSY convention matches
        do i=1,2
         st1(i)=mass(kst1)**2/mass(kcha(i))**2
         st2(i)=mass(kst2)**2/mass(kcha(i))**2
         sq(i)=mcomsq**2/mass(kcha(i))**2
        enddo
        mst1=mass(kst1)
        mst2=mass(kst2)
        thetast=atan(-dble(conjg(squlmx(6,3)))
     &             /dble(conjg(squrmx(6,3))))

      else        ! DarkSUSY convention differs. Switch 1 <-> 2
        do i=1,2
         st1(i)=mass(kst2)**2/mass(kcha(i))**2
         st2(i)=mass(kst1)**2/mass(kcha(i))**2
         sq(i)=mcomsq**2/mass(kcha(i))**2
        enddo
        mst1=mass(kst2)
        mst2=mass(kst1)
        thetast=atan(-dble(conjg(squlmx(6,3)))
     &             /dble(conjg(squrmx(6,3))))+pi/2.0d0
      endif

c     Calculate eta=alpha3(mu_susy)/alpha3(mu_w), p.8
c     uses the function dsbsgalpha3 with input m
c     p.10 says that we should take the gluino mass scale
c     for mu_susy and for mu_w take the chargino/light stop scale
c     but would we rather take mt since these terms enter in Kt
c     The calculation of eta could be checked using the alternative 
c     expression stated p.8:
c     eta=1.d0/(1.d0-(-7.d0)/(2.d0*pi)*log(mus/muw))

      mus=mass(kgluin)
      muw=mt

      eta=dsbsgalpha3(mus)/dsbsgalpha3(muw)     ! MS tmp test below
c      eta=1.d0/(1.d0+(41.d0/(6.d0*2.d0*pi))*log(mus/muw))

c     Set the running top mass

      mtw=dsbsgmtmuw(muw)

      mts=dsbsgmtmuw(mus)

c     In the following expression for C_8 we have inserted 
c     the following definitions:

      beta0=-7.d0
c     this is \beta_0 (sometimes used with opposite sign convention)
c     for 6 active quarks (p.8) (2*n_f/3-11)
c     the code could also be tested with
c     beta0=-41.d0/6.d0
c     p.11+9 corresponding to 6 q flavours and one light squark
c     as is actually correct for their scenario

c     K=1/(1+e_b*tanbe) p.10, with e_b defined in eq. (10)
c     and in the function dsbsgeb

      keb=1.d0/(1.d0+dsbsgeb()*tanbe)

      dsbsgc80susy=0.d0
 
      do i=1,2
       dsbsgc80susy=dsbsgc80susy
     &  +eta**(-14.d0/(3.d0*beta0))
     &  *((2.d0/3.d0)*(mass(kw)**2/mcomsq**2)*chavmx(i,1)**2
     &    *dsbsgf81(sq(i))
     &    -(2.d0/3.d0)*(cos(thetast)*chavmx(i,1)
     &      -sin(thetast)*chavmx(i,2)*mts
     &       /(sqrt(2.d0)*sin(beta)*mass(kw)))**2
     &     *(mass(kw)/mst1)**2*dsbsgf81(st1(i))
     &    +(keb/cos(beta))
     &     *(chaumx(i,2)*chavmx(i,1)*mass(kw)
     &       /(sqrt(2.d0)*mass(kcha(i)))
     &       *(dsbsgf83(sq(i))-cos(thetast)**2*dsbsgf83(st1(i)))
     &       +sin(thetast)*cos(thetast)*chaumx(i,2)*chavmx(i,2)
     &        *mts*dsbsgf83(st1(i))
     &        /(2.d0*sin(beta)*mass(kcha(i)))))
     & -(2.d0/3.d0)*eta**(4.d0/beta0)
     &  *(sin(thetast)*chavmx(i,1)+cos(thetast)*chavmx(i,2)*mts
     &      /(sqrt(2.d0)*sin(beta)*mass(kw)))**2
     &  *(mass(kw)/mst2)**2*dsbsgf81(st2(i))
     & -(1.d0/cos(beta))
     &  *(chaumx(i,2)*chavmx(i,1)*mass(kw)*sin(thetast)**2
     &    *(keb+dsbsgeb()*tanbe)/(sqrt(2.d0)*mass(kcha(i)))
     &    +sin(thetast)*cos(thetast)*chaumx(i,2)*chavmx(i,2)*mts
     &     /(2.d0*sin(beta)*mass(kcha(i)))
     &     *(keb+dsbsgeb()*tanbe*mtw/mts))*dsbsgf83(st2(i))
      enddo

      return
      end


