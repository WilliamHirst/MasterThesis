      function dsbsgeb()

***********************************************************************
* Function \epsilon_b in (10) of Degrassi et al.,                     *
* hep-ph/0009337                                                      *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-03                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      integer i,j
      real*8 alpha3mu,mt
      real*8 yt,a,b,x(2),y(2)
      real*8 dsbsgh2xy
      real*8 dsbsgeb,dsbsgalpha3
      real*8 pi
      parameter (pi=3.141592653589793238d0)

c     We here set mt=m_t(mt), i.e. at the scale mt which should be the
c     correct scale for terms in Kt

      mt=mtmt  ! Use DarkSUSY value at scale mt
 
c     value of alpha3mu=alpha3(\mu_0)
c     Since this term enters in Kt we should use the scale mt

      alpha3mu=dsbsgalpha3(mt)
 
c     Then set the mass ratio y_t of eq. (4)
      
      yt=mt**2/mass(khc)**2


c     mass-squared ratios for the function H_2

      a=mass(ksb1)**2/mass(kgluin)**2
      b=mass(ksb2)**2/mass(kgluin)**2

      do i=1,2
       x(i)=mass(kst1)**2/mass(kcha(i))**2
       y(i)=mass(kst2)**2/mass(kcha(i))**2
      enddo

c     In the following we use asoftu, defined at
c     the weak scale


      dsbsgeb=-(2.d0/3.d0)*(alpha3mu/pi)
     &     *(mu/mass(kgluin))*dsbsgh2xy(a,b)
       
      do j=1,2
       dsbsgeb=dsbsgeb-(yt**2/(16.d0*pi**2))*
     &    chaumx(j,2)*(asoftu(3)/mass(kcha(j)))
     &    *dsbsgh2xy(x(j),y(j))*chavmx(j,2)    
      enddo

      return
      end


