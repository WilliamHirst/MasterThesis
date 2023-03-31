
      function dsbsgkt(flag)
c      program dsackt
***********************************************************************
* Program that calculates K_t in (3.8) of Gambino and Misiak,         *
* Nucl. Phys. B611 (2001) 338                                         *
* for the calculation of b --> s gamma                                *
* Input: flag:  0 = only standard model                               *
*               1 = standard model plus SUSY corrections              *
* author:Mia Schelke, schelke@physto.se, 2003-03-24                   *
***********************************************************************
      implicit none
      include 'dsmssm.h'
      integer i,j,k,flag
      real*8 mu0,mt,mub,mb
      real*8 eta,x,alpha3mb,alpha3mu
      real*8 deriv1,deriv2,deriv,sum,temp
      real*8 ak(8),ek(8)
      real*8 dsbsgat0,dsbsgft0,dsbsget0,dsbsgat1,dsbsgft1
      real*8 dsbsgalpha3
      real*8 bsgktreal
      complex*16 dsbsgkt
      real*8 pi
      parameter (pi=3.141592653589793238d0)

c     We here set mu0=\mu_0=mt
c     \mu_0=mt is the central value of \mu_0, app.A

      mt=mtmt  ! Use DarkSUSY value at scale mt
      mu0=mt


c     We here set mub=\mu_b=mb and mb=m_b^1S=4.69
c     \mu_b=mb is the central value of \mu_b, app.A
c     For mb we use the value stated in app.A

      mb=4.69d0
      mub=mb


c     value for alpha3mb calculated at scale mb

      alpha3mb=dsbsgalpha3(mb)
      eta=dsbsgalpha3(mt)/alpha3mb

c     value of alpha3mu=alpha3(\mu_0)
      alpha3mu=dsbsgalpha3(mt)


c     x=(m_t(mu_0)/m_W)^2
      x=(mt/mass(kw))**2


c The so-called 'Magic numbers',table 1. 
      data (ak(i),i=1,8)/0.6086956521739131d+00,
     &  0.6956521739130435d+00,0.2608695652173913d+00,
     &  -0.5217391304347826d+00,
     &  0.4086d0,-0.4230d0,-0.8994d0,0.1456d0/

      data (ek(j),j=1,8)/5.2620d0,-3.8412d0,0.d0,0.d0,-1.9043d0,
     & -0.1008d0,0.1216d0,0.0183d0/


c     The following expression contains the derivatives of A^t_0 
c     and F^t_0 -- derived by hand.
c     It also contains the multiplicative factor ln(\mu_0/m_t)
c     which is zero for \mu_0=m_t
 
      deriv1=log(x)*((-9.d0*x**2+4.d0*x)/(2.d0*(x-1.d0)**4)
     &       -2.d0*(-3.d0*x**3+2.d0*x**2)/(x-1.d0)**5)
     &      +(-66.d0*x**2+306.d0*x-159.d0)/(36.d0*(x-1.d0)**3)
     &      +(22.d0*x**3-171.d0*x**2+171.d0*x-46.d0)/
     &       (12.d0*(x-1.d0)**4)
      deriv2=log(x)*(3.d0*x/(x-1.d0)**4-6.d0*x**2/(x-1.d0)**5)
     &      +(-15.d0*x**2+18.d0*x-30.d0)/(12.d0*(x-1.d0)**3)
     &      +(5.d0*x**3-9.d0*x**2+36.d0*x-8.d0)/(4.d0*(x-1.d0)**4)

      deriv=(alpha3mu/pi)*log(mu0/mt)*4.d0*x*
     & (-0.5d0*eta**(4.d0/23.d0)*deriv1
     &  +(4.d0/3.d0)*(eta**(4.d0/23.d0)-eta**(2.d0/23.d0))*deriv2)


      sum=0.d0      

      do k=1,8
        temp=ek(k)*eta**(ak(k)+11.d0/23.d0)
        sum=sum+temp
      enddo

c     In kt log terms and the deriv term vanish when we uses 
c     the central value of the \mu's  

      bsgktreal=deriv+
     &   (1.d0-(2.d0/9.d0)*alpha3mb**2)*
     &   (-0.5d0*eta**(4.d0/23.d0)*dsbsgat0(x,flag)
     &    +(4.d0/3.d0)*(eta**(4.d0/23.d0)-eta**(2.d0/23.d0))
     &     *dsbsgft0(x,flag))
     &  +alpha3mb/(4.d0*pi)*
     &   (dsbsget0(x,flag)*sum
     &    +eta**(4.d0/23.d0)*
     &     (-0.5d0*eta*dsbsgat1(x,flag)
     &      +dsbsgat0(x,flag)*(12523.d0/3174.d0-7411.d0*eta/4761.d0
     &       -2.d0*pi**2/9.d0
     &       -(4.d0/3.d0)*(log(mb/mub)+eta*log(mu0/mt)))
     &      +(4.d0/3.d0)*eta*dsbsgft1(x,flag)
     &      +dsbsgft0(x,flag)*(-50092.d0/4761.d0
     &       +1110842.d0*eta/357075.d0+16.d0*pi**2/27.d0
     &       +(32.d0/9.d0)*(log(mb/mub)+eta*log(mu0/mt))))
     &    +eta**(2.d0/23.d0)*
     &     (-(4.d0/3.d0)*eta*dsbsgft1(x,flag)
     &     +dsbsgft0(x,flag)*(2745458.d0/357075.d0
     &       -38890.d0*eta/14283.d0-(4.d0/9.d0)*pi**2
     &       -(16.d0/9.d0)*(log(mb/mub)+eta*log(mu0/mt)))))

      dsbsgkt=dcmplx(bsgktreal,
     &    -alpha3mb/(4.d0*pi)*eta**(2.d0/23.d0)*(4.d0/9.d0)*pi
     &     *dsbsgft0(x,flag))

c      write (*,*) 'K_t',dsbsgkt
c      write (*,*) 'alpha_s(mb)',alpha3mb
c      write (*,*) 'alpha_s(mt)',alpha3mu
c      write (*,*) 'eta',eta
c      write (*,*) 'K_t_r',bsgktreal
c      write (*,*) 'deriv',deriv
c      write (*,*) 'sum',sum
c      write (*,*) 'logt',log(mu0/mt)
c      write (*,*) 'mu0',mu0 
c      write (*,*) 'mt',mt 
c      write (*,*) 'logb',log(mub/mb) 
c      write (*,*) 'mub',mub 
c      write (*,*) 'mb',mb 
  
      return
      end



