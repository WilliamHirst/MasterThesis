      function dsbsgkc()
***********************************************************************
* Program that calculates K_c in (3.7) of Gambino and Misiak,         *
* Nucl. Phys. B611 (2001) 338                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-03-25                   *
***********************************************************************
      implicit none
      include 'dsmssm.h'
      integer i,j
      real*8 mu0,mub,mb
      real*8 eta,alpha3mb
      real*8 ak(8),dk(8),dtk(8),dtek(8),dtak(8),dtbk(8),dtipik(8)
      real*8 sumim,tempim,sumreal,tempreal
      real*8 bsgkcim,bsgkcreal,dsbsgalpha3
      complex*16 dsbsgkc
      real*8 pi
      parameter (pi=3.141592653589793238d0)

c     We here set mu0=\mu_0=mass(kw)
c     \mu_0=mw is the central value of \mu_0, app.A
c     For mw we use the ds value instead of the value stated 
c     in app.A

      mu0=mass(kw)

c     We here set mub=\mu_b=mb and mb=m_b^1S=4.69
c     \mu_b=mb is the central value of \mu_b, app.A
c     For mb we use the value stated in app.A
c     Note: We chose to use this hardcoded value instead of the pole
c     mass obtained in DarkSUSY (4.57 GeV) for consistency

      mb=4.69d0
      mub=mb
 
c     value for alpha3mb calculated at scale mb

      alpha3mb=dsbsgalpha3(mb)
      eta=dsbsgalpha3(mass(kw))/alpha3mb

c The so-called 'Magic numbers',table 1. 
      data (ak(i),i=1,8)/0.6086956521739131d+00,
     &  0.6956521739130435d+00,0.2608695652173913d+00,
     &  -0.5217391304347826d+00,
     &  0.4086d0,-0.4230d0,-0.8994d0,0.1456d0/


      data (dk(i),i=1,8)/1.4107d0,-0.8380d0,-0.4286d0,
     & -0.0714d0,-0.6494d0,-0.0380d0,-0.0185d0,-0.0057d0/
 
c magic numbers from Nucl. Phys. B611 (2001) 338,
c for higher order see below
c      data (dtk(i),i=1,8)/-17.6507d0,11.3460d0,3.5762d0,
c     & -2.2672d0,3.9267d0,1.1366d0,-0.5445d0,0.1653d0/

c magic numbers from Buras et al. hep-ph/0203135
      data (dtk(i),i=1,8)/-17.6507d0,11.3460d0,2.4692d0,
     & -0.8056d0,4.8898d0,-0.2308d0,-0.5290d0,0.1994d0/


      data (dtek(i),i=1,8)/9.2746d0,-6.9366d0,-0.8740d0,
     & 0.4218d0,-2.7231d0,0.4083d0,0.1465d0,0.0205d0/

c magic numbers from Nucl. Phys. B611 (2001) 338,
c for higher order see below
c      data (dtak(i),i=1,8)/0.d0,0.d0,1.d0,
c     & 1.d0,0.d0,0.d0,0.d0,0.d0/

c magic numbers from Buras et al. hep-ph/0203135
      data (dtak(i),i=1,8)/0.d0,0.d0,0.8571d0,
     & 0.6667d0,0.1298d0,0.1951d0,0.1236d0,0.0276d0/

c magic numbers from Nucl. Phys. B611 (2001) 338,
c for higher order see below
c      data (dtbk(i),i=1,8)/0.d0,0.d0,1.d0,
c     & 1.d0,0.d0,0.d0,0.d0,0.d0/

c magic numbers from Buras et al. hep-ph/0203135
      data (dtbk(i),i=1,8)/0.d0,0.d0,0.8571d0,
     & 0.6667d0,0.2637d0,0.2906d0,-0.0611d0,-0.0171d0/


c magic numbers from Nucl. Phys. B611 (2001) 338,
c for higher order see below
c      data (dtipik(i),i=1,8)/0.4702d0,0.d0,-0.4938d0,
c     & -0.4938d0,-0.8120d0,0.0776d0,-0.0507d0,0.0186d0/

c magic numbers from Buras et al. hep-ph/0203135
      data (dtipik(i),i=1,8)/0.4702d0,0.d0,-0.4268d0,
     & -0.2222d0,-0.9042d0,-0.1150d0,-0.0975d0,0.0115d0/



      sumim=0.d0

      do j=1,8
        tempim=eta**ak(j)*alpha3mb/(4.d0*pi)*
     &         (dtak(j)*1.01d0+dtbk(j)*0.09d0
     &          +dtipik(j)*pi)
        sumim=sumim+tempim
      enddo

c the scale of alpha3 that enters here is not specified in eq. (3.7),
c but it seems reasonable that it shoulbe be mb, hence we use
c slpha3mb below

      bsgkcim=sumim
     & +dble(conjg(ckm(1,2))*ckm(1,3)/(conjg(ckm(3,2))*ckm(3,3)))
     &  *alpha3mb/(4.d0*pi)*(eta**ak(3)+eta**ak(4))
     &  *(1.01d0+0.09d0)
     & +dimag(conjg(ckm(1,2))*ckm(1,3)/(conjg(ckm(3,2))*ckm(3,3)))
     &  *alpha3mb/(4.d0*pi)*(eta**ak(3)+eta**ak(4))
     &  *(0.97d0-0.04d0)


c The log terms below vanish at the central values of the \mu's
      sumreal=0.d0

      do j=1,8
        tempreal=eta**ak(j)*dk(j)
     &    +eta**ak(j)*alpha3mb/(4.d0*pi)*
     &     (2.d0*(23.d0/3.d0)*ak(j)*dk(j)*
     &       (log(mb/mub)+eta*log(mu0/mass(kw)))
     &      +dtk(j)+dtek(j)*eta
     &      +dtak(j)*0.97d0+dtbk(j)*(-0.04d0))
        sumreal=sumreal+tempreal
      enddo

c the scale of alpha3 that enters here is not specified in eq. (3.7),
c but it seems reasonable that it shoulbe be mb, hence we use
c slpha3mb below

      bsgkcreal=sumreal
     & +dble(conjg(ckm(1,2))*ckm(1,3)/(conjg(ckm(3,2))*ckm(3,3)))
     &  *alpha3mb/(4.d0*pi)*(eta**ak(3)+eta**ak(4))
     &  *(0.97d0-0.04d0)
     & -dimag(conjg(ckm(1,2))*ckm(1,3)/(conjg(ckm(3,2))*ckm(3,3)))
     &  *alpha3mb/(4.d0*pi)*(eta**ak(3)+eta**ak(4))
     &  *(1.01d0+0.09d0)

      dsbsgkc=dcmplx(bsgkcreal,bsgkcim)

c      write (*,*) 'K_c',dsbsgkc
c      write (*,*) 'Im(K_c)',bsgkcim
c      write (*,*) 'sumim',sumim  
c      write (*,*) 'Im(K_c)',dimag(dsbsgkc)
c      write (*,*) 'Re(K_c)',dble(dsbsgkc)
c      write (*,*) 'sumreal',sumreal
c      write (*,*) 'ckm(3,2)',ckm(3,2)     
c      write (*,*) 'ckm(3,3)',ckm(3,3)     
c      write (*,*) 'ckm(1,2)',ckm(1,2)     
c      write (*,*) 'ckm(1,3)',ckm(1,3)              

      return
      end



