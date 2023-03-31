      function dsbsgbofe(flag)
c      program dsacbofe     
***********************************************************************
* Program that calculates B(E_0) in (E.1) of Gambino and Misiak,      *
* Nucl. Phys. B611 (2001) 338                                         *
* for the calculation of b --> s gamma                                *
* Input: flag:  0 = only standard model                               *
*               1 = standard model plus SUSY corrections              *
* author:Mia Schelke, schelke@physto.se, 2003-03-12                   *
***********************************************************************
      implicit none
      include 'dsmssm.h'
      integer i,n,m,k,l,flag
      real*8 dsbsgbofe,mt
      real*8 z,delta,limit,bs,eta,x,alpha3mb
      real*8 sum,temp
      real*8 phi(8,8),ak(8),c0effmb(8)
      real*8 dsphi22a,dsphi22b,dsphi27a,dsphi27b
      real*8 int22a,int22b,int27a,int27b
      real*8 dsbsgat0,dsbsgft0
      real*8 dsf_int,ddilog,dsbsgalpha3,mb
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      external dsphi22a,dsphi22b,dsphi27a,dsphi27b

c     z=m_c^2/m_b^2 number from app.E hep-ph/0104034
      z=0.22d0**2
c     bs=m_b/m_s   number from app.E hep-ph/0104034
      bs=50.d0

      delta=0.9d0
      
      limit=(1.d0-delta)/z


c     For mb we use the value stated in app.A
c     Note: We chose to use this hardcoded value instead of the pole
c     mass obtained in DarkSUSY (4.57 GeV) for consistency

      mb=4.69d0

c     eta=alph3(mw)/alph3(mb)
c     value for alpha3mb calculated at scale mb
      alpha3mb=dsbsgalpha3(mb)
      eta=dsbsgalpha3(mass(kw))/alpha3mb


      mt=mtmt
c     x=(m_t(m_t)/m_W)^2
      x=(mt/mass(kw))**2         


c The so-called 'Magic numbers',table 1. 
      data (ak(i),i=1,8)/0.6086956521739131d+00,
     &  0.6956521739130435d+00,0.2608695652173913d+00,
     &  -0.5217391304347826d+00,
     &  0.4086d0,-0.4230d0,-0.8994d0,0.1456d0/

      int22a=dsf_int(dsphi22a,1.d-10,limit,1.d-4)
      int22b=dsf_int(dsphi22b,limit,1.d0/z,1.d-4)
      int27a=dsf_int(dsphi27a,0.d0,limit,1.d-4)
      int27b=dsf_int(dsphi27b,limit,1.d0/z,1.d-4)

      do n=1,8
       do m=1,8
        phi(n,m)=0.d0
       enddo
      enddo

      phi(2,2)=(16.d0*z/27.d0)*(delta*int22a+int22b)
      phi(2,7)=-(8.d0*z**2/9.d0)*(delta*int27a+int27b)
      phi(7,7)=-(2.d0/3.d0)*(log(delta))**2
     &         -(7.d0/3.d0)*log(delta)-31.d0/9.d0+delta*10.d0/3.d0
     &         +delta**2/3.d0-(2.d0/9.d0)*delta**3
     &         +(delta/3.d0)*(delta-4.d0)*log(delta)
      phi(7,8)=(8.d0/9.d0)*(ddilog(1.d0-delta)-pi**2/6.d0
     &         -delta*log(delta)+9.d0*delta/4.d0
     &         -delta**2/4.d0+delta**3/12.d0)
      
      phi(8,8)=(1.d0/27.d0)*(-2.d0*log(bs)*
     &         (delta**2+2.d0*delta+4.d0*log(1.d0-delta))
     &         +4.d0*ddilog(1.d0-delta)-2.d0*pi**2/3.d0
     &         -delta*(2.d0+delta)*log(delta)+8.d0*log(1.d0-delta)
     &         -2.d0*delta**3/3.d0+3.d0*delta**2+7.d0*delta)

      phi(1,1)=phi(2,2)/36.d0
      phi(1,2)=-phi(2,2)/3.d0
      phi(1,7)=-phi(2,7)/6.d0
      phi(1,8)=phi(2,7)/18.d0
      phi(2,8)=-phi(2,7)/3.d0

      c0effmb(1)=eta**ak(3)-eta**ak(4)
      c0effmb(2)=(2.d0/3.d0)*eta**ak(3)+(1.d0/3.d0)*eta**ak(4)
      c0effmb(3)=(2.d0/63.d0)*eta**ak(3)-(1.d0/27.d0)*eta*ak(4)
     &           -0.0659d0*eta**ak(5)+0.595d0*eta**ak(6)
     &           -0.0218d0*eta**ak(7)+0.0335d0*eta**ak(8)      
      c0effmb(4)=(1.d0/21.d0)*eta**ak(3)+(1.d0/9.d0)*eta*ak(4)
     &           +0.0237d0*eta**ak(5)-0.0173d0*eta**ak(6)
     &           -0.1336d0*eta**ak(7)-0.0316d0*eta**ak(8)     
      c0effmb(5)=(-1.d0/126.d0)*eta**ak(3)+(1.d0/108.d0)*eta*ak(4)
     &           +0.0094d0*eta**ak(5)-0.01d0*eta**ak(6)
     &           -0.001d0*eta**ak(7)-0.0017d0*eta**ak(8)    
      c0effmb(6)=(-1.d0/84.d0)*eta**ak(3)+(1.d0/36.d0)*eta*ak(4)
     &           +0.0108d0*eta**ak(5)+0.0163d0*eta**ak(6)
     &           +0.0103d0*eta**ak(7)+0.0023d0*eta**ak(8)
      c0effmb(7)=(-4.d0*dsbsgft0(x,flag)/3.d0+42678d0/30253d0)
     &           *eta**ak(1) 
     &           +(-dsbsgat0(x,flag)/2.d0+4.d0*dsbsgft0(x,flag)/3.d0
     &             -86697d0/103460d0)*eta**ak(2)
     &           -3.d0*eta**ak(3)/7.d0-eta**ak(4)/14.d0
     &           -0.6494*eta**ak(5)-0.0380*eta**ak(6)   
     &           -0.0185*eta**ak(7)-0.0057*eta**ak(8)        
      c0effmb(8)=(-dsbsgft0(x,flag)/2.d0+64017d0/121012d0)*eta**ak(1) 
     &           -0.9135*eta**ak(5)+0.0873*eta**ak(6)   
     &           -0.0571*eta**ak(7)+0.0209*eta**ak(8)        

      sum =0.d0      

      do k=1,8
       do l=k,8
        temp=c0effmb(k)*c0effmb(l)*phi(k,l)
        sum=sum+temp
       enddo
      enddo

      dsbsgbofe=alpha3mb*sum/pi
    
c      write (*,*) 'B(E_0)',dsbsgbofe
   
      return
      end



