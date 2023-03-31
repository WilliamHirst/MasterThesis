      real*8 function dsntlitlf_s(mx,vbar)
c----------------------------------------------------------------------
c       dsntlitlf_s used by capsun
c       mx: neutralino mass
c       lars bergstrom 1995-12-12
c----------------------------------------------------------------------
      implicit none
      include 'dsmpconst.h'
      real*8 mx,vbar,aa(10),dsntss
      real*8 m(10),f(10),phi(10),mn,dsff(10)
      integer i
      mn=(m_p+m_n)/2.0d0
      m(1)=1.d0*mn
      m(2)=4.d0*mn
      m(3)=12.d0*mn
      m(4)=14.d0*mn
      m(5)=16.d0*mn
      m(6)=20.d0*mn
      m(7)=24.d0*mn
      m(8)=28.d0*mn
      m(9)=32.d0*mn
      m(10)=56.d0*mn

c...Old jkg values
      f(1)=0.772d0
      f(2)=0.209d0
      f(3)=3.87d-3
      f(4)=9.4d-4
      f(5)=8.55d-3
      f(6)=1.51d-3
      f(7)=7.39d-4
      f(8)=8.13d-4
      f(9)=4.65d-4
      f(10)=1.46d-3

c...New BP2000 values
      f(1)=0.670d0
      f(2)=0.311d0
      f(3)=2.37d-3
      f(4)=1.88d-3
      f(5)=8.78d-3
      f(6)=1.93d-3
      f(7)=7.33d-4
      f(8)=7.98d-4
      f(9)=5.50d-4
      f(10)=1.42d-3


c..Old jkg values
      phi(1)=3.16d0
      phi(2)=3.4d0
      phi(3)=3.23d0
      do i=4,10
        phi(i)=3.23d0
      enddo

c...New averages based on integration of BP2000 solar model
      phi(1)=3.15d0
      phi(2)=3.40d0
      phi(3)=2.85d0
      phi(4)=3.83d0
      phi(5)=3.25d0
      do i=6,10
        phi(i)=3.22d0
      enddo

      dsff(1)=1.d0
      dsff(2)=0.986d0
      dsff(2)=dsff(2)+(1.d0-dsff(2))*
     &  dexp(-(dlog(mx)/dlog(18.2d0))**1.58d0)
      dsff(3)=0.788d0
      dsff(3)=dsff(3)+(1.d0-dsff(3))*
     &  dexp(-(dlog(mx)/dlog(61.6d0))**2.69d0)
      dsff(4)=0.613d0
      dsff(4)=dsff(4)+(1.d0-dsff(4))*
     &  dexp(-(dlog(mx)/dlog(75.2d0))**2.69d0)
      dsff(5)=0.613d0
      dsff(5)=dsff(5)+(1.d0-dsff(5))*
     &  dexp(-(dlog(mx)/dlog(75.2d0))**2.69d0)
      dsff(6)=0.613d0
      dsff(6)=dsff(6)+(1.d0-dsff(6))*
     &  dexp(-(dlog(mx)/dlog(75.2d0))**2.69d0)
      dsff(7)=0.281d0
      dsff(7)=dsff(7)+(1.d0-dsff(7))*
     &  dexp(-(dlog(mx)/dlog(71.7d0))**2.97d0)
      dsff(8)=0.281d0
      dsff(8)=dsff(8)+(1.d0-dsff(8))*
     &  dexp(-(dlog(mx)/dlog(71.7d0))**2.97d0)
      dsff(9)=0.101d0
      dsff(9)=dsff(9)+(1.d0-dsff(9))*
     &  dexp(-(dlog(mx)/dlog(57.0d0))**3.1d0)
      dsff(10)=0.00677d0
      dsff(10)=dsff(10)+(1.d0-dsff(10))*
     &  dexp(-(dlog(mx)/dlog(29.3d0))**3.36d0)
      aa(1)=1.d0
      aa(2)=4.d0
      aa(3)=12.d0
      aa(4)=14.d0
      aa(5)=16.d0
      aa(6)=20.d0
      aa(7)=24.d0
      aa(8)=28.d0
      aa(9)=32.d0
      aa(10)=56.d0
      dsntlitlf_s=0.d0
      do 20 i=1,10
        dsntlitlf_s=dsntlitlf_s+
     &  dsff(i)*f(i)*phi(i)*dsntss(mx/m(i),vbar)
     &  *m(i)**3*mx/(mx+m(i))**2
 20   continue
      return
      end
