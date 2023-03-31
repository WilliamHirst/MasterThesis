      real*8 function dsntdkgtot10(mx,sigsi,sigsd)
c----------------------------------------------------------------------
c       gtot used in damour-krauss calculation
c       units of 10**(-10) gev**(-2)
c       mx: neutralino mass
c       sigsi: spin-indep cross section in cm**2
c       sigsd: spin-dep   cross section in cm**2
c       lars bergstrom 1998-09-15
c----------------------------------------------------------------------
      implicit none
      include 'dsmpconst.h'

      real*8 mx,aa(10),fact,ntdkgtot
      real*8 m(10),f(10),mn,sig_a,dsntdkka,sigsi,sigsd
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
      f(1)=0.7095d0  ! from dk-paper
      f(2)=0.2715d0  ! from dk-paper
      f(3)=3.87d-3   ! the following from jkg review
      f(4)=9.4d-4
      f(5)=8.55d-3
      f(6)=1.51d-3
      f(7)=7.39d-4
      f(8)=8.13d-4
      f(9)=4.65d-4
      f(10)=1.46d-3
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
      ntdkgtot=0.d0
      do 20 i=1,10
      sig_a=1.*aa(i)**2*(mx*m(i))**2/(mx+m(i))**2/(mx*mn)**2*(mx+mn)**2
              ! to go from sigsi to sigma_a
              !  see d-k (7.1)
        ntdkgtot=ntdkgtot+
     &  f(i)*sig_a/m(i)*dsntdkka(aa(i),mx)
 20   continue
	  fact=1./0.389d-27  ! cm**2 to gev**-2
      dsntdkgtot10=sigsi*fact*ntdkgtot*1.d10  !see disc. after (5.16)
c  add spin-dependent piece:
      dsntdkgtot10=dsntdkgtot10+sigsd*fact*1.d10*f(1)/m(1)*
     &  dsntdkka(aa(1),mx)
      return
      end
