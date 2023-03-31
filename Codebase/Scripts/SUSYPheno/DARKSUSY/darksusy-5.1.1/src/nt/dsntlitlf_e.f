       real*8 function dsntlitlf_e(mx,vbar)
c_______________________________________________________________________
c
c    dsntlitlf_e is used by capearth to calculate the capture rate
c    in the earth.
c    mx is the neutralino mass in gev
c    written by l. bergstrom 1995-12-12
c
c=======================================================================

       implicit none
       include 'dsmpconst.h'
       real*8 mx,aa(11),dsntse,vbar,x,a
       real*8 m(11),f(11),phi(11),mn,dsff(11)
       integer i
       mn=(m_p+m_n)/2.0d0
       m(1)=16.d0*mn
       m(2)=28.d0*mn
       m(3)=24.d0*mn
       m(4)=56.d0*mn
       m(5)=40.d0*mn
       m(6)=30.d0*mn
       m(7)=23.d0*mn
       m(8)=32.d0*mn
       m(9)=59.d0*mn
       m(10)=27.d0*mn
       m(11)=52.d0*mn

c...These are the values in jkg
       f(1)=0.3d0
       f(2)=0.15d0
       f(3)=0.14d0
       f(4)=0.3d0
       f(5)=0.015d0
       f(6)=0.011d0
       f(7)=0.004d0
       f(8)=0.05d0
       f(9)=0.03d0
       f(10)=0.0d0
       f(11)=0.0d0

       phi(1)=1.2d0
       phi(2)=1.2d0
       phi(3)=1.2d0
       phi(4)=1.6d0
       phi(5)=1.2d0
       phi(6)=1.2d0
       phi(7)=1.2d0
       phi(8)=1.6d0
       phi(9)=1.6d0
       phi(10)=1.0d0
       phi(11)=1.0d0

c...These are the new values based on McDonough, Treatise on Geochemistry,
c...Vol 2, Elsevier, 2003.
       f(1)=0.298d0
       f(2)=0.162d0
       f(3)=0.154d0
       f(4)=0.319d0
       f(5)=0.0171d0
       f(6)=0.00071d0
       f(7)=0.00183d0
       f(8)=0.0063d0
       f(9)=0.0181d0
       f(10)=0.0159d0
       f(11)=0.0047d0

c...old 1.6 core, 1.2 mantle
       phi(1)=1.2d0
       phi(2)=1.24d0
       phi(3)=1.2d0
       phi(4)=1.546d0
       phi(5)=1.2d0
       phi(6)=1.56d0
       phi(7)=1.2d0
       phi(8)=1.59d0
       phi(9)=1.57d0
       phi(10)=1.20d0
       phi(11)=1.44d0

c...new full average
       phi(1)=1.20d0
       phi(2)=1.25d0
       phi(3)=1.20d0
       phi(4)=1.54d0
       phi(5)=1.2d0
       phi(6)=1.56d0
       phi(7)=1.20d0
       phi(8)=1.58d0
       phi(9)=1.56d0
       phi(10)=1.20d0
       phi(11)=1.45d0
       

       do 10 i=1,11
         dsff(i)=1.d0
 10    continue
       x=mx/m(4)
       a=1.5d0*x/(x-1.d0)**2*(13.2d0**2/vbar**2)
       dsff(4)=1.d0-0.26d0*a/(1.d0+a)

       aa(1)=16.d0
       aa(2)=28.d0
       aa(3)=24.d0
       aa(4)=56.d0
       aa(5)=40.d0
       aa(6)=30.d0
       aa(7)=23.d0
       aa(8)=32.d0
       aa(9)=59.d0
       aa(10)=27.d0
       aa(11)=52.d0

       dsntlitlf_e=0.d0

c...Now calculate f, Eq. (9-28) in jkg
       do 20 i=1,11
         dsntlitlf_e=dsntlitlf_e+
     &   dsff(i)*f(i)*phi(i)*dsntse(mx/m(i),vbar)
     &   *m(i)**3*mx/(mx+m(i))**2
 20    continue
       return
       end
