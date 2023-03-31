*****************************************************************************
*** auxiliary routine called by dsIBwwdxdy
*** author: Torsten Bringmann, 2007-02-16
*****************************************************************************

      real*8 function dsIBwwdxdy_2(x,y,m0,mw,mc1,mc2)
      implicit none
      real*8 x,y,m0,mw,mc1,mc2

      dsIBwwdxdy_2 =
     -   (-128*m0**2*(mw**2 + 2*m0**2*(-1 + x))*
     -    (mw**4*(-1 + x) + 8*m0**4*(-1 + x)*(x**2 - 2*x*y + 2*y**2) + 
     -      4*m0**2*mw**2*(2*x**2 + 2*y - x*(1 + 2*y))))/
     -  (mw**2*(-2*mc1**2 + 3*mw**2 + m0**2*(-2 + 4*x - 4*y))*
     -    (-2*mc2**2 + 3*mw**2 + m0**2*(-2 + 4*x - 4*y))*
     -    (-2*mc1**2 + mw**2 + m0**2*(-2 + 4*y))*
     -    (-2*mc2**2 + mw**2 + m0**2*(-2 + 4*y)))/2
      return
      end   ! dsIBwwdxdy_2

