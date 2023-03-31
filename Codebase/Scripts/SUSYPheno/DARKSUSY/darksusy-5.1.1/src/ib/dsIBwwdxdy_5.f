*****************************************************************************
*** auxiliary routine called by dsIBwwdxdy
*** author: Torsten Bringmann, 2007-02-16
*****************************************************************************

      real*8 function dsIBwwdxdy_5(x,y,m0,mw,mc1,mc2)
      implicit none
      real*8 x,y,m0,mw,mc1,mc2

      dsIBwwdxdy_5 = 
     -   (256*m0**2*mc1*mc2*(mw**6 + m0**2*mw**4*(-1 + 5*x - 8*y) - 
     -      16*m0**6*(-1 + x)*(x - y)*y + 
     -      4*m0**4*mw**2*(2*x**2 + 2*y*(1 + 2*y) - x*(1 + 6*y))))/
     -  (mw**4*(-2*mc1**2 + 3*mw**2 + m0**2*(-2 + 4*x - 4*y))*
     -    (-2*mc2**2 + 3*mw**2 + m0**2*(-2 + 4*x - 4*y))*
     -    (-2*mc1**2 + mw**2 + m0**2*(-2 + 4*y))*
     -    (-2*mc2**2 + mw**2 + m0**2*(-2 + 4*y)))
      return
      end   ! dsIBwwdxdy_5

