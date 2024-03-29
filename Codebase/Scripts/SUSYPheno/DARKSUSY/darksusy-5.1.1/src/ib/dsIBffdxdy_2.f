*****************************************************************************
*** auxiliary routine called by dsIBffdxdy
*** author: Torsten Bringmann, 2007-07-05
*****************************************************************************

      real*8 function dsIBffdxdy_2(x,y,m0,ml,msl1,msl2)
      implicit none
      real*8 x,y,m0,ml,msl1,msl2

      dsIBffdxdy_2 = 
     - (-8*m0**2*(ml**12*(-1 + x) + 
     -      2048*m0**12*(-1 + x)*(x - y)**2*y**2*
     -       (x**2 - 2*x*y + 2*y**2) - 
     -      512*m0**10*ml**2*y*
     -       (2*x**5 - 2*x**4*(6 + 5*y) + x**3*(17 + 22*y + 24*y**2) - 
     -         x**2*(10 + 17*y + 28*y**2 + 28*y**3) - 2*(y + 6*y**4) + 
     -         2*x*(1 + 5*y + 15*y**3 + 6*y**4)) + 
     -      4*m0**2*ml**6*(8*msl1**2*msl2**2*(2 - 2*x + x**2) - 
     -         8*ml**2*(msl1**2 + msl2**2)*(2 - 2*x + x**2) + 
     -         ml**4*(16 + 10*x**2 + 6*y - x*(19 + 6*y))) + 
     -      8*m0**4*ml**4*(16*msl1**2*msl2**2*
     -          (x**3 - 4*y - 2*x**2*y + x*(2 + 4*y)) - 
     -         4*ml**2*(msl1**2 + msl2**2)*
     -          (5*x**3 - 2*(1 + 8*y) - x**2*(5+8*y) + 2*x*(7+8*y))
     -          + ml**4*(27*x**3 - x**2*(47 + 54*y) - 
     -            2*(8 + 32*y + 15*y**2) + x*(80 + 94*y + 30*y**2))) + 
     -      32*m0**6*ml**2*(-16*msl1**2*msl2**2*(2 - 2*x + x**2)*
     -          (x - y)*y - 4*ml**2*(msl1**2 + msl2**2)*
     -          (x**4 + 4*y*(1 + 2*y) - x**3*(1 + 6*y) + 
     -            2*x**2*(2 + 9*y + 2*y**2) - 2*x*(1+10*y+4*y**2))
     -          + ml**4*(2 + 10*x**4 + 32*y + 32*y**2 + 40*y**3 - 
     -            2*x**3*(11 + 24*y) + x**2*(49 + 140*y + 64*y**2) - 
     -            2*x*(13 + 64*y + 46*y**2 + 20*y**3))) + 
     -      128*m0**8*ml**2*(4*(msl1**2 + msl2**2)*
     -          (-2 + 6*x - 5*x**2 + x**3)*(x - y)*y + 
     -         ml**2*(x**5 - x**4*(3 + 16*y) + 
     -            x**3*(9 + 72*y + 38*y**2) - 2*y*(2 + 8*y + 15*y**3) - 
     -            2*x**2*(4 + 41*y + 41*y**2 + 26*y**3) + 
     -            x*(2 + 36*y + 48*y**2 + 60*y**3 + 30*y**4)))))/
     -  ((3*ml**2 - 2*msl1**2 + m0**2*(-2 + 4*x - 4*y))*
     -    (3*ml**2 - 2*msl2**2 + m0**2*(-2 + 4*x - 4*y))*
     -    (ml**2 + 4*m0**2*(x - y))**2*(ml**2 - 4*m0**2*y)**2*
     -    (ml**2 - 2*msl1**2 + m0**2*(-2 + 4*y))*
     -    (ml**2 - 2*msl2**2 + m0**2*(-2 + 4*y)))

      return
      end   ! dsIBffdxdy_2
