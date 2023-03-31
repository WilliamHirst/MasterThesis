      function dsbsgat0(x,flag)

***********************************************************************
* Function A^t_0(x) in (2.7) of Gambino and Misiak,                   *
* Nucl. Phys. B611 (2001) 338                                         *
* x must be a positive number                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-03-12                   *
* updated by Mia Schelke 2003-04-10 to include the susy contribution  *
* as explained in eq. (5.1)                                           *
* (the references used for the susy contributions can be found in     *
*  the different fortran codes)                                       *
* Input: flag:  0 = only standard model                               *
*               1 = standard model plus SUSY corrections              *
***********************************************************************

      implicit none
      integer flag
      real*8 x
      real*8 dsbsgc70h2,dsbsgc70susy
      real*8 dsbsgat0,sm

c The Standard Model contribution, i.e. eq.(2.7) of Nucl. Phys. B611 
      sm=log(x)*(-3.d0*x**3+2.d0*x**2)/(2.d0*(x-1)**4)
     &   +(-22.d0*x**3+153.d0*x**2-159.d0*x+46.d0)/
     &   (36.d0*(x-1.d0)**3)

      dsbsgat0=sm

c Include the susy contribution

      if (flag.eq.1) then
        dsbsgat0=sm-2.d0*(dsbsgc70h2()+dsbsgc70susy())
      endif

      return
      end


