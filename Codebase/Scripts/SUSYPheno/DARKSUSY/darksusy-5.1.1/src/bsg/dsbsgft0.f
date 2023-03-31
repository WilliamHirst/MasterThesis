      function dsbsgft0(x,flag)

***********************************************************************
* Function F^t_0(x) in (2.7) of Gambino and Misiak,                   *
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
      real*8 dsbsgc80h2,dsbsgc80susy
      real*8 dsbsgft0,sm

c The Standard Model contribution, i.e. eq.(2.7) of Nucl. Phys. B611 
      sm=log(x)*3.d0*x**2/(2.d0*(x-1)**4)
     &   +(-5.d0*x**3+9.d0*x**2-30.d0*x+8.d0)/
     &   (12.d0*(x-1.d0)**3)

      dsbsgft0=sm

c Include the susy contribution

      if (flag.eq.1) then
        dsbsgft0=sm-2.d0*(dsbsgc80h2()+dsbsgc80susy())
      endif

      return
      end


