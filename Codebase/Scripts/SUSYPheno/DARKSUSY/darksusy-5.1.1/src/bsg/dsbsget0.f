      function dsbsget0(x,flag)

***********************************************************************
* Function E^t_0(x) p. 11 in Gambino and Misiak,                      *
* Nucl. Phys. B611 (2001) 338                                         *
* x must be a positive number                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-03-13                   *
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
      real*8 dsbsgc41h2,dsbsgc41susy
      real*8 dsbsget0,sm

c The Standard Model contribution, i.e. p. 11 in Nucl. Phys. B611 
      sm=log(x)*(-9.d0*x**2+16.d0*x-4.d0)/(6.d0*(x-1)**4)
     &   +(7.d0*x**3+21.d0*x**2-42.d0*x-4.d0)/
     &   (36.d0*(x-1.d0)**3)

      dsbsget0=sm

c Include the susy contribution
      
      if (flag.eq.1) then
        dsbsget0=sm+dsbsgc41h2()+dsbsgc41susy()
      endif

      return
      end


