      function dsbsgft1(x,flag)

***********************************************************************
* Function F^t_1(x) p. 11 in Gambino and Misiak,                      *
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
      real*8 dsbsgc81h2,dsbsgc81chisusy,dsbsgc81wsusy
      real*8 dsbsgc81phi1susy,dsbsgc81phi2susy
      real*8 dsbsgft1,sm
      real*8 ddilog     


c The Standard Model contribution, i.e. p.11 in Nucl. Phys. B611 
      sm=(ddilog(1.d0-1.d0/x))*
     &   (4.d0*x**4-40.d0*x**3-41.d0*x**2-x)/
     &   (3.d0*(1.d0-x)**4)
     &   +log(x)*
     &   (-144.d0*x**4+3177.d0*x**3+3661.d0*x**2+250.d0*x-32.d0)/
     &   (108.d0*(1.d0-x)**5)
     &   +(-247.d0*x**4+11890.d0*x**3+31779.d0*x**2-2966.d0*x+1016.d0)/
     &   (648.d0*(1.d0-x)**4)

      dsbsgft1=sm

c Include the susy contribution

      if (flag.eq.1) then
        dsbsgft1=sm-2.d0*(dsbsgc81h2()+dsbsgc81chisusy()
     &              +dsbsgc81wsusy()+dsbsgc81phi1susy()
     &              +dsbsgc81phi2susy())
      endif

      return
      end


