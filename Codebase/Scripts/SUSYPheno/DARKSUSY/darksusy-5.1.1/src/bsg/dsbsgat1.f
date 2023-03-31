      function dsbsgat1(x,flag)

***********************************************************************
* Function A^t_1(x) p. 11 in Gambino and Misiak,                      *
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
      real*8 dsbsgc71h2,dsbsgc71chisusy,dsbsgc71wsusy
      real*8 dsbsgc71phi1susy,dsbsgc71phi2susy
      real*8 dsbsgat1,sm
      real*8 ddilog     


c      write(*,*) 'dsbsgat1: x=',x,1.-1./x
c The Standard Model contribution, i.e. p.11 in Nucl. Phys. B611 
      sm=(ddilog(1.d0-1.d0/x))*
     &   (32.d0*x**4+244.d0*x**3-160.d0*x**2+16.d0*x)/
     &   (9.d0*(1.d0-x)**4)
     &   +log(x)*
     &   (-774.d0*x**4-2826.d0*x**3+1994.d0*x**2-130.d0*x+8.d0)/
     &   (81.d0*(1.d0-x)**5)
     &   +(-94.d0*x**4-18665.d0*x**3+20682.d0*x**2-9113.d0*x+2006.d0)/
     &   (243.d0*(1.d0-x)**4)

      dsbsgat1=sm

c Include the susy contribution

      if (flag.eq.1) then
        dsbsgat1=sm-2.d0*(dsbsgc71h2()+dsbsgc71chisusy()
     &              +dsbsgc71wsusy()+dsbsgc71phi1susy()
     &              +dsbsgc71phi2susy())
      endif

      return
      end


