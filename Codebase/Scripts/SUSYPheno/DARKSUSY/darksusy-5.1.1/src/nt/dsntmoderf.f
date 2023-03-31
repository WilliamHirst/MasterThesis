      real*8 function dsntmoderf(x)
c =================================================================
c error function 
c modified by l. bergstrom 98-09-15
c modified by p. gondolo 2000-07-19
c see test output below
c used for damour-krauss calculations
c =================================================================
      implicit none
      real*8 erf,x
      dsntmoderf=dsqrt(3.141592d0)/2.d0*erf(x)
      return
      end
