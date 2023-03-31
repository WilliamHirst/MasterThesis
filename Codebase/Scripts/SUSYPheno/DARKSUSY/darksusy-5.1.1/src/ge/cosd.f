      function cosd(x)
      implicit none
      real*8 cosd,x,pi
      parameter (pi=3.141592653589793238d0)
      cosd = cos(x*pi/180.d0)
      return
      end
