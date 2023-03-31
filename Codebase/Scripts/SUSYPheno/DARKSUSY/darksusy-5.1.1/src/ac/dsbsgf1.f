      function dsbsgf1(x)
c_______________________________________________________________________
c  function in bertolini et al, nucl phys b353 (1991) 591
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c=======================================================================
      implicit none
      real*8 dsbsgf1,x,y,res
      y = x-1.d0
      if (abs(y).lt.0.01) then
         res = (1.d0-0.6*y+0.4*y**2)/24.d0
      else if (x.lt.1.d-4) then
         res = 1.d0/6.d0
      else
         res = (x**3-6*x**2+3*x+2+6*x*log(x))/(12.*(x-1)**4)
      endif
      dsbsgf1 = res
      end
