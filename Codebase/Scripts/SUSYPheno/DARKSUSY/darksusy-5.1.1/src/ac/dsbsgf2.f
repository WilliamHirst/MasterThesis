      function dsbsgf2(x)
c_______________________________________________________________________
c  function in bertolini et al, nucl phys b353 (1991) 591
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c=======================================================================
      implicit none
      real*8 dsbsgf2,x,y,res
      y = x-1.d0
      if (abs(y).lt.0.01) then
         res = (1.d0-0.4*y+0.2*y**2)/24.d0
      else if (x.lt.1.d-4) then
         res = 1.d0/12.d0-x/6.d0
      else
         res = (2*x**3+3*x**2-6*x+1-6*x**2*log(x))/(12.*(x-1)**4)
      endif
      dsbsgf2 = res
      end
