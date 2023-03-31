      function dsbsgf3(x)
c_______________________________________________________________________
c  function in bertolini et al, nucl phys b353 (1991) 591
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c=======================================================================
      implicit none
      real*8 dsbsgf3,x,y,res
      y = x-1.d0
      if (abs(y).lt.0.01) then
         res = (1.d0-0.75d0*y+0.6d0*y**2)/3.d0
      else if (x.lt.1.d-4) then
         res = -1.5d0-log(x)
      else
         res = (x**2-4*x+3+2*log(x))/(2.*(x-1)**3)
      endif
      dsbsgf3 = res
      end
