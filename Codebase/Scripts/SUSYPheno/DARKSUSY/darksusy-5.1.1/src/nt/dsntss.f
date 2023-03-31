       real*8 function dsntss(x,vbar)
c----------------------------------------------------------------------
c       dsntss used by capsun and litlf_s
c       x=mx/m(i)
c       lars bergstrom 1995-12-12
c----------------------------------------------------------------------

       implicit none
       real*8 x,vbar,vesc,a
       vesc=1156.d0
       a=1.5d0*x/(x-1.d0)**2*(vesc**2/vbar**2)
       dsntss=(a**1.5d0/(1+a**1.5d0))**(1.d0/1.5d0)
       return
       end
