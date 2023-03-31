
*************************
      real*8 function dsntepfunc(r)
      implicit none

      real*8 r,dsntearthmass,gn
      parameter(gn=6.67259d-11) ! m^3 kg^-1 s^-1

      dsntepfunc=dsntearthmass(r)*gn/max(r,100.0d0)**2

      return
      end
