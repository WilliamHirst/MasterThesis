      real*8 function dsntnusun(enu)
c... with enu in gev, the flux is returned in units of gev^-1 cm^-2 s^-1
      implicit none

      real*8 enu
      real*8 n0,gamma,a,e0,gammap,n0p
      parameter(n0=1.3d-5,
     &          gamma=1.98d0,
     &          a=8.5d-6,
     &          e0=3.0d6,
     &          gammap=2.38d0,
     &          n0p=5.1d-3)

      if (enu.lt.e0) then
        dsntnusun=n0*enu**(-gamma-1.0d0)/(1.0d0+a*enu)
      else
        dsntnusun=n0p*enu**(-gammap-1.0d0)/(1.0d0+a*enu)
      endif

      return
      end
