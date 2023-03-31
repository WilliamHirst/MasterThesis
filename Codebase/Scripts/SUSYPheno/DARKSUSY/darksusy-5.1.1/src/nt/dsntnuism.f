      real*8 function dsntnuism(enu)
c...  with enu in gev, the flux is returned in units o
c...  gev^-1 cm^-2 s^-1 sr^-1
      implicit none

      real*8 enu
      real*8 n0,gamma,e0,gammap,n0p
      parameter(n0=3.0d-6,
     &          gamma=1.63d0,
     &          e0=4.7d5,
     &          gammap=1.95d0,
     &          n0p=1.9d-4)

      real*8 e_mu,rd
      integer fltype
      common/lbe_int3/e_mu,rd,fltype

      if (enu.lt.e0) then
        dsntnuism=rd*n0*enu**(-gamma-1.0d0)
      else
        dsntnuism=rd*n0p*enu**(-gammap-1.0d0)
      endif

      return
      end
