**********************************************************************
*** subroutine dshrgalsusy gives the susy dependent term in the
*** flux of monoenergetic gamma-rays from neutralino annihilation 
*** in the halo.
***
*** dshrgalsusy is dimensionless
*** 
*** the flux in a solid angle delta in the direction psi0 is given by:
***   cm^-2 s^-1 sr^-1 * gagarate (or gazrate)
***   * dshmjave(cospsi0,delta) * delta (in sr)
***
*** the flux per solid angle in the direction psi0 is given by:
***   cm^-2 s^-1 sr^-1 * gagarate (or gazrate)
***   * dshmj(cospsi0)
*** 
*** in case of a clumpy halo the factor fdelta has to be added 
***
*** author: joakim edsjo, edsjo@physto.se
*** modified: piero ullio (piero@tapir.caltech.edu) 00-07-13
***           Joakim Edsjo (edsjo@physto.se) 03-01-21, factor of 1/2
***           in annihilation rate added
**********************************************************************

      subroutine dshrgalsusy(gagarate,gazrate)
      implicit none
      real*8 gagarate,gazrate,dnsgaga,dnsgaz,dssigmav
      include 'dshmcom.h'
      include 'dshacom.h'
      include 'dsidtag.h'
      include 'dsprep.h'

c-----------------------------------------------------------------------

      if (.not.dshasetupcalled) then
        write(*,*) 'DS error in dshrgalsusy: dshasetup must be called',
     &    ' before any halo rate calculation'
        write(*,*) 'begins for every new model. program stopping.'
        stop
      endif

c... gamma-gamma
c      dnsgaga=2.0d0*(dssigmav(28)/1.0d-29)
      dnsgaga=2.0d0*(habr(28)*hasv/1.0d-29)
      gagarate=1.879d-11*dnsgaga*(10.0d0/hamwimp)**2*0.5d0 ! JE Corr 03-01-21
c... z gamma
      dnsgaz=(habr(29)*hasv/1.0d-29)
      gazrate=1.879d-11*dnsgaz*(10.0d0/hamwimp)**2*0.5d0 ! JE Corr 03-01-21
      end





