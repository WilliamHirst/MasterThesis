**********************************************************************
*** function dshrgacdiffsusy gives the susy dependent term in the
*** flux of gamma-rays with continuum energy spectrum per gev
*** at the energy egam (gev) from neutralino annihilation in the halo.
***
*** dshrgacdiffsusy in unit of gev^-1
*** 
*** the flux in a solid angle delta in the direction psi0 is given by:
***   cm^-2 s^-1 sr^-1 * dshrgacdiffsusy(egam,istat)
***   * dshmjave(cospsi0,delta) * delta (in sr)
***
*** the flux per solid angle in the direction psi0 is given by:
***   cm^-2 s^-1 sr^-1 * dshrgacdiffsusy(egam,istat)
***   * dshmj(cospsi0)
*** 
*** in case of a clumpy halo the factor fdelta has to be added 
***
*** author: joakim edsjo, edsjo@physto.se
*** modified: piero ullio (piero@tapir.caltech.edu) 00-07-13
***           Joakim Edsjo (edsjo@physto.se) 03-01-21, factor of 1/2
***           in annihilation rate added
**********************************************************************

      real*8 function dshrgacdiffsusy(egam,istat)
      implicit none
      include 'dshmcom.h'
      include 'dshacom.h'
      include 'dsprep.h'
      real*8 phiga,dshaloyield,egam,dnsigma
      integer istat

      phiga=dshaloyield(egam,152,istat)
      dnsigma=phiga*(hasv/1.0d-29)*0.5d0 ! JE corr 03-01-21
      dshrgacdiffsusy=1.879d-11*dnsigma*(10.0d0/hamwimp)**2
      end







