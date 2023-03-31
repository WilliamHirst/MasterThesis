**********************************************************************
*** function dshrgacsusy gives the susy dependent term in the
*** flux of gamma-rays with continuum energy spectrum above the 
*** threshold egath (gev) from neutralino annihilation in the halo.
***
*** dshrgacsusy is dimensionless
*** 
*** the flux in a solid angle delta in the direction psi0 is given by:
***   cm^-2 s^-1 sr^-1 * dshrgacsusy(egath,istat)
***   * dshmjave(cospsi0,delta) * delta (in sr)
***
*** the flux per solid angle in the direction psi0 is given by:
***   cm^-2 s^-1 sr^-1 * dshrgacsusy(egath,istat)
***   * dshmj(cospsi0)
*** 
*** in case of a clumpy halo the factor fdelta has to be added 
***
*** author: joakim edsjo, edsjo@physto.se
*** modified: piero ullio (piero@tapir.caltech.edu) 00-07-13
***           Joakim Edsjo (edsjo@physto.se) 03-01-21, factor of 1/2
***           in annihilation rate added
**********************************************************************

      real*8 function dshrgacsusy(egath,istat)
      implicit none
      include 'dshmcom.h'
      include 'dshacom.h'
      include 'dsprep.h'
      real*8 phiga,dshaloyield,egath,dnsigma
      integer istat
      phiga=dshaloyield(egath,52,istat)
      dnsigma=phiga*(hasv/1.0d-29)*0.5d0 ! JE Corr 03-01-21
      dshrgacsusy=1.879d-11*dnsigma*(10.0d0/hamwimp)**2
      end
