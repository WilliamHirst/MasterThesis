************************************************************************
*** positron propagation routines.
*** author: e.a. baltz (eabaltz@astron.berkeley.edu)
*** modified slightly by joakim edsjo (edsjo@physto.se)
*** date: jun-02-98
*** modified: jun-09-98
***           jul-06-99 paolo gondolo - calls to dshunt, ee, vv
*** Modified: Joakim Edsjo (edsjo@physto.se) 03-01-21, factor of 1/2
***           in annihilation rate added
************************************************************************



************************************************************************
*** real*8 function dsepspec calculates the differential positron
*** spectrum from neutralino annihilation in the halo.
*** input:  e - positron kinetic energy in gev
***         mchi - neutralino mass in gev
***         sigvline - <sigma v>_e+e- in cm^3 s^-1
***         sigvdnde(eep) - <sigma v> dn/de in cm^3 s^-1 gev^-1
***         ee - r8 function that gives energy as a fcn of v
***         vv - r8 function that gives v as function of energy
***         metric - r8 function that gives metric as fcn of v
*** output: spectrum in units of cm^-2 s^-1 sr^-1 gev^-1
************************************************************************

      real*8 function dsepspec(e,mchi,sigvline,sigvdnde,ee,vv,metric)
      implicit none

      include 'dsge.h'
      include 'dsepcom.h'
      include 'dsprep.h'
      real*8 e,v,vmin,mchi,dv,vprime,deltav,vv,metric,ee
      real*8 sigvline,nsq,sigvdnde,vep
      integer tabindx
      data tabindx/0/
      external sigvdnde,ee,vv,metric
      integer i
c-----------------------------------------------------------------------

      vmin=vv(mchi)
      if(e.ge.mchi) then
         dsepspec=0.0d0
         return
      endif

      nsq=(n_c/mchi)**2*0.5d0 ! JE Corr 03-01-21
      v=vv(e)
      dsepspec=0.0d0
      dv=(v-vmin)*vtol
      vprime=vmin+dv*0.5d0
 10   deltav=v-vprime
      if (dsepdiffmethod.eq.1) then
         call dshunt(vtable,13001,vprime,tabindx)
      endif
      i=max(1+int(1000.0d0*(log(deltav)+19.0d0)),1)
      dsepspec=dsepspec+table(i)*dv*
     +     sigvdnde(ee(vprime,tabindx))*metric(vprime,tabindx)
      vprime=vprime+dv
      if (vprime.lt.v) goto 10
c     line contribution
      i=max(1+int(1000.0d0*(log(v-vmin)+19.0d0)),1)
      dsepspec=dsepspec+sigvline*table(i)
      dsepspec=dsepspec*nsq*tau16*1.0d16/e**2


c...convert density to flux
      vep=2.99792458d10
      dsepspec=dsepspec*vep/(4.0d0*pi)

      return
      end
