***********************************************************************
*** dsmhcfull returns the integrand for the integration over the SM 
*** momenta that appears in the collision term. Called by dsmhboltz.
*** itegrand is multiplied by 1d20 to work better with the integration
*** routine
***
*** author: torsten bringmann (troms@physto.se), 2010-01-23
***********************************************************************

      real*8 function dsmhcint(kint)
      implicit none

      include 'dsge.h'
      include 'dsdirver.h'
      include 'dsio.h'
 

      include 'dsmssm.h'
      include 'dsmhcom.h'
      include 'dsidtag.h'

      real*8 m0,mf,kint,omega,tmpres
      real*8 ewt,dsmhm2,cn
      integer SMtype,n
c      real*8 pi
c      parameter (pi=3.141592653589793238d0)

      SMtype = 1
      tmpres = 0d0

      m0=mass(kn(1))

 10   if (SMtype.le.3) mf=0d0
      if (SMtype.eq.4) mf=mass(ke)
      if (SMtype.eq.5) mf=mass(kmu)
      if (SMtype.eq.6) mf=mass(ktau)
      if (SMtype.eq.7) mf=mass(ku)
      if (SMtype.eq.8) mf=mass(kd)
      if (SMtype.eq.9) mf=mass(ks)
      if (kint.lt.mf/10.) goto 20   
      if (SMtype.ge.7.and.Tint.lt.tqcdmax) goto 20
      omega = sqrt(mf**2+kint**2)
      ewt = exp(omega/Tint)

c... this would be the expression for the simplified amplitude
c      call dsmhm2simp(m0,SMtype,cn,n)
c      tmpres=tmpres+cn*(omega/m0)**2/omega*(ewt/(ewt+1.))/(ewt+1.)

      tmpres = tmpres+
     -    dsmhm2(omega,m0,SMtype)*ewt/(omega*(ewt+1.)**2)
  
 20   SMtype = SMtype+1
      if (SMtype.le.9) goto 10


      dsmhcint = tmpres*kint**5/(6*(2*pi)**3*m0**4*Tint)
      dsmhcint=dsmhcint*1.d20 ! rescale for integration to work fine

      return
      end

