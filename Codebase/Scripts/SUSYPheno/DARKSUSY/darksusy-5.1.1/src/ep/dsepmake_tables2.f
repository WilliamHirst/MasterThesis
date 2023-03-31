************************************************************************
*** positron propagation routines.
*** author: e.a. baltz (eabaltz@astron.berkeley.edu)
*** modified slightly by joakim edsjo (edsjo@physto.se)
*** date: jun-02-98
*** modified: jun-09-98, 2002-11-19
***   2004-01-26 (better r integration)
************************************************************************


************************************************************************
      subroutine dsepmake_tables2
*** creates auxiliary tables
************************************************************************
      implicit none
      include 'dsepcom.h'
      real*8 lne,egev,u
      integer i

      k0tau = k27*tau16
      ametric=1.0d0/3.0d0**alphexp
      av=ametric/(1.0d0-alphexp)

c...Changed to allow for lower energies 2002-11-19
      if (dsepdiffmethod.eq.1) then
         do 60 i=1,13001
            lne=real(i-1)*0.001d0-3.0d0  ! down to e-3 GeV (~0.05 GeV)
            egev=exp(lne)
            u=1.0d0/egev
            etable(i)=egev
            vtable(i)=u+av*u**(1.0d0-alphexp)
            wtable(i)=egev**2/(1.0d0+ametric*egev**alphexp)
 60      continue
      endif
      return
      end
