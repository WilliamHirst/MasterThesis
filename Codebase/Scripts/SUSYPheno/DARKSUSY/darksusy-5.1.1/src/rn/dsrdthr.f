***********************************************************************
*** subroutine dsrdthr sets up the thresholds before calling
*** dsrdens.
*** author: joakim edsjo, (edsjo@physto.se)
*** date: 98-03-03
***********************************************************************

      subroutine dsrdthr(npart,mgev,nth,thgev)
      implicit none
      include 'dsmssm.h'

c----------------------------------------------------------------------

      integer npart,nth
      real*8 mgev(20),thgev(20)

c----------------------------------------------------------------------

      nth=0

c...ww-threshold
      if (mass(kw).gt.mgev(1)) then
        nth=nth+1
        thgev(nth)=2.0d0*mass(kw)
      endif

c...zz-threshold
      if (mass(kz).gt.mgev(1)) then
        nth=nth+1
        thgev(nth)=2.0d0*mass(kz)
      endif

c...tt-bar-threshold
      if (mass(kt).gt.mgev(1)) then
        nth=nth+1
        thgev(nth)=2.0d0*mass(kt)
      endif

c...note that coannihilation thresholds are automatically added in
c...dsrdstart (called from dsrdens).

      end
