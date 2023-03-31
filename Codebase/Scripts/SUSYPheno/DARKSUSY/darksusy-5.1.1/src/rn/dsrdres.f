***********************************************************************
*** subroutine dsrdres sets up the resonances before calling
*** dsrdens.
*** author: joakim edsjo, (edsjo@physto.se)
*** date: 98-03-03
***********************************************************************

      subroutine dsrdres(npart,mgev,nres,rgev,rwid)
      implicit none
      include 'dsmssm.h'

c----------------------------------------------------------------------

      integer npart,nres
      real*8 mgev(20),rgev(20),rwid(20)

c----------------------------------------------------------------------

      nres=0

c...z resonance
      if (mass(kz).gt.mgev(1)*2.0d0) then
        nres=nres+1
        rgev(nres)=mass(kz)
        rwid(nres)=width(kz)
      endif

c...h1 resonance
      if (mass(kh1).gt.mgev(1)*2.0d0) then
        nres=nres+1
        rgev(nres)=mass(kh1)
        rwid(nres)=width(kh1)
      endif

c...h2 resonance
      if (mass(kh2).gt.mgev(1)*2.0d0) then
        nres=nres+1
        rgev(nres)=mass(kh2)
        rwid(nres)=width(kh2)
      endif

c...h3 resonance
      if (mass(kh3).gt.mgev(1)*2.0d0) then
        nres=nres+1
        rgev(nres)=mass(kh3)
        rwid(nres)=width(kh3)
      endif

c...coannihilation-resonances
      if (npart.gt.1) then
        if (mass(khc).gt.mgev(1)*2.0d0) then   ! sufficient condition
          nres=nres+1
          rgev(nres)=mass(khc)
          rwid(nres)=width(khc)
        endif
        
        if (mass(kw).gt.mgev(1)*2.0d0) then   ! sufficient condition
          nres=nres+1
          rgev(nres)=mass(kw)
          rwid(nres)=width(kw)
        endif
        
      endif


      end
