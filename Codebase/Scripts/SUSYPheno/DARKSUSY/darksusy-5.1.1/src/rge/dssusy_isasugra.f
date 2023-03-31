      subroutine dssusy_isasugra(unphys,valid)
c--------------------------------------------------------------------
c     replacement for dssusy for using ISASUGRA RGE evolution
c     author: E.A. Baltz, 2001 eabaltz@alum.mit.edu
c====================================================================
      implicit none
      integer unphys,valid
            
c...  run RGEs and set mass spectrum and mixing
      call dsrge_isasugra(unphys,valid)
      if (valid.gt.0) return
      if (unphys.ne.0) return

c...  set up models with masses and couplings
      call dsmodelsetup_isasugra

c...  Prepare for rate calculations
      if (unphys.eq.0.and.valid.le.0) then
        call dsprep
      endif

      end
