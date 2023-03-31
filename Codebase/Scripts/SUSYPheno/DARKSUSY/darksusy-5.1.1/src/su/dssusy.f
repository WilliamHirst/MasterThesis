      subroutine dssusy(unphys,hwarning)
c_______________________________________________________________________
c  set up global variables for the supersymmetric model routines
c  and prepares the rate routines for a new model.
c  author: Joakim Edsjo, edsjo@physto.se, 2001
c=======================================================================
      implicit none
      integer unphys,hwarning

c-----------------------------------------------------------------------

c... Calculate masses, couplings, etc
      call dsmodelsetup(unphys,hwarning)

c... Prepare for rate calculations
      if (unphys.eq.0.and.hwarning.eq.0) then
        call dsprep
      endif

      end


