      subroutine dsvertx
c_______________________________________________________________________
c  some couplings used in DarkSUSY
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994--
c  history:
c    951110 complex vertex constants
c    970213 joakim edsjo
c    990724 paolo gondolo trilinear higgs and goldstone couplings
c=======================================================================
c
c  vertices included:
c     see individual routines dsvertx1 and dsvertx3
c
      implicit none
c  scattering, annihilation, and chargino coannihilation
      call dsvertx1
c  sfermion coannihilation
      call dsvertx3
      return
      end
