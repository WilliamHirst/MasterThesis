      subroutine dssuconst
c_______________________________________________________________________
c  useful constants
c  author: paolo gondolo 1994-1999
c  modified: 031105 neutrino's yukawa fixed (pg)
c=======================================================================
      implicit none
      real*8 s12,s23,s13,c12,c13,c23,d13,aux
      complex*16 ed13
      integer i

      call dssuconst_couplings
      call dssuconst_ckm
      call dssuconst_yukawa
      call dssuconst_higgs

      end
