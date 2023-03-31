***********************************************************************
*** Wrapper routine kept for backwards compatibility. This routine
*** calls the latest expressions for b->s gamma
***********************************************************************
      subroutine dsbsgammafull(ratio,flag)

      implicit none

      real*8 ratio
      integer flag

c...Old (pre-2007) routines
c      call dsbsgpre2007full(ratio,flag)

c...New (2007) routines
      call dsbsg2007full(ratio,flag)
 
      return
      end

