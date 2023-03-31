!
      SUBROUTINE SYMMTEST(IN,OUT)
!
!Purpose: Test if real matrix IN(3,3) is hermitian.
!         Return 1 if not, 0 otherwise.
!
      IMPLICIT NONE
!
      DOUBLE PRECISION IN(3,3)
      INTEGER OUT
!
      OUT=0
      IF(IN(1,2).NE.IN(2,1))OUT=1
      IF(IN(1,3).NE.IN(3,1))OUT=1
      IF(IN(2,3).NE.IN(3,2))OUT=1
!
      RETURN
      END
