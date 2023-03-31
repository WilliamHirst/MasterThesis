!
      SUBROUTINE HERMTEST(IN,OUT)
!
!Purpose: Test if complex matrix IN(3,3) is hermitian.
!         Return 1 if not, 0 otherwise.
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX IN(3,3)
      INTEGER OUT
!
      OUT=0
      IF(DIMAG(IN(1,1)).NE.0.D0)OUT=1
      IF(DIMAG(IN(2,2)).NE.0.D0)OUT=1
      IF(DIMAG(IN(3,3)).NE.0.D0)OUT=1
      IF(IN(1,2).NE.CONJG(IN(2,1)))OUT=1
      IF(IN(1,3).NE.CONJG(IN(3,1)))OUT=1
      IF(IN(2,3).NE.CONJG(IN(3,2)))OUT=1
!
      RETURN
      END
