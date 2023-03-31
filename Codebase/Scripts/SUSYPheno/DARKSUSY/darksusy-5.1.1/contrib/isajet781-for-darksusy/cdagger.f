!
      SUBROUTINE CDAGGER(A,DAGA)
!
!Purpose: To compute the dagger of matrix A
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX A(3,3),DAGA(3,3)
      INTEGER I,J
!
      DO I=1,3
        DO J=1,3
          DAGA(I,J)=CONJG(A(J,I))
        END DO
      END DO
!
      RETURN
      END
