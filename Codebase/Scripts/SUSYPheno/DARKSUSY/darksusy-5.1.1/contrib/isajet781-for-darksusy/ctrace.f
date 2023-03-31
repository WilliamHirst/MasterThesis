!
      FUNCTION CTRACE(A)
!
!Purpose: To compute the trace of matrix A
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX A(3,3),CTRACE
      INTEGER I
!
      CTRACE=(0.D0,0.D0)
      DO I=1,3
        CTRACE=CTRACE+A(I,I)
      END DO
!
      RETURN
      END
