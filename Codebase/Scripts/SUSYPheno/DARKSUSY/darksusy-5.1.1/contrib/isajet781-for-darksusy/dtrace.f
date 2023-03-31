!
      FUNCTION DTRACE(A)
!
!Purpose: To compute the trace of matrix A
!
      IMPLICIT NONE
!
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION DTRACE
      INTEGER I
!
      DTRACE=0.D0
      DO I=1,3
        DTRACE=DTRACE+A(I,I)
      END DO
!
      RETURN
      END
