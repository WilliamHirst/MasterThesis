C----------------------------------------------------------------------
      REAL*8 FUNCTION TRACE(X,I)
C----------------------------------------------------------------------
C
C     Takes the trace of 250x3x3 matrix X with respect to last two indices.
C
      IMPLICIT NONE
      REAL*8 X(250,3,3)
      INTEGER I
      
      TRACE=X(I,1,1)+X(I,2,2)+X(I,3,3)
      
      RETURN
      END      
