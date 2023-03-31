C----------------------------------------------------------------------
      REAL*8 FUNCTION TR3X3(Y)
C----------------------------------------------------------------------
C
C     Takes the trace of 3x3 matrix Y.
C
      IMPLICIT NONE
      REAL*8 Y(3,3)
      
      TR3X3=Y(1,1)+Y(2,2)+Y(3,3)
      
      RETURN
      END      
