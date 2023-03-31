C---------------------------------------------------------------------- 
      REAL*8 FUNCTION GES(X)            
C----------------------------------------------------------------------
C    Phase space function given by Eq(5.10) of Greub.
C
      IMPLICIT NONE
      REAL*8 X
      
      GES=1.d0-8.d0*X**2+8.d0*X**6-X**8-24.d0*X**4*LOG(X)
      
      RETURN 
      END
