C---------------------------------------------------------------------- 
      REAL*8 FUNCTION XKAIM(X)            
C----------------------------------------------------------------------
C     Imaginary part of function \kappa defined in Eq(B6) of Greub.
C
      IMPLICIT NONE
      REAL*8 X
      REAL*8 PI,GE
C
      PI=4.d0*ATAN(1.d0)
C
      IF(X.LT.4.d0) THEN
       GE=0.d0
      ELSE
       GE=-2.d0*PI*LOG(SQRT(X)/2.d0+SQRT(X-4.d0)/2.d0)
      ENDIF
C 
      XKAIM=8.d0*GE/X
C
      RETURN 
      END
