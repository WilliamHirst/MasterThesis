C---------------------------------------------------------------------- 
      REAL*8 FUNCTION XKARE(X)            
C----------------------------------------------------------------------
C     Real part of function \kappa defined in Eq(B6) of Greub.
C
      IMPLICIT NONE
      REAL*8 X
      REAL*8 PI,GE
C
      PI=4.d0*ATAN(1.d0)
C
      IF(X.LT.4.d0) THEN
       GE=-2.d0*ATAN(SQRT(X/(4.d0-X)))**2
      ELSE
       GE=-PI**2/2.d0+2.d0*(LOG(SQRT(X)/2.d0+SQRT(X-4.d0)/2.d0))**2
      ENDIF
      XKARE=4.d0*(2.d0*GE+X)/X
C 
      RETURN 
      END
