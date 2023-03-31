C---------------------------------------------------------------------- 
      REAL*8 FUNCTION FES(X)            
C----------------------------------------------------------------------
C     Pole to MSbar conversion factor given by Eq(5.7) of Greub.
C
      IMPLICIT NONE
      REAL*8 X,PI
C
      PI=4.d0*DATAN(1.d0)       
      FES=(PI**2-31./4.d0)*(1.d0-X)**2+3./2.d0
C 
      RETURN 
      END
