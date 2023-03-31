C----------------------------------------------------------------------
      REAL*8 FUNCTION BI(X)            
C----------------------------------------------------------------------
C     Inami-Lim B function given in Table 1 of M.Ciuchini et al.
C
      IMPLICIT NONE
      REAL*8 X
C
      BI=1./4.d0*(X/(1.d0-X)+X*LOG(X)/(X-1.d0)**2)
C 
      RETURN 
      END
