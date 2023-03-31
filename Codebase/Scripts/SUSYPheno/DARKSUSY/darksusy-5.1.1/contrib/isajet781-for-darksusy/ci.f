C---------------------------------------------------------------------- 
      REAL*8 FUNCTION CI(X)            
C----------------------------------------------------------------------
C     Inami-Lim C function given in Table 1 of M.Ciuchini et al.
C
      IMPLICIT NONE
      REAL*8 X

      CI=1./8.d0*X*((X-6.d0)/(X-1.d0)+(3.d0*X+2.d0)/(X-1.d0)**2*LOG(X))
C 
      RETURN 
      END
