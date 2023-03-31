C---------------------------------------------------------------------- 
      REAL*8 FUNCTION DI(X)            
C----------------------------------------------------------------------
C     Inami-Lim D function given in Table 1 of M.Ciuchini et al.
C
      IMPLICIT NONE
      REAL*8 X

      DI=-4./9.d0*LOG(X)+(-19.d0*X**3+25.d0*X**2)/36./(X-1.d0)**3
     $   +X**2*(5.d0*X**2-2.d0*X-6.d0)/18.d0/(X-1.d0)**4*LOG(X)
C 
      RETURN 
      END
