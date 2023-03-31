C---------------------------------------------------------------------- 
      REAL*8 FUNCTION EI(X)            
C----------------------------------------------------------------------
C     Inami-Lim E function given in Table 1 of M.Ciuchini et al.
C
      IMPLICIT NONE
      REAL*8 X
C
      EI=-2./3.d0*LOG(X)
     &   +X**2*(15.d0-16.d0*X+4.*X**2)/6.d0/(X-1.d0)**4*LOG(X)
     $   +X*(18.d0-11.d0*X-X**2)/12.d0/(1.-X)**3
C 
      RETURN 
      END
