C---------------------------------------------------------------------
      SUBROUTINE FUNS(X,FN1,FN2,FN3,FN4)            
C----------------------------------------------------------------------
C
C    Computes convenient functions given by Eq(75) of Anlauf.
C
      IMPLICIT NONE
      REAL*8 X
      REAL*8 FN1,FN2,FN3,FN4
C
      IF(ABS(X-1.d0).GT.1.d-2) THEN
       FN1=(X**2-5.*X-2.d0)/12.d0/(X-1.d0)**3
     $     +X*LOG(X)/2.d0/(X-1.d0)**4
       FN2=(2.*X**2+5.*X-1.d0)/12.d0/(X-1.d0)**3
     $    -X**2*LOG(X)/2.d0/(X-1.d0)**4
       FN3=(X-3.d0)/2.d0/(X-1.d0)**2
     $    +LOG(X)/(X-1.d0)**3
       FN4=(X+1.d0)/2.d0/(X-1.d0)**2
     $    -X*LOG(X)/(X-1.d0)**3           
      ELSE
       FN1=1./24.d0-(X-1.d0)/40.d0+(X-1.d0)**2/60.d0
       FN2=1./24.d0-(X-1.d0)/60.d0+(X-1.d0)**2/120.d0
       FN3=1./3.d0-(X-1.d0)/4.d0+(X-1.d0)**2/5.d0
       FN4=1./6.d0-(X-1.d0)/12.d0+(X-1.d0)**2/20.d0
      ENDIF
C
      RETURN 
      END
