!
      FUNCTION DRE(STAR,A,B)
!
!Purpose: To find the real part of A*B with:
!         If STAR=1 compute conjugate of A
!         If STAR=2 compute conjugate of B
!
      IMPLICIT NONE
!
      DOUBLE PRECISION A,B
      DOUBLE PRECISION CON
      DOUBLE PRECISION DRE
      INTEGER STAR
!
      DRE=0.D0
      IF(STAR.EQ.1)THEN
        DRE=DBLE(CON(A)*B)
      ELSE IF(STAR.EQ.2)THEN
        DRE=DBLE(A*CON(B))
      ELSE IF(STAR.EQ.0)THEN
        DRE=DBLE(A*B)
      ELSE
        WRITE(*,*)'WRONG STAR'
      END IF
!
      RETURN
      END
