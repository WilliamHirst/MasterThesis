!
      FUNCTION CRE(STAR,A,B)
!
!Purpose: To find the real part of A*B with:
!         If STAR=1 compute conjugate of A
!         If STAR=2 compute conjugate of B
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX A,B
      DOUBLE COMPLEX CCON,CRE
      INTEGER STAR
!
      CRE=0.D0
      IF(STAR.EQ.1)THEN
        CRE=DBLE(CCON(A)*B)
      ELSE IF(STAR.EQ.2)THEN
        CRE=DBLE(A*CCON(B))
      ELSE IF(STAR.EQ.0)THEN
        CRE=DBLE(A*B)
      ELSE
        WRITE(*,*)'WRONG STAR'
      END IF
!
      RETURN
      END
