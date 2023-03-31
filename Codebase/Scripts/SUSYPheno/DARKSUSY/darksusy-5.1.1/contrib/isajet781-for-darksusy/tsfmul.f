!
      FUNCTION TSFMUL(TH,A)
!
!Purpose: To trace the matrix A but taking into account
!         restrictions on the multiple from thetas.
!
!Expects: TH input as (3) integer (1 or 0)
!         A input as (3x3) complex
!
      IMPLICIT NONE
!
      DOUBLE PRECISION A(3,3),TSFMUL
      INTEGER LOOP,TH(3),I,J
!
      TSFMUL=0.D0
      DO LOOP=1,3
        TSFMUL=TSFMUL+A(LOOP,LOOP)*DBLE(TH(LOOP))
      END DO
!
      RETURN
      END
