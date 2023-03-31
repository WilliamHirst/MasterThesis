!
      FUNCTION TCSFMUL(TH,A)
!
!Purpose: To trace the matrix A but taking into account
!         restrictions on the multiple from thetas.
!
!Expects: TH input as (3) integer (1 or 0)
!         A input as (3x3) complex
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX A(3,3),TCSFMUL
      INTEGER LOOP,TH(3),I,J
!
      TCSFMUL=(0.D0,0.D0)
      DO LOOP=1,3
        TCSFMUL=TCSFMUL+A(LOOP,LOOP)*DBLE(TH(LOOP))
      END DO
!
      RETURN
      END
