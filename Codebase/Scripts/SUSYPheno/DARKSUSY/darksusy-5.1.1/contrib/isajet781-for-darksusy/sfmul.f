!
      FUNCTION SFMUL(TH,A,B,I,J)
!
!Purpose: To multiply the two matrices A and B but taking
!         into account restrictions on the multiple from
!         thetas.
!
!Expects: TH input as (3) integer (1 or 0)
!         A, B input as (3x3) complex
!         I, J integer
!
      IMPLICIT NONE
!
      DOUBLE PRECISION A(3,3),B(3,3),SFMUL
      INTEGER LOOP,TH(3),I,J
!
      SFMUL=0.D0
      DO LOOP=1,3
        SFMUL=SFMUL+A(I,LOOP)*B(LOOP,J)*DBLE(TH(LOOP))
      END DO
!
      RETURN
      END
