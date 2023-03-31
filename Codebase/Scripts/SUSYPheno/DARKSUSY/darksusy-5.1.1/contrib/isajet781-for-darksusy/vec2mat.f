C----------------------------------------------------------------------
      SUBROUTINE VEC2MAT(G,INIT,MAT,IORD)
C----------------------------------------------------------------------
C
C    Assembles 3x3 matrix from vector G(i) starting from element i=INIT.
C    Matrix can be filled in row- or column-vise order.
C
C    Created: 5/03/07 by Azar Mustafayev.
C
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      IMPLICIT NONE
      REAL*8 G(157),MAT(3,3)
      INTEGER INIT,IORD
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/

      IF (IORD.EQ.1) THEN        ! row-vise
	MAT(1,1)=G(INIT)
	MAT(1,2)=G(INIT+1)
	MAT(1,3)=G(INIT+2)
	MAT(2,1)=G(INIT+3)
	MAT(2,2)=G(INIT+4)
	MAT(2,3)=G(INIT+5)
	MAT(3,1)=G(INIT+6)
	MAT(3,2)=G(INIT+7)
	MAT(3,3)=G(INIT+8)
      ELSEIF (IORD.EQ.-1) THEN    ! column-vise
	MAT(1,1)=G(INIT)
	MAT(2,1)=G(INIT+1)
	MAT(3,1)=G(INIT+2)
	MAT(1,2)=G(INIT+3)
	MAT(2,2)=G(INIT+4)
	MAT(3,2)=G(INIT+5)
	MAT(1,3)=G(INIT+6)
	MAT(2,3)=G(INIT+7)
	MAT(3,3)=G(INIT+8)
      ELSE
        WRITE(LOUT,*) 'VEC2MAT: unknown assembly order'
        STOP99
      ENDIF  

      RETURN
      END
