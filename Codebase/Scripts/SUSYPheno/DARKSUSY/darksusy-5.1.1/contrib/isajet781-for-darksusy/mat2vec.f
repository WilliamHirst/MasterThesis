C----------------------------------------------------------------------
      SUBROUTINE MAT2VEC(G,INIT,MAT,IORD)
C----------------------------------------------------------------------
C
C   This is reciprocal to VEC2MAT.
C  "Disassembles" 3x3 matrix into vector G(i) starting from element i=INIT.
C   Matrix can be broken in row- or column-vise order.
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
	G(INIT)  =MAT(1,1)
	G(INIT+1)=MAT(1,2)
	G(INIT+2)=MAT(1,3)
	G(INIT+3)=MAT(2,1)
	G(INIT+4)=MAT(2,2)
	G(INIT+5)=MAT(2,3)
	G(INIT+6)=MAT(3,1)
	G(INIT+7)=MAT(3,2)
	G(INIT+8)=MAT(3,3)
      ELSEIF (IORD.EQ.-1) THEN    ! column-vise
	G(INIT)  =MAT(1,1)
	G(INIT+1)=MAT(2,1)
	G(INIT+2)=MAT(3,1)
	G(INIT+3)=MAT(1,2)
	G(INIT+4)=MAT(2,2)
	G(INIT+5)=MAT(3,2)
	G(INIT+6)=MAT(1,3)
	G(INIT+7)=MAT(2,3)
	G(INIT+8)=MAT(3,3)
      ELSE
        WRITE(LOUT,*) 'MAT2VEC: unknown assembly order'
        STOP99
      ENDIF  

      RETURN
      END
