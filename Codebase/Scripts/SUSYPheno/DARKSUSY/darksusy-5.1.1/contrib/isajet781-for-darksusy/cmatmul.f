!
      FUNCTION CMATMUL(DAG,A,B,I,J)
!
!Purpose: To multiply the two matrices A and B
!         If DAG=1 compute dagger of A
!         If DAG=2 compute dagger of B
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX A(3,3),B(3,3)
      DOUBLE COMPLEX CMATMUL,DAGA(3,3),DAGB(3,3)
      INTEGER LOOP,DAG,I,J
!
      CMATMUL=(0.D0,0.D0)
      IF(DAG.EQ.1)THEN
        CALL CDAGGER(A,DAGA)
        DO LOOP=1,3
          CMATMUL=CMATMUL+DAGA(I,LOOP)*B(LOOP,J)
        END DO
      ELSE IF(DAG.EQ.2)THEN
        CALL CDAGGER(B,DAGB)
        DO LOOP=1,3
          CMATMUL=CMATMUL+A(I,LOOP)*DAGB(LOOP,J)
        END DO
      ELSE IF(DAG.EQ.0)THEN
        DO LOOP=1,3
          CMATMUL=CMATMUL+A(I,LOOP)*B(LOOP,J)
        END DO
      ELSE
        WRITE(*,*)'WRONG DAG IN MATMUL'
      END IF
!
      RETURN
      END
