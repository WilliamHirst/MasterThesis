!
      FUNCTION DMATMUL(DAG,A,B,I,J)
!
!Purpose: To multiply the two matrices A and B
!         If DAG=1 compute dagger of A
!         If DAG=2 compute dagger of B
!
      IMPLICIT NONE
!
      DOUBLE PRECISION A(3,3),B(3,3)
      DOUBLE PRECISION DMATMUL,DAGA(3,3),DAGB(3,3)
      INTEGER LOOP,DAG,I,J
!
      DMATMUL=0.D0
      IF(DAG.EQ.1)THEN
        CALL DAGGER(A,DAGA)
        DO LOOP=1,3
          DMATMUL=DMATMUL+DAGA(I,LOOP)*B(LOOP,J)
        END DO
      ELSE IF(DAG.EQ.2)THEN
        CALL DAGGER(B,DAGB)
        DO LOOP=1,3
          DMATMUL=DMATMUL+A(I,LOOP)*DAGB(LOOP,J)
        END DO
      ELSE IF(DAG.EQ.0)THEN
        DO LOOP=1,3
          DMATMUL=DMATMUL+A(I,LOOP)*B(LOOP,J)
        END DO
      ELSE
        WRITE(*,*)'WRONG DAG IN DMATMUL',DAG
      END IF
!
      RETURN
      END
