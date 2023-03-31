!
      SUBROUTINE DECOUTCOUP(Q,G)
!
!Purpose: To output m^2_U in the mass basis as it runs down.
!
      IMPLICIT NONE
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      DOUBLE COMPLEX GTMP(601),G(601)
      DOUBLE PRECISION Q,COSB,SINB
      INTEGER I,SWABS,SWROT
!
      SINB=DSQRT(TANB**2/(1+TANB**2))
      COSB=SINB/TANB
!
!These switches can be used to print all complex entries or
!turn off the rotation.
!
      SWABS=1
      SWROT=1
!
!Remember the old gs
!
      DO I=1,601
        GTMP(I)=G(I)
      END DO
!
!Perform the rotation to the mass basis if nec.
!
      IF(SWROT.EQ.1)CALL STROTBACK(GTMP,G,0)
!
!Make some adjustments to facilitate comparisons
!
      IF(SWABS.EQ.1)THEN
        DO I=1,601
          IF(ABS(G(I)).LT.1.D-30)THEN
            G(I)=(1.D-30,1.D-30)
          END IF
        END DO
      END IF
!
!Print out the terms - absolute squares first
!
      IF(SWABS.EQ.1)THEN
        I=81
        WRITE(63,52)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
!
!Now twice the size with real and imaginary parts
!
      ELSE
        I=81
        WRITE(77,54)Q,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5)
     $    ,G(I+6),G(I+7),G(I+8)
      END IF
!
!Return the gs to their original values
!
      DO I=1,601
        G(I)=GTMP(I)
      END DO
!
   52 FORMAT(SP,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10
     $,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10)
   54 FORMAT(SP,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10
     $,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X
     $,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X
     $,D17.10)
!
      RETURN
      END
