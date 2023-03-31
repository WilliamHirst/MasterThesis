!
      SUBROUTINE UOUTCOUP(G215,Q)
!
!Purpose: To output the up-type Yukawas, in either
!         the mass or current basis
!
      IMPLICIT NONE
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
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
      DOUBLE COMPLEX GTMP(215),G215(215)
      DOUBLE PRECISION Q,COSB,SINB
      INTEGER I,SWROT
!
      SINB=DSQRT(TANB**2/(1+TANB**2))
      COSB=SINB/TANB
!
!This switch can be used to choose to write the rotated values.
!
      SWROT=1
!
!Remember the old gs
!
      DO I=1,215
        GTMP(I)=G215(I)
      END DO
!
!Perform the rotation to the mass basis if nec.
!
      IF(SWROT.EQ.1)CALL ROTBACK215(G215)
!
!Make some adjustments to facilitate comparisons
!
      DO I=1,215
        IF(ABS(G215(I)).LT.1.D-30)THEN
          G215(I)=(1.D-30,1.D-30)
        END IF
      END DO
!
!Print the set of terms out.
!
      IF(Q.GE.QNH-EPS)THEN
        WRITE(35,52)Q,ABS(G215(4)),ABS(G215(5)),ABS(G215(6)),
     $                ABS(G215(7)),ABS(G215(8)),ABS(G215(9)),
     $                ABS(G215(10)),ABS(G215(11)),ABS(G215(12))
      ELSE
        WRITE(35,52)Q,ABS(G215(34))/SINB,ABS(G215(35))/SINB,
     $   ABS(G215(36))/SINB,ABS(G215(37))/SINB,ABS(G215(38))/SINB,
     $   ABS(G215(39))/SINB,ABS(G215(40))/SINB,ABS(G215(41))/SINB,
     $   ABS(G215(42))/SINB
      END IF
!
!Return the gs to their original values
!
      DO I=1,215
        G215(I)=GTMP(I)
      END DO
!
   52 FORMAT(SP,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10
     $,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10)
!
      RETURN
      END
