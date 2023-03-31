!
      SUBROUTINE SORTTH
!
!Purpose: Sorts the thresholds into ascending order
!
      IMPLICIT NONE
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      DOUBLE PRECISION TEMP
      INTEGER I,J,K
!
      DO I=1,3
        QTHSORT(I)=QTHQL(I)
        QTHSORT(I+3)=QTHUR(I)
        QTHSORT(I+6)=QTHDR(I)
        QTHSORT(I+9)=QTHLL(I)
        QTHSORT(I+12)=QTHER(I)
      END DO
      QTHSORT(16)=QNSH
      QTHSORT(17)=QNSG
      QTHSORT(18)=QNH
      QTHSORT(19)=QTHSB
      QTHSORT(20)=QTHSW
      LOCMH=18
      DO J=2,20
        TEMP=QTHSORT(J)
        DO I=J-1,1,-1
          IF(QTHSORT(I).LT.TEMP)GOTO 10
          QTHSORT(I+1)=QTHSORT(I)
!
          IF(I.EQ.LOCMH)LOCMH=I+1
!
        END DO
        I=0
  10    QTHSORT(I+1)=TEMP
!
        IF(J.EQ.18)LOCMH=I+1
!
      END DO
!
      RETURN
      END
