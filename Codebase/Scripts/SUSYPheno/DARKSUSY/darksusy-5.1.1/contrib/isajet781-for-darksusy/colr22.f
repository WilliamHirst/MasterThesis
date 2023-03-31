CDECK  ID>, COLR22.
      SUBROUTINE COLR22(I1,I2,I3,I4,ICOLOR)
C
C     Set color flow lines for 2-> 2 subprocesses
C
      REAL X
      INTEGER I1,I2,I3,I4,IC12,IC34,ICOLOR(2,100)
      INTEGER J(4),IC(4)
      DO I=1,100
       ICOLOR(1,I)=0
       ICOLOR(2,I)=0
      END DO
      J(1)=I1
      J(2)=I2
      J(3)=I3
      J(4)=I4
      DO I=1,4
       IC(I)=1
       IF (J(I).EQ.9.OR.J(I).EQ.29) THEN
         IC(I)=8
       END IF
       IF ((J(I).GE.1.AND.J(I).LE.6).OR.(J(I).GE.21.AND.J(I).LE.26).OR.
     $(J(I).GE.41.AND.J(I).LE.46)) THEN
         IC(I)=3
       END IF
       IF ((-J(I).GE.1.AND.-J(I).LE.6).OR.(-J(I).GE.21.AND.-J(I).LE.26)
     $.OR.(-J(I).GE.41.AND.-J(I).LE.46)) THEN
         IC(I)=-3
       END IF
      END DO
C     Do nothing for case of 1 1 -> 1 1 
C     For now, Select random number to determine color flow
      X=RANF()
      IC12=IC(1)*IC(2)
      IC34=IC(3)*IC(4)
C     1 1 -> 3 3*
      IF (IC12.EQ.1.AND.IC34.EQ.-9) THEN
        IF (IC(3).EQ.3) THEN
          ICOLOR(1,3)=101
          ICOLOR(2,3)=0
          ICOLOR(1,4)=0
          ICOLOR(2,4)=101
        ELSE IF (IC(3).EQ.-3) THEN
          ICOLOR(1,3)=0
          ICOLOR(2,3)=101
          ICOLOR(1,4)=101
          ICOLOR(2,4)=0
        END IF
      END IF
C     3 3 -> 3 3
      IF (IC12.EQ.9.AND.IC34.EQ.9.AND.IC(3).EQ.3) THEN
        IF (X.LT..5) THEN
          ICOLOR(1,1)=101
          ICOLOR(2,1)=0
          ICOLOR(1,2)=102
          ICOLOR(2,2)=0
          ICOLOR(1,3)=101
          ICOLOR(2,3)=0
          ICOLOR(1,4)=102
          ICOLOR(2,4)=0
        ELSE
          ICOLOR(1,1)=101
          ICOLOR(2,1)=0
          ICOLOR(1,2)=102
          ICOLOR(2,2)=0
          ICOLOR(1,3)=102
          ICOLOR(2,3)=0
          ICOLOR(1,4)=101
          ICOLOR(2,4)=0
        END IF
      END IF 
C     3* 3* -> 3* 3*
      IF (IC12.EQ.9.AND.IC34.EQ.9.AND.IC(3).EQ.-3) THEN
        IF (X.LT..5) THEN
          ICOLOR(1,1)=0
          ICOLOR(2,1)=101
          ICOLOR(1,2)=0
          ICOLOR(2,2)=102
          ICOLOR(1,3)=0
          ICOLOR(2,3)=101
          ICOLOR(1,4)=0
          ICOLOR(2,4)=102
        ELSE
          ICOLOR(1,1)=0
          ICOLOR(2,1)=101
          ICOLOR(1,2)=0
          ICOLOR(2,2)=102
          ICOLOR(1,3)=0
          ICOLOR(2,3)=102
          ICOLOR(1,4)=0
          ICOLOR(2,4)=101
        END IF
      END IF 
C     3 3* -> 3 3*
      IF (IC12.EQ.-9.AND.IC34.EQ.-9) THEN
        IF (IC(1).EQ.3.AND.IC(3).EQ.3) THEN
          IF (X.LT..5) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=0
            ICOLOR(1,2)=0
            ICOLOR(2,2)=101
            ICOLOR(1,3)=102
            ICOLOR(2,3)=0
            ICOLOR(1,4)=0
            ICOLOR(2,4)=102
          ELSE
            ICOLOR(1,1)=101
            ICOLOR(2,1)=0
            ICOLOR(1,2)=0
            ICOLOR(2,2)=102
            ICOLOR(1,3)=101
            ICOLOR(2,3)=0
            ICOLOR(1,4)=0
            ICOLOR(2,4)=102
          END IF
        END IF
        IF (IC(1).EQ.3.AND.IC(4).EQ.3) THEN
          IF (X.LT..5) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=0
            ICOLOR(1,2)=0
            ICOLOR(2,2)=101
            ICOLOR(1,3)=0
            ICOLOR(2,3)=102
            ICOLOR(1,4)=102
            ICOLOR(2,4)=0
          ELSE
            ICOLOR(1,1)=101
            ICOLOR(2,1)=0
            ICOLOR(1,2)=0
            ICOLOR(2,2)=102
            ICOLOR(1,3)=0
            ICOLOR(2,3)=102
            ICOLOR(1,4)=101
            ICOLOR(2,4)=0
          END IF
        END IF
        IF (IC(2).EQ.3.AND.IC(3).EQ.3) THEN
          IF (X.LT..5) THEN
            ICOLOR(1,1)=0
            ICOLOR(2,1)=101
            ICOLOR(1,2)=101
            ICOLOR(2,2)=0
            ICOLOR(1,3)=102
            ICOLOR(2,3)=0
            ICOLOR(1,4)=0
            ICOLOR(2,4)=102
          ELSE
            ICOLOR(1,1)=0
            ICOLOR(2,1)=101
            ICOLOR(1,2)=102
            ICOLOR(2,2)=0
            ICOLOR(1,3)=102
            ICOLOR(2,3)=0
            ICOLOR(1,4)=0
            ICOLOR(2,4)=101
          END IF
        END IF
        IF (IC(2).EQ.3.AND.IC(4).EQ.3) THEN
          IF (X.LT..5) THEN
            ICOLOR(1,1)=0
            ICOLOR(2,1)=101
            ICOLOR(1,2)=101
            ICOLOR(2,2)=0
            ICOLOR(1,3)=0
            ICOLOR(2,3)=102
            ICOLOR(1,4)=102
            ICOLOR(2,4)=0
          ELSE
            ICOLOR(1,1)=0
            ICOLOR(2,1)=101
            ICOLOR(1,2)=102
            ICOLOR(2,2)=0
            ICOLOR(1,3)=0
            ICOLOR(2,3)=101
            ICOLOR(1,4)=102
            ICOLOR(2,4)=0
          END IF
        END IF 
      END IF
C     3 3* -> 1 1
      IF (IC12.EQ.-9.AND.IC34.EQ.1) THEN
        IF (IC(1).EQ.3) THEN
          ICOLOR(1,1)=101
          ICOLOR(2,1)=0
          ICOLOR(1,2)=0
          ICOLOR(2,2)=101
        ELSE IF (IC(1).EQ.-3) THEN
          ICOLOR(1,1)=0
          ICOLOR(2,1)=101
          ICOLOR(1,2)=101
          ICOLOR(2,2)=0
        END IF        
      END IF
C     3 3* -> 8 8
      IF (IC12.EQ.-9.AND.IC34.EQ.64) THEN
        IF (IC(1).EQ.3) THEN
          IF (X.LT..5) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=0
            ICOLOR(1,2)=0
            ICOLOR(2,2)=102
            ICOLOR(1,3)=101
            ICOLOR(2,3)=103
            ICOLOR(1,4)=103
            ICOLOR(2,4)=102
          ELSE
            ICOLOR(1,1)=101
            ICOLOR(2,1)=0
            ICOLOR(1,2)=0
            ICOLOR(2,2)=102
            ICOLOR(1,3)=103
            ICOLOR(2,3)=102
            ICOLOR(1,4)=101
            ICOLOR(2,4)=103
          END IF
        END IF 
        IF (IC(2).EQ.3) THEN
          IF (X.LT..5) THEN
            ICOLOR(1,1)=0
            ICOLOR(2,1)=101
            ICOLOR(1,2)=102
            ICOLOR(2,2)=0
            ICOLOR(1,3)=103
            ICOLOR(2,3)=101
            ICOLOR(1,4)=102
            ICOLOR(2,4)=103
          ELSE
            ICOLOR(1,1)=0
            ICOLOR(2,1)=101
            ICOLOR(1,2)=102
            ICOLOR(2,2)=0
            ICOLOR(1,3)=102
            ICOLOR(2,3)=103
            ICOLOR(1,4)=103
            ICOLOR(2,4)=101
          END IF
        END IF 
      END IF
C     3 3* -> 1 8
      IF (IC12.EQ.-9.AND.IC34.EQ.8) THEN
        IF (IC(1).EQ.3.AND.IC(3).EQ.1) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=0
            ICOLOR(1,2)=0
            ICOLOR(2,2)=102
            ICOLOR(1,3)=0
            ICOLOR(2,3)=0
            ICOLOR(1,4)=101
            ICOLOR(2,4)=102
        ELSE IF (IC(1).EQ.3.AND.IC(4).EQ.1) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=0
            ICOLOR(1,2)=0
            ICOLOR(2,2)=102
            ICOLOR(1,3)=101
            ICOLOR(2,3)=102
            ICOLOR(1,4)=0
            ICOLOR(2,4)=0
        ELSE IF (IC(1).EQ.-3.AND.IC(3).EQ.1) THEN
            ICOLOR(1,1)=0
            ICOLOR(2,1)=101
            ICOLOR(1,2)=102
            ICOLOR(2,2)=0
            ICOLOR(1,3)=0
            ICOLOR(2,3)=0
            ICOLOR(1,4)=102
            ICOLOR(2,4)=101
        ELSE IF (IC(1).EQ.-3.AND.IC(4).EQ.1) THEN
            ICOLOR(1,1)=0
            ICOLOR(2,1)=101
            ICOLOR(1,2)=102
            ICOLOR(2,2)=0
            ICOLOR(1,3)=102
            ICOLOR(2,3)=101
            ICOLOR(1,4)=0
            ICOLOR(2,4)=0
        END IF 
      END IF
C     3 8 -> 1 3
      IF (IC12.EQ.24.AND.IC34.EQ.3) THEN
        IF (IC(1).EQ.3.AND.IC(3).EQ.1) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=0
            ICOLOR(1,2)=102
            ICOLOR(2,2)=101
            ICOLOR(1,3)=0
            ICOLOR(2,3)=0
            ICOLOR(1,4)=102
            ICOLOR(2,4)=0
        ELSE IF (IC(1).EQ.3.AND.IC(4).EQ.1) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=0
            ICOLOR(1,2)=102
            ICOLOR(2,2)=101
            ICOLOR(1,3)=102
            ICOLOR(2,3)=0
            ICOLOR(1,4)=0
            ICOLOR(2,4)=0
        ELSE IF (IC(2).EQ.3.AND.IC(3).EQ.1) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=102
            ICOLOR(1,2)=102
            ICOLOR(2,2)=0
            ICOLOR(1,3)=0
            ICOLOR(2,3)=0
            ICOLOR(1,4)=101
            ICOLOR(2,4)=0
        ELSE IF (IC(2).EQ.3.AND.IC(4).EQ.1) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=102
            ICOLOR(1,2)=102
            ICOLOR(2,2)=0
            ICOLOR(1,3)=101
            ICOLOR(2,3)=0
            ICOLOR(1,4)=0
            ICOLOR(2,4)=0
        END IF 
      END IF
C     3* 8 -> 1 3*
      IF (IC12.EQ.-24.AND.IC34.EQ.-3) THEN
        IF (IC(1).EQ.-3.AND.IC(3).EQ.1) THEN
            ICOLOR(1,1)=0
            ICOLOR(2,1)=101
            ICOLOR(1,2)=101
            ICOLOR(2,2)=102
            ICOLOR(1,3)=0
            ICOLOR(2,3)=0
            ICOLOR(1,4)=0
            ICOLOR(2,4)=102
        ELSE IF (IC(1).EQ.-3.AND.IC(4).EQ.1) THEN
            ICOLOR(1,1)=0
            ICOLOR(2,1)=101
            ICOLOR(1,2)=101
            ICOLOR(2,2)=102
            ICOLOR(1,3)=0
            ICOLOR(2,3)=102
            ICOLOR(1,4)=0
            ICOLOR(2,4)=0
        ELSE IF (IC(2).EQ.-3.AND.IC(3).EQ.1) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=102
            ICOLOR(1,2)=0
            ICOLOR(2,2)=101
            ICOLOR(1,3)=0
            ICOLOR(2,3)=0
            ICOLOR(1,4)=0
            ICOLOR(2,4)=102
        ELSE IF (IC(2).EQ.-3.AND.IC(4).EQ.1) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=102
            ICOLOR(1,2)=0
            ICOLOR(2,2)=101
            ICOLOR(1,3)=0
            ICOLOR(2,3)=102
            ICOLOR(1,4)=0
            ICOLOR(2,4)=0
        END IF 
      END IF
C     3 8 -> 3 8
      IF (IC12.EQ.24.AND.IC34.EQ.24) THEN
        IF (IC(1).EQ.3.AND.IC(3).EQ.3) THEN
          IF (X.LT..5) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=0
            ICOLOR(1,2)=102
            ICOLOR(2,2)=101
            ICOLOR(1,3)=103
            ICOLOR(2,3)=0
            ICOLOR(1,4)=102
            ICOLOR(2,4)=103
          ELSE
            ICOLOR(1,1)=101
            ICOLOR(2,1)=0
            ICOLOR(1,2)=102
            ICOLOR(2,2)=103
            ICOLOR(1,3)=102
            ICOLOR(2,3)=0
            ICOLOR(1,4)=101
            ICOLOR(2,4)=103
          END IF
        END IF 
        IF (IC(1).EQ.3.AND.IC(4).EQ.3) THEN
          IF (X.LT..5) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=0
            ICOLOR(1,2)=102
            ICOLOR(2,2)=101
            ICOLOR(1,3)=102
            ICOLOR(2,3)=103
            ICOLOR(1,4)=103
            ICOLOR(2,4)=0
          ELSE
            ICOLOR(1,1)=101
            ICOLOR(2,1)=0
            ICOLOR(1,2)=102
            ICOLOR(2,2)=103
            ICOLOR(1,3)=101
            ICOLOR(2,3)=103
            ICOLOR(1,4)=102
            ICOLOR(2,4)=0
          END IF
        END IF 
        IF (IC(2).EQ.3.AND.IC(3).EQ.3) THEN
          IF (X.LT..5) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=102
            ICOLOR(1,2)=102
            ICOLOR(2,2)=0
            ICOLOR(1,3)=103
            ICOLOR(2,3)=0
            ICOLOR(1,4)=101
            ICOLOR(2,4)=103
          ELSE
            ICOLOR(1,1)=102
            ICOLOR(2,1)=103
            ICOLOR(1,2)=101
            ICOLOR(2,2)=0
            ICOLOR(1,3)=102
            ICOLOR(2,3)=0
            ICOLOR(1,4)=101
            ICOLOR(2,4)=103
          END IF
        END IF 
        IF (IC(2).EQ.3.AND.IC(4).EQ.3) THEN
          IF (X.LT..5) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=102
            ICOLOR(1,2)=102
            ICOLOR(2,2)=0
            ICOLOR(1,3)=101
            ICOLOR(2,3)=103
            ICOLOR(1,4)=103
            ICOLOR(2,4)=0
          ELSE
            ICOLOR(1,1)=102
            ICOLOR(2,1)=103
            ICOLOR(1,2)=101
            ICOLOR(2,2)=0
            ICOLOR(1,3)=101
            ICOLOR(2,3)=103
            ICOLOR(1,4)=102
            ICOLOR(2,4)=0
          END IF
        END IF 
      END IF
C     3* 8 -> 3* 8
      IF (IC12.EQ.-24.AND.IC34.EQ.-24) THEN
        IF (IC(1).EQ.-3.AND.IC(3).EQ.-3) THEN
          IF (X.LT..5) THEN
            ICOLOR(1,1)=0
            ICOLOR(2,1)=101
            ICOLOR(1,2)=101
            ICOLOR(2,2)=102
            ICOLOR(1,3)=0
            ICOLOR(2,3)=103
            ICOLOR(1,4)=103
            ICOLOR(2,4)=102
          ELSE
            ICOLOR(1,1)=0
            ICOLOR(2,1)=101
            ICOLOR(1,2)=102
            ICOLOR(2,2)=103
            ICOLOR(1,3)=0
            ICOLOR(2,3)=103
            ICOLOR(1,4)=102
            ICOLOR(2,4)=101
          END IF
        END IF 
        IF (IC(1).EQ.-3.AND.IC(4).EQ.-3) THEN
          IF (X.LT..5) THEN
            ICOLOR(1,1)=0
            ICOLOR(2,1)=101
            ICOLOR(1,2)=101
            ICOLOR(2,2)=102
            ICOLOR(1,3)=103
            ICOLOR(2,3)=102
            ICOLOR(1,4)=0
            ICOLOR(2,4)=103
          ELSE
            ICOLOR(1,1)=0
            ICOLOR(2,1)=101
            ICOLOR(1,2)=102
            ICOLOR(2,2)=103
            ICOLOR(1,3)=102
            ICOLOR(2,3)=101
            ICOLOR(1,4)=0
            ICOLOR(2,4)=103
          END IF
        END IF 
        IF (IC(2).EQ.-3.AND.IC(3).EQ.-3) THEN
          IF (X.LT..5) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=102
            ICOLOR(1,2)=0
            ICOLOR(2,2)=101
            ICOLOR(1,3)=0
            ICOLOR(2,3)=103
            ICOLOR(1,4)=103
            ICOLOR(2,4)=102
          ELSE
            ICOLOR(1,1)=101
            ICOLOR(2,1)=102
            ICOLOR(1,2)=0
            ICOLOR(2,2)=103
            ICOLOR(1,3)=0
            ICOLOR(2,3)=102
            ICOLOR(1,4)=101
            ICOLOR(2,4)=103
          END IF
        END IF 
        IF (IC(2).EQ.-3.AND.IC(4).EQ.-3) THEN
          IF (X.LT..5) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=102
            ICOLOR(1,2)=0
            ICOLOR(2,2)=101
            ICOLOR(1,3)=103
            ICOLOR(2,3)=102
            ICOLOR(1,4)=0
            ICOLOR(2,4)=103
          ELSE
            ICOLOR(1,1)=101
            ICOLOR(2,1)=102
            ICOLOR(1,2)=0
            ICOLOR(2,2)=103
            ICOLOR(1,3)=101
            ICOLOR(2,3)=103
            ICOLOR(1,4)=0
            ICOLOR(2,4)=102
          END IF
        END IF 
      END IF
C     8 8 -> 3 3*
      IF (IC12.EQ.64.AND.IC34.EQ.-9) THEN
        IF (IC(3).EQ.3) THEN
          IF (X.LT..5) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=102
            ICOLOR(1,2)=102
            ICOLOR(2,2)=103
            ICOLOR(1,3)=101
            ICOLOR(2,3)=0
            ICOLOR(1,4)=0
            ICOLOR(2,4)=103
          ELSE
            ICOLOR(1,1)=101
            ICOLOR(2,1)=102
            ICOLOR(1,2)=103
            ICOLOR(2,2)=101
            ICOLOR(1,3)=103
            ICOLOR(2,3)=0
            ICOLOR(1,4)=0
            ICOLOR(2,4)=102
          END IF
        END IF 
        IF (IC(4).EQ.3) THEN
          IF (X.LT..5) THEN
            ICOLOR(1,1)=101
            ICOLOR(2,1)=102
            ICOLOR(1,2)=102
            ICOLOR(2,2)=103
            ICOLOR(1,3)=0
            ICOLOR(2,3)=103
            ICOLOR(1,4)=101
            ICOLOR(2,4)=0
          ELSE
            ICOLOR(1,1)=101
            ICOLOR(2,1)=102
            ICOLOR(1,2)=103
            ICOLOR(2,2)=101
            ICOLOR(1,3)=0
            ICOLOR(2,3)=102
            ICOLOR(1,4)=103
            ICOLOR(2,4)=0
          END IF
        END IF 
      END IF
C     8 8 -> 8 8
      IF (IC12.EQ.64.AND.IC34.EQ.64) THEN
        IF (X.LT..167) THEN
          ICOLOR(1,1)=101
          ICOLOR(2,1)=102
          ICOLOR(1,2)=103
          ICOLOR(2,2)=101
          ICOLOR(1,3)=104
          ICOLOR(2,3)=102
          ICOLOR(1,4)=103
          ICOLOR(2,4)=104
        ELSE IF (X.GE..167.AND.X.LT..334) THEN
          ICOLOR(1,1)=101
          ICOLOR(2,1)=102
          ICOLOR(1,2)=102
          ICOLOR(2,2)=103
          ICOLOR(1,3)=101
          ICOLOR(2,3)=104
          ICOLOR(1,4)=104
          ICOLOR(2,4)=103
        ELSE IF (X.GE..334.AND.X.LT..501) THEN
          ICOLOR(1,1)=101
          ICOLOR(2,1)=102
          ICOLOR(1,2)=103
          ICOLOR(2,2)=104
          ICOLOR(1,3)=103
          ICOLOR(2,3)=102
          ICOLOR(1,4)=101
          ICOLOR(2,4)=104
        ELSE IF (X.GE..501.AND.X.LT..668) THEN
          ICOLOR(1,1)=101
          ICOLOR(2,1)=102
          ICOLOR(1,2)=103
          ICOLOR(2,2)=101
          ICOLOR(1,3)=103
          ICOLOR(2,3)=104
          ICOLOR(1,4)=104
          ICOLOR(2,4)=102
        ELSE IF (X.GE..668.AND.X.LT..825) THEN
          ICOLOR(1,1)=101
          ICOLOR(2,1)=102
          ICOLOR(1,2)=104
          ICOLOR(2,2)=103
          ICOLOR(1,3)=101
          ICOLOR(2,3)=103
          ICOLOR(1,4)=104
          ICOLOR(2,4)=102
        ELSE IF (X.GE..825.AND.X.LE.1.) THEN
          ICOLOR(1,1)=101
          ICOLOR(2,1)=102
          ICOLOR(1,2)=102
          ICOLOR(2,2)=103
          ICOLOR(1,3)=104
          ICOLOR(2,3)=103
          ICOLOR(1,4)=101
          ICOLOR(2,4)=104
        END IF
      END IF
      RETURN
      END
