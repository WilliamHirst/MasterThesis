CDECK  ID>, COLR12.
      SUBROUTINE COLR12(I1,L1,I2,L2,I3,L3,ICOLOR)
C
C     Set color flow lines for 1-> 2 decay
C     I1,I2,I3 = particle IDs
C     L1,L2,L3 = line number in PARTCL where they occur
      REAL X
      INTEGER I1,I2,I3,L1,L2,L3,ICOLOR(2,100)
      INTEGER J(3),IC(3),IC23
      J(1)=I1
      J(2)=I2
      J(3)=I3
C     Set QCD color labels
      DO I=1,3
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
      IC23=IC(2)*IC(3)
C     Do nothing for case of 1 -> 1 1 
C     1 -> 3 -3
C     Note: lines in ICOLOR offset from lines in PPTCL by 2
      L1=L1+2
      L2=L2+2
      L3=L3+2
      IF (IC(1).EQ.1.AND.IC(2).EQ.3) THEN
        ICOLOR(1,L2)=200+L1
        ICOLOR(2,L2)=0
        ICOLOR(1,L3)=0
        ICOLOR(2,L3)=200+L1
      ELSE IF (IC(1).EQ.1.AND.IC(2).EQ.-3) THEN
        ICOLOR(1,L2)=0
        ICOLOR(2,L2)=200+L1
        ICOLOR(1,L3)=200+L1
        ICOLOR(2,L3)=0
      END IF
C     3 -> 3 1
      IF (IC(1).EQ.3.AND.IC23.EQ.3) THEN
        IF (IC(2).EQ.3) THEN
          ICOLOR(1,L2)=ICOLOR(1,L1)
          ICOLOR(2,L2)=0
          ICOLOR(1,L3)=0
          ICOLOR(2,L3)=0
        ELSE IF (IC(2).EQ.1) THEN
          ICOLOR(1,L2)=0
          ICOLOR(2,L2)=0
          ICOLOR(1,L3)=ICOLOR(1,L1)
          ICOLOR(2,L3)=0
        END IF
      END IF
C     3* -> 3* 1
      IF (IC(1).EQ.-3.AND.IC23.EQ.-3) THEN
        IF (IC(2).EQ.-3) THEN
          ICOLOR(1,L2)=0
          ICOLOR(2,L2)=ICOLOR(2,L1)
          ICOLOR(1,L3)=0
          ICOLOR(2,L3)=0
        ELSE IF (IC(2).EQ.1) THEN
          ICOLOR(1,L2)=0
          ICOLOR(2,L2)=0
          ICOLOR(1,L3)=0
          ICOLOR(2,L3)=ICOLOR(2,L1)
        END IF
      END IF
C     3 -> 3 8
      IF (IC(1).EQ.3.AND.IC23.EQ.24) THEN
        IF (IC(2).EQ.3) THEN
            ICOLOR(1,L2)=200+L1
            ICOLOR(2,L2)=0
            ICOLOR(1,L3)=ICOLOR(1,L1)
            ICOLOR(2,L3)=200+L1
        ELSE IF (IC(2).EQ.8) THEN
            ICOLOR(1,L2)=ICOLOR(1,L1)
            ICOLOR(2,L2)=200+L1
            ICOLOR(1,L3)=200+L1
            ICOLOR(2,L3)=0
        END IF 
      END IF
C     3* -> 3* 8
      IF (IC(1).EQ.-3.AND.IC23.EQ.-24) THEN
        IF (IC(2).EQ.-3) THEN
          ICOLOR(1,L2)=0
          ICOLOR(2,L2)=200+L1
          ICOLOR(1,L3)=200+L1
          ICOLOR(2,L3)=ICOLOR(2,L1)
        ELSE IF (IC(2).EQ.8) THEN
          ICOLOR(1,L2)=200+L1
          ICOLOR(2,L2)=ICOLOR(2,L1)
          ICOLOR(1,L3)=0
          ICOLOR(2,L3)=200+L1
        END IF 
      END IF
C     8 -> 3 3*
      IF (IC(1).EQ.8.AND.IC(2).EQ.3) THEN
        ICOLOR(1,L2)=ICOLOR(1,L1)
        ICOLOR(2,L2)=0
        ICOLOR(1,L3)=0
        ICOLOR(2,L3)=ICOLOR(2,L1)
      ELSE IF (IC(1).EQ.8.AND.IC(2).EQ.-3) THEN
        ICOLOR(1,L2)=O
        ICOLOR(2,L2)=ICOLOR(2,L1)
        ICOLOR(1,L3)=ICOLOR(1,L1)
        ICOLOR(2,L3)=0
      END IF
C     8 -> 8 1
      IF (IC(1).EQ.8.AND.IC(2).EQ.8) THEN
        ICOLOR(1,L2)=ICOLOR(1,L1)
        ICOLOR(2,L2)=ICOLOR(2,L1)
        ICOLOR(1,L3)=0
        ICOLOR(2,L3)=0
      ELSE IF (IC(1).EQ.8.AND.IC(2).EQ.1) THEN
        ICOLOR(1,L2)=O
        ICOLOR(2,L2)=0
        ICOLOR(1,L3)=ICOLOR(1,L1)
        ICOLOR(2,L3)=ICOLOR(2,L1)
      END IF
      RETURN
      END
