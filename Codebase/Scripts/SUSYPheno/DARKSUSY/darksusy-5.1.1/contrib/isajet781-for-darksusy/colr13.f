CDECK  ID>, COLR13.
      SUBROUTINE COLR13(I1,L1,I2,L2,I3,L3,I4,L4,ICOLOR)
C
C     Set color flow lines for 1-> 3 decay
C     I1,I2,I3,I4 = particle IDs
C     L1,L2,L3,L4 = line number in PARTCL where they occur
C     Isajet convention is that colored particles occur first
C       in any decay string
      INTEGER I1,I2,I3,I4,L1,L2,L3,L4,ICOLOR(2,100)
      INTEGER J(4),IC(4),IC12,IC34,IC23
      J(1)=I1
      J(2)=I2
      J(3)=I3
      J(4)=I4
C     Set QCD color labels
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
C     Do nothing for case of 1 -> 1 1 1 
C     1 -> 3 3* 1
C     Note lines in ICOLOR offset from PPTCL lines by 2
      L1=L1+2
      L2=L2+2
      L3=L3+2
      L4=L4+2
      IC12=IC(1)*IC(2)
      IC34=IC(3)*IC(4)
      IC23=IC(2)*IC(3)
      IF (IC(1).EQ.1.AND.IC(2).EQ.3) THEN
          ICOLOR(1,L2)=300+L1
          ICOLOR(2,L2)=0
          ICOLOR(1,L3)=0
          ICOLOR(2,L3)=300+L1
          ICOLOR(1,L4)=0
          ICOLOR(2,L4)=0
      ELSE IF (IC(1).EQ.1.AND.IC(2).EQ.-3) THEN
          ICOLOR(1,L2)=0
          ICOLOR(2,L2)=300+L1
          ICOLOR(1,L3)=300+L1
          ICOLOR(2,L3)=0
          ICOLOR(1,L4)=0
          ICOLOR(2,L4)=0
      END IF
C     1 -> 1 3 3*
      IF (IC12.EQ.1.AND.IC(3).EQ.3) THEN
          ICOLOR(1,L2)=0
          ICOLOR(2,L2)=0
          ICOLOR(1,L3)=300+L1
          ICOLOR(2,L3)=0
          ICOLOR(1,L4)=0
          ICOLOR(2,L4)=300+L1
      ELSE IF (IC12.EQ.1.AND.IC(3).EQ.-3) THEN
          ICOLOR(1,L2)=0
          ICOLOR(2,L2)=0
          ICOLOR(1,L3)=0
          ICOLOR(2,L3)=300+L1
          ICOLOR(1,L4)=300+L1
          ICOLOR(2,L4)=0
      END IF
C     3 -> 3 1 1
      IF (IC(1).EQ.3.AND.IC(2).EQ.3.AND.IC34.EQ.1) THEN
          ICOLOR(1,L2)=ICOLOR(1,L1)
          ICOLOR(2,L2)=0
          ICOLOR(1,L3)=0
          ICOLOR(2,L3)=0
          ICOLOR(1,L4)=0
          ICOLOR(2,L4)=0
      END IF
C     3* -> 3* 1 1
      IF (IC(1).EQ.-3.AND.IC(2).EQ.-3.AND.IC34.EQ.1) THEN
          ICOLOR(1,L2)=0
          ICOLOR(2,L2)=ICOLOR(2,L1)
          ICOLOR(1,L3)=0
          ICOLOR(2,L3)=0
          ICOLOR(1,L4)=0
          ICOLOR(2,L4)=0
      END IF
C     3 -> 3 3* 3
C     These next two decays seem only necessary for top in isalhe, 
C     which goes t-> q+qb+b in the Isajet decay table ISADECAY.DAT
      IF (IC(1).EQ.3.AND.IC(2).EQ.3.AND.IC34.EQ.-9) THEN
          ICOLOR(1,L2)=300+L1
          ICOLOR(2,L2)=0
          ICOLOR(1,L3)=0
          ICOLOR(2,L3)=300+L1
          ICOLOR(1,L4)=ICOLOR(1,L1)
          ICOLOR(2,L4)=0
      END IF
C     3* -> 3* 3 3*
      IF (IC(1).EQ.-3.AND.IC(2).EQ.-3.AND.IC34.EQ.-9) THEN
          ICOLOR(1,L2)=0
          ICOLOR(2,L2)=300+L1
          ICOLOR(1,L3)=300+L1
          ICOLOR(2,L3)=0
          ICOLOR(1,L4)=0
          ICOLOR(2,L4)=ICOLOR(2,L1)
      END IF
C     3 -> 1 1 3
      IF (IC(1).EQ.3.AND.IC(4).EQ.3.AND.IC23.EQ.1) THEN
          ICOLOR(1,L2)=0
          ICOLOR(2,L2)=0
          ICOLOR(1,L3)=0
          ICOLOR(2,L3)=0
          ICOLOR(1,L4)=ICOLOR(1,L1)
          ICOLOR(2,L4)=0
      END IF
C     3* -> 1 1 3*
      IF (IC(1).EQ.-3.AND.IC(4).EQ.-3.AND.IC23.EQ.1) THEN
          ICOLOR(1,L2)=0
          ICOLOR(2,L2)=0
          ICOLOR(1,L3)=0
          ICOLOR(2,L3)=0
          ICOLOR(1,L4)=0
          ICOLOR(2,L4)=ICOLOR(2,L1)
      END IF
C     8 -> 3 3* 1
      IF (IC(1).EQ.8.AND.IC(2).EQ.3) THEN
          ICOLOR(1,L2)=ICOLOR(1,L1)
          ICOLOR(2,L2)=0
          ICOLOR(1,L3)=0
          ICOLOR(2,L3)=ICOLOR(2,L1)
          ICOLOR(1,L4)=0
          ICOLOR(2,L4)=0
      END IF
C     8 -> 3* 3 1
      IF (IC(1).EQ.8.AND.IC(2).EQ.-3) THEN
          ICOLOR(1,L2)=0
          ICOLOR(2,L2)=ICOLOR(2,L1)
          ICOLOR(1,L3)=ICOLOR(1,L1)
          ICOLOR(2,L3)=0
          ICOLOR(1,L4)=0
          ICOLOR(2,L4)=0
      END IF
C     8 -> 1 3 3*
      IF (IC12.EQ.8.AND.IC(3).EQ.3) THEN
          ICOLOR(1,L2)=0
          ICOLOR(2,L2)=0
          ICOLOR(1,L3)=ICOLOR(1,L1)
          ICOLOR(2,L3)=0
          ICOLOR(1,L4)=0
          ICOLOR(2,L4)=ICOLOR(2,L1)
      END IF
C     8 -> 1 3* 3
      IF (IC12.EQ.8.AND.IC(3).EQ.-3) THEN
          ICOLOR(1,L2)=0
          ICOLOR(2,L2)=0
          ICOLOR(1,L3)=0
          ICOLOR(2,L3)=ICOLOR(2,L1)
          ICOLOR(1,L4)=ICOLOR(1,L1)
          ICOLOR(2,L4)=0
      END IF
      RETURN
      END
