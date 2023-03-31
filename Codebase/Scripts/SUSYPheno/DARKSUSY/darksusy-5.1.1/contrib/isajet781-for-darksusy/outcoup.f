!
      SUBROUTINE OUTCOUP(Q,UD)
!
!Purpose: To output each matrix to a different file, in either
!         the 'mass' or current basis
!
!         UD=1 if we are running up, =0 if we are running down.
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
      COMMON/EWSBSAV/CSHSH,CSHLH,CLHLH,CULSHUR,CDLSHDR,CELSHER,CULLHUR,
     $               CDLLHDR,CELLHER
      DOUBLE COMPLEX CSHSH,CSHLH,CLHLH,CULSHUR(3,3),CDLSHDR(3,3),
     $               CELSHER(3,3),CULLHUR(3,3),CDLLHDR(3,3),CELLHER(3,3)
      SAVE/EWSBSAV/
!
      COMMON/MYDECAY/MQQMASS,MUQMASS,MDQMASS,MLQMASS,MEQMASS,
     $             OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $             OFFMAXEVAL,OFFMAXQ,OFFMAXU,OFFMAXD,OFFMAXL,OFFMAXE
      DOUBLE COMPLEX MQQMASS(3,3),MUQMASS(3,3),MDQMASS(3,3),
     $               MLQMASS(3,3),MEQMASS(3,3)
      DOUBLE COMPLEX OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $               OFFMAXEVAL
      INTEGER OFFMAXQ(2),OFFMAXU(2),OFFMAXD(2),OFFMAXL(2),OFFMAXE(2)
      SAVE/MYDECAY/
!
      DOUBLE COMPLEX GTMP(601)
      DOUBLE PRECISION Q,COSB,SINB
      INTEGER I,J,SWROT,UD
!
      SINB=DSQRT(TANB**2/(1+TANB**2))
      COSB=SINB/TANB
!
!Switch which allows the user to turn off the rotation.
!
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
!      IF(SWROT.EQ.1)CALL ROTBACK(0)
      IF(SWROT.EQ.1)CALL ROTBACK(1) !This line is used to go back to the 
                                    !true mass basis.
!
!Insert the soft masses
!
      DO I=1,3
        DO J=1,3
          G(62+(I-1)*3+J)=MQQMASS(I,J)
          G(71+(I-1)*3+J)=MLQMASS(I,J)
          G(80+(I-1)*3+J)=MUQMASS(I,J)
          G(89+(I-1)*3+J)=MDQMASS(I,J)
          G(98+(I-1)*3+J)=MEQMASS(I,J)
        END DO
      END DO
!
!Make some adjustments to facilitate comparisons
!
      G(1)=G(1)
      DO I=4,12
        G(I+108)=G(I+108)/SINB
      END DO
      DO I=13,30
        G(I+108)=G(I+108)/COSB
      END DO
      G(287)=G(287)/SINB
      G(289)=G(289)/SINB
      G(288)=G(288)/COSB
      G(290)=G(290)/COSB
      IF(Q.LE.QNH-EPS)THEN
        DO I=1,27
          IF(I.LT.10)THEN
            G(I+323)=SINB*G(I+323)-COSB*CONJG(G(398))*G(I+293)
          ELSE
            G(I+323)=COSB*G(I+323)-SINB*CONJG(G(398))*G(I+293)
          END IF
        END DO
      END IF
      G(351)=G(351)+ABS(G(398))**2
      G(352)=G(352)+ABS(G(398))**2
      DO I=1,601
        IF(ABS(G(I)).LT.1.D-30)THEN
          G(I)=(1.D-30,1.D-30)
        END IF
      END DO
!
!Print each set of terms out - absolute squares first
!
      WRITE(68,56)Q,ABS(G(1)),DBLE(G(1))
      WRITE(67,56)Q,ABS(G(2)),DBLE(G(2))
      WRITE(66,56)Q,ABS(G(3)),DBLE(G(3))
      WRITE(62,56)Q,ABS(G(31)),DBLE(G(31))
      WRITE(52,55)Q,ABS(G(31)),ABS(G(32)),ABS(G(32))/ABS(G(31))
     $             ,ABS(G(321)),ABS(G(322)),ABS(G(322))/ABS(G(321))
     $             ,abs(g(2))**2/(abs(g(1))/dsqrt(3.d0/5.d0))**2
     $             ,4.D0*ABS(G(2))**2*ABS(G(108))*sinb*cosb/
     $                              (12.D0*ABS(G(32))*ABS(G(2))**2)
      WRITE(61,56)Q,ABS(G(32)),DBLE(G(32))
      WRITE(60,56)Q,ABS(G(33)),DBLE(G(33))
      WRITE(55,56)Q,ABS(G(110)),G(110)
      WRITE(54,56)Q,ABS(G(111)),G(111)
      IF((Q.GE.QNH-ABS(EPS).AND.UD.EQ.0).OR.
     $                  (Q.GE.QNH+ABS(EPS).AND.UD.EQ.1))THEN
        I=4
        WRITE(63,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=13
        WRITE(64,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=22
        WRITE(65,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=34
        WRITE(71,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=43
        WRITE(72,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=52
        WRITE(73,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=61
        WRITE(74,56)Q,ABS(G(I)),G(I)
        I=62
        WRITE(75,56)Q,ABS(G(I)),G(I)
        I=430
        WRITE(83,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=439
        WRITE(84,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=448
        WRITE(85,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=109
        WRITE(82,56)Q,ABS(G(I)),G(I)
      ELSE
        I=112
        WRITE(63,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=121
        WRITE(64,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=130
        WRITE(65,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=400
        WRITE(86,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=409
        WRITE(87,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=418
        WRITE(88,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=427
        WRITE(89,56)Q,ABS(G(I)),G(I)
      END IF
      I=63
      WRITE(76,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=72
      WRITE(79,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=81
      WRITE(77,57)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
     $,abs(2.d0*(g(81)-g(85))/(g(81)+g(85)))
      IF((G(89)-G(85)).NE.0)WRITE(53,56)Q,ABS(G(86)/(G(89)-G(85)))
     $,G(86)/(G(89)-G(85))
      I=90
      WRITE(78,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=99
      WRITE(80,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=108
      WRITE(81,56)Q,ABS(G(I)),G(I)
      I=139
      WRITE(59,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=157
      WRITE(58,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=186
      WRITE(57,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=233
      WRITE(56,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=324
      WRITE(90,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=333
      WRITE(91,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=342
      WRITE(92,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=351
      WRITE(93,55)Q,ABS(G(I)),G(I)
      I=352
      WRITE(94,55)Q,ABS(G(I)),G(I)
      I=353
      WRITE(95,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=362
      WRITE(98,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=371
      WRITE(96,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=380
      WRITE(97,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=389
      WRITE(99,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=398
      WRITE(70,56)Q,ABS(G(I)),G(I)
      I=399
      WRITE(69,56)Q,ABS(G(I)),G(I)
!
!Return the gs to their original values
!
      DO I=1,601
        G(I)=GTMP(I)
      END DO
!
   51 FORMAT(1X,A3,3X,D27.20)
   55 FORMAT(SP,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10
     $,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X
     $,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X
     $,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X
     $,D17.10,1X,D17.10,1X,D17.10,1X,D17.10)
   57 FORMAT(SP,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10
     $,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X
     $,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X
     $,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X
     $,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10)
   56 FORMAT(SP,D17.10,1X,D17.10,1X,D17.10,1X,D17.10)
!
      RETURN
      END
