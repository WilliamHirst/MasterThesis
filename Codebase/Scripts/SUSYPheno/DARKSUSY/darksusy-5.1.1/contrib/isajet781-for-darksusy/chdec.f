!
      SUBROUTINE CHDEC(Q,G,LAST)
!
!Purpose: To check if we are at the correct scale for decoupling
!         and if so, fill the vector which will be used to change
!         the threshold locations.
!         If called with LAST=1, the thresholds from the previous
!         run will be used.
!
      IMPLICIT NONE
!
      COMMON/DEC/NEWTH,ISADEC,BELOW,NSTEPTHRESH,NLTMT,
     $           THSQ,THSU,THSD,THSL,THSE
      DOUBLE PRECISION NEWTH(20)
      INTEGER ISADEC,BELOW(20),NSTEPTHRESH(19),NLTMT
      INTEGER THSQ(3),THSU(3),THSD(3),THSL(3),THSE(3)
      SAVE/DEC/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/SQEIG/MQVE,MQVA,MUPVE,MUPVA,MDVE,MDVA,MLVE,MLVA,MEVE,MEVA
      DOUBLE COMPLEX MQVE(3,3),MUPVE(3,3),MDVE(3,3),MLVE(3,3),MEVE(3,3)
     $               ,MQVA(3),MUPVA(3),MDVA(3),MLVA(3),MEVA(3)
      SAVE/SQEIG/
!
      COMMON /SQROT/ RQTOT,RUPTOT,RDTOT,RLTOT,RETOT
     $               ,RQSAV,RUPSAV,RDSAV,RLSAV,RESAV
     $               ,OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      DOUBLE COMPLEX RQTOT(3,3),RUPTOT(3,3),RDTOT(3,3)
      DOUBLE COMPLEX RLTOT(3,3),RETOT(3,3)
      DOUBLE COMPLEX RQSAV(2,3,3),RUPSAV(2,3,3),RDSAV(2,3,3)
      DOUBLE COMPLEX RLSAV(2,3,3),RESAV(2,3,3)
      INTEGER OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      SAVE /SQROT/
!
      COMMON/SFMFRZ/MQSAV,MUPSAV,MDSAV,MLSAV,MESAV
      DOUBLE COMPLEX MQSAV(3,4,3),MUPSAV(3,4,3),MDSAV(3,4,3)
      DOUBLE COMPLEX MLSAV(3,4,3),MESAV(3,4,3)
      SAVE/SFMFRZ/
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
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
      DOUBLE COMPLEX MQCURR(3,3),MUCURR(3,3),MDCURR(3,3)
      DOUBLE COMPLEX MQTMP(3,3),MUTMP(3,3),MDTMP(3,3)
      DOUBLE COMPLEX MQQB(3,3),MUQB(3,3),MDQB(3,3)
      DOUBLE COMPLEX MLQB(3,3),MEQB(3,3),CMATMUL
      DOUBLE COMPLEX CID(3,3),GMASS(601),SUM,G(601),COMBR(3,3)
      DOUBLE COMPLEX VALQ(3),VALU(3),VALD(3),VALL(3),VALE(3)
      DOUBLE COMPLEX VLQ(3,3)
      DOUBLE PRECISION Q,MAX
      INTEGER I,J,K,NSQ,NSU,NSD,NSL,NSE,LAST
!
      DATA CID(1,1)/(1.D0,0.D0)/,CID(1,2)/(0.D0,0.D0)/
     $    ,CID(1,3)/(0.D0,0.D0)/
      DATA CID(2,1)/(0.D0,0.D0)/,CID(2,2)/(1.D0,0.D0)/
     $    ,CID(2,3)/(0.D0,0.D0)/
      DATA CID(3,1)/(0.D0,0.D0)/,CID(3,2)/(0.D0,0.D0)/
     $    ,CID(3,3)/(1.D0,0.D0)/
!
!Check if user wanted to keep isajet's thresholds
!
      IF(ISADEC.EQ.1.OR.LAST.EQ.1)THEN
        DO I=1,20
          BELOW(I)=1
        END DO
!
!I can only use this method of checking for the thresholds
!if the location of the thresholds has been set in advance.
!The running will hit each threshold exactly, and EPS will be
!set in the calling routine DOWNMHIGHMZ.
!
        DO I=1,3
          THSQ(I)=1
          IF((Q-QTHQL(I)).LT.-ABS(EPS).OR.
     $           (ABS(Q-QTHQL(I)).LT.ABS(EPS).AND.EPS.LT.0))THSQ(I)=0
          THSU(I)=1
          IF((Q-QTHUR(I)).LT.-ABS(EPS).OR.
     $           (ABS(Q-QTHUR(I)).LT.ABS(EPS).AND.EPS.LT.0))THSU(I)=0
          THSD(I)=1
          IF((Q-QTHDR(I)).LT.-ABS(EPS).OR.
     $           (ABS(Q-QTHDR(I)).LT.ABS(EPS).AND.EPS.LT.0))THSD(I)=0
          THSL(I)=1
          IF((Q-QTHLL(I)).LT.-ABS(EPS).OR.
     $           (ABS(Q-QTHLL(I)).LT.ABS(EPS).AND.EPS.LT.0))THSL(I)=0
          THSE(I)=1
          IF((Q-QTHER(I)).LT.-ABS(EPS).OR.
     $           (ABS(Q-QTHER(I)).LT.ABS(EPS).AND.EPS.LT.0))THSE(I)=0
        END DO
!
!Otherwise, this is the way each threshold is checked
!
      ELSE
        DO I=1,3
          IF(BELOW(I).EQ.0.AND.Q.LT.DSQRT(ABS(MQVA(I)))
     $                           .AND.DBLE(MQVA(I)).GT.0.D0)THEN
            BELOW(I)=1
            NEWTH(I)=Q
            THSQ(I)=0
          END IF
          IF(BELOW(I+3).EQ.0.AND.Q.LT.DSQRT(ABS(MUPVA(I)))
     $                           .AND.DBLE(MUPVA(I)).GT.0.D0)THEN
            BELOW(I+3)=1
            NEWTH(I+3)=Q
            THSU(I)=0
          END IF
          IF(BELOW(I+6).EQ.0.AND.Q.LT.DSQRT(ABS(MDVA(I)))
     $                           .AND.DBLE(MDVA(I)).GT.0.D0)THEN
            BELOW(I+6)=1
            NEWTH(I+6)=Q
            THSD(I)=0
          END IF
          IF(BELOW(I+9).EQ.0.AND.Q.LT.DSQRT(ABS(MLVA(I)))
     $                           .AND.DBLE(MLVA(I)).GT.0.D0)THEN
            BELOW(I+9)=1
            NEWTH(I+9)=Q
            THSL(I)=0
          END IF
          IF(BELOW(I+12).EQ.0.AND.Q.LT.DSQRT(ABS(MEVA(I)))
     $                           .AND.DBLE(MEVA(I)).GT.0.D0)THEN
            BELOW(I+12)=1
            NEWTH(I+12)=Q
            THSE(I)=0
          END IF
        END DO
!
        IF(BELOW(16).EQ.0.AND.Q.LT.ABS(G(108)))THEN
          BELOW(16)=1
          NEWTH(16)=Q
        END IF
        IF(Q.LT.QNSG)BELOW(17)=1
        IF(Q.LT.QNH)BELOW(18)=1
        IF(BELOW(19).EQ.0.AND.Q.LT.ABS(G(31)))THEN
          BELOW(19)=1
          NEWTH(19)=Q
        END IF
        IF(BELOW(20).EQ.0.AND.Q.LT.ABS(G(32)))THEN
          BELOW(20)=1
          NEWTH(20)=Q
        END IF
      END IF
!
!Now I can set the new Ns and save the frozen eigenvector
!
      NSQ=0
      NSU=0
      NSD=0
      NSL=0
      NSE=0
      DO I=1,3
        NSQ=NSQ+THSQ(I)
        NSU=NSU+THSU(I)
        NSD=NSD+THSD(I)
        NSL=NSL+THSL(I)
        NSE=NSE+THSE(I)
      END DO
!
!Find the quark basis squark mass matrices for use later
!Note that for mq we need to choose which rotation
!we are using, and we use the choice made in the input
!file.
!
      DO I=1,3
        DO J=1,3
          IF(SVLQ.EQ.1)THEN
            VLQ(I,J)=VLU(I,J)
          ELSE
            VLQ(I,J)=VLD(I,J)
          END IF
        END DO
      END DO
      IF(NSQ+NSU+NSD+NSL+NSE.LT.15.AND.
     $   (OLDNSQ-NSQ.NE.0.OR.OLDNSU-NSU.NE.0.OR.OLDNSD-NSD.NE.0.OR.
     $    OLDNSL-NSL.NE.0.OR.OLDNSE-NSE.NE.0))THEN
        DO I=1,3
          DO J=1,3
            MQCURR(I,J)=G(62+(I-1)*3+J)
            MUCURR(I,J)=G(80+(I-1)*3+J)
            MDCURR(I,J)=G(89+(I-1)*3+J)
            MLQB(I,J)=G(71+(I-1)*3+J)
            MEQB(I,J)=G(98+(I-1)*3+J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            MQTMP(I,J)=CMATMUL(0,MQCURR,VLQ,I,J)
            MUTMP(I,J)=CMATMUL(0,MUCURR,VRU,I,J)
            MDTMP(I,J)=CMATMUL(0,MDCURR,VRD,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            MQQB(I,J)=CMATMUL(1,VLQ,MQTMP,I,J)
            MUQB(I,J)=CMATMUL(1,VRU,MUTMP,I,J)
            MDQB(I,J)=CMATMUL(1,VRD,MDTMP,I,J)
          END DO
        END DO
      END IF
!
!Find the largest off-diagonal entry for use by ORTH
!
      IF(OLDNSQ-NSQ.NE.0)THEN
        OFFMAXQVAL=(0.D0,0.D0)
        DO I=1,3
          DO J=1,3
            IF(I.NE.J)THEN
              IF(ABS(OFFMAXQVAL).LE.ABS(MQQB(I,J)))THEN
                OFFMAXQ(1)=I
                OFFMAXQ(2)=J
                OFFMAXQVAL=MQQB(I,J)
              END IF
            END IF
          END DO
        END DO
      END IF
      IF(OLDNSU-NSU.NE.0)THEN
        OFFMAXUVAL=(0.D0,0.D0)
        DO I=1,3
          DO J=1,3
            IF(I.NE.J)THEN
              IF(ABS(OFFMAXUVAL).LE.ABS(MUQB(I,J)))THEN
                OFFMAXU(1)=I
                OFFMAXU(2)=J
                OFFMAXUVAL=MUQB(I,J)
              END IF
            END IF
          END DO
        END DO
      END IF
      IF(OLDNSD-NSD.NE.0)THEN
        OFFMAXDVAL=(0.D0,0.D0)
        DO I=1,3
          DO J=1,3
            IF(I.NE.J)THEN
              IF(ABS(OFFMAXDVAL).LE.ABS(MDQB(I,J)))THEN
                OFFMAXD(1)=I
                OFFMAXD(2)=J
                OFFMAXDVAL=MDQB(I,J)
              END IF
            END IF
          END DO
        END DO
      END IF
      IF(OLDNSL-NSL.NE.0)THEN
        OFFMAXLVAL=(0.D0,0.D0)
        DO I=1,3
          DO J=1,3
            IF(I.NE.J)THEN
              IF(ABS(OFFMAXLVAL).LE.ABS(MLQB(I,J)))THEN
                OFFMAXL(1)=I
                OFFMAXL(2)=J
                OFFMAXLVAL=MLQB(I,J)
              END IF
            END IF
          END DO
        END DO
      END IF
      IF(OLDNSE-NSE.NE.0)THEN
        OFFMAXEVAL=(0.D0,0.D0)
        DO I=1,3
          DO J=1,3
            IF(I.NE.J)THEN
              IF(ABS(OFFMAXEVAL).LE.ABS(MEQB(I,J)))THEN
                OFFMAXE(1)=I
                OFFMAXE(2)=J
                OFFMAXEVAL=MEQB(I,J)
              END IF
            END IF
          END DO
        END DO
      END IF
!
!Now store the eigenvalues for use by ORTH
!
      IF(OLDNSQ.EQ.3)THEN
        DO I=1,3
          VALQ(I)=MQVA(I)
        END DO
      ELSE IF(OLDNSQ.EQ.2)THEN
        DO I=1,2
          VALQ(I)=MQVA(I)
        END DO
        VALQ(3)=MQSAV(3,4,3)
      ELSE IF(OLDNSQ.EQ.1)THEN
        VALQ(1)=MQVA(1)
        VALQ(2)=MQSAV(2,4,2)
        VALQ(3)=MQSAV(3,4,3)
      END IF
!
      IF(OLDNSU.EQ.3)THEN
        DO I=1,3
          VALU(I)=MUPVA(I)
        END DO
      ELSE IF(OLDNSU.EQ.2)THEN
        DO I=1,2
          VALU(I)=MUPVA(I)
        END DO
        VALU(3)=MUPSAV(3,4,3)
      ELSE IF(OLDNSU.EQ.1)THEN
        VALU(1)=MUPVA(1)
        VALU(2)=MUPSAV(2,4,2)
        VALU(3)=MUPSAV(3,4,3)
      END IF
!
      IF(OLDNSD.EQ.3)THEN
        DO I=1,3
          VALD(I)=MDVA(I)
        END DO
      ELSE IF(OLDNSD.EQ.2)THEN
        DO I=1,2
          VALD(I)=MDVA(I)
        END DO
        VALD(3)=MDSAV(3,4,3)
      ELSE IF(OLDNSD.EQ.1)THEN
        VALD(1)=MDVA(1)
        VALD(2)=MDSAV(2,4,2)
        VALD(3)=MDSAV(3,4,3)
      END IF
!
      IF(OLDNSL.EQ.3)THEN
        DO I=1,3
          VALL(I)=MLVA(I)
        END DO
      ELSE IF(OLDNSL.EQ.2)THEN
        DO I=1,2
          VALL(I)=MLVA(I)
        END DO
        VALL(3)=MLSAV(3,4,3)
      ELSE IF(OLDNSL.EQ.1)THEN
        VALL(1)=MLVA(1)
        VALL(2)=MLSAV(2,4,2)
        VALL(3)=MLSAV(3,4,3)
      END IF
!
      IF(OLDNSE.EQ.3)THEN
        DO I=1,3
          VALE(I)=MEVA(I)
        END DO
      ELSE IF(OLDNSE.EQ.2)THEN
        DO I=1,2
          VALE(I)=MEVA(I)
        END DO
        VALE(3)=MESAV(3,4,3)
      ELSE IF(OLDNSE.EQ.1)THEN
        VALE(1)=MEVA(1)
        VALE(2)=MESAV(2,4,2)
        VALE(3)=MESAV(3,4,3)
      END IF
!
!Only set the Rs to be different from identity if we are below
!a threshold.
!
      DO I=1,3
        DO J=1,3
          RQTOT(I,J)=CID(I,J)
          RUPTOT(I,J)=CID(I,J)
          RDTOT(I,J)=CID(I,J)
          RLTOT(I,J)=CID(I,J)
          RETOT(I,J)=CID(I,J)
        END DO
      END DO
!
!Set the rotation if we are below at least one threshold. The first
!dimension in the saved matrix is the number of active squarks in
!the region we use the rotation.
!For only 1 active squark, we need to find the compound rotation.
!MATMUL cannot be used due to the nature of the matrix with three
!indices.
!Each time I save a rotation, I need to carry out my fix using ORTH
!
      IF(NSQ+NSU+NSD+NSL+NSE.LT.15)THEN
        IF(NSQ.LT.3.AND.OLDNSQ.NE.0)THEN
          IF(OLDNSQ.EQ.3)THEN
            CALL ORTH(VALQ,MQVE,OFFMAXQ,OFFMAXQVAL,VLQ,0)
          ELSE
            DO I=1,3
              DO J=1,3
                SUM=(0.D0,0.D0)
                DO K=1,3
                  SUM=SUM+RQSAV(2,I,K)*MQVE(K,J)
                END DO
                COMBR(I,J)=SUM
              END DO
            END DO
            CALL ORTH(VALQ,COMBR,OFFMAXQ,OFFMAXQVAL,VLQ,0)
          END IF
          IF(OLDNSQ.EQ.3.AND.OLDNSQ-NSQ.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RQSAV(2,I,J)=MQVE(I,J)
                IF(NSQ.LT.2)RQSAV(1,I,J)=MQVE(I,J)
              END DO
            END DO
          ELSE IF(OLDNSQ.EQ.2.AND.OLDNSQ-NSQ.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RQSAV(1,I,J)=COMBR(I,J)
              END DO
            END DO
          END IF
          IF(NSQ.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RQTOT(I,J)=RQSAV(NSQ,I,J)
              END DO
            END DO
          END IF
        END IF
!
        IF(NSU.LT.3.AND.OLDNSU.NE.0)THEN
          IF(OLDNSU.EQ.3)THEN
            CALL ORTH(VALU,MUPVE,OFFMAXU,OFFMAXUVAL,VRU,0)
          ELSE
            DO I=1,3
              DO J=1,3
                SUM=(0.D0,0.D0)
                DO K=1,3
                  SUM=SUM+RUPSAV(2,I,K)*MUPVE(K,J)
                END DO
                COMBR(I,J)=SUM
              END DO
            END DO
            CALL ORTH(VALU,COMBR,OFFMAXU,OFFMAXUVAL,VRU,0)
          END IF
          IF(OLDNSU.EQ.3.AND.OLDNSU-NSU.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RUPSAV(2,I,J)=MUPVE(I,J)
                IF(NSU.LT.2)RUPSAV(1,I,J)=MUPVE(I,J)
              END DO
            END DO
          ELSE IF(OLDNSU.EQ.2.AND.OLDNSU-NSU.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RUPSAV(1,I,J)=COMBR(I,J)
              END DO
            END DO
          END IF
          IF(NSU.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RUPTOT(I,J)=RUPSAV(NSU,I,J)
              END DO
            END DO
          END IF
        END IF
!
        IF(NSD.LT.3.AND.OLDNSD.NE.0)THEN
          IF(OLDNSD.EQ.3)THEN
            CALL ORTH(VALD,MDVE,OFFMAXD,OFFMAXDVAL,VRD,0)
          ELSE
            DO I=1,3
              DO J=1,3
                SUM=(0.D0,0.D0)
                DO K=1,3
                  SUM=SUM+RDSAV(2,I,K)*MDVE(K,J)
                END DO
                COMBR(I,J)=SUM
              END DO
            END DO
            CALL ORTH(VALD,COMBR,OFFMAXD,OFFMAXDVAL,VRD,0)
          END IF
          IF(OLDNSD.EQ.3.AND.OLDNSD-NSD.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RDSAV(2,I,J)=MDVE(I,J)
                IF(NSD.LT.2)RDSAV(1,I,J)=MDVE(I,J)
              END DO
            END DO
          ELSE IF(OLDNSD.EQ.2.AND.OLDNSD-NSD.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RDSAV(1,I,J)=COMBR(I,J)
              END DO
            END DO
          END IF
          IF(NSD.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RDTOT(I,J)=RDSAV(NSD,I,J)
              END DO
            END DO
          END IF
        END IF
!
        IF(NSL.LT.3.AND.OLDNSL.NE.0)THEN
          IF(OLDNSL.EQ.3)THEN
            CALL ORTH(VALL,MLVE,OFFMAXL,OFFMAXLVAL,CID,0)
          ELSE
            DO I=1,3
              DO J=1,3
                SUM=(0.D0,0.D0)
                DO K=1,3
                  SUM=SUM+RLSAV(2,I,K)*MLVE(K,J)
                END DO
                COMBR(I,J)=SUM
              END DO
            END DO
            CALL ORTH(VALL,COMBR,OFFMAXL,OFFMAXLVAL,CID,0)
          END IF
          IF(OLDNSL.EQ.3.AND.OLDNSL-NSL.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RLSAV(2,I,J)=MLVE(I,J)
                IF(NSL.LT.2)RLSAV(1,I,J)=MLVE(I,J)
              END DO
            END DO
          ELSE IF(OLDNSL.EQ.2.AND.OLDNSL-NSL.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RLSAV(1,I,J)=COMBR(I,J)
              END DO
            END DO
          END IF
          IF(NSL.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RLTOT(I,J)=RLSAV(NSL,I,J)
              END DO
            END DO
          END IF
        END IF
!
        IF(NSE.LT.3.AND.OLDNSE.NE.0)THEN
          IF(OLDNSE.EQ.3)THEN
            CALL ORTH(VALE,MEVE,OFFMAXE,OFFMAXEVAL,CID,0)
          ELSE
            DO I=1,3
              DO J=1,3
                SUM=(0.D0,0.D0)
                DO K=1,3
                  SUM=SUM+RESAV(2,I,K)*MEVE(K,J)
                END DO
                COMBR(I,J)=SUM
              END DO
            END DO
            CALL ORTH(VALE,COMBR,OFFMAXE,OFFMAXEVAL,CID,0)
          END IF
          IF(OLDNSE.EQ.3.AND.OLDNSE-NSE.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RESAV(2,I,J)=MEVE(I,J)
                IF(NSE.LT.2)RESAV(1,I,J)=MEVE(I,J)
              END DO
            END DO
          ELSE IF(OLDNSE.EQ.2.AND.OLDNSE-NSE.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RESAV(1,I,J)=COMBR(I,J)
              END DO
            END DO
          END IF
          IF(NSE.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RETOT(I,J)=RESAV(NSE,I,J)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!Save the values of the sfermion matrices at each scale
!to use when running up. The first dimension of each array
!is the number of active sfermions above the threshold.
!The fourth entry in the second index is for the eigenvalues.
!
      DO I=1,3
        IF(OLDNSQ-NSQ.NE.0)THEN
          DO K=OLDNSQ,NSQ+1,-1
            DO J=1,3
              MQSAV(K,I,J)=G(62+(I-1)*3+J)
            END DO
            MQSAV(K,4,I)=VALQ(I)
          END DO
        END IF
        IF(OLDNSU-NSU.NE.0)THEN
          DO K=OLDNSU,NSU+1,-1
            DO J=1,3
              MUPSAV(K,I,J)=G(80+(I-1)*3+J)
            END DO
            MUPSAV(K,4,I)=VALU(I)
          END DO
        END IF
        IF(OLDNSD-NSD.NE.0)THEN
          DO K=OLDNSD,NSD+1,-1
            DO J=1,3
              MDSAV(K,I,J)=G(89+(I-1)*3+J)
            END DO
            MDSAV(K,4,I)=VALD(I)
          END DO
        END IF
        IF(OLDNSL-NSL.NE.0)THEN
          DO K=OLDNSL,NSL+1,-1
            DO J=1,3
              MLSAV(K,I,J)=G(71+(I-1)*3+J)
            END DO
            MLSAV(K,4,I)=VALL(I)
          END DO
        END IF
        IF(OLDNSE-NSE.NE.0)THEN
          DO K=OLDNSE,NSE+1,-1
            DO J=1,3
              MESAV(K,I,J)=G(98+(I-1)*3+J)
            END DO
            MESAV(K,4,I)=VALE(I)
          END DO
        END IF
      END DO
!
!Construct the squark matrices with saved entries when some squarks
!have decoupled.
!
      CALL MASSSQM(G)
!
!Rotate into the squark mass basis, remove the appropriate terms and
!rotate back.
!
      IF(NSQ-OLDNSQ.NE.0.OR.NSU-OLDNSU.NE.0.OR.NSD-OLDNSD.NE.0
     $     .OR.NSL-OLDNSL.NE.0.OR.NSE-OLDNSE.NE.0)THEN
        CALL REMSF(G)
      END IF
!
      RETURN
      END
