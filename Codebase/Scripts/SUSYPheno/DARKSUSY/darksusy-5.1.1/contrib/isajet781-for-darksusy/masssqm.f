!
      SUBROUTINE MASSSQM(GRUN)
!
!Purpose: To construct the quark mass basis soft mass matrices.
!
      IMPLICIT NONE
!
!Needed to save the values of mq and mu for use in decay calculations
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
      COMMON/SQEIG/MQVE,MQVA,MUPVE,MUPVA,MDVE,MDVA,MLVE,MLVA,MEVE,MEVA
      DOUBLE COMPLEX MQVE(3,3),MUPVE(3,3),MDVE(3,3),MLVE(3,3),MEVE(3,3)
     $               ,MQVA(3),MUPVA(3),MDVA(3),MLVA(3),MEVA(3)
      SAVE/SQEIG/
!
      COMMON/SFMFRZ/MQSAV,MUPSAV,MDSAV,MLSAV,MESAV
      DOUBLE COMPLEX MQSAV(3,4,3),MUPSAV(3,4,3),MDSAV(3,4,3)
      DOUBLE COMPLEX MLSAV(3,4,3),MESAV(3,4,3)
      SAVE/SFMFRZ/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/DEC/NEWTH,ISADEC,BELOW,NSTEPTHRESH,NLTMT,
     $           THSQ,THSU,THSD,THSL,THSE
      DOUBLE PRECISION NEWTH(20)
      INTEGER ISADEC,BELOW(20),NSTEPTHRESH(19),NLTMT
      INTEGER THSQ(3),THSU(3),THSD(3),THSL(3),THSE(3)
      SAVE/DEC/
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
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
      DOUBLE COMPLEX EQ(3,3),EU(3,3),ED(3,3),EL(3,3),EE(3,3)
      DOUBLE COMPLEX MQMASS(3,3),MUMASS(3,3),MDMASS(3,3)
      DOUBLE COMPLEX MLMASS(3,3),MEMASS(3,3)
      DOUBLE COMPLEX MQTMP(3,3),MUTMP(3,3),MDTMP(3,3)
      DOUBLE COMPLEX MLTMP(3,3),METMP(3,3)
      DOUBLE COMPLEX MQCURR(3,3),MUCURR(3,3),MDCURR(3,3)
      DOUBLE COMPLEX MLCURR(3,3),MECURR(3,3)
      DOUBLE COMPLEX SUM,RQTMP(3,3),RUTMP(3,3),RDTMP(3,3)
      DOUBLE COMPLEX RLTMP(3,3),RETMP(3,3)
      DOUBLE COMPLEX GRUN(601),CMATMUL,CID(3,3),VLQ(3,3)
      DOUBLE COMPLEX EVAQ(3),EVAU(3),EVAD(3),EVAL(3),EVAE(3)
      INTEGER I,J,K
!
      DATA CID(1,1)/(1.D0,0.D0)/,CID(1,2)/(0.D0,0.D0)/
     $    ,CID(1,3)/(0.D0,0.D0)/
      DATA CID(2,1)/(0.D0,0.D0)/,CID(2,2)/(1.D0,0.D0)/
     $    ,CID(2,3)/(0.D0,0.D0)/
      DATA CID(3,1)/(0.D0,0.D0)/,CID(3,2)/(0.D0,0.D0)/
     $    ,CID(3,3)/(1.D0,0.D0)/
!
      DO I=1,3
        DO J=1,3
          EQ(I,J)=CID(I,J)
          EU(I,J)=CID(I,J)
          ED(I,J)=CID(I,J)
          EL(I,J)=CID(I,J)
          EE(I,J)=CID(I,J)
          MQMASS(I,J)=CID(I,J)
          MUMASS(I,J)=CID(I,J)
          MDMASS(I,J)=CID(I,J)
          MLMASS(I,J)=CID(I,J)
          MEMASS(I,J)=CID(I,J)
          MQTMP(I,J)=CID(I,J)
          MUTMP(I,J)=CID(I,J)
          MDTMP(I,J)=CID(I,J)
          MLTMP(I,J)=CID(I,J)
          METMP(I,J)=CID(I,J)
          RQTMP(I,J)=CID(I,J)
          RUTMP(I,J)=CID(I,J)
          RDTMP(I,J)=CID(I,J)
          RLTMP(I,J)=CID(I,J)
          RETMP(I,J)=CID(I,J)
        END DO
        EVAQ(I)=(0.D0,0.D0)
        EVAU(I)=(0.D0,0.D0)
        EVAD(I)=(0.D0,0.D0)
        EVAL(I)=(0.D0,0.D0)
        EVAE(I)=(0.D0,0.D0)
      END DO
!
!Use the familiar rotation for the doublet.
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
!
!Use OLD* for the number of active squarks since we still haven't
!taken a step at this new scale
!
!
!If there are two active squarks, I need to find the combined rotation
!first from the squark mass basis to the basis where the squarks are
!diagonal at the heaviest decoupling and then from this basis to
!the original current basis. If the second squark has also decoupled
!this combined rotation is saved elsewhere.
!
      IF(OLDNSQ.LT.3.OR.OLDNSU.LT.3.OR.OLDNSD.LT.3
     $                     .OR.OLDNSL.LT.3.OR.OLDNSE.LT.3)THEN
        DO I=1,3
          DO J=1,3
            IF(OLDNSQ.EQ.2)THEN
              SUM=(0.D0,0.D0)
              DO K=1,3
                SUM=SUM+RQSAV(2,I,K)*MQVE(K,J)
              END DO
              RQTMP(I,J)=SUM
            END IF
            IF(OLDNSU.EQ.2)THEN
              SUM=(0.D0,0.D0)
              DO K=1,3
                SUM=SUM+RUPSAV(2,I,K)*MUPVE(K,J)
              END DO
              RUTMP(I,J)=SUM
            END IF          
            IF(OLDNSD.EQ.2)THEN
              SUM=(0.D0,0.D0)
              DO K=1,3
                SUM=SUM+RDSAV(2,I,K)*MDVE(K,J)
              END DO
              RDTMP(I,J)=SUM
            END IF          
            IF(OLDNSL.EQ.2)THEN
              SUM=(0.D0,0.D0)
              DO K=1,3
                SUM=SUM+RLSAV(2,I,K)*MLVE(K,J)
              END DO
              RLTMP(I,J)=SUM
            END IF          
            IF(OLDNSE.EQ.2)THEN
              SUM=(0.D0,0.D0)
              DO K=1,3
                SUM=SUM+RESAV(2,I,K)*MEVE(K,J)
              END DO
              RETMP(I,J)=SUM
            END IF          
          END DO
        END DO
!
!Work out the eigenvectors in terms of the current basis
!states. If there has been no decoupling, then we can just use
!the current basis soft mass matrices as they are.
!
!Also in this loop, work out the mass matrix in the squark mass
!basis.
!
        DO I=1,3
          IF(OLDNSQ.LT.3)THEN
            IF(OLDNSQ.EQ.2)THEN
              DO J=1,3
                EQ(I,J)=RQTMP(I,J)
              END DO
              IF(I.LT.3)EVAQ(I)=MQVA(I)
            ELSE IF(OLDNSQ.LT.2)THEN
              DO J=1,3
                EQ(I,J)=RQSAV(1,I,J)
              END DO
              IF(OLDNSQ.EQ.1.AND.I.EQ.1)EVAQ(I)=MQVA(I)
              IF(OLDNSQ.EQ.0.AND.I.EQ.1)EVAQ(I)=MQSAV(I,4,I)
              IF(I.EQ.2)EVAQ(I)=MQSAV(I,4,I)
            END IF
            IF(I.EQ.3)EVAQ(I)=MQSAV(I,4,I)
          END IF
!
          IF(OLDNSU.LT.3)THEN
            IF(OLDNSU.EQ.2)THEN
              DO J=1,3
                EU(I,J)=RUTMP(I,J)
              END DO
              IF(I.LT.3)EVAU(I)=MUPVA(I)
            ELSE IF(OLDNSU.LT.2)THEN
              DO J=1,3
                EU(I,J)=RUPSAV(1,I,J)
              END DO
              IF(OLDNSU.EQ.1.AND.I.EQ.1)EVAU(I)=MUPVA(I)
              IF(OLDNSU.EQ.0.AND.I.EQ.1)EVAU(I)=MUPSAV(I,4,I)
              IF(I.EQ.2)EVAU(I)=MUPSAV(I,4,I)
            END IF
            IF(I.EQ.3)EVAU(I)=MUPSAV(I,4,I)
          END IF
!
          IF(OLDNSD.LT.3)THEN
            IF(OLDNSD.EQ.2)THEN
              DO J=1,3
                ED(I,J)=RDTMP(I,J)
              END DO
              IF(I.LT.3)EVAD(I)=MDVA(I)
            ELSE IF(OLDNSD.LT.2)THEN
              DO J=1,3
                ED(I,J)=RDSAV(1,I,J)
              END DO
              IF(OLDNSD.EQ.1.AND.I.EQ.1)EVAD(I)=MDVA(I)
              IF(OLDNSD.EQ.0.AND.I.EQ.1)EVAD(I)=MDSAV(I,4,I)
              IF(I.EQ.2)EVAD(I)=MDSAV(I,4,I)
            END IF
            IF(I.EQ.3)EVAD(I)=MDSAV(I,4,I)
          END IF
!
          IF(OLDNSL.LT.3)THEN
            IF(OLDNSL.EQ.2)THEN
              DO J=1,3
                EL(I,J)=RLTMP(I,J)
              END DO
              IF(I.LT.3)EVAL(I)=MLVA(I)
            ELSE IF(OLDNSL.LT.2)THEN
              DO J=1,3
                EL(I,J)=RLSAV(1,I,J)
              END DO
              IF(OLDNSL.EQ.1.AND.I.EQ.1)EVAL(I)=MLVA(I)
              IF(OLDNSL.EQ.0.AND.I.EQ.1)EVAL(I)=MLSAV(I,4,I)
              IF(I.EQ.2)EVAL(I)=MLSAV(I,4,I)
            END IF
            IF(I.EQ.3)EVAL(I)=MLSAV(I,4,I)
          END IF
!
          IF(OLDNSE.LT.3)THEN
            IF(OLDNSE.EQ.2)THEN
              DO J=1,3
                EE(I,J)=RETMP(I,J)
              END DO
              IF(I.LT.3)EVAE(I)=MEVA(I)
            ELSE IF(OLDNSE.LT.2)THEN
              DO J=1,3
                EE(I,J)=RESAV(1,I,J)
              END DO
              IF(OLDNSE.EQ.1.AND.I.EQ.1)EVAE(I)=MEVA(I)
              IF(OLDNSE.EQ.0.AND.I.EQ.1)EVAE(I)=MESAV(I,4,I)
              IF(I.EQ.2)EVAE(I)=MESAV(I,4,I)
            END IF
            IF(I.EQ.3)EVAE(I)=MESAV(I,4,I)
          END IF
        END DO
!
!Fix the rotation matrices using ORTH as before so that 
!the rotation error is in the least damaging entries.
!
        IF(OLDNSQ.LT.3)CALL ORTH(EVAQ,EQ,OFFMAXQ,OFFMAXQVAL,VLQ,0)
        IF(OLDNSU.LT.3)CALL ORTH(EVAU,EU,OFFMAXU,OFFMAXUVAL,VRU,0)
        IF(OLDNSD.LT.3)CALL ORTH(EVAD,ED,OFFMAXD,OFFMAXDVAL,VRD,0)
        IF(OLDNSL.LT.3)CALL ORTH(EVAL,EL,OFFMAXL,OFFMAXLVAL,CID,0)
        IF(OLDNSE.LT.3)CALL ORTH(EVAE,EE,OFFMAXE,OFFMAXEVAL,CID,0)
!
        DO I=1,3
          MQMASS(I,I)=EVAQ(I)
          MUMASS(I,I)=EVAU(I)
          MDMASS(I,I)=EVAD(I)
          MLMASS(I,I)=EVAL(I)
          MEMASS(I,I)=EVAE(I)
        END DO
!
!I must now rotate from the squark mass basis to
!the original current basis and then to the quark mass basis.
!
        DO I=1,3
          DO J=1,3
            MQTMP(I,J)=CMATMUL(2,MQMASS,EQ,I,J)
            MUTMP(I,J)=CMATMUL(2,MUMASS,EU,I,J)
            MDTMP(I,J)=CMATMUL(2,MDMASS,ED,I,J)
            MLTMP(I,J)=CMATMUL(2,MLMASS,EL,I,J)
            METMP(I,J)=CMATMUL(2,MEMASS,EE,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            MQCURR(I,J)=CMATMUL(0,EQ,MQTMP,I,J)
            MUCURR(I,J)=CMATMUL(0,EU,MUTMP,I,J)
            MDCURR(I,J)=CMATMUL(0,ED,MDTMP,I,J)
            MLCURR(I,J)=CMATMUL(0,EL,MLTMP,I,J)
            MECURR(I,J)=CMATMUL(0,EE,METMP,I,J)
          END DO
        END DO
      END IF
!
      IF(OLDNSQ.EQ.3)THEN
        DO I=1,3
          DO J=1,3
            MQCURR(I,J)=GRUN(62+(I-1)*3+J)
          END DO
        END DO
      END IF
      IF(OLDNSU.EQ.3)THEN
        DO I=1,3
          DO J=1,3
            MUCURR(I,J)=GRUN(80+(I-1)*3+J)
          END DO
        END DO
      END IF
      IF(OLDNSD.EQ.3)THEN
        DO I=1,3
          DO J=1,3
            MDCURR(I,J)=GRUN(89+(I-1)*3+J)
          END DO
        END DO
      END IF
      IF(OLDNSL.EQ.3)THEN
        DO I=1,3
          DO J=1,3
            MLCURR(I,J)=GRUN(71+(I-1)*3+J)
          END DO
        END DO
      END IF
      IF(OLDNSE.EQ.3)THEN
        DO I=1,3
          DO J=1,3
            MECURR(I,J)=GRUN(98+(I-1)*3+J)
          END DO
        END DO
      END IF
!
      DO I=1,3
        DO J=1,3
          MQTMP(I,J)=CMATMUL(0,MQCURR,VLQ,I,J)
          MUTMP(I,J)=CMATMUL(0,MUCURR,VRU,I,J)
          MDTMP(I,J)=CMATMUL(0,MDCURR,VRD,I,J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          MQQMASS(I,J)=CMATMUL(1,VLQ,MQTMP,I,J)
          MUQMASS(I,J)=CMATMUL(1,VRU,MUTMP,I,J)
          MDQMASS(I,J)=CMATMUL(1,VRD,MDTMP,I,J)
          MLQMASS(I,J)=MLCURR(I,J)
          MEQMASS(I,J)=MECURR(I,J)
        END DO
      END DO
!
   51 FORMAT(1X,A3,3X,D27.20,3x,D27.20)
   52 FORMAT(1X,I1,I1,3X,D27.20,3x,D27.20)
   53 FORMAT(1X,D27.20,3x,D27.20)
!
      RETURN
      END
