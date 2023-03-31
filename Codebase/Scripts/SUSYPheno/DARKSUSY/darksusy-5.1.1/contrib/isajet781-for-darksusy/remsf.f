!
      SUBROUTINE REMSF(G)
!
!Purpose: To remove the entry in the matter sfermion mass matrix
!         which corresponds to the decoupled particle, and
!         return new entries for current basis and mass basis g's
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
      DOUBLE COMPLEX G(601),CMATMUL
      DOUBLE COMPLEX MQ(3,3),MUP(3,3),MD(3,3),ML(3,3),ME(3,3)
      DOUBLE COMPLEX MQT(3,3),MUPT(3,3),MDT(3,3),MLT(3,3),MET(3,3)
      DOUBLE COMPLEX MQMASS(3,3),MUPMASS(3,3),MDMASS(3,3)
      DOUBLE COMPLEX MLMASS(3,3),MEMASS(3,3)
      INTEGER I,J,NSQ,NSU,NSD,NSL,NSE
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
!THE NEXT LINE SHOULD BE UNNECESSARY
!It just checks to see that we are running down, not up.
!
      IF(OLDNSQ.LT.NSQ.OR.OLDNSU.LT.NSU.OR.OLDNSD.LT.NSD.OR.
     $   OLDNSL.LT.NSL.OR.OLDNSE.LT.NSE)RETURN
!
!Rotate the relevant matrices
!
      DO I=1,3
        DO J=1,3
          MQ(I,J)=G(62+(I-1)*3+J)
          ML(I,J)=G(71+(I-1)*3+J)
          MUP(I,J)=G(80+(I-1)*3+J)
          MD(I,J)=G(89+(I-1)*3+J)
          ME(I,J)=G(98+(I-1)*3+J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          MQT(I,J)=CMATMUL(0,MQ,RQTOT,I,J)
          MLT(I,J)=CMATMUL(0,ML,RLTOT,I,J)
          MUPT(I,J)=CMATMUL(0,MUP,RUPTOT,I,J)
          MDT(I,J)=CMATMUL(0,MD,RDTOT,I,J)
          MET(I,J)=CMATMUL(0,ME,RETOT,I,J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          MQMASS(I,J)=CMATMUL(1,RQTOT,MQT,I,J)
          MLMASS(I,J)=CMATMUL(1,RLTOT,MLT,I,J)
          MUPMASS(I,J)=CMATMUL(1,RUPTOT,MUPT,I,J)
          MDMASS(I,J)=CMATMUL(1,RDTOT,MDT,I,J)
          MEMASS(I,J)=CMATMUL(1,RETOT,MET,I,J)
        END DO
      END DO
!
!Remove the entries from the mass basis matrices
!
      DO I=1,3
        DO J=1,3
          IF(NSQ.LT.3.AND.OLDNSQ-NSQ.NE.0)THEN
            IF(I.EQ.J)THEN
              MQMASS(I,J)=MQMASS(I,J)*THSQ(I)
            ELSE
              MQMASS(I,J)=(0.D0,0.D0)
            END IF
          END IF
          IF(NSU.LT.3.AND.OLDNSU-NSU.NE.0)THEN
            IF(I.EQ.J)THEN
              MUPMASS(I,I)=MUPMASS(I,J)*THSU(I)
            ELSE
              MUPMASS(I,J)=(0.D0,0.D0)
            END IF
          END IF
          IF(NSD.LT.3.AND.OLDNSD-NSD.NE.0)THEN
            IF(I.EQ.J)THEN
              MDMASS(I,J)=MDMASS(I,J)*THSD(I)
            ELSE
              MDMASS(I,J)=(0.D0,0.D0)
            END IF
          END IF
          IF(NSL.LT.3.AND.OLDNSL-NSL.NE.0)THEN
            IF(I.EQ.J)THEN
              MLMASS(I,J)=MLMASS(I,J)*THSL(I)
            ELSE
              MLMASS(I,J)=(0.D0,0.D0)
            END IF
          END IF
          IF(NSE.LT.3.AND.OLDNSE-NSE.NE.0)THEN
            IF(I.EQ.J)THEN
              MEMASS(I,J)=MEMASS(I,J)*THSE(I)
            ELSE
              MEMASS(I,J)=(0.D0,0.D0)
            END IF
          END IF
        END DO
      END DO
!
!Finally, rotate these matrices back
!
      DO I=1,3
        DO J=1,3
          MQT(I,J)=CMATMUL(2,MQMASS,RQTOT,I,J)
          MLT(I,J)=CMATMUL(2,MLMASS,RLTOT,I,J)
          MUPT(I,J)=CMATMUL(2,MUPMASS,RUPTOT,I,J)
          MDT(I,J)=CMATMUL(2,MDMASS,RDTOT,I,J)
          MET(I,J)=CMATMUL(2,MEMASS,RETOT,I,J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          G(62+(I-1)*3+J)=CMATMUL(0,RQTOT,MQT,I,J)
          G(71+(I-1)*3+J)=CMATMUL(0,RLTOT,MLT,I,J)
          G(80+(I-1)*3+J)=CMATMUL(0,RUPTOT,MUPT,I,J)
          G(89+(I-1)*3+J)=CMATMUL(0,RDTOT,MDT,I,J)
          G(98+(I-1)*3+J)=CMATMUL(0,RETOT,MET,I,J)
        END DO
      END DO
!
      RETURN
      END
