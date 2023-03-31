!
      SUBROUTINE CHINT(Q,G)
!
!Purpose: To check if the sfermions must be reintroduced
!         into the running. If so, we must use the saved values
!         of the sfermion mass matrices as a boundary
!         condition on the upwards running.
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
      DOUBLE COMPLEX CID(3,3),G(601)
      DOUBLE PRECISION Q
      INTEGER I,J,K,NSQ,NSU,NSD,NSL,NSE
!
      DATA CID(1,1)/(1.D0,0.D0)/,CID(1,2)/(0.D0,0.D0)/
     $    ,CID(1,3)/(0.D0,0.D0)/
      DATA CID(2,1)/(0.D0,0.D0)/,CID(2,2)/(1.D0,0.D0)/
     $    ,CID(2,3)/(0.D0,0.D0)/
      DATA CID(3,1)/(0.D0,0.D0)/,CID(3,2)/(0.D0,0.D0)/
     $    ,CID(3,3)/(1.D0,0.D0)/
!
!We can use the same statements as before to set the various thetas.
!This is because we know the position of each threshold and the
!calling routine is hitting each one exactly.
!
      DO I=1,3
        THSQ(I)=1
        IF((Q-QTHQL(I)).LT.-ABS(EPS).OR.
     $         (ABS(Q-QTHQL(I)).LT.ABS(EPS).AND.EPS.LT.0))THSQ(I)=0
        THSU(I)=1
        IF((Q-QTHUR(I)).LT.-ABS(EPS).OR.
     $         (ABS(Q-QTHUR(I)).LT.ABS(EPS).AND.EPS.LT.0))THSU(I)=0
        THSD(I)=1
        IF((Q-QTHDR(I)).LT.-ABS(EPS).OR.
     $         (ABS(Q-QTHDR(I)).LT.ABS(EPS).AND.EPS.LT.0))THSD(I)=0
        THSL(I)=1
        IF((Q-QTHLL(I)).LT.-ABS(EPS).OR.
     $         (ABS(Q-QTHLL(I)).LT.ABS(EPS).AND.EPS.LT.0))THSL(I)=0
        THSE(I)=1
        IF((Q-QTHER(I)).LT.-ABS(EPS).OR.
     $         (ABS(Q-QTHER(I)).LT.ABS(EPS).AND.EPS.LT.0))THSE(I)=0
      END DO
!
!Next set the values of the number of active sfermions
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
!Only set the Rs to be different from identity if we are below
!a threshold. I do not need to reset the entries of the saved
!matrices as in CHDEC...
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
      IF(NSQ+NSU+NSD+NSL+NSE.LT.15.AND.NSQ+NSU+NSD+NSL+NSE.GT.0)THEN
        DO I=1,3
          DO J=1,3
            IF(NSQ.LT.3.AND.NSQ.GT.0)RQTOT(I,J)=RQSAV(NSQ,I,J)
            IF(NSU.LT.3.AND.NSU.GT.0)RUPTOT(I,J)=RUPSAV(NSU,I,J)
            IF(NSD.LT.3.AND.NSD.GT.0)RDTOT(I,J)=RDSAV(NSD,I,J)
            IF(NSL.LT.3.AND.NSL.GT.0)RLTOT(I,J)=RLSAV(NSL,I,J)
            IF(NSE.LT.3.AND.NSE.GT.0)RETOT(I,J)=RESAV(NSE,I,J)
          END DO
        END DO
      END IF
!
!Now insert the values of the sfermion matrices at each scale.
!
      DO I=1,3
        DO J=1,3
          IF(OLDNSQ-NSQ.NE.0)G(62+(I-1)*3+J)=MQSAV(NSQ,I,J)
          IF(OLDNSU-NSU.NE.0)G(80+(I-1)*3+J)=MUPSAV(NSU,I,J)
          IF(OLDNSD-NSD.NE.0)G(89+(I-1)*3+J)=MDSAV(NSD,I,J)
          IF(OLDNSL-NSL.NE.0)G(71+(I-1)*3+J)=MLSAV(NSL,I,J)
          IF(OLDNSE-NSE.NE.0)G(98+(I-1)*3+J)=MESAV(NSE,I,J)
        END DO
      END DO
!
      RETURN
      END
