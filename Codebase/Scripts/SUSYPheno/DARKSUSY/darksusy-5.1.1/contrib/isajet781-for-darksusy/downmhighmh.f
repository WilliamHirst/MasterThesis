!
      SUBROUTINE DOWNMHIGHMH(QEND,NSTEP)
!
!Purpose: The final run down before stop decay routine is called.
!         NOTE: If m_H is less than M_t, only run to M_t and 
!               make sure that this value is passed to SQSIX
!
      IMPLICIT NONE
!
      COMMON /BSG/GISA(31),MSQISA(3),MSLISA(3),MSUISA(3),MSDISA(3),
     &            MSEISA(3),MRNISA(3),YNFRZ(3,3),MNFRZ(3,3),TNFRZ(3,3),
     &            RTISA,RBISA,RLISA
c     MSxDEC(i) - decoupling scale of i-th generation of type x sfermion
c     MRNDEC(i) - decoupling scale of i-th RH neutrino
      REAL*8 GISA,MSQISA,MSLISA,MSUISA,MSDISA,MSEISA,MRNISA,
     &       YNFRZ,MNFRZ,TNFRZ
      REAL RTISA,RBISA,RLISA
      SAVE /BSG/
!
      COMMON/LOOPS/SSQSTEP,SW2LP
      DOUBLE PRECISION SSQSTEP
      INTEGER SW2LP
      SAVE/LOOPS/
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
!
      COMMON/RGEMS/VEVMH,RGEMS,RGEMU
      DOUBLE COMPLEX VEVMH
      DOUBLE PRECISION RGEMS,RGEMU
      SAVE/RGEMS/
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
      COMMON/DEC/NEWTH,ISADEC,BELOW,NSTEPTHRESH,NLTMT,
     $           THSQ,THSU,THSD,THSL,THSE
      DOUBLE PRECISION NEWTH(20)
      INTEGER ISADEC,BELOW(20),NSTEPTHRESH(19),NLTMT
      INTEGER THSQ(3),THSU(3),THSD(3),THSL(3),THSE(3)
      SAVE/DEC/
!
      DOUBLE PRECISION TTH(20),THIGH,TMT,T,QEND,SINB,COSB
      DOUBLE PRECISION DW(1803),THSP,PI,DT
      DOUBLE COMPLEX B,VU,VD
      DOUBLE COMPLEX W(1803),ID(3,3)
      INTEGER NSTEP,I,J,II,NHTMT,BELOWMS
      EXTERNAL CRGE601,DRGE601
!
      DATA ID(1,1)/(1.D0,0.D0)/,ID(1,2)/(0.D0,0.D0)/
     $ ,ID(1,3)/(0.D0,0.D0)/
      DATA ID(2,1)/(0.D0,0.D0)/,ID(2,2)/(1.D0,0.D0)/
     $ ,ID(2,3)/(0.D0,0.D0)/
      DATA ID(3,1)/(0.D0,0.D0)/,ID(3,2)/(0.D0,0.D0)/
     $ ,ID(3,3)/(1.D0,0.D0)/
!
      BELOWMS=0
      PI=4.D0*DATAN(1.D0)
!
      THIGH=0.D0
      DO I=1,20
        TTH(I)=LOG(QTHSORT(I)/MHIGH)
      END DO
      TMT=LOG(MT/MHIGH)
!
!Check location of lowest threshold above m_t
!
      NHTMT=NLTMT+1
!
      IF(NHTMT.NE.21)THEN
        DT=(TTH(20)-THIGH)/FLOAT(NSTEP)
      ELSE
        DT=(TMT-THIGH)/FLOAT(NSTEP)
      END IF
!
      DO II=1,NSTEP
        IF(NHTMT.NE.21)THEN
          T=THIGH+(TTH(20)-THIGH)*FLOAT(II-1)/FLOAT(NSTEP)
        ELSE
          T=THIGH+(TMT-THIGH)*FLOAT(II-1)/FLOAT(NSTEP)
        END IF
        SSQSTEP=MHIGH*EXP(T)
!
        IF(BELOWMS.EQ.0.AND.SSQSTEP.LT.RGEMS)THEN
          BELOWMS=1
          CALL ROTBACK(1)
          G(12)=G(12)*(1.D0-DBLE(RTISA))
          G(21)=G(21)*(1.D0-DBLE(RBISA))
          G(30)=G(30)*(1.D0-DBLE(RLISA))
          IF(RGEMS.LE.QNH)THEN
            G(120)=G(120)*(1.D0-DBLE(RTISA))
            G(129)=G(129)*(1.D0-DBLE(RBISA))
            G(138)=G(138)*(1.D0-DBLE(RLISA))
          END IF
          CALL ROTATE(1)
        END IF
!
        IF(II.EQ.1)EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
!
        CALL SQDIAG(G)
        CALL CHDEC(SSQSTEP,G,1)
!
        IF(COMP.EQ.0)THEN
          DO I=1,601
            DG(I)=DBLE(G(I))
          END DO
          CALL DRKSTP(601,DT,T,DG,DRGE601,DW)
          DO I=1,601
            G(I)=DCMPLX(DG(I))
          END DO
        ELSE
          CALL CRKSTP(601,DT,T,G,CRGE601,W)
        END IF
!
        EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
      END DO
!
      IF(NHTMT.NE.21)THEN
        DO I=19,1,-1
          IF(NSTEPTHRESH(I).EQ.0)GOTO 40
          IF(I.LT.LOCMH)GOTO 50
          DT=(TTH(I)-TTH(I+1))/FLOAT(NSTEPTHRESH(I))
!
          DO II=1,NSTEPTHRESH(I)
            T=TTH(I+1)+(TTH(I)-TTH(I+1))*FLOAT(II-1)
     $                                  /FLOAT(NSTEPTHRESH(I))
            SSQSTEP=MHIGH*EXP(T)
!
            IF(BELOWMS.EQ.0.AND.SSQSTEP.LT.RGEMS)THEN
              BELOWMS=1
              CALL ROTBACK(1)
              G(12)=G(12)*(1.D0-DBLE(RTISA))
              G(21)=G(21)*(1.D0-DBLE(RBISA))
              G(30)=G(30)*(1.D0-DBLE(RLISA))
              IF(RGEMS.LE.QNH)THEN
                G(120)=G(120)*(1.D0-DBLE(RTISA))
                G(129)=G(129)*(1.D0-DBLE(RBISA))
                G(138)=G(138)*(1.D0-DBLE(RLISA))
              END IF
              CALL ROTATE(1)
            END IF
!
            IF(II.EQ.1)EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
            CALL SQDIAG(G)
            CALL CHDEC(SSQSTEP,G,1)
            IF(COMP.EQ.0)THEN
              DO J=1,601
                DG(J)=DBLE(G(J))
              END DO
              CALL DRKSTP(601,DT,T,DG,DRGE601,DW)
              DO J=1,601
                G(J)=DCMPLX(DG(J))
              END DO
            ELSE
              CALL CRKSTP(601,DT,T,G,CRGE601,W)
            END IF
            EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
          END DO
!
  40    END DO
  50    CONTINUE
      END IF
      IF(NHTMT.EQ.21)THEN
        WRITE(*,*)'WARNING: M_H IS LESS THAN M_T'
      END IF
!
!SM assignments at m_H
!
      SINB=DSQRT(TANB**2/(1+TANB**2))
      COSB=SINB/TANB
      B=G(109)
      VU=G(110)
      VD=G(111)
!
      G(429)=(3.D0/5.D0*(G(1))**2+(G(2))**2)/4.D0
     $                                  *((1+TANB**2)/(1-TANB**2))**2
      G(428)=DSQRT(DBLE(VU*CONJG(VU))+DBLE(VD*CONJG(VD)))
      DO I=1,3
        DO J=1,3
          G(111+(I-1)*3+J)=SINB*G(3+(I-1)*3+J)
          G(120+(I-1)*3+J)=COSB*G(12+(I-1)*3+J)
          G(129+(I-1)*3+J)=COSB*G(21+(I-1)*3+J)
          G(399+(I-1)*3+J)=SINB*G(33+(I-1)*3+J)-COSB*G(429+(I-1)*3+J)
          G(408+(I-1)*3+J)=COSB*G(42+(I-1)*3+J)-SINB*G(438+(I-1)*3+J)
          G(417+(I-1)*3+J)=COSB*G(51+(I-1)*3+J)-SINB*G(447+(I-1)*3+J)
        END DO
      END DO
      G(287)=SINB*G(184)
      G(288)=COSB*G(185)
      G(289)=SINB*G(204)
      G(290)=COSB*G(205)
      G(427)=SINB**2*G(61)+COSB**2*G(62)-2.D0*SINB*COSB*B
!
      QEND=MHIGH*EXP(T)
!
      RETURN
      END
