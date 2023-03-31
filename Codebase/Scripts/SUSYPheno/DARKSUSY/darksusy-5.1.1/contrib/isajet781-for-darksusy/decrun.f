!
      SUBROUTINE DECRUN(GIN,GOUT,QIN,SVLQ)
!
!Purpose: To run the inputs to from m_H to the required scale.
!
!         Uses knowledge of the decoupling from RGEFLAV such as
!         decoupling points and frozen values of couplings, masses
!         and squark mixings.
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
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/LOOPS/SSQSTEP,SW2LP
      DOUBLE PRECISION SSQSTEP
      INTEGER SW2LP
      SAVE/LOOPS/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/RGEMS/VEVMH,RGEMS,RGEMU
      DOUBLE COMPLEX VEVMH
      DOUBLE PRECISION RGEMS,RGEMU
      SAVE/RGEMS/
!
      COMMON/DECCALC/T1EVE,T1EVA,USQM,COSTHT,SINTHT,GHIK,MST1,MST2,GAMMA
      DOUBLE COMPLEX T1EVE(6,6),T1EVA(6),USQM(6,6),COSTHT,SINTHT
     $              ,GHIK(601)
      DOUBLE PRECISION MST1,MST2,GAMMA
      SAVE/DECCALC/
!
      DOUBLE COMPLEX GIN(601),GRUN(601),GOUT(601),GFIX(601),W(1803)
      DOUBLE PRECISION DGRUN(601),DW(1803),QIN,TQIN
      DOUBLE PRECISION DT,T,TMH,TTH(20),TMDEC,Q,MDEC,PI
      INTEGER I,J,II,NSTEP,STTH,FINTH,NSTEPTH,BELOWMS,SVLQ
      EXTERNAL CRGE601,DRGE601
!
!Reinstate this (and certain lines below) to write out the running
!of m^2_U according to DECRUN using DECOUTCOUP.
!
!      OPEN(63,FILE='out/stdown.dat',STATUS='UNKNOWN')
      IF(QIN.LT.RGEMS)THEN
        BELOWMS=1
      ELSE
        BELOWMS=0
      END IF
!
      MDEC=QTHQL(1)
      IF(QTHUR(1).LT.QTHQL(1).AND.SVLQ.EQ.1)MDEC=QTHUR(1)
      IF(QTHDR(1).LT.QTHQL(1).AND.SVLQ.EQ.0)MDEC=QTHDR(1)
!
!Prevent running from going below mt
!
      IF(MDEC.LT.MT)THEN
        WRITE(33,*)'ATTEMPTED TO RUN BELOW MT IN DECRUN'
        MDEC=MT
      END IF
!
      PI=4.D0*DATAN(1.D0)
!
      CALL STROTATE(GIN,GRUN,0)
!
!The number of steps is fixed and large, but could be allowed
!to float or be equal to the numbers in RGEFLAV.
!
      NSTEP=1000
      NSTEPTH=500
      TMH=0.D0
      TQIN=LOG(QIN/QNH)
      TMDEC=LOG(MDEC/QNH)
      DO I=1,20
        TTH(I)=LOG(QTHSORT(I)/QNH)
      END DO
!
!Find start and end points. FINTH is threshold below MDEC
!
      STTH=LOCMH
      FINTH=20
      DO I=20,1,-1
        IF(MDEC.LT.QTHSORT(I))FINTH=FINTH-1
      END DO
!
      IF(MDEC.LT.QIN)THEN
!
!Running down past the relevant thresholds
!
        IF(STTH.NE.FINTH+1)THEN
          DO I=STTH,FINTH+2,-1
            IF(QTHSORT(I).EQ.QTHSORT(I-1))GOTO 41
            DT=(TTH(I-1)-TTH(I))/FLOAT(NSTEPTH)
            DO II=1,NSTEPTH
              T=TTH(I)+(TTH(I-1)-TTH(I))*FLOAT(II-1)/FLOAT(NSTEPTH)
              SSQSTEP=QNH*EXP(T)
!
!Finite Pierce correction
!
              IF(BELOWMS.EQ.0.AND.SSQSTEP.LT.RGEMS)THEN
                BELOWMS=1
                CALL STROTBACK(GRUN,GFIX,1)
                GFIX(12)=GFIX(12)*(1.D0-DBLE(RTISA))
                GFIX(21)=GFIX(21)*(1.D0-DBLE(RBISA))
                GFIX(30)=GFIX(30)*(1.D0-DBLE(RLISA))
                IF(RGEMS.LT.QNH)THEN
                  GFIX(120)=GFIX(120)*(1.D0-DBLE(RTISA))
                  GFIX(129)=GFIX(129)*(1.D0-DBLE(RBISA))
                  GFIX(138)=GFIX(138)*(1.D0-DBLE(RLISA))
                END IF
                CALL STROTATE(GFIX,GRUN,1)
              END IF
!
              IF(II.EQ.1)EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
              CALL SQDIAG(GRUN)
              CALL CHDEC(SSQSTEP,GRUN,1)
              IF(COMP.EQ.0)THEN
                DO J=1,601
                  DGRUN(J)=DBLE(GRUN(J))
                END DO
                CALL DRKSTP(601,DT,T,DGRUN,DRGE601,DW)
                DO J=1,601
                  GRUN(J)=DCMPLX(DGRUN(J))
                END DO
              ELSE
                CALL CRKSTP(601,DT,T,GRUN,CRGE601,W)
              END IF
              EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
!              CALL DECOUTCOUP(QNH*EXP(T),GRUN)
            END DO
  41      END DO
        END IF
!
!Now the final run to MDEC
!
        DT=(TMDEC-TTH(FINTH+1))/FLOAT(NSTEP)
        DO II=1,NSTEP
          T=TTH(FINTH+1)+(TMDEC-TTH(FINTH+1))*FLOAT(II-1)/FLOAT(NSTEP)
          SSQSTEP=QNH*EXP(T)
          Q=QNH*EXP(T)
!
          IF(BELOWMS.EQ.0.AND.SSQSTEP.LT.RGEMS)THEN
            BELOWMS=1
            CALL STROTBACK(GRUN,GFIX,1)
            GFIX(12)=GFIX(12)*(1.D0-DBLE(RTISA))
            GFIX(21)=GFIX(21)*(1.D0-DBLE(RBISA))
            GFIX(30)=GFIX(30)*(1.D0-DBLE(RLISA))
            IF(RGEMS.LT.QNH)THEN
              GFIX(120)=GFIX(120)*(1.D0-DBLE(RTISA))
              GFIX(129)=GFIX(129)*(1.D0-DBLE(RBISA))
              GFIX(138)=GFIX(138)*(1.D0-DBLE(RLISA))
            END IF
            CALL STROTATE(GFIX,GRUN,1)
          END IF
!
          IF(II.EQ.1)EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
          CALL SQDIAG(GRUN)
          CALL CHDEC(SSQSTEP,GRUN,1)
          IF(COMP.EQ.0)THEN
            DO J=1,601
              DGRUN(J)=DBLE(GRUN(J))
            END DO
            CALL DRKSTP(601,DT,T,DGRUN,DRGE601,DW)
            DO J=1,601
              GRUN(J)=DCMPLX(DGRUN(J))
            END DO
          ELSE
            CALL CRKSTP(601,DT,T,GRUN,CRGE601,W)
          END IF
          EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
!          CALL DECOUTCOUP(QNH*EXP(T),GRUN)
        END DO
      ELSE IF(MDEC.GT.QIN)THEN
!
!Running up past the relevant thresholds
!
        IF(STTH.NE.FINTH)THEN
          DO I=STTH,FINTH-1
            IF(QTHSORT(I).EQ.QTHSORT(I+1))GOTO 42
            DT=(TTH(I+1)-TTH(I))/FLOAT(NSTEPTH)
            DO II=1,NSTEPTH
              T=TTH(I)+(TTH(I+1)-TTH(I))*FLOAT(II-1)/FLOAT(NSTEPTH)
              SSQSTEP=QNH*EXP(T)
!
              IF(BELOWMS.EQ.1.AND.QNH*EXP(T+DT).GT.RGEMS)THEN
                BELOWMS=0
                CALL STROTBACK(GRUN,GFIX,1)
                GFIX(12)=GFIX(12)/(1.D0-DBLE(RTISA))
                GFIX(21)=GFIX(21)/(1.D0-DBLE(RBISA))
                GFIX(30)=GFIX(30)/(1.D0-DBLE(RLISA))
                IF(RGEMS.LT.QNH)THEN
                  GFIX(120)=GFIX(120)/(1.D0-DBLE(RTISA))
                  GFIX(129)=GFIX(129)/(1.D0-DBLE(RBISA))
                  GFIX(138)=GFIX(138)/(1.D0-DBLE(RLISA))
                END IF
                CALL STROTATE(GFIX,GRUN,1)
              END IF
!
              IF(II.EQ.1)EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
              CALL CHINT(SSQSTEP,GRUN)
              IF(COMP.EQ.0)THEN
                DO J=1,601
                  DGRUN(J)=DBLE(GRUN(J))
                END DO
                CALL DRKSTP(601,DT,T,DGRUN,DRGE601,DW)
                DO J=1,601
                  GRUN(J)=DCMPLX(DGRUN(J))
                END DO
              ELSE
                CALL CRKSTP(601,DT,T,GRUN,CRGE601,W)
              END IF
              EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
            END DO
  42      END DO
        END IF
!
!Now the final run to MDEC
!
        DT=(TMDEC-TTH(FINTH))/FLOAT(NSTEP)
        DO II=1,NSTEP
          T=TTH(FINTH)+(TMDEC-TTH(FINTH))*FLOAT(II-1)/FLOAT(NSTEP)
          SSQSTEP=QNH*EXP(T)
          Q=QNH*EXP(T)
!
          IF(BELOWMS.EQ.1.AND.QNH*EXP(T+DT).GT.RGEMS)THEN
            BELOWMS=0
            CALL STROTBACK(GRUN,GFIX,1)
            GFIX(12)=GFIX(12)/(1.D0-DBLE(RTISA))
            GFIX(21)=GFIX(21)/(1.D0-DBLE(RBISA))
            GFIX(30)=GFIX(30)/(1.D0-DBLE(RLISA))
            IF(RGEMS.LT.QNH)THEN
              GFIX(120)=GFIX(120)/(1.D0-DBLE(RTISA))
              GFIX(129)=GFIX(129)/(1.D0-DBLE(RBISA))
              GFIX(138)=GFIX(138)/(1.D0-DBLE(RLISA))
            END IF
            CALL STROTATE(GFIX,GRUN,1)
          END IF
!
          IF(II.EQ.1)EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
          CALL CHINT(SSQSTEP,GRUN)
          IF(COMP.EQ.0)THEN
            DO J=1,601
              DGRUN(J)=DBLE(GRUN(J))
            END DO
            CALL DRKSTP(601,DT,T,DGRUN,DRGE601,DW)
            DO J=1,601
              GRUN(J)=DCMPLX(DGRUN(J))
            END DO
          ELSE
            CALL CRKSTP(601,DT,T,GRUN,CRGE601,W)
          END IF
          EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
        END DO
      ELSE
        T=TQIN
      END IF
!
      QIN=QNH*EXP(T)
!
!Work out the final squark mass matrices in the quark mass basis
!
      CALL SQDIAG(GRUN)
      CALL MASSSQM(GRUN)
!
      CALL STROTBACK(GRUN,GOUT,0)
!
!Save the new values of GHIK if m_H<MDEC
!
      IF(QNH.LT.QIN)THEN
        DO I=1,601
          GHIK(I)=GOUT(I)
        END DO
      END IF
!
      CLOSE(63)
!
      RETURN
      END
