!
      SUBROUTINE DOWNMHIGHMZ(QEND,NSTEP,THCNG)
!
!Purpose: To run the couplings down from M_HIGH to M_H and
!         rotate back into the mass basis
!
!         THCNG is used to switch off the changing of the
!         thresholds before the final run to allow the Yukawas
!         to become diagonal at M_t.
!
      IMPLICIT NONE
!
      COMMON/LOOPS/SSQSTEP,SW2LP
      DOUBLE PRECISION SSQSTEP
      INTEGER SW2LP
      SAVE/LOOPS/
!
      COMMON/SMRGE/SMRGEMH,SMQSTEP,NU,SMDR2LP
      DOUBLE PRECISION SMRGEMH,SMQSTEP
      INTEGER NU,SMDR2LP
      SAVE/SMRGE/
!
      COMMON/ATMZ/G1MZ,G2MZ,G3MZ,VSMMZ,LAMBDAMZ,LAMTMZ
      DOUBLE COMPLEX G1MZ,G2MZ,G3MZ,VSMMZ,LAMBDAMZ,LAMTMZ
      SAVE/ATMZ/
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
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
      COMMON/RGEMS/VEVMH,RGEMS,RGEMU
      DOUBLE COMPLEX VEVMH
      DOUBLE PRECISION RGEMS,RGEMU
      SAVE/RGEMS/
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
!DEC:
!The common block for changing the decoupling points during each run.
!The order of the entries in NEWTH and BELOW are: mq(1-3), mu(1-3),
!md(1-3), ml(1-3), me(1-3), higgsino (mu), m_gluino, m_A0, M_1, M_2.
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/SFMFRZ/MQSAV,MUPSAV,MDSAV,MLSAV,MESAV
      DOUBLE COMPLEX MQSAV(3,4,3),MUPSAV(3,4,3),MDSAV(3,4,3)
      DOUBLE COMPLEX MLSAV(3,4,3),MESAV(3,4,3)
      SAVE/SFMFRZ/
!
      COMMON/SQEIG/MQVE,MQVA,MUPVE,MUPVA,MDVE,MDVA,MLVE,MLVA,MEVE,MEVA
      DOUBLE COMPLEX MQVE(3,3),MUPVE(3,3),MDVE(3,3),MLVE(3,3),MEVE(3,3)
     $               ,MQVA(3),MUPVA(3),MDVA(3),MLVA(3),MEVA(3)
      SAVE/SQEIG/
!
      DOUBLE PRECISION TTH(20),THIGH,DT,T,QEND,DW(1803)
      DOUBLE PRECISION TMT,TZ,PI,DGSM(32),DWSM(96)
      DOUBLE COMPLEX W(1803),ID(3,3),GSM(32),WSM(96)
      INTEGER NSTEP,I,J,K,II,NHTMT,NSTEPMT,LOOPNSTEP,BELOWMS,THCNG
      INTEGER NSTEPSM
      EXTERNAL CRGE601,DRGE601,CSMRGEDR,DSMRGEDR
!
      DATA ID(1,1)/(1.D0,0.D0)/,ID(1,2)/(0.D0,0.D0)/
     $ ,ID(1,3)/(0.D0,0.D0)/
      DATA ID(2,1)/(0.D0,0.D0)/,ID(2,2)/(1.D0,0.D0)/
     $ ,ID(2,3)/(0.D0,0.D0)/
      DATA ID(3,1)/(0.D0,0.D0)/,ID(3,2)/(0.D0,0.D0)/
     $ ,ID(3,3)/(1.D0,0.D0)/
!
!USED BY OUTCOUP. CLOSE AT THE END OF THE PROG.
!
!      OPEN(71,FILE='out/dau.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(72,FILE='out/dad.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(73,FILE='out/dae.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(74,FILE='out/dmhumu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(75,FILE='out/dmhdmu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(76,FILE='out/dmq.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(77,FILE='out/dmup.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(78,FILE='out/dmd.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(79,FILE='out/dml.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(80,FILE='out/dme.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(81,FILE='out/dmu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(82,FILE='out/db.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(83,FILE='out/dmtsfu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(84,FILE='out/dmtsfd.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(85,FILE='out/dmtsfe.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(86,FILE='out/dtriu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(87,FILE='out/dtrid.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(88,FILE='out/dtrie.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(89,FILE='out/dmhud.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(90,FILE='out/daum.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(91,FILE='out/dadm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(92,FILE='out/daem.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(93,FILE='out/dmhumum.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(94,FILE='out/dmhdmum.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(95,FILE='out/dmqm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(96,FILE='out/dmupm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(97,FILE='out/dmdm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(98,FILE='out/dmlm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(99,FILE='out/dmem.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(70,FILE='out/dmum.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(69,FILE='out/dbm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(68,FILE='out/dg1.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(67,FILE='out/dg2.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(66,FILE='out/dg3.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(65,FILE='out/dytau.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(64,FILE='out/dyb.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(63,FILE='out/dyu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(62,FILE='out/dm1.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(61,FILE='out/dm2.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(60,FILE='out/dm3.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(59,FILE='out/dgtpq.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(58,FILE='out/dgtpu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(57,FILE='out/dgtq.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(56,FILE='out/dftuq.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(55,FILE='out/dvu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(54,FILE='out/dvd.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(53,FILE='out/dmuflavcomp.dat',STATUS='UNKNOWN'
!     $ ,FORM='FORMATTED')
!      OPEN(52,FILE='out/dm12comp.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!
      BELOWMS=0
      PI=4.D0*DATAN(1.D0)
!
!Set the initial values for changing the thresholds.
!The new thresholds are set to the old ones in case they
!are not changed.
!
      DO I=1,20
        BELOW(I)=0
      END DO
      DO I=1,3
        NEWTH(I)=QTHQL(I)
        NEWTH(I+3)=QTHUR(I)
        NEWTH(I+6)=QTHDR(I)
        NEWTH(I+9)=QTHLL(I)
        NEWTH(I+12)=QTHER(I)
      END DO
      NEWTH(16)=QNSH
      NEWTH(17)=QNSG
      NEWTH(18)=QNH
      NEWTH(19)=QTHSB
      NEWTH(20)=QTHSW
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
        IF(II.EQ.1)EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
!
!Carry out the squark rotation business to make sure that we are
!in the correct basis for the squark thresholds
!
        CALL SQDIAG(G)
        CALL CHDEC(SSQSTEP,G,THCNG)
!        CALL OUTCOUP(MHIGH*EXP(T),0)
!
!Finite corrections to the Yukawas from Isajet and other
!m_SUSY conditions.
!
        IF(BELOWMS.EQ.0.AND.SSQSTEP.LT.RGEMS)THEN
          BELOWMS=1
          CALL DOWNMSCOND
        END IF
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
        EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
      END DO
!
      IF(NHTMT.NE.21)THEN
!
        IF(LOCMH.EQ.20)CALL DOWNMHCOND
!
        DO I=19,NHTMT,-1
          IF(NSTEPTHRESH(I).EQ.0)GOTO 40
          DT=(TTH(I)-TTH(I+1))/FLOAT(NSTEPTHRESH(I))
!
          DO II=1,NSTEPTHRESH(I)
            T=TTH(I+1)+(TTH(I)-TTH(I+1))*FLOAT(II-1)
     $                                  /FLOAT(NSTEPTHRESH(I))
            SSQSTEP=MHIGH*EXP(T)
            IF(II.EQ.1)EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
            CALL SQDIAG(G)
            CALL CHDEC(SSQSTEP,G,THCNG)
!            CALL OUTCOUP(MHIGH*EXP(T),0)
            IF(BELOWMS.EQ.0.AND.SSQSTEP.LT.RGEMS)THEN
              BELOWMS=1
              CALL DOWNMSCOND
            END IF
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
   40     IF(LOCMH.EQ.I)CALL DOWNMHCOND
        END DO
!
!Continue running to m_t
!
        IF(NHTMT.NE.1)THEN
          DT=(TMT-TTH(NHTMT))/FLOAT(NSTEPTHRESH(NHTMT-1))
          LOOPNSTEP=NSTEPTHRESH(NHTMT-1)
        ELSE
          NSTEPMT=INT(ABS(DLOG(MT/QTHSORT(1)))*NSTEP/DLOG(MHIGH/MT)*25)
          DT=(TMT-TTH(NHTMT))/FLOAT(NSTEPMT)
          LOOPNSTEP=NSTEPMT
        END IF
!
        DO II=1,LOOPNSTEP
          T=TTH(NHTMT)+(TMT-TTH(NHTMT))*FLOAT(II-1)/FLOAT(LOOPNSTEP)
          SSQSTEP=MHIGH*EXP(T)
          IF(II.EQ.1)EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
          CALL SQDIAG(G)
          CALL CHDEC(SSQSTEP,G,THCNG)
!          CALL OUTCOUP(MHIGH*EXP(T),0)
          IF(BELOWMS.EQ.0.AND.SSQSTEP.LT.RGEMS)THEN
            BELOWMS=1
            CALL DOWNMSCOND
          END IF
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
          EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
        END DO
!
      ELSE
        WRITE(*,*)'WARNING: ALL THRESHOLDS ARE LESS THAN MT'
        CALL DOWNMHCOND
      END IF
!
!Finish the running by continuing to M_Z. This is primarily done
!since on the first run up, the value of the quartic coupling in the
!SM was zero, as was lambda_t.
!
      NSTEPSM=50
      TZ=DLOG(MZ/MHIGH)
      DT=(TZ-TMT)/FLOAT(NSTEPSM)
      NU=3 !Decouple the top at M_Z
      DO I=1,3
        GSM(I)=G(I)
      END DO
      DO I=4,30
        GSM(I)=G(108+I)
      END DO
      GSM(31)=G(429)
      GSM(32)=G(428)
!
!Rotate into the mass basis - we will not run the full flavour
!SM RGEs
!
      CALL ROTBACKSM(GSM)
!
      IF(COMP.EQ.0)THEN
        DO I=1,32
          DGSM(I)=DBLE(GSM(I))
        END DO
      END IF
!
      DO II=1,NSTEPSM
        T=TMT+(TZ-TMT)*FLOAT(II-1)/FLOAT(NSTEPSM)
        SMQSTEP=MHIGH*EXP(T)
        IF(II.EQ.1)EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
        IF(COMP.EQ.0)THEN
          CALL DRKSTP(32,DT,T,DGSM,DSMRGEDR,DWSM)
        ELSE
          CALL CRKSTP(32,DT,T,GSM,CSMRGEDR,WSM)
        END IF
        EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
      END DO
!
      IF(COMP.EQ.0)THEN
        DO I=1,32
          GSM(I)=DCMPLX(DGSM(I))
        END DO
      END IF
!
      QEND=MHIGH*EXP(T)
!
!If any of the entries in BELOW are still 0, then
!inform the user.
!
!      DO I=1,20
!        IF(BELOW(I).EQ.0)WRITE(*,*)'THRESH ',I,' IS BELOW MT'
!      END DO
!
!Save the weak scale SM values of the quartic coupling and up Yukawa
!
      LAMBDAMZ=GSM(31)
      LAMTMZ=GSM(12)
!
!Now replace the old thresholds with the new ones
!
      DO I=1,3
        QTHQL(I)=NEWTH(I)
        QTHUR(I)=NEWTH(I+3)
        QTHDR(I)=NEWTH(I+6)
        QTHLL(I)=NEWTH(I+9)
        QTHER(I)=NEWTH(I+12)
      END DO
      QNSH=NEWTH(16)
      QNSG=NEWTH(17)
      QNH=NEWTH(18)
      QTHSB=NEWTH(19)
      QTHSW=NEWTH(20)
      CALL SORTTH
!
!Recheck the location of the first threshold above m_t
!
      NLTMT=0
      DO I=1,20
        IF(QTHSORT(I).LE.MT)NLTMT=I
      END DO
!
!Put out warning for all sfermion masses less than 0.D0
!
      DO I=1,3
        IF(QTHQL(I).LT.MT+1.D0.AND.DBLE(MQVA(I)).LE.0.D0)
     $       WRITE(*,*)'WARNING: NEGATIVE UP-LEFT SQUARK EIGENVALUE'
        IF(QTHUR(I).LT.MT+1.D0.AND.DBLE(MUPVA(I)).LE.0.d0)
     $       WRITE(*,*)'WARNING: NEGATIVE UP-RIGHT SQUARK EIGENVALUE'
        IF(QTHDR(I).LT.MT+1.D0.AND.DBLE(MDVA(I)).LE.0.D0)
     $       WRITE(*,*)'WARNING: NEGATIVE DOWN-RIGHT SQUARK EIGENVALUE'
        IF(QTHLL(I).LT.MT+1.D0.AND.DBLE(MLVA(I)).LE.0.D0)
     $       WRITE(*,*)'WARNING: NEGATIVE L-SLEPTON EIGENVALUE'
        IF(QTHER(I).LT.MT+1.D0.AND.DBLE(MEVA(I)).LE.0.D0)
     $       WRITE(*,*)'WARNING: NEGATIVE E-SLEPTON EIGENVALUE'
      END DO
!
!Fix the number of steps in case the thresholds have changed.
!
      IF(NLTMT.LT.19)THEN
        IF(NLTMT.NE.0)THEN
          NSTEPTHRESH(NLTMT)=INT(ABS(DLOG(MT/QTHSORT(NLTMT+1)))
     $                                     *NSTEP/DLOG(MHIGH/MT)*25)
        END IF
        DO I=NLTMT+1,19
          NSTEPTHRESH(I)=INT(ABS(DLOG(QTHSORT(I)/QTHSORT(I+1)))
     $                                     *NSTEP/DLOG(MHIGH/MT)*25)
          IF(NSTEPTHRESH(I).EQ.0.AND.
     $                            ABS(QTHSORT(I)-QTHSORT(I+1)).GT.1D-10)
     $                                              NSTEPTHRESH(I)=10
        END DO
      END IF
!
!If the location of the thresholds is changing, it may be useful to the
!user to know their values. The following lines can be reinstated to
!write them to the screen, or changed to print to a file.
!
!      IF(THCNG.EQ.0)THEN
!        WRITE(*,*)
!        WRITE(*,*)'RGEMS',RGEMS
!        WRITE(*,*)'MU (16)',QNSH
!        WRITE(*,*)'GLUINO (17)',QNSG
!        WRITE(*,*)'HIGGS (18)',QNH
!        WRITE(*,*)'BINO (19)',QTHSB
!        WRITE(*,*)'WINO (20)',QTHSW
!        DO I=1,3
!        WRITE(*,*)'SQUARKS (Q,UR,DR) (1,4,7)',I,' '
!     $                                    ,QTHQL(I),QTHUR(I),QTHDR(I)
!        END DO
!        DO I=1,3
!          WRITE(*,*)'SLEPTONS (L,ER) (10,13)',I,' ',QTHLL(I),QTHER(I)
!        END DO
!      END IF
!
      IF(NHTMT.GT.LOCMH)CALL DOWNMHCOND
!
!      CLOSE(52)
!      CLOSE(53)
!      CLOSE(54)
!      CLOSE(55)
!      CLOSE(56)
!      CLOSE(57)
!      CLOSE(58)
!      CLOSE(59)
!      CLOSE(60)
!      CLOSE(61)
!      CLOSE(62)
!      CLOSE(63)
!      CLOSE(64)
!      CLOSE(65)
!      CLOSE(66)
!      CLOSE(67)
!      CLOSE(68)
!      CLOSE(69)
!      CLOSE(70)
!      CLOSE(71)
!      CLOSE(72)
!      CLOSE(73)
!      CLOSE(74)
!      CLOSE(75)
!      CLOSE(76)
!      CLOSE(77)
!      CLOSE(78)
!      CLOSE(79)
!      CLOSE(80)
!      CLOSE(81)
!      CLOSE(82)
!      CLOSE(83)
!      CLOSE(84)
!      CLOSE(85)
!      CLOSE(86)
!      CLOSE(87)
!      CLOSE(88)
!      CLOSE(89)
!      CLOSE(90)
!      CLOSE(91)
!      CLOSE(92)
!      CLOSE(93)
!      CLOSE(94)
!      CLOSE(95)
!      CLOSE(96)
!      CLOSE(97)
!      CLOSE(98)
!      CLOSE(99)
!
      RETURN
      END
