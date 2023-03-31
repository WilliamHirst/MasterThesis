!
      SUBROUTINE UPMZMHIGH2(NSTEP)
!
!Purpose: To run back up to M_GUT as part of the iterative running.
!         Although different in name to distinguish from UPMZMHIGH,
!         this subroutine does almost exactly the same thing but
!         starts with saved values of GSM at M_T
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
      COMMON/WKYUK/LAMTMT,LAMBMZ,LAMTAMZ
      DOUBLE PRECISION LAMTMT,LAMBMZ,LAMTAMZ
      SAVE/WKYUK/
!
      COMMON/ATMZ/G1MZ,G2MZ,G3MZ,VSMMZ,LAMBDAMZ,LAMTMZ
      DOUBLE COMPLEX G1MZ,G2MZ,G3MZ,VSMMZ,LAMBDAMZ,LAMTMZ
      SAVE/ATMZ/
!
      COMMON/SMRGE/SMRGEMH,SMQSTEP,NU,SMDR2LP
      DOUBLE PRECISION SMRGEMH,SMQSTEP
      INTEGER NU,SMDR2LP
      SAVE/SMRGE/
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
      COMMON/LOOPS/SSQSTEP,SW2LP
      DOUBLE PRECISION SSQSTEP
      INTEGER SW2LP
      SAVE/LOOPS/
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
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
      DOUBLE COMPLEX W(1803),GSM(32),WSM(96),CID(3,3)
      DOUBLE PRECISION DW(1803),DGSM(32),DWSM(96)
      DOUBLE PRECISION PI,TMT,TZ,TTH(20),THIGH,DT,T,A1I,A2I,Q
      INTEGER NSTEP,NSTEPSM,I,J,II,LOOPNSTEP,NSTEPMT,BELOWMS
      EXTERNAL CRGE601,DRGE601,CSMRGEDR,DSMRGEDR
!
      DATA CID(1,1)/(1.D0,0.D0)/,CID(1,2)/(0.D0,0.D0)/
     $    ,CID(1,3)/(0.D0,0.D0)/
      DATA CID(2,1)/(0.D0,0.D0)/,CID(2,2)/(1.D0,0.D0)/
     $    ,CID(2,3)/(0.D0,0.D0)/
      DATA CID(3,1)/(0.D0,0.D0)/,CID(3,2)/(0.D0,0.D0)/
     $    ,CID(3,3)/(1.D0,0.D0)/
!
!Files for outcoup
!
!      OPEN(71,FILE='out/uau.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(72,FILE='out/uad.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(73,FILE='out/uae.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(74,FILE='out/umhumu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(75,FILE='out/umhdmu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(76,FILE='out/umq.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(77,FILE='out/umup.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(78,FILE='out/umd.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(79,FILE='out/uml.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(80,FILE='out/ume.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(81,FILE='out/umu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(82,FILE='out/ub.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(83,FILE='out/umtsfu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(84,FILE='out/umtsfd.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(85,FILE='out/umtsfe.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(86,FILE='out/utriu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(87,FILE='out/utrid.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(88,FILE='out/utrie.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(89,FILE='out/umhud.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(90,FILE='out/uaum.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(91,FILE='out/uadm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(92,FILE='out/uaem.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(93,FILE='out/umhumum.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(94,FILE='out/umhdmum.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(95,FILE='out/umqm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(96,FILE='out/umupm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(97,FILE='out/umdm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(98,FILE='out/umlm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(99,FILE='out/umem.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(70,FILE='out/umum.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(69,FILE='out/ubm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(68,FILE='out/ug1.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(67,FILE='out/ug2.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(66,FILE='out/ug3.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(65,FILE='out/uytau.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(64,FILE='out/uyb.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(63,FILE='out/uyu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(62,FILE='out/um1.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(61,FILE='out/um2.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(60,FILE='out/um3.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(59,FILE='out/ugtpq.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(58,FILE='out/ugtpu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(57,FILE='out/ugtq.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(56,FILE='out/uftuq.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(55,FILE='out/dvu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(54,FILE='out/dvd.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(53,FILE='out/dmuflavcomp.dat',STATUS='UNKNOWN'
!     $ ,FORM='FORMATTED')
!      OPEN(52,FILE='out/dm12comp.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!
      DO I=1,32
        GSM(I)=(0.D0,0.D0)
        DGSM(I)=0.D0
      END DO
!
      BELOWMS=1
      PI=4.D0*DATAN(1.D0)
!
!Set NSTEP to a larger number than the downward running so that I
!do not go past the unification point with MHIGH(UP)>MHIGH(DOWN)
!Fixed so that spacing (DeltaT) between points is constant.
!
      IF(UNI.EQ.1)THEN
        IF(NLTMT.NE.20)THEN
          NSTEP=NSTEP*LOG(QTHSORT(20)/1.D19)/LOG(QTHSORT(20)/MHIGH)
        ELSE
          NSTEP=NSTEP*LOG(MT/1.D19)/LOG(MT/MHIGH)
        END IF
      END IF
!
!The next line sets the upper scale so high that the couplings will
!unify before the limit is reached
!
      IF(UNI.EQ.1)MHIGH=1.D19
      TZ=LOG(MZ/MHIGH)
      TMT=LOG(MT/MHIGH)
      DO I=1,20
        TTH(I)=LOG(QTHSORT(I)/MHIGH)
      END DO
      THIGH=0.D0
!
!Run from M_Z to M_t
!
      NSTEPSM=50
      DT=(TMT-TZ)/FLOAT(NSTEPSM)
      NU=3 !top decouples at M_Z
      GSM(1)=G1MZ
      GSM(2)=G2MZ
      GSM(3)=G3MZ
      GSM(4)=DCMPLX(MWEAK(1))/VSMMZ
      GSM(8)=DCMPLX(MWEAK(2))/VSMMZ
      GSM(12)=LAMTMZ
      GSM(13)=DCMPLX(MWEAK(3))/VSMMZ
      GSM(17)=DCMPLX(MWEAK(4))/VSMMZ
      GSM(21)=DCMPLX(LAMBMZ)
      GSM(22)=DCMPLX(MWEAK(5))/VSMMZ
      GSM(26)=DCMPLX(MWEAK(6))/VSMMZ
      GSM(30)=DCMPLX(LAMTAMZ)
      GSM(31)=LAMBDAMZ
      GSM(32)=VSMMZ
!
      IF(COMP.EQ.0)THEN
        DO I=1,32
          DGSM(I)=DBLE(GSM(I))
        END DO
      END IF
!
      DO II=1,NSTEPSM
        T=TZ+(TMT-TZ)*FLOAT(II-1)/FLOAT(NSTEPSM)
        SMQSTEP=MHIGH*EXP(T)
        IF(II.EQ.1)EPS=ABS(SMQSTEP*(EXP(DT)-1)/(6.D0*PI))
        IF(COMP.EQ.0)THEN
          CALL DRKSTP(32,DT,T,DGSM,DSMRGEDR,DWSM)
        ELSE
          CALL CRKSTP(32,DT,T,GSM,CSMRGEDR,WSM)
        END IF
        EPS=-ABS(SMQSTEP*(EXP(DT)-1)/(6.D0*PI))
      END DO
!
      IF(COMP.EQ.0)THEN
        DO I=1,32
          GSM(I)=DCMPLX(DGSM(I))
        END DO
      END IF
!
!Now introduce the top yukawa and rotate
!
      GSM(12)=DCMPLX(LAMTMT)
      CALL ROTATESM(GSM)
!
!Place SM couplings at M_T into G
!
      DO I=4,30
        IF(I.LT.7)G(I-3)=GSM(I-3)
        G(I+108)=GSM(I)
      END DO
      G(428)=GSM(32)
      G(429)=GSM(31)
!
!Now run the rest of the way past the thresholds with LAMBDA_t
!
      IF(NLTMT.GE.LOCMH)CALL UPMHCOND2
!
!First threshold
!
      IF(NLTMT.NE.20)THEN
        IF(NLTMT.NE.0)THEN
          DT=(TTH(NLTMT+1)-TMT)/FLOAT(NSTEPTHRESH(NLTMT))
          LOOPNSTEP=NSTEPTHRESH(NLTMT)
        ELSE
          NSTEPMT=INT(ABS(DLOG(MT/QTHSORT(1)))*NSTEP/DLOG(MHIGH/MT)*25)
          DT=(TTH(NLTMT+1)-TMT)/FLOAT(NSTEPMT)
          LOOPNSTEP=NSTEPMT
        END IF
!
        DO II=1,LOOPNSTEP
          T=TMT+(TTH(NLTMT+1)-TMT)*FLOAT(II-1)/FLOAT(LOOPNSTEP)
          SSQSTEP=MHIGH*EXP(T)
!
          IF(BELOWMS.EQ.1.AND.MHIGH*EXP(T+DT).GE.RGEMS)THEN
            BELOWMS=0
            CALL ROTBACK(1)
            G(12)=G(12)/(1.D0-DBLE(RTISA))
            G(21)=G(21)/(1.D0-DBLE(RBISA))
            G(30)=G(30)/(1.D0-DBLE(RLISA))
            IF(RGEMS.LE.QNH)THEN
              G(120)=G(120)/(1.D0-DBLE(RTISA))
              G(129)=G(129)/(1.D0-DBLE(RBISA))
              G(138)=G(138)/(1.D0-DBLE(RLISA))
            END IF
            CALL ROTATE(1)
          END IF
!
          IF(II.EQ.1)EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
!
!Check the sfermion thresholds
!
          CALL CHINT(SSQSTEP,G)
!          CALL OUTCOUP(MHIGH*EXP(T),1)
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
          EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
        END DO
!
!Matching conditions at M_H
!
        IF(NLTMT+1.EQ.LOCMH)CALL UPMHCOND2
!
        DO I=NLTMT+2,20
!
!Don't run between degenerate thresholds
!
          IF(NSTEPTHRESH(I-1).EQ.0)GOTO 20
!
          DT=(TTH(I)-TTH(I-1))/FLOAT(NSTEPTHRESH(I-1))
!          WRITE(*,*)'RUNNING BETWEEN',QTHSORT(I-1),QTHSORT(I)
!
          DO II=1,NSTEPTHRESH(I-1)
            T=TTH(I-1)+(TTH(I)-TTH(I-1))
     $                         *FLOAT(II-1)/FLOAT(NSTEPTHRESH(I-1))
            SSQSTEP=MHIGH*EXP(T)
!
            IF(BELOWMS.EQ.1.AND.MHIGH*EXP(T+DT).GE.RGEMS)THEN
              BELOWMS=0
              CALL ROTBACK(1)
              G(12)=G(12)/(1.D0-DBLE(RTISA))
              G(21)=G(21)/(1.D0-DBLE(RBISA))
              G(30)=G(30)/(1.D0-DBLE(RLISA))
              IF(RGEMS.LE.QNH)THEN
                G(120)=G(120)/(1.D0-DBLE(RTISA))
                G(129)=G(129)/(1.D0-DBLE(RBISA))
                G(138)=G(138)/(1.D0-DBLE(RLISA))
              END IF
              CALL ROTATE(1)
            END IF
!
            IF(II.EQ.1)EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
            CALL CHINT(SSQSTEP,G)
!            CALL OUTCOUP(MHIGH*EXP(T),1)
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
            EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
          END DO
!
   20     IF(LOCMH.EQ.I)CALL UPMHCOND2
        END DO
      END IF
!
!Set the g-tilde terms equal to ( the usual g x identity )
!This prevents numerical accidents which may cause spurious
!results.
!
      DO I=1,3
        DO J=1,3
          G(138+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(147+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(156+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(165+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(174+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(185+(I-1)*3+J)=G(2)*CID(I,J)
          G(194+(I-1)*3+J)=G(2)*CID(I,J)
          G(205+(I-1)*3+J)=G(3)*CID(I,J)
          G(214+(I-1)*3+J)=G(3)*CID(I,J)
          G(223+(I-1)*3+J)=G(3)*CID(I,J)
          G(429+(I-1)*3+J)=G(429+(I-1)*3+J)*DSQRT(
     $                     ABS(G(61)+G(62)-G(351)-G(352))
     $                     /(2.D0*ABS(G(108)**2)))
          G(438+(I-1)*3+J)=G(438+(I-1)*3+J)*DSQRT(
     $                     ABS(G(61)+G(62)-G(351)-G(352))
     $                     /(2.D0*ABS(G(108)**2)))
          G(447+(I-1)*3+J)=G(447+(I-1)*3+J)*DSQRT(
     $                     ABS(G(61)+G(62)-G(351)-G(352))
     $                     /(2.D0*ABS(G(108)**2)))
        END DO
      END DO
      IF(ABS(G(61)+G(62)-G(351)-G(352)).GT.0.D0)THEN
        G(108)=G(108)*DSQRT(ABS(G(61)+G(62)-G(351)-G(352))
     $               /(2.D0*ABS(G(108)**2)))
      ELSE
!
!Fix mu to Mz if its new squared value would have been negative
!
        G(108)=G(108)*MZ/ABS(G(108))
      END IF
!
!Continue running to m_high. I need to be careful about the case that
!all thresholds are below m_t.
!
      IF(NLTMT.NE.20)THEN
        DT=(THIGH-TTH(20))/FLOAT(NSTEP)
      ELSE
        DT=(THIGH-TMT)/FLOAT(NSTEP)
      END IF
!
      DO II=1,NSTEP
        IF(NLTMT.NE.20)THEN
          T=TTH(20)+(THIGH-TTH(20))*FLOAT(II-1)/FLOAT(NSTEP)
        ELSE
          T=TMT+(THIGH-TMT)*FLOAT(II-1)/FLOAT(NSTEP)
        END IF
        SSQSTEP=MHIGH*EXP(T)
        Q=MHIGH*EXP(T)
        IF(COMP.EQ.0)THEN
          A1I=4.D0*PI/DG(1)**2
          A2I=4.D0*PI/DG(2)**2
        ELSE
          A1I=4.D0*PI/DBLE(G(1)**2)
          A2I=4.D0*PI/DBLE(G(2)**2)
        END IF
        IF(A1I.LT.A2I.AND.UNI.EQ.1)THEN
          MHIGH=Q
          GO TO 30
        END IF
!
        IF(BELOWMS.EQ.1.AND.MHIGH*EXP(T+DT).GE.RGEMS)THEN
          BELOWMS=0
          CALL ROTBACK(1)
          G(12)=G(12)/(1.D0-DBLE(RTISA))
          G(21)=G(21)/(1.D0-DBLE(RBISA))
          G(30)=G(30)/(1.D0-DBLE(RLISA))
          IF(RGEMS.LE.QNH)THEN
            G(120)=G(120)/(1.D0-DBLE(RTISA))
            G(129)=G(129)/(1.D0-DBLE(RBISA))
            G(138)=G(138)/(1.D0-DBLE(RLISA))
          END IF
          CALL ROTATE(1)
        END IF
!
        IF(II.EQ.1)EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
        CALL CHINT(SSQSTEP,G)
!        CALL OUTCOUP(MHIGH*EXP(T),1)
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
        EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
      END DO
!
!NB: If this message is written, the programme will probably
!    not succeed. It may be better to put a stop statement here.
!
      IF(UNI.EQ.1)WRITE(*,*)'ERROR: UNIFICATION NOT FOUND', 
     $' IN ITERATIVE SECTION'
      MHIGH=EXP(T)*MHIGH
 30   CONTINUE
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
