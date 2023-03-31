CDECK  ID>, QFUNC.
      SUBROUTINE QFUNC
C
C          Find approximate QMW and QTW dependence for DRELLYAN.
C          Set up /WGEN/ to generate QMW and QTW. Fit is
C          Non-resonant:
C            SIGMA=ANORM/(Q2/QMAX**2)**QPOW/(PT**2+RNU2)**PTPOW
C          Resonant:
C            SIGMA=ANORM/((Q**2-M**2)**2+M**2*GAM**2)
C          with appropriate M and GAM.
C
C          Ver. 6.23: Remove extension of region 1 under region 2
C                     to avoid discontinuity in d(sigma)/d(M)
C          Ver. 6.40: Scale Q**2 fit by QMAX**2 to avoid underflow
C                     problems. Must also change DRLLYN
C
      IMPLICIT NONE
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/DYPAR/FLW,RNU2(3),ANORM(3),QPOW(3),PTPOW(3)
      SAVE /DYPAR/
      LOGICAL FLW
      REAL      RNU2,ANORM,QPOW,PTPOW
      COMMON/DYLIM/QMIN,QMAX,QTMIN,QTMAX,YWMIN,YWMAX,XWMIN,XWMAX,THWMIN,
     2  THWMAX,PHWMIN,PHWMAX
     3  ,SETLMQ(12)
      SAVE /DYLIM/
      LOGICAL SETLMQ
      EQUIVALENCE(BLIM1(1),QMIN)
      REAL      QMIN,QMAX,QTMIN,QTMAX,YWMIN,YWMAX,XWMIN,XWMAX,THWMIN,
     +          THWMAX,PHWMIN,PHWMAX,BLIM1(12)
      COMMON/JETPAR/P(3),PT(3),YJ(3),PHI(3),XJ(3),TH(3),CTH(3),STH(3)
     1 ,JETTYP(3),SHAT,THAT,UHAT,QSQ,X1,X2,PBEAM(2)
     2 ,QMW,QW,QTW,YW,XW,THW,QTMW,PHIW,SHAT1,THAT1,UHAT1,JWTYP
     3 ,ALFQSQ,CTHW,STHW,Q0W
     4 ,INITYP(2),ISIGS,PBEAMS(5)
      SAVE /JETPAR/
      INTEGER   JETTYP,JWTYP,INITYP,ISIGS
      REAL      P,PT,YJ,PHI,XJ,TH,CTH,STH,SHAT,THAT,UHAT,QSQ,X1,X2,
     +          PBEAM,QMW,QW,QTW,YW,XW,THW,QTMW,PHIW,SHAT1,THAT1,UHAT1,
     +          ALFQSQ,CTHW,STHW,Q0W,PBEAMS
C          Jet limits
      INTEGER MXLIM
      PARAMETER (MXLIM=8)
      INTEGER MXLX12
      PARAMETER (MXLX12=12*MXLIM)
      COMMON/JETLIM/PMIN(MXLIM),PMAX(MXLIM),PTMIN(MXLIM),PTMAX(MXLIM),
     $YJMIN(MXLIM),YJMAX(MXLIM),PHIMIN(MXLIM),PHIMAX(MXLIM),
     $XJMIN(MXLIM),XJMAX(MXLIM),THMIN(MXLIM),THMAX(MXLIM),
     $SETLMJ(12*MXLIM)
      SAVE /JETLIM/
      COMMON/FIXPAR/FIXP(MXLIM),FIXPT(MXLIM),FIXYJ(MXLIM),
     $FIXPHI(MXLIM),FIXXJ(MXLIM),FIXQM,FIXQT,FIXYW,FIXXW,FIXPHW
      SAVE /FIXPAR/
      COMMON/SGNPAR/CTHS(2,MXLIM),THS(2,MXLIM),YJS(2,MXLIM),XJS(2,MXLIM)  
      SAVE /SGNPAR/
      REAL      PMIN,PMAX,PTMIN,PTMAX,YJMIN,YJMAX,PHIMIN,PHIMAX,XJMIN,
     +          XJMAX,THMIN,THMAX,BLIMS(12*MXLIM),CTHS,THS,YJS,XJS
      LOGICAL SETLMJ
      LOGICAL FIXQM,FIXQT,FIXYW,FIXXW,FIXPHW
      LOGICAL FIXP,FIXPT,FIXYJ,FIXPHI,FIXXJ
      EQUIVALENCE(BLIMS(1),PMIN(1))
      INTEGER MXGOQ,MXGOJ
      PARAMETER (MXGOQ=85,MXGOJ=8)
      COMMON/Q1Q2/GOQ(MXGOQ,MXGOJ),GOALL(MXGOJ),GODY(4),STDDY,
     $GOWW(25,2),ALLWW(2),GOWMOD(25,MXGOJ)
      SAVE /Q1Q2/
      LOGICAL GOQ,GOALL,GODY,STDDY,GOWW,ALLWW,GOWMOD
      COMMON/WCON/SIN2W,WMASS(4),WGAM(4),AQ(12,4),BQ(12,4),COUT(4),
     1MATCH(25,4),WCBR(25,4),CUTOFF,CUTPOW,TBRWW(4,2),RBRWW(12,4,2),EZ,
     2AQDP(12,4),BQDP(12,4),EZDP,WFUDGE
      SAVE /WCON/
      DOUBLE PRECISION AQDP,BQDP,EZDP
      INTEGER   MATCH
      REAL      SIN2W,WMASS,WGAM,AQ,BQ,COUT,WCBR,CUTOFF,CUTPOW,TBRWW,
     +          RBRWW,EZ,WFUDGE
      COMMON/WCON2/CUMWBR(25,3)
      REAL CUMWBR
      COMMON/WGEN/PTGN(3,3),QGEN(3,3),PTSEL(3),QSEL(3),SIGSL(3),NKL,NKH
     1,EMSQ,EMGAM,KSEL,QSELWT(3)
      SAVE /WGEN/
      INTEGER   NKL,NKH,KSEL
      REAL      PTGN,QGEN,PTSEL,QSEL,SIGSL,EMSQ,EMGAM,QSELWT
      INTEGER   MXSIGS,IOPAK
      PARAMETER (MXSIGS=3000,IOPAK=100)
      COMMON/JETSIG/SIGMA,SIGS(MXSIGS),NSIGS,INOUT(MXSIGS),SIGEVT
      SAVE /JETSIG/
      INTEGER   NSIGS,INOUT
      REAL      SIGMA,SIGS,SIGEVT
      INTEGER MXKEYS
      PARAMETER (MXKEYS=20)
      COMMON/KEYS/IKEYS,KEYON,KEYS(MXKEYS)
      COMMON/XKEYS/REAC
      SAVE /KEYS/,/XKEYS/
      LOGICAL KEYS
      LOGICAL KEYON
      CHARACTER*8 REAC
      INTEGER   IKEYS
      COMMON/HCON/ANWWWW(4,4,4),ADWWWW(2,4),AIWWWW(4)
     $,HMASS,HGAM,HGAMS(29),ETAHGG,MATCHH(29),ZSTARS(4,2)
     $,IHTYPE,HGAMSS(85,85)
      SAVE /HCON/
      DOUBLE PRECISION ANWWWW,ADWWWW,AIWWWW
      INTEGER   MATCHH,IHTYPE
      REAL      HMASS,HGAM,HGAMS,ETAHGG,ZSTARS,HGAMSS
      COMMON/TCPAR/TCMRHO,TCGRHO
      SAVE /TCPAR/
      REAL TCMRHO,TCGRHO
      COMMON/XMSSM/GOMSSM,GOSUG,GOGMSB,GOAMSB,AL3UNI,GOMMAM,GOHCAM
     $,XGLSS,XMUSS,XHASS,XTBSS
     $,XQ1SS,XDRSS,XURSS,XL1SS,XERSS
     $,XQ2SS,XSRSS,XCRSS,XL2SS,XMRSS
     $,XQ3SS,XBRSS,XTRSS,XL3SS,XTARSS,XATSS,XABSS,XATASS
     $,XM1SS,XM2SS,XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU
     $,XLAMGM,XMESGM,XN5GM,XCMGV,XMGVTO
     $,XRSLGM,XDHDGM,XDHUGM,XDYGM,XN51GM,XN52GM,XN53GM
     $,XMN3NR,XMAJNR,XANSS,XNRSS,XSBCS,
     $XCQAM,XCDAM,XCUAM,XCLAM,XCEAM,XCHDAM,XCHUAM,
     $XL1AM,XL2AM,XL3AM
      SAVE /XMSSM/
      REAL XGLSS,XMUSS,XHASS,XTBSS
     $,XQ1SS,XDRSS,XURSS,XL1SS,XERSS
     $,XQ2SS,XSRSS,XCRSS,XL2SS,XMRSS
     $,XQ3SS,XBRSS,XTRSS,XL3SS,XTARSS,XATSS,XABSS,XATASS
     $,XM1SS,XM2SS
     $,XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU
     $,XLAMGM,XMESGM,XN5GM,XCMGV,XMGVTO
     $,XRSLGM,XDHDGM,XDHUGM,XDYGM,XN51GM,XN52GM,XN53GM
     $,XMN3NR,XMAJNR,XANSS,XNRSS,XSBCS,
     $XCQAM,XCDAM,XCUAM,XCLAM,XCEAM,XCHDAM,XCHUAM,
     $XL1AM,XL2AM,XL3AM
      LOGICAL GOMSSM,GOSUG,GOGMSB,GOAMSB,AL3UNI,GOMMAM,GOHCAM
C
      REAL QT2CUT,DPT,QMN,QMX,EM,GAM,DELM,QSTOR,SUMS,DQ,ETAX,ETA,
     $Q2,XI,ALI,SIGSAV,T1,T2,T3,T4,T5,DET,DEVMAX,PTNU,ALPTNU,ALQ2,FIT,
     $DEV,DY3,DYW,SIG00,FACTOR,FAC1,C1,B1,SUM,AL1,QMAX2
      INTEGER NDIV1,NDIV2,K,I,NQS,J,N,NDIV3,NDIV4,IW,I3,II
      DIMENSION SUMS(9)
      DIMENSION QMN(3),QMX(3)
      DIMENSION SIGSAV(20,20)
C
C          QT cutoff function
      QT2CUT(QMW)=CUTOFF*QMW**CUTPOW
C
C          Entry
C
      IF(FIXQM) THEN
        NDIV1=1
      ELSE
        NDIV1=20
      ENDIF
      IF(FIXQT) THEN
        NDIV2=1
      ELSE
        NDIV2=20
      ENDIF
C
      DPT=(PTMAX(3)-PTMIN(3))/NDIV2
      YJ(3)=0
      YW=0.
      CTH(3)=0.
      STH(3)=1.
      IF(GODY(4)) JWTYP=4
      NKL=1
      NKH=1
      QMN(1)=QMIN
      QMX(1)=QMAX
      QMAX2=QMAX**2
C
C          Define resonance region
C
      IF(KEYS(3)) THEN
        IF(JWTYP.EQ.1) GO TO 99
        EM=WMASS(JWTYP)
        GAM=WGAM(JWTYP)
        DELM=20.
      ELSEIF(KEYS(7)) THEN
        EM=HMASS
        GAM=HGAM
        DELM=.201357*EM
        DELM=AMIN1(DELM,1.5*HGAM)
        DELM=AMAX1(DELM,.1*EM)
      ELSEIF(KEYS(9)) THEN
        EM=TCMRHO
        GAM=TCGRHO
        DELM=.201357*EM
        DELM=AMIN1(DELM,1.5*TCGRHO)
        DELM=AMAX1(DELM,.1*EM)
C          No resonance region for KKG
      ELSEIF(KEYS(11)) THEN
        EM=QMAX
        GAM=0.
        DELM=0.
      ENDIF
      EMGAM=EM*GAM
      EMSQ=EM**2
C          Region limits
      QMN(2)=EM-DELM
      QMN(3)=EM+DELM
      QMX(1)=QMN(2)
      QMX(2)=QMN(3)
      NKL=1
      NKH=3
      IF(QMAX.LE.QMN(3)) NKH=2
      IF(QMAX.LE.QMN(2)) NKH=1
      IF(QMIN.GE.QMN(2)) NKL=2
      IF(QMIN.GE.QMN(3)) NKL=3
      QMX(NKH)=QMAX
      QMN(NKL)=QMIN
   99 CONTINUE
C
C            Fit over regions NKL to NKH
C            Region 1 is below resonance
C            Region 2 is inside resonance
C            Region 3 is above resonance
C            FIT=ANORM/(Q2/QMAX**2)**QPOW/(PT**2+RNU2)**PTPOW
C
      DO 100 K=1,3
        ANORM(K)=0.
        PTPOW(K)=0.
        QPOW(K)=0.
        RNU2(K)=QT2CUT(QMIN)
100   CONTINUE
C
C          Loop over regions
C
      DO 200 K=NKL,NKH
        DO 210 I=1,9
210     SUMS(I)=0
        DQ=(QMX(K)-QMN(K))/NDIV1
        NQS=NDIV1
        DO 220 I=1,NDIV2
          PT(3)=PTMIN(3)+(I-1)*DPT
          QTW=PT(3)
          P(3)=PT(3)
          RNU2(K)=QT2CUT(QMN(K))
          ETAX=PT(3)**2+RNU2(K)
          ETA=ALOG(ETAX)
          DO 230 J=1,NQS
            QMW=QMN(K)+(J-1)*DQ
            Q2=QMW*QMW
            XI=ALOG(Q2/QMAX2)
            SUMS(1)=SUMS(1)+XI
            SUMS(2)=SUMS(2)+ETA
            SUMS(5)=SUMS(5)+ETA*ETA
            SUMS(4)=SUMS(4)+XI**2
            SUMS(7)=SUMS(7)+XI*ETA
C            Cross section
            IF(KEYS(3)) THEN
              CALL SIGDY
            ELSEIF(KEYS(7).AND..NOT.GOMSSM) THEN
              CALL SIGH
            ELSEIF(KEYS(7).AND.GOMSSM) THEN
              CALL SIGHSS
            ELSEIF(KEYS(9)) THEN
              CALL SIGTC
            ELSEIF(KEYS(11)) THEN
              CALL SIGKKG
            ENDIF
            IF(SIGMA.EQ.0.) GO TO 999
            AL1=ALOG(SIGMA)
            SIGSAV(I,J)=AL1
            IF(K.EQ.2) AL1=AL1+ALOG((Q2-EM**2)**2+EMGAM**2)
            SUMS(3)=SUMS(3)+AL1
            SUMS(8)=SUMS(8)+AL1*XI
            SUMS(9)=SUMS(9)+AL1*ETA
230         CONTINUE
220      CONTINUE
C
C            Find coefficients minimizing chisq
C
        N=NQS*NDIV2
        T1=N*SUMS(7)-SUMS(1)*SUMS(2)
        T2=N*SUMS(5)-SUMS(2)**2
        T3=N*SUMS(4)-SUMS(1)**2
        T4=N*SUMS(8)-SUMS(1)*SUMS(3)
        T5=N*SUMS(9)-SUMS(2)*SUMS(3)
        IF((FIXQM.OR.K.EQ.2).AND.FIXQT) THEN
          PTPOW(K)=0.
          QPOW(K)=0.
        ELSEIF(FIXQT) THEN
          PTPOW(K)=0.
          QPOW(K)=-T4/T3
        ELSEIF(FIXQM.OR.K.EQ.2) THEN
          PTPOW(K)=-T5/T2
          QPOW(K)=0.
        ELSE
          DET=T1**2-T2*T3
          PTPOW(K)=(T5*T3-T4*T1)/DET
          QPOW(K)=(T4*T2-T1*T5)/DET
        ENDIF
        ANORM(K)=(QPOW(K)*SUMS(1)+PTPOW(K)*SUMS(2)+SUMS(3))/N
C
C          Shift fit to obtain envelope for SIGDY
C
        DEVMAX=0.
        DO 240 I=1,NDIV2
          PT(3)=PTMIN(3)+(I-1)*DPT
          PTNU=PT(3)**2+RNU2(K)
          DO 250 J=1,NDIV1
            QMW=QMN(K)+(J-1)*DQ
            Q2=QMW**2
            ALPTNU=ALOG(PTNU)
            ALQ2=ALOG(Q2/QMAX2)
            IF(K.EQ.2) THEN
              FIT=EXP(ANORM(K)-PTPOW(K)*ALPTNU
     $        -ALOG((Q2-EM**2)**2+EMGAM**2))
            ELSE
              FIT=EXP(ANORM(K)-PTPOW(K)*ALPTNU-QPOW(K)*ALQ2)
            ENDIF
            DEV=SIGSAV(I,J)-ALOG(FIT)
            IF(DEV.GT.DEVMAX) DEVMAX=DEV
250       CONTINUE
240     CONTINUE
        ANORM(K)=ANORM(K)+DEVMAX
200   CONTINUE
C
C          Shift fit to obtain envelope in YW
      NDIV3=20
      IF(STDDY) THEN
        NDIV4=1
        DY3=0.
      ELSE
        NDIV4=20
        DY3=(YJMAX(3)-YJMIN(3))/(NDIV4-1)
      ENDIF
      DYW=(YWMAX-YWMIN)/(NDIV3-1)
C
      DO 300 K=NKL,NKH
        QMW=QMN(K)
        Q2=QMW**2
        QTW=QTMIN
        PT(3)=QTW
        P(3)=PT(3)
        YW=0.
        YJ(3)=0.
        CTH(3)=0.
        STH(3)=1.
        IF(KEYS(3)) THEN
          CALL SIGDY
        ELSEIF(KEYS(7).AND..NOT.GOMSSM) THEN
          CALL SIGH
        ELSEIF(KEYS(7).AND.GOMSSM) THEN
          CALL SIGHSS
        ELSEIF(KEYS(9)) THEN
          CALL SIGTC
        ELSEIF(KEYS(11)) THEN
          CALL SIGKKG
        ENDIF
        SIG00=SIGMA
        FACTOR=1.
        DO 310 IW=1,NDIV3
          YW=YWMIN+(IW-1)*DYW
          DO 320 I3=1,NDIV4
            IF(.NOT.STDDY) THEN
              YJ(3)=YJMIN(3)+(I3-1)*DY3
              CTH(3)=TANH(YJ(3))
              STH(3)=SQRT(1.-CTH(3)**2)
              IF(STH(3).EQ.0.) GO TO 320
              TH(3)=ACOS(CTH(3))
              P(3)=PT(3)/STH(3)
            ENDIF
            IF(KEYS(3)) THEN
              CALL SIGDY
            ELSEIF(KEYS(7).AND..NOT.GOMSSM) THEN
              CALL SIGH
            ELSEIF(KEYS(7).AND.GOMSSM) THEN
              CALL SIGHSS
            ELSEIF(KEYS(9)) THEN
              CALL SIGTC
            ELSEIF(KEYS(11)) THEN
              CALL SIGKKG
            ENDIF
            FAC1=SIGMA/SIG00
            FACTOR=AMAX1(FACTOR,FAC1)
320       CONTINUE
310     CONTINUE
        ANORM(K)=ALOG(FACTOR)+ANORM(K)
300   CONTINUE
C
C          Set up generating constants for PT**2 and QMW**2
C
      DO 400 K=NKL,NKH
        C1=1.-PTPOW(K)
        PTGN(1,K)=(PTMIN(3)**2+RNU2(K))**C1
        PTGN(2,K)=(PTMAX(3)**2+RNU2(K))**C1-PTGN(1,K)
        PTGN(3,K)=1./C1
        IF(K.EQ.2) THEN
          QGEN(1,2)=ATAN((QMN(2)**2-EMSQ)/EMGAM)
          QGEN(2,2)=ATAN((QMX(2)**2-EMSQ)/EMGAM)-QGEN(1,2)
          QGEN(3,2)=EMGAM
        ELSE
          B1=1.-QPOW(K)
          QGEN(1,K)=(QMN(K)/QMAX)**(2.*B1)
          QGEN(2,K)=(QMX(K)/QMAX)**(2.*B1)-QGEN(1,K)
          QGEN(3,K)=1./B1
        ENDIF
400   CONTINUE
C
      DO 410 K=1,3
410   QSELWT(K)=0.
      SUM=0.
C
      DO 420 K=NKL,NKH
        QSELWT(K)=1.
        IF(.NOT.FIXQT) QSELWT(K)=QSELWT(K)*PTGN(2,K)*PTGN(3,K)
        IF(.NOT.FIXQM) THEN
          IF(K.EQ.2) THEN
            QSELWT(K)=QSELWT(K)*QGEN(2,K)/EMGAM
          ELSE
            QSELWT(K)=QMAX**2*QSELWT(K)*QGEN(2,K)*QGEN(3,K)
          ENDIF
        ENDIF
        QSELWT(K)=EXP(ALOG(QSELWT(K))+ANORM(K))
        SUM=SUM+QSELWT(K)
420   CONTINUE
C
      DO 430 K=1,3
        QSELWT(K)=QSELWT(K)/SUM
430   CONTINUE
C
C          Write fit to output
C
      WRITE(ITLIS,4301)
4301  FORMAT(//10X,' QT AND Q FIRST GENERATED BY--'/)
      DO 440 K=NKL,NKH
        WRITE(ITLIS,4402) K,QMN(K),QMX(K)
4402    FORMAT(//5X,' REGION',I2,5X,E11.4,' < Q < ',E11.5)
        WRITE(ITLIS,4403) (PTGN(II,K),II=1,3),RNU2(K)
4403    FORMAT(/' QT**2 = (',E11.4,' + ',E11.4,' * RANF) ** ',E11.4,
     $  ' - ',E11.4)
        IF(K.NE.2) THEN
          WRITE(ITLIS,4404) QMAX2,(QGEN(II,K),II=1,3)
4404      FORMAT(/' Q**2  = ',E11.4,' * (',E11.4,' + ',E11.4,
     $    ' * RANF) ** ',E11.4)
        ELSE
          WRITE(ITLIS,4505) QGEN(3,K),QGEN(1,K),QGEN(2,K),EMSQ
4505       FORMAT(/' Q**2  = ',E11.4,' * TAN(',E11.4,' + ',E11.4,
     $    ' * RANF) + ',E11.4)
        ENDIF
        WRITE(ITLIS,4506) QSELWT(K)
4506    FORMAT(/' WEIGHT = ',E11.4)
440   CONTINUE
C
C          Set fixed limits if any
C
      IF(FIXQT) THEN
        PTMAX(3)=PTMIN(3)
        PT(3)=PTMIN(3)
        QTW=PT(3)
      ENDIF
      IF(FIXQM) THEN
        QMAX=QMIN
        QMW=QMIN
      ENDIF
      RETURN
C
C          Fit fails if SIGMA=0 in allowed range
C
999   WRITE(ITLIS,9990) QMW,QTW
9990  FORMAT(//' ERROR IN QFUNC...SIGMA=0 FOR QMW = ',E12.4,' , QTW = ',
     1E12.4/' CHECK YOUR LIMITS')
      STOP 99
      END
