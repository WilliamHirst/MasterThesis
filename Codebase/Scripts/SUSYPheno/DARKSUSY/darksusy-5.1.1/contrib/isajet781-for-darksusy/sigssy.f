CDECK  ID>, SIGSSY.
      SUBROUTINE SIGSSY
C
C          Calculate d(sigma)/d(pt**2)d(y1)d(y2) for supersymmetric
C          particle pairs, including gluinos, gauginos, and squarks.
C
C          SIGMA    = cross section summed over types allowed by
C                     JETTYPE cards (with natural equivalence.)
C          SIGS(I)  = partial cross section for I1 + I2 --> I3 + I4
C          INOUT(I) = IOPAK**3*I4 + IOPAK**2*I3 + IOPAK*I2 +I1
C
C          Extra factor of 1/2 needed for nonidentical final jets.
C          Y=-log(tan(theta/2)) gives jacobean P1*P2/E1*E2
C
C          Dec. 1992: Use cross sections from Baer and Tata, Phys. 
C          Lett. 160B, 159; Phys. Rev. D42, 2259. These papers
C          separate L and R squarks.
C
C          Gauginos are included only for MSSM. The cross sections are
C          calculated in SIGSSZ, which is called from here.
C
      IMPLICIT NONE
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/QCDPAR/ALAM,ALAM2,CUTJET,ISTRUC
      SAVE /QCDPAR/
      INTEGER   ISTRUC
      REAL      ALAM,ALAM2,CUTJET
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
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
      INTEGER MXGOQ,MXGOJ
      PARAMETER (MXGOQ=85,MXGOJ=8)
      COMMON/Q1Q2/GOQ(MXGOQ,MXGOJ),GOALL(MXGOJ),GODY(4),STDDY,
     $GOWW(25,2),ALLWW(2),GOWMOD(25,MXGOJ)
      SAVE /Q1Q2/
      LOGICAL GOQ,GOALL,GODY,STDDY,GOWW,ALLWW,GOWMOD
      INTEGER   MXSIGS,IOPAK
      PARAMETER (MXSIGS=3000,IOPAK=100)
      COMMON/JETSIG/SIGMA,SIGS(MXSIGS),NSIGS,INOUT(MXSIGS),SIGEVT
      SAVE /JETSIG/
      INTEGER   NSIGS,INOUT
      REAL      SIGMA,SIGS,SIGEVT
      COMMON/CONST/PI,SQRT2,ALFA,GF,UNITS
      SAVE /CONST/
      REAL      PI,SQRT2,ALFA,GF,UNITS
      COMMON/QSAVE/QSAVE(29,2)
      SAVE /QSAVE/
      REAL      QSAVE
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
C          SM ident code definitions. These are standard ISAJET but
C          can be changed.
      INTEGER IDUP,IDDN,IDST,IDCH,IDBT,IDTP
      INTEGER IDNE,IDE,IDNM,IDMU,IDNT,IDTAU
      INTEGER IDGL,IDGM,IDW,IDZ,IDH
      PARAMETER (IDUP=1,IDDN=2,IDST=3,IDCH=4,IDBT=5,IDTP=6)
      PARAMETER (IDNE=11,IDE=12,IDNM=13,IDMU=14,IDNT=15,IDTAU=16)
      PARAMETER (IDGL=9,IDGM=10,IDW=80,IDZ=90,IDH=81)
C          SUSY ident code definitions. They are chosen to be similar
C          to those in versions < 6.50 but may be changed.
      INTEGER ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1
      INTEGER ISNEL,ISEL,ISNML,ISMUL,ISNTL,ISTAU1
      INTEGER ISUPR,ISDNR,ISSTR,ISCHR,ISBT2,ISTP2
      INTEGER ISNER,ISER,ISNMR,ISMUR,ISNTR,ISTAU2
      INTEGER ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,ISGL
      INTEGER ISHL,ISHH,ISHA,ISHC
      INTEGER ISGRAV
      INTEGER IDTAUL,IDTAUR
      PARAMETER (ISUPL=21,ISDNL=22,ISSTL=23,ISCHL=24,ISBT1=25,ISTP1=26)
      PARAMETER (ISNEL=31,ISEL=32,ISNML=33,ISMUL=34,ISNTL=35,ISTAU1=36)
      PARAMETER (ISUPR=41,ISDNR=42,ISSTR=43,ISCHR=44,ISBT2=45,ISTP2=46)
      PARAMETER (ISNER=51,ISER=52,ISNMR=53,ISMUR=54,ISNTR=55,ISTAU2=56)
      PARAMETER (ISGL=29)
      PARAMETER (ISZ1=30,ISZ2=40,ISZ3=50,ISZ4=60,ISW1=39,ISW2=49)
      PARAMETER (ISHL=82,ISHH=83,ISHA=84,ISHC=86)
      PARAMETER (ISGRAV=91)
      PARAMETER (IDTAUL=10016,IDTAUR=20016)
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
      REAL X(2)
      INTEGER IDQ(13),IDQSS(25),JS2JT(25)
      EQUIVALENCE (X(1),X1)
      LOGICAL LLRR
      REAL QFCN,STRUC,AMASS,FQG
      REAL AMG,SIG0,SIGR,AM1,SIG,FAC,AMQ,AM,AM2,AMQ2,S,T,U,AMG2,E1,E2,
     $AMSQ,AM1SQ,AM2SQ,SIGL
      INTEGER IFL1,IFL2,IQ1,IQ2,JQ1,JQ2,I,IFLQ1,IFLQ2,IH,IQ,     
     $JQ,JQIN1,JQIN2
C
C          IDENT codes from /SSTYPE/. (Fortran 77 allows - signs in
C          parameter statements but not data statements.)
      INTEGER MSUPL,MSDNL,MSSTL,MSCHL,MSBT1,MSTP1,
     $MSUPR,MSDNR,MSSTR,MSCHR,MSBT2,MSTP2,
     $MDUP,MDDN,MDST,MDCH,MDBT,MDTP
      PARAMETER (MSUPL=-ISUPL)
      PARAMETER (MSDNL=-ISDNL)
      PARAMETER (MSSTL=-ISSTL)
      PARAMETER (MSCHL=-ISCHL)
      PARAMETER (MSBT1=-ISBT1)
      PARAMETER (MSTP1=-ISTP1)
      PARAMETER (MSUPR=-ISUPR)
      PARAMETER (MSDNR=-ISDNR)
      PARAMETER (MSSTR=-ISSTR)
      PARAMETER (MSCHR=-ISCHR)
      PARAMETER (MSBT2=-ISBT2)
      PARAMETER (MSTP2=-ISTP2)
      PARAMETER (MDUP=-IDUP)
      PARAMETER (MDDN=-IDDN)
      PARAMETER (MDST=-IDST)
      PARAMETER (MDCH=-IDCH)
      PARAMETER (MDBT=-IDBT)
      PARAMETER (MDTP=-IDTP)
      DATA IDQSS/0,
     $ISUPL,MSUPL,ISDNL,MSDNL,ISSTL,MSSTL,ISCHL,MSCHL,ISBT1,MSBT1,
     $ISTP1,MSTP1,
     $ISUPR,MSUPR,ISDNR,MSDNR,ISSTR,MSSTR,ISCHR,MSCHR,ISBT2,MSBT2,
     $ISTP2,MSTP2/
      DATA IDQ/IDGL,IDUP,MDUP,IDDN,MDDN,IDST,MDST,IDCH,MDCH,
     $IDBT,MDBT,IDTP,MDTP/
C          JS2JT: Susy jettype -> normal jettype
      DATA JS2JT/1,
     $2,3,4,5,6,7,8,9,10,11,12,13,2,3,4,5,6,7,8,9,10,11,12,13/
C
C          Functions
      QFCN(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
      FQG(S,T,U)=((16./3.)*(1./(U*T)**2+1./(S*U)**2)
     $+2.*(-2./3.)/(S*T*U**2))*(-U*S*T**2+2.*U*S*T*(AMG2-AMQ2)
     $-2.*U*S*(AMG2-AMQ2)**2-2.*S**2*AMG2*(AMG2-AMQ2))
C
C          Initialize
C
      SIGMA=0.
      NSIGS=0
      DO 100 I=1,MXSIGS
        SIGS(I)=0.
100   CONTINUE
C
C          Gluino + gluino
C
      IF(.NOT.(GOQ(1,1).AND.GOQ(1,2))) GO TO 300
      AM=AMASS(ISGL)
      CALL TWOKIN(0.,0.,AM,AM)
      IF(X1.GE.1..OR.X2.GE.1.) GO TO 300
      AM2=AM**2
      S=SHAT
      T=THAT
      U=UHAT
      E1=SQRT(P(1)**2+AM2)
      E2=SQRT(P(2)**2+AM2)
      FAC=PI*ALFQSQ**2/S**2
      FAC=FAC*(S/SCM)*(P(1)*P(2)/(E1*E2))*UNITS
C
C          gl gl ---> glss glss
      SIG=9./4.*(2.*(T-AM2)*(U-AM2)/S**2
     $+((T-AM2)*(U-AM2)-2.*AM2*(T+AM2))/(T-AM2)**2
     $+((U-AM2)*(T-AM2)-2.*AM2*(U+AM2))/(U-AM2)**2
     $+((T-AM2)*(U-AM2)+AM2*(U-T))/(S*(T-AM2))
     $+((U-AM2)*(T-AM2)+AM2*(T-U))/(S*(U-AM2))
     $+AM2*(S-4*AM2)/((T-AM2)*(U-AM2)))
      SIG=.5*FAC*SIG*QFCN(1,1)*QFCN(1,2)
      CALL SIGFIL(SIG,1,1,1,1)
C
C          qk qb ---> glss glss
      DO 220 IQ=1,5
        IQ1=2*IQ
        IQ2=IQ1+1
C          Left squark exchange
        AMQ=AMASS(IDQSS(IQ1))
        AMQ2=AMQ**2
        SIGL=(8./3.)*((T-AM2)**2+(U-AM2)**2+2.*AM2*S)/(S**2)
     $  +(32./27.)*(T-AM2)**2/(T-AMQ2)**2
     $  +(32./27.)*(U-AM2)**2/(U-AMQ2)**2
     $  +(8./3.)*((T-AM2)**2+AM2*S)/(S*(T-AMQ2))
     $  +(8./3.)*((U-AM2)**2+AM2*S)/(S*(U-AMQ2))
     $  +(8./27.)*AM2*S/((T-AMQ2)*(U-AMQ2))
        SIGL=.5*FAC*SIGL
C          Right squark exchange
        AMQ=AMASS(IDQSS(IQ1+12))
        AMQ2=AMQ**2
        SIGR=(8./3.)*((T-AM2)**2+(U-AM2)**2+2.*AM2*S)/(S**2)
     $  +(32./27.)*(T-AM2)**2/(T-AMQ2)**2
     $  +(32./27.)*(U-AM2)**2/(U-AMQ2)**2
     $  +(8./3.)*((T-AM2)**2+AM2*S)/(S*(T-AMQ2))
     $  +(8./3.)*((U-AM2)**2+AM2*S)/(S*(U-AMQ2))
     $  +(8./27.)*AM2*S/((T-AMQ2)*(U-AMQ2))
        SIGR=.5*FAC*SIGR
        SIG0=.5*(SIGL+SIGR)
C          Total
        SIG=SIG0*QFCN(IQ1,1)*QFCN(IQ2,2)
        CALL SIGFIL(SIG,IQ1,IQ2,1,1)
        SIG=SIG0*QFCN(IQ2,1)*QFCN(IQ1,2)
        CALL SIGFIL(SIG,IQ2,IQ1,1,1)
220   CONTINUE
C
C          Scalar quark + scalar (anti)quark
C
300   CONTINUE
      AMG=AMASS(ISGL)
      AMG2=AMG**2
C          IQ1 and IQ2 loop over left and right (anti)squarks
      DO 310 IQ1=2,25
      DO 320 IQ2=2,25
        IF(.NOT.(GOQ(IQ1,1).AND.GOQ(IQ2,2))) GO TO 320
        JQ1=JS2JT(IQ1)
        JQ2=JS2JT(IQ2)
C        IF(JQ1.GE.12.OR.JQ2.GE.12) GO TO 320
        IFL1=IDQSS(IQ1)
        IFL2=IDQSS(IQ2)
        IFLQ1=IDQ(JQ1)
        IFLQ2=IDQ(JQ2)
C          LLRR is true for left-left or right-right
        IF((IQ1.LE.13.AND.IQ2.LE.13).OR.(IQ1.GT.13.AND.IQ2.GT.13))
     $  THEN
          LLRR=.TRUE.
        ELSE
          LLRR=.FALSE.
        ENDIF
C          Kinematics
        AM1=AMASS(IFL1)
        AM2=AMASS(IFL2)
        AM=AM1
        CALL TWOKIN(0.,0.,AM1,AM2)
        IF(X1.GE.1..OR.X2.GE.1.) GO TO 320
        AMSQ=AM**2
        AM1SQ=AM1**2
        AM2SQ=AM2**2
        S=SHAT
        T=THAT
        U=UHAT
        E1=SQRT(P(1)**2+AM1SQ)
        E2=SQRT(P(2)**2+AM2SQ)
        FAC=PI*ALFQSQ**2/S**2
        FAC=FAC*(S/SCM)*(P(1)*P(2)/(E1*E2))*UNITS
C
C          gl gl ---> qkss qbss
C
        IF(IFL1.EQ.-IFL2) THEN
          SIG=(7./48.+3.*(U-T)**2/(16.*S**2))
     $    *(1.+2.*AMSQ*T/(T-AMSQ)**2+2.*AMSQ*U/(U-AMSQ)**2
     $    +4.*AMSQ**2/((T-AMSQ)*(U-AMSQ)))
          SIG=SIG*FAC*QFCN(1,1)*QFCN(1,2)
          SIG=.5*SIG
C          Another .5 to sum over L and R
          SIG=.5*SIG
          CALL SIGFIL(SIG,1,1,IQ1,IQ2)
        ENDIF
C
C          qk qb ---> qkss qbss
C
        IF(IFLQ1.EQ.-IFLQ2.AND.LLRR) THEN
C          Identical squark-antisquark, LL or RR
          SIG=(2./9.)*(1/(T-AMG2)**2+2/S**2-2/(3*S*(T-AMG2)))
     $     *(-S*T-(T-AMSQ)**2)*FAC*QFCN(JQ1,1)*QFCN(JQ2,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,JQ1,JQ2,IQ1,IQ2)
          SIG=(2./9.)*(1/(U-AMG2)**2+2/S**2-2/(3*S*(U-AMG2)))
     $     *(-S*U-(U-AMSQ)**2)*FAC*QFCN(JQ2,1)*QFCN(JQ1,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,JQ2,JQ1,IQ1,IQ2)
        ELSEIF(IFLQ1.EQ.-IFLQ2.AND..NOT.LLRR) THEN
C          Identical squark-antisquark, LR or RL
          SIG=(2./9.)*AMG2*S/(T-AMG2)**2*FAC*QFCN(JQ1,1)*QFCN(JQ2,2)
        SIG=.5*SIG
          CALL SIGFIL(SIG,JQ1,JQ2,IQ1,IQ2)
          SIG=(2./9.)*AMG2*S/(U-AMG2)**2*FAC*QFCN(JQ2,1)*QFCN(JQ1,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,JQ2,JQ1,IQ1,IQ2)
        ELSEIF(IFLQ1.EQ.IFLQ2.AND.LLRR) THEN
C          Identical squark-squark, LL or RR
          SIG=(1./9.)*AMG2*S*(1/(T-AMG2)**2+1/(U-AMG2)**2
     $    -(2./3.)/((T-AMG2)*(U-AMG2)))*FAC*QFCN(JQ1,1)*QFCN(JQ2,2)
          CALL SIGFIL(SIG,JQ1,JQ2,IQ1,IQ2)
        ELSEIF(IFLQ1.EQ.IFLQ2.AND..NOT.LLRR) THEN
C          Identical squark-squark, LR or RL
          SIG=(2./9.)*(1/(T-AMG2)**2*(-S*T-(T-AM1SQ)*(T-AM2SQ))
     $    +1/(U-AMG2)**2*(-S*U-(U-AM1SQ)*(U-AM2SQ)))
     $    *FAC*QFCN(JQ1,1)*QFCN(JQ2,2)
          CALL SIGFIL(SIG,JQ1,JQ2,IQ1,IQ2)
        ELSEIF(IFL1*IFL2.LT.0.AND.LLRR) THEN
C          Nonidentical squark-antisquark, LL or RR
          SIG=(2./9.)*(-S*T-(T-AM1SQ)*(T-AM2SQ))/(T-AMG2)**2*FAC
     $    *QFCN(JQ1,1)*QFCN(JQ2,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,JQ1,JQ2,IQ1,IQ2)
          SIG=(2./9.)*(-S*U-(U-AM1SQ)*(U-AM2SQ))/(U-AMG2)**2*FAC
     $    *QFCN(JQ2,1)*QFCN(JQ1,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,JQ2,JQ1,IQ1,IQ2)
        ELSEIF(IFL1*IFL2.LT.0.AND..NOT.LLRR) THEN
C          Nonidentical squark-antisquark, LR or RL
          SIG=(2./9.)*AMG2*S/(T-AMG2)**2*FAC*QFCN(JQ1,1)*QFCN(JQ2,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,JQ1,JQ2,IQ1,IQ2)
          SIG=(2./9.)*AMG2*S/(U-AMG2)**2*FAC*QFCN(JQ2,1)*QFCN(JQ1,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,JQ2,JQ1,IQ1,IQ2)
        ELSEIF(IFL1*IFL2.GT.0.AND.LLRR) THEN
C          Nonidentical squark-squark, LL or RR
          SIG=(2./9.)*AMG2*S/(T-AMG2)**2*FAC*QFCN(JQ1,1)*QFCN(JQ2,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,JQ1,JQ2,IQ1,IQ2)
          SIG=(2./9.)*AMG2*S/(U-AMG2)**2*FAC*QFCN(JQ2,1)*QFCN(JQ1,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,JQ2,JQ1,IQ1,IQ2)
        ELSEIF(IFL1*IFL2.GT.0.AND..NOT.LLRR) THEN
C          Nonidentical squark-squark, LR or RL
          SIG=(2./9.)*(-S*T-(T-AM1SQ)*(T-AM2SQ))/(T-AMG2)**2*FAC
     $    *QFCN(JQ1,1)*QFCN(JQ2,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,JQ1,JQ2,IQ1,IQ2)
          SIG=(2./9.)*(-S*U-(U-AM1SQ)*(U-AM2SQ))/(U-AMG2)**2*FAC
     $    *QFCN(JQ2,1)*QFCN(JQ1,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,JQ2,JQ1,IQ1,IQ2)
        ELSE
          STOP99
        ENDIF
C
C          q1 + q1bar --> q2ss + q2ssbar
C
        IF(IFLQ1.EQ.-IFLQ2.AND.LLRR) THEN
          DO 330 JQIN1=2,10,2
            IF(JQIN1.EQ.JQ1.OR.JQIN1.EQ.JQ2) GO TO 330
            JQIN2=MATCH(JQIN1,4)
            SIG=(4./9.)*(-S*T-(T-AM1SQ)**2)/S**2*FAC
     $      *QFCN(JQIN1,1)*QFCN(JQIN2,2)
            SIG=.5*SIG
            CALL SIGFIL(SIG,JQIN1,JQIN2,IQ1,IQ2)
            SIG=(4./9.)*(-S*U-(U-AM1SQ)**2)/S**2*FAC
     $      *QFCN(JQIN2,1)*QFCN(JQIN1,2)
            SIG=.5*SIG
            CALL SIGFIL(SIG,JQIN2,JQIN1,IQ1,IQ2)
330       CONTINUE
        ENDIF
320   CONTINUE
310   CONTINUE
C
C          Scalar quark + gluino
C
      AMG=AMASS(ISGL)
      AMG2=AMG**2
      DO 400 IQ=2,25
        AMQ=AMASS(IDQSS(IQ))
        AMQ2=AMQ**2
        JQ=JS2JT(IQ)        
C
C          Jet 1 = scalar quark
        IF(.NOT.(GOQ(JQ,1).AND.GOQ(1,2))) GO TO 410
        CALL TWOKIN(0.,0.,AMQ,AMG)
        IF(X1.GE.1..OR.X2.GE.1.) GO TO 410
        S=SHAT
        E1=SQRT(P(1)**2+AMQ2)
        E2=SQRT(P(2)**2+AMG2)
        FAC=PI*ALFQSQ**2/S**2
        FAC=FAC*S/SCM*P(1)*P(2)/(E1*E2)*UNITS
C
        T=THAT-AMQ2
        U=UHAT-AMG2
        SIG=FQG(S,T,U)*FAC/12.*QFCN(JQ,1)*QFCN(1,2)
        SIG=.5*SIG
        SIG=.5*SIG
        CALL SIGFIL(SIG,JQ,1,IQ,1)
C
        T=UHAT-AMQ2
        U=THAT-AMG2
        SIG=FQG(S,T,U)*FAC/12.*QFCN(1,1)*QFCN(JQ,2)
        SIG=.5*SIG
        SIG=.5*SIG
        CALL SIGFIL(SIG,1,JQ,IQ,1)
C
C          Jet 2 = scalar quark
410     IF(.NOT.(GOQ(1,1).AND.GOQ(JQ,2))) GO TO 400
        CALL TWOKIN(0.,0.,AMG,AMQ)
        IF(X1.GE.1..OR.X2.GE.1.) GO TO 400
        S=SHAT
        E1=SQRT(P(1)**2+AMG2)
        E2=SQRT(P(2)**2+AMQ2)
        FAC=PI*ALFQSQ**2/S**2
        FAC=FAC*S/SCM*P(1)*P(2)/(E1*E2)*UNITS
C 
        T=UHAT-AMQ2
        U=THAT-AMG2
        SIG=FQG(S,T,U)*FAC/12.*QFCN(1,1)*QFCN(JQ,2)
        SIG=.5*SIG
        SIG=.5*SIG
        CALL SIGFIL(SIG,1,JQ,1,IQ)
C
        T=THAT-AMQ2
        U=UHAT-AMG2
        SIG=FQG(S,T,U)*FAC/12.*QFCN(JQ,1)*QFCN(1,2)
        SIG=.5*SIG
        SIG=.5*SIG
        CALL SIGFIL(SIG,JQ,1,1,IQ)
400   CONTINUE
C
C          Calculate gaugino AND slepton cross sections only for MSSM
C
      IF(GOMSSM) CALL SIGSSZ
      IF(GOMSSM) CALL SIGSSL
C
      RETURN
      END
