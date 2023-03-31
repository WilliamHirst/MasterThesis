CDECK  ID>, EVOL03.
      SUBROUTINE EVOL03
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-        Setup for process 3 (DRELLYAN)
C-        Lorentz frames and perform initial and final QCD jet
C-        evolution in leading-log approximation.
C-
C-   Created  13-AUG-1991   Frank E. Paige,Serban D. Protopopescu
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
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
      INTEGER MXJETS
      PARAMETER (MXJETS=10)
      COMMON/PJETS/PJETS(5,MXJETS),IDJETS(MXJETS),QWJET(5),IDENTW 
     $,PPAIR(5,4),IDPAIR(4),JPAIR(4),NPAIR,IFRAME(MXJETS)
      SAVE /PJETS/
      INTEGER   IDJETS,IDENTW,IDPAIR,JPAIR,NPAIR,IFRAME
      REAL      PJETS,QWJET,PPAIR
      INTEGER   MXJSET,JPACK
      PARAMETER (MXJSET=400,JPACK=1000)
      COMMON/JETSET/NJSET,PJSET(5,MXJSET),JORIG(MXJSET),JTYPE(MXJSET),
     $JDCAY(MXJSET)
      SAVE /JETSET/
      INTEGER   NJSET,JORIG,JTYPE,JDCAY
      REAL      PJSET
      COMMON/JWORK/ZZC(MXJSET),JMATCH(MXJSET),TNEW,P1CM(4),
     1J1,J2,J3,J4,J5,E1CM,E2CM,E3CM,E4CM,E5CM
      SAVE /JWORK/
      LOGICAL TNEW
      EQUIVALENCE (J1,JJ(1)),(E1CM,EE(1))
      INTEGER   JMATCH,J1,J2,J3,J4,J5,JJ(5)
      REAL      ZZC,P1CM,E1CM,E2CM,E3CM,E4CM,E5CM,EE(5)
      COMMON/JWORK2/JVIR(2),PFINAL(5),SGN,ZMIN,ZMAX,DZMAX,JET,GLFORC(2),
     $ZGOOD,JIN(400),FXTEST(MXJSET)
      SAVE /JWORK2/
      LOGICAL GLFORC,ZGOOD
      INTEGER   JVIR,JET,JIN
      REAL      PFINAL,SGN,ZMIN,ZMAX,DZMAX,FXTEST
      INTEGER MXGOQ,MXGOJ
      PARAMETER (MXGOQ=85,MXGOJ=8)
      COMMON/Q1Q2/GOQ(MXGOQ,MXGOJ),GOALL(MXGOJ),GODY(4),STDDY,
     $GOWW(25,2),ALLWW(2),GOWMOD(25,MXGOJ)
      SAVE /Q1Q2/
      LOGICAL GOQ,GOALL,GODY,STDDY,GOWW,ALLWW,GOWMOD
      COMMON/FRAME/FRAME(5,3),N0JETS,N0W,N0PAIR
      SAVE /FRAME/
      INTEGER   N0JETS,N0W,N0PAIR
      REAL      FRAME
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
      REAL    EVOLMS,BP
      INTEGER I,K,J,NJFINL
C----------------------------------------------------------------------
C
C          Add W momentum and recoil jets
      N0JETS=NJSET+1
      IF(.NOT.STDDY) THEN
        DO 101 I=3,NJET
          NJSET=NJSET+1
          JORIG(NJSET)=JPACK*I
          JTYPE(NJSET)=IDJETS(I)
          JDCAY(NJSET)=0
          DO 105 K=1,5
105       PJSET(K,NJSET)=PJETS(K,I)
          IFRAME(I)=1
101     CONTINUE
        NJSET=NJSET+1
        N0W=NJSET
        JORIG(NJSET)=0
        JTYPE(NJSET)=IDENTW
        JDCAY(NJSET)=(N0W+1)*JPACK+N0W+2
        DO 120 K=1,5
120     PJSET(K,NJSET)=QWJET(K)
      ENDIF
C
C          Add W decays
      DO 110 I=1,2
        NJSET=NJSET+1
        JORIG(NJSET)=JPACK*I
        JTYPE(NJSET)=IDJETS(I)
        JDCAY(NJSET)=0
        DO 115 K=1,5
115     PJSET(K,NJSET)=PJETS(K,I)
        IFRAME(I)=2
        IF(STDDY) IFRAME(I)=1
110   CONTINUE
C
C          Set flags and maximum off-shell masses and generate
C          initial QCD parton shower.
C
      CALL ISTRAD(WFUDGE)
C
      IF(NJSET.LT.0) RETURN
C
C          Final state evolution.
C          Define Lorentz frames and JMATCH pointers for jet evolution
C          and fragmentation.
C
      IF(STDDY) THEN
        CALL IFRAMS(3,4,1,.FALSE.)
      ELSE
        CALL IFRAMS(N0W+1,N0W+2,2,.FALSE.)
        CALL IFRAMS(N0JETS,N0W,1,.FALSE.)          
      ENDIF
C
C          Set maximum off-shell masses and JDCAY flags.
C
      IF(STDDY) THEN
        NJFINL=3
        DO 310 J=3,4
          IF(IABS(JTYPE(J)).LT.10) THEN
            PJSET(5,J)=QMW
            JDCAY(J)=-1
          ENDIF
310     CONTINUE
      ELSE
        NJFINL=N0JETS
        DO 320 J=N0W+1,N0W+2
          IF(IABS(JTYPE(J)).LT.10) THEN
            PJSET(5,J)=QMW
            JDCAY(J)=-1
          ENDIF
320     CONTINUE
C          Need fudge factor for DRELLYAN
        DO 321 J=N0JETS,N0W
          IF(IABS(JTYPE(J)).LT.10) THEN
            PJSET(5,J)=EVOLMS(J,WFUDGE)
            JDCAY(J)=-1
          ENDIF
321     CONTINUE
      ENDIF
C
C          Produce final-state QCD parton cascade
C
      CALL QCDJET(NJFINL)
C
C          Reset FRAME using W momentum modified by evolution
      IF(.NOT.STDDY) THEN
        BP=0.
        DO 400 K=1,3
400     BP=BP+FRAME(K,1)*PJSET(K,N0W)
        BP=BP/FRAME(5,1)
        DO 410 K=1,3
          FRAME(K,2)=PJSET(K,N0W)+FRAME(K,1)*PJSET(4,N0W)/FRAME(5,1)
     $    +FRAME(K,1)*BP/(FRAME(4,1)+FRAME(5,1))
410     CONTINUE
        FRAME(4,2)=FRAME(4,1)*PJSET(4,N0W)/FRAME(5,1)+BP
      ENDIF
C
      RETURN
      END
