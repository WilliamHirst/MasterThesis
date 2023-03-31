CDECK  ID>, QCDT.
      SUBROUTINE QCDT(J)
C
C          Auxiliary routine for QCDJET. Calculate ZC and store in 
C          ZZC(J).  Generate new mass with ZC and store in PJSET(5,J).
C
C          Must include 1/2 symmetry factor in GAMGG. No fix is needed
C          in QCDZ since GAMGG+2*GAMQQ is used as the normalization.
C
C          Include GM, W+, W-, and Z0 radiation.
C
C          Ver 7.20: Anomalous dimensions were coded incorrectly!
C
      IMPLICIT NONE
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
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
      COMMON/QCDPAR/ALAM,ALAM2,CUTJET,ISTRUC
      SAVE /QCDPAR/
      INTEGER   ISTRUC
      REAL      ALAM,ALAM2,CUTJET
      COMMON/CONST/PI,SQRT2,ALFA,GF,UNITS
      SAVE /CONST/
      REAL      PI,SQRT2,ALFA,GF,UNITS
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
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
C
      REAL AM0,AM1,AM2,AMASS,T0,T1,T2,ZC,B0,GAMEW,GAMQQ,GAMGG,GAM,GAMZC
      REAL AM1W,AM2W,T1W,T2W,TERM,GB,PROB,RANF,RND,POW,AMNEW,AMOLD
      REAL POWEW
      INTEGER J,JTLV1,NF,IQ,JTABS,IW,JT0,JT1,IFL1,I
      INTEGER JWTYPE(4)
      DATA JWTYPE/10,80,-80,90/
C
C          Set ZC = 0 and return for W+- or Z0
C
      JTABS=IABS(JTYPE(J))
      IF(JTABS.GE.80.AND.JTABS.LE.90) THEN
        ZZC(J)=0.
        RETURN
      ENDIF
C
C          Calculate ZC
C
      AM0=PJSET(5,J)
      JTLV1=JTYPE(J)
      AM1=AMASS(JTLV1)+CUTJET
      AM2=CUTJET
      IF(AM1+AM2.GE.AM0) GO TO 300
      T0=AM0**2
      T1=AM1**2
      T2=AM2**2
C          Fix floating point problem
C     ZC=(T0-T1+T2-SQRT((T0-T1-T2)**2-4*T1*T2))/(2*T0)
      ZC=2*T2/(T0-T1+T2+SQRT((T0-T1-T2)**2-4*T1*T2))
      ZZC(J)=ZC
C          Count light fermions
      NF=3
      DO 110 IQ=4,6
      IF(AM0.LT.2*AMASS(IQ)) GO TO 120
      NF=NF+1
110   CONTINUE
120   B0=11.-2.*NF/3.
C
C          Calculate GAMMA(ZC) and GAMEW for quarks
C
      GAMEW=0.
C
C          Initial gluon
      IF(JTABS.EQ.9) THEN
        GAMQQ=(1.-2.*ZC)*(1.-ZC*(1.-ZC))/3.
        GAMGG=12.*ALOG((1.-ZC)/ZC)-9.*(1.-2.*ZC)-6.*GAMQQ
        GAMGG=0.5*GAMGG
        GAM=GAMGG+NF*GAMQQ
C
C          Initial quark
      ELSEIF(JTABS.LT.9) THEN
        GAMZC=2.*ALOG((1-ZC)/ZC)-1.5*(1.-2.*ZC)
        GAM=4./3.*GAMZC
        GAMEW=ALFA/(2.*PI)*AQ(JTABS,1)**2*GAMZC
        IF(AM0.GT.WMASS(4)) THEN
          DO 130 IW=2,4
            JT0=2*IABS(JTYPE(J))
            IF(JTYPE(J).LT.0) JT0=JT0+1
            JT1=MATCH(JT0,IW)
            IF(JT1.EQ.0) GO TO 130
            JT1=MATCH(JT1,4)
            IFL1=JT1/2
            AM1W=AMASS(IFL1)
            AM2W=AMASS(JWTYPE(IW))
            IF(AM1W+AM2W.GE.AM0) GO TO 130
            T1W=AM1W**2
            T2W=AM2W**2
C          Fix floating underflow
C           ZC=(T0-T1W+T2W-SQRT((T0-T1W-T2W)**2-4*T1W*T2W))/(2*T0)
            ZC=2*T2W/(T0-T1W+T2W+SQRT((T0-T1W-T2W)**2-4*T1W*T2W))
            GAMZC=2.*ALOG((1-ZC)/ZC)-1.5*(1.-2.*ZC)
            TERM=(AQ(JTABS,IW)**2+BQ(JTABS,IW)**2)*GAMZC
            GAMEW=GAMEW+ALFA/(2.*PI)*TERM
130       CONTINUE
        ENDIF
C
C          Initial diquark
      ELSEIF(MOD(JTABS,100).EQ.0) THEN
        GAM=8./3.*ALOG((1-ZC)/ZC)-2.*(1.-2.*ZC)
C
C          Initial gluino
      ELSEIF(JTABS.EQ.29) THEN
        GAM=6.*ALOG((1.-ZC)/ZC)-9./2.*(1.-2.*ZC)
C
C          Initial squark
      ELSEIF(JTABS.GT.20.AND.JTABS.LT.29) THEN
        GAM = 8./3.*(ALOG((1.-ZC)/ZC)-(1.-2.*ZC))
      ENDIF
C
C          Generate new mass
C
      GB=2*GAM/B0
      PROB=(ALOG(AM1/ALAM)/ALOG(AM0/ALAM))**GB
      PROB=PROB*(AM1/AM0)**(2.*GAMEW)
      IF(PROB.GT.RANF()) GO TO 300
      RND=RANF()
      POW=(1.-(1.-PROB)*RND)**(1./GB)
      AMNEW=ALAM*(AM0/ALAM)**POW
C          For quark, add effect of GM, W+-, Z0 radiation
      IF(IABS(JTYPE(J)).LT.9) THEN
        DO 200 I=1,NTRIES
          AMOLD=AMNEW
          POWEW=POW/((AMOLD/AM0)**(2.*GAMEW))**(1./GB)
          AMNEW=ALAM*(AM0/ALAM)**POWEW
          IF(ABS(AMNEW-AMOLD).LT.0.001*AMOLD) GO TO 210
200    CONTINUE
      ENDIF
210   IF(AMNEW.LE.AM1) GO TO 300
      PJSET(5,J)=AMNEW
      RETURN
C
C          Final parton -- set mass to physical value
C
300   PJSET(5,J)=AM1-CUTJET
      JDCAY(J)=0
      RETURN
      END
