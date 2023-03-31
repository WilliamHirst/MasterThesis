CDECK  ID>, ISTRAD.
      SUBROUTINE ISTRAD(FUDGE)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-      Set parameters and call QCDINI to generate initial
C-      state radiation
C-   Inputs  : 
C-     FUDGE= fudge factor
C-
C-   Created  16-AUG-1991   Serban D. Protopopescu
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
      REAL    FUDGE
      COMMON /LIMEVL/ ETTHRS,CONCUT,USELIM
      SAVE /LIMEVL/
      REAL ETTHRS,CONCUT
      LOGICAL USELIM
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
      REAL    OFF
      INTEGER I
C----------------------------------------------------------------------
C
      IF ( USELIM.AND.CONCUT.LT.1.0 ) THEN
        OFF=ETTHRS
      ELSEIF( .NOT.USELIM) THEN
        OFF=SQRT(QSQ)*FUDGE
      ELSE
        OFF=SQRT(QSQ)
      ENDIF
      DO 150 I=1,2
        PJSET(5,I)=-OFF
150   JDCAY(I)=-2
      JMATCH(1)=0
      JMATCH(2)=0
C
      CALL QCDINI(1,2)
  999 RETURN
      END
