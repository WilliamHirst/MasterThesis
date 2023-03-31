CDECK  ID>, EVOLMS.
      FUNCTION EVOLMS(J,FUDGE)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-     Set evolution mass scale for parton J
C-
C-   Returned value  : maximum mass
C-
C-   Inputs  : 
C-     J    = index to PJSET array
C-     FUDGE= fudge factor
C-
C-   Created  16-AUG-1991   Serban D. Protopopescu
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
      REAL    EVOLMS,FUDGE
      INTEGER J
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
C----------------------------------------------------------------------
C
      IF ( USELIM ) THEN
        EVOLMS=SQRT(PJSET(1,J)**2+PJSET(2,J)**2)*CONCUT
      ELSE
        EVOLMS=FUDGE*SQRT(QSQ)
      ENDIF
  999 RETURN
      END
