CDECK  ID>, IPRTNS.
      SUBROUTINE IPRTNS(NPRTNS,PRTNS,IDQ)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-      Fill PINITS common block
C-   Inputs  : 
C-     IDQ(2)= id's of partons starting reaction
C-
C-   Created  10-OCT-1991   Serban D. Protopopescu
C-   Renamed from IPINIT to avoid name clash with Cern Library
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NPRTNS,IDQ(2)
      REAL    PRTNS(4,NPRTNS)
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
      COMMON/PINITS/PINITS(5,2),IDINIT(2)
      SAVE /PINITS/
      INTEGER   IDINIT
      REAL      PINITS
      REAL AMASS, AM1SQ,AM2SQ,ROOT,QPL,QMN,P1PL,P1MN,P2PL,P2MN
      INTEGER I
C----------------------------------------------------------------------
C          sum P+ and P-, shat
C          assumes sum of transverse momenta is zero
      QPL=0
      QMN=0
      DO 1 I=1,NPRTNS
        QPL=QPL+PRTNS(4,I)+PRTNS(3,I)
        QMN=QMN+PRTNS(4,I)-PRTNS(3,I)
    1 CONTINUE
      SHAT=QPL*QMN
C
C          fill PINITS
      DO 2 I=1,2
        IDINIT(I)=IDQ(I)
        PINITS(5,I)=AMASS(IDQ(I))
        PINITS(1,I)=0.
        PINITS(2,I)=0.
    2 CONTINUE
C          and solve initial kinematics
      AM1SQ=PINITS(5,1)**2
      AM2SQ=PINITS(5,2)**2
      ROOT=SQRT((QPL*QMN-AM1SQ-AM2SQ)**2-4.*AM1SQ*AM2SQ)
      P1PL=(QPL*QMN+AM1SQ-AM2SQ+ROOT)/(2.*QMN)
      P1MN=AM1SQ/P1PL
      P2MN=(QPL*QMN+AM2SQ-AM1SQ+ROOT)/(2.*QPL)
      P2PL=AM2SQ/P2MN
      PINITS(3,1)=.5*(P1PL-P1MN)
      PINITS(4,1)=.5*(P1PL+P1MN)
      PINITS(3,2)=.5*(P2PL-P2MN)
      PINITS(4,2)=.5*(P2PL+P2MN)
  999 RETURN
      END
