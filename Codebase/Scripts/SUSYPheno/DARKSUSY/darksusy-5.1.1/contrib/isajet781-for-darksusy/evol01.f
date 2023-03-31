CDECK  ID>, EVOL01.
      SUBROUTINE EVOL01
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-        Setup for process 1 (TWOJET)
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
      COMMON/FRAME/FRAME(5,3),N0JETS,N0W,N0PAIR
      SAVE /FRAME/
      INTEGER   N0JETS,N0W,N0PAIR
      REAL      FRAME
      REAL    EVOLMS
      INTEGER I,K,J,NJSAVE,NJFINL
C----------------------------------------------------------------------
C
C          Copy momenta from /PJETS/ to /JETSET/
      N0JETS=NJSET+1
      CALL IPJSET
      NJSAVE=NJSET
C
C          Set flags and maximum off-shell masses and generate
C          initial QCD parton shower.
C
      CALL ISTRAD(1.0)
C
      IF(NJSET.LT.0) RETURN
C
C          Final state evolution.
C          Define Lorentz frames and JMATCH pointers for jet evolution
C          and fragmentation.
C
      CALL IFRAMS(N0JETS,NJSAVE,1,.FALSE.)
C
C          Set maximum off-shell masses and JDCAY flags.
C
      NJFINL=N0JETS
      DO 310 J=N0JETS,NJSAVE
        IF(IABS(JTYPE(J)).LT.10) THEN
          PJSET(5,J)=EVOLMS(J,1.0)
          JDCAY(J)=-1
        ENDIF
310   CONTINUE
C
C          Produce final-state QCD parton cascade
C
      CALL QCDJET(NJFINL)
C
      RETURN
      END
