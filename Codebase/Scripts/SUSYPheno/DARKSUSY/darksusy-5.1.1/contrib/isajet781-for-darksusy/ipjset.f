CDECK  ID>, IPJSET.
      SUBROUTINE IPJSET
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-       Initialize PJSET starting from PJETS
C-
C-   Created  14-AUG-1991   Serban D. Protopopescu
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
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
      INTEGER I,K
C----------------------------------------------------------------------
      DO 110 I=1,NJET
        NJSET=NJSET+1
        JORIG(NJSET)=JPACK*I
        JTYPE(NJSET)=IDJETS(I)
        JDCAY(NJSET)=0
        DO 115 K=1,5
115     PJSET(K,NJSET)=PJETS(K,I)
        IFRAME(I)=1
110   CONTINUE
  999 RETURN
      END
