CDECK  ID>, IFRAMS.
      SUBROUTINE IFRAMS(N1,N2,IFR,PAIR)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-     Initialize a center of mass frame for partons N1 to N2
C-     partons must be consecutive unless PAIR is true
C-     
C-   Inputs  : 
C-   N1  = first parton 
C-   N2  = last parton
C-   IFR = index of frame
C-   PAIR= if false N1, N2 denote a range
C-      if true N1 and N2 form a pair
C-
C-   Created  14-AUG-1991   Serban D. Protopopescu
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
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
      COMMON/FRAME/FRAME(5,3),N0JETS,N0W,N0PAIR
      SAVE /FRAME/
      INTEGER   N0JETS,N0W,N0PAIR
      REAL      FRAME
      INTEGER I,J,K,JADD,N1,N2,IFR
      DOUBLE PRECISION DPASS(5),DSUM(5)
      LOGICAL PAIR
C----------------------------------------------------------------------
C
      IF ( N2-N1.EQ.1.OR.PAIR ) THEN
        JMATCH(N1)=N2
        JMATCH(N2)=N1
        JADD=N2-N1
      ELSE
        JADD=1
        DO 201 I=N1,N2
          JMATCH(I)=JPACK*N1+N2
201     CONTINUE
      ENDIF
C          Need double precision boosts
      CALL DBLVEC(PJSET(1,N1),DSUM)
      DO 211 I=N1+JADD,N2
        CALL DBLVEC(PJSET(1,I),DPASS)
        DO 210 K=1,4
210     DSUM(K)=DSUM(K)+DPASS(K)
        DSUM(5)=DSQRT(DSUM(4)**2-DSUM(1)**2-DSUM(2)**2-DSUM(3)**2)
211   CONTINUE
      DO 212 K=1,5
        FRAME(K,IFR)=DSUM(K)
212   CONTINUE
C
C          Set up and generate final state QCD parton shower.
C          Boost PJSET with -FRAME.
C
      DO 240 J=N1,N2,JADD
        CALL DBOOST(-1,FRAME(1,IFR),PJSET(1,J))
240   CONTINUE
C
999   RETURN
      END
