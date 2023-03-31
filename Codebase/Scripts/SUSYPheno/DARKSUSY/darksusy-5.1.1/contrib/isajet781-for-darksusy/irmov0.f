CDECK  ID>, IRMOV0.
      SUBROUTINE IRMOV0
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-      remove 0's from PJSET
C-
C-   Created  15-OCT-1991   Serban D. Protopopescu
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
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
      INTEGER NCOUNT,I,J,K
C----------------------------------------------------------------------
C
C         remove zeroes
      NCOUNT=NJSET
      DO 160 I=3,NJSET  
  151   IF (PJSET(4,I).EQ.0.AND.I.LT.NCOUNT) THEN
          DO 155 K=I+1,NCOUNT
            DO 154 J=1,5    
              PJSET(J,K-1)=PJSET(J,K)   
  154       CONTINUE    
            JORIG(K-1)=JORIG(K) 
            JTYPE(K-1)=JTYPE(K) 
            JDCAY(K-1)=JDCAY(K) 
            ZZC(K-1)=ZZC(K)
            JMATCH(K-1)=JMATCH(K)
            IF(JMATCH(K-1).GT.I) JMATCH(K-1)=JMATCH(K-1)-1
  155     CONTINUE  
          NCOUNT=NCOUNT-1
          GOTO 151
        ENDIF
  160 CONTINUE  
      NJSET=NCOUNT  
C          remove last one if 0
      IF(PJSET(4,NJSET).EQ.0) NJSET=NJSET-1  
  999 RETURN
      END
