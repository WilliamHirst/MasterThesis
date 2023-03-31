c------------------------------------------------------------c


c------------------------------------------------------------c
C============================================================================
C Interface to CTEQ5L parton distribution
C Uses CTEQ5L from main ISAJET
C============================================================================
      DOUBLE PRECISION FUNCTION CTQ5PDF (IPARTON, X, Q, IPRT)
      IMPLICIT NONE
      INTEGER IPARTON,IPRT,IFL
      DOUBLE PRECISION X,Q,CTEQ5L,SUM5L,RAT5L
      IFL=IPARTON
      IF(IFL.GE.3) IFL=-IFL
      IF(IFL.EQ.-1) THEN
        SUM5L=CTEQ5L(-1,X,Q)
        RAT5L=CTEQ5L(-2,X,Q)
        CTQ5PDF=SUM5L/(1.D0+RAT5L)
      ELSEIF(IFL.EQ.-2) THEN
        SUM5L=CTEQ5L(-1,X,Q)
        RAT5L=CTEQ5L(-2,X,Q)
        CTQ5PDF=SUM5L*RAT5L/(1.D0+RAT5L)
      ELSEIF(IFL.GE.-5.AND.IFL.LE.5) THEN
        CTQ5PDF=CTEQ5L(IFL,X,Q)
      ELSE
        CTQ5PDF=0
      ENDIF
      RETURN
      END
