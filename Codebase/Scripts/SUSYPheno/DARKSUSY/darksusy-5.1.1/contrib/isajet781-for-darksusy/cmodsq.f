!
      FUNCTION CMODSQ(A)
!
!Purpose: Calculates the mod squared of a value. This will only
!         be a necessary function if A is complex
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX A,CMODSQ
!
      CMODSQ=CONJG(A)*A
!
      RETURN
      END
