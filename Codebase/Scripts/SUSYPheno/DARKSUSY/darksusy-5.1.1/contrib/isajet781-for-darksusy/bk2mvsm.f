CDECK  ID>, RGETOOLS.
      SUBROUTINE BK2MVSM(GIN,GOUT)
!
!Purpose: To convert the book notation G's to Martin and Vaughn notation.
!         The same subroutine will convert back as well.
!         For SMRGE it is sufficient to convert the Yukawa matrices.
!
!With f_{u,d,e} in Baer/Tata notation and Y_{u,d,e} in Martin/Vaughn
!notation, the conversion is:
!
!     Y_{u,d,e}= {f_{u,d,e}}^T
!
      IMPLICIT NONE
!
      DOUBLE PRECISION GIN(32),GOUT(32),YU(3,3),YD(3,3),YE(3,3)
      INTEGER I,J
!
!First set GOUT to GIN
!
     	DO I=1,32
     	  GOUT(I)=GIN(I)
     	END DO
!
!Convert input into 3x3 matrices
!
      DO I=1,3
        DO J=1,3
          YU(I,J)=GIN(3+(I-1)*3+J)
          YD(I,J)=GIN(12+(I-1)*3+J)
          YE(I,J)=GIN(21+(I-1)*3+J)
        END DO
      END DO
!
!Reset relevant parts of GOUT
!
      DO I=1,3
        DO J=1,3
          GOUT(3+(I-1)*3+J)=YU(J,I)
          GOUT(12+(I-1)*3+J)=YD(J,I)
          GOUT(21+(I-1)*3+J)=YE(J,I)
        END DO
      END DO
!
      RETURN
      END
