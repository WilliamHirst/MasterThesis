!
      SUBROUTINE ROTATESM(GROT)
!
!Purpose: Rotate the Yukawa matrices from the basis in which the
!         chosen Yukawa matrix is diagonal to the general current
!         basis.
!
      IMPLICIT NONE
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
      DOUBLE COMPLEX FU(3,3),FD(3,3),FUP(3,3),FDP(3,3)
      DOUBLE COMPLEX FUDUM(3,3),FDDUM(3,3)
      DOUBLE COMPLEX GROT(32)
!
      DOUBLE COMPLEX VLUT(3,3),VRUT(3,3),VLDT(3,3),VRDT(3,3)
      DOUBLE COMPLEX CMATMUL
      INTEGER I,J
!
!Convert the G's to matrices
!
      DO I=1,3
        DO J=1,3
          FUP(I,J)=GROT(3+(I-1)*3+J)
          FDP(I,J)=GROT(12+(I-1)*3+J)
        END DO
      END DO
!
!Calculate the transposes
!
      DO I=1,3
        DO J=1,3
          VLUT(I,J)=VLU(J,I)
          VRUT(I,J)=VRU(J,I)
          VLDT(I,J)=VLD(J,I)
          VRDT(I,J)=VRD(J,I)
        END DO
      END DO
!
!Now rotate the matrices
!
      DO I=1,3
        DO J=1,3
          FUDUM(I,J)=CMATMUL(0,FUP,VRUT,I,J)
          FDDUM(I,J)=CMATMUL(0,FDP,VRDT,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          GROT(3+(I-1)*3+J)=CMATMUL(1,VLUT,FUDUM,I,J)
          GROT(12+(I-1)*3+J)=CMATMUL(1,VLDT,FDDUM,I,J)
        END DO
      END DO
!
      RETURN
      END
