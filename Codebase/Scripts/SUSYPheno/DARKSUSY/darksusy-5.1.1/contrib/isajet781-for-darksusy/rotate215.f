!
      SUBROUTINE ROTATE215(G215)
!
!Purpose: Rotate all the matrices from the quark mass
!         basis to the current basis using the 
!         previously derived V's
!
      IMPLICIT NONE
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
      DOUBLE COMPLEX FUP(3,3),FDP(3,3),LUP(3,3),LDP(3,3)
!
      DOUBLE COMPLEX FUDUM(3,3),FDDUM(3,3),LUDUM(3,3),LDDUM(3,3)
!
      DOUBLE COMPLEX VLUT(3,3),VRUT(3,3),VLDT(3,3),VRDT(3,3)
      DOUBLE COMPLEX CMATMUL,G215(215)
      INTEGER I,J
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
      DO I=1,3
        DO J=1,3
          FUP(I,J)=G215(3+(I-1)*3+J)
          FDP(I,J)=G215(12+(I-1)*3+J)
          LUP(I,J)=G215(33+(I-1)*3+J)
          LDP(I,J)=G215(42+(I-1)*3+J)
        END DO
      END DO
!
!Now rotate the matrices back
!
      DO I=1,3 
        DO J=1,3
          FUDUM(I,J)=CMATMUL(0,FUP,VRUT,I,J)
          FDDUM(I,J)=CMATMUL(0,FDP,VRDT,I,J)
          LUDUM(I,J)=CMATMUL(0,LUP,VRUT,I,J)
          LDDUM(I,J)=CMATMUL(0,LDP,VRDT,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          FUP(I,J)=CMATMUL(0,VLUT,FUDUM,I,J)
          LUP(I,J)=CMATMUL(0,VLUT,LUDUM,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          G215(3+(I-1)*3+J)=CMATMUL(1,VLUT,FUDUM,I,J)
          G215(12+(I-1)*3+J)=CMATMUL(1,VLDT,FDDUM,I,J)
          G215(33+(I-1)*3+J)=CMATMUL(1,VLUT,LUDUM,I,J)
          G215(42+(I-1)*3+J)=CMATMUL(1,VLDT,LDDUM,I,J)
        END DO
      END DO
!
      RETURN
      END
