C----------------------------------------------------------------------
      SUBROUTINE YUKDIAG(G,INIT,YDIAG,VL,VR)
C----------------------------------------------------------------------
C
C  Assembles Yukawa matrix in weak eigenbasis from RGE running vector G(i)
C  starting from i=INIT and diagonalizes it.
C
C   NOTE: In our notation fermion masses are proportional to transposed
C         Yukawas in gauge eigenbasis, i.e.  m ~ VR Y^T VL^dagger,
C         while SURG111 follows m ~ Y convention.
C         One has to transpose Yukasas when passing
C         into or from those subroutines.
C
C   NOTE: VL and VR contain eigenvectors as rows, while Z - as columns. 
c
c  Ref: Misak, Pokorski & Rosiek  hep-ph/9703442
c
c  Created: 01/19/06 by Azar Mustafayev.
C 
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C   Input:
C     G(4) = Y_u(1,1)
C     G(5) = Y_u(1,2)
C     G(6) = Y_u(1,3)
C     ...    ...
C     G(12) = Y_u(3,3)
C     G(13)-G(21) = Y_d
C     G(22)-G(30) = Y_e
c   Output:
c    VL,VR - Left and right rotation matrices: YDIAG = VL^* YUK VR^T
C    YDIAG - diagonal Yukawa matrix
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv            
      IMPLICIT NONE
      REAL*8 G(157),YDIAG(3,3),VL(3,3),VR(3,3)
      INTEGER INIT
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
      REAL*8 YUK(3,3),YUKH(3,3),A(3,3),W(3),Z(3,3),WORK(3),I3(3,3)
c  YUK(3,3)  - Yukawa matrix in weak eigenbasis
c  YUKH(3,3) - hermitian conjugated YUK
c  A(3,3)  - axiliary matrix to be diagonalized
c  W(3)    - vector of eigenvalues for AR
c  Z(3,3)  - matrix of eigenvalues of AR as columns
      INTEGER I,J,K,IERR
      DATA I3/1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0/   ! 3x3 identity matrix
      
c...assemble Yukawa matrix in weak eigenbasis
      CALL VEC2MAT(G,INIT,YUK,-1)
      YUKH(1,1)=YUK(1,1)
      YUKH(2,1)=YUK(1,2)
      YUKH(3,1)=YUK(1,3)
      YUKH(1,2)=YUK(2,1)
      YUKH(2,2)=YUK(2,2)
      YUKH(3,2)=YUK(2,3)
      YUKH(1,3)=YUK(3,1)
      YUKH(2,3)=YUK(3,2)
      YUKH(3,3)=YUK(3,3)

c...Compute left rotation matrix VL by diagonalizing  YUK*YUKH     
      CALL MPROD2(A,YUK,YUKH)
      CALL EIGSYS(3,3,A,W,Z,IERR,WORK)
      IF (IERR.NE.0) THEN
        WRITE(LOUT,*) 'EISRS1 ERROR IN YUKDIAG, IERR=',IERR
        STOP99
      END IF
      DO I=1,3
        DO J=1,3
	  VL(I,J)=Z(J,I)
        ENDDO
      ENDDO
c...check if VL has all diagonal elements positive
      DO I=1,3
       IF (VL(I,I).LT.0.d0) THEN      
	 DO J=1,3
	  VL(I,J)=-VL(I,J)   ! sign change in the row
	 ENDDO
       ENDIF
      ENDDO
c...Compute right rotation matrix VR by diagonalizing  YUKH*YUK     
      CALL MPROD2(A,YUKH,YUK)
      CALL EIGSYS(3,3,A,W,Z,IERR,WORK)
      IF (IERR.NE.0) THEN
        WRITE(LOUT,*) 'EISRS1 ERROR IN YUKDIAG, IERR=',IERR
        STOP99
      END IF
      DO I=1,3
        DO J=1,3
	  VR(I,J)=Z(J,I)
        ENDDO
      ENDDO
c...check if VR has all diagonal elements positive
      DO I=1,3
       IF (VR(I,I).LT.0.d0) THEN      
	 DO J=1,3
	  VR(I,J)=-VR(I,J)   ! sign change in the row
	 ENDDO
       ENDIF
      ENDDO
c...Build diagonal Yukawa matrix
      DO I=1,3
        DO J=1,3
	  YDIAG(I,J)=SQRT(W(I))*I3(I,J)
        ENDDO
      ENDDO

      RETURN
      END
