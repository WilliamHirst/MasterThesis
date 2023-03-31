C----------------------------------------------------------------------
      SUBROUTINE AYUKDIAG(G,INIT,YDIAG,VL,VR)
C----------------------------------------------------------------------
C
C  This is inversed of YUKDIAG.
c  It rotates diagonal Yukawa matrix to weak eigenbasis and 
c  "diassembles" it filling corresponding elements of G(i).
c
C   NOTE: In our notation fermion masses are proportional to transposed
C         Yukawas in gauge eigenbasis, i.e.  m ~ VR Y^T VL^dagger,
C         while SURG111 follows m ~ Y convention.
C         One has to transpose Yukasas when passing
C         into or from those subroutines.
C
c  Ref: Misak, Pokorski & Rosiek  hep-ph/9703442
c
c  Created: 02/13/07 by Azar Mustafayev.
C 
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C   Input:
c    VL,VR - Left and right rotation matrices: YDIAG=VL^* YUK VR^T
C    YDIAG - diagonal Yukawa matrix
c   Output:
C     G(4) = Y_u(1,1)
C     G(5) = Y_u(1,2)
C     G(6) = Y_u(1,3)
C     ...    ...
C     G(12) = Y_u(3,3)
C     G(13)-G(21) = Y_d
C     G(22)-G(30) = Y_e
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv            
      IMPLICIT NONE
      REAL*8 G(157),YDIAG(3,3),VL(3,3),VR(3,3)
      INTEGER INIT
c      
      REAL*8 YUK(3,3),VLT(3,3)

      VLT(1,1)=VL(1,1)
      VLT(2,1)=VL(1,2)
      VLT(3,1)=VL(1,3)
      VLT(1,2)=VL(2,1)
      VLT(2,2)=VL(2,2)
      VLT(3,2)=VL(2,3)
      VLT(1,3)=VL(3,1)
      VLT(2,3)=VL(3,2)
      VLT(3,3)=VL(3,3)
      
      CALL MPROD3(YUK,VLT,YDIAG,VR)
 
      CALL MAT2VEC(G,INIT,YUK,-1)

      RETURN
      END
