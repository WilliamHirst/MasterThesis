C----------------------------------------------------------------------
      SUBROUTINE TACTIV(IQROT,FTMT,G)
C----------------------------------------------------------------------
C
C    Activates 3rd generation up quark in Yukawa matrix and adjusts
C    RGE-evolving vector G.
C
C    Created: 5/24/07 by Azar Mustafayev 
C
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IQROT
      REAL FTMT
      REAL*8 G(157)
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
      REAL*8 YU(3,3),YUH(3,3),xYU(3,3),WORK(3),A(3,3),W(3),Z(3,3),TEMP,
     &       ZMAX,I3(3,3),VL(3,3),VR(3,3)
c  A(3,3)  - axiliary matrix to be diagonalized
c  W(3)    - vector of eigenvalues for AR
c  Z(3,3)  - matrix of eigenvalues of AR as columns
      INTEGER I,J,K,L,IERR
      DATA I3/1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0/   ! 3x3 identity matrix
      
      IF (IQROT.EQ.0.OR.IQROT.EQ.1) THEN
        G(12)=DBLE(FTMT)
      ELSE 
c...assemble up Yukawa matrix
        CALL VEC2MAT(G,4,YU,-1)
        DO I=1,3
          DO J=1,3
	    YUH(I,J)=YU(J,I)
          ENDDO
        ENDDO
 
c...Compute left rotation matrix VL by diagonalizing  YUK*YUKH     
        DO I=1,3
        DO J=1,3
	  A(I,J)=0.d0
          DO K=1,3
	    A(I,J)=A(I,J)+YU(I,K)*YUH(K,J)
          ENDDO
        ENDDO
        ENDDO
        CALL EIGSYS(3,3,A,W,Z,IERR,WORK)
        IF (IERR.NE.0) THEN
          WRITE(LOUT,*) 'EISRS1 ERROR IN TACTIV, IERR=',IERR
          STOP99
        ENDIF
c...Restore generation structure in W and Z
        DO I=1,3
	  DO J=I,3
	    IF (ABS(Z(I,I)).LT.ABS(Z(I,J)).AND.J.NE.I) THEN
	      DO K=1,3
		TEMP=Z(K,J)	! swap columns
		Z(K,J)=Z(K,I)
		Z(K,I)=TEMP
	      ENDDO
	      TEMP=W(J)        ! swap eigenvalues
	      W(J)=W(I)
	      W(I)=TEMP
	    ENDIF
	  ENDDO
        ENDDO
c...save left rotation matrix
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
        DO I=1,3
          DO J=1,3
	    A(I,J)=0.d0
            DO K=1,3
	      A(I,J)=A(I,J)+YUH(I,K)*YU(K,J)
            ENDDO
          ENDDO
        ENDDO
        CALL EIGSYS(3,3,A,W,Z,IERR,WORK)
        IF (IERR.NE.0) THEN
          WRITE(LOUT,*) 'EISRS1 ERROR IN TACTIV, IERR=',IERR
          STOP99
        ENDIF
c...Restore generation structure in W and Z
        DO I=1,3
	  DO J=I,3
	    IF (ABS(Z(I,I)).LT.ABS(Z(I,J)).AND.J.NE.I) THEN
	      DO K=1,3
		TEMP=Z(K,J)	! swap columns
		Z(K,J)=Z(K,I)
		Z(K,I)=TEMP
	      ENDDO
	      TEMP=W(J)        ! swap eigenvalues
	      W(J)=W(I)
	      W(I)=TEMP
	    ENDIF
	  ENDDO
        ENDDO
c...save right rotation matrix
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

c...build diagonal Yukawa matrix
        DO I=1,3
          DO J=1,3
	    xYU(I,J)=SQRT(DABS(W(I)))*I3(I,J)
          ENDDO
        ENDDO

c...make the adjustment 
        xYU(3,3)=xYU(3,3)+DBLE(FTMT)
c...rotate back to original basis
        DO I=1,3
        DO J=1,3
          YU(I,J)=0.d0
          DO K=1,3
            DO L=1,3
              YU(I,J)=YU(I,J)+VL(K,I)*xYU(K,L)*VR(L,J)
            ENDDO
          ENDDO
        ENDDO
        ENDDO

        CALL MAT2VEC(G,4,YU,-1)
      ENDIF
      
      
      RETURN
      END
