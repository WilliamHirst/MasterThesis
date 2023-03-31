!
      SUBROUTINE ORTH(VA,VE,QBLARG,QBLARGVAL,QBROT,P)
!
!Purpose: To make sure that the eigenvectors are orthogonal in the
!         case that the eigenvectors are close to one another.
!         Resets one entry in one eigenvector. Chooses entry
!         to reset by finding the largest multiple in the
!         orthogonality sum.
!         QBLARG is the location of the largest off-diagonal entry
!                of M in quark basis.
!         QBLARGVAL is the value of the largest off-diagonal entry
!         QBROT is rotation from current to quark basis.
!
!The diagonalisation routine is not perfectly accurate. Just as the
!Vs come with intrinsic error in unitarity of 10^-16, the diagonalisation
!has an intrinsic deviation from unitarity, which is particulary
!acute if the eigenvectors are close, when the error can be
!as large as 10^-11. This is a problem when we are trying to calculate
!small numbers, so we must shift the error into the most convenient
!location. ORTH will use information about the eigenvectors, eigenvalues,
!and current->quark mass basis rotation to shift this error so that it
!affects the largest of the off-diagonal entries of the squark mass matrix.
!By shifting the error to affect the largest off-diagonal entry, we ensure
!that it has the smallest possible effect. For example, if we use mSUGRA
!and m_0=300, we find that the effect on the off-diagonal entries is
!of the order of 10^-6 (~10^-11*300^2). If the off-diagonal entry in the
!squared squark mass matrix is larger than 10^-6, the numerical fluctuations
!are swamped by the actual value. See Appendix A of PRD 79, (2009) 035004
!
      IMPLICIT NONE
!
      COMMON/ORTHWARNING/ORTHFLAG
      INTEGER ORTHFLAG
      SAVE/ORTHWARNING/
!
      DOUBLE COMPLEX VE(3,3),VA(3),QBROT(3,3),COMBROT(3,3)
      DOUBLE COMPLEX SUM,CMATMUL,QBLARGVAL,MAXCOMBOFF
      DOUBLE PRECISION MULT(3)
      INTEGER I,J,K,A,B,QBLARG(2),FIXE(2),FIXQB(2),P
!
      IF(P.EQ.1)WRITE(*,*)'BEFORE'
      DO I=1,3
        DO J=1,3
          COMBROT(I,J)=CMATMUL(1,VE,QBROT,I,J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          IF(P.EQ.1)WRITE(*,*)I,J,ABS(CMATMUL(1,COMBROT,COMBROT,I,J))
        END DO
      END DO
      IF(P.EQ.1)WRITE(*,*)
!
!Find two closest eigenvalues
!
      DO I=1,2
        IF(ABS(VA(1)-VA(2)).LT.ABS(VA(1)-VA(3)).AND.
     $     ABS(VA(1)-VA(2)).LT.ABS(VA(2)-VA(3)))THEN
          FIXE(1)=1
          FIXE(2)=2
        ELSE IF(ABS(VA(2)-VA(3)).LT.ABS(VA(1)-VA(3)))THEN
          FIXE(1)=2
          FIXE(2)=3
        ELSE
          FIXE(1)=1
          FIXE(2)=3
        END IF
      END DO
!
!Find the diagonal quark basis entry which corresponds to these 
!eigenvalues.
!The first index on combrot is the squark mass basis and the second
!is the quark mass basis.
!
      DO J=1,2
        IF(ABS(COMBROT(FIXE(J),1)).GT.ABS(COMBROT(FIXE(J),2)).AND.
     $     ABS(COMBROT(FIXE(J),1)).GT.ABS(COMBROT(FIXE(J),3)))THEN
          FIXQB(J)=1
        ELSE IF(ABS(COMBROT(FIXE(J),2)).GT.ABS(COMBROT(FIXE(J),3)))THEN
          FIXQB(J)=2
        ELSE
          FIXQB(J)=3
        END IF
      END DO
      IF(FIXQB(1).EQ.FIXQB(2))THEN
        WRITE(*,*)'ERROR IN FIXQB IN ORTH - EIGENVALUES ARE'
        WRITE(*,*)'PROBABLY VERY CLOSE. CHECK THE SQUARK'
        WRITE(*,*)'DIAGONALIZATION'
!
!Note: An error here is probably due combrot being far away from
!      the unit matrix. In this case, there is no clear separation
!      of entries in one basis to be identified with the eigenvalues.
!      This will probably not be a problem, since in this case it
!      is unlikely that the eigenvalues will be particularly close
!      and the fix in this subroutine will not be necessary.
!
        RETURN
!
      END IF
!
!Now fix the eigenvectors specified by FIXE to be orthogonal by
!changing one entry in the eigenvector. The eigenvector to alter
!corresponds to the contributor in the largest off-diagonal entry
!in the quark mass basis.
!
      IF((QBLARG(1)+QBLARG(2)).EQ.(FIXQB(1)+FIXQB(2)))THEN
!Don't fix anything if the error is already in the largest
!off-diagonal element
        A=0
        B=0
      ELSE IF(QBLARG(1).EQ.FIXQB(1).OR.QBLARG(2).EQ.FIXQB(1))THEN
!Fix the eigenvector corresponding to FIXQB(1)
        A=1 !A is entry of FIXE that is being fixed
        B=2
      ELSE
!Fix the eigenvector corresponding to FIXQB(2)
        A=2
        B=1
      END IF
!
      IF(A.NE.B)THEN
        DO K=1,3
          MULT(K)=ABS(CONJG(VE(K,FIXE(B)))*VE(K,FIXE(A)))
        END DO
        IF(MULT(1).GT.MULT(2).AND.MULT(1).GT.MULT(3))THEN
          VE(1,FIXE(A))=-(VE(2,FIXE(A))*CONJG(VE(2,FIXE(B)))
     $                   +VE(3,FIXE(A))*CONJG(VE(3,FIXE(B))))
     $                                          /CONJG(VE(1,FIXE(B)))
        ELSE IF(MULT(2).GT.MULT(1).AND.MULT(2).GT.MULT(3))THEN
          VE(2,FIXE(A))=-(VE(1,FIXE(A))*CONJG(VE(1,FIXE(B)))
     $                   +VE(3,FIXE(A))*CONJG(VE(3,FIXE(B))))
     $                                          /CONJG(VE(2,FIXE(B)))
        ELSE IF(MULT(3).GT.MULT(1).AND.MULT(3).GT.MULT(2))THEN
          VE(3,FIXE(A))=-(VE(1,FIXE(A))*CONJG(VE(1,FIXE(B)))
     $                   +VE(2,FIXE(A))*CONJG(VE(2,FIXE(B))))
     $                                          /CONJG(VE(3,FIXE(B)))
        ELSE IF(MULT(1).NE.MULT(2).AND.MULT(1).NE.MULT(3)
     $                            .AND.MULT(2).NE.MULT(3))THEN
          WRITE(*,*)'ERROR IN MULT OF ORTH'
        END IF 
      END IF
      IF(P.EQ.1)WRITE(*,*)'AFTER'
      DO I=1,3
        DO J=1,3
          COMBROT(I,J)=CMATMUL(1,VE,QBROT,I,J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          IF(P.EQ.1)WRITE(*,*)I,J,ABS(CMATMUL(1,COMBROT,COMBROT,I,J))
        END DO
      END DO
!
!Find largest off-diagonal value of COMBROT^dagger x COMBROT for
!for error checking purposes.
!
      MAXCOMBOFF=(0.D0,0.D0)
      DO I=1,3
        DO J=1,3
          IF(I.NE.J)THEN
            IF(ABS(MAXCOMBOFF).LT.ABS(CMATMUL(1,COMBROT,COMBROT,I,J)))
     $                         MAXCOMBOFF=CMATMUL(1,COMBROT,COMBROT,I,J)
          END IF
        END DO
      END DO
!
!If the smallest eigenvalue multiplied by the largest off-diagonal
!entry in COMBROT is greater than the largest off-diagonal entry of
!the squark mass matrix in chosen quark 'mass' basis then add 1 to
!ORTHFLAG. See the main routine for RGEFLAV.
!
      IF(ABS(VA(1)*MAXCOMBOFF).GT.ABS(QBLARGVAL))ORTHFLAG=ORTHFLAG+1
!
      RETURN
      END
