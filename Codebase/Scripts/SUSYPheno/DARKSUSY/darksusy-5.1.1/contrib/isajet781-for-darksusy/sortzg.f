!
      SUBROUTINE SORTZG(VA,VE,OLDNACTIVE)
!
!Purpose: To sort the eigenvectors into order of increasing
!         eigenvalue after normalising them.
!         If the heaviest sfermion has decoupled, place
!         the smallest eigenvalue in the place where the largest
!         would be.
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX VA(3),VE(3,3),SWAPVA,SWAPVE(3)
      DOUBLE PRECISION SUM,REZ,IMZ,THETA
      INTEGER I,J,K,OLDNACTIVE
!
      DO J=1,3
        SUM=0.D0
        DO I=1,3
          SUM=SUM+DBLE(VE(I,J)*CONJG(VE(I,J)))
        END DO
        DO I=1,3
          VE(I,J)=VE(I,J)/DSQRT(SUM)
        END DO
      END DO
!
!Fix the largest entry to be real and positive.
!
      DO I=1,3
        IF(ABS(VE(1,I)).GT.ABS(VE(2,I))
     $              .AND.ABS(VE(1,I)).GT.ABS(VE(3,I)))THEN
          REZ=DBLE(VE(1,I))
          IMZ=DIMAG(VE(1,I))
        ELSE IF(ABS(VE(2,I)).GT.ABS(VE(3,I)))THEN
          REZ=DBLE(VE(2,I))
          IMZ=DIMAG(VE(2,I))
        ELSE
          REZ=DBLE(VE(3,I))
          IMZ=DIMAG(VE(3,I))
        END IF
        THETA=DATAN2(IMZ,REZ)
        DO J=1,3
          VE(J,I)=VE(J,I)*DCMPLX(DCOS(THETA),-DSIN(THETA))
        END DO
      END DO
!
!Now sort the eigenvectors
!
      DO I=1,2
        DO J=I+1,3
          IF(DBLE(VA(I)).GT.DBLE(VA(J)))THEN
            DO K=1,3
              SWAPVE(K)=VE(K,J)
              VE(K,J)=VE(K,I)
              VE(K,I)=SWAPVE(K)
            END DO
            SWAPVA=VA(J)
            VA(J)=VA(I)
            VA(I)=SWAPVA
          END IF
        END DO
      END DO
!
!Redo the sorting if we do not have 3 active sfermions
!Need to be careful about negative eigenvalues, which are
!OK as long as the value of the eigenvalue is positive
!at the weak scale.
!
      IF(OLDNACTIVE.EQ.2)THEN
        IF(DBLE(VA(1)).LT.0.AND.DBLE(VA(2)).GE.0)THEN
          DO K=1,3
            SWAPVE(K)=VE(K,2)
            VE(K,2)=VE(K,3)
            VE(K,3)=SWAPVE(K)
          END DO
          SWAPVA=VA(2)
          VA(2)=VA(3)
          VA(3)=SWAPVA
        ELSE IF(DBLE(VA(1)).GE.0)THEN
          DO I=2,1,-1
            DO K=1,3
              SWAPVE(K)=VE(K,I)
              VE(K,I)=VE(K,3)
              VE(K,3)=SWAPVE(K)
            END DO
            SWAPVA=VA(I)
            VA(I)=VA(3)
            VA(3)=SWAPVA
          END DO
        END IF
      END IF
      IF(OLDNACTIVE.EQ.1.AND.DBLE(VA(1)).GE.0)THEN
        DO I=2,1,-1
          DO K=1,3
            SWAPVE(K)=VE(K,I)
            VE(K,I)=VE(K,I+1)
            VE(K,I+1)=SWAPVE(K)
          END DO
          SWAPVA=VA(I)
          VA(I)=VA(I+1)
          VA(I+1)=SWAPVA
        END DO
      END IF
!
      RETURN
      END
