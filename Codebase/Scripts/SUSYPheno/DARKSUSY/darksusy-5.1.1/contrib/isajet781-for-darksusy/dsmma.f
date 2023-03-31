!
      SUBROUTINE DSMMA(G,Q,SVLQ)
!
!Purpose: Find the down squark mass matrix, sort the eigenvectors and
!         eigenvalues and then (if required) write them to a file
!         called sqm2d.dat
!
      IMPLICIT NONE
!
      COMMON /RGEFNM/ FNRGE
      CHARACTER*128 FNRGE,STRADD,SQM2D
      DOUBLE COMPLEX B1EVE(6,6),B1EVA(6),DSQM(6,6),COSTHB,SINTHB
      DOUBLE PRECISION MSB1,MSB2
      DOUBLE PRECISION SUM,Q
      DOUBLE COMPLEX G(601),CDSQM(6,6),TEMP(7)
      DOUBLE COMPLEX EVERTMP(6,6),CWORK1(99),CWORK2(12)
      INTEGER I,J,K,CIERR,SVLQ
!
      DO I=1,3
        DO J=1,3
          B1EVE(I,J)=(0.D0,0.D0)
        END DO
        B1EVA(I)=(0.D0,0.D0)
      END DO
      COSTHB=(0.D0,0.D0)
      SINTHB=(0.D0,0.D0)
      MSB1=0.D0
      MSB2=0.D0
!
      CALL DOWNSQM(G,Q,DSQM)
      DO I=1,6
        DO J=1,6
          CDSQM(I,J)=DSQM(I,J)
        END DO
      END DO
      SQM2D=STRADD(FNRGE,'.sqm2d')
      OPEN(25,FILE=SQM2D,STATUS='UNKNOWN',FORM='FORMATTED')
!
      WRITE(25,*)
      WRITE(25,*)'THE 6X6 DOWN-TYPE MATRIX AT THE SCALE: ',Q
      WRITE(25,*)'IN THE CHOSEN BASIS, IS:'
!
      WRITE(25,*)
      WRITE(25,15)'d_l','s_l','b_l','d_r','s_r','b_r'
      WRITE(25,*)
!
      WRITE(25,11)'d_l',DSQM(1,1),DSQM(1,2),DSQM(1,3),
     $            DSQM(1,4),DSQM(1,5),DSQM(1,6)
      WRITE(25,11)'s_l',DSQM(2,1),DSQM(2,2),DSQM(2,3),
     $            DSQM(2,4),DSQM(2,5),DSQM(2,6)
      WRITE(25,11)'b_l',DSQM(3,1),DSQM(3,2),DSQM(3,3),
     $            DSQM(3,4),DSQM(3,5),DSQM(3,6)
      WRITE(25,11)'d_r',DSQM(4,1),DSQM(4,2),DSQM(4,3),
     $            DSQM(4,4),DSQM(4,5),DSQM(4,6)
      WRITE(25,11)'s_r',DSQM(5,1),DSQM(5,2),DSQM(5,3),
     $            DSQM(5,4),DSQM(5,5),DSQM(5,6)
      WRITE(25,11)'b_r',DSQM(6,1),DSQM(6,2),DSQM(6,3),
     $            DSQM(6,4),DSQM(6,5),DSQM(6,6)
!
      WRITE(25,*)
!
      IF(SVLQ.EQ.0)THEN
!
!Diagonalise the matrix using the LAPACK routine ZGEEV.
!
        CALL ZGEEV('V','V',6,CDSQM,6,B1EVA,B1EVE,6,EVERTMP,6,CWORK1,99
     $             ,CWORK2,CIERR)
!
        DO J=1,6
          SUM=0.D0
          DO I=1,6
            SUM=SUM+DBLE(B1EVE(I,J)*CONJG(B1EVE(I,J)))
          END DO
          DO I=1,6
            B1EVE(I,J)=B1EVE(I,J)/DSQRT(SUM)
          END DO
        END DO
!
!Find the eigenvectors with the most amount of sbottom
!
        DO I=1,5
          DO J=I+1,6
            IF(DSQRT(DBLE(B1EVE(3,I)*CONJG(B1EVE(3,I)))+DBLE(B1EVE(6,I)
     $              *CONJG(B1EVE(6,I)))).LT.
     $         DSQRT(DBLE(B1EVE(3,J)*CONJG(B1EVE(3,J)))+DBLE(B1EVE(6,J)
     $              *CONJG(B1EVE(6,J)))))THEN
              DO K=1,6
                TEMP(K)=B1EVE(K,I)
                B1EVE(K,I)=B1EVE(K,J)
                B1EVE(K,J)=TEMP(K)
              END DO
              TEMP(7)=B1EVA(I)
              B1EVA(I)=B1EVA(J)
              B1EVA(J)=TEMP(7)
            END IF
          END DO
        END DO
!
!Now find next two with most amount of scalar strange
!
        DO I=3,5
          DO J=I+1,6
            IF(DSQRT(DBLE(B1EVE(2,I)*CONJG(B1EVE(2,I)))+DBLE(B1EVE(5,I)
     $              *CONJG(B1EVE(5,I)))).LT.
     $         DSQRT(DBLE(B1EVE(2,J)*CONJG(B1EVE(2,J)))+DBLE(B1EVE(5,J)
     $              *CONJG(B1EVE(5,J)))))THEN
              DO K=1,6
                TEMP(K)=B1EVE(K,I)
                B1EVE(K,I)=B1EVE(K,J)
                B1EVE(K,J)=TEMP(K)
              END DO
              TEMP(7)=B1EVA(I)
              B1EVA(I)=B1EVA(J)
              B1EVA(J)=TEMP(7)
            END IF
          END DO
        END DO
!
!Sort the first two eigenvectors by mass
!
        IF(ABS(B1EVA(1)).GT.ABS(B1EVA(2)))THEN
          DO I=1,6
            TEMP(I)=B1EVE(I,1)
            B1EVE(I,1)=B1EVE(I,2)
            B1EVE(I,2)=TEMP(I)
          END DO
          TEMP(7)=B1EVA(1)
          B1EVA(1)=B1EVA(2)
          B1EVA(2)=TEMP(7)
        END IF
!
!Sort the next two by mass
!
        IF(ABS(B1EVA(3)).GT.ABS(B1EVA(4)))THEN
          DO I=1,6
            TEMP(I)=B1EVE(I,3)
            B1EVE(I,3)=B1EVE(I,4)
            B1EVE(I,4)=TEMP(I)
          END DO
          TEMP(7)=B1EVA(3)
          B1EVA(3)=B1EVA(4)
          B1EVA(4)=TEMP(7)
        END IF
!
!And the final two
!
        IF(ABS(B1EVA(5)).GT.ABS(B1EVA(6)))THEN
          DO I=1,6
            TEMP(I)=B1EVE(I,5)
            B1EVE(I,5)=B1EVE(I,6)
            B1EVE(I,6)=TEMP(I)
          END DO
          TEMP(7)=B1EVA(5)
          B1EVA(5)=B1EVA(6)
          B1EVA(6)=TEMP(7)
        END IF
!
!Note about sin{\theta_b}: If I define sin{\theta_b} the same way as
!BT, then -sin{\theta_b}=conjg(b1eve(6,1))
!
        COSTHB=CONJG(B1EVE(3,1))
        SINTHB=-CONJG(B1EVE(6,1))
        MSB1=DSQRT(ABS(B1EVA(1)))
        MSB2=DSQRT(ABS(B1EVA(2)))
!
        WRITE(25,*)'EIGENVECTORS ARE:'
        WRITE(25,*)
        WRITE(25,15)'b_1','b_2','s_1','s_2','d_1','d_2'
        WRITE(25,*)
!     
        WRITE(25,11)'d_l',B1EVE(1,1),B1EVE(1,2),B1EVE(1,3),
     $              B1EVE(1,4),B1EVE(1,5),B1EVE(1,6)
        WRITE(25,11)'s_l',B1EVE(2,1),B1EVE(2,2),B1EVE(2,3),
     $              B1EVE(2,4),B1EVE(2,5),B1EVE(2,6)
        WRITE(25,11)'b_l',B1EVE(3,1),B1EVE(3,2),B1EVE(3,3),
     $              B1EVE(3,4),B1EVE(3,5),B1EVE(3,6)
        WRITE(25,11)'d_r',B1EVE(4,1),B1EVE(4,2),B1EVE(4,3),
     $              B1EVE(4,4),B1EVE(4,5),B1EVE(4,6)
        WRITE(25,11)'s_r',B1EVE(5,1),B1EVE(5,2),B1EVE(5,3),
     $              B1EVE(5,4),B1EVE(5,5),B1EVE(5,6)
        WRITE(25,11)'b_r',B1EVE(6,1),B1EVE(6,2),B1EVE(6,3),
     $              B1EVE(6,4),B1EVE(6,5),B1EVE(6,6)
!
        WRITE(25,*)
        WRITE(25,*)'EIGENVALUES ARE:'
        WRITE(25,16)B1EVA(1),B1EVA(2),B1EVA(3),B1EVA(4),B1EVA(5)
     $             ,B1EVA(6)
        WRITE(25,*)
        WRITE(25,*)'MASSES ARE:'
        WRITE(25,17)DSQRT(DSQRT(DBLE(B1EVA(1)*CONJG(B1EVA(1)))))
     $             ,DSQRT(DSQRT(DBLE(B1EVA(2)*CONJG(B1EVA(2)))))
     $             ,DSQRT(DSQRT(DBLE(B1EVA(3)*CONJG(B1EVA(3)))))
     $             ,DSQRT(DSQRT(DBLE(B1EVA(4)*CONJG(B1EVA(4)))))
     $             ,DSQRT(DSQRT(DBLE(B1EVA(5)*CONJG(B1EVA(5)))))
     $             ,DSQRT(DSQRT(DBLE(B1EVA(6)*CONJG(B1EVA(6)))))
      END IF
!
!
 11   FORMAT(SP,1P,A,3X,'(',D10.3,',',D10.3,')',1X,'(',D10.3,',',
     $       D10.3,')',1X,'(',D10.3,',',D10.3,')',1X,'(',
     $       D10.3,',',D10.3,')',1X,'(',D10.3,',',D10.3,')'
     $      ,1X,'(',D10.3,',',D10.3,')')
 15   FORMAT(15X,A,19X,A,19X,A,19X,A,19X,A,19X,A)
 16   FORMAT(SP,1P,6X,'(',D10.3,',',D10.3,')',1X,'(',D10.3,',',
     $       D10.3,')',1X,'(',D10.3,',',D10.3,')',1X,'(',
     $       D10.3,',',D10.3,')',1X,'(',D10.3,',',D10.3,')'
     $      ,1X,'(',D10.3,',',D10.3,')')
 17   FORMAT(10X,F11.4,11X,F11.4,11X,F11.4,11X,F11.4,11X
     $      ,F11.4,11X,F11.4)
      WRITE(25,*)
      CLOSE(25)
!
      RETURN
      END
