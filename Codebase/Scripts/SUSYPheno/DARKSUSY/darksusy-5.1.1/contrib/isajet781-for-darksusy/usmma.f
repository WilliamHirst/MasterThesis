!
      SUBROUTINE USMMA(G,Q,SVLQ)
!
!Purpose: Find the up squark mass matrix, sort the eigenvectors and
!         eigenvalues and then (if required) write them to a file
!         called sqm2u.dat
!
      IMPLICIT NONE
!
      COMMON /RGEFNM/ FNRGE
      CHARACTER*128 FNRGE,STRADD,SQM2U
      COMMON/DECCALC/T1EVE,T1EVA,USQM,COSTHT,SINTHT,GHIK,MST1,MST2,GAMMA
      DOUBLE COMPLEX T1EVE(6,6),T1EVA(6),USQM(6,6),COSTHT,SINTHT
     $              ,GHIK(601)
      DOUBLE PRECISION MST1,MST2,GAMMA
      SAVE/DECCALC/
!
      DOUBLE PRECISION SUM,Q
      DOUBLE COMPLEX G(601),CUSQM(6,6),TEMP(7)
      DOUBLE COMPLEX EVERTMP(6,6),CWORK1(99),CWORK2(12)
      INTEGER I,J,K,CIERR,SVLQ
!
      DO I=1,3
        DO J=1,3
          T1EVE(I,J)=(0.D0,0.D0)
        END DO
        T1EVA(I)=(0.D0,0.D0)
      END DO
      COSTHT=(0.D0,0.D0)
      SINTHT=(0.D0,0.D0)
      MST1=0.D0
      MST2=0.D0
!
      CALL UPSQM(G,Q,USQM)
      DO I=1,6
        DO J=1,6
          CUSQM(I,J)=USQM(I,J)
        END DO
      END DO
      SQM2U=STRADD(FNRGE,'.sqm2u')
      OPEN(25,FILE=SQM2U,STATUS='UNKNOWN',FORM='FORMATTED')
!
      WRITE(25,*)
      WRITE(25,*)'THE 6X6 UP-TYPE MATRIX AT THE SCALE: ',Q
      WRITE(25,*)'IN THE CHOSEN BASIS, IS:'
!
      WRITE(25,*)
      WRITE(25,15)'u_l','c_l','t_l','u_r','c_r','t_r'
      WRITE(25,*)
!
      WRITE(25,11)'u_l',USQM(1,1),USQM(1,2),USQM(1,3),
     $            USQM(1,4),USQM(1,5),USQM(1,6)
      WRITE(25,11)'c_l',USQM(2,1),USQM(2,2),USQM(2,3),
     $            USQM(2,4),USQM(2,5),USQM(2,6)
      WRITE(25,11)'t_l',USQM(3,1),USQM(3,2),USQM(3,3),
     $            USQM(3,4),USQM(3,5),USQM(3,6)
      WRITE(25,11)'u_r',USQM(4,1),USQM(4,2),USQM(4,3),
     $            USQM(4,4),USQM(4,5),USQM(4,6)
      WRITE(25,11)'c_r',USQM(5,1),USQM(5,2),USQM(5,3),
     $            USQM(5,4),USQM(5,5),USQM(5,6)
      WRITE(25,11)'t_r',USQM(6,1),USQM(6,2),USQM(6,3),
     $            USQM(6,4),USQM(6,5),USQM(6,6)
!
      WRITE(25,*)
!
      IF(SVLQ.EQ.1)THEN
!
!Diagonalise the matrix using the LAPACK routine ZGEEV.
!
        CALL ZGEEV('V','V',6,CUSQM,6,T1EVA,T1EVE,6,EVERTMP,6,CWORK1,99
     $             ,CWORK2,CIERR)
!
        DO J=1,6
          SUM=0.D0
          DO I=1,6
            SUM=SUM+DBLE(T1EVE(I,J)*CONJG(T1EVE(I,J)))
          END DO
          DO I=1,6
            T1EVE(I,J)=T1EVE(I,J)/DSQRT(SUM)
          END DO
        END DO
!
!Find the eigenvectors with the most amount of stop
!
        DO I=1,5
          DO J=I+1,6
            IF(DSQRT(DBLE(T1EVE(3,I)*CONJG(T1EVE(3,I)))+DBLE(T1EVE(6,I)
     $              *CONJG(T1EVE(6,I)))).LT.
     $         DSQRT(DBLE(T1EVE(3,J)*CONJG(T1EVE(3,J)))+DBLE(T1EVE(6,J)
     $              *CONJG(T1EVE(6,J)))))THEN
              DO K=1,6
                TEMP(K)=T1EVE(K,I)
                T1EVE(K,I)=T1EVE(K,J)
                T1EVE(K,J)=TEMP(K)
              END DO
              TEMP(7)=T1EVA(I)
              T1EVA(I)=T1EVA(J)
              T1EVA(J)=TEMP(7)
            END IF
          END DO
        END DO
!
!Now find next two with most amount of scalar charm
!
        DO I=3,5
          DO J=I+1,6
            IF(DSQRT(DBLE(T1EVE(2,I)*CONJG(T1EVE(2,I)))+DBLE(T1EVE(5,I)
     $              *CONJG(T1EVE(5,I)))).LT.
     $         DSQRT(DBLE(T1EVE(2,J)*CONJG(T1EVE(2,J)))+DBLE(T1EVE(5,J)
     $              *CONJG(T1EVE(5,J)))))THEN
              DO K=1,6
                TEMP(K)=T1EVE(K,I)
                T1EVE(K,I)=T1EVE(K,J)
                T1EVE(K,J)=TEMP(K)
              END DO
              TEMP(7)=T1EVA(I)
              T1EVA(I)=T1EVA(J)
              T1EVA(J)=TEMP(7)
            END IF
          END DO
        END DO
!
!Sort the first two eigenvectors by mass
!
        IF(ABS(T1EVA(1)).GT.ABS(T1EVA(2)))THEN
          DO I=1,6
            TEMP(I)=T1EVE(I,1)
            T1EVE(I,1)=T1EVE(I,2)
            T1EVE(I,2)=TEMP(I)
          END DO
          TEMP(7)=T1EVA(1)
          T1EVA(1)=T1EVA(2)
          T1EVA(2)=TEMP(7)
        END IF
!
!Sort the next two by mass
!
        IF(ABS(T1EVA(3)).GT.ABS(T1EVA(4)))THEN
          DO I=1,6
            TEMP(I)=T1EVE(I,3)
            T1EVE(I,3)=T1EVE(I,4)
            T1EVE(I,4)=TEMP(I)
          END DO
          TEMP(7)=T1EVA(3)
          T1EVA(3)=T1EVA(4)
          T1EVA(4)=TEMP(7)
        END IF
!
!And the final two
!
        IF(ABS(T1EVA(5)).GT.ABS(T1EVA(6)))THEN
          DO I=1,6
            TEMP(I)=T1EVE(I,5)
            T1EVE(I,5)=T1EVE(I,6)
            T1EVE(I,6)=TEMP(I)
          END DO
          TEMP(7)=T1EVA(5)
          T1EVA(5)=T1EVA(6)
          T1EVA(6)=TEMP(7)
        END IF
!
!Note about sin{\theta_t}: If I define sin{\theta_t} the same way as
!BT, then -sin{\theta_t}=conjg(t1eve(6,1))
!
        COSTHT=CONJG(T1EVE(3,1))
        SINTHT=-CONJG(T1EVE(6,1))
        MST1=DSQRT(ABS(T1EVA(1)))
        MST2=DSQRT(ABS(T1EVA(2)))
!
        WRITE(25,*)'EIGENVECTORS ARE:'
        WRITE(25,*)
        WRITE(25,15)'t_1','t_2','c_1','c_2','u_1','u_2'
        WRITE(25,*)
!     
        WRITE(25,11)'u_l',T1EVE(1,1),T1EVE(1,2),T1EVE(1,3),
     $              T1EVE(1,4),T1EVE(1,5),T1EVE(1,6)
        WRITE(25,11)'c_l',T1EVE(2,1),T1EVE(2,2),T1EVE(2,3),
     $              T1EVE(2,4),T1EVE(2,5),T1EVE(2,6)
        WRITE(25,11)'t_l',T1EVE(3,1),T1EVE(3,2),T1EVE(3,3),
     $              T1EVE(3,4),T1EVE(3,5),T1EVE(3,6)
        WRITE(25,11)'u_r',T1EVE(4,1),T1EVE(4,2),T1EVE(4,3),
     $              T1EVE(4,4),T1EVE(4,5),T1EVE(4,6)
        WRITE(25,11)'c_r',T1EVE(5,1),T1EVE(5,2),T1EVE(5,3),
     $              T1EVE(5,4),T1EVE(5,5),T1EVE(5,6)
        WRITE(25,11)'t_r',T1EVE(6,1),T1EVE(6,2),T1EVE(6,3),
     $              T1EVE(6,4),T1EVE(6,5),T1EVE(6,6)
!
        WRITE(25,*)
        WRITE(25,*)'EIGENVALUES ARE:'
        WRITE(25,16)T1EVA(1),T1EVA(2),T1EVA(3),T1EVA(4),T1EVA(5)
     $             ,T1EVA(6)
        WRITE(25,*)
        WRITE(25,*)'MASSES ARE:'
        WRITE(25,17)DSQRT(DSQRT(DBLE(T1EVA(1)*CONJG(T1EVA(1)))))
     $             ,DSQRT(DSQRT(DBLE(T1EVA(2)*CONJG(T1EVA(2)))))
     $             ,DSQRT(DSQRT(DBLE(T1EVA(3)*CONJG(T1EVA(3)))))
     $             ,DSQRT(DSQRT(DBLE(T1EVA(4)*CONJG(T1EVA(4)))))
     $             ,DSQRT(DSQRT(DBLE(T1EVA(5)*CONJG(T1EVA(5)))))
     $             ,DSQRT(DSQRT(DBLE(T1EVA(6)*CONJG(T1EVA(6)))))
      END IF
!
!
 11   FORMAT(SP,1P,A,3X,'(',D10.3,',',D10.3,')',1X,'(',D10.3,',',
     $       D10.3,')',1X,'(',D10.3,',',D10.3,')',1X,'(',
     $       D10.3,',',D10.3,')',1X,'(',D10.3,',',D10.3,')'
     $      ,1X,'(',D10.3,',',D10.3,')')
 15   FORMAT(16X,A,21X,A,21X,A,21X,A,21X,A,21X,A)
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
