CDECK  ID>, RPMAIN
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C-----------------------------------------------------------------------
      SUBROUTINE RPMAIN 
C-----------------------------------------------------------------------
C     MAIN R-PARITY SUBROUTINE
C-----------------------------------------------------------------------
C     RPARTY   = logical, .TRUE. is consereved/.FALSE. is violated
C     LAMDA1   = LLE couplings
C     LAMDA2   = LUD couplings
C     LAMDA3   = UDD couplings
      COMMON/RSLASH/LAMDA1(3,3,3),LAMDA2(3,3,3),LAMDA3(3,3,3),RPARTY
      LOGICAL RPARTY
      REAL LAMDA1,LAMDA2,LAMDA3
      SAVE /RSLASH/
      CHARACTER*2 RYORN
C--Decide if want Rparity violating couplings
C-FP  Assume no Rparity violation, just return
      RPARTY=.TRUE.
      RETURN
      PRINT 100
 100  FORMAT(' R PARITY VIOLATION (Y/N)')
      READ(*,'(A1)')RYORN
      RPARTY = RYORN.NE.'Y'.AND.RYORN.NE.'y'
      IF(RPARTY) RETURN
C--If decide want rparity violation input the couplings
      PRINT 150
 150  FORMAT('Do you want lambda couplings(y/n)?')
      READ (*,'(A1)') RYORN
      IF(RYORN.EQ.'Y'.OR.RYORN.EQ.'y') THEN
        PRINT 200
 200    FORMAT(
     &      'ENTER LAMDBA IN ORDER 121,122,123,131,132,133,231,232,233')
        READ*,LAMDA1(1,2,1),LAMDA1(1,2,2),LAMDA1(1,2,3),LAMDA1(1,3,1),
     &        LAMDA1(1,3,2),LAMDA1(1,3,3),LAMDA1(2,3,1),LAMDA1(2,3,2),
     &        LAMDA1(2,3,3)
C--use the antisymmetry to find the rest
        DO K=1,3
          DO I=1,3
            DO J=1,3
              IF(I.EQ.J) LAMDA1(I,J,K) = 0.0E0
              IF(I.GT.J) LAMDA1(I,J,K) = -LAMDA1(J,I,K)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO I=1,3
          DO J=1,3
            DO K=1,3
              LAMDA1(I,J,K) = 0.
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      PRINT 250
 250  FORMAT('Do you want lambda` couplings(y/n)?')
      READ(*,'(A1)') RYORN
      IF(RYORN.EQ.'Y'.OR.RYORN.EQ.'y') THEN
        PRINT 300
        PRINT 310
        PRINT 320
        PRINT 330
 300    FORMAT('ENTER LAMBDA` IN THE ORDER')
 310    FORMAT('111,112,113,121,122,123,131,132,133,')
 320    FORMAT('211,212,213,221,222,223,231,232,233,')
 330    FORMAT('311,312,313,321,322,323,331,332,333,')
        READ*,(((LAMDA2(I,J,K) ,K=1,3),J=1,3),I=1,3)
      ELSE
        DO I=1,3
          DO J=1,3
            DO K=1,3
              LAMDA2(I,J,K) = 0.
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      PRINT 350
 350  FORMAT('Do you want lambda`` couplings(y/n)?')
      READ (*,'(A1)') RYORN
      IF(RYORN.EQ.'Y'.OR.RYORN.EQ.'y') THEN
        PRINT 400
 400    FORMAT(
     &      'ENTER LAMDBA`` ORDER 112,113,123,212,213,223,312,313,323')
        READ*,LAMDA3(1,1,2),LAMDA3(1,1,3),LAMDA3(1,2,3),LAMDA3(2,1,2),
     &        LAMDA3(2,1,3),LAMDA3(2,2,3),LAMDA3(3,1,2),LAMDA3(3,1,3),
     &        LAMDA3(3,2,3)
C--use the antisymmetry to find the rest
        DO I=1,3
          DO J=1,3
            DO K=1,3
              IF(J.EQ.K) LAMDA3(I,J,K) = 0.0E0
              IF(J.GT.K) LAMDA3(I,J,K) = -LAMDA3(I,K,J)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO I=1,3
          DO J=1,3
            DO K=1,3
              LAMDA3(I,J,K) = 0.
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      CALL RPDECY
      END
