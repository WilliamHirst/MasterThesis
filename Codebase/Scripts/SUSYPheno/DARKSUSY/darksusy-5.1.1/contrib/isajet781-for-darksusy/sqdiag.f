!
      SUBROUTINE SQDIAG(GIN)
!
!Purpose: To find the matrices which will rotate to the sfermion
!         mass basis by diagonalising.
!
!         Contains a number of conditional statements to account
!         for the decoupling of squarks.
!
      IMPLICIT NONE
!
      COMMON/SQEIG/MQVE,MQVA,MUPVE,MUPVA,MDVE,MDVA,MLVE,MLVA,MEVE,MEVA
      DOUBLE COMPLEX MQVE(3,3),MUPVE(3,3),MDVE(3,3),MLVE(3,3),MEVE(3,3)
     $               ,MQVA(3),MUPVA(3),MDVA(3),MLVA(3),MEVA(3)
      SAVE/SQEIG/
!
      COMMON /SQROT/ RQTOT,RUPTOT,RDTOT,RLTOT,RETOT
     $               ,RQSAV,RUPSAV,RDSAV,RLSAV,RESAV
     $               ,OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      DOUBLE COMPLEX RQTOT(3,3),RUPTOT(3,3),RDTOT(3,3)
      DOUBLE COMPLEX RLTOT(3,3),RETOT(3,3)
      DOUBLE COMPLEX RQSAV(2,3,3),RUPSAV(2,3,3),RDSAV(2,3,3)
      DOUBLE COMPLEX RLSAV(2,3,3),RESAV(2,3,3)
      INTEGER OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      SAVE /SQROT/
!
      DOUBLE COMPLEX GIN(601),G(601)
      DOUBLE COMPLEX MQ(3,3),MUP(3,3),MD(3,3),ML(3,3),ME(3,3)
      DOUBLE COMPLEX MQ2(2,2),MQVA2(2),MQVE2(2,2)
      DOUBLE COMPLEX MUP2(2,2),MUPVA2(2),MUPVE2(2,2)
      DOUBLE COMPLEX MD2(2,2),MDVA2(2),MDVE2(2,2)
      DOUBLE COMPLEX ML2(2,2),MLVA2(2),MLVE2(2,2)
      DOUBLE COMPLEX ME2(2,2),MEVA2(2),MEVE2(2,2)
      DOUBLE COMPLEX EVERTMP(3,3),CWORK1(99),CWORK2(6)
      DOUBLE COMPLEX EVERTMP2(2,2),CWORK12(99),CWORK22(4)
      INTEGER I,J,DIAGERR,CIERR1,CIERR2,CIERR3,CIERR4,CIERR5
!
      CIERR1=0
      CIERR2=0
      CIERR3=0
      CIERR4=0
      CIERR5=0
      DIAGERR=0
      DO I=1,3
        MQVA(I)=(0.D0,0.D0)
        MUPVA(I)=(0.D0,0.D0)
        MDVA(I)=(0.D0,0.D0)
        MLVA(I)=(0.D0,0.D0)
        MEVA(I)=(0.D0,0.D0)
      END DO
!
!First rotate to the correct basis for the scale we are
!considering.
!
      CALL ROTSQ(GIN,G)
!
      DO I=1,3
        DO J=1,3
          MQ(I,J)=G(62+(I-1)*3+J)
          MUP(I,J)=G(80+(I-1)*3+J)
          MD(I,J)=G(89+(I-1)*3+J)
          ML(I,J)=G(71+(I-1)*3+J)
          ME(I,J)=G(98+(I-1)*3+J)
        END DO
      END DO
!
!Separate out the (2x2) entries for use if the heaviest sfermion
!has already decoupled.
!
      DO I=1,2
        DO J=1,2
          IF(OLDNSQ.EQ.2)THEN
            MQ2(I,J)=MQ(I,J)
          END IF
          IF(OLDNSU.EQ.2)THEN
            MUP2(I,J)=MUP(I,J)
          END IF
          IF(OLDNSD.EQ.2)THEN
            MD2(I,J)=MD(I,J)
          END IF
          IF(OLDNSL.EQ.2)THEN
            ML2(I,J)=ML(I,J)
          END IF
          IF(OLDNSE.EQ.2)THEN
            ME2(I,J)=ME(I,J)
          END IF
        END DO
      END DO
!
!Find eigenvectors and eigenvalues. If there are only two active
!squarks in each group, then only diagonalise the (2x2).
!
      IF(OLDNSQ.EQ.3)THEN
        CALL ZGEEV('V','V',3,MQ,3,MQVA,MQVE,3,EVERTMP,3,CWORK1,99
     $             ,CWORK2,CIERR1)
      ELSE IF(OLDNSQ.EQ.2)THEN
        CALL ZGEEV('V','V',2,MQ2,2,MQVA2,MQVE2,2,EVERTMP2,2,CWORK12,99
     $             ,CWORK22,CIERR1)
      END IF
      IF(OLDNSU.EQ.3)THEN
        CALL ZGEEV('V','V',3,MUP,3,MUPVA,MUPVE,3,EVERTMP,3,CWORK1,99
     $             ,CWORK2,CIERR2)
      ELSE IF(OLDNSU.EQ.2)THEN
        CALL ZGEEV('V','V',2,MUP2,2,MUPVA2,MUPVE2,2,EVERTMP2,2,CWORK12
     $             ,99,CWORK22,CIERR2)
      END IF
      IF(OLDNSD.EQ.3)THEN
        CALL ZGEEV('V','V',3,MD,3,MDVA,MDVE,3,EVERTMP,3,CWORK1,99
     $             ,CWORK2,CIERR3)
      ELSE IF(OLDNSD.EQ.2)THEN
        CALL ZGEEV('V','V',2,MD2,2,MDVA2,MDVE2,2,EVERTMP2,2,CWORK12,99
     $             ,CWORK22,CIERR3)
      END IF
      IF(OLDNSL.EQ.3)THEN
        CALL ZGEEV('V','V',3,ML,3,MLVA,MLVE,3,EVERTMP,3,CWORK1,99
     $             ,CWORK2,CIERR4)
      ELSE IF(OLDNSL.EQ.2)THEN
        CALL ZGEEV('V','V',2,ML2,2,MLVA2,MLVE2,2,EVERTMP2,2,CWORK12,99
     $             ,CWORK22,CIERR4)
      END IF
      IF(OLDNSE.EQ.3)THEN
        CALL ZGEEV('V','V',3,ME,3,MEVA,MEVE,3,EVERTMP,3,CWORK1,99
     $             ,CWORK2,CIERR5)
      ELSE IF(OLDNSE.EQ.2)THEN
        CALL ZGEEV('V','V',2,ME2,2,MEVA2,MEVE2,2,EVERTMP2,2,CWORK12,99
     $             ,CWORK22,CIERR5)
      END IF
!
!If only two active sfermions, set the third eigenvalue to zero
!and the (3,3) entry in the rotation to (1.d0,0.d0)
!
      DO I=1,3
        DO J=1,3
          IF(OLDNSQ.EQ.2)THEN
            IF(I.LT.3.AND.J.LT.3)THEN
              MQVE(I,J)=MQVE2(I,J)
            ELSE
              MQVE(I,J)=(0.D0,0.D0)
              IF(I.EQ.J)MQVE(3,3)=(1.D0,0.D0)
            END IF
          ELSE IF(OLDNSQ.LT.2)THEN
            MQVE(I,J)=(0.D0,0.D0)
            IF(I.EQ.J)MQVE(I,J)=(1.D0,0.D0)
          END IF
          IF(OLDNSU.EQ.2)THEN
            IF(I.LT.3.AND.J.LT.3)THEN
              MUPVE(I,J)=MUPVE2(I,J)
            ELSE
              MUPVE(I,J)=(0.D0,0.D0)
              IF(I.EQ.J)MUPVE(3,3)=(1.D0,0.D0)
            END IF
          ELSE IF(OLDNSU.LT.2)THEN
            MUPVE(I,J)=(0.D0,0.D0)
            IF(I.EQ.J)MUPVE(I,J)=(1.D0,0.D0)
          END IF
          IF(OLDNSD.EQ.2)THEN
            IF(I.LT.3.AND.J.LT.3)THEN
              MDVE(I,J)=MDVE2(I,J)
            ELSE
              MDVE(I,J)=(0.D0,0.D0)
              IF(I.EQ.J)MDVE(3,3)=(1.D0,0.D0)
            END IF
          ELSE IF(OLDNSD.LT.2)THEN
            MDVE(I,J)=(0.D0,0.D0)
            IF(I.EQ.J)MDVE(I,J)=(1.D0,0.D0)
          END IF
          IF(OLDNSL.EQ.2)THEN
            IF(I.LT.3.AND.J.LT.3)THEN
              MLVE(I,J)=MLVE2(I,J)
            ELSE
              MLVE(I,J)=(0.D0,0.D0)
              IF(I.EQ.J)MLVE(3,3)=(1.D0,0.D0)
            END IF
          ELSE IF(OLDNSL.LT.2)THEN
            MLVE(I,J)=(0.D0,0.D0)
            IF(I.EQ.J)MLVE(I,J)=(1.D0,0.D0)
          END IF
          IF(OLDNSE.EQ.2)THEN
            IF(I.LT.3.AND.J.LT.3)THEN
              MEVE(I,J)=MEVE2(I,J)
            ELSE
              MEVE(I,J)=(0.D0,0.D0)
              IF(I.EQ.J)MEVE(3,3)=(1.D0,0.D0)
            END IF
          ELSE IF(OLDNSE.LT.2)THEN
            MEVE(I,J)=(0.D0,0.D0)
            IF(I.EQ.J)MEVE(I,J)=(1.D0,0.D0)
          END IF
        END DO
        IF(OLDNSQ.EQ.2)THEN
          IF(I.LT.3)THEN
            MQVA(I)=MQVA2(I)
            MQVA(3)=(0.D0,0.D0)
          END IF
        END IF
        IF(OLDNSU.EQ.2)THEN
          IF(I.LT.3)THEN
           MUPVA(I)=MUPVA2(I)
           MUPVA(3)=(0.D0,0.D0)
          END IF
        END IF
        IF(OLDNSD.EQ.2)THEN
          IF(I.LT.3)THEN
            MDVA(I)=MDVA2(I)
            MDVA(3)=(0.D0,0.D0)
          END IF
        END IF
        IF(OLDNSL.EQ.2)THEN
          IF(I.LT.3)THEN
            MLVA(I)=MLVA2(I)
            MLVA(3)=(0.D0,0.D0)
          END IF
        END IF
        IF(OLDNSE.EQ.2)THEN
          IF(I.LT.3)THEN
            MEVA(I)=MEVA2(I)
            MEVA(3)=(0.D0,0.D0)
          END IF
        END IF
      END DO
      IF(OLDNSQ.EQ.1)THEN
        MQVA(1)=MQ(1,1)
        DO I=2,3
          MQVA(I)=(0.D0,0.D0)
        END DO
      END IF
      IF(OLDNSU.EQ.1)THEN
        MUPVA(1)=MUP(1,1)
        DO I=2,3
          MUPVA(I)=(0.D0,0.D0)
        END DO
      END IF
      IF(OLDNSD.EQ.1)THEN
        MDVA(1)=MD(1,1)
        DO I=2,3
          MDVA(I)=(0.D0,0.D0)
        END DO
      END IF
      IF(OLDNSL.EQ.1)THEN
        MLVA(1)=ML(1,1)
        DO I=2,3
          MLVA(I)=(0.D0,0.D0)
        END DO
      END IF
      IF(OLDNSE.EQ.1)THEN
        MEVA(1)=ME(1,1)
        DO I=2,3
          MEVA(I)=(0.D0,0.D0)
        END DO
      END IF
!
!Sort the eigenvectors into order and normalise them
!
      CALL SORTZG(MQVA,MQVE,OLDNSQ)
      CALL SORTZG(MUPVA,MUPVE,OLDNSU)
      CALL SORTZG(MDVA,MDVE,OLDNSD)
      CALL SORTZG(MLVA,MLVE,OLDNSL)
      CALL SORTZG(MEVA,MEVE,OLDNSE)
!
!Error checking
!
      DO I=1,3
        IF(MQVE(1,I).EQ.(0.D0,0.D0).AND.MQVE(2,I).EQ.(0.D0,0.D0)
     $    .AND.MQVE(3,I).EQ.(0.D0,0.D0))DIAGERR=1
        IF(MUPVE(1,I).EQ.(0.D0,0.D0).AND.MUPVE(2,I).EQ.(0.D0,0.D0)
     $    .AND.MUPVE(3,I).EQ.(0.D0,0.D0))DIAGERR=1
        IF(MDVE(1,I).EQ.(0.D0,0.D0).AND.MDVE(2,I).EQ.(0.D0,0.D0)
     $    .AND.MDVE(3,I).EQ.(0.D0,0.D0))DIAGERR=1
        IF(MLVE(1,I).EQ.(0.D0,0.D0).AND.MLVE(2,I).EQ.(0.D0,0.D0)
     $    .AND.MLVE(3,I).EQ.(0.D0,0.D0))DIAGERR=1
        IF(MEVE(1,I).EQ.(0.D0,0.D0).AND.MEVE(2,I).EQ.(0.D0,0.D0)
     $    .AND.MEVE(3,I).EQ.(0.D0,0.D0))DIAGERR=1
      END DO
      IF(DIAGERR.EQ.1)WRITE(*,*)
     $           'ERROR IN SQDIAG'
!
      IF(CIERR1.NE.0.OR.CIERR2.NE.0.OR.CIERR3.NE.0.OR.CIERR4.NE.0
     $  .OR.CIERR5.NE.0)WRITE(*,*)'ERROR WITH DIAGONALISATION IN SQDIAG'
!
      RETURN
      END
