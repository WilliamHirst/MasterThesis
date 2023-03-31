!
      SUBROUTINE UPSQM(G,Q,OUT)
!
!Purpose: Construct the 6x6 mass matrix for up-type squarks.
!         The basis is u_l,c_l,t_l,u_r,c_r,t_r.
!
      IMPLICIT NONE
!
      COMMON/MYDECAY/MQQMASS,MUQMASS,MDQMASS,MLQMASS,MEQMASS,
     $             OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $             OFFMAXEVAL,OFFMAXQ,OFFMAXU,OFFMAXD,OFFMAXL,OFFMAXE
      DOUBLE COMPLEX MQQMASS(3,3),MUQMASS(3,3),MDQMASS(3,3),
     $               MLQMASS(3,3),MEQMASS(3,3)
      DOUBLE COMPLEX OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $               OFFMAXEVAL
      INTEGER OFFMAXQ(2),OFFMAXU(2),OFFMAXD(2),OFFMAXL(2),OFFMAXE(2)
      SAVE/MYDECAY/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      DOUBLE COMPLEX G(601),OUT(6,6)
!
      DOUBLE COMPLEX AU(3,3),MQ(3,3),MU(3,3),MTSFU(3,3)
      DOUBLE COMPLEX FUQ(3,3),FUUR(3,3)
      DOUBLE COMPLEX GLP,G2L
      DOUBLE COMPLEX AUT(3,3),MTSFUT(3,3),FUQT(3,3),FUURT(3,3)
      DOUBLE COMPLEX AUS(3,3),MTSFUS(3,3),FUQS(3,3),FUURS(3,3)
      DOUBLE COMPLEX TRIU(3,3),SFUQ(3,3),SFUUR(3,3)
      DOUBLE COMPLEX TRIUT(3,3),SFUQT(3,3),SFUURT(3,3)
      DOUBLE COMPLEX TRIUS(3,3),SFUQS(3,3),SFUURS(3,3)
      DOUBLE COMPLEX ID(3,3),CMATMUL
      DOUBLE PRECISION TANB,SINB,COSB,Q,VU,VD,VSM,VUSQ,VDSQ,VSMSQ
      INTEGER I,J
!
      DATA ID(1,1)/(1.D0,0.D0)/,ID(1,2)/(0.D0,0.D0)/
     $     ,ID(1,3)/(0.D0,0.D0)/
      DATA ID(2,1)/(0.D0,0.D0)/,ID(2,2)/(1.D0,0.D0)/
     $     ,ID(2,3)/(0.D0,0.D0)/
      DATA ID(3,1)/(0.D0,0.D0)/,ID(3,2)/(0.D0,0.D0)/
     $     ,ID(3,3)/(1.D0,0.D0)/
!
      DO I=1,6
        DO J=1,6
          OUT(I,J)=(0.D0,0.D0)
        END DO
      END DO
!
!Insert the soft mass matrices in the quark mass basis
!
      DO I=1,3
        DO J=1,3
          MQ(I,J)=MQQMASS(I,J)
          MU(I,J)=MUQMASS(I,J)
        END DO
      END DO
!
!Calculate cos (beta) and sin (beta)
!
      VU=DBLE(G(110))
      VD=DBLE(G(111))
      TANB=VU/VD
      COSB=DSQRT(1.D0/(1.D0+TANB**2))
      SINB=TANB*COSB
      VUSQ=VU**2
      VDSQ=VD**2
!
      VSM=DBLE(G(428))
      VSMSQ=VSM**2
!
!Convert the other entries in G(601) to the matrices needed
!
      DO I=1,3
        DO J=1,3
          AU(I,J)=G(33+(I-1)*3+J)
          MTSFU(I,J)=G(429+(I-1)*3+J)
          FUQ(I,J)=G(3+(I-1)*3+J) !Quartics are not run independently
          FUUR(I,J)=G(3+(I-1)*3+J)
!
          TRIU(I,J)=G(399+(I-1)*3+J)
          SFUQ(I,J)=G(111+(I-1)*3+J)
          SFUUR(I,J)=G(111+(I-1)*3+J)
        END DO
      END DO
!
!Quartics are not run independently at this time.
!Set other quartics to be their non-tilde counterparts.
!This can be fixed if the quartic running is introduced.
!*****NB: IF THE QUARTIC RUNNING IS INTRODUCED CARE MUST
!         BE TAKEN WITH G(541) AND G(542) SINCE THEY ARE
!         DEFINED DIFFERENTLY IN THE COMPLEX VERSION AS
!         OPPOSED TO THE REAL VERSION. SEE THE NOTE AT
!         THE BEGINNING OF drge601.f
!
      GLP=DSQRT(3.D0/5.D0)*G(1)
      G2L=G(2)
      G(541)=SQRT(DCMPLX(COSB**2)-DCMPLX(SINB**2))*DSQRT(3.D0/5.D0)*G(1)
      G(542)=SQRT(DCMPLX(COSB**2)-DCMPLX(SINB**2))*G(2)
!
      DO I=1,3
        DO J=1,3
          AUT(I,J)=AU(J,I)
          MTSFUT(I,J)=MTSFU(J,I)
          FUQT(I,J)=FUQ(J,I)
          FUURT(I,J)=FUUR(J,I)
          AUS(I,J)=CONJG(AU(I,J))
          MTSFUS(I,J)=CONJG(MTSFU(I,J))
          FUQS(I,J)=CONJG(FUQ(I,J))
          FUURS(I,J)=CONJG(FUUR(I,J))
!
          TRIUT(I,J)=TRIU(J,I)
          SFUQT(I,J)=SFUQ(J,I)
          SFUURT(I,J)=SFUUR(J,I)
          TRIUS(I,J)=CONJG(TRIU(I,J))
          SFUQS(I,J)=CONJG(SFUQ(I,J))
          SFUURS(I,J)=CONJG(SFUUR(I,J))
        END DO
      END DO
!
!Split matrix into 3x3 blocks.
!First, the top left block
!
      DO I=1,3
        DO J=1,3
          IF(Q.GT.QNH)THEN
            OUT(I,J)=VUSQ*CMATMUL(0,FUURS,FUURT,I,J)+MQ(I,J)
     $               +(GLP**2/12.D0-G2L**2/4.D0)*(VUSQ-VDSQ)*ID(I,J)
          ELSE
            OUT(I,J)=VSMSQ*CMATMUL(0,SFUURS,SFUURT,I,J)+MQ(I,J)
     $               -VSMSQ*(G(541)**2/12.D0-G(542)**2/4.D0)*ID(I,J)
          END IF
        END DO
      END DO
!
!Next the bottom left
!
      DO I=1,3
        DO J=1,3
          IF(Q.GT.QNH)THEN
            OUT(I+3,J)=VD*MTSFUT(I,J)-VU*AUT(I,J)
          ELSE
            OUT(I+3,J)=-(VSM*TRIUT(I,J))
          END IF
        END DO
      END DO
!
!Top right is the dagger of bottom left
!
      DO I=1,3
        DO J=1,3
          IF(Q.GT.QNH)THEN
            OUT(I,J+3)=VD*MTSFUS(I,J)-VU*AUS(I,J)
          ELSE
            OUT(I,J+3)=-(VSM*TRIUS(I,J))
          END IF
        END DO
      END DO
!
!Finally bottom right
!
      DO I=1,3
        DO J=1,3
          IF(Q.GT.QNH)THEN
            OUT(I+3,J+3)=VUSQ*CMATMUL(0,FUQT,FUQS,I,J)+MU(I,J)
     $                   -GLP**2/3.D0*(VUSQ-VDSQ)*ID(I,J)
          ELSE
            OUT(I+3,J+3)=VSMSQ*CMATMUL(0,SFUQT,SFUQS,I,J)+MU(I,J)
     $                   +VSMSQ*G(541)**2/3.D0*ID(I,J)
          END IF
        END DO
      END DO
!
      RETURN
      END
