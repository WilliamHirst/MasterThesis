!
      SUBROUTINE DOWNSQM(G,Q,OUT)
!
!Purpose: Construct the 6x6 mass matrix for down-type squarks.
!         The basis is d_l,s_l,b_l,d_r,s_r,b_r.
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
      DOUBLE COMPLEX AD(3,3),MQ(3,3),MD(3,3),MTSFD(3,3)
      DOUBLE COMPLEX FDQ(3,3),FDDR(3,3)
      DOUBLE COMPLEX GLP,G2L
      DOUBLE COMPLEX ADT(3,3),MTSFDT(3,3),FDQT(3,3),FDDRT(3,3)
      DOUBLE COMPLEX ADS(3,3),MTSFDS(3,3),FDQS(3,3),FDDRS(3,3)
      DOUBLE COMPLEX TRID(3,3),CFDQ(3,3),CFDDR(3,3)
      DOUBLE COMPLEX TRIDT(3,3),CFDQT(3,3),CFDDRT(3,3)
      DOUBLE COMPLEX TRIDS(3,3),CFDQS(3,3),CFDDRS(3,3)
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
          MD(I,J)=MDQMASS(I,J)
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
          AD(I,J)=G(42+(I-1)*3+J)
          MTSFD(I,J)=G(438+(I-1)*3+J)
          FDQ(I,J)=G(12+(I-1)*3+J) !Quartics are not run independently
          FDDR(I,J)=G(12+(I-1)*3+J)
!
          TRID(I,J)=G(408+(I-1)*3+J)
          CFDQ(I,J)=G(120+(I-1)*3+J)
          CFDDR(I,J)=G(120+(I-1)*3+J)
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
          ADT(I,J)=AD(J,I)
          MTSFDT(I,J)=MTSFD(J,I)
          FDQT(I,J)=FDQ(J,I)
          FDDRT(I,J)=FDDR(J,I)
          ADS(I,J)=CONJG(AD(I,J))
          MTSFDS(I,J)=CONJG(MTSFD(I,J))
          FDQS(I,J)=CONJG(FDQ(I,J))
          FDDRS(I,J)=CONJG(FDDR(I,J))
!
          TRIDT(I,J)=TRID(J,I)
          CFDQT(I,J)=CFDQ(J,I)
          CFDDRT(I,J)=CFDDR(J,I)
          TRIDS(I,J)=CONJG(TRID(I,J))
          CFDQS(I,J)=CONJG(CFDQ(I,J))
          CFDDRS(I,J)=CONJG(CFDDR(I,J))
        END DO
      END DO
!
!Split matrix into 3x3 blocks.
!First, the top left block
!
      DO I=1,3
        DO J=1,3
          IF(Q.GT.QNH)THEN
            OUT(I,J)=VDSQ*CMATMUL(0,FDDRS,FDDRT,I,J)+MQ(I,J)
     $               +(GLP**2/12.D0+G2L**2/4.D0)*(VUSQ-VDSQ)*ID(I,J)
          ELSE
            OUT(I,J)=VSMSQ*CMATMUL(0,CFDDRS,CFDDRT,I,J)+MQ(I,J)
     $               -VSMSQ*(G(541)**2/12.D0+G(542)**2/4.D0)*ID(I,J)
          END IF
        END DO
      END DO
!
!Next the bottom left
!
      DO I=1,3
        DO J=1,3
          IF(Q.GT.QNH)THEN
            OUT(I+3,J)=VU*MTSFDT(I,J)-VD*ADT(I,J)
          ELSE
            OUT(I+3,J)=-(VSM*TRIDT(I,J))
          END IF
        END DO
      END DO
!
!Top right is the dagger of bottom left
!
      DO I=1,3
        DO J=1,3
          IF(Q.GT.QNH)THEN
            OUT(I,J+3)=VU*MTSFDS(I,J)-VD*ADS(I,J)
          ELSE
            OUT(I,J+3)=-(VSM*TRIDS(I,J))
          END IF
        END DO
      END DO
!
!Finally bottom right
!
      DO I=1,3
        DO J=1,3
          IF(Q.GT.QNH)THEN
            OUT(I+3,J+3)=VDSQ*CMATMUL(0,FDQT,FDQS,I,J)+MD(I,J)
     $                   +GLP**2/6.D0*(VUSQ-VDSQ)*ID(I,J)
          ELSE
            OUT(I+3,J+3)=VSMSQ*CMATMUL(0,CFDQT,CFDQS,I,J)+MD(I,J)
     $                   -VSMSQ*G(541)**2/6.D0*ID(I,J)
          END IF
        END DO
      END DO
!
      RETURN
      END
