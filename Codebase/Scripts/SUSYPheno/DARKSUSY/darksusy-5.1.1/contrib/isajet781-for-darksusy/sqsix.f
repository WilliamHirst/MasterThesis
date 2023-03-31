!
      SUBROUTINE SQSIX(GIN,QEND)
!
!Purpose: Receives the weak scale parameters.
!         Computes both the up-type and down-type mass matrices,
!         and (if SVLQ=1) subsequently the flavour changing stop
!         decay rate.
!         Also carries out a calculation based upon the one-step
!         formula by Hikasa-Kobayashi for comparison with the full
!         RGE method.
!
!         Must be used in conjunction with isajet since it uses
!         isajet results for neutralino eigenvectors.
!
!         The couplings are run to the scale of the lightest decoupling
!         chosen from both left- and right-handed squarks suitable for
!         the basis chosen with SVLQ, before the mass matrices are
!         calculated.
!
      IMPLICIT NONE
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/DECCALC/T1EVE,T1EVA,USQM,COSTHT,SINTHT,GHIK,MST1,MST2,GAMMA
      DOUBLE COMPLEX T1EVE(6,6),T1EVA(6),USQM(6,6),COSTHT,SINTHT
     $              ,GHIK(601)
      DOUBLE PRECISION MST1,MST2,GAMMA
      SAVE/DECCALC/
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
      DOUBLE COMPLEX GIN(601),GDEC(601)
      DOUBLE PRECISION Q,QEND
      INTEGER I
!
      Q=QEND !QEND is the scale at which RGEFLAV called SQSIX
!
!Isajet was called with the opposite signs of mu and A_0
!since it also flips M_0. To make sure everything is
!consistent, I calculate the rest using the isajet
!numerical inputs, ie I swap the signs of M_{1,2,3}, A_0
!and mu.
!
      DO I=31,60
        GIN(I)=-GIN(I)
        GIN(I+290)=-GIN(I+290)
        IF(I.LT.58)GIN(I+369)=-GIN(I+369)
        IF(I.LT.58)GIN(I+399)=-GIN(I+399)
        IF(I.LT.34)GIN(I+568)=-GIN(I+568)
      END DO
      GIN(108)=-GIN(108)
      GIN(398)=-GIN(398)
!
      DO I=1,601
        GDEC(I)=GIN(I)
      END DO
!
!Save the couplings at m_H for use by the one-step calculation.
!GIN is passed to SQSIX in the current basis specified by SVLQ.
!
      IF(ABS(Q-QNH).LT.1.D-5)THEN
        DO I=1,601
          GHIK(I)=GIN(I)
        END DO
      END IF
!
!Input is g at m_H. Run to place where we will calculate the
!squark mass matrices, Q_DEC.
!
      CALL DECRUN(GIN,GDEC,Q,SVLQ)
!
!Calculate the (6x6) matrix for both the up-type and down-type squarks
!
      CALL USMMA(GDEC,Q,SVLQ)
      CALL DSMMA(GDEC,Q,SVLQ)
!
!Finally calculate the flavour changing decay rate, stop1-->c Neutralino
!
      IF(SVLQ.EQ.1)THEN
        CALL ST1CNEU(GDEC,Q)
      END IF
!
      RETURN
      END
