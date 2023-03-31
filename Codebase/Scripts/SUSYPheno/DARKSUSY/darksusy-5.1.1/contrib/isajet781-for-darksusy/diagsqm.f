!
      SUBROUTINE DIAGSQM(QMSQ,UMSQ,DMSQ,LMSQ,EMSQ)
!
!Purpose: To use the results of MASSSQM for the squark mass
!         squared matrices in the quark mass basis to find the
!         eigenvalues correspoding to u,c,t ; d,s,b ; e,nu,tau
!         in that order. Used by EWSB routine.
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
      DOUBLE COMPLEX QVA(3),UVA(3),DVA(3),LVA(3),EVA(3)
      DOUBLE COMPLEX QVE(3,3),UVE(3,3),DVE(3,3),LVE(3,3),EVE(3,3)
      DOUBLE COMPLEX QIN(3,3),UIN(3,3),DIN(3,3),EIN(3,3),LIN(3,3)
      DOUBLE COMPLEX EVERTMP(3,3),CWORK1(99),CWORK2(6)
      DOUBLE PRECISION VALUEQ,VALUEU,VALUED,VALUEL,VALUEE
      DOUBLE PRECISION QMSQ(3),UMSQ(3),DMSQ(3),LMSQ(3),EMSQ(3)
      INTEGER I,J
      INTEGER ENTRYQ,ENTRYU,ENTRYD,ENTRYL,ENTRYE
      INTEGER CIERR1,CIERR2,CIERR3,CIERR4,CIERR5
!
      DO I=1,3
        DO J=1,3
          QIN(I,J)=MQQMASS(I,J)
          UIN(I,J)=MUQMASS(I,J)
          DIN(I,J)=MDQMASS(I,J)
          LIN(I,J)=MLQMASS(I,J)
          EIN(I,J)=MEQMASS(I,J)
        END DO
      END DO
!
!Diagonalise and recombine into complex quantities
!
      CALL ZGEEV('V','V',3,QIN,3,QVA,QVE,3,EVERTMP,3,CWORK1,99
     $           ,CWORK2,CIERR1)
      CALL ZGEEV('V','V',3,UIN,3,UVA,UVE,3,EVERTMP,3,CWORK1,99
     $           ,CWORK2,CIERR2)
      CALL ZGEEV('V','V',3,DIN,3,DVA,DVE,3,EVERTMP,3,CWORK1,99
     $           ,CWORK2,CIERR3)
      CALL ZGEEV('V','V',3,LIN,3,LVA,LVE,3,EVERTMP,3,CWORK1,99
     $           ,CWORK2,CIERR4)
      CALL ZGEEV('V','V',3,EIN,3,EVA,EVE,3,EVERTMP,3,CWORK1,99
     $           ,CWORK2,CIERR5)
!
!Now work out the correct order of the eigenvalues
!depending on the eigenvector.
!
      DO J=1,3
        ENTRYQ=1
        VALUEQ=ABS(QVE(1,J))
        ENTRYU=1
        VALUEU=ABS(UVE(1,J))
        ENTRYD=1
        VALUED=ABS(DVE(1,J))
        ENTRYL=1
        VALUEL=ABS(LVE(1,J))
        ENTRYE=1
        VALUEE=ABS(EVE(1,J))
        DO I=2,3
          IF(ABS(QVE(I,J)).GT.VALUEQ)THEN
            VALUEQ=ABS(QVE(I,J))
            ENTRYQ=I
          END IF
          IF(ABS(UVE(I,J)).GT.VALUEU)THEN
            VALUEU=ABS(UVE(I,J))
            ENTRYU=I
          END IF
          IF(ABS(DVE(I,J)).GT.VALUED)THEN
            VALUED=ABS(DVE(I,J))
            ENTRYD=I
          END IF
          IF(ABS(LVE(I,J)).GT.VALUEL)THEN
            VALUEL=ABS(LVE(I,J))
            ENTRYL=I
          END IF
          IF(ABS(EVE(I,J)).GT.VALUEE)THEN
            VALUEE=ABS(EVE(I,J))
            ENTRYE=I
          END IF
        END DO
        QMSQ(ENTRYQ)=DSQRT(ABS(QVA(J))**2)
        UMSQ(ENTRYU)=DSQRT(ABS(UVA(J))**2)
        DMSQ(ENTRYD)=DSQRT(ABS(DVA(J))**2)
        LMSQ(ENTRYL)=DSQRT(ABS(LVA(J))**2)
        EMSQ(ENTRYE)=DSQRT(ABS(EVA(J))**2)
      END DO
!
      RETURN
      END
