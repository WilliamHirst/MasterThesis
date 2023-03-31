CDECK  ID>, EBEAM.
      REAL FUNCTION EBEAM(X,E)
C***********************************************************************
C* Computes the effective single electron spectrum from beamstrahlung at *
C* e+e- colliders, using Peskin expression Eq. 6 of SLAC-TN-04-032.    *
C* Use beamstrahlung  parameter Y; is supposed to work for Y <= 10     *
C* The quantities in the COMMON block are the beamstrahlungs parameter *
C* Y, the bunch length XL in GeV, the number of photons NGAM, and the  *
C* parameters NUCL, NUGAM, W, XKAPPA defined by Chen, as well as the   *
C* pre-factor FAC. Y, E and XLMM are read in by BEAM when it is called *
C* for the first time, with INIT=1; in this first run the other para-  *
C* meters are then computed, and simply used in later calls with       *
C* INIT = 0. This COMMON block should be present in the main program   *
C* in order to guarantee the survival of these parameters. Finally, X  *
C* is the electron energy in units of the nominal beam energy. Notice  *
C* that EBEAM  is only the part which is NOT proportional to           *
C* delta(1-X); the coefficient of the delta-function is simply         *
C* exp(-N).                                                            *
C***********************************************************************
      IMPLICIT NONE
      COMMON/EEPAR/SGMXEE,PLEP,PLEM,RSHMIN,RSHMAX,
     $UPSLON,SIGZ,IBREM,IBEAM,GAMGAM
      SAVE /EEPAR/
      REAL SGMXEE,PLEP,PLEM,RSHMIN,RSHMAX,UPSLON,SIGZ
      LOGICAL IBREM,IBEAM,GAMGAM
C
      REAL X,E,GAMMA
      REAL*8 DX,DE,Z,Y,XL,GAM,RE,XKAPPA,NUCL,NUGAM,NGAM,NGAM2,
     $YY,HTOT,H,HTOT0,DEBEAM,DNFAC
      INTEGER N,NFAC,II
      IF (X.LT.1.E-5) X=1.E-5
      IF (X.GT..99999) X=.99999
      DX=X
      DE=E
      Z=1.D0-DX
      Y=UPSLON
      XL = SIGZ*1.D12/.197327D0
      GAM = DE/5.11D-4
      RE = 1.D0/(137.D0*5.11D-4)
      XKAPPA = 2.D0/(3.D0*Y)
      NUCL = 2.5D0*Y/(SQRT(3.D0)*137.D0**2*GAM*RE)
      NUGAM = NUCL/SQRT(1.D0+Y**.6666666D0)
      NGAM = DSQRT(3.D0)*NUGAM*XL
      NGAM2=NGAM/2.D0
      YY=NGAM2*(XKAPPA*Z/X)**.333333D0
      HTOT=0.D0
      DO N=1,20
        NFAC=1
        DO II=1,N
          NFAC=NFAC*II
        END DO
        DNFAC=DBLE(FLOAT(NFAC))
        H=YY**N/DNFAC/DBLE(GAMMA(FLOAT(N)/3.))
        HTOT0=HTOT
        HTOT=HTOT+H
        IF (ABS(HTOT0-HTOT)/HTOT.LT..001D0.AND.N.GT.3) GO TO 101
      END DO
101   CONTINUE
      DEBEAM=DEXP(-XKAPPA*Z/DX)/Z/DX*HTOT*DEXP(-NGAM2)
      EBEAM=DEBEAM
      IF(EBEAM.LT.0.) EBEAM = 1.E-20
      RETURN
      END
