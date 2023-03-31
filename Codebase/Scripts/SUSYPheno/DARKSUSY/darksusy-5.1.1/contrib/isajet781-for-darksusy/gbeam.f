CDECK  ID>, GBEAM.
      REAL FUNCTION GBEAM(X,E)
C***********************************************************************
C* Computes the effective single photon spectrum from beamstrahlung at *
C* e+e- colliders, using Peskin's approximate expression, 
C* Eq. 17 of SLAC-TN-04-032.
C* For a given
C* beamstrahlungs parameter Y; is supposed to work for Y <= 10 or so.  *
C***********************************************************************
      IMPLICIT NONE
      COMMON/EEPAR/SGMXEE,PLEP,PLEM,RSHMIN,RSHMAX,
     $UPSLON,SIGZ,IBREM,IBEAM,GAMGAM
      SAVE /EEPAR/
      REAL SGMXEE,PLEP,PLEM,RSHMIN,RSHMAX,UPSLON,SIGZ
      LOGICAL IBREM,IBEAM,GAMGAM
C
      REAL X,E,GAMMA
      DOUBLE PRECISION NGAM,NUCL,NUGAM,Y,GAM,GAM13,
     ,RE,XKAPPA,Z,E2,DE,DX,DGBEAM,XL
      IF (X.LT.1.E-5) X=1.E-5
      DE=E
      DX=X
      Y=UPSLON
      GAM = DE/5.11D-4
      RE = 1.D0/(137.D0*5.11D-4)
      XL=SIGZ*1.D12/.197327D0
      XKAPPA = 2.D0/(3.D0*Y)
      NUCL = 2.5D0*Y/(SQRT(3.D0)*137.D0**2*GAM*RE)
      NUGAM = NUCL/SQRT(1.D0+Y**.6666666D0)
      NGAM = DSQRT(3.D0)*XL*NUGAM
      Z = 1.D0-DX
      E2 = DX*XKAPPA/Z
      GAM13=DBLE(GAMMA(.333333))
      DGBEAM = (NGAM/2.D0)*XKAPPA**.333333D0/DX**.666666D0/
     ,Z**1.333333D0*DEXP(-E2)/GAM13
      GBEAM=DGBEAM
      IF (GBEAM.LT.0.) GBEAM=0.
      RETURN
      END
