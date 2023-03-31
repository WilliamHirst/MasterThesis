CDECK  ID>, FBRBM.
      REAL FUNCTION FBRBM(X)
C
C     Integrand for convolution of 
C     bremsstrahlung with beamstrahlung spectra
C
      IMPLICIT NONE
      COMMON/BREMBM/QSQBM,EB,XMIN
      REAL QSQBM,EB,XMIN
      SAVE /BREMBM/
C
      REAL EBEAM,ESTRUC,X
C
      FBRBM=EBEAM(X,EB)*ESTRUC(XMIN/X,QSQBM)/X
      RETURN
      END
