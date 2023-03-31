 
      FUNCTION GAMTOT(ID)
C-----------------------------------------------------------------------
C          Calculate total width (in GeV) for ID
C-----------------------------------------------------------------------
      IMPLICIT NONE
      
C          MXSS         =  maximum number of modes
C          NSSMOD       = number of modes
C          ISSMOD       = initial particle
C          JSSMOD       = final particles
C          GSSMOD       = width
C          BSSMOD       = branching ratio
C          MSSMOD       = decay matrix element pointer
C          LSSMOD       = logical flag used internally by SSME3
      INTEGER MXSS
      PARAMETER (MXSS=1000)
      COMMON/SSMODE/NSSMOD,ISSMOD(MXSS),JSSMOD(5,MXSS),GSSMOD(MXSS)
     $,BSSMOD(MXSS),MSSMOD(MXSS),LSSMOD
      INTEGER NSSMOD,ISSMOD,JSSMOD,MSSMOD
      REAL GSSMOD,BSSMOD
      LOGICAL LSSMOD
      SAVE /SSMODE/
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
C
      INTEGER ID,I
      REAL GAMTOT
C
      GAMTOT=0
      DO 100 I=1,NSSMOD
        IF(ISSMOD(I).EQ.ID) GAMTOT=GAMTOT+GSSMOD(I)
100   CONTINUE
      RETURN
      END 
