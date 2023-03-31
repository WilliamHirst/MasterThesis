CDECK  ID>, RPRTCH
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
       FUNCTION RPRTCH(MASS,ENERGY)
C-----------------------------------------------------------------------
C      AVOID ERRORS DUE TO TAKING SQRT OF SMALL NEGATIVE NUMBERS
C      OCCURS DUE ROUNDING ERRORS
C-----------------------------------------------------------------------
       DOUBLE PRECISION RPRTCH,ENERGY,MASS
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
       RPRTCH =ENERGY**2-MASS**2
       IF(RPRTCH.LT.0) THEN
         IF(RPRTCH/(ENERGY**2+MASS**2).LT.1D-20) THEN
           RPRTCH = 0.0D0
         ELSE
           WRITE(LOUT,*) 'WARNING SQRT OF NEGATIVE NUMBER',RPRTCH
           RPRTCH = 0.0D0
         ENDIF
       ENDIF
       RPRTCH = SQRT(RPRTCH)
       END
