CDECK  ID>, RPRATE
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
      FUNCTION RPRATE(LAMCOL,M1,M2,M3)
C-----------------------------------------------------------------------
C     FUNCTION TO CALCULATE A 2 BODY R-PARITY VIOLATING DECAY RATE
C-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL LAMCOL,M1,M2,M3,RPRATE,PCM,PI
      PI = 3.1415926E0
      RPRATE = 0
      IF(M1.LT.(M2+M3)) RETURN
      RPRATE = LAMCOL*(M1**2-M2**2-M3**2)
      PCM = SQRT((M1**2-(M2+M3)**2)*(M1**2-(M2-M3)**2))/(2*M1)
      RPRATE = RPRATE*PCM/(8*PI*M1**2)
      END
