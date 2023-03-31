CDECK  ID>, RPMODA
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
      SUBROUTINE RPMODA(RATE,PIN,POUT1,POUT2,POUT3)
C-----------------------------------------------------------------------
C     SUBROUTINE TO ADDED MODE TO DECAY TABLE
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NRPMD
      PARAMETER (NRPMD=5000)
      COMMON/RPARRT/NSSMD2,ISSMD2(NRPMD),JSSMD2(5,NRPMD),
     &              GSSMD2(NRPMD),BSSMD2(NRPMD)
      INTEGER NSSMD2,ISSMD2,JSSMD2,PIN,POUT1,POUT2,POUT3
      REAL    GSSMD2,BSSMD2,RATE
      SAVE /RPARRT/
      IF(RATE.GT.1E-40.AND.NSSMD2+1.LE.NRPMD) THEN
        NSSMD2=NSSMD2+1
        ISSMD2(NSSMD2)   = PIN
        JSSMD2(1,NSSMD2) = POUT1
        JSSMD2(2,NSSMD2) = POUT2
        JSSMD2(3,NSSMD2) = POUT3 
        JSSMD2(4,NSSMD2) = 0 
        JSSMD2(5,NSSMD2) = 0
        GSSMD2(NSSMD2) = RATE
        BSSMD2(NSSMD2) = 0
      ELSEIF(NSSMD2+1.GT.NRPMD) THEN
        print *,'TOO MANY R-PARITY VIOLATING MODES'
        print *,'INCREASE NRPMD AND RERUN'
        STOP
      ENDIF
      END
