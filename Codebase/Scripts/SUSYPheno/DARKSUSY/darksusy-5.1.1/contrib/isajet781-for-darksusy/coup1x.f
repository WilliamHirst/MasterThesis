C
C **********************************************************************
C
      SUBROUTINE COUP1X(SW2 , GW,GWWA,GWWZ)
C
C this subroutine sets up the coupling constants of the gauge bosons in 
C the standard model.                                                   
C                                                                       
C input:                                                                
C       real    sw2            : square of sine of the weak angle       
C                                                                       
C output:                                                               
C       real    gw             : weak coupling constant                 
C       real    gwwa           : dimensionless coupling of w-,w+,a      
C       real    gwwz           : dimensionless coupling of w-,w+,z      
C
      IMPLICIT NONE
      REAL*8    SW2,GW,GWWA,GWWZ,ALPHA,FOURPI,EE,SW,CW
      REAL*8 RXONE, RXFOUR, RXOTE, RXPI, RIALPH
      PARAMETER( RXONE=1.0D0, RXFOUR=4.0D0, RXOTE=128.0D0 )
      PARAMETER( RXPI=3.14159265358979323846D0, RIALPH=137.0359895D0 )
C
      ALPHA = RXONE / RXOTE
C      alpha = r_one / r_ialph
      FOURPI = RXFOUR * RXPI
      EE=SQRT( ALPHA * FOURPI )
      SW=SQRT( SW2 )
      CW=SQRT( RXONE - SW2 )
C
      GW    =  EE/SW
      GWWA  =  EE
      GWWZ  =  EE*CW/SW
C
      RETURN
      END
