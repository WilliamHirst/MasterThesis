C
C ----------------------------------------------------------------------
C
      SUBROUTINE COUP2X(SW2 , GAL,GAU,GAD,GWF,GZN,GZL,GZU,GZD,G1)
C
C this subroutine sets up the coupling constants for the fermion-       
C fermion-vector vertices in the standard model.  the array of the      
C couplings specifies the chirality of the flowing-in fermion.  g??(1)  
C denotes a left-handed coupling, and g??(2) a right-handed coupling.   
C                                                                       
C input:                                                                
C       real    sw2            : square of sine of the weak angle       
C                                                                       
C output:                                                               
C       real    gal(2)         : coupling with a of charged leptons     
C       real    gau(2)         : coupling with a of up-type quarks      
C       real    gad(2)         : coupling with a of down-type quarks    
C       real    gwf(2)         : coupling with w-,w+ of fermions        
C       real    gzn(2)         : coupling with z of neutrinos           
C       real    gzl(2)         : coupling with z of charged leptons     
C       real    gzu(2)         : coupling with z of up-type quarks      
C       real    gzd(2)         : coupling with z of down-type quarks    
C       real    g1(2)          : unit coupling of fermions              
C
      IMPLICIT NONE
      REAL*8 GAL(2),GAU(2),GAD(2),GWF(2),GZN(2),GZL(2),GZU(2),GZD(2),
     &     G1(2),SW2,ALPHA,FOURPI,EE,SW,CW,EZ,EY
C
      REAL*8 RXZERO, RXHALF, RXONE, RXTWO, RTHREE, RXFOUR, RXOTE
      REAL*8 RXPI, RIALPH
      PARAMETER( RXZERO=0.0D0, RXHALF=0.5D0, RXONE=1.0D0, RXTWO=2.0D0,
     $     RTHREE=3.0D0 )
      PARAMETER( RXFOUR=4.0D0, RXOTE=128.0D0 )
      PARAMETER( RXPI=3.14159265358979323846D0, RIALPH=137.0359895D0 )
C
      ALPHA = RXONE / RXOTE
C      alpha = r_one / r_ialph
      FOURPI = RXFOUR * RXPI
      EE=SQRT( ALPHA * FOURPI )
      SW=SQRT( SW2 )
      CW=SQRT( RXONE - SW2 )
      EZ=EE/(SW*CW)
      EY=EE*(SW/CW)
C
      GAL(1) =  EE
      GAL(2) =  EE
      GAU(1) = -EE*RXTWO/RTHREE
      GAU(2) = -EE*RXTWO/RTHREE
      GAD(1) =  EE   /RTHREE
      GAD(2) =  EE   /RTHREE
      GWF(1) = -EE/SQRT(RXTWO*SW2)
      GWF(2) =  RXZERO
      GZN(1) = -EZ*  RXHALF
      GZN(2) =  RXZERO
      GZL(1) = -EZ*(-RXHALF+SW2)
      GZL(2) = -EY
      GZU(1) = -EZ*( RXHALF-SW2*RXTWO/RTHREE)
      GZU(2) =  EY*          RXTWO/RTHREE
      GZD(1) = -EZ*(-RXHALF+SW2   /RTHREE)
      GZD(2) = -EY             /RTHREE
      G1(1)  =  RXONE
      G1(2)  =  RXONE
C
      RETURN
      END
