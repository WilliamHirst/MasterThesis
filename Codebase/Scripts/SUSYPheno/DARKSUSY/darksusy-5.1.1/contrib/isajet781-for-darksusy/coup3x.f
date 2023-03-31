C
C ----------------------------------------------------------------------
C
      SUBROUTINE COUP3X(SW2,ZMASS,HMASS , 
     &                  GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH)
C
C this subroutine sets up the coupling constants of the gauge bosons and
C higgs boson in the standard model.                                    
C                                                                       
C input:                                                                
C       real    sw2            : square of sine of the weak angle       
C       real    zmass          : mass of z                              
C       real    hmass          : mass of higgs                          
C                                                                       
C output:                                                               
C       real    gwwh           : dimensionful  coupling of w-,w+,h      
C       real    gzzh           : dimensionful  coupling of z, z, h      
C       real    ghhh           : dimensionful  coupling of h, h, h      
C       real    gwwhh          : dimensionful  coupling of w-,w+,h, h   
C       real    gzzhh          : dimensionful  coupling of z, z, h, h   
C       real    ghhhh          : dimensionless coupling of h, h, h, h   
C
      IMPLICIT NONE
      REAL*8    SW2,ZMASS,HMASS,GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH,
     &        ALPHA,FOURPI,EE2,SC2,V
C
      REAL*8 RXHALF, RXONE, RXTWO, RTHREE, RXFOUR, RXOTE
      REAL*8 RXPI, RIALPH
      PARAMETER( RXHALF=0.5D0, RXONE=1.0D0, RXTWO=2.0D0, RTHREE=3.0D0 )
      PARAMETER( RXFOUR=4.0D0, RXOTE=128.0D0 )
      PARAMETER( RXPI=3.14159265358979323846D0, RIALPH=137.0359895D0 )
C
      ALPHA = RXONE / RXOTE
C      alpha = r_one / r_ialph
      FOURPI = RXFOUR * RXPI
      EE2=ALPHA*FOURPI
      SC2=SW2*( RXONE - SW2 )
      V = RXTWO * ZMASS*SQRT(SC2)/SQRT(EE2)
C
      GWWH  =   EE2/SW2*RXHALF*V
      GZZH  =   EE2/SC2*RXHALF*V
      GHHH  =  -HMASS**2/V*RTHREE
      GWWHH =   EE2/SW2*RXHALF
      GZZHH =   EE2/SC2*RXHALF
      GHHHH = -(HMASS/V)**2*RTHREE
C
      RETURN
      END
