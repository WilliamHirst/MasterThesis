C **********************************************************************
C
      SUBROUTINE MOMNTX(ENERGY,MASS,COSTH,PHI , P)
C
C This subroutine sets up a four-momentum from the four inputs.         
C                                                                       
C INPUT:                                                                
C       real    ENERGY         : energy                                 
C       real    MASS           : mass                                   
C       real    COSTH          : cos(theta)                             
C       real    PHI            : azimuthal angle                        
C                                                                       
C OUTPUT:                                                               
C       real    P(0:3)         : four-momentum                          
C
      IMPLICIT NONE
      REAL*8    P(0:3),ENERGY,MASS,COSTH,PHI,PP,SINTH
C
      P(0) = ENERGY
      IF (ENERGY.EQ.MASS) THEN
         P(1) = 0.
         P(2) = 0.
         P(3) = 0.
      ELSE
         PP=SQRT((ENERGY-MASS)*(ENERGY+MASS))
         SINTH=SQRT((1.-COSTH)*(1.+COSTH))
         P(3) = PP*COSTH
         IF (PHI.EQ.0.) THEN
            P(1) = PP*SINTH
            P(2) = 0.
         ELSE
            P(1) = PP*SINTH*COS(PHI)
            P(2) = PP*SINTH*SIN(PHI)
         ENDIF
      ENDIF
      RETURN
      END
