C
C ----------------------------------------------------------------------
C
      SUBROUTINE MOM2CX(ESUM,MASS1,MASS2,COSTH1,PHI1 , P1,P2)
C
C This subroutine sets up two four-momenta in the two particle rest     
C frame.                                                                
C                                                                       
C INPUT:                                                                
C       real    ESUM           : energy sum of particle 1 and 2         
C       real    MASS1          : mass            of particle 1          
C       real    MASS2          : mass            of particle 2          
C       real    COSTH1         : cos(theta)      of particle 1          
C       real    PHI1           : azimuthal angle of particle 1          
C                                                                       
C OUTPUT:                                                               
C       real    P1(0:3)        : four-momentum of particle 1            
C       real    P2(0:3)        : four-momentum of particle 2            
C
      IMPLICIT NONE
      REAL*8    P1(0:3),P2(0:3),
     &        ESUM,MASS1,MASS2,COSTH1,PHI1,MD2,ED,PP,SINTH1
C
      MD2=(MASS1-MASS2)*(MASS1+MASS2)
      ED=MD2/ESUM
      IF (MASS1*MASS2.EQ.0.) THEN
      PP=(ESUM-ABS(ED))*0.5D0
C
      ELSE
      PP=SQRT((MD2/ESUM)**2-2.0D0*(MASS1**2+MASS2**2)+ESUM**2)*0.5D0
      ENDIF
      SINTH1=SQRT((1.0D0-COSTH1)*(1.0D0+COSTH1))
C
      P1(0) = MAX((ESUM+ED)*0.5D0,0.D0)
      P1(1) = PP*SINTH1*COS(PHI1)
      P1(2) = PP*SINTH1*SIN(PHI1)
      P1(3) = PP*COSTH1
C
      P2(0) = MAX((ESUM-ED)*0.5D0,0.D0)
      P2(1) = -P1(1)
      P2(2) = -P1(2)
      P2(3) = -P1(3)
C
      RETURN
      END
