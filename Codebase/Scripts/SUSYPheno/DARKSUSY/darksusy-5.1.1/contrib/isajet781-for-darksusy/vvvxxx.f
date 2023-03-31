C
C ======================================================================
C
      SUBROUTINE VVVXXX(WM,WP,W3,G , VERTEX)
C
C this subroutine computes an amplitude of the three-point coupling of  
C the gauge bosons.                                                     
C                                                                       
C input:                                                                
C       complex wm(6)          : vector               flow-out w-       
C       complex wp(6)          : vector               flow-out w+       
C       complex w3(6)          : vector               j3 or a    or z   
C       real    g              : coupling constant    gw or gwwa or gwwz
C                                                                       
C output:                                                               
C       complex vertex         : amplitude               gamma(wm,wp,w3)
C
      IMPLICIT NONE
      COMPLEX*16 WM(6),WP(6),W3(6),VERTEX, 
     &        XV1,XV2,XV3,V12,V23,V31,P12,P13,P21,P23,P31,P32
      REAL*8    PWM(0:3),PWP(0:3),PW3(0:3),G
C
      REAL*8 RXZERO, RTENTH
      PARAMETER( RXZERO=0.0D0, RTENTH=0.1D0 )
C
      PWM(0)=DBLE( WM(5))
      PWM(1)=DBLE( WM(6))
      PWM(2)=DIMAG(WM(6))
      PWM(3)=DIMAG(WM(5))
      PWP(0)=DBLE( WP(5))
      PWP(1)=DBLE( WP(6))
      PWP(2)=DIMAG(WP(6))
      PWP(3)=DIMAG(WP(5))
      PW3(0)=DBLE( W3(5))
      PW3(1)=DBLE( W3(6))
      PW3(2)=DIMAG(W3(6))
      PW3(3)=DIMAG(W3(5))
C
      V12=WM(1)*WP(1)-WM(2)*WP(2)-WM(3)*WP(3)-WM(4)*WP(4)
      V23=WP(1)*W3(1)-WP(2)*W3(2)-WP(3)*W3(3)-WP(4)*W3(4)
      V31=W3(1)*WM(1)-W3(2)*WM(2)-W3(3)*WM(3)-W3(4)*WM(4)
      XV1=RXZERO
      XV2=RXZERO
      XV3=RXZERO
      IF ( ABS(WM(1)) .NE. RXZERO ) THEN
         IF (ABS(WM(1)).GE.MAX(ABS(WM(2)),ABS(WM(3)),ABS(WM(4)))
     $        *RTENTH)
     &      XV1=PWM(0)/WM(1)
      ENDIF
      IF ( ABS(WP(1)) .NE. RXZERO) THEN
         IF (ABS(WP(1)).GE.MAX(ABS(WP(2)),ABS(WP(3)),ABS(WP(4)))
     $        *RTENTH)
     &      XV2=PWP(0)/WP(1)
      ENDIF
      IF ( ABS(W3(1)) .NE. RXZERO) THEN
         IF ( ABS(W3(1)).GE.MAX(ABS(W3(2)),ABS(W3(3)),ABS(W3(4)))
     $        *RTENTH)
     &      XV3=PW3(0)/W3(1)
      ENDIF
      P12= (PWM(0)-XV1*WM(1))*WP(1)-(PWM(1)-XV1*WM(2))*WP(2)
     &    -(PWM(2)-XV1*WM(3))*WP(3)-(PWM(3)-XV1*WM(4))*WP(4)
      P13= (PWM(0)-XV1*WM(1))*W3(1)-(PWM(1)-XV1*WM(2))*W3(2)
     &    -(PWM(2)-XV1*WM(3))*W3(3)-(PWM(3)-XV1*WM(4))*W3(4)
      P21= (PWP(0)-XV2*WP(1))*WM(1)-(PWP(1)-XV2*WP(2))*WM(2)
     &    -(PWP(2)-XV2*WP(3))*WM(3)-(PWP(3)-XV2*WP(4))*WM(4)
      P23= (PWP(0)-XV2*WP(1))*W3(1)-(PWP(1)-XV2*WP(2))*W3(2)
     &    -(PWP(2)-XV2*WP(3))*W3(3)-(PWP(3)-XV2*WP(4))*W3(4)
      P31= (PW3(0)-XV3*W3(1))*WM(1)-(PW3(1)-XV3*W3(2))*WM(2)
     &    -(PW3(2)-XV3*W3(3))*WM(3)-(PW3(3)-XV3*W3(4))*WM(4)
      P32= (PW3(0)-XV3*W3(1))*WP(1)-(PW3(1)-XV3*W3(2))*WP(2)
     &    -(PW3(2)-XV3*W3(3))*WP(3)-(PW3(3)-XV3*W3(4))*WP(4)
C
      VERTEX = -(V12*(P13-P23)+V23*(P21-P31)+V31*(P32-P12))*G
C
      RETURN
      END
