C
C ----------------------------------------------------------------------
C
      SUBROUTINE JGGXXX(V1,V2,G, JVV)
C
C this subroutine computes an off-shell vector current from the three-  
C point gauge boson coupling.  the vector propagator is given in feynman
C gauge for a massless vector and in unitary gauge for a massive vector.
C                                                                       
C input:                                                                
C       complex v1(6)          : first  vector                        v1
C       complex v2(6)          : second vector                        v2
C       real    g              : coupling constant (see the table below)
C                                                                       
C output:                                                               
C       complex jvv(6)         : vector current            j^mu(v:v1,v2)
C
      IMPLICIT NONE
      COMPLEX*16 V1(6),V2(6),JVV(6),J12(0:3),
     &        SV1,SV2,V12
      REAL*8    P1(0:3),P2(0:3),Q(0:3),G,GS,S
C
      REAL*8 RXZERO
      PARAMETER( RXZERO=0.0D0 )
C
      JVV(5) = V1(5)+V2(5)
      JVV(6) = V1(6)+V2(6)
C
      P1(0)=DBLE( V1(5))
      P1(1)=DBLE( V1(6))
      P1(2)=DIMAG(V1(6))
      P1(3)=DIMAG(V1(5))
      P2(0)=DBLE( V2(5))
      P2(1)=DBLE( V2(6))
      P2(2)=DIMAG(V2(6))
      P2(3)=DIMAG(V2(5))
      Q(0)=-DBLE( JVV(5))
      Q(1)=-DBLE( JVV(6))
      Q(2)=-DIMAG(JVV(6))
      Q(3)=-DIMAG(JVV(5))
      S=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
C
      V12=V1(1)*V2(1)-V1(2)*V2(2)-V1(3)*V2(3)-V1(4)*V2(4)
      SV1= (P2(0)-Q(0))*V1(1) -(P2(1)-Q(1))*V1(2)
     &    -(P2(2)-Q(2))*V1(3) -(P2(3)-Q(3))*V1(4)
      SV2=-(P1(0)-Q(0))*V2(1) +(P1(1)-Q(1))*V2(2)
     &    +(P1(2)-Q(2))*V2(3) +(P1(3)-Q(3))*V2(4)
      J12(0)=(P1(0)-P2(0))*V12 +SV1*V2(1) +SV2*V1(1)
      J12(1)=(P1(1)-P2(1))*V12 +SV1*V2(2) +SV2*V1(2)
      J12(2)=(P1(2)-P2(2))*V12 +SV1*V2(3) +SV2*V1(3)
      J12(3)=(P1(3)-P2(3))*V12 +SV1*V2(4) +SV2*V1(4)
C
      GS=-G/S
C
      JVV(1) = GS*J12(0)
      JVV(2) = GS*J12(1)
      JVV(3) = GS*J12(2)
      JVV(4) = GS*J12(3)
C
      RETURN
      END
