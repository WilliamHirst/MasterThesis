C
C ----------------------------------------------------------------------
C
      SUBROUTINE JVVXXX(V1,V2,G,VMASS,VWIDTH , JVV)
C
C this subroutine computes an off-shell vector current from the three-  
C point gauge boson coupling.  the vector propagator is given in feynman
C gauge for a massless vector and in unitary gauge for a massive vector.
C                                                                       
C input:                                                                
C       complex v1(6)          : first  vector                        v1
C       complex v2(6)          : second vector                        v2
C       real    g              : coupling constant (see the table below)
C       real    vmass          : mass  of output vector v               
C       real    vwidth         : width of output vector v               
C                                                                       
C the possible sets of the inputs are as follows:                       
C    ------------------------------------------------------------------ 
C    |   v1   |   v2   |  jvv   |      g       |   vmass  |  vwidth   | 
C    ------------------------------------------------------------------ 
C    |   w-   |   w+   |  a/z   |  gwwa/gwwz   | 0./zmass | 0./zwidth | 
C    | w3/a/z |   w-   |  w+    | gw/gwwa/gwwz |   wmass  |  wwidth   | 
C    |   w+   | w3/a/z |  w-    | gw/gwwa/gwwz |   wmass  |  wwidth   | 
C    ------------------------------------------------------------------ 
C where all the bosons are defined by the flowing-out quantum number.   
C                                                                       
C output:                                                               
C       complex jvv(6)         : vector current            j^mu(v:v1,v2)
C
      IMPLICIT NONE
      COMPLEX*16 V1(6),V2(6),JVV(6),J12(0:3),JS,DG,
     &        SV1,SV2,S11,S12,S21,S22,V12
      REAL*8    P1(0:3),P2(0:3),Q(0:3),G,VMASS,VWIDTH,GS,S,VM2,M1,M2
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
      IF ( VMASS .NE. RXZERO ) THEN
         VM2=VMASS**2
         M1=P1(0)**2-(P1(1)**2+P1(2)**2+P1(3)**2)
         M2=P2(0)**2-(P2(1)**2+P2(2)**2+P2(3)**2)
         S11=P1(0)*V1(1)-P1(1)*V1(2)-P1(2)*V1(3)-P1(3)*V1(4)
         S12=P1(0)*V2(1)-P1(1)*V2(2)-P1(2)*V2(3)-P1(3)*V2(4)
         S21=P2(0)*V1(1)-P2(1)*V1(2)-P2(2)*V1(3)-P2(3)*V1(4)
         S22=P2(0)*V2(1)-P2(1)*V2(2)-P2(2)*V2(3)-P2(3)*V2(4)
         JS=(V12*(-M1+M2) +S11*S12 -S21*S22)/VM2
C
         DG=-G/DCMPLX( S-VM2 , MAX(SIGN( VMASS*VWIDTH ,S),RXZERO) )
C
C  for the running width, use below instead of the above dg.
C         dg=-g/dcmplx( s-vm2 , max( vwidth*s/vmass ,r_zero) )
C
         JVV(1) = DG*(J12(0)-Q(0)*JS)
         JVV(2) = DG*(J12(1)-Q(1)*JS)
         JVV(3) = DG*(J12(2)-Q(2)*JS)
         JVV(4) = DG*(J12(3)-Q(3)*JS)
C
      ELSE
         GS=-G/S
C
         JVV(1) = GS*J12(0)
         JVV(2) = GS*J12(1)
         JVV(3) = GS*J12(2)
         JVV(4) = GS*J12(3)
      END IF
C
      RETURN
      END
