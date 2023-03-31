C
C ----------------------------------------------------------------------
C
      SUBROUTINE JIOXXX(FI,FO,G,VMASS,VWIDTH , JIO)
C
C this subroutine computes an off-shell vector current from an external 
C fermion pair.  the vector boson propagator is given in feynman gauge  
C for a massless vector and in unitary gauge for a massive vector.      
C                                                                       
C input:                                                                
C       complex fi(6)          : flow-in  fermion                   |fi>
C       complex fo(6)          : flow-out fermion                   <fo|
C       real    g(2)           : coupling constants                  gvf
C       real    vmass          : mass  of output vector v               
C       real    vwidth         : width of output vector v               
C                                                                       
C output:                                                               
C       complex jio(6)         : vector current          j^mu(<fo|v|fi>)
C
      IMPLICIT NONE
      COMPLEX*16 FI(6),FO(6),JIO(6),C0,C1,C2,C3,CS,D
      REAL*8    G(2),Q(0:3),VMASS,VWIDTH,Q2,VM2,DD
C
      REAL*8 RXZERO, RXONE
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0 )
      COMPLEX*16 CXIMAG
      LOGICAL FIRST
      SAVE CXIMAG,FIRST
      DATA FIRST/.TRUE./
C
C          Fix compilation with g77
      IF(FIRST) THEN
        FIRST=.FALSE.
        CXIMAG=DCMPLX( RXZERO, RXONE )
      ENDIF
C
      JIO(5) = FO(5)-FI(5)
      JIO(6) = FO(6)-FI(6)
C
      Q(0)=DBLE( JIO(5))
      Q(1)=DBLE( JIO(6))
      Q(2)=DIMAG(JIO(6))
      Q(3)=DIMAG(JIO(5))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
      VM2=VMASS**2
C
      IF (VMASS.NE.RXZERO) THEN
C
         D=RXONE/DCMPLX( Q2-VM2 , MAX(SIGN( VMASS*VWIDTH ,Q2),RXZERO) )
C  for the running width, use below instead of the above d.
C      d=r_one/dcmplx( q2-vm2 , max( vwidth*q2/vmass ,r_zero) )
C
         IF (G(2).NE.RXZERO) THEN
C
            C0=  G(1)*( FO(3)*FI(1)+FO(4)*FI(2))
     &          +G(2)*( FO(1)*FI(3)+FO(2)*FI(4))
            C1= -G(1)*( FO(3)*FI(2)+FO(4)*FI(1))
     &          +G(2)*( FO(1)*FI(4)+FO(2)*FI(3))
            C2=( G(1)*( FO(3)*FI(2)-FO(4)*FI(1)) 
     &          +G(2)*(-FO(1)*FI(4)+FO(2)*FI(3)))*CXIMAG
            C3=  G(1)*(-FO(3)*FI(1)+FO(4)*FI(2))
     &          +G(2)*( FO(1)*FI(3)-FO(2)*FI(4))
         ELSE
C
            D=D*G(1)
            C0=  FO(3)*FI(1)+FO(4)*FI(2)
            C1= -FO(3)*FI(2)-FO(4)*FI(1)
            C2=( FO(3)*FI(2)-FO(4)*FI(1))*CXIMAG
            C3= -FO(3)*FI(1)+FO(4)*FI(2)
         END IF
C
         CS=(Q(0)*C0-Q(1)*C1-Q(2)*C2-Q(3)*C3)/VM2
C
         JIO(1) = (C0-CS*Q(0))*D
         JIO(2) = (C1-CS*Q(1))*D
         JIO(3) = (C2-CS*Q(2))*D
         JIO(4) = (C3-CS*Q(3))*D
C
      ELSE
         DD=RXONE/Q2
C
         IF (G(2).NE.RXZERO) THEN
            JIO(1) = ( G(1)*( FO(3)*FI(1)+FO(4)*FI(2))
     &                +G(2)*( FO(1)*FI(3)+FO(2)*FI(4)) )*DD
            JIO(2) = (-G(1)*( FO(3)*FI(2)+FO(4)*FI(1))
     &                +G(2)*( FO(1)*FI(4)+FO(2)*FI(3)) )*DD
            JIO(3) = ( G(1)*( FO(3)*FI(2)-FO(4)*FI(1))
     &                +G(2)*(-FO(1)*FI(4)+FO(2)*FI(3)))
     $           *DCMPLX(RXZERO,DD)
            JIO(4) = ( G(1)*(-FO(3)*FI(1)+FO(4)*FI(2))
     &                +G(2)*( FO(1)*FI(3)-FO(2)*FI(4)) )*DD
C
         ELSE
            DD=DD*G(1)
C
            JIO(1) =  ( FO(3)*FI(1)+FO(4)*FI(2))*DD
            JIO(2) = -( FO(3)*FI(2)+FO(4)*FI(1))*DD
            JIO(3) =  ( FO(3)*FI(2)-FO(4)*FI(1))*DCMPLX(RXZERO,DD)
            JIO(4) =  (-FO(3)*FI(1)+FO(4)*FI(2))*DD
         END IF
      END IF
C
      RETURN
      END
