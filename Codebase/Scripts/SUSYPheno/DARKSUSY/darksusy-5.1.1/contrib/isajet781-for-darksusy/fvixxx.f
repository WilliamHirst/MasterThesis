C
C ----------------------------------------------------------------------
C
      SUBROUTINE FVIXXX(FI,VC,G,FMASS,FWIDTH , FVI)
C
C this subroutine computes an off-shell fermion wavefunction from a     
C flowing-in external fermion and a vector boson.                       
C                                                                       
C input:                                                                
C       complex fi(6)          : flow-in  fermion                   |fi>
C       complex vc(6)          : input    vector                      v 
C       real    g(2)           : coupling constants                  gvf
C       real    fmass          : mass  of output fermion f'             
C       real    fwidth         : width of output fermion f'             
C                                                                       
C output:                                                               
C       complex fvi(6)         : off-shell fermion             |f',v,fi>
C
      IMPLICIT NONE
      COMPLEX*16 FI(6),VC(6),FVI(6),SL1,SL2,SR1,SR2,D
      REAL*8    G(2),PF(0:3),FMASS,FWIDTH,PF2
      REAL*8 RXZERO, RXONE
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0 )
      COMPLEX*16 CXIMAG
C
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
      FVI(5) = FI(5)-VC(5)
      FVI(6) = FI(6)-VC(6)
C
      PF(0)=DBLE( FVI(5))
      PF(1)=DBLE( FVI(6))
      PF(2)=DIMAG(FVI(6))
      PF(3)=DIMAG(FVI(5))
      PF2=PF(0)**2-(PF(1)**2+PF(2)**2+PF(3)**2)
C
      D=-RXONE/DCMPLX( PF2-FMASS**2,MAX(SIGN(FMASS*FWIDTH,PF2),RXZERO))
      SL1= (VC(1)+       VC(4))*FI(1)
     &    +(VC(2)-CXIMAG*VC(3))*FI(2)
      SL2= (VC(2)+CXIMAG*VC(3))*FI(1)
     &    +(VC(1)-       VC(4))*FI(2)
C
      IF ( G(2) .NE. RXZERO ) THEN
         SR1= (VC(1)-       VC(4))*FI(3)
     &       -(VC(2)-CXIMAG*VC(3))*FI(4)
         SR2=-(VC(2)+CXIMAG*VC(3))*FI(3)
     &       +(VC(1)+       VC(4))*FI(4)
C
         FVI(1) = ( G(1)*((PF(0)-PF(3))*SL1 -DCONJG(FVI(6))*SL2)
     &             +G(2)*FMASS*SR1)*D
         FVI(2) = ( G(1)*(      -FVI(6)*SL1 +(PF(0)+PF(3))*SL2)
     &             +G(2)*FMASS*SR2)*D
         FVI(3) = ( G(2)*((PF(0)+PF(3))*SR1 +DCONJG(FVI(6))*SR2)
     &             +G(1)*FMASS*SL1)*D
         FVI(4) = ( G(2)*(       FVI(6)*SR1 +(PF(0)-PF(3))*SR2)
     &             +G(1)*FMASS*SL2)*D
C
      ELSE          
         FVI(1) = G(1)*((PF(0)-PF(3))*SL1 -DCONJG(FVI(6))*SL2)*D
         FVI(2) = G(1)*(      -FVI(6)*SL1 +(PF(0)+PF(3))*SL2)*D
         FVI(3) = G(1)*FMASS*SL1*D
         FVI(4) = G(1)*FMASS*SL2*D
      END IF
C
      RETURN          
      END
