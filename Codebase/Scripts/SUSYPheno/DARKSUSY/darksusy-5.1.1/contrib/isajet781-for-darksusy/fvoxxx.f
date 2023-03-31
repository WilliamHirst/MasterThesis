C
C ----------------------------------------------------------------------
C
      SUBROUTINE FVOXXX(FO,VC,G,FMASS,FWIDTH , FVO)
C
C this subroutine computes an off-shell fermion wavefunction from a     
C flowing-out external fermion and a vector boson.                      
C                                                                       
C input:                                                                
C       complex fo(6)          : flow-out fermion                   <fo|
C       complex vc(6)          : input    vector                      v 
C       real    g(2)           : coupling constants                  gvf
C       real    fmass          : mass  of output fermion f'             
C       real    fwidth         : width of output fermion f'             
C                                                                       
C output:                                                               
C       complex fvo(6)         : off-shell fermion             <fo,v,f'|
C
      IMPLICIT NONE
      COMPLEX*16 FO(6),VC(6),FVO(6),SL1,SL2,SR1,SR2,D
      REAL*8    G(2),PF(0:3),FMASS,FWIDTH,PF2
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
      FVO(5) = FO(5)+VC(5)
      FVO(6) = FO(6)+VC(6)
C
      PF(0)=DBLE( FVO(5))
      PF(1)=DBLE( FVO(6))
      PF(2)=DIMAG(FVO(6))
      PF(3)=DIMAG(FVO(5))
      PF2=PF(0)**2-(PF(1)**2+PF(2)**2+PF(3)**2)
C
      D=-RXONE/DCMPLX( PF2-FMASS**2,MAX(SIGN(FMASS*FWIDTH,PF2),RXZERO))
      SL1= (VC(1)+       VC(4))*FO(3)
     &    +(VC(2)+CXIMAG*VC(3))*FO(4)
      SL2= (VC(2)-CXIMAG*VC(3))*FO(3)
     &    +(VC(1)-       VC(4))*FO(4)
C
      IF ( G(2) .NE. RXZERO ) THEN
         SR1= (VC(1)-       VC(4))*FO(1)
     &       -(VC(2)+CXIMAG*VC(3))*FO(2)
         SR2=-(VC(2)-CXIMAG*VC(3))*FO(1)
     &       +(VC(1)+       VC(4))*FO(2)
C
         FVO(1) = ( G(2)*( (PF(0)+PF(3))*SR1        +FVO(6)*SR2)
     &             +G(1)*FMASS*SL1)*D
         FVO(2) = ( G(2)*( DCONJG(FVO(6))*SR1 +(PF(0)-PF(3))*SR2)
     &             +G(1)*FMASS*SL2)*D
         FVO(3) = ( G(1)*( (PF(0)-PF(3))*SL1        -FVO(6)*SL2)
     &             +G(2)*FMASS*SR1)*D
         FVO(4) = ( G(1)*(-DCONJG(FVO(6))*SL1 +(PF(0)+PF(3))*SL2)
     &             +G(2)*FMASS*SR2)*D
C
      ELSE
         FVO(1) = G(1)*FMASS*SL1*D
         FVO(2) = G(1)*FMASS*SL2*D
         FVO(3) = G(1)*( (PF(0)-PF(3))*SL1        -FVO(6)*SL2)*D
         FVO(4) = G(1)*(-DCONJG(FVO(6))*SL1 +(PF(0)+PF(3))*SL2)*D
      END IF
C
      RETURN          
      END
