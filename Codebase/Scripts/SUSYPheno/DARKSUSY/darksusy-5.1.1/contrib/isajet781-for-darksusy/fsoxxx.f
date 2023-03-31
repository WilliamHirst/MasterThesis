C
C ----------------------------------------------------------------------
C
      SUBROUTINE FSOXXX(FO,SC,GC,FMASS,FWIDTH , FSO)
C
C this subroutine computes an off-shell fermion wavefunction from a     
C flowing-out external fermion and a vector boson.                      
C                                                                       
C input:                                                                
C       complex*16 fo(6)          : flow-out fermion                   <fo|
C       complex*16 sc(6)          : input    scalar                      s 
C       complex*16 gc(2)          : coupling constants                 gchf
C       real*8     fmass          : mass  of output fermion f'             
C       real*8     fwidth         : width of output fermion f'             
C                                                                       
C output:                                                               
C       complex fso(6)         : off-shell fermion             <fo,s,f'|
C
      IMPLICIT NONE
      COMPLEX*16 FO(6),SC(6),FSO(6),GC(2),SL1,SL2,SR1,SR2,DS
      REAL*8     PF(0:3),FMASS,FWIDTH,PF2,P0P3,P0M3
C
      FSO(5) = FO(5)+SC(2)
      FSO(6) = FO(6)+SC(3)
C
      PF(0)=DBLE( FSO(5))
      PF(1)=DBLE( FSO(6))
      PF(2)=DIMAG(FSO(6))
      PF(3)=DIMAG(FSO(5))
      PF2=PF(0)**2-(PF(1)**2+PF(2)**2+PF(3)**2)
C
      DS=-SC(1)/DCMPLX(PF2-FMASS**2,MAX(DSIGN(FMASS*FWIDTH ,PF2),0D0))
      P0P3=PF(0)+PF(3)
      P0M3=PF(0)-PF(3)
      SL1=GC(2)*(P0P3*FO(3)      +FSO(6) *FO(4))
      SL2=GC(2)*(P0M3*FO(4)+DCONJG(FSO(6))*FO(3))
      SR1=GC(1)*(P0M3*FO(1)      -FSO(6) *FO(2))
      SR2=GC(1)*(P0P3*FO(2)-DCONJG(FSO(6))*FO(1))
C
      FSO(1) = ( GC(1)*FMASS*FO(1) + SL1 )*DS
      FSO(2) = ( GC(1)*FMASS*FO(2) + SL2 )*DS
      FSO(3) = ( GC(2)*FMASS*FO(3) + SR1 )*DS
      FSO(4) = ( GC(2)*FMASS*FO(4) + SR2 )*DS
C
      RETURN          
      END
