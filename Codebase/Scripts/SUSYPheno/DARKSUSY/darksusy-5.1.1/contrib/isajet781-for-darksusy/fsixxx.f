C
C ----------------------------------------------------------------------
C
      SUBROUTINE FSIXXX(FI,SC,GC,FMASS,FWIDTH , FSI)
C
C this subroutine computes an off-shell fermion wavefunction from a     
C flowing-in external fermion and a vector boson.                       
C                                                                       
C input:                                                                
C       complex*16 fi(6)          : flow-in  fermion                   |fi>
C       complex*16 sc(3)          : input    scalar                      s 
C       complex*16 gc(2)          : coupling constants                 gchf
C       real*8    fmass          : mass  of output fermion f'             
C       real*8    fwidth         : width of output fermion f'             
C                                                                       
C output:                                                               
C       complex fsi(6)         : off-shell fermion             |f',s,fi>
C
      IMPLICIT NONE
      COMPLEX*16 FI(6),SC(3),FSI(6),GC(2),SL1,SL2,SR1,SR2,DS
      REAL*8     PF(0:3),FMASS,FWIDTH,PF2,P0P3,P0M3
C
      FSI(5) = FI(5)-SC(2)
      FSI(6) = FI(6)-SC(3)
C
      PF(0)=DBLE( FSI(5))
      PF(1)=DBLE( FSI(6))
      PF(2)=DIMAG(FSI(6))
      PF(3)=DIMAG(FSI(5))
      PF2=PF(0)**2-(PF(1)**2+PF(2)**2+PF(3)**2)
C
      DS=-SC(1)/DCMPLX(PF2-FMASS**2,MAX(DSIGN(FMASS*FWIDTH ,PF2),0D0))
      P0P3=PF(0)+PF(3)
      P0M3=PF(0)-PF(3)
      SL1=GC(1)*(P0P3*FI(1)+DCONJG(FSI(6))*FI(2))
      SL2=GC(1)*(P0M3*FI(2)      +FSI(6) *FI(1))
      SR1=GC(2)*(P0M3*FI(3)-DCONJG(FSI(6))*FI(4))
      SR2=GC(2)*(P0P3*FI(4)      -FSI(6) *FI(3))
C
      FSI(1) = ( GC(1)*FMASS*FI(1) + SR1 )*DS
      FSI(2) = ( GC(1)*FMASS*FI(2) + SR2 )*DS
      FSI(3) = ( GC(2)*FMASS*FI(3) + SL1 )*DS
      FSI(4) = ( GC(2)*FMASS*FI(4) + SL2 )*DS
C
      RETURN          
      END
