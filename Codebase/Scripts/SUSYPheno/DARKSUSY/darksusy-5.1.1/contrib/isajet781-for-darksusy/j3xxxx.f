C
C ----------------------------------------------------------------------
C
      SUBROUTINE J3XXXX(FI,FO,GAF,GZF,ZMASS,ZWIDTH , J3)
C
C this subroutine computes the sum of photon and z currents with the    
C suitable weights ( j(w3) = cos(theta_w) j(z) + sin(theta_w) j(a) ).   
C the output j3 is useful as an input of vvvxxx, jvvxxx or w3w3xx.      
C the photon propagator is given in feynman gauge, and the z propagator 
C is given in unitary gauge.                                            
C                                                                       
C input:                                                                
C       complex fi(6)          : flow-in  fermion                   |fi>
C       complex fo(6)          : flow-out fermion                   <fo|
C       real    gaf(2)         : fi couplings with a                 gaf
C       real    gzf(2)         : fi couplings with z                 gzf
C       real    zmass          : mass  of z                             
C       real    zwidth         : width of z                             
C                                                                       
C output:                                                               
C       complex j3(6)          : w3 current             j^mu(<fo|w3|fi>)
C
      IMPLICIT NONE
      COMPLEX*16 FI(6),FO(6),J3(6),
     &        C0L,C1L,C2L,C3L,CSL,C0R,C1R,C2R,C3R,CSR,DZ,DDIF
      REAL*8    GAF(2),GZF(2),Q(0:3),ZMASS,ZWIDTH,ZM2,ZMW,Q2,DA,WW,
     &        CW,SW,GN,GZ3L,GA3L
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
      J3(5) = FO(5)-FI(5)
      J3(6) = FO(6)-FI(6)
C
      Q(0)=-DBLE( J3(5))
      Q(1)=-DBLE( J3(6))
      Q(2)=-DIMAG(J3(6))
      Q(3)=-DIMAG(J3(5))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
      ZM2=ZMASS**2
      ZMW=ZMASS*ZWIDTH
C
      DA=RXONE/Q2
      WW=MAX(DSIGN( ZMW ,Q2),RXZERO)
      DZ=RXONE/DCMPLX( Q2-ZM2 , WW )
      DDIF=DCMPLX( -ZM2 , WW )*DA*DZ
C
C ddif is the difference : ddif=da-dz
C  for the running width, use below instead of the above ww,dz and ddif.
C      ww=max( zwidth*q2/zmass ,r_zero)
C      dz=r_one/dcmplx( q2-zm2 , ww )
C      ddif=dcmplx( -zm2 , ww )*da*dz
C
      CW=RXONE/SQRT(RXONE+(GZF(2)/GAF(2))**2)
      SW=SQRT((RXONE-CW)*(RXONE+CW))
      GN=GAF(2)*SW
      GZ3L=GZF(1)*CW
      GA3L=GAF(1)*SW
      C0L=  FO(3)*FI(1)+FO(4)*FI(2)
      C0R=  FO(1)*FI(3)+FO(2)*FI(4)
      C1L=-(FO(3)*FI(2)+FO(4)*FI(1))
      C1R=  FO(1)*FI(4)+FO(2)*FI(3)
      C2L= (FO(3)*FI(2)-FO(4)*FI(1))*CXIMAG
      C2R=(-FO(1)*FI(4)+FO(2)*FI(3))*CXIMAG
      C3L= -FO(3)*FI(1)+FO(4)*FI(2)
      C3R=  FO(1)*FI(3)-FO(2)*FI(4)
      CSL=(Q(0)*C0L-Q(1)*C1L-Q(2)*C2L-Q(3)*C3L)/ZM2
      CSR=(Q(0)*C0R-Q(1)*C1R-Q(2)*C2R-Q(3)*C3R)/ZM2
C
      J3(1) =  GZ3L*DZ*(C0L-CSL*Q(0))+GA3L*C0L*DA
     &       + GN*(C0R*DDIF-CSR*Q(0)*DZ)
      J3(2) =  GZ3L*DZ*(C1L-CSL*Q(1))+GA3L*C1L*DA
     &       + GN*(C1R*DDIF-CSR*Q(1)*DZ)
      J3(3) =  GZ3L*DZ*(C2L-CSL*Q(2))+GA3L*C2L*DA
     &       + GN*(C2R*DDIF-CSR*Q(2)*DZ)
      J3(4) =  GZ3L*DZ*(C3L-CSL*Q(3))+GA3L*C3L*DA
     &       + GN*(C3R*DDIF-CSR*Q(3)*DZ)
C
      RETURN
      END
