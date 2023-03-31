C
C ----------------------------------------------------------------------
C
      SUBROUTINE JW3WXX(W1,W2,W3,G1,G2,WMASS,WWIDTH,VMASS,VWIDTH , JW3W)
C
C this subroutine computes an off-shell w+, w-, w3, z or photon current 
C from the four-point gauge boson coupling, including the contributions 
C of w exchange diagrams.  the vector propagator is given in feynman    
C gauge for a photon and in unitary gauge for w and z bosons.  if one   
C sets wmass=0.0, then the ggg-->g current is given (see sect 2.9.1 of 
C the manual).                                                          
C                                                                       
C input:                                                                
C       complex w1(6)          : first  vector                        w1
C       complex w2(6)          : second vector                        w2
C       complex w3(6)          : third  vector                        w3
C       real    g1             : first  coupling constant               
C       real    g2             : second coupling constant               
C                                                  (see the table below)
C       real    wmass          : mass  of internal w                    
C       real    wwidth         : width of internal w                    
C       real    vmass          : mass  of output w'                     
C       real    vwidth         : width of output w'                     
C                                                                       
C the possible sets of the inputs are as follows:                       
C   ------------------------------------------------------------------- 
C   |  w1  |  w2  |  w3  | g1 | g2 |wmass|wwidth|vmass|vwidth || jw3w | 
C   ------------------------------------------------------------------- 
C   |  w-  |  w3  |  w+  | gw |gwwz|wmass|wwidth|zmass|zwidth ||  z   | 
C   |  w-  |  w3  |  w+  | gw |gwwa|wmass|wwidth|  0. |  0.   ||  a   | 
C   |  w-  |  z   |  w+  |gwwz|gwwz|wmass|wwidth|zmass|zwidth ||  z   | 
C   |  w-  |  z   |  w+  |gwwz|gwwa|wmass|wwidth|  0. |  0.   ||  a   | 
C   |  w-  |  a   |  w+  |gwwa|gwwz|wmass|wwidth|zmass|zwidth ||  z   | 
C   |  w-  |  a   |  w+  |gwwa|gwwa|wmass|wwidth|  0. |  0.   ||  a   | 
C   ------------------------------------------------------------------- 
C   |  w3  |  w-  |  w3  | gw | gw |wmass|wwidth|wmass|wwidth ||  w+  | 
C   |  w3  |  w+  |  w3  | gw | gw |wmass|wwidth|wmass|wwidth ||  w-  | 
C   |  w3  |  w-  |  z   | gw |gwwz|wmass|wwidth|wmass|wwidth ||  w+  | 
C   |  w3  |  w+  |  z   | gw |gwwz|wmass|wwidth|wmass|wwidth ||  w-  | 
C   |  w3  |  w-  |  a   | gw |gwwa|wmass|wwidth|wmass|wwidth ||  w+  | 
C   |  w3  |  w+  |  a   | gw |gwwa|wmass|wwidth|wmass|wwidth ||  w-  | 
C   |  z   |  w-  |  z   |gwwz|gwwz|wmass|wwidth|wmass|wwidth ||  w+  | 
C   |  z   |  w+  |  z   |gwwz|gwwz|wmass|wwidth|wmass|wwidth ||  w-  | 
C   |  z   |  w-  |  a   |gwwz|gwwa|wmass|wwidth|wmass|wwidth ||  w+  | 
C   |  z   |  w+  |  a   |gwwz|gwwa|wmass|wwidth|wmass|wwidth ||  w-  | 
C   |  a   |  w-  |  a   |gwwa|gwwa|wmass|wwidth|wmass|wwidth ||  w+  | 
C   |  a   |  w+  |  a   |gwwa|gwwa|wmass|wwidth|wmass|wwidth ||  w-  | 
C   ------------------------------------------------------------------- 
C where all the bosons are defined by the flowing-out quantum number.   
C                                                                       
C output:                                                               
C       complex jw3w(6)        : w current             j^mu(w':w1,w2,w3)
C
      IMPLICIT NONE
      COMPLEX*16  W1(6),W2(6),W3(6),JW3W(6)
      COMPLEX*16 DW1(0:3),DW2(0:3),DW3(0:3),
     &           JJ(0:3),J4(0:3),
     &           DV,W12,W32,W13,
     &           JQ
      REAL*8     G1,G2,WMASS,WWIDTH,VMASS,VWIDTH
      REAL*8     P1(0:3),P2(0:3),P3(0:3),Q(0:3),
     &           DG2,DMV,DWV,MV2,Q2
C
      REAL*8 RXZERO
      PARAMETER( RXZERO=0.0D0 )
C
      JW3W(5) = W1(5)+W2(5)+W3(5)
      JW3W(6) = W1(6)+W2(6)+W3(6)
C
      DW1(0)=DCMPLX(W1(1))
      DW1(1)=DCMPLX(W1(2))
      DW1(2)=DCMPLX(W1(3))
      DW1(3)=DCMPLX(W1(4))
      DW2(0)=DCMPLX(W2(1))
      DW2(1)=DCMPLX(W2(2))
      DW2(2)=DCMPLX(W2(3))
      DW2(3)=DCMPLX(W2(4))
      DW3(0)=DCMPLX(W3(1))
      DW3(1)=DCMPLX(W3(2))
      DW3(2)=DCMPLX(W3(3))
      DW3(3)=DCMPLX(W3(4))
      P1(0)=DBLE(      W1(5))
      P1(1)=DBLE(      W1(6))
      P1(2)=DBLE(DIMAG(W1(6)))
      P1(3)=DBLE(DIMAG(W1(5)))
      P2(0)=DBLE(      W2(5))
      P2(1)=DBLE(      W2(6))
      P2(2)=DBLE(DIMAG(W2(6)))
      P2(3)=DBLE(DIMAG(W2(5)))
      P3(0)=DBLE(      W3(5))
      P3(1)=DBLE(      W3(6))
      P3(2)=DBLE(DIMAG(W3(6)))
      P3(3)=DBLE(DIMAG(W3(5)))
      Q(0)=-(P1(0)+P2(0)+P3(0))
      Q(1)=-(P1(1)+P2(1)+P3(1))
      Q(2)=-(P1(2)+P2(2)+P3(2))
      Q(3)=-(P1(3)+P2(3)+P3(3))
C
      Q2 =Q(0)**2 -(Q(1)**2 +Q(2)**2 +Q(3)**2)
      DG2=DBLE(G1)*DBLE(G2)
      DMV=DBLE(VMASS)
      DWV=DBLE(VWIDTH)
      MV2=DMV**2
      IF (VMASS.EQ. RXZERO) THEN
      DV = 1.0D0/DCMPLX( Q2 )
      ELSE
      DV = 1.0D0/DCMPLX( Q2 -MV2 , DMAX1(DSIGN(DMV*DWV,Q2 ),0.D0) )
      ENDIF
C  for the running width, use below instead of the above dv.
C      dv = 1.0d0/dcmplx( q2 -mv2 , dmax1(dwv*q2/dmv,0.d0) )
C
      W12=DW1(0)*DW2(0)-DW1(1)*DW2(1)-DW1(2)*DW2(2)-DW1(3)*DW2(3)
      W32=DW3(0)*DW2(0)-DW3(1)*DW2(1)-DW3(2)*DW2(2)-DW3(3)*DW2(3)
C
      IF ( WMASS .NE. RXZERO ) THEN
         W13=DW1(0)*DW3(0)-DW1(1)*DW3(1)-DW1(2)*DW3(2)-DW1(3)*DW3(3)
C
         J4(0)=DG2*( DW1(0)*W32 + DW3(0)*W12 - 2.D0*DW2(0)*W13 )
         J4(1)=DG2*( DW1(1)*W32 + DW3(1)*W12 - 2.D0*DW2(1)*W13 )
         J4(2)=DG2*( DW1(2)*W32 + DW3(2)*W12 - 2.D0*DW2(2)*W13 )
         J4(3)=DG2*( DW1(3)*W32 + DW3(3)*W12 - 2.D0*DW2(3)*W13 )
C
         JJ(0)=J4(0)
         JJ(1)=J4(1)
         JJ(2)=J4(2)
         JJ(3)=J4(3)
C
      ELSE
C
         W12=DW1(0)*DW2(0)-DW1(1)*DW2(1)-DW1(2)*DW2(2)-DW1(3)*DW2(3)
         W32=DW3(0)*DW2(0)-DW3(1)*DW2(1)-DW3(2)*DW2(2)-DW3(3)*DW2(3)
         W13=DW1(0)*DW3(0)-DW1(1)*DW3(1)-DW1(2)*DW3(2)-DW1(3)*DW3(3)
C
         J4(0)=DG2*( DW1(0)*W32 - DW2(0)*W13 )
         J4(1)=DG2*( DW1(1)*W32 - DW2(1)*W13 )
         J4(2)=DG2*( DW1(2)*W32 - DW2(2)*W13 )
         J4(3)=DG2*( DW1(3)*W32 - DW2(3)*W13 )
C
         JJ(0)=J4(0)
         JJ(1)=J4(1)
         JJ(2)=J4(2)
         JJ(3)=J4(3)
C
      END IF
C
      IF ( VMASS .NE. RXZERO ) THEN
C
         JQ=(JJ(0)*Q(0)-JJ(1)*Q(1)-JJ(2)*Q(2)-JJ(3)*Q(3))/MV2
C
         JW3W(1) = DCMPLX( (JJ(0)-JQ*Q(0))*DV )
         JW3W(2) = DCMPLX( (JJ(1)-JQ*Q(1))*DV )
         JW3W(3) = DCMPLX( (JJ(2)-JQ*Q(2))*DV )
         JW3W(4) = DCMPLX( (JJ(3)-JQ*Q(3))*DV )
C
      ELSE
C
         JW3W(1) = DCMPLX( JJ(0)*DV )
         JW3W(2) = DCMPLX( JJ(1)*DV )
         JW3W(3) = DCMPLX( JJ(2)*DV )
         JW3W(4) = DCMPLX( JJ(3)*DV )
      END IF
C
      RETURN
      END
