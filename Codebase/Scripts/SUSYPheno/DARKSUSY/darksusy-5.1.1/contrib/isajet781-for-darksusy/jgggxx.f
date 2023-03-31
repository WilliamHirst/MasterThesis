C
C
C ----------------------------------------------------------------------
C
      SUBROUTINE JGGGXX(W1,W2,W3,G, JW3W)
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
C       real    g             : first  coupling constant               
C                                                  (see the table below)
C                                                                       
C output:                                                               
C       complex jw3w(6)        : w current             j^mu(w':w1,w2,w3)
C
      IMPLICIT NONE
      COMPLEX*16  W1(6),W2(6),W3(6),JW3W(6)
      COMPLEX*16 DW1(0:3),DW2(0:3),DW3(0:3),
     &           JJ(0:3),DV,W32,W13
      REAL*8     P1(0:3),P2(0:3),P3(0:3),Q(0:3),G,DG2,Q2
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
C
      DG2=DBLE(G)*DBLE(G)
C
      DV = 1.0D0/DCMPLX( Q2 )
C
C  for the running width, use below instead of the above dv.
C      dv = 1.0d0/dcmplx( q2 -mv2 , dmax1(dwv*q2/dmv,0.d0) )
C
      W32=DW3(0)*DW2(0)-DW3(1)*DW2(1)-DW3(2)*DW2(2)-DW3(3)*DW2(3)
C
C     
      W13=DW1(0)*DW3(0)-DW1(1)*DW3(1)-DW1(2)*DW3(2)-DW1(3)*DW3(3)
C     
      JJ(0)=DG2*( DW1(0)*W32 - DW2(0)*W13 )
      JJ(1)=DG2*( DW1(1)*W32 - DW2(1)*W13 )
      JJ(2)=DG2*( DW1(2)*W32 - DW2(2)*W13 )
      JJ(3)=DG2*( DW1(3)*W32 - DW2(3)*W13 )
C     
      JW3W(1) = DCMPLX( JJ(0)*DV )
      JW3W(2) = DCMPLX( JJ(1)*DV )
      JW3W(3) = DCMPLX( JJ(2)*DV )
      JW3W(4) = DCMPLX( JJ(3)*DV )
C
      RETURN
      END
