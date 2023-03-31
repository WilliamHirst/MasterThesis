C
C ----------------------------------------------------------------------
C
      SUBROUTINE GGGGXX(WM,W31,WP,W32,G, VERTEX)
C
C this subroutine computes an amplitude of the four-point coupling of   
C the w-, w+ and two w3/z/a.  the amplitude includes the contributions  
C of w exchange diagrams.  the internal w propagator is given in unitary
C gauge.  if one sets wmass=0.0, then the gggg vertex is given (see sect
C 2.9.1 of the manual).
C                                                                       
C input:                                                                
C       complex wm(0:3)        : flow-out w-                         wm 
C       complex w31(0:3)       : first    w3/z/a                     w31
C       complex wp(0:3)        : flow-out w+                         wp 
C       complex w32(0:3)       : second   w3/z/a                     w32
C       real    g              : coupling of w31 with w-/w+             
C                                                  (see the table below)
C                                                                       
C the possible sets of the inputs are as follows:                       
C   -------------------------------------------                         
C   |  wm  |  w31 |  wp  |  w32 |  g31 |  g32 |                         
C   -------------------------------------------                         
C   |  w-  |  w3  |  w+  |  w3  |  gw  |  gw  |                         
C   |  w-  |  w3  |  w+  |  z   |  gw  | gwwz |                         
C   |  w-  |  w3  |  w+  |  a   |  gw  | gwwa |                         
C   |  w-  |  z   |  w+  |  z   | gwwz | gwwz |                         
C   |  w-  |  z   |  w+  |  a   | gwwz | gwwa |                         
C   |  w-  |  a   |  w+  |  a   | gwwa | gwwa |                         
C   -------------------------------------------                         
C where all the bosons are defined by the flowing-out quantum number.   
C                                                                       
C output:                                                               
C       complex vertex         : amplitude          gamma(wm,w31,wp,w32)
C
      IMPLICIT NONE
      COMPLEX*16    WM(6),W31(6),WP(6),W32(6),VERTEX
      COMPLEX*16 DV1(0:3),DV2(0:3),DV3(0:3),DV4(0:3),
     &           DVERTX,V12,V13,V14,V23,V24,V34
      REAL*8       PWM(0:3),PW31(0:3),PWP(0:3),PW32(0:3),G
      REAL*8     DP1(0:3),DP2(0:3),DP3(0:3),DP4(0:3)
      REAL*8 RXZERO, RXONE
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0 )
C
      PWM(0)=DBLE( WM(5))
      PWM(1)=DBLE( WM(6))
      PWM(2)=DIMAG(WM(6))
      PWM(3)=DIMAG(WM(5))
      PWP(0)=DBLE( WP(5))
      PWP(1)=DBLE( WP(6))
      PWP(2)=DIMAG(WP(6))
      PWP(3)=DIMAG(WP(5))
      PW31(0)=DBLE( W31(5))
      PW31(1)=DBLE( W31(6))
      PW31(2)=DIMAG(W31(6))
      PW31(3)=DIMAG(W31(5))
      PW32(0)=DBLE( W32(5))
      PW32(1)=DBLE( W32(6))
      PW32(2)=DIMAG(W32(6))
      PW32(3)=DIMAG(W32(5))
C
      DV1(0)=DCMPLX(WM(1))
      DV1(1)=DCMPLX(WM(2))
      DV1(2)=DCMPLX(WM(3))
      DV1(3)=DCMPLX(WM(4))
      DP1(0)=DBLE(PWM(0))
      DP1(1)=DBLE(PWM(1))
      DP1(2)=DBLE(PWM(2))
      DP1(3)=DBLE(PWM(3))
      DV2(0)=DCMPLX(W31(1))
      DV2(1)=DCMPLX(W31(2))
      DV2(2)=DCMPLX(W31(3))
      DV2(3)=DCMPLX(W31(4))
      DP2(0)=DBLE(PW31(0))
      DP2(1)=DBLE(PW31(1))
      DP2(2)=DBLE(PW31(2))
      DP2(3)=DBLE(PW31(3))
      DV3(0)=DCMPLX(WP(1))
      DV3(1)=DCMPLX(WP(2))
      DV3(2)=DCMPLX(WP(3))
      DV3(3)=DCMPLX(WP(4))
      DP3(0)=DBLE(PWP(0))
      DP3(1)=DBLE(PWP(1))
      DP3(2)=DBLE(PWP(2))
      DP3(3)=DBLE(PWP(3))
      DV4(0)=DCMPLX(W32(1))
      DV4(1)=DCMPLX(W32(2))
      DV4(2)=DCMPLX(W32(3))
      DV4(3)=DCMPLX(W32(4))
      DP4(0)=DBLE(PW32(0))
      DP4(1)=DBLE(PW32(1))
      DP4(2)=DBLE(PW32(2))
      DP4(3)=DBLE(PW32(3))
C
      V12= DV1(0)*DV2(0)-DV1(1)*DV2(1)-DV1(2)*DV2(2)-DV1(3)*DV2(3)
      V13= DV1(0)*DV3(0)-DV1(1)*DV3(1)-DV1(2)*DV3(2)-DV1(3)*DV3(3)
      V14= DV1(0)*DV4(0)-DV1(1)*DV4(1)-DV1(2)*DV4(2)-DV1(3)*DV4(3)
      V23= DV2(0)*DV3(0)-DV2(1)*DV3(1)-DV2(2)*DV3(2)-DV2(3)*DV3(3)
      V24= DV2(0)*DV4(0)-DV2(1)*DV4(1)-DV2(2)*DV4(2)-DV2(3)*DV4(3)
      V34= DV3(0)*DV4(0)-DV3(1)*DV4(1)-DV3(2)*DV4(2)-DV3(3)*DV4(3)
C
         DVERTX = V14*V23 -V13*V24
C
         VERTEX = DCMPLX( DVERTX ) * (G*G)
C
      RETURN
      END
