C
C ----------------------------------------------------------------------
C
      SUBROUTINE W3W3XX(WM,W31,WP,W32,G31,G32,WMASS,WWIDTH , VERTEX)
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
C       real    g31            : coupling of w31 with w-/w+             
C       real    g32            : coupling of w32 with w-/w+             
C                                                  (see the table below)
C       real    wmass          : mass  of w                             
C       real    wwidth         : width of w                             
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
      COMPLEX*16 DV1(0:3),DV2(0:3),DV3(0:3),DV4(0:3),DVERTX,
     &           V12,V13,V14,V23,V24,V34
      REAL*8     G31,G32,WMASS,WWIDTH
C
      REAL*8 RXZERO, RXONE
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0 )
C
      DV1(0)=DCMPLX(WM(1))
      DV1(1)=DCMPLX(WM(2))
      DV1(2)=DCMPLX(WM(3))
      DV1(3)=DCMPLX(WM(4))
      DV2(0)=DCMPLX(W31(1))
      DV2(1)=DCMPLX(W31(2))
      DV2(2)=DCMPLX(W31(3))
      DV2(3)=DCMPLX(W31(4))
      DV3(0)=DCMPLX(WP(1))
      DV3(1)=DCMPLX(WP(2))
      DV3(2)=DCMPLX(WP(3))
      DV3(3)=DCMPLX(WP(4))
      DV4(0)=DCMPLX(W32(1))
      DV4(1)=DCMPLX(W32(2))
      DV4(2)=DCMPLX(W32(3))
      DV4(3)=DCMPLX(W32(4))
C
      IF ( DBLE(WMASS) .NE. RXZERO ) THEN
C         dm2inv = r_one / dmw2
C
         V12= DV1(0)*DV2(0)-DV1(1)*DV2(1)-DV1(2)*DV2(2)-DV1(3)*DV2(3)
         V13= DV1(0)*DV3(0)-DV1(1)*DV3(1)-DV1(2)*DV3(2)-DV1(3)*DV3(3)
         V14= DV1(0)*DV4(0)-DV1(1)*DV4(1)-DV1(2)*DV4(2)-DV1(3)*DV4(3)
         V23= DV2(0)*DV3(0)-DV2(1)*DV3(1)-DV2(2)*DV3(2)-DV2(3)*DV3(3)
         V24= DV2(0)*DV4(0)-DV2(1)*DV4(1)-DV2(2)*DV4(2)-DV2(3)*DV4(3)
         V34= DV3(0)*DV4(0)-DV3(1)*DV4(1)-DV3(2)*DV4(2)-DV3(3)*DV4(3)
C
         DVERTX = V12*V34 +V14*V23 -2.D0*V13*V24
C
         VERTEX = DCMPLX( DVERTX ) * (G31*G32)
C
      ELSE
         V12= DV1(0)*DV2(0)-DV1(1)*DV2(1)-DV1(2)*DV2(2)-DV1(3)*DV2(3)
         V13= DV1(0)*DV3(0)-DV1(1)*DV3(1)-DV1(2)*DV3(2)-DV1(3)*DV3(3)
         V14= DV1(0)*DV4(0)-DV1(1)*DV4(1)-DV1(2)*DV4(2)-DV1(3)*DV4(3)
         V23= DV2(0)*DV3(0)-DV2(1)*DV3(1)-DV2(2)*DV3(2)-DV2(3)*DV3(3)
         V24= DV2(0)*DV4(0)-DV2(1)*DV4(1)-DV2(2)*DV4(2)-DV2(3)*DV4(3)
         V34= DV3(0)*DV4(0)-DV3(1)*DV4(1)-DV3(2)*DV4(2)-DV3(3)*DV4(3)
C
         DVERTX = V14*V23 -V13*V24
C
         VERTEX = DCMPLX( DVERTX ) * (G31*G32)
      END IF
C
      RETURN
      END
