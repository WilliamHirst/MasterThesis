C
C ======================================================================
C
      SUBROUTINE WWWWXX(WM1,WP1,WM2,WP2,GWWA,GWWZ,ZMASS,ZWIDTH , VERTEX)
C
C this subroutine computes an amplitude of the four-point w-/w+         
C coupling, including the contributions of photon and z exchanges.  the 
C photon propagator is given in feynman gauge and the z propagator is   
C given in unitary gauge.                                               
C                                                                       
C input:                                                                
C       complex wm1(0:3)       : first  flow-out w-                  wm1
C       complex wp1(0:3)       : first  flow-out w+                  wp1
C       complex wm2(0:3)       : second flow-out w-                  wm2
C       complex wp2(0:3)       : second flow-out w+                  wp2
C       real    gwwa           : coupling constant of w and a       gwwa
C       real    gwwz           : coupling constant of w and z       gwwz
C       real    zmass          : mass  of z                             
C       real    zwidth         : width of z                             
C                                                                       
C output:                                                               
C       complex vertex         : amplitude        gamma(wm1,wp1,wm2,wp2)
C
      IMPLICIT NONE
      COMPLEX*16    WM1(6),WP1(6),WM2(6),WP2(6),VERTEX
      COMPLEX*16 DV1(0:3),DV2(0:3),DV3(0:3),DV4(0:3),
     &           J12(0:3),J34(0:3),J14(0:3),J32(0:3),DVERTX,
     &           SV1,SV2,SV3,SV4,TV1,TV2,TV3,TV4,DZS,DZT,
     &           V12,V13,V14,V23,V24,V34,JS12,JS34,JS14,JS32,JS,JT
      REAL*8       PWM1(0:3),PWP1(0:3),PWM2(0:3),PWP2(0:3),
     &           GWWA,GWWZ,ZMASS,ZWIDTH
      REAL*8     Q(0:3),K(0:3),DP1(0:3),DP2(0:3),DP3(0:3),DP4(0:3),
     &           DGWWA2,DGWWZ2,DGW2,DMZ,DWIDTH,S,T,DAS,DAT
C
      REAL*8 RXZERO, RXONE, RXTWO
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0, RXTWO=2.0D0 )
C
      PWM1(0)=DBLE( WM1(5))
      PWM1(1)=DBLE( WM1(6))
      PWM1(2)=DIMAG(WM1(6))
      PWM1(3)=DIMAG(WM1(5))
      PWP1(0)=DBLE( WP1(5))
      PWP1(1)=DBLE( WP1(6))
      PWP1(2)=DIMAG(WP1(6))
      PWP1(3)=DIMAG(WP1(5))
      PWM2(0)=DBLE( WM2(5))
      PWM2(1)=DBLE( WM2(6))
      PWM2(2)=DIMAG(WM2(6))
      PWM2(3)=DIMAG(WM2(5))
      PWP2(0)=DBLE( WP2(5))
      PWP2(1)=DBLE( WP2(6))
      PWP2(2)=DIMAG(WP2(6))
      PWP2(3)=DIMAG(WP2(5))
C
      DV1(0)=DCMPLX(WM1(1))
      DV1(1)=DCMPLX(WM1(2))
      DV1(2)=DCMPLX(WM1(3))
      DV1(3)=DCMPLX(WM1(4))
      DP1(0)=DBLE(PWM1(0))
      DP1(1)=DBLE(PWM1(1))
      DP1(2)=DBLE(PWM1(2))
      DP1(3)=DBLE(PWM1(3))
      DV2(0)=DCMPLX(WP1(1))
      DV2(1)=DCMPLX(WP1(2))
      DV2(2)=DCMPLX(WP1(3))
      DV2(3)=DCMPLX(WP1(4))
      DP2(0)=DBLE(PWP1(0))
      DP2(1)=DBLE(PWP1(1))
      DP2(2)=DBLE(PWP1(2))
      DP2(3)=DBLE(PWP1(3))
      DV3(0)=DCMPLX(WM2(1))
      DV3(1)=DCMPLX(WM2(2))
      DV3(2)=DCMPLX(WM2(3))
      DV3(3)=DCMPLX(WM2(4))
      DP3(0)=DBLE(PWM2(0))
      DP3(1)=DBLE(PWM2(1))
      DP3(2)=DBLE(PWM2(2))
      DP3(3)=DBLE(PWM2(3))
      DV4(0)=DCMPLX(WP2(1))
      DV4(1)=DCMPLX(WP2(2))
      DV4(2)=DCMPLX(WP2(3))
      DV4(3)=DCMPLX(WP2(4))
      DP4(0)=DBLE(PWP2(0))
      DP4(1)=DBLE(PWP2(1))
      DP4(2)=DBLE(PWP2(2))
      DP4(3)=DBLE(PWP2(3))
      DGWWA2=DBLE(GWWA)**2
      DGWWZ2=DBLE(GWWZ)**2
      DGW2  =DGWWA2+DGWWZ2
      DMZ   =DBLE(ZMASS)
      DWIDTH=DBLE(ZWIDTH)
C
      V12= DV1(0)*DV2(0)-DV1(1)*DV2(1)-DV1(2)*DV2(2)-DV1(3)*DV2(3)
      V13= DV1(0)*DV3(0)-DV1(1)*DV3(1)-DV1(2)*DV3(2)-DV1(3)*DV3(3)
      V14= DV1(0)*DV4(0)-DV1(1)*DV4(1)-DV1(2)*DV4(2)-DV1(3)*DV4(3)
      V23= DV2(0)*DV3(0)-DV2(1)*DV3(1)-DV2(2)*DV3(2)-DV2(3)*DV3(3)
      V24= DV2(0)*DV4(0)-DV2(1)*DV4(1)-DV2(2)*DV4(2)-DV2(3)*DV4(3)
      V34= DV3(0)*DV4(0)-DV3(1)*DV4(1)-DV3(2)*DV4(2)-DV3(3)*DV4(3)
C
      Q(0)=DP1(0)+DP2(0)
      Q(1)=DP1(1)+DP2(1)
      Q(2)=DP1(2)+DP2(2)
      Q(3)=DP1(3)+DP2(3)
      K(0)=DP1(0)+DP4(0)
      K(1)=DP1(1)+DP4(1)
      K(2)=DP1(2)+DP4(2)
      K(3)=DP1(3)+DP4(3)
C
      S=Q(0)**2-Q(1)**2-Q(2)**2-Q(3)**2
      T=K(0)**2-K(1)**2-K(2)**2-K(3)**2
C
      DAS=-RXONE/S
      DAT=-RXONE/T
      DZS=-RXONE/DCMPLX( S-DMZ**2 , DMAX1(DSIGN(DMZ*DWIDTH,S),RXZERO) )
      DZT=-RXONE/DCMPLX( T-DMZ**2 , DMAX1(DSIGN(DMZ*DWIDTH,T),RXZERO) )
C
      SV1= (DP2(0)+Q(0))*DV1(0) -(DP2(1)+Q(1))*DV1(1)
     &    -(DP2(2)+Q(2))*DV1(2) -(DP2(3)+Q(3))*DV1(3)
      SV2=-(DP1(0)+Q(0))*DV2(0) +(DP1(1)+Q(1))*DV2(1)
     &    +(DP1(2)+Q(2))*DV2(2) +(DP1(3)+Q(3))*DV2(3)
      SV3= (DP4(0)-Q(0))*DV3(0) -(DP4(1)-Q(1))*DV3(1)
     &    -(DP4(2)-Q(2))*DV3(2) -(DP4(3)-Q(3))*DV3(3)
      SV4=-(DP3(0)-Q(0))*DV4(0) +(DP3(1)-Q(1))*DV4(1)
     &    +(DP3(2)-Q(2))*DV4(2) +(DP3(3)-Q(3))*DV4(3)
C
      TV1= (DP4(0)+K(0))*DV1(0) -(DP4(1)+K(1))*DV1(1)
     &    -(DP4(2)+K(2))*DV1(2) -(DP4(3)+K(3))*DV1(3)
      TV2=-(DP3(0)-K(0))*DV2(0) +(DP3(1)-K(1))*DV2(1)
     &    +(DP3(2)-K(2))*DV2(2) +(DP3(3)-K(3))*DV2(3)
      TV3= (DP2(0)-K(0))*DV3(0) -(DP2(1)-K(1))*DV3(1)
     &    -(DP2(2)-K(2))*DV3(2) -(DP2(3)-K(3))*DV3(3)
      TV4=-(DP1(0)+K(0))*DV4(0) +(DP1(1)+K(1))*DV4(1)
     &    +(DP1(2)+K(2))*DV4(2) +(DP1(3)+K(3))*DV4(3)
C
      J12(0)=(DP1(0)-DP2(0))*V12 +SV1*DV2(0) +SV2*DV1(0)
      J12(1)=(DP1(1)-DP2(1))*V12 +SV1*DV2(1) +SV2*DV1(1)
      J12(2)=(DP1(2)-DP2(2))*V12 +SV1*DV2(2) +SV2*DV1(2)
      J12(3)=(DP1(3)-DP2(3))*V12 +SV1*DV2(3) +SV2*DV1(3)
      J34(0)=(DP3(0)-DP4(0))*V34 +SV3*DV4(0) +SV4*DV3(0)
      J34(1)=(DP3(1)-DP4(1))*V34 +SV3*DV4(1) +SV4*DV3(1)
      J34(2)=(DP3(2)-DP4(2))*V34 +SV3*DV4(2) +SV4*DV3(2)
      J34(3)=(DP3(3)-DP4(3))*V34 +SV3*DV4(3) +SV4*DV3(3)
C
      J14(0)=(DP1(0)-DP4(0))*V14 +TV1*DV4(0) +TV4*DV1(0)
      J14(1)=(DP1(1)-DP4(1))*V14 +TV1*DV4(1) +TV4*DV1(1)
      J14(2)=(DP1(2)-DP4(2))*V14 +TV1*DV4(2) +TV4*DV1(2)
      J14(3)=(DP1(3)-DP4(3))*V14 +TV1*DV4(3) +TV4*DV1(3)
      J32(0)=(DP3(0)-DP2(0))*V23 +TV3*DV2(0) +TV2*DV3(0)
      J32(1)=(DP3(1)-DP2(1))*V23 +TV3*DV2(1) +TV2*DV3(1)
      J32(2)=(DP3(2)-DP2(2))*V23 +TV3*DV2(2) +TV2*DV3(2)
      J32(3)=(DP3(3)-DP2(3))*V23 +TV3*DV2(3) +TV2*DV3(3)
C
      JS12=Q(0)*J12(0)-Q(1)*J12(1)-Q(2)*J12(2)-Q(3)*J12(3)
      JS34=Q(0)*J34(0)-Q(1)*J34(1)-Q(2)*J34(2)-Q(3)*J34(3)
      JS14=K(0)*J14(0)-K(1)*J14(1)-K(2)*J14(2)-K(3)*J14(3)
      JS32=K(0)*J32(0)-K(1)*J32(1)-K(2)*J32(2)-K(3)*J32(3)
C
      JS=J12(0)*J34(0)-J12(1)*J34(1)-J12(2)*J34(2)-J12(3)*J34(3)
      JT=J14(0)*J32(0)-J14(1)*J32(1)-J14(2)*J32(2)-J14(3)*J32(3)
C
      DVERTX = (V12*V34 +V14*V23 -RXTWO*V13*V24)*DGW2
C     &        +(dzs*dgwwz2+das*dgwwa2)*js -dzs*dgwwz2*js12*js34/dmz**2
C     &        +(dzt*dgwwz2+dat*dgwwa2)*jt -dzt*dgwwz2*js14*js32/dmz**2
C
      VERTEX = -DCMPLX( DVERTX )
C
      RETURN
      END
