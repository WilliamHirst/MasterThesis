C
C ----------------------------------------------------------------------
C
      SUBROUTINE JWWWXX(W1,W2,W3,GWWA,GWWZ,ZMASS,ZWIDTH,WMASS,WWIDTH ,
     &                  JWWW)
C
C this subroutine computes an off-shell w+/w- current from the four-    
C point gauge boson coupling, including the contributions of photon and 
C z exchanges.  the vector propagators for the output w and the internal
C z bosons are given in unitary gauge, and that of the internal photon  
C is given in feynman gauge.                                            
C                                                                       
C input:                                                                
C       complex w1(6)          : first  vector                        w1
C       complex w2(6)          : second vector                        w2
C       complex w3(6)          : third  vector                        w3
C       real    gwwa           : coupling constant of w and a       gwwa
C       real    gwwz           : coupling constant of w and z       gwwz
C       real    zmass          : mass  of internal z                    
C       real    zwidth         : width of internal z                    
C       real    wmass          : mass  of output w                      
C       real    wwidth         : width of output w                      
C                                                                       
C the possible sets of the inputs are as follows:                       
C   ------------------------------------------------------------------- 
C   |  w1  |  w2  |  w3  |gwwa|gwwz|zmass|zwidth|wmass|wwidth || jwww | 
C   ------------------------------------------------------------------- 
C   |  w-  |  w+  |  w-  |gwwa|gwwz|zmass|zwidth|wmass|wwidth ||  w+  | 
C   |  w+  |  w-  |  w+  |gwwa|gwwz|zmass|zwidth|wmass|wwidth ||  w-  | 
C   ------------------------------------------------------------------- 
C where all the bosons are defined by the flowing-out quantum number.   
C                                                                       
C output:                                                               
C       complex jwww(6)        : w current             j^mu(w':w1,w2,w3)
C
      IMPLICIT NONE
      COMPLEX*16  W1(6),W2(6),W3(6),JWWW(6)
      COMPLEX*16 DW1(0:3),DW2(0:3),DW3(0:3),
     &           JJ(0:3),JS(0:3),JT(0:3),J4(0:3),
     &           JT12(0:3),JT32(0:3),J12(0:3),J32(0:3),
     &           DZS,DZT,DW,W12,W32,W13,P1W2,P2W1,P3W2,P2W3,
     &           JK12,JK32,JSW3,JTW1,P3JS,KSW3,P1JT,KTW1,JQ
      REAL*8     GWWA,GWWZ,ZMASS,ZWIDTH,WMASS,WWIDTH
      REAL*8     P1(0:3),P2(0:3),P3(0:3),Q(0:3),KS(0:3),KT(0:3),
     &           DGWWA2,DGWWZ2,DGW2,DMZ,DWZ,DMW,DWW,MZ2,MW2,Q2,KS2,KT2,
     &           DAS,DAT
C
      JWWW(5) = W1(5)+W2(5)+W3(5)
      JWWW(6) = W1(6)+W2(6)+W3(6)
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
      KS(0)=P1(0)+P2(0)
      KS(1)=P1(1)+P2(1)
      KS(2)=P1(2)+P2(2)
      KS(3)=P1(3)+P2(3)
      KT(0)=P2(0)+P3(0)
      KT(1)=P2(1)+P3(1)
      KT(2)=P2(2)+P3(2)
      KT(3)=P2(3)+P3(3)
      Q2 =Q(0)**2 -(Q(1)**2 +Q(2)**2 +Q(3)**2)
      KS2=KS(0)**2-(KS(1)**2+KS(2)**2+KS(3)**2)
      KT2=KT(0)**2-(KT(1)**2+KT(2)**2+KT(3)**2)
      DGWWA2=DBLE(GWWA)**2
      DGWWZ2=DBLE(GWWZ)**2
      DGW2  =DGWWA2+DGWWZ2
      DMZ=DBLE(ZMASS)
      DWZ=DBLE(ZWIDTH)
      DMW=DBLE(WMASS)
      DWW=DBLE(WWIDTH)
      MZ2=DMZ**2
      MW2=DMW**2
C
      DAS=-DGWWA2/KS2
      DAT=-DGWWA2/KT2
      DZS=-DGWWZ2/DCMPLX( KS2-MZ2 , DMAX1(DSIGN(DMZ*DWZ,KS2),0.D0) )
      DZT=-DGWWZ2/DCMPLX( KT2-MZ2 , DMAX1(DSIGN(DMZ*DWZ,KT2),0.D0) )
      DW =-1.0D0/DCMPLX( Q2 -MW2 , DMAX1(DSIGN(DMW*DWW,Q2 ),0.D0) )
C  for the running width, use below instead of the above dw.
C      dw =-1.0d0/dcmplx( q2 -mw2 , dmax1(dww*q2/dmw,0.d0) )
C
      W12=DW1(0)*DW2(0)-DW1(1)*DW2(1)-DW1(2)*DW2(2)-DW1(3)*DW2(3)
      W32=DW3(0)*DW2(0)-DW3(1)*DW2(1)-DW3(2)*DW2(2)-DW3(3)*DW2(3)
C
      P1W2= (P1(0)+KS(0))*DW2(0)-(P1(1)+KS(1))*DW2(1)
     &     -(P1(2)+KS(2))*DW2(2)-(P1(3)+KS(3))*DW2(3)
      P2W1= (P2(0)+KS(0))*DW1(0)-(P2(1)+KS(1))*DW1(1)
     &     -(P2(2)+KS(2))*DW1(2)-(P2(3)+KS(3))*DW1(3)
      P3W2= (P3(0)+KT(0))*DW2(0)-(P3(1)+KT(1))*DW2(1)
     &     -(P3(2)+KT(2))*DW2(2)-(P3(3)+KT(3))*DW2(3)
      P2W3= (P2(0)+KT(0))*DW3(0)-(P2(1)+KT(1))*DW3(1)
     &     -(P2(2)+KT(2))*DW3(2)-(P2(3)+KT(3))*DW3(3)
C
      JT12(0)= (P1(0)-P2(0))*W12 + P2W1*DW2(0) - P1W2*DW1(0)
      JT12(1)= (P1(1)-P2(1))*W12 + P2W1*DW2(1) - P1W2*DW1(1)
      JT12(2)= (P1(2)-P2(2))*W12 + P2W1*DW2(2) - P1W2*DW1(2)
      JT12(3)= (P1(3)-P2(3))*W12 + P2W1*DW2(3) - P1W2*DW1(3)
      JT32(0)= (P3(0)-P2(0))*W32 + P2W3*DW2(0) - P3W2*DW3(0)
      JT32(1)= (P3(1)-P2(1))*W32 + P2W3*DW2(1) - P3W2*DW3(1)
      JT32(2)= (P3(2)-P2(2))*W32 + P2W3*DW2(2) - P3W2*DW3(2)
      JT32(3)= (P3(3)-P2(3))*W32 + P2W3*DW2(3) - P3W2*DW3(3)
C
      JK12=(JT12(0)*KS(0)-JT12(1)*KS(1)-JT12(2)*KS(2)-JT12(3)*KS(3))/MZ2
      JK32=(JT32(0)*KT(0)-JT32(1)*KT(1)-JT32(2)*KT(2)-JT32(3)*KT(3))/MZ2
C
      J12(0)=JT12(0)*(DAS+DZS)-KS(0)*JK12*DZS
      J12(1)=JT12(1)*(DAS+DZS)-KS(1)*JK12*DZS
      J12(2)=JT12(2)*(DAS+DZS)-KS(2)*JK12*DZS
      J12(3)=JT12(3)*(DAS+DZS)-KS(3)*JK12*DZS
      J32(0)=JT32(0)*(DAT+DZT)-KT(0)*JK32*DZT
      J32(1)=JT32(1)*(DAT+DZT)-KT(1)*JK32*DZT
      J32(2)=JT32(2)*(DAT+DZT)-KT(2)*JK32*DZT
      J32(3)=JT32(3)*(DAT+DZT)-KT(3)*JK32*DZT
C
      JSW3=J12(0)*DW3(0)-J12(1)*DW3(1)-J12(2)*DW3(2)-J12(3)*DW3(3)
      JTW1=J32(0)*DW1(0)-J32(1)*DW1(1)-J32(2)*DW1(2)-J32(3)*DW1(3)
C
      P3JS= (P3(0)-Q(0))*J12(0)-(P3(1)-Q(1))*J12(1)
     &     -(P3(2)-Q(2))*J12(2)-(P3(3)-Q(3))*J12(3)
      KSW3= (KS(0)-Q(0))*DW3(0)-(KS(1)-Q(1))*DW3(1)
     &     -(KS(2)-Q(2))*DW3(2)-(KS(3)-Q(3))*DW3(3)
      P1JT= (P1(0)-Q(0))*J32(0)-(P1(1)-Q(1))*J32(1)
     &     -(P1(2)-Q(2))*J32(2)-(P1(3)-Q(3))*J32(3)
      KTW1= (KT(0)-Q(0))*DW1(0)-(KT(1)-Q(1))*DW1(1)
     &     -(KT(2)-Q(2))*DW1(2)-(KT(3)-Q(3))*DW1(3)
C
      JS(0)= (KS(0)-P3(0))*JSW3 + P3JS*DW3(0) - KSW3*J12(0)
      JS(1)= (KS(1)-P3(1))*JSW3 + P3JS*DW3(1) - KSW3*J12(1)
      JS(2)= (KS(2)-P3(2))*JSW3 + P3JS*DW3(2) - KSW3*J12(2)
      JS(3)= (KS(3)-P3(3))*JSW3 + P3JS*DW3(3) - KSW3*J12(3)
      JT(0)= (KT(0)-P1(0))*JTW1 + P1JT*DW1(0) - KTW1*J32(0)
      JT(1)= (KT(1)-P1(1))*JTW1 + P1JT*DW1(1) - KTW1*J32(1)
      JT(2)= (KT(2)-P1(2))*JTW1 + P1JT*DW1(2) - KTW1*J32(2)
      JT(3)= (KT(3)-P1(3))*JTW1 + P1JT*DW1(3) - KTW1*J32(3)
C
      W13=DW1(0)*DW3(0)-DW1(1)*DW3(1)-DW1(2)*DW3(2)-DW1(3)*DW3(3)
C
      J4(0)=DGW2*( DW1(0)*W32 + DW3(0)*W12 - 2.D0*DW2(0)*W13 )
      J4(1)=DGW2*( DW1(1)*W32 + DW3(1)*W12 - 2.D0*DW2(1)*W13 )
      J4(2)=DGW2*( DW1(2)*W32 + DW3(2)*W12 - 2.D0*DW2(2)*W13 )
      J4(3)=DGW2*( DW1(3)*W32 + DW3(3)*W12 - 2.D0*DW2(3)*W13 )
C
      JJ(0)=J4(0)
      JJ(1)=J4(1)
      JJ(2)=J4(2)
      JJ(3)=J4(3)
C
      JQ=(JJ(0)*Q(0)-JJ(1)*Q(1)-JJ(2)*Q(2)-JJ(3)*Q(3))/MW2
C
      JWWW(1) = DCMPLX( (JJ(0)-JQ*Q(0))*DW )
      JWWW(2) = DCMPLX( (JJ(1)-JQ*Q(1))*DW )
      JWWW(3) = DCMPLX( (JJ(2)-JQ*Q(2))*DW )
      JWWW(4) = DCMPLX( (JJ(3)-JQ*Q(3))*DW )
C
      RETURN
      END
