CDECK  ID>, SURG26.
C-----------------------------------------------------------------
      SUBROUTINE SURG26(T,G,F)
C-----------------------------------------------------------------
C
C     Right hand side of full renormalization group equations
C          dG_i/dT = F_i(G)
C     using the thresholds MSS for each mass calculated with the 
C     couplings G0 frozen by SUGFRZ. 
C     Added right neutrino RGE's on 9/24/99
C     Upgrade to 2-loop RGE's for MSSM on 2/11/00
C     Add RGEs for vd and vu on 2/19/03
C     Improved Yukawa thresholds to co-incide with improved A. Box
C     treatment on 8/13/09
      IMPLICIT NONE
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
C          Frozen couplings from RG equations:
C     GSS( 1) = g_1        GSS( 2) = g_2        GSS( 3) = g_3
C     GSS( 4) = y_tau      GSS( 5) = y_b        GSS( 6) = y_t
C     GSS( 7) = M_1        GSS( 8) = M_2        GSS( 9) = M_3
C     GSS(10) = A_tau      GSS(11) = A_b        GSS(12) = A_t
C     GSS(13) = M_hd^2     GSS(14) = M_hu^2     GSS(15) = M_er^2
C     GSS(16) = M_el^2     GSS(17) = M_dnr^2    GSS(18) = M_upr^2
C     GSS(19) = M_upl^2    GSS(20) = M_taur^2   GSS(21) = M_taul^2
C     GSS(22) = M_btr^2    GSS(23) = M_tpr^2    GSS(24) = M_tpl^2
C     GSS(25) = mu         GSS(26) = B          GSS(27) = Y_N
C     GSS(28) = M_nr       GSS(29) = A_n        GSS(30) = vdq
C     GSS(31) = vuq
C          Masses:
C     MSS( 1) = glss     MSS( 2) = upl      MSS( 3) = upr
C     MSS( 4) = dnl      MSS( 5) = dnr      MSS( 6) = stl
C     MSS( 7) = str      MSS( 8) = chl      MSS( 9) = chr
C     MSS(10) = b1       MSS(11) = b2       MSS(12) = t1
C     MSS(13) = t2       MSS(14) = nuel     MSS(15) = numl
C     MSS(16) = nutl     MSS(17) = el-      MSS(18) = er-
C     MSS(19) = mul-     MSS(20) = mur-     MSS(21) = tau1
C     MSS(22) = tau2     MSS(23) = z1ss     MSS(24) = z2ss
C     MSS(25) = z3ss     MSS(26) = z4ss     MSS(27) = w1ss
C     MSS(28) = w2ss     MSS(29) = hl0      MSS(30) = hh0
C     MSS(31) = ha0      MSS(32) = h+
C          Unification:
C     MGUTSS  = M_GUT    GGUTSS  = g_GUT    AGUTSS  = alpha_GUT
      COMMON /SUGMG/ MSS(32),GSS(31),MGUTSS,GGUTSS,AGUTSS,FTGUT,
     $FBGUT,FTAGUT,FNGUT
      REAL MSS,GSS,MGUTSS,GGUTSS,AGUTSS,FTGUT,FBGUT,FTAGUT,FNGUT
      SAVE /SUGMG/
      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,MHLNEG,MHCNEG,
     $IGUTST,ASM3,
     $VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,MHDSMG,MHUSMG,MUMG,BMG,
     $FT2Z1,FB2Z1,FL2Z1
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3,VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,
     $MHDSMG,MHUSMG,MUMG,BMG,FT2Z1,FB2Z1,FL2Z1
      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG,MHLNEG,MHCNEG,IGUTST
      SAVE /SUGPAS/
      REAL*8 T,G(31),F(31)
      REAL*8 FAC,COSB,B3,B2,XTAU,XT,XB,B1,B12,B11,TANB,MT,DMZ,
     $Q,B33,B32,PI,SINB,BETA,B21,B13,B22,B31,B23,XN,SMV,SPMV
      REAL*8 A11,A12,A13,A21,A22,A23,A31,A32,A33,A14,A24
      REAL*8 C11,C12,C13,C21,C22,C23,C31,C32,C33
      REAL*8 D11,D12,D13,D21,D22,D23,D31,D32,D33
      REAL*8 C1,C2,C3,CP1,CP2,CP3,CPP1,CPP2,CPP3,BSY1,BSY2,BSY3
      REAL*8 BY4,BY5,BY6,SIG1,SIG2,SIG3
      REAL*8 G72LP,G82LP,G92LP,FAC2LP,G102lP,G112lP,G122lP
      REAL*8 G132LP,G142LP,G152LP,G162LP,G172LP,G182LP,G192LP,G262LP
      REAL*8 G202LP,G212LP,G222LP,G232LP,G242LP,G252LP,G292LP,G282LP
      REAL*8 DMQ1,DMQ3,DMU1,DMU3,DMD1,DMD3,DML1,DML3,DME1,DME3
      REAL*8 MHSM,VSM,MU2SM,LAMSM,Y2S,Y4S,CHI4S
      INTEGER THHL,THHH,THMU,THQ3,THU3,THD3,THL3,THE3
      INTEGER TH1,TH2,TH3,TH2LP,THTOP
      INTEGER NSL,NSD,NSE,NSH,NSU,NSQ,NH,NSG,NU,NE,ND,NN
C
      DATA ND/3/,NE/3/,NN/3/
      DATA B11/7.96/,B12/5.4/,B13/17.6/,B21/1.8/,B22/25./,B23/24./
      DATA B31/2.2/,B32/9./,B33/14./
      DATA A11/5.2/,A12/2.8/,A13/3.6/,A21/6./,A22/6./,A23/2./
      DATA A31/4./,A32/4./,A33/0./,A14/1.2/,A24/2./
      DATA C11/1.7/,C12/.5/,C13/1.5/,C21/1.5/,C22/1.5/,C23/.5/
      DATA C31/2./,C32/2./,C33/0./
      DATA D11/3.98/,D12/2.7/,D13/8.8/,D21/.9/,D22/5.833/,D23/12./
      DATA D31/1.1/,D32/4.5/,D33/-26./
      DATA C1/.8666667/,C2/3./,C3/5.333333/
      DATA CP1/.4666667/,CP2/3./,CP3/5.333333/
      DATA CPP1/1.8/,CPP2/3./,CPP3/0./
      DATA BSY1/6.6/,BSY2/1./,BSY3/-3./
C
C-----THESE ARE VALID FROM MZ TO MGUT --------------------------------
      PI=4.*DATAN(1.D0)
      TANB=XTANB
      MT=AMT
      Q=MGUT*EXP(T)
      BETA=ATAN(TANB)
      SINB=SIN(BETA)
      COSB=SQRT(1.-SINB**2)
      FAC=2./16.D0/PI**2
      DMZ=91.17D0
C-----HERE, PUT IN A CONSTANT VALUE OF LAMSM AND ASSUME SMALL EVOLUTION
      MHSM=120.D0
      VSM=248.D0
      MU2SM=-MHSM**2/2.D0
      LAMSM=-MU2SM/VSM**2
C-----CALCULATE 1-LOOP THRESHOLD EFFECTS -----------------------------
      TH2LP=1
      NSQ=0
      NSU=0
      NSD=0
      NSL=0
      NSE=0
      THQ3=0
      THU3=0
      THD3=0
      THL3=0
      THE3=0
      DMQ1=DBLE(SQRT(MAX(0.,GSS(19))))
      DMQ3=DBLE(SQRT(MAX(0.,GSS(24))))
      DMU1=DBLE(SQRT(MAX(0.,GSS(18))))
      DMU3=DBLE(SQRT(MAX(0.,GSS(23))))
      DMD1=DBLE(SQRT(MAX(0.,GSS(17))))
      DMD3=DBLE(SQRT(MAX(0.,GSS(22))))
      DML1=DBLE(SQRT(MAX(0.,GSS(16))))
      DML3=DBLE(SQRT(MAX(0.,GSS(21))))
      DME1=DBLE(SQRT(MAX(0.,GSS(15))))
      DME3=DBLE(SQRT(MAX(0.,GSS(20))))
C-----SET THRESHOLDS -------------------------------------------------
      IF (Q.GT.MAX(DMZ,DMQ1)) THEN
        NSQ=NSQ+2
      END IF
      IF (Q.GT.MAX(DMZ,DMQ3)) THEN
        NSQ=NSQ+1
        THQ3=1
      END IF
      IF (Q.GT.MAX(DMZ,DMU1)) THEN
        NSU=NSU+2
      END IF
      IF (Q.GT.MAX(DMZ,DMU3)) THEN
        NSU=NSU+1
        THU3=1
      END IF
      IF (Q.GT.MAX(DMZ,DMD1)) THEN
        NSD=NSD+2
      END IF
      IF (Q.GT.MAX(DMZ,DMD3)) THEN
        NSD=NSD+1
        THD3=1
      END IF
      IF (Q.GT.MAX(DMZ,DME1)) THEN
        NSE=NSE+2
      END IF
      IF (Q.GT.MAX(DMZ,DME3)) THEN
        NSE=NSE+1
        THE3=1
      END IF
      IF (Q.GT.MAX(DMZ,DML1)) THEN
        NSL=NSL+2
      END IF
      IF (Q.GT.MAX(DMZ,DML3)) THEN
        NSL=NSL+1
        THL3=1
      END IF
      IF (Q.GT.DBLE(ABS(MU))) THEN
        NSH=2
      ELSE 
        NSH=0
      END IF
      IF (Q.GT.DBLE(ABS(MSS(1)))) THEN
        NSG=1
      ELSE
        NSG=0
      END IF
      IF (Q.GT.DBLE(MSS(31))) THEN
        NH=2
      ELSE
        NH=1
      END IF
      IF (Q.GT.DBLE(ABS(GSS(7)))) THEN
        TH1=1
      ELSE
        TH1=0
      END IF
      IF (Q.GT.DBLE(ABS(GSS(8)))) THEN
        TH2=1
      ELSE
        TH2=0
      END IF
      IF (Q.GT.MT) THEN
        NU=3
        THTOP=1
      ELSE
        NU=2
        THTOP=0
      END IF
      THHL=1
      THHH=NH/2
      THMU=NSH/2
      TH3=NSG
      B1=2.*(17*NU/12.D0+5*ND/12.D0+5*NE/4.D0+NN/4.D0)/5.D0+
     $ NSQ/30.D0+4*NSU/15.D0+NSD/15.D0+NSL/10.D0+NSE/5.D0+
     $ 1.*NSH/5.D0+1.*NH/10.D0
      B2=-22./3.D0+.5D0*(NU+ND)+1.*(NE+NN)/6.D0+
     $ 1.*NSQ/2.D0+1.*NSL/6.D0+1.*NSH/3.D0+1.*NH/6.D0+4.*TH2/3.D0
      B3=2.*(NU+ND)/3.D0+1.*NSQ/3.D0+1.*NSU/6.D0+1.*NSD/6.D0+
     $2.D0*NSG-11.D0
      IF (Q.GT.DBLE(MSUSY)) THEN
      F(1)=G(1)/16.D0/PI**2*(B1*G(1)**2+TH2LP/16.D0/PI**2*G(1)**2*
     $(B11*G(1)**2+B12*G(2)**2+B13*G(3)**2-A11*G(6)**2-A12*G(5)**2
     $-A13*G(4)**2-A14*G(27)**2))
      F(2)=G(2)/16.D0/PI**2*(B2*G(2)**2+TH2LP/16.D0/PI**2*G(2)**2*
     $(B21*G(1)**2+B22*G(2)**2+B23*G(3)**2-A21*G(6)**2-A22*G(5)**2
     $-A23*G(4)**2-A24*G(27)**2))
      F(3)=G(3)/16.D0/PI**2*(B3*G(3)**2+TH2LP/16.D0/PI**2*G(3)**2*
     $(B31*G(1)**2+B32*G(2)**2+B33*G(3)**2-A31*G(6)**2-A32*G(5)**2
     $-A33*G(4)**2))
      ELSE
      F(1)=G(1)/16.D0/PI**2*(B1*G(1)**2+TH2LP/16.D0/PI**2*G(1)**2*
     $(D11*G(1)**2+D12*G(2)**2+D13*G(3)**2-C11*G(6)**2-C12*G(5)**2
     $-C13*G(4)**2))
      F(2)=G(2)/16.D0/PI**2*(B2*G(2)**2+TH2LP/16.D0/PI**2*G(2)**2*
     $(D21*G(1)**2+D22*G(2)**2+D23*G(3)**2-C21*G(6)**2-C22*G(5)**2
     $-C23*G(4)**2))
      F(3)=G(3)/16.D0/PI**2*(B3*G(3)**2+TH2LP/16.D0/PI**2*G(3)**2*
     $(D31*G(1)**2+D32*G(2)**2+D33*G(3)**2-C31*G(6)**2-C32*G(5)**2
     $-C33*G(4)**2))
      ENDIF
c     The MSSM RGEs for Yukawas include threshold theta fns, so
C     should morph smoothly into SM ones as thresholds are passed
c      IF (Q.LT.MSUSY) THEN
c        F(4)=G(4)/16.D0/PI**2*(5*G(4)**2*COSB**2/2.D0+3*G(6)**2*SINB**2*
c     $   THTOP+3*G(5)**2*COSB**2-9*G(1)**2/4.D0-9*G(2)**2/4.D0
c     $   -SINB**2*(3*G(6)**2*THTOP-3*G(5)**2-G(4)**2))
c        F(5)=G(5)/16.D0/PI**2*(9*G(5)**2*COSB**2/2.D0+
c     $3*G(6)**2*SINB**2*THTOP/2.D0+G(4)**2*COSB**2-G(1)**2/4.D0-
c     $9*G(2)**2/4.D0-8*G(3)**2-SINB**2*(3*G(6)**2*THTOP-
c     $3*G(5)**2-G(4)**2))
c        F(6)=G(6)/16.D0/PI**2*(9*G(6)**2*SINB**2/2.D0*THTOP+
c     $3*G(5)**2*COSB**2/2.D0+G(4)**2*COSB**2-17.*G(1)**2/20.D0
c     $-9*G(2)**2/4.D0-8*G(3)**2+COSB**2*
c     $(3*G(6)**2*THTOP-3*G(5)**2-G(4)**2))
c      ELSE
        BY4=2.5D0*(COSB**2*THHL+SINB**2*THHH)*G(4)**2+
     $.5D0*THMU*(2*THL3+THE3)*G(4)**2+
     $3*SINB**2*(THHL-THHH)*G(6)**2+3*(COSB**2*THHL+SINB**2*THHH)*
     $G(5)**2-
     $3*G(1)**2/5.D0*(15./4.D0-(THL3/4.D0+THE3)*TH1-
     $(THHL/2.D0+THL3-2*THE3)*THMU*TH1)-G(2)**2*(9./4.D0-
     $3*THL3*TH2/4.D0-(3*THHL/2.D0-3*THL3)*THMU*TH2)
        BY4=BY4+G(27)**2
        BY5=4.5D0*(COSB**2*THHL+SINB**2*THHH)*G(5)**2+
     $.5D0*THMU*(2*THQ3+THD3)*G(5)**2+
     $.5D0*(SINB**2*THHL+COSB**2*THHH-
     $4*SINB**2*(THHL-THHH))*G(6)**2+.5D0*THMU*THU3*
     $G(6)**2+3*SINB**2*(THHL-THHH)*G(6)**2+
     $(COSB**2*THHL+SINB**2*THHH)*G(4)**2-
     $3*G(1)**2/5.D0*(5./12.D0-(THQ3/36.D0+THD3/9.D0)*TH1
     $-(THHL/2.D0-THQ3/3.D0-2*THD3/3.D0)*THMU*TH1)
     $-G(2)**2*(9./4.D0-3*THQ3*TH2/4.D0-(3*THHL/2.D0-3*THQ3)*THMU*TH2)-
     $G(3)**2*(8.D0-4*(THQ3+THD3)*TH3/3.D0)
        BY6=4.5D0*(SINB**2*THHL+COSB**2*THHH)*G(6)**2+
     $.5D0*THMU*(2*THQ3+THD3)*G(6)**2+
     $.5D0*(COSB**2*THHL+SINB**2*THHH-
     $4*COSB**2*(THHL-THHH))*G(5)**2+.5D0*THMU*THD3*
     $G(5)**2+COSB**2*(THHL-THHH)*(3*G(5)**2+G(4)**2)-
     $3*G(1)**2/5.D0*(17./12.D0-(THQ3/36.D0+4*THU3/9.D0)*TH1-
     $(THHL/2.D0+THQ3/3.D0-4*THU3/3.D0)*THMU*TH1)
     $-G(2)**2*(9./4.D0-3*THQ3*TH2/4.D0-(3*THHL/2.D0-3*THQ3)*THMU*TH2)-
     $G(3)**2*(8.D0-4*(THQ3+THU3)*TH3/3.D0)
        BY6=BY6+G(27)**2
        F(4)=G(4)/16.D0/PI**2*(BY4+
     $       TH2LP/16.D0/PI**2*((CPP1*BSY1+CPP1**2/2.D0)*G(1)**4+
     $       (CPP2*BSY2+CPP2**2/2.D0)*G(2)**4+(CPP3*BSY3+CPP3**2/2.D0)*
     $       G(3)**4+9*G(1)**2*G(2)**2/5.D0+
     $       G(5)**2*(-.4*G(1)**2+16*G(3)**2)+
     $       G(4)**2*(1.2*G(1)**2+6*G(2)**2)-
     $       (3*G(6)**2*G(5)**2+9*G(5)**4+9*G(5)**2*G(4)**2+
     $       10*G(4)**4)-3*G(27)**4-3*G(27)**2*G(6)**2
     $       -3*G(27)**2*G(4)**2))
        F(5)=G(5)/16.D0/PI**2*(BY5+
     $       TH2LP/16.D0/PI**2*((CP1*BSY1+CP1**2/2.D0)*G(1)**4+
     $       (CP2*BSY2+CP2**2/2.D0)*G(2)**4+(CP3*BSY3+CP3**2/2.D0)*
     $       G(3)**4+G(1)**2*G(2)**2+8*G(1)**2*G(3)**2/9.D0+
     $       8*G(2)**2*G(3)**2+.8*G(6)**2*G(1)**2+
     $       G(5)**2*(.4*G(1)**2+6*G(2)**2+16*G(3)**2)+
     $       1.2*G(4)**2*G(1)**2-
     $       (22*G(5)**4+5*G(6)**2*G(5)**2+3*G(5)**2*G(4)**2+
     $       3*G(4)**4+5*G(6)**4)-G(27)**2*G(6)**2-G(27)**2*G(4)**2))
        F(6)=G(6)/16.D0/PI**2*(BY6+
     $       TH2LP/16.D0/PI**2*((C1*BSY1+C1**2/2.D0)*G(1)**4+
     $       (C2*BSY2+C2**2/2.D0)*G(2)**4+(C3*BSY3+C3**2/2.D0)*
     $       G(3)**4+G(1)**2*G(2)**2+136*G(1)**2*G(3)**2/45.D0+
     $       8*G(2)**2*G(3)**2+G(6)**2*(1.2*G(1)**2+6*G(2)**2+
     $       16*G(3)**2)+.4*G(5)**2*G(1)**2-
     $       (22*G(6)**4+5*G(6)**2*G(5)**2+5*G(5)**4+
     $       G(5)**2*G(4)**2)-3*G(27)**4-3*G(27)**2*G(6)**2-
     $       G(27)**2*G(4)**2))
c      END IF
      FAC2LP=(FAC/2.D0)**2
      G72LP=2*G(1)**2*(B11*G(1)**2*(G(7)+G(7))+
     $B12*G(2)**2*(G(7)+G(8))+B13*G(3)**2*(G(7)+G(9))+
     $A11*G(6)**2*(G(12)-G(7))+A12*G(5)**2*(G(11)-G(7))+
     $A13*G(4)**2*(G(10)-G(7))+A14*G(27)**2*(G(29)-G(7)))
      G82LP=2*G(2)**2*(B21*G(1)**2*(G(7)+G(8))+
     $B22*G(2)**2*(G(8)+G(8))+B23*G(3)**2*(G(8)+G(9))+
     $A21*G(6)**2*(G(12)-G(8))+A22*G(5)**2*(G(11)-G(8))+
     $A23*G(4)**2*(G(10)-G(8))+A24*G(27)**2*(G(29)-G(8)))
      G92LP=2*G(3)**2*(B31*G(1)**2*(G(9)+G(7))+
     $B32*G(2)**2*(G(9)+G(8))+B33*G(3)**2*(G(9)+G(9))+
     $A31*G(6)**2*(G(12)-G(9))+A32*G(5)**2*(G(11)-G(9))+
     $A33*G(4)**2*(G(10)-G(9)))
      F(7)=FAC*B1*G(1)**2*G(7)+TH2LP*FAC2LP*G72LP
      F(8)=FAC*B2*G(2)**2*G(8)+TH2LP*FAC2LP*G82LP
      F(9)=FAC*B3*G(3)**2*G(9)+TH2LP*FAC2LP*G92LP
      XTAU=G(21)+G(20)+G(13)+G(10)**2
      XB=G(24)+G(22)+G(13)+G(11)**2
      XT=G(24)+G(23)+G(14)+G(12)**2
      XN=G(21)+G(28)+G(14)+G(29)**2
      SMV=G(14)-G(13)+(2*G(19)+G(24))-(2*G(16)+G(21))
     $    -2*(2*G(18)+G(23))+(2*G(17)+G(22))+(2*G(15)+G(20))
      SPMV=-2*G(4)**2*G(20)+(6*G(1)**2*(2*G(15)+G(20)))/5.D0+ 
     $G(4)**2*(G(13)+G(21))-(3*(G(1)**2+5*G(2)**2)*
     $(G(13)-G(14)+2*G(16)+G(21)))/10.D0-2*G(5)**2*G(22)+ 
     $(2*(G(1)**2+20*G(3)**2)*(2*G(17)+G(22)))/15.D0+
     $4*G(6)**2*G(23)-(16*(G(1)**2+5*G(3)**2)*(2*G(18)+G(23)))/15.D0+
     $G(5)**2*(3*G(13)-G(24))-G(6)**2*(3*G(14)+G(24))+ 
     $((G(1)**2+45*G(2)**2+80*G(3)**2)*(2*G(19)+G(24)))/30.D0
     $+G(27)**2*(G(21)-G(14))
      SIG1=.2*G(1)**2*(3*(G(14)+G(13))+(2*G(19)+G(24))+
     $3*(2*G(16)+G(21))+8*(2*G(18)+G(23))+2*(2*G(17)+G(22))+
     $6*(2*G(15)+G(20)))
      SIG2=G(2)**2*(G(14)+G(13)+3*(2*G(19)+G(24))+(2*G(16)+G(21)))
      SIG3=G(3)**2*(2*(2*G(19)+G(24))+(2*G(18)+G(23))+
     $(2*G(17)+G(22)))
C
      G102LP=(-54*G(1)**4-(18*G(1)**2*G(2)**2)/5.D0- 
     $(12*G(1)**2*G(4)**2)/5.D0+(4*G(1)**2*G(5)**2)/5.D0)*G(7)+ 
     $((-18*G(1)**2*G(2)**2)/5.D0-30*G(2)**4- 
     $12*G(2)**2*G(4)**2)*G(8)-32*G(3)**2*G(5)**2*G(9)+ 
     $(12*G(1)**2*G(4)**2*G(10))/5.D0+12*G(2)**2*G(4)**2*G(10)-
     $40*G(4)**4*G(10)-36*G(5)**4*G(11)+G(5)**2*
     $((-4*G(1)**2)/5.D0+32*G(3)**2-18*G(4)**2-6*G(6)**2)*G(11)+ 
     $G(5)**2*(-18*G(4)**2*G(10)-6*G(6)**2*G(12))
     $-6*G(27)**2*(G(12)*G(6)**2+G(10)*G(4)**2+G(29)*
     $(G(6)**2+2*G(27)**2+G(4)**2))
      G112LP=((-574*G(1)**4)/45.D0-2*G(1)**2*G(2)**2- 
     $(16*G(1)**2*G(3)**2)/9.D0-(12*G(1)**2*G(4)**2)/5.D0- 
     $(4*G(1)**2*G(5)**2)/5.D0-(8*G(1)**2*G(6)**2)/5.D0)*G(7)+ 
     $(-2*G(1)**2*G(2)**2-30*G(2)**4-16*G(2)**2*G(3)**2-
     $12*G(2)**2*G(5)**2)*G(8)+((-16*G(1)**2*G(3)**2)/9.D0- 
     $16*G(2)**2*G(3)**2+(64*G(3)**4)/9.D0-
     $32*G(3)**2*G(5)**2)*G(9)+(12*G(1)**2*G(4)**2*G(10))/5.D0- 
     $12*G(4)**4*G(10)-88*G(5)**4*G(11)+G(5)**2*((4*G(1)**2)/5.D0+
     $12*G(2)**2+32*G(3)**2-6*G(4)**2-10*G(6)**2)*G(11)+ 
     $(8*G(1)**2*G(6)**2*G(12))/5.D0-20*G(6)**4*G(12)+ 
     $G(5)**2*(-6*G(4)**2*G(10)-10*G(6)**2*G(12))
     $-2*G(27)**2*(G(12)*G(6)**2+G(10)*G(4)**2+G(29)*
     $(G(6)**2+G(4)**2))
      G122LP=((-5486*G(1)**4)/225.D0-2*G(1)**2*G(2)**2- 
     $(272*G(1)**2*G(3)**2)/45.D0-(4*G(1)**2*G(5)**2)/5.D0- 
     $(12*G(1)**2*G(6)**2)/5.D0)*G(7)+(-2*G(1)**2*G(2)**2-
     $30*G(2)**4-16*G(2)**2*G(3)**2-12*G(2)**2*G(6)**2)*G(8)+
     $((-272*G(1)**2*G(3)**2)/45.D0-16*G(2)**2*G(3)**2+
     $(64*G(3)**4)/9.D0-32*G(3)**2*G(6)**2)*G(9)- 
     $2*G(4)**2*G(5)**2*G(10)+(4*G(1)**2*G(5)**2*G(11))/5.D0- 
     $2*G(4)**2*G(5)**2*G(11)-20*G(5)**4*G(11)- 
     $10*G(5)**2*G(6)**2*G(11)+((12*G(1)**2)/5.D0+12*G(2)**2+
     $32*G(3)**2-10*G(5)**2)*G(6)**2*G(12)-88*G(6)**4*G(12)
     $-2*G(27)**2*(3*G(12)*G(6)**2+G(10)*G(4)**2+G(29)*
     $(3*G(6)**2+6*G(27)**2+G(4)**2))
      G132LP=XB*((-4*G(1)**2)/5.D0+32*G(3)**2)*G(5)**2- 
     $(12*G(1)**2*G(4)**2*(-XTAU-2*G(7)**2+2*G(7)*G(10)))/5.D0- 
     $12*G(4)**4*(XTAU+G(10)**2)-36*G(5)**4*(XB+G(11)**2)+ 
     $G(5)**2*((8*G(1)**2*G(7)*(-G(7)+G(11)))/5.D0- 
     $64*G(3)**2*G(9)*(-G(9)+G(11)))-6*G(5)**2*G(6)**2*
     $(XB+XT+2*G(11)*G(12))-1.2*G(1)**2*SPMV+33*G(2)**4*G(8)**2+
     $3.6*G(2)**2*G(1)**2*(G(8)**2+G(7)**2+G(7)*G(8))+
     $24.84*G(1)**4*G(7)**2+3*G(2)**2*SIG2+.6*G(1)**2*SIG1
     $-2*G(4)**2*G(27)**2*(2*G(29)*G(10)+XN+XTAU)
      G142LP=XT*((8*G(1)**2)/5.D0+32*G(3)**2)*G(6)**2- 
     $6*G(5)**2*G(6)**2*(XB+XT+2*G(11)*G(12))- 
     $36*G(6)**4*(XT+G(12)**2)+G(6)**2*((-16*G(1)**2*G(7)*
     $(-G(7)+G(12)))/5.D0-64*G(3)**2*G(9)*(-G(9)+G(12)))+
     $1.2*G(1)**2*SPMV+33*G(2)**4*G(8)**2+
     $3.6*G(2)**2*G(1)**2*(G(8)**2+G(7)**2+G(7)*G(8))+
     $24.84*G(1)**4*G(7)**2+3*G(2)**2*SIG2+.6*G(1)**2*SIG1
     $-12*G(27)**4*(G(29)**2+XN)-2*G(4)**2*G(27)**2*
     $(2*G(29)*G(10)+XN+XTAU)
      G152LP=2.4*G(1)**2*SPMV+112.32*G(1)**4*G(7)**2+2.4*G(1)**2*SIG1
      G162LP=-1.2*G(1)**2*SPMV+33*G(2)**4*G(8)**2+3.6*G(2)**2*G(1)**2*
     $(G(8)**2+G(7)**2+G(8)*G(7))+24.84*G(1)**4*G(7)**2+
     $3*G(2)**2*SIG2+.6*G(1)**2*SIG1
      G172LP=.8*G(1)**2*SPMV-42.66667*G(3)**4*G(9)**2+
     $128*G(3)**2*G(1)**2*(G(9)**2+G(7)**2+G(7)*G(9))/45.D0+808*
     $G(1)**4*G(7)**2/75.D0+16*G(3)**2*SIG3/3.D0+4*G(1)**2*SIG1/15.D0
      G182LP=-1.6*G(1)**2*SPMV-42.66667*G(3)**4*G(9)**2+
     $512*G(3)**2*G(1)**2*(G(9)**2+G(7)**2+G(7)*G(9))/45.D0+
     $3424*G(1)**4*G(7)**2/75.D0+16*G(3)**2*SIG3/3.D0+
     $16*G(1)**2*SIG1/15.D0
      G192LP=.4*G(1)**2*SPMV-128.*G(3)**4*G(9)**2/3.D0+
     $32*G(3)**2*G(2)**2*(G(9)**2+G(8)**2+G(9)*G(8))+
     $32*G(3)**2*G(1)**2*(G(9)**2+G(7)**2+G(9)*G(7))/45.D0+
     $33*G(2)**4*G(8)**2+
     $.4*G(2)**2*G(1)**2*(G(8)**2+G(7)**2+G(8)*G(7))+
     $199*G(1)**4*G(7)**2/75.D0+16*G(3)**2*SIG3/3.D0+
     $3*G(2)**2*SIG2+G(1)**2*SIG1/15.D0
      G202LP=XTAU*((-12*G(1)**2)/5.D0+12*G(2)**2)*G(4)**2- 
     $16*G(4)**4*(XTAU+G(10)**2)+G(4)**2*((24*G(1)**2*G(7)*
     $(-G(7)+G(10)))/5.D0-24*G(2)**2*G(8)*(-G(8)+G(10)))- 
     $12*G(4)**2*G(5)**2*(XB+XTAU+2*G(10)*G(11))+
     $2.4*G(1)**2*SPMV+112.32*G(1)**4*G(7)**2+2.4*G(1)**2*SIG1
     $-4*G(27)**2*G(4)**2*(2*G(29)*G(10)+XN+XTAU)
      G212LP=(-12*G(1)**2*G(4)**2*(-XTAU-2*G(7)**2+2*G(7)*G(10)))/5.D0-
     $12*G(4)**4*(XTAU+G(10)**2)-6*G(4)**2*G(5)**2*
     $(XB+XTAU+2*G(10)*G(11))
     $-1.2*G(1)**2*SPMV+33*G(2)**4*G(8)**2+3.6*G(2)**2*G(1)**2*
     $(G(8)**2+G(7)**2+G(8)*G(7))+24.84*G(1)**4*G(7)**2+
     $3*G(2)**2*SIG2+.6*G(1)**2*SIG1
     $-6*G(27)**2*(2*G(12)*G(29)*G(6)**2+G(6)**2*XT+
     $G(6)**2*XN+2*G(27)**2*(XN+G(29)**2))
      G222LP=XB*((4*G(1)**2)/5.D0+12*G(2)**2)*G(5)**2- 
     $4*G(4)**2*G(5)**2*(XB+XTAU+2*G(10)*G(11))-32*G(5)**4*
     $(XB+G(11)**2)+G(5)**2*((-8*G(1)**2*G(7)*(-G(7)+G(11)))/5.D0- 
     $24*G(2)**2*G(8)*(-G(8)+G(11)))-4*G(5)**2*G(6)**2*(XB+XT+
     $2*G(11)*G(12))+.8*G(1)**2*SPMV-42.66667*G(3)**4*G(9)**2+
     $128*G(3)**2*G(1)**2*(G(9)**2+G(7)**2+G(7)*G(9))/45.D0+808*
     $G(1)**4*G(7)**2/75.D0+16*G(3)**2*SIG3/3.D0+4*G(1)**2*SIG1/15.D0
      G232LP=XT*((-4*G(1)**2)/5.D0+12*G(2)**2)*G(6)**2- 
     $4*G(5)**2*G(6)**2*(XB+XT+2*G(11)*G(12))-32*G(6)**4*
     $(XT+G(12)**2)+G(6)**2*((8*G(1)**2*G(7)*(-G(7)+G(12)))/5.D0- 
     $24*G(2)**2*G(8)*(-G(8)+G(12)))-1.6*G(1)**2*SPMV-42.66667*
     $G(3)**4*G(9)**2+
     $512*G(3)**2*G(1)**2*(G(9)**2+G(7)**2+G(7)*G(9))/45.D0+3424*
     $G(1)**4*G(7)**2/75.D0+16*G(3)**2*SIG3/3.D0+16*G(1)**2*SIG1/15.D0
     $-4*G(6)**2*G(27)**2*(2*G(12)*G(29)+XT+XN)
      G242LP=(-4*G(1)**2*G(5)**2*(-XB-2*G(7)**2+2*G(7)*G(11)))/5.D0- 
     $2*G(4)**2*G(5)**2*(XB+XTAU+2*G(10)*G(11))-20*G(5)**4*
     $(XB+G(11)**2)-(8*G(1)**2*G(6)**2*(-XT-2*G(7)**2+
     $2*G(7)*G(12)))/5.D0-20*G(6)**4*(XT+G(12)**2)
     $+.4*G(1)**2*SPMV-128.*G(3)**4*G(9)**2/3.D0+
     $32*G(3)**2*G(2)**2*(G(9)**2+G(8)**2+G(9)*G(8))+
     $32*G(3)**2*G(1)**2*(G(9)**2+G(7)**2+G(9)*G(7))/45.D0+
     $33*G(2)**4*G(8)**2+
     $.4*G(2)**2*G(1)**2*(G(8)**2+G(7)**2+G(8)*G(7))+
     $199*G(1)**4*G(7)**2/75.D0+16*G(3)**2*SIG3/3.D0+
     $3*G(2)**2*SIG2+G(1)**2*SIG1/15.D0
     $-2*G(6)**2*G(27)**2*(2*G(12)*G(29)+XT+XN)
C     ADD IN MU 2-LOOP TERM SOMETIME...
      G252LP=0.
      G282LP=-16*G(27)**4*(G(29)**2+XN)-12*G(6)**2*G(27)**2*
     $(2*G(12)*G(29)+XT+XN)-4*G(4)**2*G(27)**2*(2*G(29)*G(10)+
     $XN+XTAU)+G(27)**2*G(29)*(-24*G(7)*G(1)**2/5.D0-24*G(2)**2*G(8))
     $+G(27)**2*(2.4*G(1)**2*(2*G(7)**2+XN)+12*G(2)**2*
     $(2*G(8)**2+XN))
      G292LP=-32*G(3)**2*G(6)**2*G(9)+(12*G(1)**2*G(4)**2*G(10))/5.D0
     $-12*G(4)**4*G(10)-6*G(4)**2*G(5)**2*G(10)- 
     $6*G(4)**2*G(5)**2*G(11)-6*G(5)**2*G(6)**2*G(11)+ 
     $(8*G(1)**2*G(6)**2*G(12))/5.D0+32*G(3)**2*G(6)**2*G(12)- 
     $6*G(5)**2*G(6)**2*G(12)-36*G(6)**4*G(12)+ 
     $(-6*G(4)**2*G(10)-18*G(6)**2*G(12))*G(27)**2+
     $((12*G(1)**2)/5.D0+12*G(2)**2-6*G(4)**2-
     $18*G(6)**2)*G(29)*G(27)**2-40*G(29)*G(27)**4+ 
     $G(7)*((-414*G(1)**4)/25.D0-(18*G(1)**2*G(2)**2)/5.D0- 
     $(12*G(1)**2*G(4)**2)/5.D0-(8*G(1)**2*G(6)**2)/5.D0- 
     $(12*G(1)**2*G(27)**2)/5.D0)+G(8)*((-18*G(1)**2*G(2)**2)/5.D0- 
     $30*G(2)**4-12*G(2)**2*G(27)**2)
C
      F(10)=FAC*(9*G(1)**2*G(7)/5.D0+3*G(2)**2*G(8)+3*G(5)**2*G(11)+
     $4*G(4)**2*G(10)+G(27)**2*G(29))+FAC2LP*TH2LP*G102LP
      F(11)=FAC*(7*G(1)**2*G(7)/15.D0+3*G(2)**2*G(8)+16*G(3)**2*
     $G(9)/3.D0+6*G(5)**2*G(11)+G(6)**2*G(12)+G(4)**2*G(10))
     $+FAC2LP*TH2LP*G112LP
      F(12)=FAC*(13*G(1)**2*G(7)/15.D0+3*G(2)**2*G(8)+16*G(3)**2*
     $G(9)/3.D0+G(5)**2*G(11)+6*G(6)**2*G(12)+G(27)**2*G(29))
     $+FAC2LP*TH2LP*G122LP
      F(13)=FAC*(-.6*G(1)**2*G(7)**2-3*G(2)**2*G(8)**2+
     $3*G(5)**2*XB+G(4)**2*XTAU-.3*G(1)**2*SMV)+FAC2LP*TH2LP*G132LP
      F(14)=FAC*(-.6*G(1)**2*G(7)**2-3*G(2)**2*G(8)**2+3*G(6)**2*XT
     $+G(27)**2*XN+.3*G(1)**2*SMV)+FAC2LP*TH2LP*G142LP
      F(15)=FAC*(-2.4*G(1)**2*G(7)**2+.6*G(1)**2*SMV)
     $+FAC2LP*TH2LP*G152LP
      F(16)=FAC*(-.6*G(1)**2*G(7)**2-3*G(2)**2*G(8)**2-.3*G(1)**2*SMV)
     $+FAC2LP*TH2LP*G162LP
      F(17)=FAC*(-4*G(1)**2*G(7)**2/15.D0-16*G(3)**2*G(9)**2/3.D0+
     $.2*G(1)**2*SMV)+FAC2LP*TH2LP*G172LP
      F(18)=FAC*(-16*G(1)**2*G(7)**2/15.D0-16*G(3)**2*G(9)**2/3.D0-
     $.4*G(1)**2*SMV)+FAC2LP*TH2LP*G182LP
      F(19)=FAC*(-G(1)**2*G(7)**2/15.D0-3*G(2)**2*G(8)**2-
     $16*G(3)**2*G(9)**2/3.D0+.1*G(1)**2*SMV)+FAC2LP*TH2LP*G192LP
      F(20)=FAC*(-2.4*G(1)**2*G(7)**2+2*G(4)**2*XTAU+.6*G(1)**2*SMV)+
     $FAC2LP*TH2LP*G202LP
      F(21)=FAC*(-.6*G(1)**2*G(7)**2-3*G(2)**2*G(8)**2+G(4)**2*XTAU
     $+G(27)**2*XN-.3*G(1)**2*SMV)+FAC2LP*TH2LP*G212LP
      F(22)=FAC*(-4*G(1)**2*G(7)**2/15.D0-16*G(3)**2*G(9)**2/3.D0+
     $2*G(5)**2*XB+.2*G(1)**2*SMV)+FAC2LP*TH2LP*G222LP
      F(23)=FAC*(-16*G(1)**2*G(7)**2/15.D0-16*G(3)**2*G(9)**2/3.D0+
     $2*G(6)**2*XT-.4*G(1)**2*SMV)+FAC2LP*TH2LP*G232LP
      F(24)=FAC*(-G(1)**2*G(7)**2/15.D0-3*G(2)**2*G(8)**2-
     $16*G(3)**2*G(9)**2/3.D0+G(6)**2*XT+G(5)**2*XB+.1*G(1)**2*SMV)
     $+FAC2LP*TH2LP*G242LP
      F(25)=FAC*G(25)/2.D0*(-.6*G(1)**2-3*G(2)**2+3*G(6)**2+
     $       3*G(5)**2+G(4)**2+G(27)**2)+FAC2LP*TH2LP*G252LP
C     Version <isa7.81 had sign flip of gauge pieces of F(26)=B
      F(26)=FAC*(.6*G(1)**2*G(7)+3*G(2)**2*G(8)+3*G(6)**2*G(12)+
     $       3*G(5)**2*G(11)+G(4)**2*G(10))
      IF (Q.GT.AMNRMJ) THEN
      F(27)=G(27)/16.D0/PI**2*(3*G(6)**2+G(4)**2+4*G(27)**2-
     $       3*G(2)**2-3*G(1)**2/5.D0+TH2LP/16.D0/PI**2*(-10*G(27)**4
     $       +G(27)**2*(1.2*G(1)**2+3*(2*G(2)**2-3*G(6)**2-
     $       G(4)**2))+207*G(1)**4/50.D0+.2*G(1)**2*(9*G(2)**2+
     $       4*G(6)**2+6*G(4)**2)+.5*(15*G(2)**4-18*G(6)**4-
     $       6*G(4)**4+32*G(3)**2*G(6)**2-6*G(5)**2*G(6)**2-
     $       6*G(5)**2*G(4)**2)))
      F(28)=FAC*2*G(27)**2*XN+FAC2LP*TH2LP*G282LP
      F(29)=FAC*(3*G(1)**2*G(7)/5.D0+3*G(2)**2*G(8)+3*G(6)**2*G(12)+
     $4*G(27)**2*G(29)+G(4)**2*G(10))+FAC2LP*TH2LP*G292LP
      ELSE
      F(27)=0.
      F(28)=0.
      F(29)=0.
      END IF
      F(30)=G(30)*(.75*(G(1)**2/5.D0+G(2)**2)-3*G(5)**2-G(4)**2)/
     $16.D0/PI**2
      F(31)=G(31)*(.75*(G(1)**2/5.D0+G(2)**2)-3*G(6)**2)/16.D0/PI**2
      RETURN
      END
