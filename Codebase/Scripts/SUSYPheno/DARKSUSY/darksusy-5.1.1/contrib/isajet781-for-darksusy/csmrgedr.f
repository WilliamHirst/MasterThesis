CDECK  ID>, CSMRGEDR.
      SUBROUTINE CSMRGEDR(T,G,F)
!
!Purpose: Contains Standard Model RGE's dG_i/dT = F_i(G) and
!         dY_i/dT = F_i(Y) for Gauge and Yukawas in the DR-bar
!         renormalisation scheme.
!
!This is the double complex version
!
!The common block contains information about n_g (= 2n_f) and SMDR2LP
!which tells CSMRGEDR whether to run the full two loops.
!
!These RGE's come from Martin and Vaughn. Unless specifically
!stated otherwise, all the conventions that follow this
!comment are in Martin and Vaughn notation and are therefore
!different from those in the book by Baer and Tata.
!
!This subroutine takes input in Baer and Tata notation and converts
!to Martin and Vaughn notation. At the end, the conversion is
!run in reverse so that both the inputs and outputs are in Baer/Tata
!notation. Note that in Baer/Tata notation the actual Yukawa
!matrices are the transpose of f_{u,d,e} -- f_{u,d,e} being the inputs
!G(4-12), G(13-21) and G(22-30) respectively. The inputs/outputs are
!in the following order:
!
!     G(  1) = g_1      G(  2) = g_2      G(  3) = g_3
!     G(  4) = FU(1,1)  G(  5) = FU(1,2)  G( 12) = FU(3,3)
!     G( 13) = FD(1,1)  G( 22) = FE(1,1)  G( 30) = FE(3,3)
!     G( 31) = LAMBDA   G( 32) = VEV_SM
!
!
!For more information on the conversion between notations, see
!the subroutine CBK2MVSM in crgeutils.f
!
!Requires crgeutils.f for matrix multiplication, traces and conversions.
!
      IMPLICIT NONE
!
!Couplings and masses from isajet to be used in the thresholds for
!the one loop gauge running
!
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
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/SMRGE/SMRGEMH,SMQSTEP,NU,SMDR2LP
      DOUBLE PRECISION SMRGEMH,SMQSTEP
      INTEGER NU,SMDR2LP
      SAVE/SMRGE/
!
      DOUBLE PRECISION T
      DOUBLE COMPLEX G(32),GPR(32)
!
!     GPR(  1) = g_1      GPR(  2) = g_2      GPR(  3) = g_3
!     GPR(  4) = YU(1,1)  GPR(  5) = YU(1,2)  GPR( 12) = YU(3,3)
!     GPR( 13) = YD(1,1)  GPR( 22) = YE(1,1)  GPR( 30) = YE(3,3)
!     GPR( 31) = LAMBDA - HIGGS QUARTIC COUPLING
!     GPR( 32) = V - HIGGS VEV
!
      DOUBLE COMPLEX F(32)
!
      DOUBLE COMPLEX YU(3,3),YD(3,3),YE(3,3),DY(3,3,3)
!
!     YU        = Y_u   YD        = Y_d   YE        = Y_e
!     DY(1,I,J) = DY_u  DY(2,I,J) = DY_d  DY(3,I,J) = DY_e
!
      DOUBLE PRECISION PI,B1LP(3),B2LP(3,3),C(3,3),ID(3,3)
      DOUBLE COMPLEX Y2,Y4,CHI4,H
      DOUBLE COMPLEX CMATMUL,CTRACE,SUM
!
      DOUBLE COMPLEX DUM1U(3,3),DUM1D(3,3),DUM1E(3,3)
      DOUBLE COMPLEX DUM2U(3,3),DUM2D(3,3),DUM2E(3,3)
      DOUBLE COMPLEX YUDYU(3,3),YDDYD(3,3),YEDYE(3,3)
      DOUBLE COMPLEX YUDYU2(3,3),YDDYD2(3,3),YEDYE2(3,3)
      DOUBLE COMPLEX YUDYU3(3,3),YDDYD3(3,3),YEDYE3(3,3)
      DOUBLE COMPLEX YUDYUYDDYD(3,3),YDDYDYUDYU(3,3)
      DOUBLE COMPLEX TYUDYU,TYDDYD,TYEDYE
      DOUBLE COMPLEX TYUDYU2,TYDDYD2,TYEDYE2
      DOUBLE COMPLEX TYUDYU3,TYDDYD3,TYEDYE3
      DOUBLE COMPLEX TYUDYUYDDYD,TYDDYDYUDYU
      DOUBLE COMPLEX BETA1U(3,3),BETA1D(3,3),BETA1E(3,3)
      DOUBLE COMPLEX BETA2U(3,3),BETA2D(3,3),BETA2E(3,3)
      DOUBLE COMPLEX DUMLUD1(3,3),DUMLUD2(3,3),DUMLUD(3,3),TDUMLUD
      DOUBLE COMPLEX BETALAM1,BETALAM2,BETAVEV1,BETAVEV2
!
      INTEGER I,J,ND,NE,NNU,NSQ,NSU,NSD,NSL,NSE,NSH,NH,NSW,NSG
      DOUBLE PRECISION NG,Q
!
      DATA ID(1,1)/1.D0/,ID(1,2)/0.D0/,ID(1,3)/0.D0/
      DATA ID(2,1)/0.D0/,ID(2,2)/1.D0/,ID(2,3)/0.D0/
      DATA ID(3,1)/0.D0/,ID(3,2)/0.D0/,ID(3,3)/1.D0/
      DATA C(1,1)/1.7D0/,C(1,2)/.5D0/,C(1,3)/1.5D0/
      DATA C(2,1)/1.5D0/,C(2,2)/1.5D0/,C(2,3)/.5D0/
      DATA C(3,1)/2.D0/,C(3,2)/2.D0/,C(3,3)/0.D0/
!
!First convert notations - GPR is input in Martin/Vaughn notation
!
      CALL CBK2MVSM(G,GPR)
!
      DO I=1,31
        F(I)=(0.D0,0.D0)
      END DO
      DO I=1,3
        DO J=1,3
          BETA1U(I,J)=(0.D0,0.D0)
          BETA1D(I,J)=(0.D0,0.D0)
          BETA1E(I,J)=(0.D0,0.D0)
          BETA2U(I,J)=(0.D0,0.D0)
          BETA2D(I,J)=(0.D0,0.D0)
          BETA2E(I,J)=(0.D0,0.D0)
        END DO
      END DO
      BETALAM1=(0.D0,0.D0)
      BETALAM2=(0.D0,0.D0)
      BETAVEV1=(0.D0,0.D0)
      BETAVEV2=(0.D0,0.D0)
!
      PI=4.D0*ATAN(1.D0)
      Q=SMQSTEP
      NG=3.D0
      ND=3
      NE=3
      NNU=3
!
      IF (Q.GT.DBLE(MSS(2))) THEN
        NSQ=3
        NSU=3
        NSD=3
      ELSE
        NSQ=0
        NSU=0
        NSD=0
      END IF
      IF (Q.GT.DBLE(MSS(17))) THEN
        NSL=3
        NSE=3
      ELSE
        NSL=0
        NSE=0
      END IF
!
      IF ((Q-QNSH).GT.ABS(EPS).OR.
     $         (ABS(Q-QNSH).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSH=2
      ELSE
        NSH=0
      END IF
      IF ((Q-QNSG).GT.ABS(EPS).OR.
     $         (ABS(Q-QNSG).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSG=1
      ELSE
        NSG=0
      END IF
      IF ((Q-QNH).GT.ABS(EPS).OR.
     $         (ABS(Q-QNH).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NH=2
      ELSE
        NH=1
      END IF
      IF ((Q-QTHSW).GT.ABS(EPS).OR.
     $         (ABS(Q-QTHSW).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSW=1
      ELSE
        NSW=0
      END IF
!
!Convert input into 3x3 matrices for use in rge's
!
      DO I=1,3
        DO J=1,3
          YU(I,J)=GPR(3+(I-1)*3+J)
          YD(I,J)=GPR(12+(I-1)*3+J)
          YE(I,J)=GPR(21+(I-1)*3+J)
        END DO
      END DO
!
!The separated out contributions are from PRD 49 4882 (1194),
!Castano,Piard,Ramond
!
      B1LP(1)=-(2.D0/5.D0*(17.D0/12.D0*DBLE(NU)+5.D0/12.D0*DBLE(ND)
     $        +5.D0/4.D0*DBLE(NE)+1.D0/4.D0*DBLE(NNU))+1.D0/30.D0
     $        *DBLE(NSQ)+4.D0/15.D0*DBLE(NSU)+1.D0/15.D0*DBLE(NSD)
     $        +1.D0/10.D0*DBLE(NSL)+1.D0/5.D0*DBLE(NSE)
     $        +1.D0/5.D0*DBLE(NSH)+1.D0/10.D0*DBLE(NH))
      B1LP(2)=22.D0/3.D0-(1.D0/2.D0*(DBLE(NU)+DBLE(ND))+1.D0/6.D0
     $        *(DBLE(NE)+DBLE(NNU))+1.D0/2.D0*DBLE(NSQ)
     $        +1.D0/6.D0*DBLE(NSL)+1.D0/3.D0*DBLE(NSH)+1.D0/6.D0
     $        *DBLE(NH)+4.D0/3.D0*DBLE(NSW))
      B1LP(3)=11.D0-(2.D0/3.D0*(DBLE(NU)+DBLE(ND))+1.D0/3.D0
     $        *DBLE(NSQ)+1.D0/6.D0*DBLE(NSU)+1.D0/6.D0*DBLE(NSD)
     $        +2.D0*DBLE(NSG))
!
      DO I=1,3
        DO J=1,3
          YUDYU(I,J)=CMATMUL(1,YU,YU,I,J)
          YDDYD(I,J)=CMATMUL(1,YD,YD,I,J)
          YEDYE(I,J)=CMATMUL(1,YE,YE,I,J)
        END DO
      END DO
      TYUDYU=CTRACE(YUDYU)
      TYDDYD=CTRACE(YDDYD)
      TYEDYE=CTRACE(YEDYE)
      DO I=1,3
        DO J=1,3
          YUDYU2(I,J)=CMATMUL(0,YUDYU,YUDYU,I,J)
          YDDYD2(I,J)=CMATMUL(0,YDDYD,YDDYD,I,J)
          YEDYE2(I,J)=CMATMUL(0,YEDYE,YEDYE,I,J)
        END DO
      END DO
      TYUDYU2=CTRACE(YUDYU2)
      TYDDYD2=CTRACE(YDDYD2)
      TYEDYE2=CTRACE(YEDYE2)
!
      Y2=3.D0*TYUDYU+3.D0*TYDDYD+TYEDYE
      H=3.D0*TYUDYU2+3.D0*TYDDYD2+TYEDYE2
!
      IF(SMDR2LP.EQ.1)THEN
!
        DO I=1,3
          DO J=1,3
            YUDYUYDDYD(I,J)=CMATMUL(0,YUDYU,YDDYD,I,J)
            YDDYDYUDYU(I,J)=CMATMUL(0,YDDYD,YUDYU,I,J)
            DUMLUD1(I,J)=YUDYU(I,J)+YDDYD(I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            DUMLUD2(I,J)=CMATMUL(0,DUMLUD1,YDDYD,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            YUDYU3(I,J)=CMATMUL(0,YUDYU2,YUDYU,I,J)
            YDDYD3(I,J)=CMATMUL(0,YDDYD2,YDDYD,I,J)
            YEDYE3(I,J)=CMATMUL(0,YEDYE2,YEDYE,I,J)
            DUMLUD(I,J)=CMATMUL(0,YUDYU,DUMLUD2,I,J)
          END DO
        END DO
        TYUDYU3=CTRACE(YUDYU3)
        TYDDYD3=CTRACE(YDDYD3)
        TYEDYE3=CTRACE(YEDYE3)
        TYUDYUYDDYD=CTRACE(YUDYUYDDYD)
        TDUMLUD=CTRACE(DUMLUD)
!
!This is the transpose of the Machacek and Vaughn b_kl
!
        B2LP(1,1)=-NG*19.D0/15.D0-9.D0/50.D0
        B2LP(1,2)=-NG*3.D0/5.D0-9.D0/10.D0
        B2LP(1,3)=-NG*44.D0/15.D0
        B2LP(2,1)=-NG*1.D0/5.D0-3.D0/10.D0
        B2LP(2,2)=136.D0/3.D0-NG*49.D0/3.D0-13.D0/6.D0
        B2LP(2,3)=-NG*4.D0
        B2LP(3,1)=-NG*11.D0/30.D0
        B2LP(3,2)=-NG*3.D0/2.D0
        B2LP(3,3)=102.D0-NG*76.D0/3.D0
!
        Y4=(83.D0/40.D0*GPR(1)**2+27.D0/8.D0*GPR(2)**2+28.D0*GPR(3)**2)
     $    *TYUDYU+(-1.D0/40.D0*GPR(1)**2+27.D0/8.D0*GPR(2)**2+28.D0
     $    *GPR(3)**2)*TYDDYD+(93.D0/40.D0*GPR(1)**2+9.D0/8.D0*GPR(2)**2)
     $    *TYEDYE
        CHI4=9.D0/4.D0*(3.D0*TYUDYU2+3.D0*TYDDYD2+TYEDYE2-2.D0/3.D0
     $       *TYUDYUYDDYD)
      ELSE
        DO I=1,3
          DO J=1,3
            YUDYUYDDYD(I,J)=(0.D0,0.D0)
            YDDYDYUDYU(I,J)=(0.D0,0.D0)
            DUMLUD1(I,J)=(0.D0,0.D0)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            DUMLUD2(I,J)=(0.D0,0.D0)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            YUDYU3(I,J)=(0.D0,0.D0)
            YDDYD3(I,J)=(0.D0,0.D0)
            YEDYE3(I,J)=(0.D0,0.D0)
            DUMLUD(I,J)=(0.D0,0.D0)
          END DO
        END DO
        TYUDYU3=(0.D0,0.D0)
        TYDDYD3=(0.D0,0.D0)
        TYEDYE3=(0.D0,0.D0)
        TYUDYUYDDYD=(0.D0,0.D0)
        TDUMLUD=(0.D0,0.D0)
!
!This is the transpose of the Machacek and Vaughn b_kl
!
        DO I=1,3
          DO J=1,3
            B2LP(I,J)=(0.D0,0.D0)
          END DO
        END DO
!
        Y4=(0.D0,0.D0)
        CHI4=(0.D0,0.D0)
      END IF
!
!First, the running of the couplings. This is scheme
!independent to two loops.
!
!
!The B2LP part is different to the MSSM running since the matrix
!is defined as the transpose of the MSSM one.
!
      DO I=1,3
        SUM=0.D0
        IF(SMDR2LP.EQ.1)THEN
          DO J=1,3
            SUM=SUM+B2LP(I,J)*GPR(J)**2
          END DO
        END IF
        F(I)=-GPR(I)**3/16.D0/PI**2*(B1LP(I)+DBLE(SMDR2LP)/16.D0/PI**2
     $       *(SUM+(C(I,1)*TYUDYU+C(I,2)*TYDDYD+C(I,3)*TYEDYE)))
      END DO
!
!Now the Yukawas
!
!
!DUM1U is a dummy matrix for the calculation of BETA1U
!
      DO I=1,3
        DO J=1,3
          DUM1U(I,J)=3.D0/2.D0*(YUDYU(I,J)-YDDYD(I,J))
     $              +Y2*ID(I,J)-(17.D0/20.D0*GPR(1)**2+9.D0/4.D0
     $              *GPR(2)**2+8.D0*GPR(3)**2)*ID(I,J)
          DUM1D(I,J)=3.D0/2.D0*(YDDYD(I,J)-YUDYU(I,J))
     $              +Y2*ID(I,J)-(1.D0/4.D0*GPR(1)**2+9.D0/4.D0*GPR(2)**2
     $              +8.D0*GPR(3)**2)*ID(I,J)
          DUM1E(I,J)=3.D0/2.D0*YEDYE(I,J)+Y2*ID(I,J)-9.D0/4.D0
     $              *(GPR(1)**2+GPR(2)**2)*ID(I,J)
          IF(SMDR2LP.EQ.1)THEN
            DUM2U(I,J)=3.D0/2.D0*YUDYU2(I,J)
     $                -YUDYUYDDYD(I,J)
     $                -1.D0/4.D0*YDDYDYUDYU(I,J)
     $                +11.D0/4.D0*YDDYD2(I,J)
     $                +Y2*(5.D0/4.D0*YDDYD(I,J)-9.D0/4.D0
     $                *YUDYU(I,J))-CHI4*ID(I,J)+3.D0/2.D0*GPR(31)**2
     $                *ID(I,J)-2.D0*GPR(31)*(3.D0*YUDYU(I,J)
     $                +YDDYD(I,J))+(221.D0/80.D0*GPR(1)**2
     $                +117.D0/16.D0*GPR(2)**2+20.D0*GPR(3)**2)
     $                *YUDYU(I,J)-(17.D0/80.D0*GPR(1)**2
     $                -27.D0/16.D0*GPR(2)**2+20.D0*GPR(3)**2)
     $                *YDDYD(I,J)+Y4*ID(I,J)+((7.D0/150.D0
     $                +2.D0/3.D0*NG)*GPR(1)**4-9.D0/20.D0*GPR(1)**2
     $                *GPR(2)**2+19.D0/15.D0*GPR(1)**2*GPR(3)**2
     $                -(101.D0/8.D0-2.D0*NG)*GPR(2)**4+9.D0*GPR(2)**2
     $                *GPR(3)**2-(292.D0/3.D0-16.D0/3.D0*NG)*GPR(3)**4)
     $                *ID(I,J)
            DUM2D(I,J)=3.D0/2.D0*YDDYD2(I,J)
     $                -YDDYDYUDYU(I,J)
     $                -1.D0/4.D0*YUDYUYDDYD(I,J)
     $                +11.D0/4.D0*YUDYU2(I,J)
     $                +Y2*(5.D0/4.D0*YUDYU(I,J)-9.D0/4.D0
     $                *YDDYD(I,J))-CHI4*ID(I,J)+3.D0/2.D0*GPR(31)**2
     $                *ID(I,J)-2.D0*GPR(31)*(3.D0*YDDYD(I,J)
     $                +YUDYU(I,J))+(161.D0/80.D0*GPR(1)**2
     $                +117.D0/16.D0*GPR(2)**2+20.D0*GPR(3)**2)
     $                *YDDYD(I,J)-(77.D0/80.D0*GPR(1)**2
     $                -27.D0/16.D0*GPR(2)**2+20.D0*GPR(3)**2)
     $                *YUDYU(I,J)+Y4*ID(I,J)+(-(37.D0/300.D0
     $                -4.D0/15.D0*NG)*GPR(1)**4-27.D0/20.D0*GPR(1)**2
     $                *GPR(2)**2+31.D0/15.D0*GPR(1)**2*GPR(3)**2
     $                -(101.D0/8.D0-2.D0*NG)*GPR(2)**4+9.D0*GPR(2)**2
     $                *GPR(3)**2-(292.D0/3.D0-16.D0/3.D0*NG)*GPR(3)**4)
     $                *ID(I,J)
            DUM2E(I,J)=3.D0/2.D0*YEDYE2(I,J)
     $                -Y2*9.D0/4.D0*YEDYE(I,J)-CHI4*ID(I,J)
     $                +3.D0/2.D0*GPR(31)**2*ID(I,J)-6.D0*GPR(31)
     $                *YEDYE(I,J)+(441.D0/80.D0*GPR(1)**2+117.D0/16.D0
     $                *GPR(2)**2)*YEDYE(I,J)+Y4*ID(I,J)+((21.D0/100.D0
     $                +8.D0/5.D0*NG)*GPR(1)**4+27.D0/20.D0*GPR(1)**2
     $                *GPR(2)**2-(101.D0/8.D0-2.D0*NG)*GPR(2)**4)
     $                *ID(I,J)
          ELSE
            DUM2U(I,J)=(0.D0,0.D0)
            DUM2D(I,J)=(0.D0,0.D0)
            DUM2E(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
!
!Now calculate the beta functions for the Yukawas
!
      DO I=1,3
        DO J=1,3
          BETA1U(I,J)=CMATMUL(0,YU,DUM1U,I,J)
          BETA1D(I,J)=CMATMUL(0,YD,DUM1D,I,J)
          BETA1E(I,J)=CMATMUL(0,YE,DUM1E,I,J)
          IF(SMDR2LP.EQ.1)THEN
            BETA2U(I,J)=CMATMUL(0,YU,DUM2U,I,J)
            BETA2D(I,J)=CMATMUL(0,YD,DUM2D,I,J)
            BETA2E(I,J)=CMATMUL(0,YE,DUM2E,I,J)
          END IF
!
!Calculate the differentials DY
!
          DY(1,I,J)=1.D0/16.D0/PI**2*BETA1U(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2U(I,J)
          DY(2,I,J)=1.D0/16.D0/PI**2*BETA1D(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2D(I,J)
          DY(3,I,J)=1.D0/16.D0/PI**2*BETA1E(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2E(I,J)
!
!Convert into form readable by RKSTP
!
          F(3+(I-1)*3+J)=DY(1,I,J)
          F(12+(I-1)*3+J)=DY(2,I,J)
          F(21+(I-1)*3+J)=DY(3,I,J)
        END DO
      END DO
!
!Next we have the running of the Higgs Quartic Coupling.
!This is more complicated:
!Programmed here is the MS-bar running. It therefore needs the MS-bar
!gauge and Yukawas.
!The gauge couplings and Yukawas needed to be converted to MS-bar using
!the Martin and Vaughn conversion in hep-ph/9308222.
!The following is after the conversion, so all Yukawas and Gauge
!couplings are still in the DR-bar scheme.
!
      BETALAM1=12*GPR(31)**2-(9.D0/5.D0*GPR(1)**2+9.D0*GPR(2)**2)
     $       *GPR(31)+9.D0/4.D0*(3.D0/25.D0*GPR(1)**4+2.D0/5.D0
     $       *GPR(1)**2*GPR(2)**2+GPR(2)**4)+4.D0*Y2*GPR(31)-4*H
      IF(SMDR2LP.EQ.1)THEN
        BETALAM2=-78.D0*GPR(31)**3+18.D0*(3.D0/5.D0*GPR(1)**2
     $          +3.D0*GPR(2)**2)*GPR(31)**2-((265.D0/8.D0-10*NG)
     $          *GPR(2)**4-117.D0/20.D0*GPR(1)**2*GPR(2)**2-9.D0/25.D0
     $          *(229.D0/24.D0+50.D0/9.D0*NG)*GPR(1)**4)*GPR(31)
     $          +(473.D0/8.D0-8.D0*NG)*GPR(2)**6-3.D0/5.D0
     $          *(121.D0/24.D0+8.D0/3.D0*NG)*GPR(1)**2*GPR(2)**4
     $          -9.D0/25.D0*(239.D0/24.D0+40.D0/9.D0*NG)
     $          *GPR(1)**4*GPR(2)**2-27.D0/125.D0*(59.D0/24.D0
     $          +40.D0/9.D0*NG)*GPR(1)**6+(-14.D0/5.D0*GPR(1)**2
     $          +18.D0*GPR(2)**2-128.D0*GPR(3)**2)*TYUDYU2+(34.D0/5.D0
     $          *GPR(1)**2+18.D0*GPR(2)**2-128.D0*GPR(3)**2)
     $          *TYDDYD2+(-42.D0/5.D0*GPR(1)**2
     $          +6.D0*GPR(2)**2)*TYEDYE2-3.D0/2.D0*GPR(2)**4
     $          *Y2+GPR(31)*((83.D0/10.D0*GPR(1)**2+27.D0/2.D0*GPR(2)**2
     $          +112.D0*GPR(3)**2)*TYUDYU+(-1.D0/10.D0*GPR(1)**2
     $          +27.D0/2.D0*GPR(2)**2+112.D0*GPR(3)**2)*TYDDYD
     $          +(93.D0/10.D0*GPR(1)**2+9.D0/2.D0*GPR(2)**2)*TYEDYE)
     $          +3.D0/5.D0*GPR(1)**2*((-57.D0/10.D0*GPR(1)**2+21.D0
     $          *GPR(2)**2)*TYUDYU+(3.D0/2.D0*GPR(1)**2+9.D0*GPR(2)**2)
     $          *TYDDYD+(-15.D0/2.D0*GPR(1)**2+11.D0*GPR(2)**2)
     $          *TYEDYE)-24.D0*GPR(31)**2*Y2-GPR(31)*H+6.D0*GPR(31)
     $          *TYUDYUYDDYD+20.D0*(3.D0*TYUDYU3
     $          +3.D0*TYDDYD3+TYEDYE3)-12.D0
     $          *TDUMLUD
      END IF
!
      F(31)=1.D0/(16.D0*PI**2)*BETALAM1+1.D0/(16.D0*PI**2)**2*BETALAM2
!
!Caculate the betas for the standard model vev.
!As with lambda this is the MS-bar running with DR-bar inputs except
!v and lambda
!
      BETAVEV1=9.D0/4.D0*(1.D0/5.D0*GPR(1)**2+GPR(2)**2)-Y2
      IF(SMDR2LP.EQ.1)THEN
        BETAVEV2=-3.D0/2.D0*GPR(31)**2-(83.D0/40.D0*GPR(1)**2+27.D0/8.D0
     $           *GPR(2)**2+28.D0*GPR(3)**2)*TYUDYU-(-1.D0/40.D0
     $           *GPR(1)**2+27.D0/8.D0*GPR(2)**2+28.D0*GPR(3)**2)*TYDDYD
     $           -(93.D0/40.D0*GPR(1)**2+9.D0/8.D0*GPR(2)**2)*TYEDYE
     $           +CHI4-27.D0/80.D0*GPR(1)**2*GPR(2)**2-(93.D0/800.D0
     $           +1.D0/2.D0*NG)*GPR(1)**4+(463.D0/32.D0-5.D0/2.D0*NG)
     $           *GPR(2)**4
      END IF
!
      F(32)=GPR(32)*(1.D0/(16.D0*PI**2)*BETAVEV1+1.D0/(16.D0*PI**2)**2
     $      *BETAVEV2)
!
!Convert the output from MV notation to book notation
!
      CALL CBK2MVSM(F,F)
!
      RETURN
      END
