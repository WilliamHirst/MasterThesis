CDECK  ID>, SUGFRZ.
C------------------------------------------------------------------
      SUBROUTINE SUGFRZ(Q,G,G0,IG)
C------------------------------------------------------------------
C
C     Freeze out final soft breaking parameters
C
      IMPLICIT NONE
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
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
C          SUSY parameters
C          AMGLSS               = gluino mass
C          AMULSS               = up-left squark mass
C          AMELSS               = left-selectron mass
C          AMERSS               = right-slepton mass
C          AMNiSS               = sneutrino mass for generation i
C          TWOM1                = Higgsino mass = - mu
C          RV2V1                = ratio v2/v1 of vev's
C          AMTLSS,AMTRSS        = left,right stop masses
C          AMT1SS,AMT2SS        = light,heavy stop masses
C          AMBLSS,AMBRSS        = left,right sbottom masses
C          AMB1SS,AMB2SS        = light,heavy sbottom masses
C          AMLLSS,AMLRSS        = left,right stau masses
C          AML1SS,AML2SS        = light,heavy stau masses
C          AMZiSS               = signed mass of Zi
C          ZMIXSS               = Zi mixing matrix
C          AMWiSS               = signed Wi mass
C          GAMMAL,GAMMAR        = Wi left, right mixing angles
C          AMHL,AMHH,AMHA       = neutral Higgs h0, H0, A0 masses
C          AMHC                 = charged Higgs H+ mass
C          ALFAH                = Higgs mixing angle
C          AAT                  = stop trilinear term
C          THETAT               = stop mixing angle
C          AAB                  = sbottom trilinear term
C          THETAB               = sbottom mixing angle
C          AAL                  = stau trilinear term
C          THETAL               = stau mixing angle
C          AMGVSS               = gravitino mass
C          MTQ                  = top mass at MSUSY
C          MBQ                  = bottom mass at MSUSY
C          MLQ                  = tau mass at MSUSY
C          FBMA                 = b-Yukawa at mA scale
C          VUQ                  = Hu vev at MSUSY
C          VDQ                  = Hd vev at MSUSY
C          SGNM3                = sign of gaugino mass M3
      COMMON/SSPAR/GORGE,AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS(4,4)
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,
     $VUQ,VDQ,SGNM3
      REAL AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,VUQ,VDQ,SGNM3
      LOGICAL GORGE
      REAL AMZISS(4)
      EQUIVALENCE (AMZISS(1),AMZ1SS)
      SAVE /SSPAR/
      INTEGER IG(31)
      REAL Q,MT,PI,TMZ1
      REAL*8 G(31),DPI
      REAL G0(31)
      INTEGER I
C
      MT=AMT
      PI=4.*ATAN(1.)
      DPI=4.D0*DATAN(1.D0)
      TMZ1=2*ABS(AMZ1SS)
      DO 200 I=1,3
        G0(I)=SNGL(G(I))
200   CONTINUE
C          Freeze out Yukawa couplings and A-terms at HIGFRZ
      DO 205 I=4,6
        IF (Q.LT.HIGFRZ.AND.IG(I).EQ.0) THEN
          G0(I)=SNGL(G(I))
          IG(I)=1
          G0(I+6)=SNGL(G(I+6))
          IG(I+6)=1
        ELSE IF (IG(I).EQ.0) THEN
          G0(I)=SNGL(G(I))
          G0(I+6)=SNGL(G(I+6))
        END IF
205   CONTINUE
C         Extract b Yukawa at mA for use in Higgs decay rates
      IF (Q.GT.AMHA) THEN
        FBMA=SNGL(G(5))
      END IF
      IF (Q.GT.TMZ1) THEN
        FT2Z1=SNGL(G(6))
        FB2Z1=SNGL(G(5))
        FL2Z1=SNGL(G(4))
      END IF
C         Extract vu, vd at mt for mt calculation...
      IF (Q.GT.AMT) THEN
        VUMT=SNGL(G(30))
        VDMT=SNGL(G(31))
        ASMTP=SNGL(G(3)**2/4.D0/DPI)
      END IF
      IF (Q.GT.HIGFRZ) THEN
        ASMSS=SNGL(G(3)**2/4.D0/DPI)
        M3Q=SNGL(G(9))
      END IF
C          Freeze out running gluino mass at MGL
        IF (Q.LT.ABS(SNGL(G(9))).AND.IG(9).EQ.0) THEN
          G0(9)=SNGL(G(9))
          IG(9)=1
            ASM3=SNGL(G(3)**2/4.D0/DPI)
        ELSE IF (IG(9).EQ.0) THEN
          G0(9)=SNGL(G(9))
        END IF
C          Freeze out Higgs paremeters at HIGFRZ
      DO 211 I=13,14
        IF (Q.LT.HIGFRZ.AND.IG(I).EQ.0) THEN
          G0(I)=SNGL(G(I))
          IG(I)=1
          G0(I-6)=SNGL(G(I-6))
          IG(I-6)=1
          G0(I+12)=SNGL(G(I+12))
          IG(I+12)=1
          G0(I+17)=SNGL(G(I+17))
          IG(I+17)=1
        ELSE IF (IG(I).EQ.0) THEN
          G0(I)=SNGL(G(I))
          G0(I-6)=SNGL(G(I-6))
          G0(I+12)=SNGL(G(I+12))
          G0(I+17)=SNGL(G(I+17))
        END IF
211   CONTINUE
C          Freeze out 1st/2nd gen. soft terms at own masses
      DO 220 I=15,19
        IF (Q.LT.SQRT(ABS(SNGL(G(I)))).AND.IG(I).EQ.0) THEN
          G0(I)=SNGL(G(I))
          IG(I)=1
        ELSE IF (IG(I).EQ.0) THEN
          G0(I)=SNGL(G(I))
        END IF
220   CONTINUE
C          Freeze out third gen. soft terms at HIGFRZ
      DO 221 I=20,24
        IF (Q.LT.HIGFRZ.AND.IG(I).EQ.0) THEN
          G0(I)=SNGL(G(I))
          IG(I)=1
        ELSE IF (IG(I).EQ.0) THEN
          G0(I)=SNGL(G(I))
        END IF
221   CONTINUE
C          Freeze our N_R parameters at Majorana mass scale
      DO 230 I=27,29
        IF (G(I).NE.0.D0) G0(I)=SNGL(G(I))
        IF (Q.LT.AMNRMJ.AND.IG(I).EQ.0.) THEN
          IG(I)=1
        END IF
230   CONTINUE
100   RETURN
      END
