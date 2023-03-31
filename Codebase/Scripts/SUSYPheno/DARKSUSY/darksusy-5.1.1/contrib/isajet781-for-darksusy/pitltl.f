CDECK  ID>, PITLTL.
        REAL FUNCTION PITLTL(P2,G1,G2,G3,CTHW)
C-----------------------------------------------------------------------
C          PITLTL: t_L squark self-energy
C     Taken from Damien M. Pierce, Jonathan A. Bagger, Konstantin T. Matchev,
C     Ren-jie Zhang, Nucl.Phys.B491:3-67,1997, hep-ph/9606211
C          P2 = 4-momentum squared
C          CTHW = Cos(theta_W) in DR bar scheme
C-----------------------------------------------------------------------
      IMPLICIT NONE
C          Standard model parameters
C          AMUP,...,AMTP        = quark masses
C          AME,AMMU,AMTAU       = lepton masses
C          AMW,AMZ              = W,Z masses
C          GAMW,GAMZ            = W,Z widths
C          ALFAEM,SN2THW,ALFA3  = SM couplings
C          ALQCD4               = 4 flavor lambda
      COMMON/SSSM/AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      REAL AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      SAVE /SSSM/
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
      COMMON/SSINF/XLAM
      DOUBLE PRECISION XLAM
      COMPLEX*16 SSB0,SSG,SSF
      DOUBLE PRECISION SSA0
      COMPLEX TEMP,N(4,4),AC0TTL(4),BC0TTL(4),FTTLL(4),GTTLL(4)
      REAL P2,P,GG,GGP,CTHW,CHWW2,BE,FT,FB,C(4),DU(4),G1,G2,G3
     $,MH0(4),MHP(2),LHTLTI(4,2),LHTLBI(2,2),FBTLL(4),VMAT(2,2)
     $,UMAT(2,2),GBTLL(4),AP0TTL(4),BP0TTL(4)
     $,APPBTL(2),BPPBTL(2),ACPBTL(2),BCPBTL(2)
     $,LSTLTL(2),LSTLTR(2),LHTLTL(2),LHTLTR(2),LSTLTI(2,2)
     $,LHTLBL(2),LHTLBR(2)
      REAL COST,SINT,COST2,SINT2,COSB,SINB,COSB2,SINB2
     $,COSL,SINL,COSL2,SINL2,THX,THY
      REAL SINA2,COSA2,SINBE2,COSBE2,I3UL,I3DL,I3EL,I3NL,YUR,YUL
     $,YDR,YDL,YER,YEL,YNL,EUL,SWW2,GUL,SINA,COSA,SINBE,COSBE
      INTEGER I,II,THIW1,THIW2
      COMPLEX IMAG
      PARAMETER (IMAG=(0.,1.))
      REAL PI,SR2
      PI=4*ATAN(1.)
      SR2=SQRT(2.)
      P=SQRT(P2)
      CHWW2=CTHW**2
      COST=COS(THETAT)
      SINT=SIN(THETAT)
      COSB=COS(THETAB)
      SINB=SIN(THETAB)
      COSL=COS(THETAL)
      SINL=SIN(THETAL)
      COST2=COST**2
      SINT2=1-COST2
      COSB2=COSB**2
      SINB2=1-COSB2
      COSL2=COSL**2
      SINL2=1-COSL2
      GG=G2
      GGP=SQRT(3./5.)*G1
      BE=ATAN(VUQ/VDQ)
      FT=MTQ/VUQ
      FB=MBQ/VDQ
      SINA2=SIN(ALFAH)**2
      COSA2=COS(ALFAH)**2
      SINBE2=SIN(BE)**2
      COSBE2=COS(BE)**2
      SINA=SIN(ALFAH)
      COSA=COS(ALFAH)
      SINBE=SIN(BE)
      COSBE=COS(BE)
      C(1)=-(COSA2-SINA2)
      C(2)=COSA2-SINA2
      C(3)=-(COSBE2-SINBE2)
      C(4)=COSBE2-SINBE2
      DU(1)=SINA2
      DU(2)=COSA2
      DU(3)=SINBE2
      DU(4)=COSBE2
      MH0(1)=MSS(30)
      MH0(2)=MSS(29)
      MH0(3)=AMZ
      MH0(4)=MSS(31)
      MHP(2)=AMW
      MHP(1)=MSS(32)
      I3UL=1./2.
      I3DL=-1./2.
      I3EL=-1./2.
      I3NL=1./2.
      YUR=-4./3.
      YUL=1./3.
      YDR=2./3.
      YDL=1./3.
      YER=2.
      YEL=-1.
      YNL=-1.
      EUL=2./3.
      SWW2=1-(AMW/AMZ)**2
      GUL=I3UL-EUL*SWW2
      THX=SIGN(1.,TAN(GAMMAL))
      THY=SIGN(1.,TAN(GAMMAR))
      IF (SIGN(1.,AMW1SS).EQ.1.) THEN
         THIW1=0
      ELSE
         THIW1=1
      END IF
      IF (SIGN(1.,AMW2SS).EQ.1.) THEN
         THIW2=0
      ELSE
         THIW2=1
      END IF
      UMAT(1,1)=-(-1)**THIW1*SIN(GAMMAR)
      UMAT(1,2)=(-1)**THIW1*COS(GAMMAR)
      UMAT(2,1)=-(-1)**THIW2*THY*COS(GAMMAR)
      UMAT(2,2)=-(-1)**THIW2*THY*SIN(GAMMAR)
      VMAT(1,1)=-SIN(GAMMAL)
      VMAT(1,2)=COS(GAMMAL)
      VMAT(2,1)=-THX*COS(GAMMAL)
      VMAT(2,2)=-THX*SIN(GAMMAL)
      DO II=1,4
        IF (SIGN(1.,AMZISS(II)).EQ.1.) THEN
          I=0
        ELSE
          I=1
        END IF
        N(II,1)=-IMAG**I*ZMIXSS(4,II)
        N(II,2)=-IMAG**I*ZMIXSS(3,II)
        N(II,3)=IMAG**I*ZMIXSS(2,II)
        N(II,4)=IMAG**I*ZMIXSS(1,II)
      ENDDO
      DO I=1,4
        AP0TTL(I)=0.
        BP0TTL(I)=0.
        AC0TTL(I)=0.
        BC0TTL(I)=0.
      ENDDO
      BP0TTL(1)=GGP/SR2*YUL
      BP0TTL(2)=SR2*GG*I3UL
      AP0TTL(4)=FT
      DO I=1,4
        DO II=1,4
          AC0TTL(I)=AC0TTL(I)+CONJG(N(I,II))*AP0TTL(II)
          BC0TTL(I)=BC0TTL(I)+N(I,II)*BC0TTL(II)
        ENDDO
      ENDDO
      DO I=1,4
        FTTLL(I)=CONJG(AC0TTL(I))*AC0TTL(I)
     $+CONJG(BC0TTL(I))*BC0TTL(I)
        GTTLL(I)=CONJG(BC0TTL(I))*AC0TTL(I)
     $+CONJG(AC0TTL(I))*BC0TTL(I)
      ENDDO
      DO I=1,2
        APPBTL(I)=0.
        BPPBTL(I)=0.
        ACPBTL(I)=0.
        BCPBTL(I)=0.
      ENDDO
      APPBTL(1)=GG
      BPPBTL(2)=-FB
      DO I=1,2
        DO II=1,2
          ACPBTL(I)=ACPBTL(I)+VMAT(I,II)*APPBTL(II)
          BCPBTL(I)=BCPBTL(I)+UMAT(I,II)*BPPBTL(II)
        ENDDO
      ENDDO
      DO I=1,2
        FBTLL(I)=ACPBTL(I)*ACPBTL(I)+BCPBTL(I)*BCPBTL(I)
        GBTLL(I)=BCPBTL(I)*ACPBTL(I)+ACPBTL(I)*BCPBTL(I)
      ENDDO
      LSTLTL(1)=GG*AMZ/CTHW*GUL*COSBE
      LSTLTR(1)=-FT/SR2*TWOM1
      LSTLTL(2)=-GG*AMZ/CTHW*GUL*SINBE+SR2*FT*MTQ
      LSTLTR(2)=-FT/SR2*AAT
C     First index is for Higgs sector, second is for sfermion
      LSTLTI(1,1)=LSTLTL(1)*COST-LSTLTR(1)*SINT
      LSTLTI(1,2)=LSTLTL(1)*SINT+LSTLTR(1)*COST
      LSTLTI(2,1)=LSTLTL(2)*COST-LSTLTR(2)*SINT
      LSTLTI(2,2)=LSTLTL(2)*SINT+LSTLTR(2)*COST
      LHTLTI(1,1)=COSA*LSTLTI(1,1)-SINA*LSTLTI(2,1)
      LHTLTI(1,2)=COSA*LSTLTI(1,2)-SINA*LSTLTI(2,2)
      LHTLTI(2,1)=SINA*LSTLTI(1,1)+COSA*LSTLTI(2,1)
      LHTLTI(2,2)=SINA*LSTLTI(1,2)+COSA*LSTLTI(2,2)
      LHTLTL(1)=0.
      LHTLTR(1)=FT/SR2*(-TWOM1*COSBE-AAT*SINBE)
      LHTLTL(2)=0.
      LHTLTR(2)=-FT/SR2*(-TWOM1*SINBE+AAT*COSBE)
      LHTLTI(3,1)=LHTLTL(1)*COST-LHTLTR(1)*SINT
      LHTLTI(3,2)=LHTLTL(1)*SINT+LHTLTR(1)*COST
      LHTLTI(4,1)=LHTLTL(2)*COST-LHTLTR(2)*SINT
      LHTLTI(4,2)=LHTLTL(2)*SINT+LHTLTR(2)*COST
      LHTLBL(1)=-GG*AMW/SR2*(COSBE2-SINBE2)-FT*MTQ*SINBE
     $+FB*MBQ*COSBE
      LHTLBR(1)=FB*(-TWOM1*SINBE-AAB*COSBE)
      LHTLBL(2)=GG*AMW/SR2*2.*COSBE*SINBE-FT*MTQ*COSBE
     $-FB*MBQ*SINBE
      LHTLBR(2)=FB*(-TWOM1*COSBE+AAB*SINBE)
      LHTLBI(1,1)=LHTLBL(1)*COSB-LHTLBR(1)*SINB
      LHTLBI(1,2)=LHTLBL(1)*SINB+LHTLBR(1)*COSB
      LHTLBI(2,1)=LHTLBL(2)*COSB-LHTLBR(2)*SINB
      LHTLBI(2,2)=LHTLBL(2)*SINB+LHTLBR(2)*COSB
      TEMP=4.*G3**2/3.*(2.*SSG(P,MSS(1),AMTP)
     $+COST2*SSF(P,MSS(13),0.)+SINT2*SSF(P,MSS(12),0.)
     $+COST2*SSA0(MSS(13))+SINT2*SSA0(MSS(12)))
     $+FT**2*(SINT2*SSA0(MSS(13))+COST2*SSA0(MSS(12)))
     $+FB**2*(SINB2*SSA0(MSS(10))+COSB2*SSA0(MSS(11)))
     $+4.*GG**2/CHWW2*GUL**2*SSA0(AMZ)+2.*GG**2*SSA0(AMW)
     $+(2./3.)**2*4.*PI/137.036*(COST2*SSF(P,MSS(13),0.)
     $+SINT2*SSF(P,MSS(12),0.))+GG**2/CHWW2*GUL**2
     $*(COST2*SSF(P,MSS(13),AMZ)+SINT2*SSF(P,MSS(12),AMZ))
     $+GG**2/2.*(COSB2*SSF(P,MSS(10),AMW)+SINB2*SSF(P,MSS(11),AMW))
     $+GG**2/4.*(COST2*SSA0(MSS(13))+SINT2*SSA0(MSS(12))
     $+2.*(COSB2*SSA0(MSS(10))+SINB2*SSA0(MSS(11))))
     $+GG**2*I3UL*(3.*(I3UL*(SSA0(MSS(2))+SSA0(MSS(8))
     $+COST2*SSA0(MSS(13))+SINT2*SSA0(MSS(12)))
     $+I3DL*(SSA0(MSS(4))+SSA0(MSS(6))
     $+COSB2*SSA0(MSS(10))+SINB2*SSA0(MSS(11))))
     $+I3EL*(SSA0(MSS(17))+SSA0(MSS(19))
     $+COSL2*SSA0(MSS(21))+SINL2*SSA0(MSS(22)))
     $+I3NL*(SSA0(MSS(14))+SSA0(MSS(15))+SSA0(MSS(16))))
     $+GGP**2/4.*YUL**2*(COST2*SSA0(MSS(13))+SINT2*SSA0(MSS(12)))
     $+GGP**2/4.*YUL*(3.*(YUL*(SSA0(MSS(2))+SSA0(MSS(8))
     $+COST2*SSA0(MSS(13))+SINT2*SSA0(MSS(12)))
     $+YUR*(SSA0(MSS(3))+SSA0(MSS(9))+SINT2*SSA0(MSS(13))
     $+COST2*SSA0(MSS(12)))+YDL*(SSA0(MSS(4))+SSA0(MSS(6))
     $+COSB2*SSA0(MSS(10))+SINB2*SSA0(MSS(11)))
     $+YDR*(SSA0(MSS(5))+SSA0(MSS(7))
     $+SINB2*SSA0(MSS(10))+COSB2*SSA0(MSS(11))))
     $+YEL*(SSA0(MSS(17))+SSA0(MSS(19))
     $+COSL2*SSA0(MSS(21))+SINL2*SSA0(MSS(22)))
     $+YER*(SSA0(MSS(18))+SSA0(MSS(20))
     $+SINL2*SSA0(MSS(21))+COSL2*SSA0(MSS(22)))
     $+YNL*(SSA0(MSS(14))+SSA0(MSS(15))+SSA0(MSS(16))))
      DO I=1,4
        TEMP=TEMP+(FT**2*DU(I)-GG**2*GUL/2./CHWW2*C(I))
     $*SSA0(MH0(I))/2.
      ENDDO
      DO I=3,4
        TEMP=TEMP+(FB**2*DU(I)+GG**2*(GUL/2./CHWW2-I3UL)*C(I))
     $*SSA0(MHP(I-2))
      ENDDO
      DO I=1,4
          TEMP=TEMP+(LHTLTI(I,1))**2*SSB0(P2,MH0(I),MSS(13))
     $+(LHTLTI(I,2))**2*SSB0(P2,MH0(I),MSS(12))
      ENDDO
      DO I=1,2
        DO II=1,2
          TEMP=TEMP+(LHTLBI(I,II))**2*SSB0(P2,MSS(9+II),MHP(I))
        ENDDO
      ENDDO
      DO I=1,4
        TEMP=TEMP+FTTLL(I)*SSG(P,ABS(MSS(22+I)),AMTP)
     $-2.*GTTLL(I)*ABS(MSS(22+I))*MTQ*SSB0(P2,ABS(MSS(22+I)),AMTP)
      ENDDO
      DO I=1,2
        TEMP=TEMP+FBTLL(I)*SSG(P,ABS(MSS(26+I)),AMBT)
     $-2.*GBTLL(I)*ABS(MSS(26+I))*MBQ*SSB0(P2,ABS(MSS(26+I)),AMBT)
      ENDDO
      PITLTL=REAL(TEMP)/16./PI**2
      RETURN
      END
