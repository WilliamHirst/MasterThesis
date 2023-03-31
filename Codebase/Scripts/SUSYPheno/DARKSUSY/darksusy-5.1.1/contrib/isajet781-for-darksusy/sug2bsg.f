C---------------------------------------------------------------------
      SUBROUTINE SUG2BSG            
C----------------------------------------------------------------------
C     SUGRA - BSG interface.
C     Fills BSG common blocks from SUGRA ones.
C
C    Note: some of the parameters are transferred to the main IsaBSG
C          routine directly
C
C     Created: 02/22/07 by Azar Mustafayev
C     Modified: 6/12/07 by Azar Mustafayev - compatibility with ISAJET 7.75
C
      IMPLICIT NONE
c************************
c   ISAJET common blocks
c************************
C     XSUGIN contains the inputs to SUGRA:
C     XSUGIN(1) = M_0        XSUGIN(2) = M_(1/2)  XSUGIN(3) = A_0
C     XSUGIN(4) = tan(beta)  XSUGIN(5) = sgn(mu)  XSUGIN(6) = M_t
C     XSUGIN(7) = SUG BC scale
C     XGMIN(1) = LAM         XGMIN(2)  = M_MES    XGMIN(3)  = XN5
C     XGMIN(4) = tan(beta)   XGMIN(5)  = sgn(mu)  XGMIN(6) = M_t
C     XGMIN(7) = CGRAV       XGMIN(8)  =RSL       XGMIN(9)  = DEL_HD
C     XGMIN(10)  = DEL_HU    XGMIN(11) = DY       XGMIN(12) = N5_1
C     XGMIN(13)  = N5_2      XGMIN(14) = N5_3
C     XNRIN(1) = M_N3        XNRIN(2) = M_MAJ     XNRIN(3) = ANSS 
C     XNRIN(4) = M_N3SS
C     XISAIN contains the MSSMi inputs in natural order.
      COMMON /SUGXIN/ XISAIN(24),XSUGIN(7),XGMIN(14),XNRIN(4),
     $XAMIN(11)
      REAL XISAIN,XSUGIN,XGMIN,XNRIN,XAMIN
      SAVE /SUGXIN/
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
C     XNUSUG contains non-universal GUT scale soft terms for SUGRA:
C     XNUSUG(1)=M1 XNUSUG(2)=M2 XNUSUG(3)=M3
C     XNUSUG(4)=A_tau XNUSUG(5)=A_b XNUSUG(6)=A_t
C     XNUSUG(7)=m_Hd XNUSUG(8)=m_Hu XNUSUG(9)=m_eR XNUSUG(10)=m_eL
C     XNUSUG(11)=m_dR XNUSUG(12)=m_uR XNUSUG(13)=m_uL XNUSUG(14)=m_lR
C     XNUSUG(15)=m_lL XNUSUG(16)=m_bR XNUSUG(17)=m_tR XNUSUG(18)=m_tL
C     XNUSUG(19)=mu(Q) XNUSUG(20)=mA(Q)
      COMMON /SUGNU/ XNUSUG(20),INUHM
      REAL XNUSUG
      INTEGER INUHM
      SAVE /SUGNU/
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
c+CDE,SSNU
      COMMON /BSG/GISA(31),MSQISA(3),MSLISA(3),MSUISA(3),MSDISA(3),
     &            MSEISA(3),MRNISA(3),YNFRZ(3,3),MNFRZ(3,3),TNFRZ(3,3),
     &            RSIGT,RSIGB,RSIGL
c       GISA(i) - values of RGE parameters at MZ in DRbar:
C     GISA( 1) = g_1        GISA( 2) = g_2        GISA( 3) = g_3
C     GISA( 4) = y_tau      GISA( 5) = y_b        GISA( 6) = y_t
C     GISA( 7) = M_1        GISA( 8) = M_2        GISA( 9) = M_3
C     GISA(10) = A_tau      GISA(11) = A_b        GISA(12) = A_t
C     GISA(13) = M_hd^2     GISA(14) = M_hu^2     GISA(15) = M_er^2
C     GISA(16) = M_el^2     GISA(17) = M_dnr^2    GISA(18) = M_upr^2
C     GISA(19) = M_upl^2    GISA(20) = M_taur^2   GISA(21) = M_taul^2
C     GISA(22) = M_btr^2    GISA(23) = M_tpr^2    GISA(24) = M_tpl^2
C     GISA(25) = mu         GISA(26) = B          GISA(27) = Y_N
C     GISA(28) = M_nr       GISA(29) = A_n        GISA(30) = vdq
C     GISA(31) = vuq
c
c     MSxDEC(i) - decoupling scale of i-th generation of type x sfermion
c     MRNDEC(i) - decoupling scale of i-th RH neutrino
c     RSIGT,RSIGB,RSIGL - radiative corrections to top, bottom and tau
c                         Yukawas at MSUSY
      REAL*8 GISA,MSQISA,MSLISA,MSUISA,MSDISA,MSEISA,MRNISA,
     &       YNFRZ,MNFRZ,TNFRZ
      REAL RSIGT,RSIGB,RSIGL
      SAVE /BSG/
c************************
c   ISABSG common blocks
c************************
      COMMON /BSGSM/ MZ,MW,MB,MC,MS,MT,MTAU,BSGXW,S12,S23,S13,BSGAEM,
     &               BSGTHW
      REAL*8 MZ,MW,MB,MC,MS,MT,MTAU,BSGXW,S12,S23,S13,BSGAEM,BSGTHW
c  NB: The following name changes has to be made in order not to interfere
c      with ISAJET names, but only in this subroutine. In all other
c      subroutines original names are restored.
c       XW     -> BSGXW
c       ALFAEM -> BSGAEM
c       SN2THW -> BSGTHW
      SAVE /BSGSM/
      COMMON /GGN/ M1,M2,M3,ABOT,ATOP,ATAU
      REAL*8 M1,M2,M3,ABOT,ATOP,ATAU
      SAVE /GGN/      
      COMMON /GLNN/ MCH0(4),MDSB(6)
      REAL*8 MCH0,MDSB
      SAVE /GLNN/
      COMMON /BSGDEC/MSQDEC(3),MSLDEC(3),MSUDEC(3),MSDDEC(3),
     &               MSEDEC(3),MRNDEC(3),BSGRHN
      REAL*8 MSQDEC,MSLDEC,MSUDEC,MSDDEC,MSEDEC,MRNDEC
      INTEGER BSGRHN
      SAVE /BSGDEC/
      COMMON /BSGSUG/ TANB,BSGV,BSGVP,BSGQ,BSGMU,MSTP1,MSTP2,MSCHL,
     &                MSCHR,MSUPL,MSEL,MSW1,MGLU,MHPLUS,MHA0,BSGMZI(4),
     &                BSGZMIX(4,4),BSGMW1,BSGMW2,BSGGAML,BSGGAMR,
     &                BSGTHET,BSGMTQ,BSGMBQ,
     &                MSTLQ,MSTRQ,BSGMGUT,BSGNGUT,BSGFT,BSGNRIN(21),
     &                BSGGMIN(14),BSGGUN,BSGNUS(20),BSGAMIN(7),
     &                BSGENU,BSGFN(3),BSGMHU,BSGMHD,
     &                BSGNUHM,BSG3UN,BSG5ON
      REAL*8 TANB,BSGV,BSGVP,BSGQ,BSGMU,MSTP1,MSTP2,MSCHL,MSCHR,MSUPL,
     &       MSEL,MSW1,MGLU,MHPLUS,MHA0,BSGMZI,BSGZMIX,BSGMW1,BSGMW2,
     &       BSGGAML,BSGGAMR,BSGTHET,BSGMTQ,BSGMBQ,
     &       MSTLQ,MSTRQ,BSGMGUT,BSGNGUT,BSGFT,BSGNRIN,BSGGMIN,BSGGUN,
     &       BSGNUS,BSGAMIN,BSGENU,BSGFN,BSGMHU,BSGMHD
      INTEGER BSG3UN,BSGNUHM
      LOGICAL BSG5ON
c  NB: The following name changes has to be made in order not to interfere
c      with ISAJET names, but only in this subroutine. In all other
c      subroutines original names are restored.
c   V      -> BSGV               XRHNIN -> BSGNRIN
c   VP     -> BSGVP              XGMIN  -> BSGMIN
c   MU     -> BSGMU              IAL3UN -> BSG3UN
c   MSUSY  -> BSGQ               GGUTSS -> BSGGUN
c   MTQ    -> BSGMTQ             XNUSUG -> BSGNUS
c   MBQ    -> BSGMBQ             XAMIN  -> BSGAMIN
c   MGUT   -> BSGMGUT            AMZISS -> BSGMZI
c   FNGUT  -> BSGNGUT            ZMIXSS -> BSGZMIX
c   FTMT   -> BSGFT              AMWiSS -> BSGMWi
c   THETAT -> BSGTHET            EPSNU  -> BSGENU
c   GAMMAL -> BSGGAML            FTRHLD -> BSGFN
c   GAMMAR -> BSGGAMR            LND5ON -> BSG5ON
c   IRHN   -> BSGRHN             INUHM  -> BSGNUHM
c   MHUSMG -> BSGMHU             MHDSMG -> BSGMHD
      SAVE /BSGSUG/
      INTEGER I,J
      
      
      M1 = GISA(7)
      M2 = GISA(8)
      M3 = GISA(9)
      ABOT = GISA(11)
      ATOP = GISA(12)
      ATAU = GISA(10)
      
      MZ = DBLE(AMZ)
      MW = DBLE(AMW)
      MB = DBLE(AMBT)
      MC = 1.5d0
      MS = 0.2d0
      MT = DBLE(AMTP)
      MTAU = DBLE(AMTAU)
      
      S12 = 0.22d0
      S23 = 0.04d0
      S13 = 0.003d0
      
      BSGXW = DBLE(XW)
      BSGAEM = DBLE(ALFAEM)
      BSGTHW = DBLE(SN2THW)
      
      TANB = DBLE(XTANB)
      BSGV = DBLE(V)            ! v_u
      BSGVP = DBLE(VP)          ! v_d 
      BSGMU = DBLE(MU)
      
      BSGQ = DBLE(MSUSY)
      BSGMGUT = DBLE(MGUT)
      
      MSTP1 = DBLE(MSS(12))
      MSTP2 = DBLE(MSS(13))
      BSGTHET = DBLE(THETAT)
      
      MSCHL = DBLE(MSS(8))
      MSCHR = DBLE(MSS(9))
      MSUPL = DBLE(MSS(2))
      MSEL  = DBLE(MSS(17))
      MSW1  = DBLE(MSS(27))
      MGLU  = DBLE(MSS(1))
      MHPLUS= DBLE(MSS(32))
      MHA0  = DBLE(MSS(31))
      
      BSGMZI(1)=DBLE(AMZ1SS)
      BSGMZI(2)=DBLE(AMZ2SS)
      BSGMZI(3)=DBLE(AMZ3SS)
      BSGMZI(4)=DBLE(AMZ4SS)
      DO I=1,4
        DO J=1,4
          BSGZMIX(I,J)=DBLE(ZMIXSS(I,J))
        ENDDO
      ENDDO
      BSGGAML = DBLE(GAMMAL)
      BSGGAMR = DBLE(GAMMAR)
      BSGMW1 = DBLE(AMW1SS)
      BSGMW2 = DBLE(AMW2SS)
      
      BSGMTQ = DBLE(MTQ)
      BSGMBQ = DBLE(MBQ)
      MSTLQ = DBLE(AMTLSS)
      MSTRQ = DBLE(AMTRSS)
      
      BSGFT = DBLE(FTMT)
      BSGNGUT = DBLE(FNGUT)
      
c      DO I=1,21
c        BSGNRIN(I) = DBLE(XRHNIN(I))
c      ENDDO
c      BSGENU = EPSNU
c      DO I=1,3
c        BSGFN(I) = FTRHLD(I)
c      ENDDO
c      BSGRHN = IRHN
c      BSG5ON = LND5ON
      BSGENU = 1.d0
      IF (XNRIN(2).LT.1.E19) THEN 
        BSGRHN=1
      ELSE
        BSGRHN=0
      ENDIF
      BSG5ON = .FALSE.
      BSGFN(1)=0.d0
      BSGFN(2)=0.d0
      BSGFN(3)=XNRIN(2)
      
      
      DO I=1,14
        BSGGMIN(I) = DBLE(XGMIN(i))
      ENDDO
      DO I=1,20
        BSGNUS(I) = DBLE(XNUSUG(i))
      ENDDO
      DO I=1,7
        BSGAMIN(I) = DBLE(XAMIN(i))
      ENDDO
      
      BSG3UN = IAL3UN
      BSGGUN = DBLE(GGUTSS)
      
      BSGNUHM = INUHM
      BSGMHU = DBLE(MHUSMG)
      BSGMHD = DBLE(MHDSMG)

      DO I=1,3
        MSQDEC(I)=MSQISA(I)
	MSLDEC(I)=MSLISA(I)
	MSUDEC(I)=MSUISA(I)
	MSDDEC(I)=MSDISA(I)
	MSEDEC(I)=MSEISA(I)
	MRNDEC(I)=MRNISA(I)
      ENDDO
      
      RETURN
      END
