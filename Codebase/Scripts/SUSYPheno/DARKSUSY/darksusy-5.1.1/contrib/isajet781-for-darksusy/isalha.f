CDECK  ID>, ISALHA.
C--------------------------------------------------------------------
      SUBROUTINE ISALHA(LOUT,IMODEL,IMODIN,MT)
C--------------------------------------------------------------------
C
C     Print SUGRA parameters in 'SUSY Les Houches Accord 2' (SLHA2) format
C      C. Balazs, Apr. 21 2009, v2.0
C      C. Balazs, Jan.  5 2005, v1.0
C      C. Balazs, Jul. 24 2003, v0.1
C
C     LOUT   = Output file ID#
C     IMODEL = model type for SUGRA
C     IMODIN = input model type to control formatting
C     MT     = top mass
C
      IMPLICIT NONE
CsB   ISAJET common blocks from SUGPRT ...
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
CsB   Additional ISAJET common blocks ...
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
      COMMON/CONST/PI,SQRT2,ALFA,GF,UNITS
      SAVE /CONST/
      REAL      PI,SQRT2,ALFA,GF,UNITS
CsB   End of ISAJET common blocks
C
      REAL GPX,SIN2W,ALEMI,AS,MT,ASMB,MBMB,SUALFS
      INTEGER IMODEL,J,K,IMODIN,LOUT
      INTEGER I,I1,I2,J1,J2
      REAL RMSUSY,SG,THX,THY
C
CsB   Local ISAJET related variables
      CHARACTER*40 VERSN,VISAJE
C
CsB   Local LHA3 related variables
      Integer iPDG(33),ISA2LHA3(33),I2L3GSS(50),IModSLHA(10)
      DIMENSION CHAF(33),ModelDescr(10),SoftParaLHA(50)
      CHARACTER CHAF*16, ModelDescr*50,SoftParaLHA*16
CsB   Particle names for SLHA (in this order)
      DATA CHAF /
     $' top',' h^0',' H^0',' A^0',' H^+',
     $' dnl',' upl',' stl',' chl',' b1',' t1',
     $' el-',' nuel',' mul-',' numl',' tau1',' nutl',
     $' glss',' z1ss',' z2ss',' w1ss',' z3ss',' z4ss',' w2ss',
     $' dnr',' upr',' str',' chr',' b2',' t2',
     $' er-',' mur-',' tau2' /
CsB   PDG codes of the above
      DATA iPDG /
     &      6,     25,     35,     36,     37,
     &1000001,1000002,1000003,1000004,1000005,1000006,
     &1000011,1000012,1000013,1000014,1000015,1000016,
     &1000021,1000022,1000023,1000024,1000025,1000035,1000037,
     &2000001,2000002,2000003,2000004,2000005,2000006,
     &2000011,2000013,2000015 /
CsB   MSS indices of the above
      Data ISA2LHA3 /
     $ 0,29,30,31,32,
     $ 4, 2, 6, 8,10,12,
     $17,14,19,15,21,16,
     $ 1,23,24,27,25,26,28,
     $ 5, 3, 7, 9,11,13,
     $18,20,22/
CsB   Soft parameters for SLHA (in this order)
      DATA SoftParaLHA /
     ,'M_1(Q)','M_2(Q)','M_3(Q)','     ','     ',
     ,'      ','      ','      ','     ','     ',
     ,'      ','      ','      ','     ','     ',
     ,'      ','      ','      ','     ','     ',
     ,'      ','      ','      ','     ','     ',
     ,'      ','      ','      ','     ','     ',
     ,'MeL(Q)','MmuL(Q)','MtauL(Q)','MeR(Q)','MmuR(Q)',
     ,'MtauR(Q)','      ','      ','     ','     ',
     ,'MqL1(Q)','MqL2(Q)','MqL3(Q)','MuR(Q)','McR(Q)',
     ,'MtR(Q)','MdR(Q)','MsR(Q)','MbR(Q)','    '/
CsB   GSS indices of the above
      Data I2L3GSS /
     $  7, 8, 9, 0, 0, 0, 0, 0, 0, 0,
     $  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     $  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     $ 16,16,21,15,15,20, 0, 0, 0, 0,
     $ 19,19,24,18,18,23,17,17,22, 0/
CsB   Short model descriptions
      Data ModelDescr /
     $'Minimal supergravity (mSUGRA) model',
     $'Minimal gauge mediated (GMSB) model',
     $'Non-universal supergravity model',
     $'Supergravity model with truly unified couplings',
     $'Non-minimal gauge mediated (GMSB) model',
     $'Supergravity model with right-handed neutrinos',
     $'Minimal anomaly mediated (AMSB) model',
     $'Non-minimal anomaly mediated (AMSB) model',
     $'Mixed moduli-anomaly mediated model',
     $'Hypercharged anomaly mediation (HCAMSB) model'
     $/
CsB   ISAJET -> SLHA model numbering
      Data IModSLHA / 1,2,13,14,15,16, 3,18,19,20/
C
      LOGICAL TESTING
      TESTING = .FALSE.
C
      PI=4.*ATAN(1.)
      GPX=SQRT(.6)*GSS(1)
      SIN2W=GPX**2/(GSS(2)**2+GPX**2)
      ALEMI=4*PI/GSS(2)**2/SIN2W
      AS=.1172
C
CsB   SLHA output file is assumed to be open
C
CsB   Write LHA3 header
      WRITE(LOUT,7000)
     . ' ISAJET SUSY parameters in SUSY Les Houches Accord 2 format'
      WRITE(LOUT,7000)
     . ' Created by ISALHA 2.0 Last revision: C. Balazs 21 Apr 2009'
      VERSN=VISAJE()
      VERSN=VERSN(14:)
C
      WRITE(LOUT,7001)    'SPINFO',
     ,                    'Program information'
      WRITE(LOUT,7012) 1, 'ISASUGRA from ISAJET       ',
     ,                    'Spectrum Calculator'
      WRITE(LOUT,7012) 2,  VERSN,
     ,                    'Version number'
C
      WRITE(LOUT,7001) 'MODSEL', 'Model selection'
      WRITE(LOUT,7010) 1, IModSLHA(IMODIN), ModelDescr(IMODIN)
C
      Call SETCON
      ASMB=SUALFS(AMBT**2,.36,AMTP,3)
      MBMB=4.2
      WRITE(LOUT,7001) 'SMINPUTS', 'Standard Model inputs'
      WRITE(LOUT,7011) 1, ALEMI, 'alpha_em^(-1)' ! (MZ) SM MSbar'
      WRITE(LOUT,7011) 2,    GF, 'G_Fermi'
      WRITE(LOUT,7011) 3,    AS, 'alpha_s(M_Z)'
      WRITE(LOUT,7011) 4,   AMZ, 'm_{Z}(pole)'
C     WRITE(LOUT,7011) 5,  AMBT, 'm_{b}(pole)'
      WRITE(LOUT,7011) 5,  MBMB, 'm_{b}(m_{b})'
      WRITE(LOUT,7011) 6,  AMTP, 'm_{top}(pole)'
      WRITE(LOUT,7011) 7, AMTAU, 'm_{tau}(pole)'
C
      WRITE(LOUT,7001) 'MINPAR', 'SUSY breaking input parameters'
C     Print inputs and GUT couplings for SUGRA/AMSB models
      IF(IMODEL.EQ.1.OR.IMODEL.EQ.7.OR.IMODEL.EQ.9.OR.IMODEL.EQ.10) THEN
        IF(IMODEL.EQ.1) THEN
          WRITE(LOUT,7011) 1, XSUGIN(1), 'm_0'
          WRITE(LOUT,7011) 2, XSUGIN(2), 'm_{1/2}'
          WRITE(LOUT,7011) 3, XSUGIN(4), 'tan(beta)'
          WRITE(LOUT,7011) 4, XSUGIN(5), 'sign(mu)'
          WRITE(LOUT,7011) 5, XSUGIN(3), 'A_0'
C          WRITE(LOUT,1000) XSUGIN(1),XSUGIN(2),XSUGIN(3),XSUGIN(4),
C     $    XSUGIN(5),XSUGIN(6)
C1000      FORMAT(
C     $    ' M_0,  M_(1/2),  A_0,  tan(beta),  sgn(mu),  M_t ='
C     $    /4F10.3,2X,F6.1,F10.3)
        ELSE IF (IMODEL.EQ.7) THEN
          WRITE(LOUT,7011) 1, XSUGIN(1), 'm_0'
          WRITE(LOUT,7011) 2, XSUGIN(2), 'm_{3/2}'
          WRITE(LOUT,7011) 3, XSUGIN(4), 'tan(beta)'
          WRITE(LOUT,7011) 4, XSUGIN(5), 'sign(mu)'
C          WRITE(LOUT,1018) XSUGIN(1),XSUGIN(2),XSUGIN(4),XSUGIN(5),
C     $    XSUGIN(6)
C1018      FORMAT(
C     $    ' M_0,  M_(3/2),  tan(beta),  sgn(mu),  M_t ='
C     $    /3F10.3,2X,F6.1,2F10.3)
        ELSE IF (IMODEL.EQ.9) THEN
          WRITE(LOUT,7011) 10, XSUGIN(1), 'alpha'
          WRITE(LOUT,7011)  2, XSUGIN(2), 'm_{3/2}'
          WRITE(LOUT,7011)  3, XSUGIN(4), 'tan(beta)'
          WRITE(LOUT,7011)  4, XSUGIN(5), 'sign(mu)'
          WRITE(LOUT,7011) 11, XAMIN( 1), 'Moduli nQ'
          WRITE(LOUT,7011) 12, XAMIN( 2), '       nD'
          WRITE(LOUT,7011) 13, XAMIN( 3), '       nU'
          WRITE(LOUT,7011) 14, XAMIN( 4), '       nL'
          WRITE(LOUT,7011) 15, XAMIN( 5), '       nE'
          WRITE(LOUT,7011) 16, XAMIN( 6), '      nHd'
          WRITE(LOUT,7011) 17, XAMIN( 7), '      nHu'
          WRITE(LOUT,7011) 18, XAMIN( 8), '       L1'
          WRITE(LOUT,7011) 19, XAMIN( 9), '       L2'
          WRITE(LOUT,7011) 20, XAMIN(10), '       L3'
C          WRITE(LOUT,1019) XSUGIN(1),XSUGIN(2),XSUGIN(4),XSUGIN(5),
C     $    XSUGIN(6)
C1019      FORMAT(
C     $    ' alpha,  M_(3/2),  tan(beta),  sgn(mu),  M_t ='
C     $    /3F10.3,2X,F6.1,2F10.3)
C          WRITE(LOUT,1020) XAMIN(1),XAMIN(2),XAMIN(4),XAMIN(5),
C     $    XAMIN(6),XAMIN(7),XAMIN(8),XAMIN(9),XAMIN(10)
C1020      FORMAT(
C     $    ' Moduli nQ, nD, nU, nL, nE, nHd, nHu, L1, L2, L3 ='
C     $    /10F4.1)
        ELSE IF (IMODEL.EQ.10) THEN
          WRITE(LOUT,7011) 10, XAMIN(11), 'alpha'
          WRITE(LOUT,7011)  2, XSUGIN(2), 'm_{3/2}'
          WRITE(LOUT,7011)  3, XSUGIN(4), 'tan(beta)'
          WRITE(LOUT,7011)  4, XSUGIN(5), 'sign(mu)'
C          WRITE(LOUT,1040) XAMIN(11),XSUGIN(2),XSUGIN(4),XSUGIN(5),
C     $    XSUGIN(6)
C1040      FORMAT(
C     $    ' alpha,  M_(3/2),  tan(beta),  sgn(mu),  M_t ='
C     $    /3F10.3,2X,F6.1,2F10.3)
        END IF
C
        WRITE(LOUT,7001) 'EXTPAR', 
     $   'Non-universal SUSY breaking parameters'
C       Write out non-universal GUT scale parameters
        IF(XSUGIN(7).NE.0) THEN !!! Check this !!!
          WRITE(LOUT,7011) 0, XSUGIN(7), 'Input scale'
C          WRITE(LOUT,1026) XSUGIN(7)
C1026      FORMAT(' Q_max= ',E12.4)
        ELSE
          WRITE(LOUT,7011) 0, MGUTSS, 'Input scale'
        ENDIF
        IF (INUHM.EQ.1) THEN
C         WRITE(LOUT,7011) 21, MHDSQ, 'Down type Higgs mass squared'
C         WRITE(LOUT,7011) 22, MHUSQ, ' Up  type Higgs mass squared'
          WRITE(LOUT,7011) 21, MHDSMG, 'Down type Higgs mass squared'
          WRITE(LOUT,7011) 22, MHUSMG, ' Up  type Higgs mass squared'
C          WRITE(LOUT,1023)
C          WRITE(LOUT,1021) MHDSQ,MHUSQ
C          WRITE(LOUT,1022) MHDSMG,MHUSMG
C1021      FORMAT(/,' M_Hd^2(Q)= ',E10.3,3X,'M_Hu^2(Q)= ',E10.3)
C1022      FORMAT(' M_Hd^2(MGUT)= ',E10.3,3X,'M_Hu^2(MGUT)= ',E10.3)
C1023      FORMAT(/,' NUHM model has been selected:')
        END IF
        IF(XNUSUG(1).LT.1.E19.OR.XNUSUG(2).LT.1.E19.OR.XNUSUG(3)
     $  .LT.1.E19) THEN
          WRITE(LOUT,7011) 1, XNUSUG(1), ' U(1)_Y gaugino (Bino) mass'
          WRITE(LOUT,7011) 2, XNUSUG(2), 'SU(2)_L gaugino (Wino) mass'
          WRITE(LOUT,7011) 3, XNUSUG(3), 'SU(3)_C gaugino (gluino) mass'
C          WRITE(LOUT,1010) XNUSUG(1),XNUSUG(2),XNUSUG(3)
C1010      FORMAT(/' M_1(GUT)= ',F8.2,'    M_2(GUT)= ',F8.2,
C     $    '    M_3(GUT)= ',F8.2)
        END IF
        IF(XNUSUG(4).LT.1.E19.OR.XNUSUG(5).LT.1.E19.OR.XNUSUG(6)
     $  .LT.1.E19) THEN
          WRITE(LOUT,7011) 11, XNUSUG(6), '  Top  trilinear coupling'
          WRITE(LOUT,7011) 12, XNUSUG(5), 'Bottom trilinear coupling'
          WRITE(LOUT,7011) 13, XNUSUG(4), '  Tau  trilinear coupling'
C          WRITE(LOUT,1011) XNUSUG(4),XNUSUG(5),XNUSUG(6)
C1011      FORMAT(/' A_tau(GUT)= ',F8.2,'    A_b(GUT)= ',F8.2,
C     $    '    A_t(GUT)= ',F8.2)
        END IF
        IF(XNUSUG(7).LT.1.E19.OR.XNUSUG(8).LT.1.E19) THEN
          WRITE(LOUT,7011) 21, XNUSUG(7)**2, 'Down type Higgs mass^2'
          WRITE(LOUT,7011) 22, XNUSUG(8)**2, ' Up  type Higgs mass^2'
C          WRITE(LOUT,1012) XNUSUG(7),XNUSUG(8)
C1012      FORMAT(/' M_Hd(GUT)= ',F8.2,'    M_Hu(GUT)= ',F8.2)
        END IF
        IF (XNUSUG(9).LT.1.E19.OR.XNUSUG(10).LT.1.E19) THEN
          WRITE(LOUT,7011) 34, XNUSUG( 9), 'Right scalar electron mass'
          WRITE(LOUT,7011) 31, XNUSUG(10), 'Left 1st gen. slepton mass'
C          WRITE(LOUT,1013) XNUSUG(9),XNUSUG(10)
C1013      FORMAT(/' M_eR(GUT)= ',F8.2,'    M_eL(GUT)= ',F8.2)
        END IF
        IF(XNUSUG(11).LT.1.E19.OR.XNUSUG(12).LT.1.E19.OR.XNUSUG(13)
     $  .LT.1.E19) THEN
          WRITE(LOUT,7011) 47, XNUSUG(11), 'Right scalar down mass'
          WRITE(LOUT,7011) 44, XNUSUG(12), 'Right scalar  up  mass'
          WRITE(LOUT,7011) 41, XNUSUG(13), 'Left 1st gen. squark mass'
C          WRITE(LOUT,1014) XNUSUG(11),XNUSUG(12),XNUSUG(13)
C1014      FORMAT(' M_dR(GUT)= ',F8.2,'    M_uR(GUT)= ',F8.2,
C     $    '    M_uL(GUT)=',F8.2)
        END IF
        IF(XNUSUG(14).LT.1.E19.OR.XNUSUG(15).LT.1.E19) THEN
          WRITE(LOUT,7011) 36, XNUSUG(14), 'Right scalar tau mass'
          WRITE(LOUT,7011) 33, XNUSUG(15), 'Left 3rd gen. slepton mass'
C          WRITE(LOUT,1015) XNUSUG(14),XNUSUG(15)
C1015      FORMAT(/' M_tauR(GUT)= ',F8.2,'    M_tauL(GUT)= ',F8.2)
        END IF
        IF(XNUSUG(16).LT.1.E19.OR.XNUSUG(17).LT.1.E19.OR.XNUSUG(18)
     $  .LT.1.E19) THEN
          WRITE(LOUT,7011) 49, XNUSUG(16), 'Right scalar bottom mass'
          WRITE(LOUT,7011) 46, XNUSUG(17), 'Right scalar top mass'
          WRITE(LOUT,7011) 43, XNUSUG(18), 'Left 3rd gen. squark mass'
C          WRITE(LOUT,1016) XNUSUG(16),XNUSUG(17),XNUSUG(18)
C1016      FORMAT(' M_bR(GUT)= ',F8.2,'    M_tR(GUT)= ',F8.2,
C     $    '    M_tL(GUT)=',F8.2)
        END IF
C       Right-handed neutrino parameters
        IF (XNRIN(2).LT.1.E19) THEN
          WRITE(LOUT,7011) 101, XNRIN(1), 'M(nu_tau)'
          WRITE(LOUT,7011) 102, XNRIN(2), 'M(N_R)'
          WRITE(LOUT,7011) 103, XNRIN(3), 'A_N'
          WRITE(LOUT,7011) 104, XNRIN(4), 'M(NRSS)'
          WRITE(LOUT,7011) 105,     FNMZ, 'FN(M_Z)'
          WRITE(LOUT,7011) 106,    FNGUT, 'FN(M_{GUT})'
C          WRITE(LOUT,1017) XNRIN(1),XNRIN(2),XNRIN(3),XNRIN(4),
C     $    FNMZ,FNGUT
C1017      FORMAT(' Right-handed neutrino parameters:'/
C     $    ' M(nu_tau)=',E10.3,'   M(N_R) =',E10.3,
C     $    '   A_N=',F8.2,'   M(NRSS)=',F8.2/
C     $    ' FN(M_Z)  =',F8.4, '   FN(M_GUT) =',F8.4)
        END IF
CsB     Non-minimal parameters for AMSB 
        IF (IMODEL.EQ.7 .AND. IMODIN.EQ.8) THEN
          WRITE(LOUT,7011) 101, XAMIN( 1), 'cQ'
          WRITE(LOUT,7011) 102, XAMIN( 2), 'cD'
          WRITE(LOUT,7011) 103, XAMIN( 3), 'cU'
          WRITE(LOUT,7011) 104, XAMIN( 4), 'cL'
          WRITE(LOUT,7011) 105, XAMIN( 5), 'cE'
          WRITE(LOUT,7011) 106, XAMIN( 6), 'cHd'
          WRITE(LOUT,7011) 107, XAMIN( 7), 'cHu'
        END IF
C
CC          Unification results
C        WRITE(LOUT,1001) MGUTSS,GGUTSS,AGUTSS
C1001    FORMAT(/' ISASUGRA unification:'/' M_GUT      =',E10.3,
C     $  '   g_GUT          =',F5.3,3X,'   alpha_GUT =',F5.3)
C        WRITE(LOUT,999) FTGUT,FBGUT,FTAGUT
C999     FORMAT(' FT_GUT     =',F6.3,
C     $  '       FB_GUT         =',F6.3,3X,'  FL_GUT =',F6.3)
CC
C     Print inputs for GMSB models
      ELSE IF (IMODEL.EQ.2) THEN
        WRITE(LOUT,7011) 1, XGMIN(1), 'Lambda scale of soft SSB'
        WRITE(LOUT,7011) 2, XGMIN(2), 'M_mess overall messenger scale'
        WRITE(LOUT,7011) 3, XGMIN(4), 'tan(beta)'
        WRITE(LOUT,7011) 4, XGMIN(5), 'sign(mu)'
        WRITE(LOUT,7011) 5, XGMIN(3), 'N_5 messenger index'
        WRITE(LOUT,7011) 6, XGMIN(7), 'c_grav gravitino mass factor'
C        WRITE(LOUT,1002) (XGMIN(J),J=1,7)
C1002    FORMAT(
C     $  ' Lambda,  M_mes,  N_5,  tan(beta),  sgn(mu),  M_t,  C_grav='
C     $  /2E10.3,2F10.3,2X,F6.1,F10.3,1X,E10.3)
        WRITE(LOUT,7011)  51, XGMIN(12), 'N5_1  U(1)_Y messenger index'
        WRITE(LOUT,7011)  52, XGMIN(13), 'N5_2 SU(2)_L messenger index'
        WRITE(LOUT,7011)  53, XGMIN(14), 'N5_3 SU(3)_C messenger index'
        WRITE(LOUT,7011) 101, XGMIN( 8), 'Rsl'
        WRITE(LOUT,7011) 102, XGMIN( 9), 'dmH_d^2'
        WRITE(LOUT,7011) 103, XGMIN(10), 'dmH_u^2'
        WRITE(LOUT,7011) 104, XGMIN(11), 'd_Y'
C        WRITE(LOUT,1020) (XGMIN(J),J=8,14)
C1020    FORMAT(/' GMSB2 model input:'/
C     $  ' Rsl,    dmH_d^2,   dmH_u^2,     d_Y,     N5_1,  N5_2,  N5_3='
C     $  /F7.3,1X,E10.3,1X,E10.3,1X,E10.3,2X,3F7.3)
C        WRITE(LOUT,1003) AMGVSS
C1003    FORMAT(/' M(gravitino)=',E10.3)
      END IF
C
      Go to 1234
C     Weak scale couplings
      WRITE(LOUT,1004) ALEMI,SIN2W,AS
1004  FORMAT(/' 1/alpha_em =',F8.2,2X,
     $'   sin**2(thetaw) =',F6.4,2X,'   alpha_s   =',F5.3)
      WRITE(LOUT,1005) GSS(7),GSS(8),GSS(9)
1005  FORMAT(' M_1        =',F8.2,2X,
     $'   M_2            =',F8.2,'   M_3       =',F8.2)
      WRITE(LOUT,1006) MU,B,HIGFRZ
1006  FORMAT(' mu(Q)      =',F8.2,2X,
     $'   B(Q)           =',F8.2,'   Q         =',F8.2)
      WRITE(LOUT,1007) GSS(13),GSS(14)
1007  FORMAT(' M_H1^2     =',E10.3,'   M_H2^2         =',E10.3)
C
1234  Continue
C
C          Print mass spectrum from ISASUGRA
C
C     WRITE(LOUT,7000) ' '
C     WRITE(LOUT,6999) ' M_{GUT} =', MGUTSS
      WRITE(LOUT,7001) 'MASS', 'Scalar and gaugino mass spectrum'
      WRITE(LOUT,7000) ' PDG code   mass                 particle'
C
      If (Testing) then
        WRITE(LOUT,2000) MSS(1),MSS(2),MSS(3),MSS(4),MSS(5),MSS(10),
     $  MSS(11),MSS(12),MSS(13),MSS(14),MSS(17),MSS(18),MSS(16),
     $  MSS(21),MSS(22),MSS(23),MSS(24),MSS(25),MSS(26),MSS(27),
     $  MSS(28),MSS(29),MSS(30),MSS(31),MSS(32)
2000    FORMAT(/' ISAJET masses (with signs):'/
     $  ' M(GL)  =',F9.2/
     $  ' M(UL)  =',F9.2,'   M(UR)  =',F9.2,'   M(DL)  =',F9.2,
     $  '   M(DR) =',F9.2/
     $  ' M(B1)  =',F9.2,'   M(B2)  =',F9.2,'   M(T1)  =',F9.2,
     $  '   M(T2) =',F9.2/
     $  ' M(SN)  =',F9.2,'   M(EL)  =',F9.2,'   M(ER)  =',F9.2/
     $  ' M(NTAU)=',F9.2,'   M(TAU1)=',F9.2,'   M(TAU2)=',F9.2/
     $  ' M(Z1)  =',F9.2,'   M(Z2)  =',F9.2,'   M(Z3)  =',F9.2,
     $  '   M(Z4) =',F9.2/
     $  ' M(W1)  =',F9.2,'   M(W2)  =',F9.2/
     $  ' M(HL)  =',F9.2,'   M(HH)  =',F9.2,'   M(HA)  =',F9.2,
     $  '   M(H+) =',F9.2)
      EndIf
C
C     WRITE(LOUT,7013) iPDG(1),  MT, CHAF(1)
      WRITE(LOUT,7013)      24, AMW, ' W^+'
      DO 370 I=2,33
        sg = 1.
CsB     The signs of the (EW) gaugino masses are flipped according to ISAWIG
        If (iPDG(I).Eq.1000022 .or. iPDG(I).Eq.1000023 .or.
     .      iPDG(I).Eq.1000024 .or. iPDG(I).Eq.1000025 .or.
     .      iPDG(I).Eq.1000035 .or. iPDG(I).Eq.1000037) sg = -1.
        WRITE(LOUT,7013) iPDG(I), sg*MSS(ISA2LHA3(I)), CHAF(I)
 370  CONTINUE
C
C     SUSY scale
      RMSUSY = HIGFRZ !!! check this
C     WRITE(LOUT,7000) ' Higgs mixing'
      WRITE(LOUT,7001) 'ALPHA','Effective Higgs mixing parameter'
      WRITE(LOUT,7016) -ALFAH, 'alpha' ! Sign flips for LHA3
C
      If (Testing) then
        WRITE(LOUT,2001) THETAT,THETAB,THETAL,ALFAH
2001    FORMAT(/,' theta_t=',F9.4,'   theta_b=',F9.4,
     $  '   theta_l=',F9.4,'   alpha_h=',F9.4)
      EndIf
C
C     Write out chargino /neutralino masses/eigenvectors
C
      If (Testing) then
        WRITE(LOUT,3100) AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS
3100    FORMAT(/' NEUTRALINO MASSES (SIGNED) =',4F10.3)
        DO 100 J=1,4
          WRITE(LOUT,3200) J,(ZMIXSS(K,J),K=1,4)
3200      FORMAT(' EIGENVECTOR ',I1,'       =',4F10.5)
100     CONTINUE
        WRITE(LOUT,3300) AMW1SS,AMW2SS
3300    FORMAT(/' CHARGINO MASSES (SIGNED)  =',2F10.3)
        WRITE(LOUT,3400) GAMMAL,GAMMAR
3400    FORMAT(' GAMMAL, GAMMAR             =',2F10.5/)
      EndIf
C
CsB   Comparing Eq.(7.67a) in the book of Baer-Tata with Eq.(16) in
C     hep-ph/0311123 (the les Houches accord), knowing that the accord
C     (implicitly) uses the Haber-Kane translation between angles and
C     mixing matrices it's obvious that one has to flip the sign of the
C     sfernion mixing angles in ISAJET.
C     HB: I think instead we should dump out isajet mixing matrix since
C     SLHA doesn't ascribe to any particluar mixing scheme 10/28/05
CsB   Reverting sign of the mixing angles again on Apr 2009.
      WRITE(LOUT,7001) 'STOPMIX','stop mixing matrix'
      WRITE(LOUT,7021) 1, 1,  COS(THETAT), 'O_{11}'
      WRITE(LOUT,7021) 1, 2,  SIN(THETAT), 'O_{12}'
      WRITE(LOUT,7021) 2, 1, -SIN(THETAT), 'O_{21}'
      WRITE(LOUT,7021) 2, 2,  COS(THETAT), 'O_{22}'
      WRITE(LOUT,7001) 'SBOTMIX','sbottom mixing matrix'
      WRITE(LOUT,7021) 1, 1,  COS(THETAB), 'O_{11}'
      WRITE(LOUT,7021) 1, 2,  SIN(THETAB), 'O_{12}'
      WRITE(LOUT,7021) 2, 1, -SIN(THETAB), 'O_{21}'
      WRITE(LOUT,7021) 2, 2,  COS(THETAB), 'O_{22}'
      WRITE(LOUT,7001) 'STAUMIX','stau mixing matrix'
      WRITE(LOUT,7021) 1, 1,  COS(THETAL), 'O_{11}'
      WRITE(LOUT,7021) 1, 2,  SIN(THETAL), 'O_{12}'
      WRITE(LOUT,7021) 2, 1, -SIN(THETAL), 'O_{21}'
      WRITE(LOUT,7021) 2, 2,  COS(THETAL), 'O_{22}'
      WRITE(LOUT,7001) 'NMIX','neutralino mixing matrix'
CsB   For the gaugino mixing matrices I follow ISAWIG1200 to the letter
CsB   Ascending mass order (rows) and in the order
C     (bino, w3ino, higgs1, higgs2) (columns)
      DO I1=1,4
        DO I2=1,4
           sg = 1.
           If (I2.GT.2) sg = -1.
           J1 = 5 - I2
           J2 = I1
          WRITE(LOUT,7021) I1, I2, sg*ZMIXSS(J1,J2)
        EndDo
      EndDo
C
      THX=SIGN(1.,1./TAN(GAMMAL))
      THY=SIGN(1.,1./TAN(GAMMAR))
      WRITE(LOUT,7001) 'UMIX','chargino U mixing matrix'
      WRITE(LOUT,7021) 1, 1, -1.0*SIN(GAMMAL), 'U_{11}'
      WRITE(LOUT,7021) 1, 2, -1.0*COS(GAMMAL), 'U_{12}'
      WRITE(LOUT,7021) 2, 1, -THX*COS(GAMMAL), 'U_{21}'
      WRITE(LOUT,7021) 2, 2,  THX*SIN(GAMMAL), 'U_{22}'
      WRITE(LOUT,7001) 'VMIX','chargino V mixing matrix'
      WRITE(LOUT,7021) 1, 1, -1.0*SIN(GAMMAR), 'V_{11}'
      WRITE(LOUT,7021) 1, 2, -1.0*COS(GAMMAR), 'V_{12}'
      WRITE(LOUT,7021) 2, 1, -THY*COS(GAMMAR), 'V_{21}'
      WRITE(LOUT,7021) 2, 2,  THY*SIN(GAMMAR), 'V_{22}'
C
      WRITE(LOUT,7002) 'GAUGE',RMSUSY !!! check: are these at Q=RMSUSY?
      WRITE(LOUT,7011) 1, SQRT(.6)*GSS(1), 'g`'
      WRITE(LOUT,7011) 2, GSS(2), 'g_2'
      WRITE(LOUT,7011) 3, GSS(3), 'g_3'
C
      WRITE(LOUT,7002) 'YU',RMSUSY
      WRITE(LOUT,7021) 3, 3, GSS( 6), 'y_t' !!! check: are these at Q=RMSUSY?
C
      WRITE(LOUT,7002) 'YD',RMSUSY
      WRITE(LOUT,7021) 3, 3, GSS( 5), 'y_b'
C
      WRITE(LOUT,7002) 'YE',RMSUSY
      WRITE(LOUT,7021) 3, 3, GSS( 4), 'y_tau'
C
      WRITE(LOUT,7002) 'HMIX',RMSUSY,'Higgs mixing parameters' !!! check: are th
      WRITE(LOUT,7011) 1,           MU, 'mu(Q)'
      WRITE(LOUT,7011) 2,        XtanB, 'tan(beta)(M_GUT)'
      WRITE(LOUT,7011) 3, Sqrt(2.)*VEV, 'Higgs vev at Q'
      WRITE(LOUT,7011) 4,   MSS(31)**2, 'm_A^2(Q)'
C
      WRITE(LOUT,7002) 'MSOFT',RMSUSY,
     ,                 'DRbar SUSY breaking parameters' !!! check: are these at
      Do I=1,3
        If (I2L3GSS(I).NE.0)
     ,  WRITE(LOUT,7011) I,          GSS(I2L3GSS(I)),   SoftParaLHA(I)
      End Do
      Do I=4,50
        If (I2L3GSS(I).NE.0) !!! Fix sign - if necessary
     ,  WRITE(LOUT,7011) I, Sqrt(Abs(GSS(I2L3GSS(I)))), SoftParaLHA(I)
      End Do
C
      WRITE(LOUT,7002) 'AU',RMSUSY
      WRITE(LOUT,7021) 1, 1, GSS(12), 'A_u'
      WRITE(LOUT,7021) 2, 2, GSS(12), 'A_c'
      WRITE(LOUT,7021) 3, 3, GSS(12), 'A_t'
C
      WRITE(LOUT,7002) 'AD',RMSUSY
      WRITE(LOUT,7021) 1, 1, GSS(11), 'A_d'
      WRITE(LOUT,7021) 2, 2, GSS(11), 'A_s'
      WRITE(LOUT,7021) 3, 3, GSS(11), 'A_b'
C
      WRITE(LOUT,7002) 'AE',RMSUSY
      WRITE(LOUT,7021) 1, 1, GSS(10), 'A_e'
      WRITE(LOUT,7021) 2, 2, GSS(10), 'A_mu'
      WRITE(LOUT,7021) 3, 3, GSS(10), 'A_tau'
C
C          Print ISAJET MSSMi equivalent input
C
      If (Testing) then
        WRITE(LOUT,3000)
3000    FORMAT(/' ISAJET equivalent input:')
        WRITE(LOUT,3001) MSS(1),MU,MSS(31),XSUGIN(4)
3001    FORMAT(' MSSMA: ',4F8.2)
        WRITE(LOUT,3002) SQRT(GSS(19)),SQRT(GSS(17)),SQRT(GSS(18)),
     $  SQRT(GSS(16)),SQRT(GSS(15))
3002    FORMAT(' MSSMB: ',5F8.2)
        WRITE(LOUT,3003) SIGN(1.,GSS(24))*SQRT(ABS(GSS(24))),
     $  SQRT(GSS(22)),SIGN(1.,GSS(23))*SQRT(ABS(GSS(23))),
     $  SQRT(GSS(21)),SQRT(GSS(20)),GSS(12),GSS(11),GSS(10)
3003    FORMAT(' MSSMC: ',8F8.2)
        WRITE(LOUT,3004)
3004    FORMAT(' MSSMD: SAME AS MSSMB (DEFAULT)')
        WRITE(LOUT,3005) GSS(7),GSS(8)
3005    FORMAT(' MSSME: ',2F8.2)
      EndIf
C
      Close(91)
C
CsB LHA3 format statements
C
C     Formats for user information printout.
 5000 FORMAT(1x,17('*'),1x,'ISALHA v2.0: SUSY SPECTRUM '
     &     ,'INTERFACE',1x,17('*')/1x,'*',3x
     &     ,'ISALHA: Last Change',1x,A,1x,'-',1x,'C. Balazs')
 5001 FORMAT(1x,'*',3x,'Writing spectrum file on unit: ',I3)
 5002 FORMAT(1x,'*',3x,'Reading spectrum file on unit: ',I3)
 5003 FORMAT(1x,'*',3x,'Spectrum Calculator was: ',A,' version ',A)
 5100 FORMAT(1x,'*',1x,'Model parameters:'/1x,'*',1x,'----------------')
 5200 FORMAT(1x,'*',1x,3x,'m_0',6x,'m_{1/2}',5x,'A_0',3x,'tan(beta)',
     &     3x,'sgn(mu)',3x,'m_t'/1x,'*',1x,4(F8.2,1x),I8,2x,F8.2)
 5300 FORMAT(1x,'*'/1x,'*',1x,'Model spectrum :'/1x,'*',1x
     &     ,'----------------')
 5400 FORMAT(1x,'*',1x,A)
 5500 FORMAT(1x,'*',1x,A,':')
 5600 FORMAT(1x,'*',2x,2x,'M_GUT',2x,2x,'g_GUT',2x,1x,'alpha_GUT'/
     &       1x,'*',2x,1P,2(1x,E8.2),2x,E8.2)
 5700 FORMAT(1x,'*',4x,4x,'~d',2x,1x,4x,'~u',2x,1x,4x,'~s',2x,1x,
     &     4x,'~c',2x,1x,1x,'~b(12)',1x,1x,1x,'~t(12)'/1x,'*',2x,'L',1x
     &     ,6(F8.2,1x)/1x,'*',2x,'R',1x,6(F8.2,1x))
 5800 FORMAT(1x,'*'/1x,'*',4x,4x,'~e',2x,1x,2x,'~nu_e',2x,1x,3x,'~mu',2x
     &     ,1x,1x,'~nu_mu',1x,1x,'~tau(12)',1x,1x,'~nu_tau'/1x,'*',2x
     &     ,'L',1x,6(F8.2,1x)/1x,'*',2x,'R',1x,6(F8.2,1x))
 5900 FORMAT(1x,'*'/1x,'*',4x,4x,'~g',2x,1x,1x,'~chi_10',1x,1x,'~chi_20'
     &     ,1x,1x,'~chi_30',1x,1x,'~chi_40',1x,1x,'~chi_1+',1x
     &     ,1x,'~chi_2+'/1x,'*',3x,1x,7(F8.2,1x))
 6000 FORMAT(1x,'*'/1x,'*',4x,4x,'h0',2x,1x,4x,'H0',2x,1x,4x,'A0',2x
     &     ,1x,4x,'H+'/1x,'*',3x,1x,5(F8.2,1x))
 6100 FORMAT(1x,'*',11x,'|',3x,'~B',3x,'|',2x,'~W_3',2x,'|',2x
     &     ,'~H_1',2x,'|',2x,'~H_2',2x,'|'/1x,'*',3x,'~chi_10',1x,4('|'
     &     ,1x,F6.3,1x),'|'/1x,'*',3x,'~chi_20',1x,4('|'
     &     ,1x,F6.3,1x),'|'/1x,'*',3x,'~chi_30',1x,4('|'
     &     ,1x,F6.3,1x),'|'/1x,'*',3x,'~chi_40',1x,4('|'
     &     ,1x,F6.3,1x),'|')
 6200 FORMAT(1x,'*'/1x,'*',6x,'L',4x,'|',3x,'~W',3x,'|',3x,'~H',3x,'|'
     &     ,12x,'R',4x,'|',3x,'~W',3x,'|',3x,'~H',3x,'|'/1x,'*',3x
     &     ,'~chi_1+',1x,2('|',1x,F6.3,1x),'|',9x,'~chi_1+',1x,2('|',1x
     &     ,F6.3,1x),'|'/1x,'*',3x,'~chi_2+',1x,2('|',1x,F6.3,1x),'|',9x
     &     ,'~chi_2+',1x,2('|',1x,F6.3,1x),'|')
 6300 FORMAT(1x,'*'/1x,'*',8x,'|',2x,'~b_L',2x,'|',2x,'~b_R',2x,'|',8x
     &     ,'|',2x,'~t_L',2x,'|',2x,'~t_R',2x,'|',10x
     &     ,'|',1x,'~tau_L',1x,'|',1x,'~tau_R',1x,'|'/
     &     1x,'*',3x,'~b_1',1x,2('|',1x,F6.3,1x),'|',3x,'~t_1',1x,2('|'
     &     ,1x,F6.3,1x),'|',3x,'~tau_1',1x,2('|',1x,F6.3,1x),'|'/
     &     1x,'*',3x,'~b_2',1x,2('|',1x,F6.3,1x),'|',3x,'~t_2',1x,2('|'
     &     ,1x,F6.3,1x),'|',3x,'~tau_2',1x,2('|',1x,F6.3,1x),'|')
 6400 FORMAT(1x,'*',3x,'  A_b = ',F8.2,4x,'      A_t = ',F8.2,4x
     &     ,'A_tau = ',F8.2)
 6450 FORMAT(1x,'*',3x,'alpha = ',F8.2,4x,'tan(beta) = ',F8.2,4x
     &     ,'   mu = ',F8.2)
 6500 FORMAT(1x,32('*'),1x,'END OF ISALHA',1x,31('*'))
C
C     Format to use for comments
 6999 FORMAT('# ',A,1x,E16.8)
 7000 FORMAT('# ',A)
C     Format to use for block statements
 7001 FORMAT('Block',1x,A,3x,'#',1x,A)
 7002 FORMAT('Block',1x,A,1x,'Q=',1P,E16.8,0P,3x,'#',1x,A)
C     Indexed Int
 7010 FORMAT(1x,I5,1x,I5,3x,'#',1x,A)
C     Indexed Double
 7011 FORMAT(1x,I5,3x,1P,E16.8,0P,3x,'#',1x,A)
C     Indexed Char(12)
 7012 FORMAT(1x,I5,3x,A27,3x,'#',1x,A)
C     Long Indexed Double
 7013 FORMAT(1x,I9,3x,1P,E16.8,0P,3x,'#',1x,A)
C     Indexed Double without leading integer
 7016 FORMAT(8x,1P,E16.8,0P,3x,'#',1x,A)
C     Double Matrix
 7022 FORMAT(1x,I2,1x,I2,3x,1P,E16.8,3x,E16.8,0P,3x,'#',1x,A)
C     Single matrix
 7021 FORMAT(1x,I2,1x,I2,3x,1P,E16.8,0P,3x,'#',1x,A)
C     Write Decay Table
 7500 FORMAT('Decay',1x,I9,1x,'WIDTH=',1P,E16.8,0P,3x,'#',1x,A)
 7501 FORMAT(4x,1P,E16.8,0P,3x,I2,3x,'IDA=',1x,5(1x,I9),3x,'#',1x,A)
C
      RETURN
      END
