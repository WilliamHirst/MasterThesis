CDECK  ID>, DOMSSM.
      SUBROUTINE DOMSSM
C-----------------------------------------------------------------------
C          Initialize MSSM masses and decay modes from ISASUSY. 
C          Check for validity with ISAJET masses.
C          Decay modes are transfered to /DKYTAB/ by /SETDKY/.
C
C          F.E. Paige, November, 1992
C
C          Ver. 7.01: Add test so that AMASS is not called if ID = 0
C          Ver. 7.07: Add checking for LEP bounds.
C          Ver. 7.10: Add SUGRA interface
C          Ver. 7.32: Extend to large tanb solution
C          Ver. 7.33: Add gauge-mediated SUSY model
C          Ver. 7.38: NOGRAV turns off gravitino and weaker decays
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C          ISAJET common blocks
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/QLMASS/AMLEP(100),NQLEP,NMES,NBARY
      SAVE /QLMASS/
      INTEGER   NQLEP,NMES,NBARY
      REAL      AMLEP
      COMMON/XMSSM/GOMSSM,GOSUG,GOGMSB,GOAMSB,AL3UNI,GOMMAM,GOHCAM
     $,XGLSS,XMUSS,XHASS,XTBSS
     $,XQ1SS,XDRSS,XURSS,XL1SS,XERSS
     $,XQ2SS,XSRSS,XCRSS,XL2SS,XMRSS
     $,XQ3SS,XBRSS,XTRSS,XL3SS,XTARSS,XATSS,XABSS,XATASS
     $,XM1SS,XM2SS,XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU
     $,XLAMGM,XMESGM,XN5GM,XCMGV,XMGVTO
     $,XRSLGM,XDHDGM,XDHUGM,XDYGM,XN51GM,XN52GM,XN53GM
     $,XMN3NR,XMAJNR,XANSS,XNRSS,XSBCS,
     $XCQAM,XCDAM,XCUAM,XCLAM,XCEAM,XCHDAM,XCHUAM,
     $XL1AM,XL2AM,XL3AM
      SAVE /XMSSM/
      REAL XGLSS,XMUSS,XHASS,XTBSS
     $,XQ1SS,XDRSS,XURSS,XL1SS,XERSS
     $,XQ2SS,XSRSS,XCRSS,XL2SS,XMRSS
     $,XQ3SS,XBRSS,XTRSS,XL3SS,XTARSS,XATSS,XABSS,XATASS
     $,XM1SS,XM2SS
     $,XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU
     $,XLAMGM,XMESGM,XN5GM,XCMGV,XMGVTO
     $,XRSLGM,XDHDGM,XDHUGM,XDYGM,XN51GM,XN52GM,XN53GM
     $,XMN3NR,XMAJNR,XANSS,XNRSS,XSBCS,
     $XCQAM,XCDAM,XCUAM,XCLAM,XCEAM,XCHDAM,XCHUAM,
     $XL1AM,XL2AM,XL3AM
      LOGICAL GOMSSM,GOSUG,GOGMSB,GOAMSB,AL3UNI,GOMMAM,GOHCAM
      COMMON/NODCAY/NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOGRAV,
     $NOB,NOTAU
      LOGICAL NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOGRAV,
     $NOB,NOTAU
      SAVE /NODCAY/
C          ISASUSY common blocks
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
C          MXSS         =  maximum number of modes
C          NSSMOD       = number of modes
C          ISSMOD       = initial particle
C          JSSMOD       = final particles
C          GSSMOD       = width
C          BSSMOD       = branching ratio
C          MSSMOD       = decay matrix element pointer
C          LSSMOD       = logical flag used internally by SSME3
      INTEGER MXSS
      PARAMETER (MXSS=1000)
      COMMON/SSMODE/NSSMOD,ISSMOD(MXSS),JSSMOD(5,MXSS),GSSMOD(MXSS)
     $,BSSMOD(MXSS),MSSMOD(MXSS),LSSMOD
      INTEGER NSSMOD,ISSMOD,JSSMOD,MSSMOD
      REAL GSSMOD,BSSMOD
      LOGICAL LSSMOD
      SAVE /SSMODE/
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
C          SM ident code definitions. These are standard ISAJET but
C          can be changed.
      INTEGER IDUP,IDDN,IDST,IDCH,IDBT,IDTP
      INTEGER IDNE,IDE,IDNM,IDMU,IDNT,IDTAU
      INTEGER IDGL,IDGM,IDW,IDZ,IDH
      PARAMETER (IDUP=1,IDDN=2,IDST=3,IDCH=4,IDBT=5,IDTP=6)
      PARAMETER (IDNE=11,IDE=12,IDNM=13,IDMU=14,IDNT=15,IDTAU=16)
      PARAMETER (IDGL=9,IDGM=10,IDW=80,IDZ=90,IDH=81)
C          SUSY ident code definitions. They are chosen to be similar
C          to those in versions < 6.50 but may be changed.
      INTEGER ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1
      INTEGER ISNEL,ISEL,ISNML,ISMUL,ISNTL,ISTAU1
      INTEGER ISUPR,ISDNR,ISSTR,ISCHR,ISBT2,ISTP2
      INTEGER ISNER,ISER,ISNMR,ISMUR,ISNTR,ISTAU2
      INTEGER ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,ISGL
      INTEGER ISHL,ISHH,ISHA,ISHC
      INTEGER ISGRAV
      INTEGER IDTAUL,IDTAUR
      PARAMETER (ISUPL=21,ISDNL=22,ISSTL=23,ISCHL=24,ISBT1=25,ISTP1=26)
      PARAMETER (ISNEL=31,ISEL=32,ISNML=33,ISMUL=34,ISNTL=35,ISTAU1=36)
      PARAMETER (ISUPR=41,ISDNR=42,ISSTR=43,ISCHR=44,ISBT2=45,ISTP2=46)
      PARAMETER (ISNER=51,ISER=52,ISNMR=53,ISMUR=54,ISNTR=55,ISTAU2=56)
      PARAMETER (ISGL=29)
      PARAMETER (ISZ1=30,ISZ2=40,ISZ3=50,ISZ4=60,ISW1=39,ISW2=49)
      PARAMETER (ISHL=82,ISHH=83,ISHA=84,ISHC=86)
      PARAMETER (ISGRAV=91)
      PARAMETER (IDTAUL=10016,IDTAUR=20016)
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
C
      INTEGER NOUT
      PARAMETER (NOUT=33)
      INTEGER IDOUT(NOUT)
      REAL AMASS,AMPL
      REAL AMI,SUMGAM,SUMMJ,WIDMX
      REAL QSUSY,ASMB,MBMB,ASMT,MTMT,SUALFS,PI,GG
      DOUBLE PRECISION SSMQCD
      INTEGER I,J,K,IFL1,IFL2,IFL3,JSPIN,INDEX,IALLOW,IITEST,IMDL
      INTEGER IMHL,IMHC
C
      DATA IDOUT/
     $IDTP,ISGL,ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1,ISUPR,ISDNR,
     $ISSTR,ISCHR,ISBT2,ISTP2,ISEL,ISMUL,ISTAU1,ISNEL,ISNML,ISNTL,
     $ISER,ISMUR,ISTAU2,ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,
     $ISHL,ISHH,ISHA,ISHC/
      DATA AMPL/2.4E18/,IAL3UN/0/
C
C          Generate masses and decays
C
C     FIRST SET HIGH SCALE FOR SUSY BCs; default is M_GUT
      XSUGIN(7)=XSBCS
      IF (XMGVTO.LT.1.E19) AMGVSS=XMGVTO
      IF(.NOT.GOMSSM) RETURN
      LOUT=ITLIS
      IF (AL3UNI) IAL3UN=1
      IF (INUHM.EQ.1) THEN
        MU=XNUSUG(19)
        AMHA=XNUSUG(20)
        TWOM1=-MU
      END IF
      IF(GOSUG) THEN
C          SUGRA input
C          First solve renormalization group equations
        IF (XMAJNR.LT.1.E19) THEN
          XNRIN(1)=XMN3NR
          XNRIN(2)=XMAJNR
          XNRIN(3)=XANSS
          XNRIN(4)=XNRSS
        ELSE
          XNRIN(2)=1.E20
        END IF
        IF (GOAMSB.OR.GOMMAM.OR.GOHCAM) THEN
          XA0SU=0.
          XAMIN(1)=XCQAM
          XAMIN(2)=XCDAM
          XAMIN(3)=XCUAM
          XAMIN(4)=XCLAM
          XAMIN(5)=XCEAM
          XAMIN(6)=XCHDAM
          XAMIN(7)=XCHUAM
          XAMIN(8)=XL1AM
          XAMIN(9)=XL2AM
          XAMIN(10)=XL3AM
          XAMIN(11)=XM0SU
          IF (GOAMSB) THEN
            IMDL=7
            CALL SUGRA(XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU,AMASS(6),IMDL)
          ELSE IF (GOMMAM) THEN
            IMDL=9
            CALL SUGRA(XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU,AMASS(6),IMDL)
          ELSE IF (GOHCAM) THEN
            IMDL=10
            CALL SUGRA(0.,XMHSU,XA0SU,XTGBSU,XSMUSU,AMASS(6),IMDL)
          END IF
        ELSE
          IMDL=1
          CALL SUGRA(XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU,AMASS(6),IMDL)
        END IF
        IF (NOGOOD.EQ.1) THEN
          WRITE(LOUT,*) 'SUGRA BAD POINT: TACHYONIC PARTICLES!'
        ELSE IF (NOGOOD.EQ.2) THEN
          WRITE(LOUT,*) 'SUGRA BAD POINT: NO EW SYMMETRY BREAKING!'
        ELSE IF (NOGOOD.EQ.3) THEN
          WRITE(LOUT,*) 'SUGRA BAD POINT: M(H_P)^2<0!'
        ELSE IF (NOGOOD.EQ.4) THEN
          WRITE(LOUT,*) 'SUGRA BAD POINT: YUKAWA>10!'
        ELSE IF (NOGOOD.EQ.5) THEN
          WRITE(LOUT,*) 'SUGRA BAD POINT: Z1SS NOT LSP!'
        ELSE IF (NOGOOD.EQ.7) THEN
          WRITE(LOUT,*) 'SUGRA BAD POINT: XT EWSB IS BAD!'
        ELSE IF (NOGOOD.EQ.8) THEN
          WRITE(LOUT,*) 'SUGRA BAD POINT: MHL^2<0!'
        END IF
        IF (MHPNEG.EQ.1) THEN
          WRITE(LOUT,*) 'SUGRA BAD POINT: M(H_P)^2<0!!'
          NOGOOD=3
        END IF
        IF(NOGOOD.NE.0) STOP99
        IF(ITACHY.NE.0) THEN
          WRITE(LOUT,*) 'WARNING: TACHYONIC SLEPTONS AT GUT SCALE'
          WRITE(LOUT,*) '         POINT MAY BE INVALID'
        ENDIF
C          Then calculate masses and decays
        CALL SSMSSM(XISAIN(1),XISAIN(2),XISAIN(3),
     $ XISAIN(4),XISAIN(5),XISAIN(6),XISAIN(7),XISAIN(8),XISAIN(9),
     $ XISAIN(10),XISAIN(11),XISAIN(12),XISAIN(13),XISAIN(14),
     $ XISAIN(15),XISAIN(16),XISAIN(17),XISAIN(18),XISAIN(19),
     $ XISAIN(20),XISAIN(21),XISAIN(22),XISAIN(23),XISAIN(24),
     $ AMASS(6),IALLOW,IMDL,IMHL,IMHC)
      ELSE IF(GOGMSB) THEN
C          GMSB input
        XGMIN(8)=XRSLGM
        XGMIN(9)=XDHDGM
        XGMIN(10)=XDHUGM
        XGMIN(11)=XDYGM
        XGMIN(12)=XN51GM
        XGMIN(13)=XN52GM
        XGMIN(14)=XN53GM
C          First solve renormalization group equations
        IMDL=2
        CALL SUGRA(XLAMGM,XMESGM,XN5GM,XTGBSU,XSMUSU,AMASS(6),IMDL)
        IF (NOGOOD.EQ.1) THEN
          WRITE(LOUT,*) 'GMSB BAD POINT: TACHYONIC PARTICLES!'
        ELSE IF (NOGOOD.EQ.2) THEN
          WRITE(LOUT,*) 'GMSB BAD POINT: NO EW SYMMETRY BREAKING!'
        ELSE IF (NOGOOD.EQ.3) THEN
          WRITE(LOUT,*) 'GMSB BAD POINT: M(H_P)^2<0!'
        ELSE IF (NOGOOD.EQ.4) THEN
          WRITE(LOUT,*) 'GMSB BAD POINT: YUKAWA>100!'
        ELSE IF (NOGOOD.EQ.7) THEN
          WRITE(LOUT,*) 'GMSB BAD POINT: XT EWSB IS BAD!'
        ELSE IF (NOGOOD.EQ.8) THEN
          WRITE(LOUT,*) 'GMSB BAD POINT: MHL^2<0!'
        END IF
        IF (MHPNEG.EQ.1) THEN
          WRITE(LOUT,*) 'GMSB BAD POINT: M(H_P)^2<0!!'
          NOGOOD=3
        END IF
        IF(NOGOOD.NE.0) STOP99
        IF(ITACHY.NE.0) THEN
          WRITE(LOUT,*) 'WARNING: TACHYONIC SLEPTONS AT HIGH SCALE'
          WRITE(LOUT,*) '         POINT MAY BE INVALID'
        ENDIF
C          Then calculate masses and decays
        AMGVSS=XLAMGM*XMESGM*XCMGV/SQRT(3.)/AMPL
        CALL SSMSSM(XISAIN(1),XISAIN(2),XISAIN(3),
     $  XISAIN(4),XISAIN(5),XISAIN(6),XISAIN(7),XISAIN(8),XISAIN(9),
     $  XISAIN(10),XISAIN(11),XISAIN(12),XISAIN(13),XISAIN(14),
     $  XISAIN(15),XISAIN(16),XISAIN(17),XISAIN(18),XISAIN(19),
     $  XISAIN(20),XISAIN(21),XISAIN(22),XISAIN(23),XISAIN(24),
     $  AMASS(6),IALLOW,IMDL,IMHL,IMHC)
      ELSE
C          Weak scale input
C          Values of 1.E20 indicate that SSMASS should calculate
C          M_1 and M_2 from M_3
C          First do fermion masses at QSUSY since SUGRA is not called
        QSUSY=SQRT(XQ3SS*XTRSS)
        PI=4.*ATAN(1.)
C          Define heavy quark pole masses and LambdaQCD:
        AMBT=AMASS(5)
        AMTP=AMASS(6)
        ALQCD4=0.177
        ASMB=SUALFS(AMBT**2,.36,AMTP,3)
        MBMB=AMBT*(1.-4*ASMB/3./PI)
        MBQ=SSMQCD(DBLE(MBMB),DBLE(QSUSY))
        ASMT=SUALFS(AMTP**2,.36,AMTP,3)
        MTMT=AMTP/(1.+4*ASMT/3./PI+(16.11-1.04*(5.-6.63/AMTP))*
     $  (ASMT/PI)**2)
        MTQ=SSMQCD(DBLE(MTMT),DBLE(QSUSY))
        MLQ=1.7463
C       Define TANBQ parameters= TANB for MSSM runs, but not for SUGRA
        AMW=80.423
        ALFAEM=1./128.
        SN2THW=.232
        GG=SQRT(4*PI*ALFAEM/SN2THW)
        VUQ=SQRT(2*AMW**2/GG**2/(1.+1./XTBSS**2))
        VDQ=VUQ/XTBSS
        CALL SSMSSM(XGLSS,XMUSS,XHASS,XTBSS,XQ1SS,XDRSS,XURSS,XL1SS,
     $  XERSS,XQ2SS,XSRSS,XCRSS,XL2SS,XMRSS,XQ3SS,XBRSS,XTRSS,XL3SS,
     $  XTARSS,XATSS,XABSS,XATASS,XM1SS,XM2SS,AMASS(6),IALLOW,0,
     $IMHL,IMHC)
      ENDIF
C
C          Test parameters
C
      IF(IALLOW.NE.0) THEN
        WRITE(LOUT,1000)
1000    FORMAT(//' MSSM WARNING: Z1SS IS NOT LSP')
      ENDIF
      CALL SSTEST(IALLOW)
      IITEST=IALLOW/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,1002)
1002    FORMAT(' MSSM WARNING: Z -> Z1SS Z1SS TOO BIG')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,1004)
1004    FORMAT(' MSSM WARNING: Z -> CHARGINOS ALLOWED')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,1008)
1008    FORMAT(' MSSM WARNING: Z -> Z1SS Z2SS TOO BIG')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,1008)
1016    FORMAT(' MSSM WARNING: Z -> SQUARKS OR SLEPTONS')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,1032)
1032    FORMAT(' MSSM WARNING: Z -> Z* HL0 TOO BIG')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,1064)
1064    FORMAT(' MSSM WARNING: Z -> HL0 HA0 ALLOWED')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,1128)
1128    FORMAT(' MSSM WARNING: Z -> H+ H- ALLOWED')
      ENDIF
C
C          Store masses in /QLMASS/
C
      CALL FLAVOR(ISUPL,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMULSS
      CALL FLAVOR(ISDNL,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMDLSS
      CALL FLAVOR(ISSTL,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMSLSS
      CALL FLAVOR(ISCHL,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMCLSS
      CALL FLAVOR(ISBT1,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMB1SS
      CALL FLAVOR(ISTP1,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMT1SS
      CALL FLAVOR(ISUPR,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMURSS
      CALL FLAVOR(ISDNR,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMDRSS
      CALL FLAVOR(ISSTR,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMSRSS
      CALL FLAVOR(ISCHR,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMCRSS
      CALL FLAVOR(ISBT2,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMB2SS
      CALL FLAVOR(ISTP2,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMT2SS
C
      CALL FLAVOR(ISNEL,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMN1SS
      CALL FLAVOR(ISEL,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMELSS
      CALL FLAVOR(ISNML,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMN2SS
      CALL FLAVOR(ISMUL,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMMLSS
      CALL FLAVOR(ISNTL,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMN3SS
      CALL FLAVOR(ISTAU1,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AML1SS
      CALL FLAVOR(ISER,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMERSS
      CALL FLAVOR(ISMUR,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AMMRSS
      CALL FLAVOR(ISTAU2,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=AML2SS
C
      CALL FLAVOR(ISGL,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=ABS(AMGLSS)
      CALL FLAVOR(ISZ1,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=ABS(AMZ1SS)
      CALL FLAVOR(ISZ2,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=ABS(AMZ2SS)
      CALL FLAVOR(ISZ3,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=ABS(AMZ3SS)
      CALL FLAVOR(ISZ4,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=ABS(AMZ4SS)
      CALL FLAVOR(ISW1,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=ABS(AMW1SS)
      CALL FLAVOR(ISW2,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=ABS(AMW2SS)
C
      CALL FLAVOR(ISHL,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=ABS(AMHL)
      CALL FLAVOR(ISHH,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=ABS(AMHH)
      CALL FLAVOR(ISHA,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=ABS(AMHA)
      CALL FLAVOR(ISHC,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=ABS(AMHC)
C
C          Check decays with ISAJET masses
C          NOGRAV turns off gravitino decays and all weaker ones
C
      WIDMX=0
      IF(NOGRAV) THEN
        DO 90 J=1,NSSMOD
          DO 91 K=1,5
            IF(JSSMOD(K,J).EQ.ISGRAV) WIDMX=MAX(WIDMX,GSSMOD(J))
91        CONTINUE
90      CONTINUE
      ENDIF
      WIDMX=1.01*WIDMX
C
      DO 100 I=1,NOUT
        SUMGAM=0
        AMI=AMASS(IDOUT(I))
        DO 110 J=1,NSSMOD
          IF(IDOUT(I).NE.ISSMOD(J)) GO TO 110
          SUMMJ=0
          DO 111 K=1,5
            IF(JSSMOD(K,J).NE.0) SUMMJ=SUMMJ+AMASS(JSSMOD(K,J))
111       CONTINUE
          IF(SUMMJ.GE.AMI.OR.GSSMOD(J).LT.WIDMX) GSSMOD(J)=0
          SUMGAM=SUMGAM+GSSMOD(J)
110     CONTINUE
        DO 120 J=1,NSSMOD
          IF(IDOUT(I).NE.ISSMOD(J)) GO TO 120
          IF(SUMGAM.NE.0) THEN
            BSSMOD(J)=GSSMOD(J)/SUMGAM
          ELSE
            BSSMOD(J)=0
          ENDIF
120     CONTINUE
100   CONTINUE
C
      RETURN
      END
