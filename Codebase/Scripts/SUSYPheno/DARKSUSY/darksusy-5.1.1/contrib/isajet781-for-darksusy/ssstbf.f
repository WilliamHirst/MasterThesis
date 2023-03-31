CDECK  ID>, SSSTBF.
        SUBROUTINE SSSTBF
C-----------------------------------------------------------------------
C
C        This program gives stop squark branching fractions to gauginos
C        according to Baer and Tata.
C        If no other modes are allowed, stop -> c z_i through loops is
C        used as the default.
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
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
C
      COMPLEX ZI,ZONE,ZA,ZB,ZPP,ZPM,ZAUIZ,ZBUIZ
      DOUBLE PRECISION SSALFS,DWID
      REAL SSXLAM
      REAL WID,AWD(2),BW(2),FB,FT,XM,YM,THX,THY,AU1,MZ1,WIDC1
      REAL PI,SR2,G,GP,TANB,COTB,MPL,MMI,AH
      REAL AUIZ,MZIZ,SINT,COST,AS,BS,SNZI,THIZ
      INTEGER IZ,ISTOP,IDSTOP
      REAL AMSTOP,BWP(2),A
      REAL MW1,MW2,SNW1,SNW2,CS2THW,BETA,TN2THW,SINB,COSB
      REAL EPSILON,DELTAL,DELTAR,MSUSY,KMTB,KMCB,COS2B
      INTEGER ISZIZ(4)
      DATA ZONE/(1.,0.)/,ZI/(0.,1.)/
C
C          Partly duplicated from SSMASS
C
      CS2THW=1.-SN2THW
      TN2THW=SN2THW/CS2THW
      PI=4.*ATAN(1.)
      SR2=SQRT(2.)
      G=SQRT(4*PI*ALFAEM/SN2THW)
      GP=G*SQRT(SN2THW/(1.-SN2THW))
      TANB=1./RV2V1
      COTB=RV2V1
      BETA=ATAN(TANB)
      MSUSY=SQRT(MAX(AMZ**2,AMTLSS*AMTRSS*SIGN(1.,AMTLSS*AMTRSS)))
      KMTB=0.9991
      KMCB=0.0413
C          Reconstruct masses from SSMASS
      FB=G*MBQ/SR2/AMW/COS(BETA)
      FT=G*MTQ/SR2/AMW/SIN(BETA)
      MW1=ABS(AMW1SS)
      MW2=ABS(AMW2SS)
      SNW1=SIGN(1.,AMW1SS)
      SNW2=SIGN(1.,AMW2SS)
      XM=1./TAN(GAMMAL)
      YM=1./TAN(GAMMAR)
      THX=SIGN(1.,XM)
      THY=SIGN(1.,YM)
C
      AWD(1)=-G*SNW1*SIN(GAMMAR)
      AWD(2)=-G*SNW2*THY*COS(GAMMAR)
      BW(1)=-FT*SNW1*COS(GAMMAR)
      BW(2)=FT*SNW2*THY*SIN(GAMMAR)
      BWP(1)=-FB*COS(GAMMAL)
      BWP(2)=FB*THX*SIN(GAMMAL)
      MMI=AMW1SS
      MPL=AMW2SS
      COST=COS(THETAT)
      SINT=SIN(THETAT)
      COSB=COS(THETAB)
      SINB=SIN(THETAB)
      COS2B=COS(2*THETAB)
C
C          Compute stop_i branching fractions to charm + zi if no other
C          modes are allowed. WIDC1 is taken from Hikasa-Kobayashi
C          and subsequent decays to charm are scaled in terms of this.
c          WIDC1 updated from Andrew Box 6/14/07 for GOSUG=.true.
C
      ISZIZ(1)=ISZ1
      ISZIZ(2)=ISZ2
      ISZIZ(3)=ISZ3
      ISZIZ(4)=ISZ4
      AU1=-G/SR2*ZMIXSS(3,1)-GP/3./SR2*ZMIXSS(4,1)
      MZ1=ABS(AMZISS(1))
      DO 100 ISTOP=1,2
        IF(ISTOP.EQ.1) THEN
          AMSTOP=AMT1SS
          IDSTOP=ISTP1
        ELSE
          AMSTOP=AMT2SS
          IDSTOP=ISTP2
        ENDIF
        IF(AMSTOP.LT.(MW1+AMBT).AND.AMSTOP.GT.(AMCH+MZ1)) THEN
          IF (GORGE) THEN
            DELTAL=-LOG(MGUTSS/MSUSY)/16./PI**2*KMTB*KMCB*FB**2
     $         *(GSS(19)+GSS(24)+2*GSS(13)+2*GSS(22)
     $         +2.D0*GSS(11)**2)
            DELTAR=LOG(MGUTSS/MSUSY)/16./PI**2*KMTB*KMCB*FB**2
     $         *MTQ*2.D0*GSS(11)
            EPSILON=(DELTAL*COST-DELTAR*SINT)/(AMT1SS**2-AMCLSS**2)
            WIDC1=EPSILON**2/16./PI*AMT1SS
     $        *(1.-(MZ1**2/AMT1SS**2))**2*AU1**2
          ELSE      
            WIDC1=3.E-10*AMT1SS*(1.-MZ1**2/AMT1SS**2)**2
          END IF
          DO 110 IZ=1,4
            MZIZ=ABS(AMZISS(IZ))
            AUIZ=-G/SR2*ZMIXSS(3,IZ)-GP/3./SR2*ZMIXSS(4,IZ)
            IF (AMT1SS.GT.(AMCH+MZIZ)) THEN
              WID=AUIZ**2*(AMSTOP**2-MZIZ**2)/AU1**2
     $        /(AMSTOP**2-MZ1**2)*WIDC1
              CALL SSSAVE(IDSTOP,WID,ISZIZ(IZ),IDCH,0,0,0)
            END IF
110       CONTINUE
        ELSEIF(AMSTOP.LT.(MW1+AMBT).AND.AMSTOP.LE.(AMCH+MZ1)) THEN
          WRITE(LOUT,1000) ISTOP
1000      FORMAT(' ERROR IN SSSTBF: NO ALLOWED MODE FOR STOP',I2)
        END IF
100   CONTINUE
C
C          stop_i -> gluino + top
C
      IF (AMT1SS.GT.(AMGLSS+AMTP)) THEN
        WID=2*SSALFS(DBLE(AMT1SS**2))*AMT1SS*((1.-AMGLSS**2/AMT1SS**2-
     $  AMTP**2/AMT1SS**2)-SGNM3*2*SIN(2*THETAT)*AMTP*AMGLSS/AMT1SS**2)
     $  *SQRT(SSXLAM(1.,AMGLSS**2/AMT1SS**2,AMTP**2/AMT1SS**2))/3.
        CALL SSSAVE(ISTP1,WID,ISGL,IDTP,0,0,0)
      END IF
C
      IF (AMT2SS.GT.(AMGLSS+AMTP)) THEN
        WID=2*SSALFS(DBLE(AMT2SS**2))*AMT2SS*((1.-AMGLSS**2/AMT2SS**2-
     $  AMTP**2/AMT2SS**2)+SGNM3*2*SIN(2*THETAT)*AMTP*AMGLSS/AMT2SS**2)
     $  *SQRT(SSXLAM(1.,AMGLSS**2/AMT2SS**2,AMTP**2/AMT2SS**2))/3.
        CALL SSSAVE(ISTP2,WID,ISGL,IDTP,0,0,0)
      END IF
C
C          stop_1 -> top + zino_i
C
      DO 200 IZ=1,4
        MZIZ=ABS(AMZISS(IZ))
        SNZI=SIGN(1.,AMZISS(IZ))
        IF (SNZI.EQ.1.) THEN
           THIZ=0.
        ELSE
           THIZ=1.
        END IF
        ZAUIZ=ZI**(THIZ-1.)*SNZI
     $  *(-G/SR2*ZMIXSS(3,IZ)-GP/3./SR2*ZMIXSS(4,IZ))
        ZBUIZ=ZI**(THIZ-1.)*4*GP*ZMIXSS(4,IZ)/3./SR2
        ZPP=ZI**THIZ
        ZPM=(-ZI)**THIZ
        ZA=((ZI*ZAUIZ-ZPP*FT*ZMIXSS(1,IZ))*COST
     $  -(ZI*ZBUIZ-ZPM*FT*ZMIXSS(1,IZ))*SINT)/2.
        ZB=((-ZI*ZAUIZ-ZPP*FT*ZMIXSS(1,IZ))*COST
     $  -(ZI*ZBUIZ+ZPM*FT*ZMIXSS(1,IZ))*SINT)/2.
        AS=ZA*CONJG(ZA)
        BS=ZB*CONJG(ZB)
        IF (AMT1SS.GT.(AMTP+MZIZ)) THEN
          WID=(AS*(AMT1SS**2-(AMTP+MZIZ)**2)+BS
     $    *(AMT1SS**2-(AMTP-MZIZ)**2))/8./PI/AMT1SS
     $    *SQRT(SSXLAM(1.,AMTP**2/AMT1SS**2,MZIZ**2/AMT1SS**2))
          CALL SSSAVE(ISTP1,WID,ISZIZ(IZ),IDTP,0,0,0)
        END IF
200   CONTINUE
C
C          Wino decays
C
      IF (AMT1SS.GT.(AMBT+MW1)) THEN
        A=AWD(1)*COST-BW(1)*SINT
        AS=A*A
        WID=AMT1SS*((AS+BWP(1)**2*COST**2)*(1.-MW1**2/AMT1SS**2-
     $   AMBT**2/AMT1SS**2)-4*MW1*AMBT*BWP(1)*COST*A/AMT1SS**2)
     $   *SQRT(SSXLAM(1.,MW1**2/AMT1SS**2,AMBT**2/AMT1SS**2))/16./PI
        CALL SSSAVE(ISTP1,WID,ISW1,IDBT,0,0,0)
      END IF
      IF (AMT1SS.GT.(AMBT+MW2)) THEN
        A=AWD(2)*COST-BW(2)*SINT
        AS=A*A
        WID=AMT1SS*((AS+BWP(2)**2*COST**2)*(1.-MW2**2/AMT1SS**2-
     $   AMBT**2/AMT1SS**2)-4*MW2*AMBT*BWP(2)*COST*A/AMT1SS**2)
     $   *SQRT(SSXLAM(1.,MW2**2/AMT1SS**2,AMBT**2/AMT1SS**2))/16./PI
        CALL SSSAVE(ISTP1,WID,ISW2,IDBT,0,0,0)
      END IF
C 
      IF (AMT2SS.GT.(AMBT+MW1)) THEN
        A=AWD(1)*SINT+BW(1)*COST
        AS=A*A
        WID=AMT2SS*((AS+BWP(1)**2*SINT**2)*(1.-MW1**2/AMT2SS**2-
     $   AMBT**2/AMT2SS**2)-4*MW1*AMBT*BWP(1)*SINT*A/AMT2SS**2)
     $   *SQRT(SSXLAM(1.,MW1**2/AMT2SS**2,AMBT**2/AMT2SS**2))/16./PI
        CALL SSSAVE(ISTP2,WID,ISW1,IDBT,0,0,0)
      END IF
      IF (AMT2SS.GT.(AMBT+MW2)) THEN
        A=AWD(2)*SINT+BW(2)*COST
        AS=A*A
        WID=AMT2SS*((AS+BWP(2)**2*SINT**2)*(1.-MW2**2/AMT2SS**2-
     $   AMBT**2/AMT2SS**2)-4*MW2*AMBT*BWP(2)*SINT*A/AMT2SS**2)
     $  *SQRT(SSXLAM(1.,MW2**2/AMT2SS**2,AMBT**2/AMT2SS**2))/16./PI
        CALL SSSAVE(ISTP2,WID,ISW2,IDBT,0,0,0)
      END IF
C
C          stop_2 -> stop_1 + X modes
C
      IF (AMT2SS.GT.(AMT1SS+AMZ)) THEN
        WID=G**2*COST**2*SINT**2
     $  *(SQRT(SSXLAM(AMT2SS**2,AMZ**2,AMT1SS**2)))**3
     $  /64./PI/CS2THW/AMT2SS**3/AMZ**2
        CALL SSSAVE(ISTP2,WID,IDZ,ISTP1,0,0,0)
      END IF
C
      IF (AMT2SS.GT.(AMT1SS+AMHL)) THEN
        AH=G*AMW*SIN(BETA-ALFAH)*(1.-5.*TN2THW/3.)*SINT*COST/2.
     $  +G*AMTP*COS(2.*THETAT)*(TWOM1*SIN(ALFAH)+AAT*COS(ALFAH))/2.
     $  /AMW/SIN(BETA)
        WID=AH**2/16./PI/AMT2SS**3
     $  *SQRT(SSXLAM(AMT2SS**2,AMHL**2,AMT1SS**2))
        CALL SSSAVE(ISTP2,WID,ISHL,ISTP1,0,0,0)
      END IF
C
      IF (AMT2SS.GT.(AMT1SS+AMHH)) THEN
        AH=-G*AMW*COS(BETA-ALFAH)*(1.-5.*TN2THW/3.)*SINT*COST/2.
     $  +G*AMTP*COS(2.*THETAT)*(TWOM1*COS(ALFAH)-AAT*SIN(ALFAH))/2.
     $  /AMW/SIN(BETA)
        WID=AH**2/16./PI/AMT2SS**3
     $  *SQRT(SSXLAM(AMT2SS**2,AMHH**2,AMT1SS**2))
        CALL SSSAVE(ISTP2,WID,ISHH,ISTP1,0,0,0)
      END IF
C
      IF (AMT2SS.GT.(AMT1SS+AMHA)) THEN
        AH=G*AMTP*(TWOM1-AAT/TANB)/2./AMW
        WID=AH**2/16./PI/AMT2SS**3
     $  *SQRT(SSXLAM(AMT2SS**2,AMHA**2,AMT1SS**2))
        CALL SSSAVE(ISTP2,WID,ISHA,ISTP1,0,0,0)
      END IF
C
C          t_i --> b_i + W decays
C
      IF (AMT1SS.GT.(AMB1SS+AMW)) THEN
        WID=G**2*COST**2*COSB**2*(SSXLAM(AMT1SS**2,AMB1SS**2,
     $AMW**2))**1.5/32./PI/AMT1SS**3/AMW**2
        CALL SSSAVE(ISTP1,WID,IDW,ISBT1,0,0,0)
      END IF
C
      IF (AMT1SS.GT.(AMB2SS+AMW)) THEN
        WID=G**2*COST**2*SINB**2*(SSXLAM(AMT1SS**2,AMB2SS**2,
     $AMW**2))**1.5/32./PI/AMT1SS**3/AMW**2
        CALL SSSAVE(ISTP1,WID,IDW,ISBT2,0,0,0)
      END IF
C
      IF (AMT2SS.GT.(AMB1SS+AMW)) THEN
        WID=G**2*SINT**2*COSB**2*(SSXLAM(AMT2SS**2,AMB1SS**2,
     $AMW**2))**1.5/32./PI/AMT2SS**3/AMW**2
        CALL SSSAVE(ISTP2,WID,IDW,ISBT1,0,0,0)
      END IF
C
      IF (AMT2SS.GT.(AMB2SS+AMW)) THEN
        WID=G**2*SINT**2*SINB**2*(SSXLAM(AMT2SS**2,AMB2SS**2,
     $AMW**2))**1.5/32./PI/AMT2SS**3/AMW**2
        CALL SSSAVE(ISTP2,WID,IDW,ISBT2,0,0,0)
      END IF
C
C          t_i --> b_i + H+ decays
C
      IF (AMT1SS.GT.(AMB1SS+AMHC)) THEN
        A=G/SR2/AMW*(AMTP*AMBT*(COTB+TANB)*SINT*SINB+
     $(AMBT**2*TANB+AMTP**2*COTB-AMW**2*SIN(2*BETA))*
     $COST*COSB-AMTP*(TWOM1-AAT*COTB)*SINT*COSB-AMBT*
     $(TWOM1-AAB*TANB)*SINB*COST)
        AS=A*A
        WID=AS*SQRT(SSXLAM(AMT1SS**2,AMB1SS**2,AMHC**2))/
     $      16./PI/AMT1SS**3
        CALL SSSAVE(ISTP1,WID,ISHC,ISBT1,0,0,0)
      END IF
C
      IF (AMT1SS.GT.(AMB2SS+AMHC)) THEN
        A=G/SR2/AMW*(-AMTP*AMBT*(COTB+TANB)*SINT*COSB+
     $(AMBT**2*TANB+AMTP**2*COTB-AMW**2*SIN(2*BETA))*
     $COST*SINB-AMTP*(TWOM1-AAT*COTB)*SINT*SINB+AMBT*
     $(TWOM1-AAB*TANB)*COST*COSB)
        AS=A*A
        WID=AS*SQRT(SSXLAM(AMT1SS**2,AMB2SS**2,AMHC**2))/
     $      16./PI/AMT1SS**3
        CALL SSSAVE(ISTP1,WID,ISHC,ISBT2,0,0,0)
      END IF
C
      IF (AMT2SS.GT.(AMB1SS+AMHC)) THEN
        A=G/SR2/AMW*(-AMTP*AMBT*(COTB+TANB)*COST*SINT+
     $(AMBT**2*TANB+AMTP**2*COTB-AMW**2*SIN(2*BETA))*
     $SINT*COSB+AMTP*(TWOM1-AAT*COTB)*COST*COSB-AMBT*
     $(TWOM1-AAB*TANB)*SINT*SINB)
        AS=A*A
        WID=AS*SQRT(SSXLAM(AMT2SS**2,AMB1SS**2,AMHC**2))/
     $      16./PI/AMT2SS**3
        CALL SSSAVE(ISTP2,WID,ISHC,ISBT1,0,0,0)
      END IF
C
      IF (AMT2SS.GT.(AMB2SS+AMHC)) THEN
        A=G/SR2/AMW*(AMTP*AMBT*(COTB+TANB)*COST*COSB+
     $(AMBT**2*TANB+AMTP**2*COTB-AMW**2*SIN(2*BETA))*
     $SINT*SINB+AMTP*(TWOM1-AAT*COTB)*SINB*COST+AMBT*
     $(TWOM1-AAB*TANB)*COSB*SINT)
        AS=A*A
        WID=AS*SQRT(SSXLAM(AMT2SS**2,AMB2SS**2,AMHC**2))/
     $      16./PI/AMT2SS**3
        CALL SSSAVE(ISTP2,WID,ISHC,ISBT2,0,0,0)
      END IF
C
C
C          stop_2 -> top + zino_i
C
      DO 500 IZ=1,4
        MZIZ=ABS(AMZISS(IZ))
        SNZI=SIGN(1.,AMZISS(IZ))
        IF (SNZI.EQ.1.) THEN
           THIZ=0.
        ELSE
           THIZ=1.
        END IF
        ZAUIZ=ZI**(THIZ-1.)*SNZI
     $  *(-G/SR2*ZMIXSS(3,IZ)-GP/3./SR2*ZMIXSS(4,IZ))
        ZBUIZ=ZI**(THIZ-1.)*4*GP*ZMIXSS(4,IZ)/3./SR2
        ZPP=ZI**THIZ
        ZPM=(-ZI)**THIZ
        ZA=((ZI*ZAUIZ-ZPP*FT*ZMIXSS(1,IZ))*SINT
     $  +(ZI*ZBUIZ-ZPM*FT*ZMIXSS(1,IZ))*COST)/2.
        ZB=((-ZI*ZAUIZ-ZPP*FT*ZMIXSS(1,IZ))*SINT
     $  +(ZI*ZBUIZ+ZPM*FT*ZMIXSS(1,IZ))*COST)/2.
        AS=ZA*CONJG(ZA)
        BS=ZB*CONJG(ZB)
        IF (AMT2SS.GT.(AMTP+MZIZ)) THEN
          WID=(AS*(AMT2SS**2-(AMTP+MZIZ)**2)+BS
     $    *(AMT2SS**2-(AMTP-MZIZ)**2))/8./PI/AMT2SS
     $    *SQRT(SSXLAM(1.,AMTP**2/AMT2SS**2,MZIZ**2/AMT2SS**2))
          CALL SSSAVE(ISTP2,WID,ISZIZ(IZ),IDTP,0,0,0)
        END IF
500   CONTINUE
C     Implement Andrew Box t1-> b+W+Z1 decay
C
      IF (AMT1SS.GT.(AMBT+AMW+MZ1)) THEN
          CALL STBWZ1(DWID)
          CALL SSSAVE(ISTP1,SNGL(DWID),ISZIZ(1),IDW,IDBT,0,0)
      END IF
C 
C          Normalize branching ratios 
C
       CALL SSNORM(ISTP1)
       CALL SSNORM(ISTP2)
C
       RETURN
       END
