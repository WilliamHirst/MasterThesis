CDECK  ID>, SSGLBF.
        SUBROUTINE SSGLBF
C-----------------------------------------------------------------------
C
C       This subroutine gives gluino branching fractions to gauginos
C       according to Baer,Barger,Karatas,Tata (Phys.Rev.D36,96(1987)
C       (Now includes gluino->gluon+zino1,2,3,4 loop decays. Jan 1990)
C       Also includes contribution due to non-degenerate t1-t2 stops
C       Also includes contribution due to non-degenerate b_L and b_R
C       Updated to include mixed sbottom states b1 and b2: 10/9/96
C
C       Auxiliary functions are called SSGxyi, where normally x 
C       indicates the SUSY particle, y the SM particle(s), and i is
C       a counter.
C
C       Baer's GLUBF
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
C          Temporary parameters for functions
      COMMON/SSTMP/TMP(10),ITMP(10)
      REAL TMP
      INTEGER ITMP
      SAVE /SSTMP/
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
C
      EXTERNAL SSGWQ1,SSGWQ2,SSGZG1,SSGZG2,SSGZG3,SSGZT
      EXTERNAL SSGWT1,SSGWT2,SSGWT3,SSGWT4,SSGWT5,SSGWT6,
     $SSGWT7,SSGWT8
      EXTERNAL SSGX1,SSGX2,SSGX3,SSGX4,SSGX5,SSGX6,SSGX7,SSGX8
     $,SSGX9
      REAL WIDU,WIDD,WIDI,WIDL,WIDR
      REAL SSXINT,SSXLAM,XUPPER,XTCON,SUM,GMQK,C,G,FT,GP
     $,UPPER,GS,FB,TANB,PI,SR2,UPPR,ALF3
      REAL MW1,MW2,SNW1,SNW2,XM,YM,THX,THY,COST,SINT,COSB,SINB
      REAL MZIZ,AUIZ,ADIZ,BUIZ,BDIZ
      REAL WID,SNIZ,FACT
      REAL XT1,XT2,XT3,XT4,XT5,XT6,XT7,XT8,XT9,AL1,AL2,BE1,BE2
      REAL XLL,XRR,XL1R1,XL2R2,XL1R2,XL2R1,XL1L2,XR1R2,XLR1
      REAL XMST1,XMST2,XMST12,XLR2,XMSB1,XMSB2,XMSB12
      REAL BW1,BW2,GT1,GT2,GT1T2,GB1,GB2,GT1B1,GT1B2,GT2B1,GT2B2
      REAL KUL,KUR,KDL,KDR,KCL,KCR,KSL,KSR,KBL,KBR,KT1,KT2
      REAL XKUL,XKUR,XKDL,XKDR,XKSL,XKSR,XKCL,XKCR,XKBL,XKBR,
     $XKT1,XKT2
      REAL XI1UL,XI1UR,XI1DL,XI1DR,XI1SL,XI1SR,XI1CL,XI1CR,
     $XI1BL,XI1BR,XI1T1,XI1T2,XIT1,XIT2
      REAL ALT1,ART1,ALT2,ART2,ALB1,ARB1,ALB2,ARB2
      REAL BETA,AMPL
      INTEGER IZ,ISZI(4),THIZ
      COMPLEX ZONE,ZI,ZAT1(2),ZAT2(2),ZADW1,ZADW2,ZAUW1,ZAUW2
      COMPLEX ZAUIZ,ZADIZ,ZBUIZ,ZBDIZ,Z1(2),Z2(2)
      DOUBLE PRECISION SSALFS
      SAVE ZONE,ZI
      DATA ZONE,ZI/(1.,0.),(0.,1.)/

C
C          Partly duplicated from SSMASS.
C
      AMPL=2.4E18
      PI=4.*ATAN(1.)
      SR2=SQRT(2.)
      G=SQRT(4*PI*ALFAEM/SN2THW)
      GP=G*SQRT(SN2THW/(1.-SN2THW))
      ALF3=SSALFS(DBLE(AMGLSS**2))
      GS=SQRT(4.*PI*ALF3)
      TANB=1./RV2V1
      BETA=ATAN(TANB)
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
      COST=COS(THETAT)
      SINT=SIN(THETAT)
      COSB=COS(THETAB)
      SINB=SIN(THETAB)
C
      ZADW1=ZI*G*SNW1*SIN(GAMMAR)
      ZAUW1=ZI*G*SIN(GAMMAL)
      ZADW2=ZI*G*SNW2*COS(GAMMAR)*THY
      ZAUW2=ZI*G*COS(GAMMAL)*THX
      BW1=-FT*SNW1*COS(GAMMAR)
      BW2=FT*SNW2*SIN(GAMMAR)*THY
      ZAT1(1)=ZADW1*COST+ZI*BW1*SINT
      ZAT2(1)=ZADW1*SINT-ZI*BW1*COST
      ZAT1(2)=ZADW2*COST+ZI*BW2*SINT
      ZAT2(2)=ZADW2*SINT-ZI*BW2*COST
      C=SSALFS(DBLE(AMGLSS**2))*AMGLSS/8./PI**2
C
C          gluino --> w1 + qk + qb
C
C     Now includes sbottom as well as stop mixing/Yukawa effects
C     as of 3/31/97, thanks to M. Drees
      TMP(1)=MW1
      UPPR=(AMGLSS**2-MW1**2)/2./AMGLSS
      IF (AMGLSS.GT.(MW1+AMUP+AMDN)) THEN
        IF (AMGLSS.LT.AMULSS.AND.AMGLSS.LT.AMDLSS) THEN
          TMP(2)=AMULSS
          WIDU=ZADW1*CONJG(ZADW1)*SSXINT(0.,SSGWQ1,UPPR)
          TMP(2)=AMDLSS
          WIDD=ZAUW1*CONJG(ZAUW1)*SSXINT(0.,SSGWQ1,UPPR)
          TMP(2)=AMULSS
          TMP(3)=AMDLSS
          WIDI=-SGNM3*2*REAL(ZAUW1*ZADW1)*SSXINT(0.,SSGWQ2,UPPR)
          WID=GS**2/2./AMGLSS/(2*PI)**5*(WIDU+WIDD+WIDI)
        ELSE IF (AMGLSS.LT.AMULSS.AND.AMGLSS.GE.AMDLSS) THEN
          TMP(2)=AMULSS
          WIDU=ZADW1*CONJG(ZADW1)*SSXINT(0.,SSGWQ1,UPPR)
          WID=GS**2/2./AMGLSS/(2*PI)**5*WIDU
        ELSE IF (AMGLSS.LT.AMDLSS.AND.AMGLSS.GE.AMULSS) THEN
          TMP(2)=AMDLSS
          WIDD=ZAUW1*CONJG(ZAUW1)*SSXINT(0.,SSGWQ1,UPPR)
          WID=GS**2/2./AMGLSS/(2*PI)**5*WIDD
        ELSE
          WID=0.
        END IF
        CALL SSSAVE(ISGL,WID,+ISW1,+IDDN,-IDUP,0,0)
        Z1(1)=1.
        Z1(2)=-Z1(1)
        Z2(1)=G*SIN(GAMMAR)
        Z2(2)=Z2(1)
        CALL SSME3(2,AMULSS,Z1,Z2)
        Z1(1)=G*SIN(GAMMAL)*SNW1
        Z1(2)=-Z1(1)
        Z2(1)=1.
        Z2(2)=Z2(1)
        CALL SSME3(3,AMDLSS,Z1,Z2)
        CALL SSSAVE(ISGL,WID,-ISW1,+IDUP,-IDDN,0,0)
        Z1(1)=1.
        Z1(2)=-Z1(1)
        Z2(1)=G*SIN(GAMMAL)
        Z2(2)=Z2(1)
        CALL SSME3(2,AMDLSS,Z1,Z2)
        Z1(1)=G*SIN(GAMMAR)*SNW1
        Z1(2)=-Z1(1)
        Z2(1)=1.
        Z2(2)=Z2(1)
        CALL SSME3(3,AMULSS,Z1,Z2)
      END IF
C
      IF (AMGLSS.GT.(MW1+AMST+AMCH)) THEN
        IF (AMGLSS.LT.AMCLSS.AND.AMGLSS.LT.AMSLSS) THEN
          TMP(2)=AMCLSS
          WIDU=ZADW1*CONJG(ZADW1)*SSXINT(0.,SSGWQ1,UPPR)
          TMP(2)=AMSLSS
          WIDD=ZAUW1*CONJG(ZAUW1)*SSXINT(0.,SSGWQ1,UPPR)
          TMP(2)=AMCLSS
          TMP(3)=AMSLSS
          WIDI=-SGNM3*2*REAL(ZAUW1*ZADW1)*SSXINT(0.,SSGWQ2,UPPR)
          WID=GS**2/2./AMGLSS/(2*PI)**5*(WIDU+WIDD+WIDI)
        ELSE IF (AMGLSS.LT.AMCLSS.AND.AMGLSS.GE.AMSLSS) THEN
          TMP(2)=AMCLSS
          WIDU=ZADW1*CONJG(ZADW1)*SSXINT(0.,SSGWQ1,UPPR)
          WID=GS**2/2./AMGLSS/(2*PI)**5*WIDU
        ELSE IF (AMGLSS.LT.AMSLSS.AND.AMGLSS.GE.AMCLSS) THEN
          TMP(2)=AMSLSS
          WIDD=ZAUW1*CONJG(ZAUW1)*SSXINT(0.,SSGWQ1,UPPR)
          WID=GS**2/2./AMGLSS/(2*PI)**5*WIDD
        ELSE
          WID=0.
        END IF
        CALL SSSAVE(ISGL,WID,+ISW1,+IDST,-IDCH,0,0)
        Z1(1)=1.
        Z1(2)=-Z1(1)
        Z2(1)=G*SIN(GAMMAR)
        Z2(2)=Z2(1)
        CALL SSME3(2,AMCLSS,Z1,Z2)
        Z1(1)=G*SIN(GAMMAL)*SNW1
        Z1(2)=-Z1(1)
        Z2(1)=1.
        Z2(2)=Z2(1)
        CALL SSME3(3,AMSLSS,Z1,Z2)
        CALL SSSAVE(ISGL,WID,-ISW1,+IDCH,-IDST,0,0)
        Z1(1)=1.
        Z1(2)=-Z1(1)
        Z2(1)=G*SIN(GAMMAL)
        Z2(2)=Z2(1)
        CALL SSME3(2,AMSLSS,Z1,Z2)
        Z1(1)=G*SIN(GAMMAR)*SNW1
        Z1(2)=-Z1(1)
        Z2(1)=1.
        Z2(2)=Z2(1)
        CALL SSME3(3,AMCLSS,Z1,Z2)
      END IF
C
      IF (AMGLSS.GT.(MW1+AMTP+AMBT)) THEN
        ALT1=-G*SIN(GAMMAR)*COST+FT*COS(GAMMAR)*SINT
        ART1=-FB*COS(GAMMAL)*COST
        ALT2=-G*SIN(GAMMAR)*SINT-FT*COS(GAMMAR)*COST
        ART2=-FB*COS(GAMMAL)*SINT
        ALB1=-G*SIN(GAMMAL)*COSB+FB*COS(GAMMAL)*SINB
        ARB1=-FT*COS(GAMMAR)*COSB
        ALB2=-G*SIN(GAMMAL)*SINB-FB*COS(GAMMAL)*COSB
        ARB2=-FT*COS(GAMMAR)*SINB
        UPPER=(AMGLSS**2+AMTP**2-(MW1+AMBT)**2)/2./AMGLSS
        FACT=GS**2*PI**2/(2*PI)**5/2./AMGLSS
        TMP(1)=MW1
        TMP(2)=AMGLSS
        TMP(3)=AMTP
        IF (AMGLSS.LT.(AMTP+AMT1SS)) THEN
          TMP(6)=AMT1SS
          TMP(7)=AMT1SS
          GT1=(ALT1**2+ART1**2)*(SSXINT(AMTP,SSGWT1,UPPER)+
     ,        SIN(2*THETAT)*SSXINT(AMTP,SSGWT8,UPPER))
          TMP(7)=AMT2SS
          GT1T2=2*(ALT1*ALT2+ART1*ART2)*(SINT**2-COST**2)*
     ,          SSXINT(AMTP,SSGWT8,UPPER)
        ELSE
          GT1=0.
          GT1T2=0.
        END IF
        IF (AMGLSS.LT.(AMTP+AMT2SS)) THEN
          TMP(6)=AMT2SS
          TMP(7)=AMT2SS
          GT2=(ALT2**2+ART2**2)*(SSXINT(AMTP,SSGWT1,UPPER)-
     ,        SIN(2*THETAT)*SSXINT(AMTP,SSGWT8,UPPER))
        ELSE
          GT2=0.
        END IF
        IF (AMGLSS.LT.(AMB1SS+AMBT)) THEN
          TMP(4)=AMB1SS
          TMP(8)=SNW1
C          Rewrite UPPR=(AMGLSS**2+AMBT**2-(AMTP+MW1)**2)/2./AMGLSS
          UPPR=((AMGLSS-AMTP-MW1)*(AMGLSS+AMTP+MW1)+AMBT**2)/(2*AMGLSS)
          GB1=(ALB1**2+ARB1**2)*SSXINT(AMBT,SSGWT2,UPPR)-
     ,        ALB1*ARB1*SSXINT(AMBT,SSGWT3,UPPR)
        ELSE
          GB1=0.
        END IF
        IF (AMGLSS.LT.(AMB2SS+AMBT)) THEN
          TMP(4)=AMB2SS
          TMP(8)=SNW1
C          Rewrite UPPR=(AMGLSS**2+AMBT**2-(AMTP+MW1)**2)/2./AMGLSS
          UPPR=((AMGLSS-AMTP-MW1)*(AMGLSS+AMTP+MW1)+AMBT**2)/(2*AMGLSS)
          GB2=(ALB2**2+ARB2**2)*SSXINT(AMBT,SSGWT2,UPPR)-
     ,        ALB2*ARB2*SSXINT(AMBT,SSGWT3,UPPR)
        ELSE
          GB2=0.
        END IF
        IF (AMGLSS.LT.(AMTP+AMT1SS).AND.AMGLSS.LT.
     $    (AMB1SS+AMBT)) THEN
          TMP(4)=AMB1SS
          TMP(6)=AMT1SS
          TMP(8)=SNW1
          GT1B1=(COST*SINB*ART1*ALB1+SINT*COSB*ALT1*ARB1)*
     ,           SSXINT(AMTP,SSGWT6,UPPER)-
     ,          (COST*COSB*ALT1*ALB1+SINT*SINB*ART1*ARB1)*
     ,           SSXINT(AMTP,SSGWT4,UPPER)+
     ,          (COST*COSB*ALT1*ARB1+SINT*SINB*ART1*ALB1)*
     ,           SSXINT(AMTP,SSGWT5,UPPER)-
     ,          (COST*SINB*ART1*ARB1+SINT*COSB*ALT1*ALB1)*
     ,           SSXINT(AMTP,SSGWT7,UPPER)
        ELSE
          GT1B1=0.
        END IF
        IF (AMGLSS.LT.(AMTP+AMT1SS).AND.AMGLSS.LT.
     $    (AMB2SS+AMBT)) THEN
          TMP(4)=AMB2SS
          TMP(6)=AMT1SS
          TMP(8)=SNW1
          GT1B2=(-COST*COSB*ART1*ALB2+SINT*SINB*ALT1*ARB2)*
     ,           SSXINT(AMTP,SSGWT6,UPPER)-
     ,          (COST*SINB*ALT1*ALB2-SINT*COSB*ART1*ARB2)*
     ,           SSXINT(AMTP,SSGWT4,UPPER)+
     ,          (COST*SINB*ALT1*ARB2-SINT*COSB*ART1*ALB2)*
     ,           SSXINT(AMTP,SSGWT5,UPPER)-
     ,          (-COST*COSB*ART1*ARB2+SINT*SINB*ALT1*ALB2)*
     ,           SSXINT(AMTP,SSGWT7,UPPER)
        ELSE
          GT1B2=0.
        END IF
        IF (AMGLSS.LT.(AMTP+AMT2SS).AND.AMGLSS.LT.
     $    (AMB1SS+AMBT)) THEN
          TMP(4)=AMB1SS
          TMP(6)=AMT2SS
          TMP(8)=SNW1
          GT2B1=(SINT*SINB*ART2*ALB1-COST*COSB*ALT2*ARB1)*
     ,           SSXINT(AMTP,SSGWT6,UPPER)-
     ,          (SINT*COSB*ALT2*ALB1-COST*SINB*ART2*ARB1)*
     ,           SSXINT(AMTP,SSGWT4,UPPER)+
     ,          (SINT*COSB*ALT2*ARB1-COST*SINB*ART2*ALB1)*
     ,           SSXINT(AMTP,SSGWT5,UPPER)-
     ,          (SINT*SINB*ART2*ARB1-COST*COSB*ALT2*ALB1)*
     ,           SSXINT(AMTP,SSGWT7,UPPER)
        ELSE
          GT2B1=0.
        END IF
        IF (AMGLSS.LT.(AMTP+AMT2SS).AND.AMGLSS.LT.
     $    (AMB2SS+AMBT)) THEN
          TMP(4)=AMB2SS
          TMP(6)=AMT2SS
          TMP(8)=SNW1
          GT2B2=(-SINT*COSB*ART2*ALB2-COST*SINB*ALT2*ARB2)*
     ,           SSXINT(AMTP,SSGWT6,UPPER)-
     ,          (SINT*SINB*ALT2*ALB2+COST*COSB*ART2*ARB2)*
     ,           SSXINT(AMTP,SSGWT4,UPPER)+
     ,          (SINT*SINB*ALT2*ARB2+COST*COSB*ART2*ALB2)*
     ,           SSXINT(AMTP,SSGWT5,UPPER)-
     ,          (-SINT*COSB*ART2*ARB2-COST*SINB*ALT2*ALB2)*
     ,           SSXINT(AMTP,SSGWT7,UPPER)
        ELSE
          GT2B2=0.
        END IF
        WID=GT1+GT2+GT1T2+GB1+GB2+GT1B1+GT1B2+GT2B1+GT2B2
        WID=FACT*WID
        IF (WID.GT.0.) THEN
        CALL SSSAVE(ISGL,WID,+ISW1,+IDBT,-IDTP,0,0)
        Z1(1)=SINT-COST
        Z1(2)=SINT+COST
        Z2(1)=ALT1+SNW1*ART1
        Z2(2)=ALT1-SNW1*ART1
        CALL SSME3(2,AMT1SS,Z1,Z2)
        Z1(1)=-COST-SINT
        Z1(2)=-COST+SINT
        Z2(1)=ALT2+SNW1*ART2
        Z2(2)=ALT2-SNW1*ART2
        CALL SSME3(2,AMT2SS,Z1,Z2)
        Z1(1)=SNW1*ALB1+ARB1
        Z1(2)=-SNW1*ALB1+ARB1
        Z2(1)=-COSB+SINB
        Z2(2)=-COSB-SINB
        CALL SSME3(3,AMB1SS,Z1,Z2)
        Z1(1)=SNW1*ALB2+ARB2
        Z1(2)=-SNW1*ALB2+ARB2
        Z2(1)=-SINB-COSB
        Z2(2)=-SINB+COSB
        CALL SSME3(3,AMB2SS,Z1,Z2)
        CALL SSSAVE(ISGL,WID,-ISW1,+IDTP,-IDBT,0,0)
        Z1(1)=SINB-COSB
        Z1(2)=SINB+COSB
        Z2(1)=ALB1+SNW1*ARB1
        Z2(2)=ALB1-SNW1*ARB1
        CALL SSME3(2,AMB1SS,Z1,Z2)
        Z1(1)=-COSB-SINB
        Z1(2)=-COSB+SINB
        Z2(1)=ALB2+SNW1*ARB2
        Z2(2)=ALB2-SNW1*ARB2
        CALL SSME3(2,AMB2SS,Z1,Z2)
        Z1(1)=SNW1*ALT1+ART1
        Z1(2)=-SNW1*ALT1+ART1
        Z2(1)=-COST+SINT
        Z2(2)=-COST-SINT
        CALL SSME3(3,AMT1SS,Z1,Z2)
        Z1(1)=SNW1*ALT2+ART2
        Z1(2)=-SNW1*ALT2+ART2
        Z2(1)=-SINT-COST
        Z2(2)=-SINT+COST
        CALL SSME3(3,AMT2SS,Z1,Z2)
        END IF
      END IF
C
C          gluino --> w2 + qk + qb
C
      TMP(1)=MW2
      UPPR=(AMGLSS**2-MW2**2)/2./AMGLSS
      IF (AMGLSS.GT.(MW2+AMUP+AMDN)) THEN
        IF (AMGLSS.LT.AMULSS.AND.AMGLSS.LT.AMDLSS) THEN
          TMP(2)=AMULSS
          WIDU=ZADW2*CONJG(ZADW2)*SSXINT(0.,SSGWQ1,UPPR)
          TMP(2)=AMDLSS
          WIDD=ZAUW2*CONJG(ZAUW2)*SSXINT(0.,SSGWQ1,UPPR)
          TMP(2)=AMULSS
          TMP(3)=AMDLSS
          WIDI=-SGNM3*2*REAL(ZAUW2*ZADW2)*SSXINT(0.,SSGWQ2,UPPR)
          WID=GS**2/2./AMGLSS/(2*PI)**5*(WIDU+WIDD+WIDI)
        ELSE IF (AMGLSS.LT.AMULSS.AND.AMGLSS.GE.AMDLSS) THEN
          TMP(2)=AMULSS
          WIDU=ZADW2*CONJG(ZADW2)*SSXINT(0.,SSGWQ1,UPPR)
          WID=GS**2/2./AMGLSS/(2*PI)**5*WIDU
        ELSE IF (AMGLSS.LT.AMDLSS.AND.AMGLSS.GE.AMULSS) THEN
          TMP(2)=AMDLSS
          WIDD=ZAUW2*CONJG(ZAUW2)*SSXINT(0.,SSGWQ1,UPPR)
          WID=GS**2/2./AMGLSS/(2*PI)**5*WIDD
        ELSE
          WID=0.
        END IF
        CALL SSSAVE(ISGL,WID,+ISW2,+IDDN,-IDUP,0,0)
        Z1(1)=1.
        Z1(2)=-Z1(1)
        Z2(1)=G*THY*COS(GAMMAR)
        Z2(2)=Z2(1)
        CALL SSME3(2,AMULSS,Z1,Z2)
        Z1(1)=G*THX*COS(GAMMAL)*SNW2
        Z1(2)=-Z1(1)
        Z2(1)=1.
        Z2(2)=Z2(1)
        CALL SSME3(3,AMDLSS,Z1,Z2)
        CALL SSSAVE(ISGL,WID,-ISW2,+IDUP,-IDDN,0,0)
        Z1(1)=1.
        Z1(2)=-Z1(1)
        Z2(1)=G*THX*COS(GAMMAL)
        Z2(2)=Z2(1)
        CALL SSME3(2,AMDLSS,Z1,Z2)
        Z1(1)=G*THY*COS(GAMMAR)*SNW2
        Z1(2)=-Z1(1)
        Z2(1)=1.
        Z2(2)=Z2(1)
        CALL SSME3(3,AMULSS,Z1,Z2)
      END IF
C
      IF (AMGLSS.GT.(MW2+AMST+AMCH)) THEN
        IF (AMGLSS.LT.AMCLSS.AND.AMGLSS.LT.AMSLSS) THEN
          TMP(2)=AMCLSS
          WIDU=ZADW2*CONJG(ZADW2)*SSXINT(0.,SSGWQ1,UPPR)
          TMP(2)=AMSLSS
          WIDD=ZAUW2*CONJG(ZAUW2)*SSXINT(0.,SSGWQ1,UPPR)
          TMP(2)=AMCLSS
          TMP(3)=AMSLSS
          WIDI=-SGNM3*2*REAL(ZAUW2*ZADW2)*SSXINT(0.,SSGWQ2,UPPR)
          WID=GS**2/2./AMGLSS/(2*PI)**5*(WIDU+WIDD+WIDI)
        ELSE IF (AMGLSS.LT.AMCLSS.AND.AMGLSS.GE.AMSLSS) THEN
          TMP(2)=AMCLSS
          WIDU=ZADW2*CONJG(ZADW2)*SSXINT(0.,SSGWQ1,UPPR)
          WID=GS**2/2./AMGLSS/(2*PI)**5*WIDU
        ELSE IF (AMGLSS.LT.AMSLSS.AND.AMGLSS.GE.AMCLSS) THEN
          TMP(2)=AMSLSS
          WIDD=ZAUW2*CONJG(ZAUW2)*SSXINT(0.,SSGWQ1,UPPR)
          WID=GS**2/2./AMGLSS/(2*PI)**5*WIDD
        ELSE
          WID=0.
        END IF
        CALL SSSAVE(ISGL,WID,+ISW2,+IDST,-IDCH,0,0)
        Z1(1)=1.
        Z1(2)=-Z1(1)
        Z2(1)=G*THY*COS(GAMMAR)
        Z2(2)=Z2(1)
        CALL SSME3(2,AMCLSS,Z1,Z2)
        Z1(1)=G*THX*COS(GAMMAL)*SNW2
        Z1(2)=-Z1(1)
        Z2(1)=1.
        Z2(2)=Z2(1)
        CALL SSME3(3,AMSLSS,Z1,Z2)
        CALL SSSAVE(ISGL,WID,-ISW2,+IDCH,-IDST,0,0)
        Z1(1)=1.
        Z1(2)=-Z1(1)
        Z2(1)=G*THX*COS(GAMMAL)
        Z2(2)=Z2(1)
        CALL SSME3(2,AMSLSS,Z1,Z2)
        Z1(1)=G*THY*COS(GAMMAR)*SNW2
        Z1(2)=-Z1(1)
        Z2(1)=1.
        Z2(2)=Z2(1)
        CALL SSME3(3,AMCLSS,Z1,Z2)
      END IF
C
      IF (AMGLSS.GT.(MW2+AMTP+AMBT)) THEN
        ALT1=-G*THY*COS(GAMMAR)*COST-FT*THY*SIN(GAMMAR)*SINT
        ART1=FB*THX*SIN(GAMMAL)*COST
        ALT2=-G*THY*COS(GAMMAR)*SINT+FT*THY*SIN(GAMMAR)*COST
        ART2=FB*THX*SIN(GAMMAL)*SINT
        ALB1=-G*THX*COS(GAMMAL)*COSB-FB*THX*SIN(GAMMAL)*SINB
        ARB1=FT*THY*SIN(GAMMAR)*COSB
        ALB2=-G*THX*COS(GAMMAL)*SINB+FB*THX*SIN(GAMMAL)*COSB
        ARB2=FT*THY*SIN(GAMMAR)*SINB
        UPPER=(AMGLSS**2+AMTP**2-(MW2+AMBT)**2)/2./AMGLSS
        FACT=GS**2*PI**2/(2*PI)**5/2./AMGLSS
        TMP(1)=MW2
        TMP(2)=AMGLSS
        TMP(3)=AMTP
        IF (AMGLSS.LT.(AMTP+AMT1SS)) THEN
          TMP(6)=AMT1SS
          TMP(7)=AMT1SS
          GT1=(ALT1**2+ART1**2)*(SSXINT(AMTP,SSGWT1,UPPER)+
     ,        SIN(2*THETAT)*SSXINT(AMTP,SSGWT8,UPPER))
          TMP(7)=AMT2SS
          GT1T2=2*(ALT1*ALT2+ART1*ART2)*(SINT**2-COST**2)*
     ,          SSXINT(AMTP,SSGWT8,UPPER)
        ELSE
          GT1=0.
          GT1T2=0.
        END IF
        IF (AMGLSS.LT.(AMTP+AMT2SS)) THEN
          TMP(6)=AMT2SS
          TMP(7)=AMT2SS
          GT2=(ALT2**2+ART2**2)*(SSXINT(AMTP,SSGWT1,UPPER)-
     ,        SIN(2*THETAT)*SSXINT(AMTP,SSGWT8,UPPER))
        ELSE
          GT2=0.
        END IF
        IF (AMGLSS.LT.(AMB1SS+AMBT)) THEN
          TMP(4)=AMB1SS
          TMP(8)=SNW2
C          Rewrite UPPR=(AMGLSS**2+AMBT**2-(AMTP+MW2)**2)/2./AMGLSS
          UPPR=((AMGLSS-AMTP-MW2)*(AMGLSS+AMTP+MW2)+AMBT**2)/(2*AMGLSS)
          GB1=(ALB1**2+ARB1**2)*SSXINT(AMBT,SSGWT2,UPPR)-
     ,        ALB1*ARB1*SSXINT(AMBT,SSGWT3,UPPR)
        ELSE
          GB1=0.
        END IF
        IF (AMGLSS.LT.(AMB2SS+AMBT)) THEN
          TMP(4)=AMB2SS
          TMP(8)=SNW2
          UPPR=((AMGLSS-AMTP-MW2)*(AMGLSS+AMTP+MW2)+AMBT**2)/(2*AMGLSS)
          GB2=(ALB2**2+ARB2**2)*SSXINT(AMBT,SSGWT2,UPPR)-
     ,        ALB2*ARB2*SSXINT(AMBT,SSGWT3,UPPR)
        ELSE
          GB2=0.
        END IF
        IF (AMGLSS.LT.(AMTP+AMT1SS).AND.AMGLSS.LT.
     $    (AMB1SS+AMBT)) THEN
          TMP(4)=AMB1SS
          TMP(6)=AMT1SS
          TMP(8)=SNW2
          GT1B1=(COST*SINB*ART1*ALB1+SINT*COSB*ALT1*ARB1)*
     ,           SSXINT(AMTP,SSGWT6,UPPER)-
     ,          (COST*COSB*ALT1*ALB1+SINT*SINB*ART1*ARB1)*
     ,           SSXINT(AMTP,SSGWT4,UPPER)+
     ,          (COST*COSB*ALT1*ARB1+SINT*SINB*ART1*ALB1)*
     ,           SSXINT(AMTP,SSGWT5,UPPER)-
     ,          (COST*SINB*ART1*ARB1+SINT*COSB*ALT1*ALB1)*
     ,           SSXINT(AMTP,SSGWT7,UPPER)
        ELSE
          GT1B1=0.
        END IF
        IF (AMGLSS.LT.(AMTP+AMT1SS).AND.AMGLSS.LT.
     $    (AMB2SS+AMBT)) THEN
          TMP(4)=AMB2SS
          TMP(6)=AMT1SS
          TMP(8)=SNW2
          GT1B2=(-COST*COSB*ART1*ALB2+SINT*SINB*ALT1*ARB2)*
     ,           SSXINT(AMTP,SSGWT6,UPPER)-
     ,          (COST*SINB*ALT1*ALB2-SINT*COSB*ART1*ARB2)*
     ,           SSXINT(AMTP,SSGWT4,UPPER)+
     ,          (COST*SINB*ALT1*ARB2-SINT*COSB*ART1*ALB2)*
     ,           SSXINT(AMTP,SSGWT5,UPPER)-
     ,          (-COST*COSB*ART1*ARB2+SINT*SINB*ALT1*ALB2)*
     ,           SSXINT(AMTP,SSGWT7,UPPER)
        ELSE
          GT1B2=0.
        END IF
        IF (AMGLSS.LT.(AMTP+AMT2SS).AND.AMGLSS.LT.
     $    (AMB1SS+AMBT)) THEN
          TMP(4)=AMB1SS
          TMP(6)=AMT2SS
          TMP(8)=SNW2
          GT2B1=(SINT*SINB*ART2*ALB1-COST*COSB*ALT2*ARB1)*
     ,           SSXINT(AMTP,SSGWT6,UPPER)-
     ,          (SINT*COSB*ALT2*ALB1-COST*SINB*ART2*ARB1)*
     ,           SSXINT(AMTP,SSGWT4,UPPER)+
     ,          (SINT*COSB*ALT2*ARB1-COST*SINB*ART2*ALB1)*
     ,           SSXINT(AMTP,SSGWT5,UPPER)-
     ,          (SINT*SINB*ART2*ARB1-COST*COSB*ALT2*ALB1)*
     ,           SSXINT(AMTP,SSGWT7,UPPER)
        ELSE
          GT2B1=0.
        END IF
        IF (AMGLSS.LT.(AMTP+AMT2SS).AND.AMGLSS.LT.
     $    (AMB2SS+AMBT)) THEN
          TMP(4)=AMB2SS
          TMP(6)=AMT2SS
          TMP(8)=SNW2
          GT2B2=(-SINT*COSB*ART2*ALB2-COST*SINB*ALT2*ARB2)*
     ,           SSXINT(AMTP,SSGWT6,UPPER)-
     ,          (SINT*SINB*ALT2*ALB2+COST*COSB*ART2*ARB2)*
     ,           SSXINT(AMTP,SSGWT4,UPPER)+
     ,          (SINT*SINB*ALT2*ARB2+COST*COSB*ART2*ALB2)*
     ,           SSXINT(AMTP,SSGWT5,UPPER)-
     ,          (-SINT*COSB*ART2*ARB2-COST*SINB*ALT2*ALB2)*
     ,           SSXINT(AMTP,SSGWT7,UPPER)
        ELSE
          GT2B2=0.
        END IF
        WID=GT1+GT2+GT1T2+GB1+GB2+GT1B1+GT1B2+GT2B1+GT2B2
        WID=FACT*WID
        IF (WID.GT.0.) THEN
        CALL SSSAVE(ISGL,WID,+ISW2,+IDBT,-IDTP,0,0)
        Z1(1)=SINT-COST
        Z1(2)=SINT+COST
        Z2(1)=ALT1+SNW2*ART1
        Z2(2)=ALT1-SNW2*ART1
        CALL SSME3(2,AMT1SS,Z1,Z2)
        Z1(1)=-COST-SINT
        Z1(2)=-COST+SINT
        Z2(1)=ALT2+SNW2*ART2
        Z2(2)=ALT2-SNW2*ART2
        CALL SSME3(2,AMT2SS,Z1,Z2)
        Z1(1)=SNW2*ALB1+ARB1
        Z1(2)=-SNW2*ALB1+ARB1
        Z2(1)=-COSB+SINB
        Z2(2)=-COSB-SINB
        CALL SSME3(3,AMB1SS,Z1,Z2)
        Z1(1)=SNW2*ALB2+ARB2
        Z1(2)=-SNW2*ALB2+ARB2
        Z2(1)=-SINB-COSB
        Z2(2)=-SINB+COSB
        CALL SSME3(3,AMB2SS,Z1,Z2)
        CALL SSSAVE(ISGL,WID,-ISW2,+IDTP,-IDBT,0,0)
        Z1(1)=SINB-COSB
        Z1(2)=SINB+COSB
        Z2(1)=ALB1+SNW2*ARB1
        Z2(2)=ALB1-SNW2*ARB1
        CALL SSME3(2,AMB1SS,Z1,Z2)
        Z1(1)=-COSB-SINB
        Z1(2)=-COSB+SINB
        Z2(1)=ALB2+SNW2*ARB2
        Z2(2)=ALB2-SNW2*ARB2
        CALL SSME3(2,AMB2SS,Z1,Z2)
        Z1(1)=SNW2*ALT1+ART1
        Z1(2)=-SNW2*ALT1+ART1
        Z2(1)=-COST+SINT
        Z2(2)=-COST-SINT
        CALL SSME3(3,AMT1SS,Z1,Z2)
        Z1(1)=SNW2*ALT2+ART2
        Z1(2)=-SNW2*ALT2+ART2
        Z2(1)=-SINT-COST
        Z2(2)=-SINT+COST
        CALL SSME3(3,AMT2SS,Z1,Z2)
        END IF
      END IF
C
C       gluino --> zi decays, zi = z1, z2, z3, z4
C       the auiz etc, below are Atilde's etc. of PRD 42,1568 (1990)
C
      ISZI(1)=ISZ1
      ISZI(2)=ISZ2
      ISZI(3)=ISZ3
      ISZI(4)=ISZ4
      DO 100 IZ=1,4
        MZIZ=ABS(AMZISS(IZ))
        AUIZ=G/SR2*ZMIXSS(3,IZ)+GP/3./SR2*ZMIXSS(4,IZ)
        ADIZ=-G/SR2*ZMIXSS(3,IZ)+GP/3./SR2*ZMIXSS(4,IZ)
        BUIZ=4*GP*ZMIXSS(4,IZ)/3./SR2
        BDIZ=-2*GP/3./SR2*ZMIXSS(4,IZ)
        SNIZ=SIGN(1.,AMZISS(IZ))
        THIZ=0
        IF (AMZISS(IZ).LT.0.) THIZ=1
        ZAUIZ=ZI**(THIZ-1)*SNIZ
     $  *(-G/SR2*ZMIXSS(3,IZ)-GP/3./SR2*ZMIXSS(4,IZ))
        ZBUIZ=ZI**(THIZ-1)*4*GP*ZMIXSS(4,IZ)/3./SR2
        ZADIZ=ZI**(THIZ-1)*SNIZ
     $  *(G/SR2*ZMIXSS(3,IZ)-GP/3./SR2*ZMIXSS(4,IZ))
        ZBDIZ=-2*ZI**(THIZ-1)*GP*ZMIXSS(4,IZ)/3./SR2
C          Radiative gluino --> gluon + zi loop decay
        IF (AMGLSS.GT.MZIZ) THEN
          IF (AMGLSS.LT.(AMT1SS+AMTP)) THEN
            TMP(1)=AMTP
            TMP(2)=MZIZ
            TMP(3)=AMT1SS
            XIT1=SSXINT(0.,SSGZG1,1.)
            XI1T1=SSXINT(0.,SSGZG2,1.)
            XKT1=SSXINT(0.,SSGZG3,1.)
          ELSE
            XIT1=0.
            XI1T1=0.
            XKT1=0.
          END IF
          IF (AMGLSS.LT.(AMT2SS+AMTP)) THEN
            TMP(1)=AMTP
            TMP(2)=MZIZ
            TMP(3)=AMT2SS
            XIT2=SSXINT(0.,SSGZG1,1.)
            XI1T2=SSXINT(0.,SSGZG2,1.)
            XKT2=SSXINT(0.,SSGZG3,1.)
          ELSE
            XIT2=0.
            XI1T2=0.
            XKT2=0.
          END IF
C         !!! NEEDS UPDATE FOR SBOTTOM MIXING !!!
          IF (AMGLSS.LT.(AMB1SS+AMBT)) THEN
            TMP(1)=AMBT
            TMP(2)=MZIZ
            TMP(3)=AMB1SS
            XI1BL=SSXINT(0.,SSGZG2,1.)
            XKBL=SSXINT(0.,SSGZG3,1.)
          ELSE
            XI1BL=0.
            XKBL=0.
          END IF
          IF (AMGLSS.LT.(AMB2SS+AMBT)) THEN
            TMP(1)=AMBT
            TMP(2)=MZIZ
            TMP(3)=AMB2SS
            XI1BR=SSXINT(0.,SSGZG2,1.)
            XKBR=SSXINT(0.,SSGZG3,1.)
          ELSE
            XI1BR=0.
            XKBR=0.
          END IF
          IF (AMGLSS.LT.(AMULSS+AMUP)) THEN
            TMP(1)=AMUP
            TMP(2)=MZIZ
            TMP(3)=AMULSS
            XI1UL=SSXINT(0.,SSGZG2,1.)
            XKUL=SSXINT(0.,SSGZG3,1.)
          ELSE
            XI1UL=0.
            XKUL=0.
          END IF
          IF (AMGLSS.LT.(AMURSS+AMUP)) THEN
            TMP(1)=AMUP
            TMP(2)=MZIZ
            TMP(3)=AMURSS
            XI1UR=SSXINT(0.,SSGZG2,1.)
            XKUR=SSXINT(0.,SSGZG3,1.)
          ELSE
            XI1UR=0.
            XKUR=0.
          END IF
          IF (AMGLSS.LT.(AMDLSS+AMDN)) THEN
            TMP(1)=AMDN
            TMP(2)=MZIZ
            TMP(3)=AMDLSS
            XI1DL=SSXINT(0.,SSGZG2,1.)
            XKDL=SSXINT(0.,SSGZG3,1.)
          ELSE
            XI1DL=0.
            XKDL=0.
          END IF
          IF (AMGLSS.LT.(AMDRSS+AMDN)) THEN
            TMP(1)=AMDN
            TMP(2)=MZIZ
            TMP(3)=AMDRSS
            XI1DR=SSXINT(0.,SSGZG2,1.)
            XKDR=SSXINT(0.,SSGZG3,1.)
          ELSE
            XI1DR=0.
            XKDR=0.
          END IF
          IF (AMGLSS.LT.(AMSLSS+AMST)) THEN
            TMP(1)=AMST
            TMP(2)=MZIZ
            TMP(3)=AMSLSS
            XI1SL=SSXINT(0.,SSGZG2,1.)
            XKSL=SSXINT(0.,SSGZG3,1.)
          ELSE
            XI1SL=0.
            XKSL=0.
          END IF
          IF (AMGLSS.LT.(AMSRSS+AMST)) THEN
            TMP(1)=AMST
            TMP(2)=MZIZ
            TMP(3)=AMSRSS
            XI1SR=SSXINT(0.,SSGZG2,1.)
            XKSR=SSXINT(0.,SSGZG3,1.)
          ELSE
            XI1SR=0.
            XKSR=0.
          END IF
          IF (AMGLSS.LT.(AMCLSS+AMCH)) THEN
            TMP(1)=AMCH
            TMP(2)=MZIZ
            TMP(3)=AMCLSS
            XI1CL=SSXINT(0.,SSGZG2,1.)
            XKCL=SSXINT(0.,SSGZG3,1.)
          ELSE
            XI1CL=0.
            XKCL=0.
          END IF
          IF (AMGLSS.LT.(AMCRSS+AMCH)) THEN
            TMP(1)=AMCH
            TMP(2)=MZIZ
            TMP(3)=AMCRSS
            XI1CR=SSXINT(0.,SSGZG2,1.)
            XKCR=SSXINT(0.,SSGZG3,1.)
          ELSE
            XI1CR=0.
            XKCR=0.
          END IF
          KUL=AUIZ*(XKUL*(MZIZ+SGNM3*SNIZ*AMGLSS)+MZIZ*XI1UL)
          KUR=-BUIZ*(XKUR*(MZIZ+SGNM3*SNIZ*AMGLSS)+MZIZ*XI1UR)
          KDL=ADIZ*(XKDL*(MZIZ+SGNM3*SNIZ*AMGLSS)+MZIZ*XI1DL)
          KDR=-BDIZ*(XKDR*(MZIZ+SGNM3*SNIZ*AMGLSS)+MZIZ*XI1DR)
          KCL=AUIZ*(XKCL*(MZIZ+SGNM3*SNIZ*AMGLSS)+MZIZ*XI1CL)
          KCR=-BUIZ*(XKCR*(MZIZ+SGNM3*SNIZ*AMGLSS)+MZIZ*XI1CR)
          KSL=ADIZ*(XKSL*(MZIZ+SGNM3*SNIZ*AMGLSS)+MZIZ*XI1SL)
          KSR=-BDIZ*(XKSR*(MZIZ+SGNM3*SNIZ*AMGLSS)+MZIZ*XI1SR)
          KBL=ADIZ*(XKBL*(MZIZ+SGNM3*SNIZ*AMGLSS)+MZIZ*XI1BL)
          KBR=-BDIZ*(XKBR*(MZIZ+SGNM3*SNIZ*AMGLSS)+MZIZ*XI1BR)
          KT1=(MZIZ*(XKT1+XI1T1)*(AUIZ*COST-FT*ZMIXSS(1,IZ)*SINT)
     $      +SGNM3*SNIZ*AMGLSS*XKT1*(AUIZ*COST-FT*ZMIXSS(1,IZ)*SINT)+
     $      SNIZ*AMTP*XIT1*(BUIZ*SINT+FT*ZMIXSS(1,IZ)*COST))*COST+
     $      (MZIZ*(XKT1+XI1T1)*(-BUIZ*SINT-FT*ZMIXSS(1,IZ)*COST)
     $      +SGNM3*SNIZ*AMGLSS*XKT1*(-BUIZ*SINT-FT*ZMIXSS(1,IZ)*COST)-
     $      SNIZ*AMTP*XIT1*(AUIZ*COST-FT*ZMIXSS(1,IZ)*SINT))*SINT
          KT2=(MZIZ*(XKT2+XI1T2)*(AUIZ*SINT+FT*ZMIXSS(1,IZ)*COST)
     $      +SGNM3*SNIZ*AMGLSS*XKT2*(AUIZ*SINT+FT*ZMIXSS(1,IZ)*COST)+
     $      SNIZ*AMTP*XIT2*(-BUIZ*COST+FT*ZMIXSS(1,IZ)*SINT))*SINT+
     $      (-MZIZ*(XKT2+XI1T2)*(BUIZ*COST-FT*ZMIXSS(1,IZ)*SINT)
     $      -SGNM3*SNIZ*AMGLSS*XKT2*(BUIZ*COST-FT*ZMIXSS(1,IZ)*SINT)+
     $      SNIZ*AMTP*XIT2*(AUIZ*SINT+FT*ZMIXSS(1,IZ)*COST))*COST
          SUM=(KUL+KUR+KDL+KDR+KSL+KSR+KCL+KCR+KBL+KBR+KT1+
     $         KT2)**2/AMGLSS**2
          WID=ALF3**2*AMGLSS*(1.-MZIZ**2/AMGLSS**2)/256./PI**3*SUM
          CALL SSSAVE(ISGL,WID,ISZI(IZ),IDGL,0,0,0)
        END IF
C          3 body gluino --> q + qb + zi decay, q=u,d
        UPPR=(AMGLSS**2-MZIZ**2)/2./AMGLSS
        IF (AMGLSS.GT.(MZIZ+2*AMUP)) THEN
          IF (AMGLSS.LT.AMULSS.AND.AMGLSS.LT.AMURSS) THEN
            TMP(1)=MZIZ
            TMP(2)=AMULSS
            TMP(3)=AMULSS
            WIDL=2*AUIZ**2*(SSXINT(0.,SSGWQ1,UPPR)+SGNM3*SNIZ*
     $          SSXINT(0.,SSGWQ2,UPPR))
            TMP(2)=AMURSS
            TMP(3)=AMURSS
            WIDR=2*BUIZ**2*(SSXINT(0.,SSGWQ1,UPPR)+SGNM3*SNIZ*
     $          SSXINT(0.,SSGWQ2,UPPR))
            WID=WIDL+WIDR
          ELSE IF (AMGLSS.LT.AMULSS.AND.AMGLSS.GE.AMURSS) THEN
            TMP(2)=AMULSS
            TMP(3)=AMULSS
            WID=2*AUIZ**2*(SSXINT(0.,SSGWQ1,UPPR)+SGNM3*SNIZ*
     $          SSXINT(0.,SSGWQ2,UPPR))
          ELSE IF (AMGLSS.LT.AMURSS.AND.AMGLSS.GE.AMULSS) THEN
            TMP(2)=AMURSS
            TMP(3)=AMURSS
            WID=2*BUIZ**2*(SSXINT(0.,SSGWQ1,UPPR)+SGNM3*SNIZ*
     $          SSXINT(0.,SSGWQ2,UPPR))
          ELSE
            WID=0.
          END IF
          WID=GS**2/AMGLSS/2./(2*PI)**5*WID
          IF (WID.GT.0.) THEN
            CALL SSSAVE(ISGL,WID,ISZI(IZ),IDUP,-IDUP,0,0)
C           Enter decay matrix element info
            Z1(1)=1.
            Z1(2)=-Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*(-1.)*(THIZ+1)*AUIZ)
            Z2(2)=Z2(1)
            CALL SSME3(2,AMULSS,Z1,Z2)
            Z1(1)=1.
            Z1(2)=Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*BUIZ)
            Z2(2)=-Z2(1)
            CALL SSME3(2,AMURSS,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*(-1.)*(THIZ+1)*AUIZ
            Z1(2)=-Z1(1)
            Z2(1)=1.
            Z2(2)=Z2(1)
            CALL SSME3(3,AMULSS,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*BUIZ
            Z1(2)=Z1(1)
            Z2(1)=1.
            Z2(2)=-Z2(1)
            CALL SSME3(3,AMURSS,Z1,Z2)
          END IF
        END IF
C
        IF (AMGLSS.GT.(MZIZ+2*AMDN)) THEN
          IF (AMGLSS.LT.AMDLSS.AND.AMGLSS.LT.AMDRSS) THEN
            TMP(1)=MZIZ
            TMP(2)=AMDLSS
            TMP(3)=AMDLSS
            WIDL=2*ADIZ**2*(SSXINT(0.,SSGWQ1,UPPR)+SGNM3*SNIZ*
     $          SSXINT(0.,SSGWQ2,UPPR))
            TMP(2)=AMDRSS
            TMP(3)=AMDRSS
            WIDR=2*BDIZ**2*(SSXINT(0.,SSGWQ1,UPPR)+SGNM3*SNIZ*
     $          SSXINT(0.,SSGWQ2,UPPR))
            WID=WIDL+WIDR
          ELSE IF (AMGLSS.LT.AMDLSS.AND.AMGLSS.GE.AMDRSS) THEN
            TMP(2)=AMDLSS
            TMP(3)=AMDLSS
            WID=2*ADIZ**2*(SSXINT(0.,SSGWQ1,UPPR)+SGNM3*SNIZ*
     $          SSXINT(0.,SSGWQ2,UPPR))
          ELSE IF (AMGLSS.LT.AMDRSS.AND.AMGLSS.GE.AMDLSS) THEN
            TMP(2)=AMDRSS
            TMP(3)=AMDRSS
            WID=2*BDIZ**2*(SSXINT(0.,SSGWQ1,UPPR)+SGNM3*SNIZ*
     $          SSXINT(0.,SSGWQ2,UPPR))
          ELSE
            WID=0.
          END IF
          WID=GS**2/AMGLSS/2./(2*PI)**5*WID
          IF (WID.GT.0.) THEN
          CALL SSSAVE(ISGL,WID,ISZI(IZ),IDDN,-IDDN,0,0)
C           Enter decay matrix element info
            Z1(1)=1.
            Z1(2)=-Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*(-1.)*(THIZ+1)*ADIZ)
            Z2(2)=Z2(1)
            CALL SSME3(2,AMDLSS,Z1,Z2)
            Z1(1)=1.
            Z1(2)=Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*BDIZ)
            Z2(2)=-Z2(1)
            CALL SSME3(2,AMDRSS,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*(-1.)*(THIZ+1)*ADIZ
            Z1(2)=-Z1(1)
            Z2(1)=1.
            Z2(2)=Z2(1)
            CALL SSME3(3,AMDLSS,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*BDIZ
            Z1(2)=Z1(1)
            Z2(1)=1.
            Z2(2)=-Z2(1)
            CALL SSME3(3,AMDRSS,Z1,Z2)
          END IF
        END IF
C          3 body gluino --> q + qb + zi decay, q=s
        IF (AMGLSS.GT.(MZIZ+2*AMST)) THEN
          IF (AMGLSS.LT.AMSLSS.AND.AMGLSS.LT.AMSRSS) THEN
            TMP(1)=MZIZ
            TMP(2)=AMSLSS
            TMP(3)=AMSLSS
            WIDL=2*ADIZ**2*(SSXINT(0.,SSGWQ1,UPPR)+SGNM3*SNIZ*
     $          SSXINT(0.,SSGWQ2,UPPR))
            TMP(2)=AMSRSS
            TMP(3)=AMSRSS
            WIDR=2*BDIZ**2*(SSXINT(0.,SSGWQ1,UPPR)+SGNM3*SNIZ*
     $          SSXINT(0.,SSGWQ2,UPPR))
            WID=WIDL+WIDR
          ELSE IF (AMGLSS.LT.AMSLSS.AND.AMGLSS.GE.AMSRSS) THEN
            TMP(2)=AMSLSS
            TMP(3)=AMSLSS
            WID=2*ADIZ**2*(SSXINT(0.,SSGWQ1,UPPR)+SGNM3*SNIZ*
     $          SSXINT(0.,SSGWQ2,UPPR))
          ELSE IF (AMGLSS.LT.AMSRSS.AND.AMGLSS.GE.AMSLSS) THEN
            TMP(2)=AMSRSS
            TMP(3)=AMSRSS
            WID=2*BDIZ**2*(SSXINT(0.,SSGWQ1,UPPR)+SGNM3*SNIZ*
     $          SSXINT(0.,SSGWQ2,UPPR))
          ELSE
            WID=0.
          END IF
          WID=GS**2/AMGLSS/2./(2*PI)**5*WID
          IF (WID.GT.0.) THEN
          CALL SSSAVE(ISGL,WID,ISZI(IZ),IDST,-IDST,0,0)
C           Enter decay matrix element info
            Z1(1)=1.
            Z1(2)=-Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*(-1.)*(THIZ+1)*ADIZ)
            Z2(2)=Z2(1)
            CALL SSME3(2,AMDLSS,Z1,Z2)
            Z1(1)=1.
            Z1(2)=Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*BDIZ)
            Z2(2)=-Z2(1)
            CALL SSME3(2,AMDRSS,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*(-1.)*(THIZ+1)*ADIZ
            Z1(2)=-Z1(1)
            Z2(1)=1.
            Z2(2)=Z2(1)
            CALL SSME3(3,AMDLSS,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*BDIZ
            Z1(2)=Z1(1)
            Z2(1)=1.
            Z2(2)=-Z2(1)
            CALL SSME3(3,AMDRSS,Z1,Z2)
          END IF
        END IF
C          3 body gluino --> q + qb + zi decay, q=c
        IF (AMGLSS.GT.(MZIZ+2*AMCH)) THEN
          IF (AMGLSS.LT.AMCLSS.AND.AMGLSS.LT.AMCRSS) THEN
            TMP(1)=MZIZ
            TMP(2)=AMCLSS
            TMP(3)=AMCLSS
            WIDL=2*AUIZ**2*(SSXINT(0.,SSGWQ1,UPPR)+SGNM3*SNIZ*
     $          SSXINT(0.,SSGWQ2,UPPR))
            TMP(2)=AMCRSS
            TMP(3)=AMCRSS
            WIDR=2*BUIZ**2*(SSXINT(0.,SSGWQ1,UPPR)+SGNM3*SNIZ*
     $          SSXINT(0.,SSGWQ2,UPPR))
            WID=WIDL+WIDR
          ELSE IF (AMGLSS.LT.AMCLSS.AND.AMGLSS.GE.AMCRSS) THEN
            TMP(2)=AMCLSS
            TMP(3)=AMCLSS
            WID=2*AUIZ**2*(SSXINT(0.,SSGWQ1,UPPR)+SGNM3*SNIZ*
     $          SSXINT(0.,SSGWQ2,UPPR))
          ELSE IF (AMGLSS.LT.AMCRSS.AND.AMGLSS.GE.AMCLSS) THEN
            TMP(2)=AMCRSS
            TMP(3)=AMCRSS
            WID=2*BUIZ**2*(SSXINT(0.,SSGWQ1,UPPR)+SGNM3*SNIZ*
     $          SSXINT(0.,SSGWQ2,UPPR))
          ELSE
            WID=0.
          END IF
          WID=GS**2/AMGLSS/2./(2*PI)**5*WID
          IF (WID.GT.0.) THEN
          CALL SSSAVE(ISGL,WID,ISZI(IZ),IDCH,-IDCH,0,0)
C           Enter decay matrix element info
            Z1(1)=1.
            Z1(2)=-Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*(-1.)*(THIZ+1)*AUIZ)
            Z2(2)=Z2(1)
            CALL SSME3(2,AMULSS,Z1,Z2)
            Z1(1)=1.
            Z1(2)=Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*BUIZ)
            Z2(2)=-Z2(1)
            CALL SSME3(2,AMURSS,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*(-1.)*(THIZ+1)*AUIZ
            Z1(2)=-Z1(1)
            Z2(1)=1.
            Z2(2)=Z2(1)
            CALL SSME3(3,AMULSS,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*BUIZ
            Z1(2)=Z1(1)
            Z2(1)=1.
            Z2(2)=-Z2(1)
            CALL SSME3(3,AMURSS,Z1,Z2)
          END IF
        END IF
C          3 body gluino --> q + qb + zi decay, q=b 
        XTCON=ALF3/8./PI**4/AMGLSS
        IF (AMGLSS.GT.(MZIZ+2*AMBT).AND.AMGLSS.LT.
     $     (AMB1SS+AMBT)) THEN
          TMP(1)=AMGLSS
          TMP(2)=AMBT
          TMP(3)=MZIZ
          TMP(4)=AMB1SS
          TMP(5)=AMB1SS
          XUPPER=(AMGLSS**2+AMBT**2-(AMBT+MZIZ)**2)/2./AMGLSS
          XT1=SSXINT(AMBT,SSGX1,XUPPER)
          XT2=SSXINT(AMBT,SSGX2,XUPPER)
          XT3=SSXINT(AMBT,SSGX3,XUPPER)
          XT4=SSXINT(AMBT,SSGX4,XUPPER)
          XT5=SSXINT(AMBT,SSGX5,XUPPER)
          XT6=SSXINT(AMBT,SSGX6,XUPPER)
          XT7=SSXINT(AMBT,SSGX7,XUPPER)
          XT8=SSXINT(AMBT,SSGX8,XUPPER)
          XT9=SSXINT(AMBT,SSGX9,XUPPER)
          AL1=ADIZ*COSB-FB*ZMIXSS(2,IZ)*SINB
          BE1=FB*ZMIXSS(2,IZ)*COSB+BDIZ*SINB
C         ---- here, al2 is (-) al2 of tata notes-----
          AL2=BDIZ*SINB+FB*ZMIXSS(2,IZ)*COSB
          BE2=-FB*ZMIXSS(2,IZ)*SINB+ADIZ*COSB
          XLL=(AL1**2+BE1**2)*XT1-4*AMBT*MZIZ*SNIZ*AL1*
     $     BE1*XT3+SGNM3*AMGLSS*(SNIZ*MZIZ*(AL1**2*XT2/AMGLSS/
     $     MZIZ+BE1**2*AMBT**2*XT5)-AL1*BE1*AMBT*(XT4-
     $     MZIZ**2*XT5))
          XRR=(AL2**2+BE2**2)*XT1-4*AMBT*MZIZ*SNIZ*AL2*
     $     BE2*XT3+SGNM3*AMGLSS*(SNIZ*MZIZ*(AL2**2*XT2/AMGLSS/
     $     MZIZ+BE2**2*AMBT**2*XT5)-AL2*BE2*AMBT*(XT4-
     $     MZIZ**2*XT5))
          XL1R1=SGNM3*2*AMGLSS*AMBT*((AL1*AL2+BE1*BE2)*SNIZ*AMBT*
     $     MZIZ*XT6-(AL2*BE1+AL1*BE2)*XT7)
          XL2R2=XL1R1
          XL1R2=BE1*BE2*XT8+AL1*AL2*AMBT**2*XT4-AMBT*MZIZ*
     $     SNIZ*(AL1*BE2+AL2*BE1)*XT9
          XL2R1=XL1R2
          XMSB1=COSB**2*XLL+SINB**2*XRR-SINB*COSB*(XL1R1+XL1R2+
     $     XL2R1+XL2R2)
        ELSE
          XMSB1=0.
        END IF
        IF (AMGLSS.GT.(MZIZ+2*AMBT).AND.AMGLSS.LT.
     $     (AMB2SS+AMBT)) THEN
          TMP(1)=AMGLSS
          TMP(2)=AMBT
          TMP(3)=MZIZ
          TMP(4)=AMB2SS
          TMP(5)=AMB2SS
          XUPPER=(AMGLSS**2+AMBT**2-(AMBT+MZIZ)**2)/2./AMGLSS
          XT1=SSXINT(AMBT,SSGX1,XUPPER)
          XT2=SSXINT(AMBT,SSGX2,XUPPER)
          XT3=SSXINT(AMBT,SSGX3,XUPPER)
          XT4=SSXINT(AMBT,SSGX4,XUPPER)
          XT5=SSXINT(AMBT,SSGX5,XUPPER)
          XT6=SSXINT(AMBT,SSGX6,XUPPER)
          XT7=SSXINT(AMBT,SSGX7,XUPPER)
          XT8=SSXINT(AMBT,SSGX8,XUPPER)
          XT9=SSXINT(AMBT,SSGX9,XUPPER)
          AL1=ADIZ*SINB+FB*ZMIXSS(2,IZ)*COSB
          BE1=FB*ZMIXSS(2,IZ)*SINB-BDIZ*COSB
C         ---- here, al2 is (-) al2 of tata notes-----
          AL2=-BDIZ*COSB+FB*ZMIXSS(2,IZ)*SINB
          BE2=FB*ZMIXSS(2,IZ)*COSB+ADIZ*SINB
          XLL=(AL1**2+BE1**2)*XT1-4*AMBT*MZIZ*SNIZ*AL1*
     $     BE1*XT3+SGNM3*AMGLSS*(SNIZ*MZIZ*(AL1**2*XT2/AMGLSS/
     $     MZIZ+BE1**2*AMBT**2*XT5)-AL1*BE1*AMBT*(XT4-
     $     MZIZ**2*XT5))
          XRR=(AL2**2+BE2**2)*XT1-4*AMBT*MZIZ*SNIZ*AL2*
     $     BE2*XT3+SGNM3*AMGLSS*(SNIZ*MZIZ*(AL2**2*XT2/AMGLSS/
     $     MZIZ+BE2**2*AMBT**2*XT5)-AL2*BE2*AMBT*(XT4-
     $     MZIZ**2*XT5))
          XL1R1=SGNM3*2*AMGLSS*AMBT*((AL1*AL2+BE1*BE2)*SNIZ*AMBT*
     $     MZIZ*XT6-(AL2*BE1+AL1*BE2)*XT7)
          XL2R2=XL1R1
          XL1R2=BE1*BE2*XT8+AL1*AL2*AMBT**2*XT4-AMBT*MZIZ*
     $     SNIZ*(AL1*BE2+AL2*BE1)*XT9
          XL2R1=XL1R2
          XMSB2=SINB**2*XLL+COSB**2*XRR+SINB*COSB*(XL1R1+XL1R2+
     $     XL2R1+XL2R2)
        ELSE
          XMSB2=0.
        END IF
C       ----cross term between b_1 and b_2 graphs -----------
        IF (AMGLSS.GT.(MZIZ+2*AMBT).AND.AMGLSS.LT.
     $     (AMB1SS+AMBT)) THEN
          TMP(1)=AMGLSS
          TMP(2)=AMBT
          TMP(3)=MZIZ
          TMP(4)=AMB1SS
          TMP(5)=AMB2SS
          XUPPER=(AMGLSS**2+AMBT**2-(AMBT+MZIZ)**2)/2./AMGLSS
          XT1=SSXINT(AMBT,SSGX1,XUPPER)
          XT2=SSXINT(AMBT,SSGX2,XUPPER)
          XT3=SSXINT(AMBT,SSGX3,XUPPER)
          XT4=SSXINT(AMBT,SSGX4,XUPPER)
          XT5=SSXINT(AMBT,SSGX5,XUPPER)
          XT6=SSXINT(AMBT,SSGX6,XUPPER)
          XT7=SSXINT(AMBT,SSGX7,XUPPER)
          XT8=SSXINT(AMBT,SSGX8,XUPPER)
          XT9=SSXINT(AMBT,SSGX9,XUPPER)
          AL1=ADIZ*COSB-FB*ZMIXSS(2,IZ)*SINB
          AL2=ADIZ*SINB+FB*ZMIXSS(2,IZ)*COSB
          BE1=FB*ZMIXSS(2,IZ)*COSB+BDIZ*SINB
          BE2=FB*ZMIXSS(2,IZ)*SINB-BDIZ*COSB
          XL1L2=COSB*SINB*(2*(AL1*AL2+BE1*BE2)*XT1-4*SNIZ*AMBT*
     $     MZIZ*(AL1*BE2+AL2*BE1)*XT3+SGNM3*AMGLSS*(2*MZIZ*
     $     SNIZ*(AL1*AL2*XT2/AMGLSS/MZIZ+BE1*BE2*AMBT**2*
     $     XT5)-(AL1*BE2+AL2*BE1)*AMBT*(XT4-MZIZ**2*XT5)))
          AL1=-BDIZ*SINB-FB*ZMIXSS(2,IZ)*COSB
          AL2=BDIZ*COSB-FB*ZMIXSS(2,IZ)*SINB
          BE1=-FB*ZMIXSS(2,IZ)*SINB+ADIZ*COSB
          BE2=FB*ZMIXSS(2,IZ)*COSB+ADIZ*SINB
          XR1R2=-COSB*SINB*(2*(AL1*AL2+BE1*BE2)*XT1+4*SNIZ*AMBT*
     $     MZIZ*(AL1*BE2+AL2*BE1)*XT3+SGNM3*AMGLSS*(2*MZIZ*
     $     SNIZ*(AL1*AL2*XT2/AMGLSS/MZIZ+BE1*BE2*AMBT**2*
     $     XT5)+(AL1*BE2+AL2*BE1)*AMBT*(XT4-MZIZ**2*XT5)))
          AL1=ADIZ*COSB-FB*ZMIXSS(2,IZ)*SINB
          AL2=BDIZ*COSB-FB*ZMIXSS(2,IZ)*SINB
          BE1=FB*ZMIXSS(2,IZ)*COSB+BDIZ*SINB
          BE2=FB*ZMIXSS(2,IZ)*COSB+ADIZ*SINB
          XL1R1=-SGNM3*2*AMGLSS*AMBT*COSB**2*(SNIZ*(AL1*AL2-BE1*BE2)*
     $     AMBT*MZIZ*XT6-(AL2*BE1-AL1*BE2)*XT7)
          XL1R2=COSB**2*(BE1*BE2*XT8-AL1*AL2*AMBT**2*XT4+AMBT*
     $     MZIZ*SNIZ*XT9*(-AL1*BE2+BE1*AL2))
          XLR1=2*(XL1R1+XL1R2)
          AL1=ADIZ*SINB+FB*ZMIXSS(2,IZ)*COSB
          AL2=-BDIZ*SINB-FB*ZMIXSS(2,IZ)*COSB
          BE1=FB*ZMIXSS(2,IZ)*SINB-BDIZ*COSB
          BE2=-FB*ZMIXSS(2,IZ)*SINB+ADIZ*COSB
          TMP(4)=AMB2SS
          TMP(5)=AMB1SS
          XT8=SSXINT(AMBT,SSGX8,XUPPER)
          XT9=SSXINT(AMBT,SSGX9,XUPPER)
          XL1R1=-SGNM3*2*AMGLSS*AMBT*SINB**2*(SNIZ*(-AL1*AL2+BE1*BE2)*
     $     AMBT*MZIZ*XT6+(AL2*BE1-AL1*BE2)*XT7)
          XL1R2=-SINB**2*(BE1*BE2*XT8-AL1*AL2*AMBT**2*XT4+AMBT*
     $     MZIZ*SNIZ*XT9*(-AL1*BE2+BE1*AL2))
          XLR2=2*(XL1R1+XL1R2)
          XMSB12=XL1L2+XR1R2+XLR1+XLR2
        ELSE
          XMSB12=0.
        END IF
        WID=XTCON*(XMSB1+XMSB2+XMSB12)
        IF (WID.GT.0.) THEN
          CALL SSSAVE(ISGL,WID,ISZI(IZ),IDBT,-IDBT,0,0)
          Z1(1)=((ZI*ZADIZ-FB*ZMIXSS(2,IZ)*ZI**THIZ)*COSB-
     $           (ZI*ZBDIZ-FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*SINB)/2.
          Z1(2)=((-ZI*ZADIZ-FB*ZMIXSS(2,IZ)*ZI**THIZ)*COSB-
     $           (ZI*ZBDIZ+FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*SINB)/2.
          Z2(1)=(COSB-SINB)/2.
          Z2(2)=-(COSB+SINB)/2.
          CALL SSME3(3,AMB1SS,Z1,Z2)
          Z1(1)=(COSB-SINB)/2.
          Z1(2)=-(COSB+SINB)/2.
          Z2(1)=CONJG((ZI*ZADIZ-FB*ZMIXSS(2,IZ)*ZI**THIZ)*COSB-
     $           (ZI*ZBDIZ-FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*SINB)/2.
          Z2(2)=-CONJG((-ZI*ZADIZ-FB*ZMIXSS(2,IZ)*ZI**THIZ)*COSB-
     $           (ZI*ZBDIZ+FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*SINB)/2.
          CALL SSME3(2,AMB1SS,Z1,Z2)
          Z1(1)=((ZI*ZADIZ-FB*ZMIXSS(2,IZ)*ZI**THIZ)*SINB+
     $           (ZI*ZBDIZ-FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*COSB)/2.
          Z1(2)=((-ZI*ZADIZ-FB*ZMIXSS(2,IZ)*ZI**THIZ)*SINB+
     $           (ZI*ZBDIZ+FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*COSB)/2.
          Z2(1)=(COSB+SINB)/2.
          Z2(2)=(COSB-SINB)/2.
          CALL SSME3(3,AMB2SS,Z1,Z2)
          Z1(1)=(COSB+SINB)/2.
          Z1(2)=(COSB-SINB)/2.
          Z2(1)=CONJG((ZI*ZADIZ-FB*ZMIXSS(2,IZ)*ZI**THIZ)*SINB+
     $           (ZI*ZBDIZ-FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*COSB)/2.
          Z2(2)=-CONJG((-ZI*ZADIZ-FB*ZMIXSS(2,IZ)*ZI**THIZ)*SINB+
     $           (ZI*ZBDIZ+FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*COSB)/2.
          CALL SSME3(2,AMB2SS,Z1,Z2)
        END IF
C          3 body gluino --> q + qb + zi decay, q=t
        IF (AMGLSS.GT.(MZIZ+2*AMTP).AND.AMGLSS.LT.
     $     (AMT1SS+AMTP)) THEN
          TMP(1)=AMGLSS
          TMP(2)=AMTP
          TMP(3)=MZIZ
          TMP(4)=AMT1SS
          TMP(5)=AMT1SS
          XUPPER=(AMGLSS**2+AMTP**2-(AMTP+MZIZ)**2)/2./AMGLSS
          XT1=SSXINT(AMTP,SSGX1,XUPPER)
          XT2=SSXINT(AMTP,SSGX2,XUPPER)
          XT3=SSXINT(AMTP,SSGX3,XUPPER)
          XT4=SSXINT(AMTP,SSGX4,XUPPER)
          XT5=SSXINT(AMTP,SSGX5,XUPPER)
          XT6=SSXINT(AMTP,SSGX6,XUPPER)
          XT7=SSXINT(AMTP,SSGX7,XUPPER)
          XT8=SSXINT(AMTP,SSGX8,XUPPER)
          XT9=SSXINT(AMTP,SSGX9,XUPPER)
          AL1=AUIZ*COST-FT*ZMIXSS(1,IZ)*SINT
          BE1=FT*ZMIXSS(1,IZ)*COST+BUIZ*SINT
C         ---- here, al2 is (-) al2 of tata notes-----
          AL2=BUIZ*SINT+FT*ZMIXSS(1,IZ)*COST
          BE2=-FT*ZMIXSS(1,IZ)*SINT+AUIZ*COST
          XLL=(AL1**2+BE1**2)*XT1-4*AMTP*MZIZ*SNIZ*AL1*
     $     BE1*XT3+SGNM3*AMGLSS*(SNIZ*MZIZ*(AL1**2*XT2/AMGLSS/
     $     MZIZ+BE1**2*AMTP**2*XT5)-AL1*BE1*AMTP*(XT4-
     $     MZIZ**2*XT5))
          XRR=(AL2**2+BE2**2)*XT1-4*AMTP*MZIZ*SNIZ*AL2*
     $     BE2*XT3+SGNM3*AMGLSS*(SNIZ*MZIZ*(AL2**2*XT2/AMGLSS/
     $     MZIZ+BE2**2*AMTP**2*XT5)-AL2*BE2*AMTP*(XT4-
     $     MZIZ**2*XT5))
          XL1R1=SGNM3*2*AMGLSS*AMTP*((AL1*AL2+BE1*BE2)*SNIZ*AMTP*
     $     MZIZ*XT6-(AL2*BE1+AL1*BE2)*XT7)
          XL2R2=XL1R1
          XL1R2=BE1*BE2*XT8+AL1*AL2*AMTP**2*XT4-AMTP*MZIZ*
     $     SNIZ*(AL1*BE2+AL2*BE1)*XT9
          XL2R1=XL1R2
          XMST1=COST**2*XLL+SINT**2*XRR-SINT*COST*(XL1R1+XL1R2+
     $     XL2R1+XL2R2)
        ELSE
          XMST1=0.
        END IF
        IF (AMGLSS.GT.(MZIZ+2*AMTP).AND.AMGLSS.LT.
     $     (AMT2SS+AMTP)) THEN
          TMP(1)=AMGLSS
          TMP(2)=AMTP
          TMP(3)=MZIZ
          TMP(4)=AMT2SS
          TMP(5)=AMT2SS
          XUPPER=(AMGLSS**2+AMTP**2-(AMTP+MZIZ)**2)/2./AMGLSS
          XT1=SSXINT(AMTP,SSGX1,XUPPER)
          XT2=SSXINT(AMTP,SSGX2,XUPPER)
          XT3=SSXINT(AMTP,SSGX3,XUPPER)
          XT4=SSXINT(AMTP,SSGX4,XUPPER)
          XT5=SSXINT(AMTP,SSGX5,XUPPER)
          XT6=SSXINT(AMTP,SSGX6,XUPPER)
          XT7=SSXINT(AMTP,SSGX7,XUPPER)
          XT8=SSXINT(AMTP,SSGX8,XUPPER)
          XT9=SSXINT(AMTP,SSGX9,XUPPER)
          AL1=AUIZ*SINT+FT*ZMIXSS(1,IZ)*COST
          BE1=FT*ZMIXSS(1,IZ)*SINT-BUIZ*COST
C         ---- here, al2 is (-) al2 of tata notes-----
          AL2=-BUIZ*COST+FT*ZMIXSS(1,IZ)*SINT
          BE2=FT*ZMIXSS(1,IZ)*COST+AUIZ*SINT
          XLL=(AL1**2+BE1**2)*XT1-4*AMTP*MZIZ*SNIZ*AL1*
     $     BE1*XT3+SGNM3*AMGLSS*(SNIZ*MZIZ*(AL1**2*XT2/AMGLSS/
     $     MZIZ+BE1**2*AMTP**2*XT5)-AL1*BE1*AMTP*(XT4-
     $     MZIZ**2*XT5))
          XRR=(AL2**2+BE2**2)*XT1-4*AMTP*MZIZ*SNIZ*AL2*
     $     BE2*XT3+SGNM3*AMGLSS*(SNIZ*MZIZ*(AL2**2*XT2/AMGLSS/
     $     MZIZ+BE2**2*AMTP**2*XT5)-AL2*BE2*AMTP*(XT4-
     $     MZIZ**2*XT5))
          XL1R1=SGNM3*2*AMGLSS*AMTP*((AL1*AL2+BE1*BE2)*SNIZ*AMTP*
     $     MZIZ*XT6-(AL2*BE1+AL1*BE2)*XT7)
          XL2R2=XL1R1
          XL1R2=BE1*BE2*XT8+AL1*AL2*AMTP**2*XT4-AMTP*MZIZ*
     $     SNIZ*(AL1*BE2+AL2*BE1)*XT9
          XL2R1=XL1R2
          XMST2=SINT**2*XLL+COST**2*XRR+SINT*COST*(XL1R1+XL1R2+
     $     XL2R1+XL2R2)
        ELSE
          XMST2=0.
        END IF
C       ----cross term between t_1 and t_2 graphs -----------
        IF (AMGLSS.GT.(MZIZ+2*AMTP).AND.AMGLSS.LT.
     $     (AMT1SS+AMTP)) THEN
          TMP(1)=AMGLSS
          TMP(2)=AMTP
          TMP(3)=MZIZ
          TMP(4)=AMT1SS
          TMP(5)=AMT2SS
          XUPPER=(AMGLSS**2+AMTP**2-(AMTP+MZIZ)**2)/2./AMGLSS
          XT1=SSXINT(AMTP,SSGX1,XUPPER)
          XT2=SSXINT(AMTP,SSGX2,XUPPER)
          XT3=SSXINT(AMTP,SSGX3,XUPPER)
          XT4=SSXINT(AMTP,SSGX4,XUPPER)
          XT5=SSXINT(AMTP,SSGX5,XUPPER)
          XT6=SSXINT(AMTP,SSGX6,XUPPER)
          XT7=SSXINT(AMTP,SSGX7,XUPPER)
          XT8=SSXINT(AMTP,SSGX8,XUPPER)
          XT9=SSXINT(AMTP,SSGX9,XUPPER)
          AL1=AUIZ*COST-FT*ZMIXSS(1,IZ)*SINT
          AL2=AUIZ*SINT+FT*ZMIXSS(1,IZ)*COST
          BE1=FT*ZMIXSS(1,IZ)*COST+BUIZ*SINT
          BE2=FT*ZMIXSS(1,IZ)*SINT-BUIZ*COST
          XL1L2=COST*SINT*(2*(AL1*AL2+BE1*BE2)*XT1-4*SNIZ*AMTP*
     $     MZIZ*(AL1*BE2+AL2*BE1)*XT3+SGNM3*AMGLSS*(2*MZIZ*
     $     SNIZ*(AL1*AL2*XT2/AMGLSS/MZIZ+BE1*BE2*AMTP**2*
     $     XT5)-(AL1*BE2+AL2*BE1)*AMTP*(XT4-MZIZ**2*XT5)))
          AL1=-BUIZ*SINT-FT*ZMIXSS(1,IZ)*COST
          AL2=BUIZ*COST-FT*ZMIXSS(1,IZ)*SINT
          BE1=-FT*ZMIXSS(1,IZ)*SINT+AUIZ*COST
          BE2=FT*ZMIXSS(1,IZ)*COST+AUIZ*SINT
          XR1R2=-COST*SINT*(2*(AL1*AL2+BE1*BE2)*XT1+4*SNIZ*AMTP*
     $     MZIZ*(AL1*BE2+AL2*BE1)*XT3+SGNM3*AMGLSS*(2*MZIZ*
     $     SNIZ*(AL1*AL2*XT2/AMGLSS/MZIZ+BE1*BE2*AMTP**2*
     $     XT5)+(AL1*BE2+AL2*BE1)*AMTP*(XT4-MZIZ**2*XT5)))
          AL1=AUIZ*COST-FT*ZMIXSS(1,IZ)*SINT
          AL2=BUIZ*COST-FT*ZMIXSS(1,IZ)*SINT
          BE1=FT*ZMIXSS(1,IZ)*COST+BUIZ*SINT
          BE2=FT*ZMIXSS(1,IZ)*COST+AUIZ*SINT
          XL1R1=-SGNM3*2*AMGLSS*AMTP*COST**2*(SNIZ*(AL1*AL2-BE1*BE2)*
     $     AMTP*MZIZ*XT6-(AL2*BE1-AL1*BE2)*XT7)
          XL1R2=COST**2*(BE1*BE2*XT8-AL1*AL2*AMTP**2*XT4+AMTP*
     $     MZIZ*SNIZ*XT9*(-AL1*BE2+BE1*AL2))
          XLR1=2*(XL1R1+XL1R2)
          AL1=AUIZ*SINT+FT*ZMIXSS(1,IZ)*COST
          AL2=-BUIZ*SINT-FT*ZMIXSS(1,IZ)*COST
          BE1=FT*ZMIXSS(1,IZ)*SINT-BUIZ*COST
          BE2=-FT*ZMIXSS(1,IZ)*SINT+AUIZ*COST
          TMP(4)=AMT2SS
          TMP(5)=AMT1SS
          XT8=SSXINT(AMTP,SSGX8,XUPPER)
          XT9=SSXINT(AMTP,SSGX9,XUPPER)
          XL1R1=-SGNM3*2*AMGLSS*AMTP*SINT**2*(SNIZ*(-AL1*AL2+BE1*BE2)*
     $     AMTP*MZIZ*XT6+(AL2*BE1-AL1*BE2)*XT7)
          XL1R2=-SINT**2*(BE1*BE2*XT8-AL1*AL2*AMTP**2*XT4+AMTP*
     $     MZIZ*SNIZ*XT9*(-AL1*BE2+BE1*AL2))
          XLR2=2*(XL1R1+XL1R2)
          XMST12=XL1L2+XR1R2+XLR1+XLR2
        ELSE
          XMST12=0.
        END IF
          WID=XTCON*(XMST1+XMST2+XMST12)
        IF (WID.GT.0.) THEN
          CALL SSSAVE(ISGL,WID,ISZI(IZ),IDTP,-IDTP,0,0)
          Z1(1)=((ZI*ZAUIZ-FT*ZMIXSS(1,IZ)*ZI**THIZ)*COST-
     $           (ZI*ZBUIZ-FT*ZMIXSS(1,IZ)*(-ZI)**THIZ)*SINT)/2.
          Z1(2)=((-ZI*ZAUIZ-FT*ZMIXSS(1,IZ)*ZI**THIZ)*COST-
     $           (ZI*ZBUIZ+FT*ZMIXSS(1,IZ)*(-ZI)**THIZ)*SINT)/2.
          Z2(1)=(COST-SINT)/2.
          Z2(2)=-(COST+SINT)/2.
          CALL SSME3(3,AMT1SS,Z1,Z2)
          Z1(1)=(COST-SINT)/2.
          Z1(2)=-(COST+SINT)/2.
          Z2(1)=CONJG((ZI*ZAUIZ-FT*ZMIXSS(1,IZ)*ZI**THIZ)*COST-
     $           (ZI*ZBUIZ-FT*ZMIXSS(1,IZ)*(-ZI)**THIZ)*SINT)/2.
          Z2(2)=-CONJG((-ZI*ZAUIZ-FT*ZMIXSS(1,IZ)*ZI**THIZ)*COST-
     $           (ZI*ZBUIZ+FT*ZMIXSS(1,IZ)*(-ZI)**THIZ)*SINT)/2.
          CALL SSME3(2,AMT1SS,Z1,Z2)
          Z1(1)=((ZI*ZAUIZ-FT*ZMIXSS(1,IZ)*ZI**THIZ)*SINT+
     $           (ZI*ZBUIZ-FT*ZMIXSS(1,IZ)*(-ZI)**THIZ)*COST)/2.
          Z1(2)=((-ZI*ZAUIZ-FT*ZMIXSS(1,IZ)*ZI**THIZ)*SINT+
     $           (ZI*ZBUIZ+FT*ZMIXSS(1,IZ)*(-ZI)**THIZ)*COST)/2.
          Z2(1)=(COST+SINT)/2.
          Z2(2)=(COST-SINT)/2.
          CALL SSME3(3,AMT2SS,Z1,Z2)
          Z1(1)=(COST+SINT)/2.
          Z1(2)=(COST-SINT)/2.
          Z2(1)=CONJG((ZI*ZAUIZ-FT*ZMIXSS(1,IZ)*ZI**THIZ)*SINT+
     $           (ZI*ZBUIZ-FT*ZMIXSS(1,IZ)*(-ZI)**THIZ)*COST)/2.
          Z2(2)=-CONJG((-ZI*ZAUIZ-FT*ZMIXSS(1,IZ)*ZI**THIZ)*SINT+
     $           (ZI*ZBUIZ+FT*ZMIXSS(1,IZ)*(-ZI)**THIZ)*COST)/2.
          CALL SSME3(2,AMT2SS,Z1,Z2)
        END IF
100   CONTINUE
C
C          gluino --> quark + squark mode
C     
      IF (AMGLSS.GT.(AMULSS+AMUP)) THEN
        GMQK=ALF3*AMGLSS*(1.+AMUP**2/AMGLSS**2-AMULSS**2/AMGLSS**2)*
     $  SQRT(SSXLAM(1.,AMUP**2/AMGLSS**2,AMULSS**2/AMGLSS**2))/8.
        CALL SSSAVE(ISGL,GMQK,-ISUPL,+IDUP,0,0,0)
        CALL SSSAVE(ISGL,GMQK,+ISUPL,-IDUP,0,0,0)
      END IF
      IF (AMGLSS.GT.(AMDLSS+AMDN)) THEN
        GMQK=ALF3*AMGLSS*(1.+AMDN**2/AMGLSS**2-AMDLSS**2/AMGLSS**2)*
     $  SQRT(SSXLAM(1.,AMDN**2/AMGLSS**2,AMDLSS**2/AMGLSS**2))/8.
        CALL SSSAVE(ISGL,GMQK,-ISDNL,+IDDN,0,0,0)
        CALL SSSAVE(ISGL,GMQK,+ISDNL,-IDDN,0,0,0)
      END IF
      IF (AMGLSS.GT.(AMURSS+AMUP)) THEN
        GMQK=ALF3*AMGLSS*(1.+AMUP**2/AMGLSS**2-AMURSS**2/AMGLSS**2)*
     $  SQRT(SSXLAM(1.,AMUP**2/AMGLSS**2,AMURSS**2/AMGLSS**2))/8.
        CALL SSSAVE(ISGL,GMQK,-ISUPR,+IDUP,0,0,0)
        CALL SSSAVE(ISGL,GMQK,+ISUPR,-IDUP,0,0,0)
      END IF
      IF (AMGLSS.GT.(AMDRSS+AMDN)) THEN
        GMQK=ALF3*AMGLSS*(1.+AMDN**2/AMGLSS**2-AMDRSS**2/AMGLSS**2)*
     $  SQRT(SSXLAM(1.,AMDN**2/AMGLSS**2,AMDRSS**2/AMGLSS**2))/8.
        CALL SSSAVE(ISGL,GMQK,-ISDNR,+IDDN,0,0,0)
        CALL SSSAVE(ISGL,GMQK,+ISDNR,-IDDN,0,0,0)
      END IF
C
      IF (AMGLSS.GT.(AMSLSS+AMST)) THEN
        GMQK=ALF3*AMGLSS*(1.+AMST**2/AMGLSS**2-AMSLSS**2/AMGLSS**2)*
     $  SQRT(SSXLAM(1.,AMST**2/AMGLSS**2,AMSLSS**2/AMGLSS**2))/8.
        CALL SSSAVE(ISGL,GMQK,-ISSTL,+IDST,0,0,0)
        CALL SSSAVE(ISGL,GMQK,+ISSTL,-IDST,0,0,0)
      END IF
      IF (AMGLSS.GT.(AMSRSS+AMST)) THEN
        GMQK=ALF3*AMGLSS*(1.+AMST**2/AMGLSS**2-AMSRSS**2/AMGLSS**2)*
     $  SQRT(SSXLAM(1.,AMST**2/AMGLSS**2,AMSRSS**2/AMGLSS**2))/8.
        CALL SSSAVE(ISGL,GMQK,-ISSTR,+IDST,0,0,0)
        CALL SSSAVE(ISGL,GMQK,+ISSTR,-IDST,0,0,0)
      END IF
C     
      IF (AMGLSS.GT.(AMCLSS+AMCH)) THEN
        GMQK=ALF3*AMGLSS*(1.+AMCH**2/AMGLSS**2-AMCLSS**2/AMGLSS**2)*
     $  SQRT(SSXLAM(1.,AMCH**2/AMGLSS**2,AMCLSS**2/AMGLSS**2))/8.
        CALL SSSAVE(ISGL,GMQK,-ISCHL,+IDCH,0,0,0)
        CALL SSSAVE(ISGL,GMQK,+ISCHL,-IDCH,0,0,0)
      END IF
      IF (AMGLSS.GT.(AMCRSS+AMCH)) THEN
        GMQK=ALF3*AMGLSS*(1.+AMCH**2/AMGLSS**2-AMCRSS**2/AMGLSS**2)*
     $  SQRT(SSXLAM(1.,AMCH**2/AMGLSS**2,AMCRSS**2/AMGLSS**2))/8.
        CALL SSSAVE(ISGL,GMQK,-ISCHR,+IDCH,0,0,0)
        CALL SSSAVE(ISGL,GMQK,+ISCHR,-IDCH,0,0,0)
      END IF
C     
      IF (AMGLSS.GT.(AMB1SS+AMBT)) THEN
        GMQK=ALF3*AMGLSS*(1.+AMBT**2/AMGLSS**2-AMB1SS**2/AMGLSS**2+
     $   SGNM3*2*SIN(2*THETAB)*AMBT/AMGLSS)*
     $  SQRT(SSXLAM(1.,AMBT**2/AMGLSS**2,AMB1SS**2/AMGLSS**2))/8.
        CALL SSSAVE(ISGL,GMQK,-ISBT1,+IDBT,0,0,0)
        CALL SSSAVE(ISGL,GMQK,+ISBT1,-IDBT,0,0,0)
      END IF
C     
      IF (AMGLSS.GT.(AMB2SS+AMBT)) THEN
        GMQK=ALF3*AMGLSS*(1.+AMBT**2/AMGLSS**2-AMB2SS**2/AMGLSS**2-
     $   SGNM3*2*SIN(2*THETAB)*AMBT/AMGLSS)*
     $  SQRT(SSXLAM(1.,AMBT**2/AMGLSS**2,AMB2SS**2/AMGLSS**2))/8.
        CALL SSSAVE(ISGL,GMQK,-ISBT2,+IDBT,0,0,0)
        CALL SSSAVE(ISGL,GMQK,+ISBT2,-IDBT,0,0,0)
      END IF
C     
      IF (AMGLSS.GT.(AMT1SS+AMTP)) THEN
        GMQK=ALF3*AMGLSS*(1.+AMTP**2/AMGLSS**2-AMT1SS**2/AMGLSS**2+
     $   SGNM3*2*SIN(2*THETAT)*AMTP/AMGLSS)*
     $  SQRT(SSXLAM(1.,AMTP**2/AMGLSS**2,AMT1SS**2/AMGLSS**2))/8.
        CALL SSSAVE(ISGL,GMQK,-ISTP1,+IDTP,0,0,0)
        CALL SSSAVE(ISGL,GMQK,+ISTP1,-IDTP,0,0,0)
      END IF
C     
      IF (AMGLSS.GT.(AMT2SS+AMTP)) THEN
        GMQK=ALF3*AMGLSS*(1.+AMTP**2/AMGLSS**2-AMT2SS**2/AMGLSS**2-
     $   SGNM3*2*SIN(2*THETAT)*AMTP/AMGLSS)*
     $  SQRT(SSXLAM(1.,AMTP**2/AMGLSS**2,AMT2SS**2/AMGLSS**2))/8.
        CALL SSSAVE(ISGL,GMQK,-ISTP2,+IDTP,0,0,0)
        CALL SSSAVE(ISGL,GMQK,+ISTP2,-IDTP,0,0,0)
      END IF
C
C     Decay to gravitino
C
      IF (AMGLSS.GT.AMGVSS) THEN
        WID=AMGLSS**5/48./PI/(AMGVSS*AMPL)**2
        CALL SSSAVE(ISGL,WID,91,IDGL,0,0,0)
      END IF
C
C          Normalize branching ratios
C
      CALL SSNORM(ISGL)
C
      RETURN
      END
