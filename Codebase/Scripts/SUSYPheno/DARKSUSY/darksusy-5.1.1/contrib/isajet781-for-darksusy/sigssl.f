CDECK  ID>, SIGSSL.
      SUBROUTINE SIGSSL
C
C          Calculate d(sigma)/d(pt**2)d(y1)d(y2) for supersymmetric
C          sleptons and sneutrinos in MSSM using cross
C          sections from Baer and Tata.
C
C          SIGMA    = cross section summed over types allowed by
C                     JETTYPE cards.
C          SIGS(I)  = partial cross section for I1 + I2 --> I3 + I4
C          INOUT(I) = IOPAK**3*I4 + IOPAK**2*I3 + IOPAK*I2 +I1
C          JETTYP -> IDENT mapping:
C          GLSS, UPSSL, UBSSL, ..., UPSSR, UBSSR, ...,
C          W1SS+, W1SS-, WS22+, W2SS-, Z1SS, Z2SS, Z3SS, Z4SS
C          NUEL, ANUEL, EL-, ..., TAUL+
C
C          Extra factor of 1/2 needed for nonidentical final jets.
C          Y=-log(tan(theta/2)) gives jacobean P1*P2/E1*E2
C
C          Called from SIGSSY and so does not reinitialize /JETSIG/.
C
C
      IMPLICIT NONE
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/CONST/PI,SQRT2,ALFA,GF,UNITS
      SAVE /CONST/
      REAL      PI,SQRT2,ALFA,GF,UNITS
      COMMON/JETPAR/P(3),PT(3),YJ(3),PHI(3),XJ(3),TH(3),CTH(3),STH(3)
     1 ,JETTYP(3),SHAT,THAT,UHAT,QSQ,X1,X2,PBEAM(2)
     2 ,QMW,QW,QTW,YW,XW,THW,QTMW,PHIW,SHAT1,THAT1,UHAT1,JWTYP
     3 ,ALFQSQ,CTHW,STHW,Q0W
     4 ,INITYP(2),ISIGS,PBEAMS(5)
      SAVE /JETPAR/
      INTEGER   JETTYP,JWTYP,INITYP,ISIGS
      REAL      P,PT,YJ,PHI,XJ,TH,CTH,STH,SHAT,THAT,UHAT,QSQ,X1,X2,
     +          PBEAM,QMW,QW,QTW,YW,XW,THW,QTMW,PHIW,SHAT1,THAT1,UHAT1,
     +          ALFQSQ,CTHW,STHW,Q0W,PBEAMS
      INTEGER   MXSIGS,IOPAK
      PARAMETER (MXSIGS=3000,IOPAK=100)
      COMMON/JETSIG/SIGMA,SIGS(MXSIGS),NSIGS,INOUT(MXSIGS),SIGEVT
      SAVE /JETSIG/
      INTEGER   NSIGS,INOUT
      REAL      SIGMA,SIGS,SIGEVT
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
      INTEGER MXGOQ,MXGOJ
      PARAMETER (MXGOQ=85,MXGOJ=8)
      COMMON/Q1Q2/GOQ(MXGOQ,MXGOJ),GOALL(MXGOJ),GODY(4),STDDY,
     $GOWW(25,2),ALLWW(2),GOWMOD(25,MXGOJ)
      SAVE /Q1Q2/
      LOGICAL GOQ,GOALL,GODY,STDDY,GOWW,ALLWW,GOWMOD
      COMMON/QCDPAR/ALAM,ALAM2,CUTJET,ISTRUC
      SAVE /QCDPAR/
      INTEGER   ISTRUC
      REAL      ALAM,ALAM2,CUTJET
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
      COMMON/WCON/SIN2W,WMASS(4),WGAM(4),AQ(12,4),BQ(12,4),COUT(4),
     1MATCH(25,4),WCBR(25,4),CUTOFF,CUTPOW,TBRWW(4,2),RBRWW(12,4,2),EZ,
     2AQDP(12,4),BQDP(12,4),EZDP,WFUDGE
      SAVE /WCON/
      DOUBLE PRECISION AQDP,BQDP,EZDP
      INTEGER   MATCH
      REAL      SIN2W,WMASS,WGAM,AQ,BQ,COUT,WCBR,CUTOFF,CUTPOW,TBRWW,
     +          RBRWW,EZ,WFUDGE
      COMMON/WCON2/CUMWBR(25,3)
      REAL CUMWBR
C
      REAL X(2)
      EQUIVALENCE (X(1),X1)
      EQUIVALENCE (S,SHAT),(T,THAT),(U,UHAT)
      INTEGER JS2JT(25),IW2JS(4),IW2IM(4),IZ2JS(4),IS2UD(25)
      SAVE JS2JT,IW2JS,IW2IM,IZ2JS,IS2UD
      INTEGER IDLSS(18)
      SAVE IDLSS
      INTEGER IL2JS(18),IS2LN(18),II
      SAVE IL2JS,IS2LN
      REAL SIG,S,T,U,FAC,AM22,AM12,TT,GP,G,
     $E1,E2
      INTEGER IQ,IQ1,IQ2,IH
      REAL QFCN,STRUC,PSIFCN,AMASS
      REAL SR2,AML,AMN,SIGW,PROPZ
      REAL CS2THW,TNTHW,CTTHW,AL(2),BE(2),ESQ,XWI(2),YWI(2)
      REAL ALL(2),BEL(2),EL1
      REAL EQ1,XMGG,XMZZ,XMGZ,XM,CTH2L
      REAL SIGUT,SIGTU,EHAT,PHAT,EBM,TPP,AMWI,AMQ,PROPW
      REAL A,B,ASPBS,ASMBS,TM1,TM2,TM3,COTB,TANB
      INTEGER JTYP1,JTYP2,IFLQ,IUD(13)
      INTEGER IFLL,IL,IN,IDL,IDN,IL1,IL2,JTYPL1,JTYPL2,IDL1,IDL2
C
C          IDENT codes from /SSTYPE/. (Fortran 77 allows - signs in
C          parameter statements but not data statements.)
      INTEGER MSUPL,MSDNL,MSSTL,MSCHL,MSBT1,MSTP1,
     $MSUPR,MSDNR,MSSTR,MSCHR,MSBT2,MSTP2,MSW1,MSW2,
     $MSNEL,MSEL,MSNML,MSMUL,MSNTL,MSTAU1,MSER,MSMUR,MSTAU2
      PARAMETER (MSUPL=-ISUPL)
      PARAMETER (MSDNL=-ISDNL)
      PARAMETER (MSSTL=-ISSTL)
      PARAMETER (MSCHL=-ISCHL)
      PARAMETER (MSBT1=-ISBT1)
      PARAMETER (MSTP1=-ISTP1)
      PARAMETER (MSUPR=-ISUPR)
      PARAMETER (MSDNR=-ISDNR)
      PARAMETER (MSSTR=-ISSTR)
      PARAMETER (MSCHR=-ISCHR)
      PARAMETER (MSBT2=-ISBT2)
      PARAMETER (MSTP2=-ISTP2)
      PARAMETER (MSW1=-ISW1)
      PARAMETER (MSW2=-ISW2)
      PARAMETER (MSNEL=-ISNEL)
      PARAMETER (MSEL=-ISEL)
      PARAMETER (MSNML=-ISNML)
      PARAMETER (MSMUL=-ISMUL)
      PARAMETER (MSNTL=-ISNTL)
      PARAMETER (MSTAU1=-ISTAU1)
      PARAMETER (MSER=-ISER)
      PARAMETER (MSMUR=-ISMUR)
      PARAMETER (MSTAU2=-ISTAU2)
      DATA IDLSS/ISNEL,MSNEL,ISEL,MSEL,ISNML,MSNML,ISMUL,MSMUL,
     $ISNTL,MSNTL,ISTAU1,MSTAU1,ISER,MSER,ISMUR,MSMUR,
     $ISTAU2,MSTAU2/
      DATA IUD/0,1,-1,2,-2,2,-2,1,-1,2,-2,1,-1/
C
C          JS2JT: Susy jettype -> normal jettype
      DATA JS2JT/1,
     $2,3,4,5,6,7,8,9,10,11,12,13,2,3,4,5,6,7,8,9,10,11,12,13/
C          IW2JS: Wino index -> susy jettype
      DATA IW2JS/26,27,28,29/
C          IW2IM: Wino index -> match code
      DATA IW2IM/2,3,2,3/
C          IZ2JS: Zino index -> susy jettype
      DATA IZ2JS/30,31,32,33/
C          IS2UD: Susy jettype -> u/d code
      DATA IS2UD/0,1,1,2,2,2,2,1,1,2,2,1,1,1,1,2,2,2,2,1,1,2,2,1,1/
      DATA IS2LN/1,1,2,2,1,1,2,2,1,1,2,2,2,2,2,2,2,2/
      DATA IL2JS/34,35,36,37,38,39,40,41,42,43,44,45,46,47,
     $48,49,50,51/
C
C          Functions
      QFCN(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
      PSIFCN(AM12,AM22,TT)=((S+TT-AM12)/(2*S)
     $-AM12*(AM22-TT)/(AM12-TT)**2
     $+(TT*(AM22-AM12)+AM22*(S-AM22+AM12))/(S*(AM12-TT)))
C
C          Constants from Baer and Tata,
C
      G=SQRT(4*PI*ALFAEM/SN2THW)
      GP=G*SQRT(SN2THW/(1.-SN2THW))
C          Quark couplings to Z
      CS2THW=1.-SN2THW
      TNTHW=SQRT(SN2THW/CS2THW)
      CTTHW=1./TNTHW
      AL(1)=(CTTHW/4.-5*TNTHW/12.)
      AL(2)=(TNTHW/12.-CTTHW/4.)
      BE(1)=-(CTTHW+TNTHW)/4.
      BE(2)=-BE(1)
      ALL(1)=(CTTHW+TNTHW)/4.
      ALL(2)=(-CTTHW+3*TNTHW)/4.
      BEL(1)=-(CTTHW+TNTHW)/4.
      BEL(2)=-BEL(1)
      ESQ=4*PI*ALFAEM
      SR2=SQRT(2.)
      COTB=RV2V1
      TANB=1./COTB
C
C         qk qb --> slss slbss
C
C
C     Left-leftbar slepton pair production
C
      DO 200 IL=1,6
        IL1=2*IL-1
        IL2=IL1+1
        AML=AMASS(IDLSS(IL1))
        JTYPL1=IL2JS(IL1)
        JTYPL2=IL2JS(IL2)
        IDL1=IDLSS(IL1)
        IDL2=IDLSS(IL2)
          IF (.NOT.(GOQ(JTYPL1,1).AND.GOQ(JTYPL2,2))) GO TO 210
          CALL TWOKIN(0.,0.,AML,AML)
          IF (X1.GE.1..OR.X2.GE.1.) GO TO 210
          E1=SQRT(P(1)**2+AML**2)
          E2=SQRT(P(2)**2+AML**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          DO 220 IQ1=2,11
            IFLQ=IS2UD(IQ1)
            IFLL=IS2LN(IL1)
            IF (IFLQ.EQ.1) THEN
              EQ1=2./3.
            ELSE
              EQ1=-1./3.
            END IF
            IF (IFLL.EQ.1) THEN
              EL1=0.
            ELSE
              EL1=-1.
            END IF
            IQ2=MATCH(IQ1,4)
            PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
            IF (IQ2.EQ.0.OR.IQ2.GE.12) GO TO 220
              XMGG=EL1**2*EQ1**2/S/S
              CTH2L=1.
              IF (JTYPL1.EQ.44) CTH2L=COS(2*THETAL)
              XMZZ=(AL(IFLQ)**2+BE(IFLQ)**2)*(ALL(IFLL)-BEL(IFLL)*
     $              CTH2L)**2/PROPZ
              XMGZ=2*EL1*EQ1*AL(IFLQ)*(ALL(IFLL)-BEL(IFLL)*CTH2L)*
     $             (S-AMZ**2)/S/PROPZ
              XM=2*ESQ*ESQ*(U*T-AML**4)/3.
              SIG=XM*(XMGG+XMZZ+XMGZ)
              SIG=SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
              SIG=.5*SIG
              CALL SIGFIL(SIG,IQ1,IQ2,JTYPL1,JTYPL2)
220       CONTINUE
210     CONTINUE
200   CONTINUE
C          stau_1 + stau_2 bar
      IF (GOQ(44,1).AND.GOQ(51,2)) THEN
          CALL TWOKIN(0.,0.,AML1SS,AML2SS)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 231
          E1=SQRT(P(1)**2+AML1SS**2)
          E2=SQRT(P(2)**2+AML2SS**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
          DO 230 IQ1=2,11
            IFLQ=IS2UD(IQ1)
            IQ2=MATCH(IQ1,4)
            IF (IQ2.EQ.0.OR.IQ2.GE.12) GO TO 230
            SIG=2*ESQ**2*(AL(IFLQ)**2+BE(IFLQ)**2)*BEL(2)**2*
     $     SIN(2*THETAL)**2*(U*T-AML1SS**2*AML2SS**2)/3./PROPZ
            SIG=.5*SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
            CALL SIGFIL(SIG,IQ1,IQ2,44,51)
230         CONTINUE
231       CONTINUE
      END IF
C
C
C     Right-rightbar slepton pair production
C
      DO 300 IL=1,3
        IL1=11+2*IL
        IL2=IL1+1
        AML=AMASS(IDLSS(IL1))
        JTYPL1=IL2JS(IL1)
        JTYPL2=IL2JS(IL2)
        IDL1=IDLSS(IL1)
        IDL2=IDLSS(IL2)
          IF (.NOT.(GOQ(JTYPL1,1).AND.GOQ(JTYPL2,2))) GO TO 310
          CALL TWOKIN(0.,0.,AML,AML)
          IF (X1.GE.1..OR.X2.GE.1.) GO TO 310
          E1=SQRT(P(1)**2+AML**2)
          E2=SQRT(P(2)**2+AML**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          DO 320 IQ1=2,11
            IFLQ=IS2UD(IQ1)
            IFLL=IS2LN(IL1)
            IF (IFLQ.EQ.1) THEN
              EQ1=2./3.
            ELSE
              EQ1=-1./3.
            END IF
            IF (IFLL.EQ.1) THEN
              EL1=0.
            ELSE
              EL1=-1.
            END IF
            IQ2=MATCH(IQ1,4)
            PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
            IF (IQ2.EQ.0.OR.IQ2.GE.12) GO TO 320
              XMGG=EL1**2*EQ1**2/S/S
              CTH2L=1.
              IF (JTYPL1.EQ.50) CTH2L=COS(2*THETAL)
              XMZZ=(AL(IFLQ)**2+BE(IFLQ)**2)*(ALL(IFLL)+BEL(IFLL)*
     $              CTH2L)**2/PROPZ
              XMGZ=2*EL1*EQ1*AL(IFLQ)*(ALL(IFLL)+BEL(IFLL)*CTH2L)*
     $             (S-AMZ**2)/S/PROPZ
              XM=2*ESQ*ESQ*(U*T-AML**4)/3.
              SIG=XM*(XMGG+XMZZ+XMGZ)
              SIG=SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
              SIG=.5*SIG
              CALL SIGFIL(SIG,IQ1,IQ2,JTYPL1,JTYPL2)
320       CONTINUE
310     CONTINUE
300   CONTINUE
C          stau_2 bar + stau_1
      IF (GOQ(51,1).AND.GOQ(44,2)) THEN
          CALL TWOKIN(0.,0.,AML2SS,AML1SS)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 331
          E1=SQRT(P(1)**2+AML2SS**2)
          E2=SQRT(P(2)**2+AML1SS**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
          DO 330 IQ1=2,11
            IFLQ=IS2UD(IQ1)
            IQ2=MATCH(IQ1,4)
            IF (IQ2.EQ.0.OR.IQ2.GE.12) GO TO 330
            SIG=2*ESQ**2*(AL(IFLQ)**2+BE(IFLQ)**2)*BEL(2)**2*
     $     SIN(2*THETAL)**2*(U*T-AML1SS**2*AML2SS**2)/3./PROPZ
            SIG=.5*SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
            CALL SIGFIL(SIG,IQ1,IQ2,51,44)
330         CONTINUE
331       CONTINUE
      END IF
C
C
C     Leftbar-left slepton pair production
C
      DO 400 IL=1,6
        IL1=2*IL
        IL2=IL1-1
        AML=AMASS(IDLSS(IL1))
        JTYPL1=IL2JS(IL1)
        JTYPL2=IL2JS(IL2)
        IDL1=IDLSS(IL1)
        IDL2=IDLSS(IL2)
          IF (.NOT.(GOQ(JTYPL1,1).AND.GOQ(JTYPL2,2))) GO TO 410
          CALL TWOKIN(0.,0.,AML,AML)
          IF (X1.GE.1..OR.X2.GE.1.) GO TO 410
          E1=SQRT(P(1)**2+AML**2)
          E2=SQRT(P(2)**2+AML**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          DO 420 IQ1=2,11
            IFLQ=IS2UD(IQ1)
            IFLL=IS2LN(IL1)
            IF (IFLQ.EQ.1) THEN
              EQ1=2./3.
            ELSE
              EQ1=-1./3.
            END IF
            IF (IFLL.EQ.1) THEN
              EL1=0.
            ELSE
              EL1=-1.
            END IF
            IQ2=MATCH(IQ1,4)
            PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
            IF (IQ2.EQ.0.OR.IQ2.GE.12) GO TO 420
              XMGG=EL1**2*EQ1**2/S/S
              CTH2L=1.
              IF (JTYPL1.EQ.45) CTH2L=COS(2*THETAL)
              XMZZ=(AL(IFLQ)**2+BE(IFLQ)**2)*(ALL(IFLL)-BEL(IFLL)*
     $              CTH2L)**2/PROPZ
              XMGZ=2*EL1*EQ1*AL(IFLQ)*(ALL(IFLL)-BEL(IFLL)*CTH2L)*
     $             (S-AMZ**2)/S/PROPZ
              XM=2*ESQ*ESQ*(U*T-AML**4)/3.
              SIG=XM*(XMGG+XMZZ+XMGZ)
              SIG=SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
              SIG=.5*SIG
              CALL SIGFIL(SIG,IQ1,IQ2,JTYPL1,JTYPL2)
420       CONTINUE
410     CONTINUE
400   CONTINUE
C          stau_1 bar + stau_2
      IF (GOQ(45,1).AND.GOQ(50,2)) THEN
          CALL TWOKIN(0.,0.,AML1SS,AML2SS)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 431
          E1=SQRT(P(1)**2+AML1SS**2)
          E2=SQRT(P(2)**2+AML2SS**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
          DO 430 IQ1=2,11
            IFLQ=IS2UD(IQ1)
            IQ2=MATCH(IQ1,4)
            IF (IQ2.EQ.0.OR.IQ2.GE.12) GO TO 430
            SIG=2*ESQ**2*(AL(IFLQ)**2+BE(IFLQ)**2)*BEL(2)**2*
     $     SIN(2*THETAL)**2*(U*T-AML1SS**2*AML2SS**2)/3./PROPZ
            SIG=.5*SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
            CALL SIGFIL(SIG,IQ1,IQ2,45,50)
430         CONTINUE
431       CONTINUE
      END IF
C
C
C     Rightbar-right slepton pair production
C
      DO 500 IL=1,3
        IL1=12+2*IL
        IL2=IL1-1
        AML=AMASS(IDLSS(IL1))
        JTYPL1=IL2JS(IL1)
        JTYPL2=IL2JS(IL2)
        IDL1=IDLSS(IL1)
        IDL2=IDLSS(IL2)
          IF (.NOT.(GOQ(JTYPL1,1).AND.GOQ(JTYPL2,2))) GO TO 510
          CALL TWOKIN(0.,0.,AML,AML)
          IF (X1.GE.1..OR.X2.GE.1.) GO TO 510
          E1=SQRT(P(1)**2+AML**2)
          E2=SQRT(P(2)**2+AML**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          DO 520 IQ1=2,11
            IFLQ=IS2UD(IQ1)
            IFLL=IS2LN(IL1)
            IF (IFLQ.EQ.1) THEN
              EQ1=2./3.
            ELSE
              EQ1=-1./3.
            END IF
            IF (IFLL.EQ.1) THEN
              EL1=0.
            ELSE
              EL1=-1.
            END IF
            IQ2=MATCH(IQ1,4)
            PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
            IF (IQ2.EQ.0.OR.IQ2.GE.12) GO TO 520
              XMGG=EL1**2*EQ1**2/S/S
              CTH2L=1.
              IF (JTYPL1.EQ.51) CTH2L=COS(2*THETAL)
              XMZZ=(AL(IFLQ)**2+BE(IFLQ)**2)*(ALL(IFLL)+BEL(IFLL)*
     $              CTH2L)**2/PROPZ
              XMGZ=2*EL1*EQ1*AL(IFLQ)*(ALL(IFLL)+BEL(IFLL)*CTH2L)*
     $             (S-AMZ**2)/S/PROPZ
              XM=2*ESQ*ESQ*(U*T-AML**4)/3.
              SIG=XM*(XMGG+XMZZ+XMGZ)
              SIG=SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
              SIG=.5*SIG
              CALL SIGFIL(SIG,IQ1,IQ2,JTYPL1,JTYPL2)
520       CONTINUE
510     CONTINUE
500   CONTINUE
C          stau_2 + stau_1 bar
      IF (GOQ(50,1).AND.GOQ(45,2)) THEN
          CALL TWOKIN(0.,0.,AML2SS,AML1SS)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 531
          E1=SQRT(P(1)**2+AML2SS**2)
          E2=SQRT(P(2)**2+AML1SS**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
          DO 530 IQ1=2,11
            IFLQ=IS2UD(IQ1)
            IQ2=MATCH(IQ1,4)
            IF (IQ2.EQ.0.OR.IQ2.GE.12) GO TO 530
            SIG=2*ESQ**2*(AL(IFLQ)**2+BE(IFLQ)**2)*BEL(2)**2*
     $     SIN(2*THETAL)**2*(U*T-AML1SS**2*AML2SS**2)/3./PROPZ
            SIG=.5*SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
            CALL SIGFIL(SIG,IQ1,IQ2,50,45)
530         CONTINUE
531       CONTINUE
      END IF
C
C          slepton+sneutrino-bar via W-*
C
      DO 600 II=1,3
        IL=4*II-1
        IN=IL-1
        IDL=IDLSS(IL)
        IDN=IDLSS(IN)
        AML=AMASS(IDL)
        AMN=AMASS(IDN)
        JTYP1=IL2JS(IL)
        JTYP2=IL2JS(IN)
          IF(.NOT.(GOQ(JTYP1,1).AND.GOQ(JTYP2,2))) GO TO 610
          CALL TWOKIN(0.,0.,AML,AMN)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 610
          E1=SQRT(P(1)**2+AML**2)
          E2=SQRT(P(2)**2+AMN**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGW=G**4*(U*T-AML**2*AMN**2)/12./PROPW
          IF (JTYP1.EQ.44) SIGW=SIGW*COS(THETAL)**2
            SIG=.5*SIGW*FAC*QFCN(3,1)*QFCN(4,2)
            CALL SIGFIL(SIG,3,4,JTYP1,JTYP2)
            SIG=.5*SIGW*FAC*QFCN(4,1)*QFCN(3,2)
            CALL SIGFIL(SIG,4,3,JTYP1,JTYP2)
            SIG=.5*SIGW*FAC*QFCN(9,1)*QFCN(6,2)
            CALL SIGFIL(SIG,9,6,JTYP1,JTYP2)
            SIG=.5*SIGW*FAC*QFCN(6,1)*QFCN(9,2)
            CALL SIGFIL(SIG,6,9,JTYP1,JTYP2)
610     CONTINUE
600   CONTINUE
C          stau_2 +nu_tau bar
      IF (GOQ(50,1).AND.GOQ(43,2)) THEN
          CALL TWOKIN(0.,0.,AML2SS,AMN3SS)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 620
          E1=SQRT(P(1)**2+AML2SS**2)
          E2=SQRT(P(2)**2+AMN3SS**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGW=G**4*(U*T-AML2SS**2*AMN3SS**2)/12./PROPW
          SIGW=SIGW*SIN(THETAL)**2
          SIG=.5*SIGW*FAC*QFCN(3,1)*QFCN(4,2)
          CALL SIGFIL(SIG,3,4,50,43)
          SIG=.5*SIGW*FAC*QFCN(4,1)*QFCN(3,2)
          CALL SIGFIL(SIG,4,3,50,43)
          SIG=.5*SIGW*FAC*QFCN(9,1)*QFCN(6,2)
          CALL SIGFIL(SIG,9,6,50,43)
          SIG=.5*SIGW*FAC*QFCN(6,1)*QFCN(9,2)
          CALL SIGFIL(SIG,6,9,50,43)
620       CONTINUE
      END IF
C
C          sneutrino-bar+slepton via W-*
C
      DO 700 II=1,3
        IN=4*II-2
        IL=IN+1
        IDL=IDLSS(IL)
        IDN=IDLSS(IN)
        AML=AMASS(IDL)
        AMN=AMASS(IDN)
        JTYP1=IL2JS(IN)
        JTYP2=IL2JS(IL)
          IF(.NOT.(GOQ(JTYP1,1).AND.GOQ(JTYP2,2))) GO TO 710
          CALL TWOKIN(0.,0.,AMN,AML)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 710
          E1=SQRT(P(1)**2+AMN**2)
          E2=SQRT(P(2)**2+AML**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGW=G**4*(U*T-AML**2*AMN**2)/12./PROPW
          IF (JTYP2.EQ.44) SIGW=SIGW*COS(THETAL)**2
            SIG=.5*SIGW*FAC*QFCN(3,1)*QFCN(4,2)
            CALL SIGFIL(SIG,3,4,JTYP1,JTYP2)
            SIG=.5*SIGW*FAC*QFCN(4,1)*QFCN(3,2)
            CALL SIGFIL(SIG,4,3,JTYP1,JTYP2)
            SIG=.5*SIGW*FAC*QFCN(9,1)*QFCN(6,2)
            CALL SIGFIL(SIG,9,6,JTYP1,JTYP2)
            SIG=.5*SIGW*FAC*QFCN(6,1)*QFCN(9,2)
            CALL SIGFIL(SIG,6,9,JTYP1,JTYP2)
710     CONTINUE
700   CONTINUE
C          nu_tau bar + STAU_2
      IF (GOQ(43,1).AND.GOQ(50,2)) THEN
          CALL TWOKIN(0.,0.,AMN3SS,AML2SS)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 720
          E1=SQRT(P(1)**2+AMN3SS**2)
          E2=SQRT(P(2)**2+AML2SS**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGW=G**4*(U*T-AML2SS**2*AMN3SS**2)/12./PROPW
          SIGW=SIGW*SIN(THETAL)**2
          SIG=.5*SIGW*FAC*QFCN(3,1)*QFCN(4,2)
          CALL SIGFIL(SIG,3,4,43,50)
          SIG=.5*SIGW*FAC*QFCN(4,1)*QFCN(3,2)
          CALL SIGFIL(SIG,4,3,43,50)
          SIG=.5*SIGW*FAC*QFCN(9,1)*QFCN(6,2)
          CALL SIGFIL(SIG,9,6,43,50)
          SIG=.5*SIGW*FAC*QFCN(6,1)*QFCN(9,2)
          CALL SIGFIL(SIG,6,9,43,50)
720       CONTINUE
      END IF
C
C          slepton-bar+sneutrino via W+*
C
      DO 800 II=1,3
        IL=4*II
        IN=IL-3
        IDL=IDLSS(IL)
        IDN=IDLSS(IN)
        AML=AMASS(IDL)
        AMN=AMASS(IDN)
        JTYP1=IL2JS(IL)
        JTYP2=IL2JS(IN)
          IF(.NOT.(GOQ(JTYP1,1).AND.GOQ(JTYP2,2))) GO TO 810
          CALL TWOKIN(0.,0.,AML,AMN)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 810
          E1=SQRT(P(1)**2+AML**2)
          E2=SQRT(P(2)**2+AMN**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGW=G**4*(U*T-AML**2*AMN**2)/12./PROPW
          IF (JTYP1.EQ.45) SIGW=SIGW*COS(THETAL)**2
            SIG=.5*SIGW*FAC*QFCN(2,1)*QFCN(5,2)
            CALL SIGFIL(SIG,2,5,JTYP1,JTYP2)
            SIG=.5*SIGW*FAC*QFCN(5,1)*QFCN(2,2)
            CALL SIGFIL(SIG,5,2,JTYP1,JTYP2)
            SIG=.5*SIGW*FAC*QFCN(8,1)*QFCN(7,2)
            CALL SIGFIL(SIG,8,7,JTYP1,JTYP2)
            SIG=.5*SIGW*FAC*QFCN(7,1)*QFCN(8,2)
            CALL SIGFIL(SIG,7,8,JTYP1,JTYP2)
810     CONTINUE
800   CONTINUE
C          stau_2 bar+nu_tau
      IF (GOQ(51,1).AND.GOQ(42,2)) THEN
          CALL TWOKIN(0.,0.,AML2SS,AMN3SS)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 820
          E1=SQRT(P(1)**2+AML2SS**2)
          E2=SQRT(P(2)**2+AMN3SS**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGW=G**4*(U*T-AML2SS**2*AMN3SS**2)/12./PROPW
          SIGW=SIGW*SIN(THETAL)**2
          SIG=.5*SIGW*FAC*QFCN(2,1)*QFCN(5,2)
          CALL SIGFIL(SIG,2,5,51,42)
          SIG=.5*SIGW*FAC*QFCN(5,1)*QFCN(2,2)
          CALL SIGFIL(SIG,5,2,51,42)
          SIG=.5*SIGW*FAC*QFCN(8,1)*QFCN(7,2)
          CALL SIGFIL(SIG,8,7,51,42)
          SIG=.5*SIGW*FAC*QFCN(7,1)*QFCN(8,2)
          CALL SIGFIL(SIG,7,8,51,42)
820       CONTINUE
      END IF
C
C          sneutrino+slepton-bar via W+*
C
      DO 900 II=1,3
        IN=4*II-3
        IL=IN+3
        IDL=IDLSS(IL)
        IDN=IDLSS(IN)
        AML=AMASS(IDL)
        AMN=AMASS(IDN)
        JTYP1=IL2JS(IN)
        JTYP2=IL2JS(IL)
          IF(.NOT.(GOQ(JTYP1,1).AND.GOQ(JTYP2,2))) GO TO 910
          CALL TWOKIN(0.,0.,AMN,AML)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 910
          E1=SQRT(P(1)**2+AMN**2)
          E2=SQRT(P(2)**2+AML**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGW=G**4*(U*T-AML**2*AMN**2)/12./PROPW
          IF (JTYP2.EQ.45) SIGW=SIGW*COS(THETAL)**2
            SIG=.5*SIGW*FAC*QFCN(2,1)*QFCN(5,2)
            CALL SIGFIL(SIG,2,5,JTYP1,JTYP2)
            SIG=.5*SIGW*FAC*QFCN(5,1)*QFCN(2,2)
            CALL SIGFIL(SIG,5,2,JTYP1,JTYP2)
            SIG=.5*SIGW*FAC*QFCN(8,1)*QFCN(7,2)
            CALL SIGFIL(SIG,8,7,JTYP1,JTYP2)
            SIG=.5*SIGW*FAC*QFCN(7,1)*QFCN(8,2)
            CALL SIGFIL(SIG,7,8,JTYP1,JTYP2)
910     CONTINUE
900   CONTINUE
C          nu_tau + stau_2 bar
      IF (GOQ(42,1).AND.GOQ(51,2)) THEN
          CALL TWOKIN(0.,0.,AMN3SS,AML2SS)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 920
          E1=SQRT(P(1)**2+AMN3SS**2)
          E2=SQRT(P(2)**2+AML2SS**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGW=G**4*(U*T-AML2SS**2*AMN3SS**2)/12./PROPW
          SIGW=SIGW*SIN(THETAL)**2
          SIG=.5*SIGW*FAC*QFCN(2,1)*QFCN(5,2)
          CALL SIGFIL(SIG,2,5,42,51)
          SIG=.5*SIGW*FAC*QFCN(5,1)*QFCN(2,2)
          CALL SIGFIL(SIG,5,2,42,51)
          SIG=.5*SIGW*FAC*QFCN(8,1)*QFCN(7,2)
          CALL SIGFIL(SIG,8,7,42,51)
          SIG=.5*SIGW*FAC*QFCN(7,1)*QFCN(8,2)
          CALL SIGFIL(SIG,7,8,42,51)
920       CONTINUE
      END IF
C
      RETURN
      END
