CDECK  ID>, STBWZ1.
      SUBROUTINE STBWZ1(GAMMA)
C
CPurpose: To carry out all the necessary calculations and calls to
C         find the decay width of
C                  stop --> bottom + W + lightest neutralino
C
C         Returns GAMMA, the width of the decay, in GeV.        
C
      IMPLICIT NONE
C
CIsajet common blocks
C
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
C
C
CMy common blocks
C
      COMMON/ST3WID/SBWIDTH,CHWIDTH,TWIDTH
      DOUBLE PRECISION SBWIDTH(2),CHWIDTH(2),TWIDTH
      SAVE/ST3WID/
C
      COMMON/ST3MAS/MSB,MST,MCH,MNEU,MW,MB,MT
      DOUBLE PRECISION MSB(2),MST(1),MCH(2),MNEU(1),MW,MB,MT
      SAVE/ST3MAS/
C
      COMMON/ST3COUP/CA1,CB1,CA2,CB2,CC2,CD2,CE2,CF2,CA3,CB3,CC3,
     $               CA12,CB12,CC12,CD12,CA13,CB13,CC13,CD13,
     $               CA23,CB23,CC23,CD23,CE23,CF23,CG23,CH23,
     $               A11,A32
      DOUBLE PRECISION A11(2),A32
      DOUBLE COMPLEX CA1(2,2),CB1(2,2),CA2(2,2),CB2(2,2),CC2(2,2),
     $               CD2(2,2),CE2(2,2),CF2(2,2),CA3,CB3,CC3,
     $               CA12(2,2),CB12(2,2),CC12(2,2),CD12(2,2),
     $               CA13(2),CB13(2),CC13(2),CD13(2),
     $               CA23(2),CB23(2),CC23(2),CD23(2),CE23(2),
     $               CF23(2),CG23(2),CH23(2)
      SAVE/ST3COUP/
C
      DOUBLE PRECISION GAMMA,INTMSQ,PI,GPME,COST,SINT,SNW1,SNW2,
     $                 AWD(2),BW(2),BWP(2),XM,YM,THX,THY,
     $                 FT,FB,TANB,MZIZ,SNZI,THIZ,YIP,XIP,YIM,XIM,
     $                 SN1,X1(2),Y1(2),LT1I(2),KT1I(2),COSB,SINB,
     $                 THT,THB,GME,A21(2),B21(2)
      DOUBLE COMPLEX OL1I(2),OR1I(2),ZAUIZ,ZBUIZ,ZADIZ,ZBDIZ,
     $               ZPP,ZPM,ZA,ZB,A12(2),A22(2),A31,B12(2),
     $               B22(2),B31
      INTEGER I,J
C
      GAMMA=0.D0
      PI=4.D0*DATAN(1.D0)
C
CEnter the masses into the common blocks
C
      MSB(1)=DBLE(AMB1SS)
      MSB(2)=DBLE(AMB2SS)
      MST(1)=DBLE(AMT1SS)
      MCH(1)=DABS(DBLE(AMW1SS))
      MCH(2)=DABS(DBLE(AMW2SS))
      MNEU(1)=DABS(DBLE(AMZ1SS))
      MW=DBLE(AMW)
      MB=DBLE(AMBT)
      MT=DBLE(AMTP)
C
CCheck the kinematics. If the stop is too light for the three body
Cdecay or the two body decays are available, then do not calculate
Cthe width and return a value of 0.D0.
C
      IF(MST(1).LT.(MB+MW+MNEU(1)))THEN
        WRITE(*,*)
        WRITE(*,*)'STOP TOO LIGHT FOR THREE BODY DECAY'
        RETURN
      ELSE IF(MST(1).GT.(MNEU(1)+MT)
     .   .OR.MST(1).GT.(DABS(DBLE(AMZ2SS))+MT)
     .   .OR.MST(1).GT.(DABS(DBLE(AMZ3SS))+MT)
     .   .OR.MST(1).GT.(DABS(DBLE(AMZ4SS))+MT)
     .   .OR.MST(1).GT.(DABS(DBLE(AMW1SS))+MB)
     .   .OR.MST(1).GT.(DABS(DBLE(AMW2SS))+MB)
     .   .OR.MST(1).GT.(MSB(1)+MW).OR.MST(1).GT.(MSB(2)+MW))THEN
C       WRITE(*,*)
C       WRITE(*,*)'TWO BODY STOP DECAY OPEN'
        RETURN
      END IF
C
CDeal with the widths. The ident codes (25,45,39,49,6) are in the
Cisajet manual.
C
      SBWIDTH(1)=0.D0
      SBWIDTH(2)=0.D0
      CHWIDTH(1)=0.D0
      CHWIDTH(2)=0.D0
      TWIDTH=0.D0
C
      DO I=1,NSSMOD
        IF(ISSMOD(I).EQ.25)SBWIDTH(1)=SBWIDTH(1)+DBLE(GSSMOD(I))
        IF(ISSMOD(I).EQ.45)SBWIDTH(2)=SBWIDTH(2)+DBLE(GSSMOD(I))
        IF(ISSMOD(I).EQ.39)CHWIDTH(1)=CHWIDTH(1)+DBLE(GSSMOD(I))
        IF(ISSMOD(I).EQ.49)CHWIDTH(2)=CHWIDTH(2)+DBLE(GSSMOD(I))
        IF(ISSMOD(I).EQ.6)TWIDTH=TWIDTH+DBLE(GSSMOD(I))
      END DO
C
CNow take care of the couplings
C
      THB=DBLE(THETAB) !These angles are still BT convention.
      THT=DBLE(THETAT) !The following defs. rely on this fact.
      COST=DCOS(THT)
      SINT=DSIN(THT)
      COSB=DCOS(THB)
      SINB=DSIN(THB)
C
CNOTE: This mostly follows isajet, but many variables are now dble
CThe values of g and g' are not precise since they should be found
Cby running up to m_s using DR-bar SM RGEs.
C
      GME=DSQRT(4.D0*PI*DBLE(ALFAEM/SN2THW))
      GPME=GME*DSQRT(DBLE(SN2THW)/(1.D0-DBLE(SN2THW)))
      SNW1=DSIGN(1.D0,DBLE(AMW1SS))
      SNW2=DSIGN(1.D0,DBLE(AMW2SS))
      XM=1.D0/DTAN(DBLE(GAMMAL))
      YM=1.D0/DTAN(DBLE(GAMMAR))
      THX=DSIGN(1.D0,XM)
      THY=DSIGN(1.D0,YM)
      TANB=1.D0/DBLE(RV2V1)
      FT=GME*DBLE(MTQ)/DSQRT(2.D0)/DBLE(AMW)/DSIN(DATAN(TANB))
      FB=GME*DBLE(MBQ)/DSQRT(2.D0)/DBLE(AMW)/DCOS(DATAN(TANB))
      AWD(1)=-GME*SNW1*DSIN(DBLE(GAMMAR))
      AWD(2)=-GME*SNW2*THY*DCOS(DBLE(GAMMAR))
      BW(1)=-FT*SNW1*DCOS(DBLE(GAMMAR))
      BW(2)=FT*SNW2*THY*DSIN(DBLE(GAMMAR))
      BWP(1)=-FB*DCOS(DBLE(GAMMAL))
      BWP(2)=FB*THX*DSIN(DBLE(GAMMAL))
C
      LT1I(1)=AWD(1)*COST-BW(1)*SINT
      KT1I(1)=BWP(1)*COST
      LT1I(2)=AWD(2)*COST-BW(2)*SINT
      KT1I(2)=BWP(2)*COST
C
      MZIZ=DABS(DBLE(AMZISS(1)))
      SNZI=DSIGN(1.D0,DBLE(AMZISS(1)))
      IF(SNZI.EQ.1.D0)THEN
        THIZ=0.D0
      ELSE
        THIZ=1.D0
      END IF
      ZAUIZ=(0.D0,1.D0)**(THIZ-1.D0)*SNZI
     $  *(-GME/DSQRT(2.D0)*DBLE(ZMIXSS(3,1))-GPME/3.D0/DSQRT(2.D0)
     $  *DBLE(ZMIXSS(4,1)))
      ZBUIZ=(0.D0,1.D0)**(THIZ-1.D0)*4.D0*GPME*DBLE(ZMIXSS(4,1))
     $  /3.D0/DSQRT(2.D0)
      ZPP=(0.D0,1.D0)**THIZ
      ZPM=(-(0.D0,1.D0))**THIZ
      ZA=(((0.D0,1.D0)*ZAUIZ-ZPP*FT*DBLE(ZMIXSS(1,1)))*COST
     $  -((0.D0,1.D0)*ZBUIZ-ZPM*FT*DBLE(ZMIXSS(1,1)))*SINT
     $   )/2.D0
      ZB=((-(0.D0,1.D0)*ZAUIZ-ZPP*FT*DBLE(ZMIXSS(1,1)))*COST
     $  -((0.D0,1.D0)*ZBUIZ+ZPM*FT*DBLE(ZMIXSS(1,1)))*SINT
     &   )/2.D0
C
      A31=DCONJG(ZA)
      B31=DCONJG(-ZB)
C
      XIM=.5D0*(SNW1*DSIGN(1.D0,DBLE(AMZISS(1)))*(DCOS(DBLE(GAMMAR))
     $  *DBLE(ZMIXSS(1,1))/DSQRT(2.D0)+DSIN(DBLE(GAMMAR))
     $  *DBLE(ZMIXSS(3,1)))-DCOS(DBLE(GAMMAL))
     $  *DBLE(ZMIXSS(2,1))/DSQRT(2.D0)+DSIN(DBLE(GAMMAL))
     $  *DBLE(ZMIXSS(3,1)))
      YIM=.5D0*(-SNW1*DSIGN(1.D0,DBLE(AMZISS(1)))*(DCOS(DBLE(GAMMAR))
     $  *DBLE(ZMIXSS(1,1))/DSQRT(2.D0)+DSIN(DBLE(GAMMAR))
     $  *DBLE(ZMIXSS(3,1)))-DCOS(DBLE(GAMMAL))
     $  *DBLE(ZMIXSS(2,1))/DSQRT(2.D0)+DSIN(DBLE(GAMMAL))
     $  *DBLE(ZMIXSS(3,1)))
      XIP=.5D0*(SNW2*DSIGN(1.D0,DBLE(AMZISS(1)))*THY
     $  *(-DSIN(DBLE(GAMMAR))
     $  *DBLE(ZMIXSS(1,1))/DSQRT(2.D0)+DCOS(DBLE(GAMMAR))
     $  *DBLE(ZMIXSS(3,1)))+THX*(DSIN(DBLE(GAMMAL))
     $  *DBLE(ZMIXSS(2,1))/DSQRT(2.D0)+DCOS(DBLE(GAMMAL))
     $  *DBLE(ZMIXSS(3,1))))
      YIP=.5D0*(-SNW2*DSIGN(1.D0,DBLE(AMZISS(1)))*THY
     $  *(-DSIN(DBLE(GAMMAR))
     $  *DBLE(ZMIXSS(1,1))/DSQRT(2.D0)+DCOS(DBLE(GAMMAR))
     $  *DBLE(ZMIXSS(3,1)))+THX*(DSIN(DBLE(GAMMAL))
     $  *DBLE(ZMIXSS(2,1))/DSQRT(2.D0)+DCOS(DBLE(GAMMAL))
     $  *DBLE(ZMIXSS(3,1))))
      SN1=DSIGN(1.D0,DBLE(AMZISS(1)))
      X1(1)=XIM
      Y1(1)=YIM
      X1(2)=XIP
      Y1(2)=YIP
C
      OL1I(1)=(-(0.D0,1.D0))**((1.D0-SN1)/2.D0)*(X1(1)-Y1(1))
      OR1I(1)=(-(0.D0,1.D0))**((1.D0-SN1)/2.D0)*(X1(1)+Y1(1))
      OL1I(2)=(-(0.D0,1.D0))**((1.D0-SN1)/2.D0)*(X1(2)-Y1(2))
      OR1I(2)=(-(0.D0,1.D0))**((1.D0-SN1)/2.D0)*(X1(2)+Y1(2))
C
      ZADIZ=(0.D0,1.D0)**(THIZ-1.D0)*(-1.D0)*SNZI
     $      *(-GME/DSQRT(2.D0)*DBLE(ZMIXSS(3,1))+GPME/3.D0/DSQRT(2.D0)
     $      *DBLE(ZMIXSS(4,1)))
      ZBDIZ=(0.D0,1.D0)**(THIZ-1.D0)*(-2.D0)*GPME
     $      *DBLE(ZMIXSS(4,1))/3.D0/DSQRT(2.D0)
      ZA=((0.D0,1.D0)*ZADIZ-FB*DBLE(ZMIXSS(2,1))*(0.D0,1.D0)**THIZ)
     $   *COSB/2.D0-
     $       ((0.D0,1.D0)*ZBDIZ-FB*DBLE(ZMIXSS(2,1))
     $   *(-(0.D0,1.D0))**THIZ)*SINB/2.D0
      ZB=(-(0.D0,1.D0)*ZADIZ-FB*DBLE(ZMIXSS(2,1))*(0.D0,1.D0)**THIZ)
     $   *COSB/2.D0-
     $       ((0.D0,1.D0)*ZBDIZ+FB*DBLE(ZMIXSS(2,1))
     $    *(-(0.D0,1.D0))**THIZ)*SINB/2.D0
C
      A12(1)=CONJG(ZA)
      B12(1)=CONJG(-ZB)
C
      ZA=((0.D0,1.D0)*ZADIZ-FB*DBLE(ZMIXSS(2,1))*(0.D0,1.D0)**THIZ)
     $   *SINB/2.D0+
     $       ((0.D0,1.D0)*ZBDIZ-FB*DBLE(ZMIXSS(2,1))
     $   *(-(0.D0,1.D0))**THIZ)*COSB/2.D0
      ZB=(-(0.D0,1.D0)*ZADIZ-FB*DBLE(ZMIXSS(2,1))*(0.D0,1.D0)**THIZ)
     $   *SINB/2.D0+
     $       ((0.D0,1.D0)*ZBDIZ+FB*DBLE(ZMIXSS(2,1))
     $    *(-(0.D0,1.D0))**THIZ)*COSB/2.D0
C
      A12(2)=CONJG(ZA)
      B12(2)=CONJG(-ZB)
C
      A11(1)=-GME*COSB*COST/DSQRT(2.D0)
      A11(2)=-GME*SINB*COST/DSQRT(2.D0)
C
      DO I=1,2
        A21(I)=1.D0/2.D0*(LT1I(I)+KT1I(I))
        B21(I)=1.D0/2.D0*(LT1I(I)-KT1I(I))
        A22(I)=-GME/2.D0*(OL1I(I)+OR1I(I))
        B22(I)=-GME/2.D0*(OL1I(I)-OR1I(I))
      END DO
C
      A32=-GME/(2.D0*DSQRT(2.D0))
C
CThe following are used directly in the calculation of the terms
Cin the matrix element. They are defined here to increase efficiency
Cof the integration.
C
      DO I=1,2
        DO J=1,2
          CA1(I,J)=CONJG(A12(I))*A12(J)+CONJG(B12(I))*B12(J)
          CB1(I,J)=CONJG(A12(I))*A12(J)-CONJG(B12(I))*B12(J)
          CA2(I,J)=A21(I)*CONJG(A22(I))*A21(J)*A22(J)
     .             +B21(I)*CONJG(B22(I))*B21(J)*B22(J)
     .             +A21(I)*CONJG(A22(I))*B21(J)*B22(J)
     .             +A21(I)*A22(J)*B21(J)*CONJG(B22(I))
     .             +A21(I)*A21(J)*CONJG(B22(I))*B22(J)
     .             +B21(I)*CONJG(B22(I))*A21(J)*A22(J)
     .             +B21(I)*B22(J)*A21(J)*CONJG(A22(I))
     .             +B21(I)*B21(J)*CONJG(A22(I))*A22(J)
          CB2(I,J)=A21(I)*CONJG(A22(I))*A21(J)*A22(J)
     .             +B21(I)*CONJG(B22(I))*B21(J)*B22(J)
     .             -A21(I)*CONJG(A22(I))*B21(J)*B22(J)
     .             +A21(I)*A22(J)*B21(J)*CONJG(B22(I))
     .             -A21(I)*A21(J)*CONJG(B22(I))*B22(J)
     .             -B21(I)*CONJG(B22(I))*A21(J)*A22(J)
     .             +B21(I)*B22(J)*A21(J)*CONJG(A22(I))
     .             -B21(I)*B21(J)*CONJG(A22(I))*A22(J)
          CC2(I,J)=(A21(I)*CONJG(A22(I))*A21(J)*A22(J)
     .             -B21(I)*CONJG(B22(I))*B21(J)*B22(J)
     .             -A21(I)*CONJG(A22(I))*B21(J)*B22(J)
     .             +A21(I)*A22(J)*B21(J)*CONJG(B22(I))
     .             -A21(I)*A21(J)*CONJG(B22(I))*B22(J)
     .             +B21(I)*CONJG(B22(I))*A21(J)*A22(J)
     .             -B21(I)*B22(J)*A21(J)*CONJG(A22(I))
     .             +B21(I)*B21(J)*CONJG(A22(I))*A22(J))*MCH(J)+
     .             (A21(J)*CONJG(A22(J))*A21(I)*A22(I)
     .             -B21(J)*CONJG(B22(J))*B21(I)*B22(I)
     .             -A21(J)*CONJG(A22(J))*B21(I)*B22(I)
     .             +A21(J)*A22(I)*B21(I)*CONJG(B22(J))
     .             -A21(J)*A21(I)*CONJG(B22(J))*B22(I)
     .             +B21(J)*CONJG(B22(J))*A21(I)*A22(I)
     .             -B21(J)*B22(I)*A21(I)*CONJG(A22(J))
     .             +B21(J)*B21(I)*CONJG(A22(J))*A22(I))*MCH(I)
          CD2(I,J)=A21(I)*CONJG(A22(I))*A21(J)*A22(J)
     .             +B21(I)*CONJG(B22(I))*B21(J)*B22(J)
     .             -A21(I)*CONJG(A22(I))*B21(J)*B22(J)
     .             -A21(I)*A22(J)*B21(J)*CONJG(B22(I))
     .             +A21(I)*A21(J)*CONJG(B22(I))*B22(J)
     .             -B21(I)*CONJG(B22(I))*A21(J)*A22(J)
     .             -B21(I)*B22(J)*A21(J)*CONJG(A22(I))
     .             +B21(I)*B21(J)*CONJG(A22(I))*A22(J)
          CE2(I,J)=A21(I)*CONJG(A22(I))*A21(J)*A22(J)
     .             +B21(I)*CONJG(B22(I))*B21(J)*B22(J)
     .             +A21(I)*CONJG(A22(I))*B21(J)*B22(J)
     .             -A21(I)*A22(J)*B21(J)*CONJG(B22(I))
     .             -A21(I)*A21(J)*CONJG(B22(I))*B22(J)
     .             +B21(I)*CONJG(B22(I))*A21(J)*A22(J)
     .             -B21(I)*B22(J)*A21(J)*CONJG(A22(I))
     .             -B21(I)*B21(J)*CONJG(A22(I))*A22(J)
          CF2(I,J)=(A21(I)*CONJG(A22(I))*A21(J)*A22(J)
     .             -B21(I)*CONJG(B22(I))*B21(J)*B22(J)
     .             -A21(I)*CONJG(A22(I))*B21(J)*B22(J)
     .             -A21(I)*A22(J)*B21(J)*CONJG(B22(I))
     .             +A21(I)*A21(J)*CONJG(B22(I))*B22(J)
     .             +B21(I)*CONJG(B22(I))*A21(J)*A22(J)
     .             +B21(I)*B22(J)*A21(J)*CONJG(A22(I))
     .             -B21(I)*B21(J)*CONJG(A22(I))*A22(J))*MCH(J)+
     .             (A21(J)*CONJG(A22(J))*A21(I)*A22(I)
     .             -B21(J)*CONJG(B22(J))*B21(I)*B22(I)
     .             -A21(J)*CONJG(A22(J))*B21(I)*B22(I)
     .             -A21(J)*A22(I)*B21(I)*CONJG(B22(J))
     .             +A21(J)*A21(I)*CONJG(B22(J))*B22(I)
     .             +B21(J)*CONJG(B22(J))*A21(I)*A22(I)
     .             +B21(J)*B22(I)*A21(I)*CONJG(A22(J))
     .             -B21(J)*B21(I)*CONJG(A22(J))*A22(I))*MCH(I)
          CA12(I,J)=A21(I)*CONJG(A22(I))*A12(J)
     .              +A21(I)*CONJG(B22(I))*B12(J)
     .              +CONJG(A22(I))*B21(I)*B12(J)
     .              +B21(I)*CONJG(B22(I))*A12(J)
          CB12(I,J)=A21(I)*CONJG(A22(I))*A12(J)
     .              -A21(I)*CONJG(B22(I))*B12(J)
     .              +CONJG(A22(I))*B21(I)*B12(J)
     .              -B21(I)*CONJG(B22(I))*A12(J)
          CC12(I,J)=A21(I)*CONJG(A22(I))*A12(J)
     .              -A21(I)*CONJG(B22(I))*B12(J)
     .              -CONJG(A22(I))*B21(I)*B12(J)
     .              +B21(I)*CONJG(B22(I))*A12(J)
          CD12(I,J)=A21(I)*CONJG(A22(I))*A12(J)
     .              +A21(I)*CONJG(B22(I))*B12(J)
     .              -CONJG(A22(I))*B21(I)*B12(J)
     .              -B21(I)*CONJG(B22(I))*A12(J)
        END DO
        CA13(I)=A31*CONJG(A12(I))+B31*CONJG(B12(I))+
     .          CONJG(A12(I))*B31+A31*CONJG(B12(I))
        CB13(I)=A31*CONJG(A12(I))+B31*CONJG(B12(I))-
     .          CONJG(A12(I))*B31-A31*CONJG(B12(I))
        CC13(I)=A31*CONJG(A12(I))-B31*CONJG(B12(I))-
     .          CONJG(A12(I))*B31+A31*CONJG(B12(I))
        CD13(I)=A31*CONJG(A12(I))-B31*CONJG(B12(I))+
     .          CONJG(A12(I))*B31-A31*CONJG(B12(I))
        CA23(I)=A21(I)*CONJG(A22(I))*A31+B21(I)*CONJG(A22(I))*A31+
     .          A21(I)*CONJG(B22(I))*B31+B21(I)*CONJG(B22(I))*B31+
     .          A21(I)*CONJG(B22(I))*A31+B21(I)*CONJG(B22(I))*A31+
     .          A21(I)*CONJG(A22(I))*B31+B21(I)*CONJG(A22(I))*B31
        CB23(I)=A21(I)*CONJG(A22(I))*A31+B21(I)*CONJG(A22(I))*A31+
     .          A21(I)*CONJG(B22(I))*B31+B21(I)*CONJG(B22(I))*B31-
     .          A21(I)*CONJG(B22(I))*A31-B21(I)*CONJG(B22(I))*A31-
     .          A21(I)*CONJG(A22(I))*B31-B21(I)*CONJG(A22(I))*B31
        CC23(I)=A21(I)*CONJG(A22(I))*A31+B21(I)*CONJG(A22(I))*A31-
     .          A21(I)*CONJG(B22(I))*B31-B21(I)*CONJG(B22(I))*B31+
     .          A21(I)*CONJG(B22(I))*A31+B21(I)*CONJG(B22(I))*A31-
     .          A21(I)*CONJG(A22(I))*B31-B21(I)*CONJG(A22(I))*B31
        CD23(I)=A21(I)*CONJG(A22(I))*A31+B21(I)*CONJG(A22(I))*A31-
     .          A21(I)*CONJG(B22(I))*B31-B21(I)*CONJG(B22(I))*B31-
     .          A21(I)*CONJG(B22(I))*A31-B21(I)*CONJG(B22(I))*A31+
     .          A21(I)*CONJG(A22(I))*B31+B21(I)*CONJG(A22(I))*B31
        CE23(I)=A21(I)*CONJG(A22(I))*A31-B21(I)*CONJG(A22(I))*A31+
     .          A21(I)*CONJG(B22(I))*B31-B21(I)*CONJG(B22(I))*B31-
     .          A21(I)*CONJG(B22(I))*A31+B21(I)*CONJG(B22(I))*A31-
     .          A21(I)*CONJG(A22(I))*B31+B21(I)*CONJG(A22(I))*B31
        CF23(I)=A21(I)*CONJG(A22(I))*A31-B21(I)*CONJG(A22(I))*A31+
     .          A21(I)*CONJG(B22(I))*B31-B21(I)*CONJG(B22(I))*B31+
     .          A21(I)*CONJG(B22(I))*A31-B21(I)*CONJG(B22(I))*A31+
     .          A21(I)*CONJG(A22(I))*B31-B21(I)*CONJG(A22(I))*B31
        CG23(I)=A21(I)*CONJG(A22(I))*A31-B21(I)*CONJG(A22(I))*A31-
     .          A21(I)*CONJG(B22(I))*B31+B21(I)*CONJG(B22(I))*B31-
     .          A21(I)*CONJG(B22(I))*A31+B21(I)*CONJG(B22(I))*A31+
     .          A21(I)*CONJG(A22(I))*B31-B21(I)*CONJG(A22(I))*B31
        CH23(I)=A21(I)*CONJG(A22(I))*A31-B21(I)*CONJG(A22(I))*A31-
     .          A21(I)*CONJG(B22(I))*B31+B21(I)*CONJG(B22(I))*B31+
     .          A21(I)*CONJG(B22(I))*A31-B21(I)*CONJG(B22(I))*A31-
     .          A21(I)*CONJG(A22(I))*B31+B21(I)*CONJG(A22(I))*B31
      END DO
      CA3=A32**2*(ABS(A31+B31))**2
      CB3=A32**2*(ABS(A31-B31))**2
      CC3=A32**2*(ABS(A31)**2-ABS(B31)**2)
C
CWith the couplings complete, start the integration
C
      CALL ST3INT(INTMSQ)
C
CFinally calculate the width in GeV
C
      GAMMA=INTMSQ/(64.D0*PI**3*MST(1))
C
      RETURN
      END
