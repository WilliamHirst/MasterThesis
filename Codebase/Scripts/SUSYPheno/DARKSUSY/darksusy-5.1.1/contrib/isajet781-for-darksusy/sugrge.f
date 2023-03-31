CDECK  ID>, SUGRGE.
      SUBROUTINE SUGRGE(M0,MHF,A0,TANB,SGNMU,MT,G,G0,IG,W2
     $,NSTEP,IMODEL,BADMU)
C
C          Make one complete iteration of the renormalization group
C          equations from MZ to MGUT and back, setting the boundary
C          conditions on each end.
C
      IMPLICIT NONE
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
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
      COMMON/SSINF/XLAM
      DOUBLE PRECISION XLAM
C
      COMMON /BSG/GISA(31),MSQISA(3),MSLISA(3),MSUISA(3),MSDISA(3),
     &            MSEISA(3),MRNISA(3),YNFRZ(3,3),MNFRZ(3,3),TNFRZ(3,3),
     &            RTISA,RBISA,RLISA
c     MSxDEC(i) - decoupling scale of i-th generation of type x sfermion
c     MRNDEC(i) - decoupling scale of i-th RH neutrino
      REAL*8 GISA,MSQISA,MSLISA,MSUISA,MSDISA,MSEISA,MRNISA,
     &       YNFRZ,MNFRZ,TNFRZ
      REAL RTISA,RBISA,RLISA
      SAVE /BSG/
      EXTERNAL SURG26
      DOUBLE PRECISION DDILOG,XLM
      COMPLEX*16 SSB0,SSB1
      REAL*8 G(31),W2(93),T,DT,DY,DPI,DAS
      REAL*8 BTHAT,BBHAT,BLHAT
      REAL M0,MHF,A0,TANB,SGNMU,MT,G0(31)
      INTEGER IG(31),NSTEP,IMODEL
      REAL PI,TZ,A1I,A2I,A3I,GGUT,AGUTI,SIG1,SIG2,
     $MH1S,MH2S,MUS,MZ,TGUT,AGUT,Q,ASMT,MTMT,
     $QNEW,XLAMGM,XMESGM,XN5GM,XC,G3GUT,THRF,THRG,MBMZ,
     $M2,AM2,MSN,MG,MT1,MT2,MB1,MB2,MW1,MW2,AMU,
     $RSIGT,RSIGL,RSIGB,DAEM,ALEMDR,ALEM,MTAMZ,
     $TANBQ,SIN2BQ,SINB,COSB,XWMSB
      REAL SSRSGT,SSRSGB,SSRSGL,SUALFS,MZQ,SIGA
      REAL CF,CA,ZETA2,ZETA3,ST2LP,COTB
      INTEGER I,II
      LOGICAL BADMU
C
      DATA ZETA2/1.644934/,ZETA3/1.202057/
      DATA MZ/91.187/
C
C          Re-initialize weak scale parameters
C
      XLAMGM=M0
      XMESGM=MHF
      XN5GM=A0
      PI=4.*ATAN(1.)
      DPI=4.D0*DATAN(1.D0)
      CF=4./3.
      CA=3.
C     Here we input alpha_s^MSbar(MZ)
      ASMZ=0.1172
      MTAMZ=1.7463
C     Value of mb(MZ)^DRbar taken from PRD66, 074007 (2002).
      MBMZ=2.83
      SINB=SIN(ATAN(TANB))
      COSB=COS(ATAN(TANB))
      COTB=1./TANB
C
C     Calculate fermion masses including loop corrections
C
      M2=G0(8)
      AM2=ABS(M2)
      MSN=MSS(16)
      MG=ABS(MSS(1))
      MT1=MSS(12)
      MT2=MSS(13)
      MB1=MSS(10)
      MB2=MSS(11)
      MW1=ABS(MSS(27))
      MW2=ABS(MSS(28))
      AMU=ABS(MU)
      MUS=MU**2
      XLAM=DLOG(DBLE(HIGFRZ**2))
C      FTMT=MTMT/V
C      FBMZ=MBMZ/VP
C      FTAMZ=MTAMZ/VP
C          Be careful in using our convention vs Pierce et al.
C          cos(th)>1/sqrt(2):eigenstates same; cs-> -cs
C          cos(th)<1/sqrt(2);flip mass eigenstates; c <-> s interchange
C          Formula remains invariant under these switches
C          Use negative gaugino masses for consistency
C          Now input self-energies for consistency: RSIG*
C
C     here, add in 2 loop QCD correction to mt(DRbar) from Bednyakov
C     et al. Eq. 61.
      ST2LP=CF*(ASMTP/4./PI)**2*(-43.-12*ZETA2+CF*(-59./8.+30*ZETA2-
     ,48*LOG(2.)*ZETA2+12*ZETA3)+
     ,CA*(1093./24.-8*ZETA2+24*LOG(2.)*ZETA2-6*ZETA3))
      MTMT=MT/(1.+5*ASMTP/3./PI+ST2LP)
      FTMT=MTMT/SINB/VEV
      FBMZ=MBMZ/COSB/VEV
      FTAMZ=MTAMZ/COSB/VEV
      RSIGT=SSRSGT(MT**2)
      RSIGB=SSRSGB(MBQ**2)
      RSIGL=SSRSGL(MLQ**2)
C
C     Here, conversion from MSbar to DRbar is done at MZ.
C     Effect of sparticles is included by decoupling
C     beta functions in RGEs
      DAS=DBLE(ASMZ)/2.D0/DPI*(.5)
      ALEM=1./137.036
      DAEM=0.0682-ALEM/2./PI*(-7*LOG(AMW/AMZ))
      ALEMDR=ALEM/(1.-DAEM)
      XWMSB=.23113
C      XW=.2324-1.03E-7*(AMT**2-138.**2)
C     Convert XW to DRbar
      XW=1.-(1.-XWMSB)/(1.-ALEMDR/12./PI)      
      A1MZ=5*ALEMDR/3./(1.-XW)
      A2MZ=ALEMDR/XW
C      ALEM=1./128.
C      A1MZ=5*ALEM/3./(1.-XW)
C      A2MZ=ALEM/XW
      G(1)=DSQRT(4*DPI*A1MZ)
      G(2)=DSQRT(4*DPI*A2MZ)
      G(3)=DSQRT(4*DPI*ASMZ/(1.D0-DAS))
      G(4)=DBLE(FTAMZ)
      G(5)=DBLE(FBMZ)
      G(6)=G(6)
      G(25)=DBLE(MU)
      G(26)=DBLE(B)
      G(27)=0.D0
      G(28)=0.D0
      G(29)=0.D0
      G(30)=DBLE(VP)
      G(31)=DBLE(V)
C          Compute gauge mediated threshold functions
      IF (IMODEL.EQ.2) THEN
        XLM=XLAMGM/XMESGM
        THRF=((1.D0+XLM)*(LOG(1.D0+XLM)-2*DDILOG(XLM/(1.D0+XLM))+
     ,        .5*DDILOG(2*XLM/(1.D0+XLM)))+
     ,       (1.D0-XLM)*(LOG(1.D0-XLM)-2*DDILOG(-XLM/(1.D0-XLM))+
     ,        .5*DDILOG(-2*XLM/(1.D0-XLM))))/XLM**2
        THRG=((1.D0+XLM)*LOG(1.D0+XLM)+(1.D0-XLM)*LOG(1.D0-XLM))/XLM**2
      END IF
C
C          Run back up to mgut with approximate susy spectra
C
      IF (IMODEL.EQ.1) THEN
        IF (XSUGIN(7).EQ.0.) THEN 
          MGUT=1.E19
        ELSE
          MGUT=XSUGIN(7)
        END IF
      ELSE IF (IMODEL.EQ.2) THEN
        MGUT=XMESGM
      END IF
      TZ=DLOG(DBLE(MZ)/DBLE(MGUT))
      TGUT=0.D0
      DT=(TGUT-TZ)/DBLE(FLOAT(NSTEP))
      Q=MZ
      DO 250 II=1,NSTEP
        T=TZ+(TGUT-TZ)*FLOAT(II-1)/DBLE(FLOAT(NSTEP))
        Q=SNGL(MGUT*DEXP(T))
        QNEW=SNGL(MGUT*DEXP(T+DT))
        IF (Q.LE.MT.AND.QNEW.GT.MT) G(6)=DBLE(FTMT)
C       Implement sparticle threshold corrections at Q=HIGFRZ
        IF (Q.LE.HIGFRZ.AND.QNEW.GT.HIGFRZ) THEN
          G(6)=G(6)/(1.D0-DBLE(RSIGT))
          G(5)=G(5)/(1.D0-DBLE(RSIGB))
          G(4)=G(4)/(1.D0-DBLE(RSIGL))
          IF (INUHM.EQ.1) THEN
            G(13)=DBLE(MHDSQ)
            G(14)=DBLE(MHUSQ)
          END IF
        END IF
        IF (Q.LE.XNRIN(2).AND.QNEW.GT.XNRIN(2)) THEN
          G(27)=DBLE(FNMZ)
          G(28)=DBLE(G0(28))
          G(29)=DBLE(G0(29))
        END IF
        CALL DRKSTP(31,DT,T,G,SURG26,W2)
        A1I=SNGL(4*DPI/G(1)**2)
        A2I=SNGL(4*DPI/G(2)**2)
        A3I=SNGL(4*DPI/G(3)**2)
C       TEST YUKAWA DIVERGENCE
        IF (G(4).GT.5.D0.OR.G(5).GT.5.D0.OR.
     $G(6).GT.5.D0.OR.G(27).GT.5.D0) THEN
          NOGOOD=4
          GO TO 100
        END IF
        IF (A1I.LT.A2I.AND.XSUGIN(7).EQ.0.) GO TO 30
250   CONTINUE
      IF (IMODEL.EQ.1.AND.XSUGIN(7).EQ.0.) THEN
        WRITE(LOUT,*) 'SUGRGE ERROR: NO UNIFICATION FOUND'
        NOGOOD=1
        GO TO 100
      END IF
30    IF (XSUGIN(7).EQ.0.) THEN
        MGUT=QNEW
      ELSE
        MGUT=XSUGIN(7)
      END IF
      AGUT=SNGL((G(1)**2/4.D0/DPI+G(2)**2/4.D0/DPI)/2.D0)
      GGUT=SQRT(4*PI*AGUT)
      AGUTI=1./AGUT
      FTAGUT=SNGL(G(4))
      FBGUT=SNGL(G(5))
      FTGUT=SNGL(G(6))
      IF (INUHM.EQ.1) THEN
        MHDSMG=SNGL(G(13))
        MHUSMG=SNGL(G(14))
      END IF
      MUMG=SNGL(G(25))
      BMG=SNGL(G(26))
      IF (XNRIN(2).LT.1.E19.AND.XNRIN(1).EQ.0.) THEN
C     IMPOSE FN-FT UNIFICATION
        FNGUT=SNGL(G(6))
      ELSE
        FNGUT=SNGL(G(27))
      END IF
      G3GUT=SNGL(G(3))
      MGUTSS=MGUT
      AGUTSS=AGUT
      GGUTSS=GGUT
C
C          Set GUT boundary condition
C
      DO 260 I=1,3
        IF (IMODEL.EQ.1) THEN
          G(I)=G(I)
          G(I+6)=DBLE(MHF)
          G(I+9)=DBLE(A0)
        ELSE IF (IMODEL.EQ.2) THEN
          G(I)=G(I)
          G(I+6)=XGMIN(11+I)*XGMIN(8)*THRG*(G(I)/4.D0/DPI)**2*XLAMGM
          G(I+9)=0.D0
        END IF
      IF (XNRIN(2).LT.1.E19) THEN
        G(27)=DBLE(FNGUT)
        G(28)=DBLE(XNRIN(4))**2
        G(29)=DBLE(XNRIN(3))
      ELSE
        G(27)=0.D0
        G(28)=0.D0
        G(29)=0.D0
      END IF
260   CONTINUE
C     OVERWRITE ALFA_3 UNIFICATION TO GET ALFA_3(MZ) RIGHT
      IF (IMODEL.EQ.1.AND.IAL3UN.NE.0) G(3)=DBLE(GGUT)
      IF (IMODEL.EQ.1) THEN
        DO 270 I=13,24
          G(I)=DBLE(M0)**2
270     CONTINUE
      IF (INUHM.EQ.1) THEN
        G(13)=DBLE(MHDSMG)
        G(14)=DBLE(MHUSMG)
      END IF
C          Set possible non-universal GUT scale boundary conditions
      DO 280 I=1,6
        IF (XNUSUG(I).LT.1.E19) THEN
          G(I+6)=DBLE(XNUSUG(I))
        END IF
280   CONTINUE
      DO 281 I=7,18
        IF (XNUSUG(I).LT.1.E19) THEN
          G(I+6)=SIGN(1.,XNUSUG(I))*DBLE(XNUSUG(I))**2
        END IF
281   CONTINUE
      ELSE IF (IMODEL.EQ.2) THEN
       XC=2*THRF*XLAMGM**2
       DY=DSQRT(3.D0/5.D0)*G(1)*XGMIN(11)
       G(13)=XC*(.75*XGMIN(13)*(G(2)/4.D0/DPI)**4+.6D0*.25*
     , XGMIN(12)*(G(1)/4.D0/DPI)**4)+DBLE(XGMIN(9))-DY
       G(14)=XC*(.75*XGMIN(13)*(G(2)/4.D0/DPI)**4+.6D0*.25*
     , XGMIN(12)*(G(1)/4.D0/DPI)**4)+DBLE(XGMIN(10))+DY
       G(15)=XC*(.6*XGMIN(12)*(G(1)/4.D0/DPI)**4)+2*DY
       G(16)=XC*(.75*XGMIN(13)*(G(2)/4.D0/DPI)**4+.6D0*.25*
     , XGMIN(12)*(G(1)/4.D0/DPI)**4)-DY
       G(17)=XC*(4*XGMIN(14)*(G(3)/4.D0/DPI)**4/3.D0+.6D0*XGMIN(12)*
     , (G(1)/4.D0/DPI)**4/9.D0)+2*DY/3.D0
       G(18)=XC*(4*XGMIN(14)*(G(3)/4.D0/DPI)**4/3.D0+
     , .6D0*4*XGMIN(12)*(G(1)/4.D0/DPI)**4/9.D0)-4*DY/3.D0
       G(19)=XC*(4*XGMIN(14)*(G(3)/4.D0/DPI)**4/3.D0+.75*XGMIN(13)*
     ,(G(2)/4.D0/DPI)**4+.6*XGMIN(12)*(G(1)/4.D0/DPI)**4/36.D0)+DY/3.D0
       G(20)=G(15)
       G(21)=G(16)
       G(22)=G(17)
       G(23)=G(18)
       G(24)=G(19)
      ELSE IF (IMODEL.EQ.7.OR.IMODEL.EQ.9.OR.IMODEL.EQ.10) THEN
       G(1)=G(1)
       G(2)=G(2)
       G(3)=G(3)
       BLHAT=G(4)*(-9*G(1)**2/5.D0-3*G(2)**2+3*G(5)**2+4*G(4)**2)
       BBHAT=G(5)*(-7*G(1)**2/15.D0-3*G(2)**2-16*G(3)**2/3.D0+
     ,             G(6)**2+6*G(5)**2+G(4)**2)
       BTHAT=G(6)*(-13*G(1)**2/15.D0-3*G(2)**2-16*G(3)**2/3.D0+
     ,             6*G(6)**2+G(5)**2)
       G(7)=33*MHF*G(1)**2/5.D0/16.D0/DPI**2
       IF (IMODEL.EQ.10) THEN
         G(7)=G(7)+XAMIN(11)*MHF
       END IF
       G(8)=MHF*G(2)**2/16.D0/DPI**2
       G(9)=-3*MHF*G(3)**2/16.D0/DPI**2
       G(10)=-BLHAT*MHF/G(4)/16.D0/DPI**2
       G(11)=-BBHAT*MHF/G(5)/16.D0/DPI**2
       G(12)=-BTHAT*MHF/G(6)/16.D0/DPI**2
       G(13)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0+3*G(5)*BBHAT+
     ,G(4)*BLHAT)*MHF**2/(16*DPI**2)**2+XAMIN(6)*DBLE(M0)**2
       G(14)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0+3*G(6)*BTHAT)*
     ,        MHF**2/(16*DPI**2)**2+XAMIN(7)*DBLE(M0)**2
       G(15)=(-198*G(1)**4/25.D0)*MHF**2/(16*DPI**2)**2+
     ,XAMIN(5)*DBLE(M0)**2
       G(16)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0)*MHF**2/(16*DPI**2)**2
     ,+XAMIN(4)*DBLE(M0)**2
       G(17)=(-22*G(1)**4/25.D0+8*G(3)**4)*MHF**2/(16*DPI**2)**2+
     ,XAMIN(2)*DBLE(M0)**2
       G(18)=(-88*G(1)**4/25.D0+8*G(3)**4)*MHF**2/(16*DPI**2)**2+
     ,XAMIN(3)*DBLE(M0)**2
       G(19)=(-11*G(1)**4/50.D0-3*G(2)**4/2.D0+8*G(3)**4)*
     ,        MHF**2/(16*DPI**2)**2+XAMIN(1)*DBLE(M0)**2
       G(20)=(-198*G(1)**4/25.D0+2*G(4)*BLHAT)*MHF**2/(16*DPI**2)**2
     ,+XAMIN(5)*DBLE(M0)**2
       G(21)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0+G(4)*BLHAT)*
     ,        MHF**2/(16*DPI**2)**2+XAMIN(4)*DBLE(M0)**2
       G(22)=(-22*G(1)**4/25.D0+8*G(3)**4+2*G(5)*BBHAT)*
     , MHF**2/(16*DPI**2)**2+XAMIN(2)*DBLE(M0)**2
       G(23)=(-88*G(1)**4/25.D0+8*G(3)**4+2*G(6)*BTHAT)*
     , MHF**2/(16*DPI**2)**2+XAMIN(3)*DBLE(M0)**2
       G(24)=(-11*G(1)**4/50.D0-3*G(2)**4/2.D0+8*G(3)**4+G(5)*BBHAT+
     ,        G(6)*BTHAT)*MHF**2/(16*DPI**2)**2+XAMIN(1)*DBLE(M0)**2
      END IF
      IF (IMODEL.EQ.9) THEN
        CALL MMAMSB(M0,MHF,G)
      END IF
      DO 285 I=1,31
        IG(I)=0
285   CONTINUE
C          Check for tachyonic sleptons at GUT scale
      IF (G(15).LT.0.D0.OR.G(16).LT.0.D0) THEN
        ITACHY=2
      ELSE
        ITACHY=0
      END IF
C
C          Run back down to weak scale
C
      TZ=DLOG(DBLE(MZ)/DBLE(MGUT))
      TGUT=0.D0
      DT=(TZ-TGUT)/DBLE(FLOAT(NSTEP))
      DO 290 II=1,NSTEP+2
        T=TGUT+(TZ-TGUT)*FLOAT(II-1)/DBLE(FLOAT(NSTEP))
        Q=SNGL(MGUT*DEXP(T))
        CALL DRKSTP(31,DT,T,G,SURG26,W2)
C       Here, DRKSTP advances T by DT
        QNEW=SNGL(MGUT*DEXP(T))
C       TEST YUKAWA DIVERGENCE
        IF (G(4).GT.5.D0.OR.G(5).GT.5.D0.OR.
     $    G(6).GT.5.D0.OR.G(27).GT.5.D0) THEN
          NOGOOD=4
          GO TO 100
        END IF
        CALL SUGFRZ(QNEW,G,G0,IG)
        IF (Q.GE.AMNRMJ.AND.QNEW.LT.AMNRMJ.AND.XNRIN(1).EQ.0.) THEN
          FNMZ=SNGL(G(27))
        END IF
        IF (Q.GT.HIGFRZ.AND.QNEW.LE.HIGFRZ) THEN
          G(6)=G(6)*(1.D0-DBLE(RSIGT))
          G(5)=G(5)*(1.D0-DBLE(RSIGB))
          G(4)=G(4)*(1.D0-DBLE(RSIGL))
        END IF
        IF (QNEW.LT.AMNRMJ) THEN
          G(27)=0.D0
          G(28)=0.D0
          G(29)=0.D0
        END IF
        IF (NOGOOD.NE.0) GO TO 100
        IF (QNEW.LT.MZ) GO TO 40
290   CONTINUE
40    CONTINUE
C
C          Electroweak breaking constraints; tree level
C
      VUQ=G0(31)
      VDQ=G0(30)
      TANBQ=VUQ/VDQ
      SIN2BQ=SIN(2*ATAN(TANBQ))
      MZQ=SQRT((G0(2)**2+.6*G0(1)**2)*(VUQ**2+VDQ**2)/2.)
      BADMU=.FALSE.
      IF (INUHM.NE.1) THEN
      MUS=(G0(13)-G0(14)*TANBQ**2)/(TANBQ**2-1.)-MZQ**2/2.
C          Compute loop corrections using masses from last iteration
      CALL SUGEFF(G0,SIG1,SIG2)
      MH1S=G0(13)+SIG1
      MH2S=G0(14)+SIG2
      MUS=(MH1S-MH2S*TANBQ**2)/(TANBQ**2-1.)-MZQ**2/2.
C          If MUS<0, set it to MZ**2 and continue
      IF (MUS.LT.0.) THEN
        MUS=AMZ**2
      END IF
      MU=SQRT(MUS)*SIGN(1.,SGNMU)
      B=(G0(13)+G0(14)+2*MUS)*SIN2BQ/MU/2.
      CALL SUGMAS(G0,0,IMODEL,SIGA)
      IF (NOGOOD.NE.0) GO TO 100
C
C           Electroweak breaking constraints; loop level
C
      CALL SUGEFF(G0,SIG1,SIG2)
      MH1S=G0(13)+SIG1
      MH2S=G0(14)+SIG2
      MUS=(MH1S-MH2S*TANBQ**2)/(TANBQ**2-1.)-MZQ**2/2.
      IF (MUS.LT.0.) THEN
C        NOGOOD=2
C        GO TO 100
         MUS=MZ**2
      END IF
      MU=SQRT(MUS)*SIGN(1.,SGNMU)
      B=(MH1S+MH2S+2*MUS)*SIN2BQ/MU/2.
C
C     Once more, with feeling!
C
      CALL SUGEFF(G0,SIG1,SIG2)
      MH1S=G0(13)+SIG1
      MH2S=G0(14)+SIG2
      MUS=(MH1S-MH2S*TANBQ**2)/(TANBQ**2-1.)-MZQ**2/2.
      IF (MUS.LT.0.) THEN
C        NOGOOD=2
C        GO TO 100
         BADMU=.TRUE.
         MUS=MZ**2
      END IF
      MU=SQRT(MUS)*SIGN(1.,SGNMU)
      B=(MH1S+MH2S+2*MUS)*SIN2BQ/MU/2.
      CALL SUGMAS(G0,1,IMODEL,SIGA)
      ELSE
        MUS=MU**2
        B=AMHA**2/MU/(COTB+TANB)
        CALL SUGMAS(G0,1,IMODEL,SIGA)
        CALL SUGEFF(G0,SIG1,SIG2)
        MHDSQ=(-2*SIG1-MZQ**2-2*MUS+2*TANB**2*(AMHA**2+SIGA+MZQ**2+
     ,SIG2-SIG1)+TANB**4*(2*MUS-MZQ**2-2*SIG2-2*AMHA**2+2*SIGA))
     ,/2./(1.-TANB**4)
        MHUSQ=(2*AMHA**2+2*SIGA+MZQ**2-2*MUS-2*SIG1+2*TANB**2*(
     ,-AMHA**2+SIGA-MZQ**2-SIG1+SIG2)+TANB**4*(2*MUS+2*SIG2+MZQ**2))
     ,/2./(1.-TANB**4)
      END IF
C
C  Save radiative corrections to Yukawas for b->s gamma computation
C
      RTISA=RSIGT
      RBISA=RSIGB
      RLISA=RSIGL

100   RETURN
      END
