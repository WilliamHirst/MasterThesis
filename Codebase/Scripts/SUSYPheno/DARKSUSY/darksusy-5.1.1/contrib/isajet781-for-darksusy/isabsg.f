C
CDECK  ID>, ISABSG.
C--------------------------------------------------------------------
      subroutine ISABSG(IMODEL,M0,MHF,A0,BFBSG,IWR)
C--------------------------------------------------------------------
C    Calculates branching ratio of b->s\gamma decay.
C
C    Note: This version uses just running MBQ.
C   
C    Ref: H.Anlauf hep-ph/9406286;
C         C.Greub, T.Hurth and D.Wyler  PRD54 (1996);
C         A.Masiero et al  NPB353 (1991) 591-649;
C         M.Ciuchini et al  hep-ph/9304257
C
C    Created by H.Baer and M.Brhlik
C    Modified:
c      03/13/07 by Azar Mustafayev  - converted to double precision
c               and improved interface with ISASUGRA
c      05/29/07 by Azar Mustafayev  - neutrino sector added
c      12/13/07 by Azar Mustafayev  - fixes to speedup
c
c  NB: This version is design to work with ISAJET 7.75
c
C--------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IMODEL,IWR
      REAL M0,MHF,A0
      REAL*8 BFBSG
      COMMON /BSG/GISA(31),MSQISA(3),MSLISA(3),MSUISA(3),MSDISA(3),
     &            MSEISA(3),MRNISA(3),YNFRZ(3,3),MNFRZ(3,3),TNFRZ(3,3),
     &            RSIGT,RSIGB,RSIGL
c  GISA(i) - values of RGE parameters at MZ in DRbar:
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
c     
      REAL*8 G(157)
c   G(i) - values of RGE parameters at MZ in DRbar:
C     G(1) = g_1
C     G(2) = g_2
C     G(3) = g_3
C     G(4) = Y_u(1,1)
C     G(5) = Y_u(1,2)
C     G(6) = Y_u(1,3)
C     ...    ...
C     G(12) = Y_u(3,3)
C     G(13)-G(21) = Y_d
C     G(22)-G(30) = Y_e
C     G(31) = M_1
C     G(32) = M_2
C     G(33) = M_3
C     G(34)-G(42) = h_u
C     G(43)-G(51) = h_d
C     G(52)-G(60) = h_e
C     G(61)= mu
C     G(62)= B*mu
C     G(63)= m^2_Hu
C     G(64)= m^2_Hd
C     G(65)-G(73) = m^2_Q
C     G(74)-G(82) = m^2_L
C     G(83)-G(91) = m^2_u
C     G(92)-G(100)= m^2_d
C     G(101)-G(109)=m^2_e
C     G(110)= v_u
C     G(111)= v_d
C     G(112)-G(120) = Y_nu  
C     G(121)-G(129) = M_RHN 
C     G(130)-G(138) = h_nu
C     G(139)-G(147) = m^2_nu
C     G(148)-G(156) = \kappa
C     G(157) = \lambda   - SM quartic higgs coupling
      REAL*8 CKM(3,3),YU(3,3),YD(3,3),YE(3,3),HU(3,3),HD(3,3),HE(3,3),
     &       xYU(3,3),xYD(3,3),xHU(3,3),xHD(3,3),GF(157),xYUMT(3,3),
     &       VULMT(3,3),VURMT(3,3),
     &       YUGUT(3,3),YDGUT(3,3),YEGUT(3,3),YNGUT(3,3)
      REAL SQ,SQOLD,SFTMT
      REAL*8 HUGUT(9),HDGUT(9),HEGUT(9)
      REAL*8 CSM(9),CH1(9),CH2(19),
     $     CC111(12),CC112(12),CC113(12),
     $     CC121(12),CC122(12),CC123(12),
     $     CC211(17),CC212(17),CC213(17),
     $     CC221(17),CC222(17),CC223(17),
     $     DSSM(9),DSH1(9),DSH2(9),
     $     DSC111(9),DSC112(9),DSC113(9),
     $     DSC121(9),DSC122(9),DSC123(9),
     $     DSC211(9),DSC212(9),DSC213(9),
     $     DSC221(9),DSC222(9),DSC223(9),
     $     WSM(27),WHP1(27),WHP2(57),
     $     WSC111(36),WSC112(36),WSC113(36),
     $     WSC121(36),WSC122(36),WSC123(36),
     $     WSC211(51),WSC212(51),WSC213(51),
     $     WSC221(51),WSC222(51),WSC223(51),
     $     CI1(8),WCI1(24),CI2(8),WCI2(24),
     $     GK(2,3),HK(2,3)
      REAL*8 GAM(6,6),W2(471)

      REAL*8 AGE(4,6,6),AHE(4,6,6)
      INTEGER IGF(157)
      EXTERNAL RGE157,GAMMASM,GAMMAHP,GAMMAC1,GAMMAC2,GAMMAWB1,GAMMAWB2
      REAL SUALFE,SUALFS
      REAL*8 DDILOG,EI,BI,CI,GES,FES
C
      COMMON/BSGSM/ MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      REAL*8 MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      SAVE /BSGSM/
      COMMON /BSGDEC/MSQDEC(3),MSLDEC(3),MSUDEC(3),MSDDEC(3),
     &               MSEDEC(3),MRNDEC(3),IRHN
      REAL*8 MSQDEC,MSLDEC,MSUDEC,MSDDEC,MSEDEC,MRNDEC
      INTEGER IRHN
      SAVE /BSGDEC/
      COMMON /BSGSUG/ TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,
     &                MSCHR,MSUPL,MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS(4),
     &                ZMIXSS(4,4),AMW1SS,AMW2SS,GAMMAL,GAMMAR,THETAT,
     &                MTQ,MBQ,MSTLQ,MSTRQ,MGUT,FNGUT,FTMT,XRHNIN(21),
     &                XGMIN(14),GGUTSS,XNUSUG(20),XAMIN(7),EPSNU,
     &                FTRHLD(3),MHUSMG,MHDSMG,
     &                INUHM,IAL3UN,LND5ON
      REAL*8 TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,MSCHR,MSUPL,
     &       MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS,ZMIXSS,AMW1SS,AMW2SS,
     &       GAMMAL,GAMMAR,THETAT,MTQ,MBQ,MSTLQ,MSTRQ,MGUT,
     &       FNGUT,FTMT,XRHNIN,XGMIN,GGUTSS,XNUSUG,XAMIN,EPSNU,FTRHLD,
     &       MHUSMG,MHDSMG
      INTEGER IAL3UN,INUHM
      LOGICAL LND5ON
      SAVE /BSGSUG/
      COMMON /CHRGN/ MCHA(2),MSQU(3),MSQD(3),SMIX(6),RMIX(2),ABMIX(2)
c       MCHA - chargino masses
c       MSQU - up-squark masses
c       MSQD - down-squark masses
      REAL*8 MCHA,MSQU,MSQD,SMIX,RMIX,ABMIX
      SAVE /CHRGN/        
      COMMON /GLNN/ MCH0(4),MDSB(6)
c       MCH0 - neutralino masses
c       MDSB - 
      REAL*8 MCH0,MDSB
      COMMON /GGN/ M1,M2,M3,ABOT,ATOP,ATAU
c        M1,M2,M3  - gaugino masses at MZ
c        ABOT,ATOP,ATAU  - soft trilinear scalar couplings at MZ 
      REAL*8 M1,M2,M3,ABOT,ATOP,ATAU
      SAVE /GGN/      
      COMMON /G3/ G3
c        G3 = g_3 - strong coupling at current stage of RGE evolution
      REAL*8 G3
      SAVE /G3/
      COMMON /SCALEQ/ Q
      REAL*8 Q
      SAVE /SCALEQ/
C
      REAL*8 PI,SR2,C12,C13,C23,
     &       XLAMGM,XMESGM,XN5GM,XLM,THRF,THRG,BETA,SINB,COSB,VEV,
     &       MELMZ,MMUMZ,MUPMZ,MDOMZ,MCHMZ,MSTMZ,FELMZ,FMUMZ,FUPMZ,
     &       FDOMZ,FCHMZ,FSTMZ,TZ,TGUT,T,DT,QOLD,XC,DY,BLHAT,BBHAT,
     &       BTHAT,QNEW,COS2B
      REAL*8 VULQ(3,3),VURQ(3,3),VDLQ(3,3),VDRQ(3,3)
      REAL*8 FTAGUT,FBGUT,FTGUT,FELGUT,FMUGUT,FDOGUT,FSTGUT,FUPGUT,
     &       FCHGUT
      REAL*8 X,Y,FN1,FN2,FN3,FN4,YPS,SCALE,XWA,XCB,ALES,ALEL,G3MHP,MSQ,
     &       MCH,G3SS,QMAX,SWI1,SWI2,SWI3,CER,CEG,ALFA2,RALPH,FACG7,
     &       FACG8,FACN7,FACN8,SINTW2I,XTW,ASMB,AEMB,G2,GEF,BQLOG,ZI,
     &       LOGZI,RER2,RER7,RER8,REDE,IMR2,IMR7,IMR8,IMDE,EF,GVIRT,
     &       GAMMASG,gammasg0,GAMMASL,GBREMS,XWB
      REAL*8 c3smnai,c2smnai,c3hpnai,c3c11nai,c3c12nai,c3c13nai,
     &       c3c23nai,C2SMEFF,C3SMEFF,C2HPEFF,C3HPEFF,C2C11EFF,C3C11EFF,
     &       C2C12EFF,C3C12EFF,C2C13EFF,C3C13EFF,C2C21EFF,C3C21EFF,
     &       C2C22EFF,C3C22EFF,C2C23EFF,C3C23EFF,C3GS1EFF,C2GS1EFF,
     &       C3GS2EFF,C2GS2EFF,C3GS3EFF,C2GS3EFF,C3GS4EFF,C2GS4EFF,
     &       C3GS5EFF,C2GS5EFF,C3GS6EFF,C2GS6EFF,C3NS1EFF,C2NS1EFF,
     &       C3NS2EFF,C2NS2EFF,C3NS3EFF,C2NS3EFF,C3NS4EFF,C2NS4EFF,
     &       C3NS5EFF,C2NS5EFF,C3NS6EFF,C2NS6EFF,C3CHEFF,C3GLEFF,
     &       C3NEEFF,C2CHEFF,C2GLEFF,C2NEEFF,C3C1EFF,C3C2EFF,
     &       C3C1NAI,C3C2NAI,C3CHNAI,C3GLNAI,C3NENAI,c3c21nai,c3c22nai
      REAL*8 BRSL,COBA,BRSG,BRSG0,brup,brdo,BFBSG0
      REAL XQ2,XMT
      INTEGER NSTEP,I,J,K,L,II
      DATA NSTEP/30000/
      
      real*8 ESQ(3,3),ESD(3,3),EAD(3,3),EMD(3,3)
      
      PI=4.d0*DATAN(1.d0)
      SR2=SQRT(2.d0)
      
c...Set experimentaly measured parameters
      BRSL=0.104    ! semileptonic branching ratio
      COBA=0.95     ! CKM factor

c...Reset the branching fraction             
      BFBSG=-1.d0

c...Print logo
      IF (IWR.EQ.1) THEN
        PRINT*,'%%%%%%%%%%%%%%%%%%% ISABSG %%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      ENDIF

c...Transfer input parameters
      CALL SUG2BSG

C...DEFINE CKM MATRIX
      C12=SQRT(1.d0-S12**2)
      C13=SQRT(1.d0-S13**2)
      C23=SQRT(1.d0-S23**2)
      CKM(1,1)= C12*C13
      CKM(1,2)= S12*C13
      CKM(1,3)= S13
      CKM(2,1)=-S12*C23-C12*S23*S13
      CKM(2,2)= C12*C23-S12*S23*S13
      CKM(2,3)= S23*C13
      CKM(3,1)= S12*S23-C12*C23*S13
      CKM(3,2)=-C12*S23-S12*C23*S13
      CKM(3,3)= C23*C13
      
C---------------------------------------------
      CALL CHARGINO(GK,HK)

C...Compute gauge mediated threshold functions
      IF (IMODEL.EQ.2) THEN
        XLAMGM=M0
        XMESGM=MHF
        XN5GM=A0
        XLM=XLAMGM/XMESGM
        THRF=((1.D0+XLM)*(LOG(1.D0+XLM)-2*DDILOG(XLM/(1.D0+XLM))+
     ,        .5*DDILOG(2*XLM/(1.D0+XLM)))+
     ,       (1.D0-XLM)*(LOG(1.D0-XLM)-2*DDILOG(-XLM/(1.D0-XLM))+
     ,        .5*DDILOG(-2*XLM/(1.D0-XLM))))/XLM**2
        THRG=((1.D0+XLM)*LOG(1.D0+XLM)+(1.D0-XLM)*LOG(1.D0-XLM))/XLM**2
      END IF

C...Assemble (s)quark matrices and transform into correct basis if necessary
      BETA=ATAN(TANB)
      SINB=SIN(BETA)
      COSB=COS(BETA)
      COS2B=COS(2.d0*BETA)
      VEV=SQRT(V**2+VP**2)
      DO I=1,3
        DO J=1,3
	  xYU(I,J)=0.d0
	  xYD(I,J)=0.d0
	  xHU(I,J)=0.d0
	  xHD(I,J)=0.d0
          YU(I,J)=0.d0
          YD(I,J)=0.d0
          YE(I,J)=0.d0
          HU(I,J)=0.d0
          HD(I,J)=0.d0
          HE(I,J)=0.d0
        ENDDO
      ENDDO
      
      xYU(1,1)= 0.0027 /SINB/VEV 
      xYU(2,2)= 0.620  /SINB/VEV
C      xYU(3,3)= GISA(6)
      xYU(3,3)= 0.D0
      DO I=1,3
        DO J=1,3
          DO K=1,3
            YU(I,J)=YU(I,J)+CKM(K,I)*xYU(K,J)
          ENDDO
        ENDDO
      ENDDO
      HU(3,3)= GISA(12)*GISA(6)
      YD(1,1)= 0.0058 /COSB/VEV 
      YD(2,2)= 0.113  /COSB/VEV 
      YD(3,3)= GISA(5)
      HD(3,3)= GISA(11)*GISA(5)
C...(s)lepton sector
      YE(1,1)= 0.00051/COSB/VEV 
      YE(2,2)= 0.106  /COSB/VEV
      YE(3,3)= GISA(4)
      HE(3,3)= GISA(10)*GISA(4)

C...Set RGE parameters at M_Z in DRbar      
c   NOTE difference in Yukawa and trilinear definitions!      
      G(1)=GISA(1)
      G(2)=GISA(2)
      G(3)=GISA(3)

      K=0
      DO I=1,3
        DO J=1,3
          G(4+K) = YU(J,I)
          G(13+K)= YD(J,I)
          G(22+K)= YE(J,I)
          G(34+K)= HU(J,I)
          G(43+K)= HD(J,I)
          G(52+K)= HE(J,I)
          K=K+1
        ENDDO
      ENDDO
      
      G(31)=GISA(7)
      G(32)=GISA(8)
      G(33)=GISA(9)
      
      DO I=65,157
        G(I)=0.d0
      ENDDO
      G(61)=GISA(25)
      G(62)=GISA(26)*GISA(25)
      G(63)=GISA(14)
      G(64)=GISA(13)
      
      G(65)=GISA(19)
      G(69)=GISA(19)
      G(73)=GISA(24)
      
      G(74)=GISA(16)
      G(78)=GISA(16)
      G(82)=GISA(21)
      
      G(83)=GISA(18)
      G(87)=GISA(18)
      G(91)=GISA(23)
      
      G(92)=GISA(17)
      G(96)=GISA(17)
      G(100)=GISA(22)
      
      G(101)=GISA(15)
      G(105)=GISA(15)
      G(109)=GISA(20)
      
      G(110)=GISA(31)
      G(111)=GISA(30)
      
      G(120)=GISA(27)
      G(129)=GISA(28)
      G(138)=GISA(29)*GISA(27)
      G(147)=0.d0
      G(156)=0.d0
      G(157)=0.5
C-------------------------------------------------------------------
C----- RUN RGEs FROM MZ TO MGUT -------------------------------------
C-------------------------------------------------------------------
      TZ=DLOG(MZ/MGUT)
      TGUT=0.d0
      DT=(TGUT-TZ)/FLOAT(NSTEP)
      DO II=1,NSTEP
        T=TZ+DT*FLOAT(II-1)
        QOLD=Q
        Q=MGUT*EXP(T)
c...Correct up Yukawa matrix at M_top
        IF (QOLD.LE.MT.AND.Q.GT.MT) THEN
          SFTMT=SNGL(FTMT)
	  CALL TACTIV(4,SFTMT,G)
	ENDIF
C...Implement sparticle threshold corrections at Q=HIGFRZ
        IF (QOLD.LE.MSUSY.AND.Q.GT.MSUSY) THEN
	  CALL YUKDIAG(G,4,xYU,VULQ,VURQ)
	  xYU(3,3)=xYU(3,3)/(1.d0-DBLE(RSIGT))
	  CALL AYUKDIAG(G,4,xYU,VULQ,VURQ)
	  CALL YUKDIAG(G,13,xYD,VDLQ,VDRQ)
	  xYD(3,3)=xYD(3,3)/(1.d0-DBLE(RSIGB))
	  CALL AYUKDIAG(G,13,xYD,VDLQ,VDRQ)
	  G(30)=G(30)/(1.d0-DBLE(RSIGL))
      CALL GLUNENO(G,GAM,AGE,AHE)
        END IF
c...neutrino sector
	IF (IRHN.GT.0) THEN
	  IF (QOLD.LE.DABS(MRNDEC(3)).AND.Q.GT.DABS(MRNDEC(3))) THEN
            G(120)=GISA(27)
            G(129)=GISA(28)
            G(138)=GISA(29)*GISA(27)
	  ENDIF
        ENDIF
        CALL DRKSTP(157,DT,T,G,RGE157,W2)   
      END DO
C----------------------------------------------------------------------
C----- SET GUT-scale boundary conditions ------------------------------
C----------------------------------------------------------------------

c...Yukawa couplings      
      CALL VEC2MAT(G,4,YUGUT,-1)
      CALL VEC2MAT(G,13,YDGUT,-1)
      CALL VEC2MAT(G,22,YEGUT,-1)
      IF (IRHN.EQ.1)  CALL VEC2MAT(G,112,YNGUT,-1)
      
      CALL BSGGUT(YUGUT,YDGUT,YEGUT,YNGUT,G,IMODEL,M0,MHF,A0)
c...reset freezeout flags
      DO I=1,157
        IF(I.LT.148) THEN 
	  IGF(I)=0
	ELSE
	  IGF(I)=1
	ENDIF
	GF(I)=0.d0
      ENDDO
C----------------------------------------------------------------------
C-------- SET VALUES OF WILSON COEFFICIENTS ---------------------------
C----------------------------------------------------------------------
      CALL WILSON(GK,HK,
     &  	  CSM,CH1,CH2,CC111,CC112,CC113,CC121,CC122,CC123,
     &  	  CC211,CC212,CC213,CC221,CC222,CC223,
     &  	  DSSM,DSH1,DSH2,DSC111,DSC112,DSC113,DSC121,
     &  	  DSC122,DSC123,DSC211,DSC212,DSC213,DSC221,
     &  	  DSC222,DSC223)
C
      c3smnai=  -1./2.*(CSM(1)+DSSM(1))
     $          +(CSM(3)+DSSM(3))
     $          -1./2.*(CSM(4)+DSSM(4)+2./9.d0)
     $          -1./4.*(CSM(5)+DSSM(5)-7./9.d0)
     $          +1./4.*(CSM(7)+DSSM(7)+1.d0)
     $          -1./4.*(CSM(9)+DSSM(9)+9.d0)          
      c2smnai=  -1./2.*(CSM(1)+DSSM(1))
     $          +(CSM(2)+DSSM(2))
     $          -1./2.*(CSM(4)+DSSM(4)+2./9.d0)
     $          -1./4.*(CSM(5)+DSSM(5)-7./9.d0)
     $          +1./4.*(CSM(7)+DSSM(7)+1.d0)         
      c3hpnai=  -1./2.d0*(CH1(1)+DSH1(1))
     $          +(CH1(3)+DSH1(3))
     $          -1./2.d0*(CH1(4)+DSH1(4))
     $          -1./4.d0*(CH1(5)+DSH1(5))
     $          +1./4.d0*(CH1(7)+DSH1(7))
     $          -1./4.d0*(CH1(9)+DSH1(9))        
      c3c11nai= -1./2.d0*(CC211(1)+DSC211(1))                                
     $          +(CC211(3)+DSC211(3))
     $          -1./2.d0*(CC211(4)+DSC211(4))
     $          -1./4.d0*(CC211(5)+DSC211(5))
     $          +1./4.d0*(CC211(7)+DSC211(7))
     $          -1./4.d0*(CC211(9)+DSC211(9))                  
      c3c12nai= -1./2.d0*(CC212(1)+DSC212(1))
     $          +(CC212(3)+DSC212(3))
     $          -1./2.d0*(CC212(4)+DSC212(4))
     $          -1./4.d0*(CC212(5)+DSC212(5))
     $          +1./4.d0*(CC212(7)+DSC212(7))
     $          -1./4.d0*(CC212(9)+DSC212(9))                  
      c3c13nai= -1./2.d0*(CC213(1)+DSC213(1))
     $          +(CC213(3)+DSC213(3))
     $          -1./2.d0*(CC213(4)+DSC213(4))
     $          -1./4.d0*(CC213(5)+DSC213(5))
     $          +1./4.d0*(CC213(7)+DSC213(7))
     $          -1./4.d0*(CC213(9)+DSC213(9))                    
      c3c21nai= -1./2.d0*(CC221(1)+DSC221(1))
     $          +(CC221(3)+DSC221(3))
     $          -1./2.d0*(CC221(4)+DSC221(4))
     $          -1./4.d0*(CC221(5)+DSC221(5))
     $          +1./4.d0*(CC221(7)+DSC221(7))
     $          -1./4.d0*(CC221(9)+DSC221(9))                         
      c3c22nai= -1./2.d0*(CC222(1)+DSC222(1))
     $          +(CC222(3)+DSC222(3))
     $          -1./2.d0*(CC222(4)+DSC222(4))
     $          -1./4.d0*(CC222(5)+DSC222(5))
     $          +1./4.d0*(CC222(7)+DSC222(7))
     $          -1./4.d0*(CC222(9)+DSC222(9))                               
      c3c23nai= -1./2.d0*(CC223(1)+DSC223(1))
     $          +(CC223(3)+DSC223(3))
     $          -1./2.d0*(CC223(4)+DSC223(4))
     $          -1./4.d0*(CC223(5)+DSC223(5))
     $          +1./4.d0*(CC223(7)+DSC223(7))
     $          -1./4.d0*(CC223(9)+DSC223(9))                           
C----------------------------------------------------------------------
C------- NOW WE RUN from MGUT past MZ to Q=1 GeV ----------------------
C----------------------------------------------------------------------
      TZ=LOG(1.d0/MGUT)
      DT=TZ/FLOAT(NSTEP)
C
      DO 400 II=1,NSTEP+2   
        T=    DT*FLOAT(II-1)
        Q   = MGUT*EXP(T)
        QOLD= MGUT*EXP(T-dt)
        QNEW= MGUT*EXP(T+dt)
C
        RMIX(1)=0.d0
        RMIX(2)=0.d0 
        SMIX(1)=0.d0
        SMIX(2)=0.d0  
        SMIX(3)=0.d0
        SMIX(4)=0.d0  
        SMIX(5)=0.d0
        SMIX(6)=0.d0  
        ABMIX(1)=0.d0
        ABMIX(2)=0.d0      
C
C-------- Evolution of MSSM parameters --------------------------------
C
        CALL DRKSTP(157,DT,T,G,RGE157,W2)
C...TEST YUKAWA DIVERGENCE
        DO I=4,30
	  IF (G(i).gt.5.d0) THEN
            PRINT*,'ISABSG: NON-PERTURBATIVE YUKAWA'
	    PRINT*,'G(',i,')>= 5:',G(i),' Q=',Q
          END IF
	ENDDO
c...Decouple top quark in  up Yukawa elements at M_top
        IF (Q.LT.MT.AND.QOLD.GE.MT) THEN
	  CALL YUKDIAG(G,4,xYUMT,VULMT,VURMT)
	  xYUMT(3,3)=0.d0
	  CALL AYUKDIAG(G,4,xYUMT,VULMT,VURMT)
	ENDIF
C...Implement sparticle threshold corrections at Q=HIGFRZ
        IF (QOLD.GT.MSUSY.AND.Q.LE.MSUSY) THEN
      CALL GLUNENO(G,GAM,AGE,AHE)
	  CALL YUKDIAG(G,4,xYU,VULQ,VURQ)
	  xYU(3,3)=xYU(3,3)*(1.d0-DBLE(RSIGT))
	  CALL AYUKDIAG(G,4,xYU,VULQ,VURQ)
	  CALL YUKDIAG(G,13,xYD,VDLQ,VDRQ)
	  xYD(3,3)=xYD(3,3)*(1.d0-DBLE(RSIGB))
	  CALL AYUKDIAG(G,13,xYD,VDLQ,VDRQ)
	  G(30)=G(30)*(1.d0-DBLE(RSIGL))
        END IF
c...decoupling in neutrino sector
        IF (QNEW.LT.DABS(MRNDEC(3))) THEN
          G(120)=0.d0
          G(129)=0.d0
          G(138)=0.d0
	ENDIF
C
C------- Evolution of Wilson coefficients -----------------------------
C
        G3=G(3)
        ALES=G(3)**2/4.d0/PI
        ALEL=DBLE(SUALFE(SNGL(MW**2)))
C
c ...  SM contribution
C
       IF(Q.LT.MT.AND.Q.GT.MW) THEN
        T=DT*FLOAT(II-1)
        CALL DRKSTP(9,DT,T,CSM,GAMMASM,WSM)
       ENDIF

       IF(Q.GT.MW.AND.QNEW.LT.MW) THEN  
         CSM(1)=CSM(1)+DSSM(1)
         CSM(2)=CSM(2)+DSSM(2)
         CSM(3)=CSM(3)+DSSM(3)       
         CSM(4)=CSM(4)+2./9.d0+DSSM(4)
         CSM(5)=CSM(5)-7./9.d0+DSSM(5)
         CSM(6)=CSM(6)+2./9.d0+DSSM(6)
         CSM(7)=CSM(7)+1.d0+DSSM(7)
         CSM(8)=CSM(8)-3./2.d0+DSSM(8)
         CSM(9)=CSM(9)+9.d0+DSSM(9)             
       ENDIF
        C2SMEFF=-1./2.d0*CSM(1)+CSM(2)-1./2.d0*CSM(4)-1./4.d0*CSM(5)
     $          +1./4.d0*CSM(7)            
        C3SMEFF=-1./2.d0*CSM(1)+CSM(3)-1./2.d0*CSM(4)-1./4.d0*CSM(5)
     $          +1./4.d0*CSM(7)-1./4.d0*CSM(9) 
C
C ...  CH (charged higgs) contribution  
C
       IF(Q.GT.MHPLUS.AND.QNEW.LT.MHPLUS) THEN
        G3MHP=G(3)
        X=(MT/MHPLUS)**2
        CH2(11)=-1./2.d0*X/TANB**2
        CH2(17)= 16.*PI**2*X/G3MHP**2
       ENDIF
C       
       IF(MHPLUS.LT.MT) THEN
        IF(Q.LT.MT.AND.Q.GT.MW) THEN
         T=DT*FLOAT(II-1)
         CALL DRKSTP(9,DT,T,CH1,GAMMASM,WHP1)     
         IF(Q.GT.MHPLUS.AND.QNEW.LE.MHPLUS) THEN
          CH1(1)=CH1(1)+DSH1(1)   
          CH1(2)=CH1(2)+DSH1(2)
          CH1(3)=CH1(3)+DSH1(3)
          CH1(4)=CH1(4)+DSH1(4)
          CH1(5)=CH1(5)+DSH1(5)
          CH1(6)=CH1(6)+DSH1(6)
          CH1(7)=CH1(7)+DSH1(7)
          CH1(8)=CH1(8)+DSH1(8)
          CH1(9)=CH1(9)+DSH1(9)    
         ENDIF
         C2HPEFF=-1./2.d0*CH1(1)+CH1(2)-1./2.d0*CH1(4)-1./4.d0*CH1(5)
     $           +1./4.d0*CH1(7)      
         C3HPEFF=-1./2.d0*CH1(1)+CH1(3)-1./2.d0*CH1(4)-1./4.d0*CH1(5)
     $           +1./4.d0*CH1(7)-1./4.d0*CH1(9)       
        ENDIF
       ELSE
        IF(Q.LT.MHPLUS.AND.Q.GT.MT) THEN
         T=DT*FLOAT(II-1)
         CALL DRKSTP(19,DT,T,CH2,GAMMAHP,WHP2)     
          SMIX(1)=CH2(10)
          SMIX(2)=CH2(11)
          SMIX(3)=CH2(12)
          SMIX(4)=CH2(13)
          SMIX(5)=CH2(14)
          SMIX(6)=CH2(15)
C
         IF(Q.GT.MT.AND.QNEW.LE.MT) THEN
          CH2(1)=CH2(1)+DSH2(1) 
          CH2(2)=CH2(2)+DSH2(2)
          CH2(3)=CH2(3)+DSH2(3)
          CH2(4)=CH2(4)+DSH2(4)
          CH2(5)=CH2(5)+DSH2(5)
          CH2(6)=CH2(6)+DSH2(6)
          CH2(7)=CH2(7)+DSH2(7)
          CH2(8)=CH2(8)+DSH2(8)
          CH2(9)=CH2(9)+DSH2(9)    
         ENDIF
        ENDIF                           
C       
       IF(Q.LT.MT.AND.Q.GT.MW) THEN
         T=DT*FLOAT(II-1)
         CALL DRKSTP(9,DT,T,CH2,GAMMASM,WHP2)    
C
         C2HPEFF=-1./2.d0*CH2(1)+CH2(2)-1./2.d0*CH2(4)-1./4.d0*CH2(5)
     $           +1./4.*CH2(7)       
         C3HPEFF=-1./2.d0*CH2(1)+CH2(3)-1./2.d0*CH2(4)-1./4.d0*CH2(5)
     $           +1./4.d0*CH2(7)-1./4.d0*CH2(9) 
        
	ENDIF
       ENDIF
C
C ...  chargino(1)-squark(1) contribution  
C
       MSQ=MSQU(1)
       MCH=MCHA(1)    
C
       IF(Q.GT.MSQ.AND.QNEW.LE.MSQ) THEN
        G3SS=G(3)
        X=(MW/MSQU(1))**2
        Y=MCHA(1)/MBQ
        CC111(10)= GK(1,1)*HK(1,1)*16.*PI**2*X*Y/G3SS**2
        CC111(11)= GK(1,1)*HK(1,1)*16.*PI**2*X*Y/G3SS**2
        CC111(12)=-GK(1,1)*GK(1,1)*16.*PI**2*X/G3SS**2
       ENDIF    
C
       IF(MSQU(1).GT.MCHA(1)) THEN
         IF(MCHA(1).GT.(MW+0.5)) THEN
           QMAX=MCH
         ELSE
           QMAX=MW
         ENDIF
         IF(Q.LT.MSQ.AND.Q.GT.QMAX) THEN  
           T=DT*FLOAT(II-1)
           CALL DRKSTP(12,DT,T,CC111,GAMMAC1,WSC111)        
         ENDIF
         IF(Q.GT.QMAX.AND.QNEW.LT.QMAX) THEN
	   CC111(1)=CC111(1)+DSC111(1)        
           CC111(2)=CC111(2)+DSC111(2)
           CC111(3)=CC111(3)+DSC111(3)
           CC111(4)=CC111(4)+DSC111(4)
           CC111(5)=CC111(5)+DSC111(5)
           CC111(6)=CC111(6)+DSC111(6)
           CC111(7)=CC111(7)+DSC111(7)
           CC111(8)=CC111(8)+DSC111(8)
           CC111(9)=CC111(9)+DSC111(9)
 	 ENDIF
         IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
           T=DT*FLOAT(II-1)
           CALL DRKSTP(9,DT,T,CC111,GAMMASM,WSC111)       
         ENDIF
c
         C2C11EFF=-1./2.d0*CC111(1)+CC111(2)
     $            -1./2.d0*CC111(4)-1./4.d0*CC111(5)
     $            +1./4.d0*CC111(7)       
         C3C11EFF=-1./2.d0*CC111(1)+CC111(3)
     $            -1./2.d0*CC111(4)-1./4.d0*CC111(5)
     $            +1./4.d0*CC111(7)-1./4.d0*CC111(9)       
       ELSE 
        IF(Q.LT.MCH.AND.Q.GT.MSQ) THEN
         SWI1=1.
        ELSE
         SWI1=0.
        ENDIF              
        RMIX(1)=RMIX(1)+SWI1*CC211(13)
        RMIX(2)=RMIX(2)+SWI1*CC211(14)       
        IF(MSQU(1).GT.MT) THEN
         SMIX(1)=0.
         SMIX(2)=0.
         SMIX(3)=0.
         SMIX(4)=0.
         SMIX(5)=0.
         SMIX(6)=0.
        ENDIF
        IF(Q.LT.MCHA(1).AND.Q.GT.MSQU(2)) THEN
         SWI2=1.
        ELSE
         SWI2=0.
        ENDIF                
        IF(Q.LT.MCHA(1).AND.Q.GT.MSQU(3)) THEN
         SWI3=1.
        ELSE
         SWI3=0.
        ENDIF                
        ABMIX(1)=SWI2*CC212(13)+SWI3*CC213(13)
        ABMIX(2)=SWI2*CC212(14)+SWI3*CC213(14)
C 
        IF(MSQU(1).GT.(MW+0.5)) THEN
         QMAX=MSQ
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MCH.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(17,DT,T,CC211,GAMMAC2,WSC211)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LT.QMAX) THEN
         CC211(1)=CC211(1)+DSC211(1)        
         CC211(2)=CC211(2)+DSC211(2)
         CC211(3)=CC211(3)+DSC211(3)
         CC211(4)=CC211(4)+DSC211(4)
         CC211(5)=CC211(5)+DSC211(5)
         CC211(6)=CC211(6)+DSC211(6)
         CC211(7)=CC211(7)+DSC211(7)
         CC211(8)=CC211(8)+DSC211(8)
         CC211(9)=CC211(9)+DSC211(9)        
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC211,GAMMASM,WSC211)      
        ENDIF                  
         C2C11EFF=-1./2.d0*CC211(1)+CC211(2)
     $            -1./2.d0*CC211(4)-1./4.d0*CC211(5)
     $            +1./4.d0*CC211(7)      
         C3C11EFF=-1./2.d0*CC211(1)+CC211(3)
     $            -1./2.d0*CC211(4)-1./4.d0*CC211(5)
     $            +1./4.d0*CC211(7)-1./4.d0*CC211(9)       
       ENDIF
C
C ...  chargino(1)-squark(2) contribution  
C
       MSQ=MSQU(2)
       MCH=MCHA(1)    
c
       IF(Q.GT.MSQ.AND.QNEW.LE.MSQ) THEN
        G3SS=G(3)
        X=(MW/MSQU(2))**2
        Y=MCHA(1)/MBQ
        CC112(10)= GK(1,2)*HK(1,2)*16.*PI**2*X*Y/G3SS**2
        CC112(11)= GK(1,2)*HK(1,2)*16.*PI**2*X*Y/G3SS**2
        CC112(12)=-GK(1,2)*GK(1,2)*16.*PI**2*X/G3SS**2
       ENDIF  
C
       IF(MSQU(2).GT.MCHA(1)) THEN
        IF(MCHA(1).GT.(MW+0.5)) THEN
         QMAX=MCH
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MSQ.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(12,DT,T,CC112,GAMMAC1,WSC112)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LT.QMAX) THEN
         CC112(1)=CC112(1)+DSC112(1)        
         CC112(2)=CC112(2)+DSC112(2)
         CC112(3)=CC112(3)+DSC112(3)
         CC112(4)=CC112(4)+DSC112(4)
         CC112(5)=CC112(5)+DSC112(5)
         CC112(6)=CC112(6)+DSC112(6)
         CC112(7)=CC112(7)+DSC112(7)
         CC112(8)=CC112(8)+DSC112(8)
         CC112(9)=CC112(9)+DSC112(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC112,GAMMASM,WSC112)      
        ENDIF
         C2C12EFF=-1./2.d0*CC112(1)+CC112(2)
     $            -1./2.d0*CC112(4)-1./4.d0*CC112(5)
     $            +1./4.d0*CC112(7)       
         C3C12EFF=-1./2.d0*CC112(1)+CC112(3)
     $            -1./2.d0*CC112(4)-1./4.d0*CC112(5)
     $            +1./4.d0*CC112(7)-1./4.d0*CC112(9)
       ELSE 
        IF(Q.LT.MCH.AND.Q.GT.MSQ) THEN
         SWI1=1.
        ELSE
         SWI1=0.
        ENDIF              
        RMIX(1)=RMIX(1)+SWI1*CC212(13)
        RMIX(2)=RMIX(2)+SWI1*CC212(14)
        IF(MSQU(2).GT.MT) THEN
         SMIX(1)=0.
         SMIX(2)=0.
         SMIX(3)=0.
         SMIX(4)=0.
         SMIX(5)=0.
         SMIX(6)=0.
        ENDIF
        IF(Q.LT.MCHA(1).AND.Q.GT.MSQU(1)) THEN
         SWI2=1.
        ELSE
         SWI2=0.
        ENDIF                
        IF(Q.LT.MCHA(1).AND.Q.GT.MSQU(3)) THEN
         SWI3=1.
        ELSE
         SWI3=0.
        ENDIF                
        ABMIX(1)=SWI2*CC211(13)+SWI3*CC213(13)
        ABMIX(2)=SWI2*CC211(14)+SWI3*CC213(14)
C 
        IF(MSQU(2).GT.(MW+0.5)) THEN
         QMAX=MSQ
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MCH.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(17,DT,T,CC212,GAMMAC2,WSC212)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LE.QMAX) THEN
         CC212(1)=CC212(1)+DSC212(1)        
         CC212(2)=CC212(2)+DSC212(2)
         CC212(3)=CC212(3)+DSC212(3)
         CC212(4)=CC212(4)+DSC212(4)
         CC212(5)=CC212(5)+DSC212(5)
         CC212(6)=CC212(6)+DSC212(6)
         CC212(7)=CC212(7)+DSC212(7)
         CC212(8)=CC212(8)+DSC212(8)
         CC212(9)=CC212(9)+DSC212(9)        
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC212,GAMMASM,WSC212)      
        ENDIF                  
         C2C12EFF=-1./2.d0*CC212(1)+CC212(2)
     $            -1./2.d0*CC212(4)-1./4.d0*CC212(5)
     $            +1./4.d0*CC212(7)       
         C3C12EFF=-1./2.d0*CC212(1)+CC212(3)
     $            -1./2.d0*CC212(4)-1./4.d0*CC212(5)
     $            +1./4.d0*CC212(7)-1./4.d0*CC212(9) 
       ENDIF
C  
C ...  chargino(1)-squark(3) contribution  
C
       MSQ=MSQU(3)
       MCH=MCHA(1)    
       IF(Q.GT.MSQ.AND.QNEW.LE.MSQ) THEN
        G3SS=G(3)
        X=(MW/MSQU(3))**2
        Y=MCHA(1)/MBQ
        CC113(10)= HK(1,3)*16.*PI**2*X*Y/G3SS**2
        CC113(11)= HK(1,3)*16.*PI**2*X*Y/G3SS**2
        CC113(12)=-GK(1,3)*16.*PI**2*X/G3SS**2
       ENDIF    
C
       IF(MSQU(3).GT.MCHA(1)) THEN
        IF(MCHA(1).GT.(MW+0.5)) THEN
         QMAX=MCH
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MSQ.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(12,DT,T,CC113,GAMMAC1,WSC113)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LE.QMAX) THEN
         CC113(1)=CC113(1)+DSC113(1)        
         CC113(2)=CC113(2)+DSC113(2)
         CC113(3)=CC113(3)+DSC113(3)
         CC113(4)=CC113(4)+DSC113(4)
         CC113(5)=CC113(5)+DSC113(5)
         CC113(6)=CC113(6)+DSC113(6)
         CC113(7)=CC113(7)+DSC113(7)
         CC113(8)=CC113(8)+DSC113(8)
         CC113(9)=CC113(9)+DSC113(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC113,GAMMASM,WSC113)      
        ENDIF
         C2C13EFF=-1./2.d0*CC113(1)+CC113(2)
     $            -1./2.d0*CC113(4)-1./4.d0*CC113(5)
     $            +1./4.d0*CC113(7)       
         C3C13EFF=-1./2.d0*CC113(1)+CC113(3)
     $            -1./2.d0*CC113(4)-1./4.d0*CC113(5)
     $            +1./4.d0*CC113(7)-1./4.d0*CC113(9)       
       ELSE 
        IF(Q.LT.MCH.AND.Q.GT.MSQ) THEN
         SWI1=1.
        ELSE
         SWI1=0.
        ENDIF              
        RMIX(1)=RMIX(1)+SWI1*CC213(13)
        RMIX(2)=RMIX(2)+SWI1*CC213(14)
        IF(MSQU(3).GT.MT) THEN
         SMIX(1)=0.
         SMIX(2)=0.
         SMIX(3)=0.
         SMIX(4)=0.
         SMIX(5)=0.
         SMIX(6)=0.
        ENDIF
        IF(Q.LT.MCHA(1).AND.Q.GT.MSQU(1)) THEN
         SWI2=1.
        ELSE
         SWI2=0.
        ENDIF                
        IF(Q.LT.MCHA(1).AND.Q.GT.MSQU(2)) THEN
         SWI3=1.
        ELSE
         SWI3=0.
        ENDIF                
        ABMIX(1)=SWI2*CC211(13)+SWI3*CC212(13)
        ABMIX(2)=SWI2*CC211(14)+SWI3*CC212(14)
C 
        IF(MSQU(3).GT.(MW+0.5)) THEN
         QMAX=MSQ
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MCH.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(17,DT,T,CC213,GAMMAC2,WSC213)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LT.QMAX) THEN
         CC213(1)=CC213(1)+DSC213(1)        
         CC213(2)=CC213(2)+DSC213(2)
         CC213(3)=CC213(3)+DSC213(3)
         CC213(4)=CC213(4)+DSC213(4)
         CC213(5)=CC213(5)+DSC213(5)
         CC213(6)=CC213(6)+DSC213(6)
         CC213(7)=CC213(7)+DSC213(7)
         CC213(8)=CC213(8)+DSC213(8)
         CC213(9)=CC213(9)+DSC213(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC213,GAMMASM,WSC213)      
        ENDIF                  
         C2C13EFF=-1./2.d0*CC213(1)+CC213(2)
     $            -1./2.d0*CC213(4)-1./4.d0*CC213(5)
     $            +1./4.d0*CC213(7)       
         C3C13EFF=-1./2.d0*CC213(1)+CC213(3)
     $            -1./2.d0*CC213(4)-1./4.d0*CC213(5)
     $            +1./4.d0*CC213(7)-1./4.d0*CC213(9)       
       ENDIF
C
C ...  chargino(2)-squark(1) contribution  
C
       MSQ=MSQU(1)
       MCH=MCHA(2)    
       IF(Q.GT.MSQ.AND.QNEW.LE.MSQ) THEN
        G3SS=G(3)
        X=(MW/MSQU(1))**2
        Y=MCHA(2)/MBQ
        CC121(10)= GK(2,1)*HK(2,1)*16.*PI**2*X*Y/G3SS**2
        CC121(11)= GK(2,1)*HK(2,1)*16.*PI**2*X*Y/G3SS**2
        CC121(12)=-GK(2,1)*GK(2,1)*16.*PI**2*X/G3SS**2
       ENDIF    
C
       IF(MSQU(1).GT.MCHA(2)) THEN
        IF(MCHA(2).GT.(MW+0.5)) THEN
         QMAX=MCH
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MSQ.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(12,DT,T,CC121,GAMMAC1,WSC121)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LE.QMAX) THEN
         CC121(1)=CC121(1)+DSC121(1)        
         CC121(2)=CC121(2)+DSC121(2)
         CC121(3)=CC121(3)+DSC121(3)
         CC121(4)=CC121(4)+DSC121(4)
         CC121(5)=CC121(5)+DSC121(5)
         CC121(6)=CC121(6)+DSC121(6)
         CC121(7)=CC121(7)+DSC121(7)
         CC121(8)=CC121(8)+DSC121(8)
         CC121(9)=CC121(9)+DSC121(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC121,GAMMASM,WSC121)      
        ENDIF
         C2C21EFF=-1./2.d0*CC121(1)+CC121(2)
     $            -1./2.d0*CC121(4)-1./4.d0*CC121(5)
     $            +1./4.d0*CC121(7)       
         C3C21EFF=-1./2.d0*CC121(1)+CC121(3)
     $            -1./2.d0*CC121(4)-1./4.d0*CC121(5)
     $            +1./4.d0*CC121(7)-1./4.d0*CC121(9)       
       ELSE 
        IF(Q.LT.MCH.AND.Q.GT.MSQ) THEN
         SWI1=1.
        ELSE
         SWI1=0.
        ENDIF              
        RMIX(1)=RMIX(1)+SWI1*CC221(13)
        RMIX(2)=RMIX(2)+SWI1*CC221(14)
        IF(MSQU(1).GT.MT) THEN
         SMIX(1)=0.
         SMIX(2)=0.
         SMIX(3)=0.
         SMIX(4)=0.
         SMIX(5)=0.
         SMIX(6)=0.
        ENDIF
        IF(Q.LT.MCHA(2).AND.Q.GT.MSQU(2)) THEN
         SWI2=1.
        ELSE
         SWI2=0.
        ENDIF                
        IF(Q.LT.MCHA(2).AND.Q.GT.MSQU(3)) THEN
         SWI3=1.
        ELSE
         SWI3=0.
        ENDIF                    
        ABMIX(1)=SWI2*CC222(13)+SWI3*CC223(13)
        ABMIX(2)=SWI2*CC222(14)+SWI3*CC223(14)
C 
        IF(MSQU(1).GT.MW) THEN
         QMAX=MSQ
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MCH.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(17,DT,T,CC221,GAMMAC2,WSC221)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LE.QMAX) THEN
         CC221(1)=CC221(1)+DSC221(1)        
         CC221(2)=CC221(2)+DSC221(2)
         CC221(3)=CC221(3)+DSC221(3)
         CC221(4)=CC221(4)+DSC221(4)
         CC221(5)=CC221(5)+DSC221(5)
         CC221(6)=CC221(6)+DSC221(6)
         CC221(7)=CC221(7)+DSC221(7)
         CC221(8)=CC221(8)+DSC221(8)
         CC221(9)=CC221(9)+DSC221(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC221,GAMMASM,WSC221)      
        ENDIF                  
         C2C21EFF=-1./2.d0*CC221(1)+CC221(2)
     $            -1./2.d0*CC221(4)-1./4.d0*CC221(5)
     $            +1./4.d0*CC221(7)      
         C3C21EFF=-1./2.d0*CC221(1)+CC221(3)
     $            -1./2.d0*CC221(4)-1./4.d0*CC221(5)
     $            +1./4.d0*CC221(7)-1./4.d0*CC221(9)       
       ENDIF
C
C ...  chargino(2)-squark(2) contribution  
C
       MSQ=MSQU(2)
       MCH=MCHA(2)    
       IF(Q.GT.MSQ.AND.QNEW.LE.MSQ) THEN
        G3SS=G(3)
        X=(MW/MSQU(2))**2
        Y=MCHA(2)/MBQ
        CC122(10)= GK(2,2)*HK(2,2)*16.*PI**2*X*Y/G3SS**2
        CC122(11)= GK(2,2)*HK(2,2)*16.*PI**2*X*Y/G3SS**2
        CC122(12)=-GK(2,2)*GK(2,2)*16.*PI**2*X/G3SS**2
       ENDIF    
C
       IF(MSQU(2).GT.MCHA(2)) THEN
        IF(MCHA(2).GT.(MW+0.5)) THEN
         QMAX=MCH
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MSQ.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(12,DT,T,CC122,GAMMAC1,WSC122)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LE.QMAX) THEN
         CC122(1)=CC122(1)+DSC122(1)        
         CC122(2)=CC122(2)+DSC122(2)
         CC122(3)=CC122(3)+DSC122(3)
         CC122(4)=CC122(4)+DSC122(4)
         CC122(5)=CC122(5)+DSC122(5)
         CC122(6)=CC122(6)+DSC122(6)
         CC122(7)=CC122(7)+DSC122(7)
         CC122(8)=CC122(8)+DSC122(8)
         CC122(9)=CC122(9)+DSC122(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC122,GAMMASM,WSC122)      
C
        ENDIF
         C2C22EFF=-1./2.d0*CC122(1)+CC122(2)
     $            -1./2.d0*CC122(4)-1./4.d0*CC122(5)
     $            +1./4.d0*CC122(7)           
         C3C22EFF=-1./2.d0*CC122(1)+CC122(3)
     $            -1./2.d0*CC122(4)-1./4.d0*CC122(5)
     $            +1./4.d0*CC122(7)-1./4.d0*CC122(9)           
       ELSE 
        IF(Q.LT.MCH.AND.Q.GT.MSQ) THEN
         SWI1=1.
        ELSE
         SWI1=0.
        ENDIF              
        RMIX(1)=RMIX(1)+SWI1*CC222(13)
        RMIX(2)=RMIX(2)+SWI1*CC222(14)
        IF(MSQU(2).GT.MT) THEN
         SMIX(1)=0.
         SMIX(2)=0.
         SMIX(3)=0.
         SMIX(4)=0.
         SMIX(5)=0.
         SMIX(6)=0.
        ENDIF
        IF(Q.LT.MCHA(2).AND.Q.GT.MSQU(1)) THEN
         SWI2=1.
        ELSE
         SWI2=0.
        ENDIF                
        IF(Q.LT.MCHA(2).AND.Q.GT.MSQU(3)) THEN
         SWI3=1.
        ELSE
         SWI3=0.
        ENDIF                    
        ABMIX(1)=SWI2*CC221(13)+SWI3*CC223(13)
        ABMIX(2)=SWI2*CC221(14)+SWI3*CC223(14)
C 
        IF(MSQU(2).GT.(MW+0.5)) THEN
         QMAX=MSQ
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MCH.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(17,DT,T,CC222,GAMMAC2,WSC222)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LE.QMAX) THEN
         CC222(1)=CC222(1)+DSC222(1)        
         CC222(2)=CC222(2)+DSC222(2)
         CC222(3)=CC222(3)+DSC222(3)
         CC222(4)=CC222(4)+DSC222(4)
         CC222(5)=CC222(5)+DSC222(5)
         CC222(6)=CC222(6)+DSC222(6)
         CC222(7)=CC222(7)+DSC222(7)
         CC222(8)=CC222(8)+DSC222(8)
         CC222(9)=CC222(9)+DSC222(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC222,GAMMASM,WSC222)      
        ENDIF                  
         C2C22EFF=-1./2.d0*CC222(1)+CC222(2)
     $            -1./2.d0*CC222(4)-1./4.d0*CC222(5)
     $            +1./4.d0*CC222(7)       
         C3C22EFF=-1./2.d0*CC222(1)+CC222(3)
     $            -1./2.d0*CC222(4)-1./4.d0*CC222(5)
     $            +1./4.d0*CC222(7)-1./4.d0*CC222(9)       
       ENDIF
C
C ...  chargino(2)-squark(3) contribution  
C
       MSQ=MSQU(3)
       MCH=MCHA(2)   
       IF(Q.GT.MSQ.AND.QNEW.LE.MSQ) THEN
        G3SS=G(3)
        X=(MW/MSQU(3))**2
        Y=MCHA(2)/MBQ
        CC123(10)= HK(2,3)*16.*PI**2*X*Y/G3SS**2
        CC123(11)= HK(2,3)*16.*PI**2*X*Y/G3SS**2
        CC123(12)=-GK(2,3)*16.*PI**2*X/G3SS**2
       ENDIF    
C
       IF(MSQU(3).GT.MCHA(2)) THEN
        IF(MCHA(2).GT.(MW+0.5)) THEN
         QMAX=MCH
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MSQ.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(12,DT,T,CC123,GAMMAC1,WSC123)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LE.QMAX) THEN
         CC123(1)=CC123(1)+DSC123(1)        
         CC123(2)=CC123(2)+DSC123(2)
         CC123(3)=CC123(3)+DSC123(3)
         CC123(4)=CC123(4)+DSC123(4)
         CC123(5)=CC123(5)+DSC123(5)
         CC123(6)=CC123(6)+DSC123(6)
         CC123(7)=CC123(7)+DSC123(7)
         CC123(8)=CC123(8)+DSC123(8)
         CC123(9)=CC123(9)+DSC123(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC123,GAMMASM,WSC123)      
        ENDIF
         C2C23EFF=-1./2.d0*CC123(1)+CC123(2)
     $            -1./2.d0*CC123(4)-1./4.d0*CC123(5)
     $            +1./4.d0*CC123(7)       
         C3C23EFF=-1./2.d0*CC123(1)+CC123(3)
     $            -1./2.d0*CC123(4)-1./4.d0*CC123(5)
     $            +1./4.d0*CC123(7)-1./4.d0*CC123(9)       
       ELSE 
        IF(Q.LT.MCH.AND.Q.GT.MSQ) THEN
         SWI1=1.
        ELSE
         SWI1=0.
        ENDIF              
        RMIX(1)=RMIX(1)+SWI1*CC223(13)
        RMIX(2)=RMIX(2)+SWI1*CC223(14)
        IF(MSQU(3).GT.MT) THEN
         SMIX(1)=0.
         SMIX(2)=0.
         SMIX(3)=0.
         SMIX(4)=0.
         SMIX(5)=0.
         SMIX(6)=0.
        ENDIF
        IF(Q.LT.MCHA(2).AND.Q.GT.MSQU(1)) THEN
         SWI2=1.
        ELSE
         SWI2=0.
        ENDIF                
        IF(Q.LT.MCHA(2).AND.Q.GT.MSQU(2)) THEN
         SWI3=1.
        ELSE
         SWI3=0.
        ENDIF                    
        ABMIX(1)=SWI2*CC221(13)+SWI3*CC222(13)
        ABMIX(2)=SWI2*CC221(14)+SWI3*CC222(14)
C 
        IF(MSQU(3).GT.(MW+0.5)) THEN
         QMAX=MSQ
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MCH.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(17,DT,T,CC223,GAMMAC2,WSC223)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LE.QMAX) THEN
         CC223(1)=CC223(1)+DSC223(1)        
         CC223(2)=CC223(2)+DSC223(2)
         CC223(3)=CC223(3)+DSC223(3)
         CC223(4)=CC223(4)+DSC223(4)
         CC223(5)=CC223(5)+DSC223(5)
         CC223(6)=CC223(6)+DSC223(6)
         CC223(7)=CC223(7)+DSC223(7)
         CC223(8)=CC223(8)+DSC223(8)
         CC223(9)=CC223(9)+DSC223(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC223,GAMMASM,WSC223)      
        ENDIF                  
         C2C23EFF=-1./2.d0*CC223(1)+CC223(2)
     $            -1./2.d0*CC223(4)-1./4.d0*CC223(5)
     $            +1./4.d0*CC223(7)       
         C3C23EFF=-1./2.d0*CC223(1)+CC223(3)
     $            -1./2.d0*CC223(4)-1./4.d0*CC223(5)
     $            +1./4.d0*CC223(7)-1./4.d0*CC223(9)       
       ENDIF
c
c ... gluino and neutralinos contributions ...
c      GLUNENO must be recalled at MSUSY to set squark mass matrices
c
       IF(Q.GT.MW.AND.QNEW.LE.MW) THEN  

       CER=4./3.d0
       CEG=3.d0
       ALFA2=ALFAEM/SN2THW
       RALPH=ALES/ALFA2
       FACG7= 6.d0*(-1./3.d0)*RALPH*CER/(CKM(3,2)*CKM(3,3))
       FACG8=-RALPH/(CKM(3,2)*CKM(3,3))
       FACN7= 3.d0*(-1./3.d0)/(CKM(3,2)*CKM(3,3))
       FACN8=-1.d0/(CKM(3,2)*CKM(3,3))
C
C ...  gluino-dsquark1 ...
C              
       X=MGLU**2/MDSB(1)**2
       CALL FUNS(X,FN1,FN2,FN3,FN4) 
       C3GS1EFF=FACG7*(MW/MDSB(1))**2
     &          *(GAM(3,1)*GAM(2,1)*FN2
     $           -GAM(6,1)*GAM(2,1)*MGLU/MBQ*FN4)
       C2GS1EFF=FACG8*(MW/MDSB(1))**2
     &          *(GAM(3,1)*GAM(2,1)*(-CEG*FN1+(2*CER-CEG)*FN2)
     $           -GAM(6,1)*GAM(2,1)*MGLU/MBQ*(-CEG*FN3+(2*CER-CEG)*FN4))
C
C ...  gluino-dsquark2 ...
C             
       X=MGLU**2/MDSB(2)**2
       CALL FUNS(X,FN1,FN2,FN3,FN4) 
       C3GS2EFF=FACG7*(MW/MDSB(2))**2
     &          *(GAM(3,2)*GAM(2,2)*FN2
     $           -GAM(6,2)*GAM(2,2)*MGLU/MBQ*FN4)
       C2GS2EFF=FACG8*(MW/MDSB(2))**2
     &          *(GAM(3,2)*GAM(2,2)*(-CEG*FN1+(2*CER-CEG)*FN2)
     $           -GAM(6,2)*GAM(2,2)*MGLU/MBQ*(-CEG*FN3+(2*CER-CEG)*FN4))
C
C ...  gluino-dsquark3 ...
C       
       X=MGLU**2/MDSB(3)**2
       CALL FUNS(X,FN1,FN2,FN3,FN4) 
       C3GS3EFF=FACG7*(MW/MDSB(3))**2
     &          *(GAM(3,3)*GAM(2,3)*FN2
     $           -GAM(6,3)*GAM(2,3)*MGLU/MBQ*FN4)
       C2GS3EFF=FACG8*(MW/MDSB(3))**2
     &          *(GAM(3,3)*GAM(2,3)*(-CEG*FN1+(2*CER-CEG)*FN2)
     $           -GAM(6,3)*GAM(2,3)*MGLU/MBQ*(-CEG*FN3+(2*CER-CEG)*FN4))
C
C ...  gluino-dsquark4 ...
C       
       X=MGLU**2/MDSB(4)**2
       CALL FUNS(X,FN1,FN2,FN3,FN4) 
       C3GS4EFF=FACG7*(MW/MDSB(4))**2
     &          *(GAM(3,4)*GAM(2,4)*FN2
     $           -GAM(6,4)*GAM(2,4)*MGLU/MBQ*FN4)
       C2GS4EFF=FACG8*(MW/MDSB(4))**2
     &          *(GAM(3,4)*GAM(2,4)*(-CEG*FN1+(2*CER-CEG)*FN2)
     $           -GAM(6,4)*GAM(2,4)*MGLU/MBQ*(-CEG*FN3+(2*CER-CEG)*FN4))
C
C ...  gluino-dsquark5 ...
C       
       X=MGLU**2/MDSB(5)**2
       CALL FUNS(X,FN1,FN2,FN3,FN4) 
       C3GS5EFF=FACG7*(MW/MDSB(5))**2
     &          *(GAM(3,5)*GAM(2,5)*FN2
     $           -GAM(6,5)*GAM(2,5)*MGLU/MBQ*FN4)
       C2GS5EFF=FACG8*(MW/MDSB(5))**2
     &          *(GAM(3,5)*GAM(2,5)*(-CEG*FN1+(2*CER-CEG)*FN2)
     $           -GAM(6,5)*GAM(2,5)*MGLU/MBQ*(-CEG*FN3+(2*CER-CEG)*FN4))
C
C ...  gluino-dsquark6 ...
C       
       X=MGLU**2/MDSB(6)**2
       CALL FUNS(X,FN1,FN2,FN3,FN4) 
       C3GS6EFF=FACG7*(MW/MDSB(6))**2
     &          *(GAM(3,6)*GAM(2,6)*FN2
     $           -GAM(6,6)*GAM(2,6)*MGLU/MBQ*FN4)
       C2GS6EFF=FACG8*(MW/MDSB(6))**2
     &          *(GAM(3,6)*GAM(2,6)*(-CEG*FN1+(2*CER-CEG)*FN2)
     $           -GAM(6,6)*GAM(2,6)*MGLU/MBQ*(-CEG*FN3+(2*CER-CEG)*FN4))
C
C ...  neutralinos-dsquark1 ...
C              
       C3NS1EFF=0.
       C2NS1EFF=0.
       DO K=1,4                    ! loop through neutralinos
         X=MCH0(K)**2/MDSB(1)**2
         CALL FUNS(X,FN1,FN2,FN3,FN4) 
         C3NS1EFF=C3NS1EFF
     $            +FACN7*(MW/MDSB(1))**2
     $             *(2.*AGE(K,3,1)*AGE(K,2,1)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,1)-AHE(K,3,1))
     $               *AGE(K,2,1)*MCH0(K)/MBQ*FN4)  
         C2NS1EFF=C2NS1EFF
     $            +FACN8*(MW/MDSB(1))**2
     $             *(2.*AGE(K,3,1)*AGE(K,2,1)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,1)-AHE(K,3,1))
     $               *AGE(K,2,1)*MCH0(K)/MBQ*FN4)  
       ENDDO
C
C ...  neutralinos-dsquark2 ...
C              
       C3NS2EFF=0.
       C2NS2EFF=0.
       DO K=1,4
         X=MCH0(K)**2/MDSB(2)**2
         CALL FUNS(X,FN1,FN2,FN3,FN4) 
         C3NS2EFF=C3NS2EFF
     $            +FACN7*(MW/MDSB(2))**2
     $             *(2.*AGE(K,3,2)*AGE(K,2,2)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,2)-AHE(K,3,2))
     $               *AGE(K,2,2)*MCH0(K)/MBQ*FN4)  
         C2NS2EFF=C2NS2EFF
     $            +FACN8*(MW/MDSB(2))**2
     $             *(2.*AGE(K,3,2)*AGE(K,2,2)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,2)-AHE(K,3,2))
     $               *AGE(K,2,2)*MCH0(K)/MBQ*FN4)  
       ENDDO
C
C ...  neutralinos-dsquark3 ...
C              
       C3NS3EFF=0.
       C2NS3EFF=0.
       DO K=1,4
         X=MCH0(K)**2/MDSB(3)**2
         CALL FUNS(X,FN1,FN2,FN3,FN4) 
         C3NS3EFF=C3NS3EFF
     $            +FACN7*(MW/MDSB(3))**2
     $             *(2.*AGE(K,3,3)*AGE(K,2,3)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,3)-AHE(K,3,3))
     $               *AGE(K,2,3)*MCH0(K)/MBQ*FN4)  
         C2NS3EFF=C2NS3EFF
     $            +FACN8*(MW/MDSB(3))**2
     $             *(2.*AGE(K,3,3)*AGE(K,2,3)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,3)-AHE(K,3,3))
     $               *AGE(K,2,3)*MCH0(K)/MBQ*FN4)  
       ENDDO
C
C ...  neutralinos-dsquark4 ...
C              
       C3NS4EFF=0.
       C2NS4EFF=0.
       DO K=1,4
         X=MCH0(K)**2/MDSB(4)**2
         CALL FUNS(X,FN1,FN2,FN3,FN4) 
         C3NS4EFF=C3NS4EFF
     $            +FACN7*(MW/MDSB(4))**2
     $             *(2.*AGE(K,3,4)*AGE(K,2,4)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,4)-AHE(K,3,4))
     $               *AGE(K,2,4)*MCH0(K)/MBQ*FN4)  
         C2NS4EFF=C2NS4EFF
     $            +FACN8*(MW/MDSB(4))**2
     $             *(2.*AGE(K,3,4)*AGE(K,2,4)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,4)-AHE(K,3,4))
     $               *AGE(K,2,4)*MCH0(K)/MBQ*FN4)  
       ENDDO
c
C ...  neutralinos-dsquark5 ...
C              
       C3NS5EFF=0.
       C2NS5EFF=0.
       DO K=1,4
         X=MCH0(K)**2/MDSB(5)**2
         CALL FUNS(X,FN1,FN2,FN3,FN4) 
         C3NS5EFF=C3NS5EFF
     $            +FACN7*(MW/MDSB(5))**2
     $             *(2.*AGE(K,3,5)*AGE(K,2,5)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,5)-AHE(K,3,5))
     $               *AGE(K,2,5)*MCH0(K)/MBQ*FN4)  
         C2NS5EFF=C2NS5EFF
     $            +FACN8*(MW/MDSB(5))**2
     $             *(2.*AGE(K,3,5)*AGE(K,2,5)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,5)-AHE(K,3,5))
     $               *AGE(K,2,5)*MCH0(K)/MBQ*FN4)  
       ENDDO
C
C ...  neutralinos-dsquark6 ...
C              
       C3NS6EFF=0.
       C2NS6EFF=0.
       DO K=1,4
         X=MCH0(K)**2/MDSB(6)**2
         CALL FUNS(X,FN1,FN2,FN3,FN4) 
         C3NS6EFF=C3NS6EFF
     $            +FACN7*(MW/MDSB(6))**2
     $             *(2.*AGE(K,3,6)*AGE(K,2,6)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,6)-AHE(K,3,6))
     $               *AGE(K,2,6)*MCH0(K)/MBQ*FN4)  
         C2NS6EFF=C2NS6EFF
     $            +FACN8*(MW/MDSB(6))**2
     $             *(2.*AGE(K,3,6)*AGE(K,2,6)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,6)-AHE(K,3,6))
     $               *AGE(K,2,6)*MCH0(K)/MBQ*FN4)  
       ENDDO
      ENDIF
C
      SINTW2I=1./SN2THW
      
      XTW=MT**2/MW**2
C
c ... here we switch from Anlauf's and Ciuchini's 
c     normalization to that of Greub:
C     Ci(Greub)=Ci(Ciuchini,Buras)                for i=1,...,6
c     C7(Greub)=-Qb*C7(Anlauf)=1/3*C7(Anlauf)
c     C8(Greub)=-C8(Anlauf)
c ...
c
      IF(Q.GT.MW.AND.QNEW.LE.MW) THEN 
       CI1(1)= ALES/4./PI*11./2.d0
       CI1(2)= 1.d0-ALES/4./PI*11./6.d0
       CI1(3)=-ALES/24.d0/PI*(EI(XTW)-2./3.d0)
     $        +ALEL/6.d0/PI*SINTW2I*(2.*BI(XTW)+CI(XTW))
       CI1(4)= ALES/8.d0/PI*(EI(XTW)-2./3.d0)
       CI1(5)=-ALES/24.d0/PI*(EI(XTW)-2./3.d0)
       CI1(6)= ALES/8.d0/PI*(EI(XTW)-2./3.d0)
       CI1(7)= 1./3.d0*(C3SMEFF+C3HPEFF
     $                  +C3C11EFF+C3C12EFF+C3C13EFF
     $                  +C3C21EFF+C3C22EFF+C3C23EFF
     $                  +C3GS1EFF+C3GS2EFF+C3GS3EFF
     $                  +C3GS4EFF+C3GS5EFF+C3GS6EFF
     $                  +C3NS1EFF+C3NS2EFF+C3NS3EFF
     $                  +C3NS4EFF+C3NS5EFF+C3NS6EFF)
       CI1(8)=-1.d0*(C2SMEFF+C2HPEFF
     $               +C2C11EFF+C2C12EFF+C2C13EFF
     $               +C2C21EFF+C2C22EFF+C2C23EFF
     $               +C2GS1EFF+C2GS2EFF+C2GS3EFF
     $               +C2GS4EFF+C2GS5EFF+C2GS6EFF
     $               +C2NS1EFF+C2NS2EFF+C2NS3EFF
     $               +C2NS4EFF+C2NS5EFF+C2NS6EFF)
c
       CI2(1)= ALES/4.d0/PI*11./2.d0
       CI2(2)= 1.d0-ALES/4.d0/PI*11./6.d0
       CI2(3)=-ALES/24.d0/PI*(EI(XTW)-2./3.d0)
     $        +ALEL/6.d0/PI*SINTW2I*(2.*BI(XTW)+CI(XTW))
       CI2(4)= ALES/8.d0/PI*(EI(XTW)-2./3.d0)
       CI2(5)=-ALES/24.d0/PI*(EI(XTW)-2./3.d0)
       CI2(6)= ALES/8.d0/PI*(EI(XTW)-2./3.d0)
       CI2(7)= 1./3.d0*(C3SMEFF+C3HPEFF
     $                  +C3C11EFF+C3C12EFF+C3C13EFF
     $                  +C3C21EFF+C3C22EFF+C3C23EFF
     $                  +C3GS1EFF+C3GS2EFF+C3GS3EFF
     $                  +C3GS4EFF+C3GS5EFF+C3GS6EFF
     $                  +C3NS1EFF+C3NS2EFF+C3NS3EFF
     $                  +C3NS4EFF+C3NS5EFF+C3NS6EFF)
       CI2(8)=-1.d0*(C2SMEFF+C2HPEFF
     $               +C2C11EFF+C2C12EFF+C2C13EFF
     $               +C2C21EFF+C2C22EFF+C2C23EFF
     $               +C2GS1EFF+C2GS2EFF+C2GS3EFF
     $               +C2GS4EFF+C2GS5EFF+C2GS6EFF
     $               +C2NS1EFF+C2NS2EFF+C2NS3EFF
     $               +C2NS4EFF+C2NS5EFF+C2NS6EFF)
c
        C3CHEFF=C3C11EFF+C3C12EFF+C3C13EFF
     $         +C3C21EFF+C3C22EFF+C3C23EFF
        C3GLEFF=C3GS1EFF+C3GS2EFF+C3GS3EFF
     $         +C3GS4EFF+C3GS5EFF+C3GS6EFF
        C3NEEFF=C3NS1EFF+C3NS2EFF+C3NS3EFF
     $         +C3NS4EFF+C3NS5EFF+C3NS6EFF
c
        C2CHEFF=C2C11EFF+C2C12EFF+C2C13EFF
     $         +C2C21EFF+C2C22EFF+C2C23EFF
        C2GLEFF=C2GS1EFF+C2GS2EFF+C2GS3EFF
     $         +C2GS4EFF+C2GS5EFF+C2GS6EFF
        C2NEEFF=C2NS1EFF+C2NS2EFF+C2NS3EFF
     $         +C2NS4EFF+C2NS5EFF+C2NS6EFF
      ENDIF
C
C           
      IF(Q.LT.MW) THEN
        T=DT*FLOAT(II-1)
        CALL DRKSTP(8,DT,T,CI1,GAMMAWB1,WCI1)
        T=DT*FLOAT(II-1)
        CALL DRKSTP(8,DT,T,CI2,GAMMAWB2,WCI2)
C
C------- Computation of b->s\gamma branching ratio --------------------
C
       IF(Q.GT.(2.d0*MB)) THEN
         XQ2 = SNGL(Q**2)
	 XMT = SNGL(MT)
         ASMB= DBLE(SUALFS(XQ2,.36,XMT,3))
         AEMB= DBLE(SUALFE(XQ2))
         ALFA2=AEMB/SN2THW
         G2=SQRT(4.d0*PI*ALFA2)
         GEF=SQRT(2.d0)/8.d0*G2**2/MW**2      
       ENDIF
       BQLOG=LOG(MB/Q)
       ZI=(MC/MB)**2
       LOGZI=LOG(ZI)
       RER2=2./243.d0
     &      *(-833.d0+144.*PI**2*ZI**(3./2.)
     $        +(1728.d0-180.*PI**2-1296.d0*1.20206
     $          +(1296.d0-324.*PI**2)*LOGZI+108.d0*LOGZI**2
     $          +36.d0*LOGZI**3)*ZI
     $        +(648.d0+72.d0*PI**2+(432.d0-216.d0*PI**2)*LOGZI
     $          +36.d0*LOGZI**3)*ZI**2
     $        +(-54.d0-84.d0*PI**2+1092.d0*LOGZI
     $          -756.d0*LOGZI**2)*ZI**3)
       IMR2=16.*PI/81.d0
     &      *(-5.d0+(45.d0-3.d0*PI**2+9.*LOGZI+9.*LOGZI**2)*ZI
     $        +(-3.d0*PI**2+9.d0*LOGZI**2)*ZI**2
     $        +(28.d0-12.d0*LOGZI)*ZI**3)
       RER7= 8./9.d0*(4.d0-PI**2)
       IMR7= 0.d0
       RER8=-4./27.d0*(-33.d0+2.*PI**2)
       IMR8=-4./27.d0*(-6.d0*PI)
       REDE=ASMB/4.d0/PI
     $      *((416./81.d0*BQLOG+RER2)*CI1(2)
     $        +(32./3.d0*BQLOG+RER7)*CI1(7)
     $        -(32./9.d0*BQLOG+RER8)*CI1(8))
       IMDE=ASMB/4.d0/PI
     $      *(IMR2*CI1(2)+IMR7*CI1(7)+IMR8*CI1(8))
       EF=(1.d0-8./3.d0*ASMB/PI)
c
C...Compute bremsstrahlung corrections
c
       CALL BREMS(AEMB,ASMB,MS,MC,MB,Q,CI1,GBREMS)
c
c...Implement virtual QCD corrections according to Eq(5.6) of Greub.
c
       GVIRT=AEMB/32.d0/PI**4*EF*(CI2(7)**2+2.d0*CI2(7)*REDE)
c     
       GAMMASG=GVIRT+GBREMS
       gammasg0=AEMB/32.d0/PI**4*EF*CI1(7)**2
c
c...Evaluate semileptonic decay width as in Eq(5.9) of Greub.
c
       GAMMASL=1./192.d0/PI**3*GES(MC/MB)
     $         *(1.d0-2./3.d0/PI*ASMB*FES(MC/MB))
c
c...Compute b->s\gamma branching ratio using Eq(5.8) of Greub.
c
       BRSG = COBA*GAMMASG/GAMMASL*BRSL
       BRSG0= COBA*GAMMASG0/GAMMASL*BRSL
c
C...Evaluate theoretical uncertainties from scale variation
c
       if(q.Gt.(2.*MB).AND.QNEW.LT.(2.*MB)) then
         brup=BRSG*1.d+4
       ENDIF
       if(q.Gt.(1./2.*MB).AND.QNEW.LT.(1./2.*MB)) then
         brdo=BRSG*1.d+4
       ENDIF
c
C...Computation of the final result
c
       IF(q.Gt.MB.AND.QNEW.LT.MB) THEN
         BFBSG =BRSG*1.d+4
         BFBSG0=BRSG0*1.d+4
c       write(6,*) 'q=',q,ci1(7),ci2(7),rede,gvirt,gbrems
c       write(6,*) 'ci=',ci1(2),ci1(7),ci1(8),rer2,rer7,rer8,asmb
c       write(6,*) 'brsg,brsg0=',brsg,brsg0
        ENDIF
      ENDIF
c---------------------------------------------------------------
400   CONTINUE 
c
        C3C1NAI=C3C11NAI+C3C12NAI+C3C13NAI
        C3C2NAI=C3C21NAI+C3C22NAI+C3C23NAI
        C3CHNAI=C3C11NAI+C3C12NAI+C3C13NAI
     $         +C3C21NAI+C3C22NAI+C3C23NAI
        C3GLNAI=C3GS1EFF+C3GS2EFF+C3GS3EFF
     $         +C3GS4EFF+C3GS5EFF+C3GS6EFF
        C3NENAI=C3NS1EFF+C3NS2EFF+C3NS3EFF
     $         +C3NS4EFF+C3NS5EFF+C3NS6EFF
C
        C3C1EFF=C3C11EFF+C3C12EFF+C3C13EFF
        C3C2EFF=C3C21EFF+C3C22EFF+C3C23EFF
        C3CHEFF=C3C11EFF+C3C12EFF+C3C13EFF
     $         +C3C21EFF+C3C22EFF+C3C23EFF
        C3GLEFF=C3GS1EFF+C3GS2EFF+C3GS3EFF
     $         +C3GS4EFF+C3GS5EFF+C3GS6EFF
        C3NEEFF=C3NS1EFF+C3NS2EFF+C3NS3EFF
     $         +C3NS4EFF+C3NS5EFF+C3NS6EFF
c
      IF (IWR.EQ.1) THEN
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'SM:          '
        WRITE(6,*)'QCD corrected=',c3smeff,' naive=',C3SMNAI                
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'higgs:       '
        WRITE(6,*)'QCD corrected=',c3hpeff,' naive=',C3HPNAI      
        WRITE(6,*)'charged Higgs mass=',mhplus
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'chargino(1)-stop(1): '
        WRITE(6,*)'QCD corrected=',c3c11eff,' naive=',C3C11NAI       
        WRITE(6,*)'chargino mass=',mcha(1),' squark mass=',msqu(1)  
        WRITE(6,*)'G=',gk(1,1),' H=',hk(1,1)
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'chargino(1)-stop(2): '
        WRITE(6,*)'QCD corrected=',c3c12eff,' naive=',C3C12NAI        
        WRITE(6,*)'chargino mass=',mcha(1),' squark mass=',msqu(2)  
        WRITE(6,*)'G=',gk(1,2),' H=',hk(1,2)
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'chargino(1)-up+charm squarks: '
        WRITE(6,*)'QCD corrected=',c3c13eff,' naive=',C3C13NAI       
        WRITE(6,*)'chargino mass=',mcha(1),' squark mass=',msqu(3)  
        WRITE(6,*)'G=',gk(1,3),' H=',hk(1,3)
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'chargino(2)-stop(1): '
        WRITE(6,*)'QCD corrected=',c3c21eff,' naive=',C3C21NAI       
        WRITE(6,*)'chargino mass=',mcha(2),' squark mass=',msqu(1)  
        WRITE(6,*)'G=',gk(2,1),' H=',hk(2,1)
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'chargino(2)-stop(2): '
        WRITE(6,*)'QCD corrected=',c3c22eff,' naive=',C3C22NAI        
        WRITE(6,*)'chargino mass=',mcha(2),' squark mass=',msqu(2)  
        WRITE(6,*)'G=',gk(2,2),' H=',hk(2,2)
	WRITE(6,*)'____________________________________________'
	WRITE(6,*)'chargino(2)-up+charm squarks: '
        WRITE(6,*)'QCD corrected=',c3c23eff,' naive=',C3C23NAI       
        WRITE(6,*)'chargino mass=',mcha(2),' squark mass=',msqu(3)  
        WRITE(6,*)'G=',gk(2,3),' H=',hk(2,3)
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'sum of chargino-squarks contributions: '
        WRITE(6,*)'chargino(1)'
        WRITE(6,*)'QCD corrected=',c3c1eff,' naive=',C3C1NAI   
        WRITE(6,*)'chargino(2)'
        WRITE(6,*)'QCD corrected=',c3c2eff,' naive=',C3C2NAI    
        WRITE(6,*)'total'
        WRITE(6,*)'QCD corrected=',c3cheff,' naive=',C3ChNAI    
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'sum of gluino-squarks contributions: '
        WRITE(6,*)'QCD corrected=',1./3.*c3gleff,' naive=',
     $              1./3.*C3glNAI    
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'sum of neutralino-squarks contributions: '
        WRITE(6,*)'QCD corrected=',c3neeff,' naive=',C3neNAI    
        WRITE(6,*)'____________________________________________'
        write(6,*) 'BFBSG,BFBSG0=',BFBSG,BFBSG0,'  x 10^-4'
      END IF
C--------------------------------------------------------------
1000  continue              
      
      RETURN
      END
