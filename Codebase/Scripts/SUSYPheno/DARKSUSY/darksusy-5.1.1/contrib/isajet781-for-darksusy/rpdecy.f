CDECK  ID>, RPDECY
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
      SUBROUTINE RPDECY
C----------------------------------------------------------------------- 
C     SUBROUTINE TO CALCULATE ALL THE R-PARITY VIOLATING RATES
C     2-BODY SQUARK AND SLEPTON
C     3-BODY NEUTRALINO, GLUINO AND CHARGINO
C----------------------------------------------------------------------- 
C
      IMPLICIT NONE
C
C          ISAJET common blocks and EQUIVALENCE
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
      REAL MSS(72)
      EQUIVALENCE (AMGLSS,MSS(1))
C
C--Common block containing the Rparity violating couplings
C
      COMMON/RSLASH/LAMDA1(3,3,3),LAMDA2(3,3,3),LAMDA3(3,3,3),RPARTY
      LOGICAL RPARTY
      REAL LAMDA1,LAMDA2,LAMDA3
      SAVE /RSLASH/
C   
C--Common block containing the Standard Model parameters          
C 
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
C
C--Common Block to contain R-parity violating decay rates
C
      INTEGER NRPMD
      PARAMETER (NRPMD=5000)
      COMMON/RPARRT/NSSMD2,ISSMD2(NRPMD),JSSMD2(5,NRPMD),
     &              GSSMD2(NRPMD),BSSMD2(NRPMD)
      REAL GSSMD2,BSSMD2
      INTEGER ISSMD2,JSSMD2,NSSMD2
      SAVE /RPARRT/
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
      REAL M(4),RESM(6),WIDTH(6),A(6),B(6)
      LOGICAL CRSTRM(6)
C--External isajet gaussian quadrature routine
      REAL RPINF1,RPINF2
      EXTERNAL RPRATE,RPINF1,RPINF2
C--Local Variables
      REAL MQU(6),MLP(3),MSLLT(6),MSLRT(3),MSQLT(6),MSQRT(6),MIXING(6),
     &     RATE,RPRATE,SLLTWD(6),SLRTWD(3),SQLTWD(6),CHARWD,CHARM(2), 
     &     SQRTWD(6),PI,G,ECHAR,CBETA,SBETA,CWEAK,LAMCOL,NEUTWD,  
     &     SWEAK,NPRIME(4,4),TBETA,GLUWD,WMXVSS(2,2),WMXUSS(2,2),
     &     BMIXSS(2,2),TMIXSS(2,2),LMIXSS(2,2),LCHAR(4),RCHAR(4),
     &     MCHAR(2),THX,THY
      INTEGER I,J,K,L,PIN,POUT1,POUT2,POUT3,CHANEL,N,NMSIGN(4),MIX,
     &        CMSIGN(2)
      REAL ZERO,ONE,EPS
      PARAMETER (ZERO=0.0,ONE=1.0,EPS=1E-40)
C--Couplings, etc
      PI = 3.1415926E0
      ECHAR = SQRT(4*PI*ALFAEM) 
      G     = SQRT(4*PI*ALFA2)
      TBETA = 1./RV2V1
      CBETA = COS(ATAN(TBETA))
      SBETA = SIN(ATAN(TBETA))
      CWEAK = COS(ASIN(SQRT(SN2THW)))
      SWEAK = SQRT(SN2THW)
C--Neutralino mixing
      DO I=1,4
        NMSIGN(I) = -INT(AMZISS(I)/ABS(AMZISS(I)))
      ENDDO
      DO I=1,4
        NPRIME(I,1) = ZMIXSS(4,I)*CWEAK+ZMIXSS(3,I)*SWEAK
        NPRIME(I,2) = -ZMIXSS(4,I)*SWEAK+ZMIXSS(3,I)*CWEAK
        NPRIME(I,3) = -ZMIXSS(2,I)
        NPRIME(I,4) = -ZMIXSS(1,I)
      ENDDO
C--Chargino mixing, fixed 01/04/00 PR
      THX=SIGN(1.0D0,1.0D0/TAN(GAMMAL))
      THY=SIGN(1.0D0,1.0D0/TAN(GAMMAR))
      CMSIGN(1) = -INT(AMW1SS/ABS(AMW1SS))
      CMSIGN(2) = -INT(AMW2SS/ABS(AMW2SS))
      WMXVSS(1,1) =    -SIN(GAMMAR)
      WMXVSS(1,2) =    -COS(GAMMAR)
      WMXVSS(2,1) =-THY*COS(GAMMAR)
      WMXVSS(2,2) = THY*SIN(GAMMAR) 
      WMXUSS(1,1) =    -SIN(GAMMAL)
      WMXUSS(1,2) =    -COS(GAMMAL)
      WMXUSS(2,1) =-THX*COS(GAMMAL)
      WMXUSS(2,2) = THX*SIN(GAMMAL)
      CHARM(1) = AMW1SS
      CHARM(2) = AMW2SS
C--Number of R-parity violating modes
      NSSMD2    = 0
C--Set up local mass variables 
      MQU(1)   = AMDN 
      MQU(2)   = AMUP
      MQU(3)   = AMST
      MQU(4)   = AMCH
      MQU(5)   = AMBT
      MQU(6)   = AMTP
      MLP(1)   = AME 
      MLP(2)   = AMMU
      MLP(3)   = AMTAU
      MSLLT(1) = AMELSS
      MSLLT(2) = AMN1SS
      MSLLT(3) = AMMLSS
      MSLLT(4) = AMN2SS
      MSLLT(5) = AML1SS
      MSLLT(6) = AMN3SS
      MSLRT(1) = AMERSS
      MSLRT(2) = AMMRSS
      MSLRT(3) = AML2SS
      MSQLT(1) = AMDLSS
      MSQLT(2) = AMULSS
      MSQLT(3) = AMSLSS
      MSQLT(4) = AMCLSS
      MSQLT(5) = AMB1SS
      MSQLT(6) = AMT1SS
      MSQRT(1) = AMDRSS
      MSQRT(2) = AMURSS
      MSQRT(3) = AMSRSS
      MSQRT(4) = AMCRSS
      MSQRT(5) = AMB2SS
      MSQRT(6) = AMT2SS
C--Scalar top/bottom/tau mixing, bug fix 01/04/00 PR
      TMIXSS(1,1) =  COS(THETAT) 
      TMIXSS(1,2) = -SIN(THETAT)
      TMIXSS(2,1) =  SIN(THETAT)
      TMIXSS(2,2) =  COS(THETAT)
      BMIXSS(1,1) =  COS(THETAB)
      BMIXSS(1,2) = -SIN(THETAB)
      BMIXSS(2,1) =  SIN(THETAB)
      BMIXSS(2,2) =  COS(THETAB)
      LMIXSS(1,1) =  COS(THETAL)
      LMIXSS(1,2) = -SIN(THETAL)
      LMIXSS(2,1) =  SIN(THETAL)
      LMIXSS(2,2) =  COS(THETAL)
C--Now begin the rate calculation
C  Scalar decay rates via Rparity violation
C--CHECKED AGAINST CALCULATIONS PR 1/4/99
C--First the rates of left charged leptons
      DO I=1,3
        DO J=1,3
          DO K=1,3
C--Via LLE
            RATE = ZERO
            IF(I.NE.3) THEN 
              IF(ABS(LAMDA1(J,I,K)).GT.EPS) 
     &          RATE = RPRATE(LAMDA1(J,I,K)**2,
     &                      MSLLT(2*I-1),MLP(K),ZERO)
                CALL RPMODA(RATE,30+2*I,10+2*K,-(9+2*J),0)
            ELSE
C--New for left/right stau mixing
              IF(ABS(LAMDA1(J,I,K)).GT.EPS) THEN
                RATE = RPRATE(LAMDA1(J,I,K)**2*LMIXSS(1,1)**2,
     &                      MSLLT(2*I-1),MLP(K),ZERO) 
                CALL RPMODA(RATE,36,10+2*K,-(9+2*J),0)
                RATE = RPRATE(LAMDA1(J,I,K)**2*LMIXSS(1,2)**2,
     &                      MSLRT(I),MLP(K),ZERO)
                CALL RPMODA(RATE,56,10+2*K,-(9+2*J),0)
              ENDIF
            ENDIF
C--Via LUD
            RATE = ZERO
            IF(J.EQ.1) POUT1 = -1
            IF(J.GT.1) POUT1 = -2*J
            IF(K.EQ.1) POUT2 = 2
            IF(K.GT.1) POUT2 = 2*K-1
            IF(I.NE.3) THEN
              IF(ABS(LAMDA2(I,J,K)).GT.EPS) 
     &          RATE = RPRATE(3*LAMDA2(I,J,K)**2,MSLLT(2*I-1),
     &                        MQU(2*J),MQU(2*K-1))
                CALL RPMODA(RATE,30+2*I,POUT1,POUT2,0)
            ELSE
C--New for left/right stau mixing
              IF(ABS(LAMDA2(I,J,K)).GT.EPS) THEN
                RATE=RPRATE(3*LAMDA2(I,J,K)**2*LMIXSS(1,1)**2,
     &                      MSLLT(2*I-1), MQU(2*J),MQU(2*K-1))
                CALL RPMODA(RATE,30+2*I,POUT1,POUT2,0)
                RATE=RPRATE(3*LAMDA2(I,J,K)**2*LMIXSS(1,2)**2,
     &                      MSLRT(I), MQU(2*J),MQU(2*K-1))
                CALL RPMODA(RATE,56,POUT1,POUT2,0)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C--Now right charged leptons via LLE
      DO K=1,3
        DO I=1,3
          DO J=1,3
            RATE = ZERO
            IF(K.NE.3) THEN
              IF(ABS(LAMDA1(I,J,K)).GT.EPS) 
     &          RATE = RPRATE(LAMDA1(I,J,K)**2,MSLRT(K),
     &                     MLP(J),ZERO) 
                CALL RPMODA(RATE,50+2*K,10+2*J,9+2*I,0)
            ELSE
C--New for left/right stau mixing
              IF(ABS(LAMDA1(I,J,K)).GT.EPS) THEN
                RATE = RPRATE(LAMDA1(I,J,K)**2*LMIXSS(2,1)**2,
     &                     MSLLT(2*K-1),MLP(J),ZERO) 
                CALL RPMODA(RATE,36,10+2*J,9+2*I,0)
                RATE = RPRATE(LAMDA1(I,J,K)**2*LMIXSS(2,2)**2,
     &                     MSLRT(K),MLP(J),ZERO) 
                CALL RPMODA(RATE,56,10+2*J,9+2*I,0)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C--Now sneutrinos 
      DO I=1,3
        DO J=1,3
          DO K=1,3
C--Via LLE
            RATE = ZERO
            IF(ABS(LAMDA1(J,I,K)).GT.EPS) 
     &        RATE = RPRATE(LAMDA1(J,I,K)**2,MSLLT(2*I),
     &                      MLP(J),MLP(K))
              CALL RPMODA(RATE,29+2*I,-(10+2*J),10+2*K,0)
C--Via LUD
            RATE = ZERO
            IF(ABS(LAMDA2(I,J,K)).GT.EPS)
     &         RATE = RPRATE(3*LAMDA2(I,J,K)**2,MSLLT(2*I)
     &                      ,MQU(2*J-1),MQU(2*K-1))
              IF(J.EQ.1) POUT1 = -2
              IF(J.GT.1) POUT1 = 1-2*J
              IF(K.EQ.1) POUT2 = 2
              IF(K.GT.1) POUT2 = 2*K-1
              CALL RPMODA(RATE,29+2*I,POUT1,POUT2,0)
          ENDDO
        ENDDO
      ENDDO  
C--Now left up squarks via LUD
      DO J=1,3
        DO I=1,3
          DO K=1,3
            RATE = ZERO
            IF(J.EQ.1) PIN   = 21
            IF(J.GT.1) PIN   = 20+2*J
            IF(K.EQ.1) POUT1 = 2
            IF(K.GT.1) POUT1 = 2*K-1
            IF(J.NE.3) THEN
              IF(ABS(LAMDA2(I,J,K)).GT.EPS) 
     &           RATE = RPRATE(LAMDA2(I,J,K)**2,MSQLT(2*J),
     &                        MQU(2*K-1),MLP(I))
              CALL RPMODA(RATE,PIN,-(2*I+10),POUT1,0)
            ELSE
C-- New for left/right stop mixing
              IF(ABS(LAMDA2(I,J,K)).GT.EPS) THEN
                RATE = RPRATE(LAMDA2(I,J,K)**2*TMIXSS(1,1)**2,
     &                        MSQLT(2*J),MQU(2*K-1),MLP(I)) 
                CALL RPMODA(RATE,26,-(2*I+10),POUT1,0)
                RATE = RPRATE(LAMDA2(I,J,K)**2*TMIXSS(1,2)**2,
     &                        MSQRT(2*J),MQU(2*K-1),MLP(I)) 
                CALL RPMODA(RATE,46,-(2*I+10),POUT1,0)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO 
C--Now left down squarks via LUD
      DO J=1,3
        DO I=1,3
          DO K=1,3
            RATE = ZERO 
            IF(J.EQ.1) PIN   = 22
            IF(J.GT.1) PIN   = 19+2*J
            IF(K.EQ.1) POUT1 = 2
            IF(K.GT.1) POUT1 = 2*K-1
            IF(J.NE.3) THEN
              IF(ABS(LAMDA2(I,J,K)).GT.EPS) 
     &          RATE = RPRATE(LAMDA2(I,J,K)**2,MSQLT(2*J-1),
     &                        MQU(2*K-1),ZERO)
              CALL RPMODA(RATE,PIN,-9-2*I,POUT1,0)
            ELSE
C--New for left/right sbottom mixing
              IF(ABS(LAMDA2(I,J,K)).GT.EPS) THEN
                RATE = RPRATE(LAMDA2(I,J,K)**2*BMIXSS(1,1)**2,
     &                        MSQLT(2*J-1),MQU(2*K-1),ZERO)
                CALL RPMODA(RATE,25,-9-2*I,POUT1,0)
C--bug fix 9/04/02 by P.R.
                RATE = RPRATE(LAMDA2(I,J,K)**2*BMIXSS(2,1)**2,
     &                        MSQRT(2*J-1),MQU(2*K-1),ZERO)
                CALL RPMODA(RATE,45,-9-2*I,POUT1,0)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO 
C--Now right up squarks via UDD
      DO I=1,3
        DO J=1,3
          DO K=1,3
            RATE = ZERO 
            IF(I.EQ.1) PIN   = 41
            IF(I.GT.1) PIN   = 40+2*I
            IF(J.EQ.1) POUT1 = -2
            IF(J.GT.1) POUT1 = 1-2*J
            IF(K.EQ.1) POUT2 = -2
            IF(K.GT.1) POUT2 = 1-2*K
            IF(I.NE.3) THEN
              IF(ABS(LAMDA3(I,J,K)).GT.EPS.AND.J.LT.K) 
     &           RATE = RPRATE(2*LAMDA3(I,J,K)**2,MSQRT(2*I),
     &                       MQU(2*K-1),MQU(2*J-1))
              CALL RPMODA(RATE,PIN,POUT2,POUT1,0)
            ELSE
C--New for left right stop mixing
              IF(ABS(LAMDA3(I,J,K)).GT.EPS.AND.J.LT.K) THEN 
                RATE = RPRATE(2*LAMDA3(I,J,K)**2*TMIXSS(2,1)**2,
     &                     MSQLT(2*I),MQU(2*K-1),MQU(2*J-1))
                CALL RPMODA(RATE,26,POUT2,POUT1,0)
                RATE = RPRATE(2*LAMDA3(I,J,K)**2*TMIXSS(2,2)**2,
     &                     MSQRT(2*I),MQU(2*K-1),MQU(2*J-1))
                CALL RPMODA(RATE,46,POUT2,POUT1,0)
              ENDIF       
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C--Now right down squarks 
      DO I=1,3
        DO J=1,3
          DO K=1,3
C--via LUD 1st rate
            RATE = ZERO 
            IF(I.EQ.1) PIN   = 42
            IF(I.GT.1) PIN   = 39+2*I
            IF(K.EQ.1) POUT1 = 2
            IF(K.GT.1) POUT1 = 2*K-1
            IF(I.NE.3) THEN
              IF(ABS(LAMDA2(J,K,I)).GT.EPS) 
     &          RATE = RPRATE(LAMDA2(J,K,I)**2,MSQRT(2*I-1),
     &                        ZERO,MQU(2*K-1))
              CALL RPMODA(RATE,PIN,9+2*J,POUT1,0)
            ELSE
C--New for left/right sbottom mixing
              IF(ABS(LAMDA2(J,K,I)).GT.EPS) THEN
C--bug fix 09/04/02 by P.R.
                RATE =  RPRATE(LAMDA2(J,K,I)**2*BMIXSS(2,1)**2,
     &                         MSQLT(2*I-1),ZERO,MQU(2*K-1))
                CALL RPMODA(RATE,25,9+2*J,POUT1,0)
                RATE =  RPRATE(LAMDA2(J,K,I)**2*BMIXSS(2,2)**2,
     &                         MSQRT(2*I-1),ZERO,MQU(2*K-1))
                CALL RPMODA(RATE,45,9+2*J,POUT1,0)
              ENDIF
            ENDIF
C--via LUD 2nd rate
              RATE = ZERO
              IF(K.EQ.1) PIN = 42
              IF(K.GT.1) PIN = 39+2*K
              IF(J.EQ.1) POUT1 = 1
              IF(J.GT.1) POUT1 = 2*J 
              IF(K.NE.3) THEN
                IF(ABS(LAMDA2(I,J,K)).GT.EPS)
     &            RATE = RPRATE(LAMDA2(I,J,K)**2,MSQRT(2*K-1),
     &                          MLP(I),MQU(2*J))
                CALL RPMODA(RATE,PIN,10+2*I,POUT1,0) 
              ELSE
C-- New for left/right sbottom mixing
                IF(ABS(LAMDA2(I,J,K)).GT.EPS) THEN
                  RATE = RPRATE(LAMDA2(I,J,K)**2*BMIXSS(2,1)**2,
     &                           MSQLT(2*K-1),MLP(I),MQU(2*J))
                  CALL RPMODA(RATE,25,10+2*I,POUT1,0) 
C--bug fix 09/04/02 by P.R.
                  RATE = RPRATE(LAMDA2(I,J,K)**2*BMIXSS(2,2)**2,
     &                           MSQRT(2*K-1),MLP(I),MQU(2*J))
                  CALL RPMODA(RATE,45,10+2*I,POUT1,0) 
                ENDIF
              ENDIF
C--via UDD
              RATE = ZERO
              IF(I.EQ.1) PIN   = 42
              IF(I.GT.1) PIN   = 39+2*I
              IF(J.EQ.1) POUT1 = -1
              IF(J.GT.1) POUT1 = -2*J
              IF(K.EQ.1) POUT2 = -2
              IF(K.GT.1) POUT2 = 1-2*K
              IF(I.NE.3) THEN
                IF(ABS(LAMDA3(J,K,I)).GT.EPS) 
     &            RATE = RPRATE(2*LAMDA3(J,K,I)**2,MSQRT(2*I-1),
     &                            MQU(2*K-1),MQU(2*J))
                CALL RPMODA(RATE,PIN,POUT1,POUT2,0)
              ELSE
C--New for left right sbottom mixing
                IF(ABS(LAMDA3(J,K,I)).GT.EPS) THEN
                  RATE = RPRATE(2*LAMDA3(J,K,I)**2*BMIXSS(2,1)**2,
     &                            MSQLT(2*I-1),MQU(2*K-1),MQU(2*J))
                  CALL RPMODA(RATE,25,POUT1,POUT2,0)
                  RATE = RPRATE(2*LAMDA3(J,K,I)**2*BMIXSS(2,2)**2,
     &                            MSQRT(2*I-1),MQU(2*K-1),MQU(2*J))
                  CALL RPMODA(RATE,45,POUT1,POUT2,0)
                ENDIF
              ENDIF
          ENDDO
        ENDDO
      ENDDO
C------------END OF SCALAR DECAY RATES
C--Masses and Widths for the 3-body decays
C--First obtain all the widths we need and put in the arrays
      DO I=1,6
        IF(I.LE.3) SLRTWD(I) =  ZERO
        SLLTWD(I) =  ZERO
        SQLTWD(I) =  ZERO
        SQRTWD(I) =  ZERO
      ENDDO
      DO J=1,NSSMOD
        DO I=1,6
          IF(ISSMOD(J).EQ.(20+I)) SQLTWD(I) = SQLTWD(I)+GSSMOD(J)
          IF(ISSMOD(J).EQ.(30+I)) SLLTWD(I) = SLLTWD(I)+GSSMOD(J)
          IF(ISSMOD(J).EQ.(40+I)) SQRTWD(I) = SQRTWD(I)+GSSMOD(J)
          IF(ISSMOD(J).EQ.(50+2*I).AND.I.LE.3) SLRTWD(I) = 
     &                                        SLRTWD(I)+GSSMOD(J)
        ENDDO
      ENDDO
      DO J=1,NSSMD2
        DO I=1,6
          IF(ISSMD2(J).EQ.(20+I)) SQLTWD(I) = SQLTWD(I)+GSSMD2(J)
          IF(ISSMD2(J).EQ.(30+I)) SLLTWD(I) = SLLTWD(I)+GSSMD2(J)
          IF(ISSMD2(J).EQ.(40+I)) SQRTWD(I) = SQRTWD(I)+GSSMD2(J)
          IF(ISSMD2(J).EQ.(50+2*I).AND.I.LE.3) SLRTWD(I) =
     &                                        SLRTWD(I)+GSSMD2(J)
        ENDDO
      ENDDO
C--Change order of down and up to make loops easier
      RATE = SQLTWD(2)
      SQLTWD(2) = SQLTWD(1)
      SQLTWD(1) = RATE
      RATE = SQRTWD(2)
      SQRTWD(2) = SQRTWD(1)
      SQRTWD(1) = RATE
C--Now calculate the rates
C--Neutralino Decay Rates via R-parity violation
      DO CHANEL=1,4
        DO L=1,4
C--First calculate the charges we will need everywhere
          LCHAR(1) = -(ECHAR*NPRIME(L,1)/3.
     &                 +G*(.5-SN2THW/3)*NPRIME(L,2)/CWEAK)
          LCHAR(2) = ECHAR*NPRIME(L,1)*2/3.
     &                 +G*(.5-2*SN2THW/3)*NPRIME(L,2)/CWEAK
          LCHAR(3) = -(ECHAR*NPRIME(L,1)
     &                 +G*(.5-SN2THW)*NPRIME(L,2)/CWEAK)
          LCHAR(4) = G*NPRIME(L,2)/(2*CWEAK)
          RCHAR(1) = (ECHAR*NPRIME(L,1)
     &                -G*SN2THW*NPRIME(L,2)/CWEAK)/3.
          RCHAR(2) = -(ECHAR*NPRIME(L,1)-
     &                  G*SN2THW*NPRIME(L,2)/CWEAK)*2/3.
          RCHAR(3) =  ECHAR*NPRIME(L,1)
     &                -G*SN2THW*NPRIME(L,2)/CWEAK
          RCHAR(4) = ZERO
          MCHAR(1) = G*NPRIME(L,3)/(2*AMW*CBETA)
          MCHAR(2) = G*NPRIME(L,4)/(2*AMW*SBETA)
          DO I=1,3
            DO J=1,3
              DO K=1,3
                DO MIX=1,3
                  MIXING(2*MIX-1) = ONE
                  MIXING(2*MIX)   = ZERO
                  RESM(2*MIX-1)   = ZERO
                  RESM(2*MIX)     = ZERO
                  WIDTH(2*MIX-1)  = ZERO
                  WIDTH(2*MIX)    = ZERO
                ENDDO   
                NEUTWD = ZERO
                M(4) = ABS(AMZISS(L))
C--Charged lepton LQD mode
                IF(CHANEL.EQ.1) THEN
                  WIDTH(1)=SLLTWD(2*I)
                  WIDTH(3)=SQLTWD(2*J)
                  WIDTH(5)=SQRTWD(2*K-1)      
                  M(1) = MLP(I)
                  M(2) = MQU(2*J)
                  M(3) = MQU(2*K-1)
                  RESM(1) = MSLLT(2*I-1)
                  RESM(3) = MSQLT(2*J)
                  RESM(5) = MSQRT(2*K-1)
                  IF(I.NE.3) THEN
                    A(1) = MCHAR(1)*M(1)
                    B(1) = LCHAR(3)
                  ELSE
C--left/right stau mixing
                    DO MIX=1,2  
                      MIXING(MIX) = LMIXSS(1,MIX)
                      A(MIX) = MCHAR(1)*M(1)*LMIXSS(1,MIX)+
     &                         RCHAR(3)*LMIXSS(2,MIX)
                      B(MIX) = MCHAR(1)*M(1)*LMIXSS(2,MIX)+
     &                         LCHAR(3)*LMIXSS(1,MIX)
                    ENDDO         
                    RESM(2) = MSLRT(3)
                    WIDTH(2)= SLRTWD(3)
                  ENDIF
                  IF(J.NE.3) THEN
                    A(3) = MCHAR(2)*M(2)
                    B(3) = LCHAR(2)
                  ELSE
C--left/right stop mixing
                    DO MIX=1,2
                      MIXING(2+MIX) = TMIXSS(1,MIX) 
                      A(2+MIX) = MCHAR(2)*M(2)*TMIXSS(1,MIX)+
     &                           RCHAR(2)*TMIXSS(2,MIX)
                      B(2+MIX) = MCHAR(2)*M(2)*TMIXSS(2,MIX)+
     &                           LCHAR(2)*TMIXSS(1,MIX)     
                    ENDDO      
                    RESM(4) = MSQRT(2*J)
                    WIDTH(4) = SQRTWD(2*J)
                  ENDIF
                  IF(K.NE.3) THEN
                    A(5) = MCHAR(1)*M(3)
                    B(5) = RCHAR(1)
                  ELSE
C--left/right sbottom mixing
                    DO MIX=1,2  
                      MIXING(4+MIX) = BMIXSS(2,MIX)
                      A(4+MIX) = MCHAR(1)*M(3)*BMIXSS(2,MIX)+
     &                           LCHAR(1)*BMIXSS(1,MIX)
                      B(4+MIX) = MCHAR(1)*M(3)*BMIXSS(1,MIX)+
     &                           RCHAR(1)*BMIXSS(2,MIX)     
                    ENDDO  
                    RESM(5) = MSQLT(2*K-1)
                    RESM(6) = MSQRT(2*K-1)    
                    WIDTH(5)=SQLTWD(2*K-1)
                    WIDTH(6)=SQRTWD(2*K-1)      
                  ENDIF
                  LAMCOL = 6*LAMDA2(I,J,K)**2 
                  POUT1 = -10-2*I
                  IF(J.EQ.1) POUT2 = -1
                  IF(J.GT.1) POUT2 = -2*J
                  IF(K.EQ.1) POUT3 = 2
                  IF(K.GT.1) POUT3 = 2*K-1
C--neutrino LQD mode
                ELSEIF(CHANEL.EQ.2) THEN
                  WIDTH(1) = SLLTWD(2*I-1)
                  WIDTH(3) = SQLTWD(2*J-1)
                  WIDTH(5) = SQRTWD(2*K-1) 
                  M(1) = ZERO
                  M(2) = MQU(2*J-1)
                  M(3) = MQU(2*K-1)
                  RESM(1) = MSLLT(2*I)
                  RESM(3) = MSQLT(2*J-1)
C--bug fix 09/04/02 by P.R.
                  RESM(5) = MSQRT(2*K-1)
                  A(1) = ZERO
                  B(1) = LCHAR(4)
                  IF(J.NE.3) THEN
                    A(3) = MCHAR(1)*M(2)
                    B(3) = LCHAR(1)
                  ELSE
C--left/right sbottom mixing
                    DO MIX=1,2 
                      MIXING(2+MIX) = BMIXSS(1,MIX)
                      A(2+MIX) = MCHAR(1)*M(2)*BMIXSS(1,MIX)+
     &                           RCHAR(1)*BMIXSS(2,MIX)
                      B(2+MIX) = MCHAR(1)*M(2)*BMIXSS(2,MIX)+
     &                           LCHAR(1)*BMIXSS(1,MIX)
                    ENDDO      
                    RESM(4) = MSQRT(2*J-1)
                    WIDTH(4) = SQRTWD(2*J-1)
                  ENDIF
                  IF(K.NE.3) THEN
                    A(5) = MCHAR(1)*M(3)
                    B(5) = RCHAR(1)
                  ELSE
C--left/right sbottom mixing
                    DO MIX=1,2  
                      MIXING(4+MIX) = BMIXSS(2,MIX) 
                      A(4+MIX) = MCHAR(1)*M(3)*BMIXSS(2,MIX)+
     &                           LCHAR(1)*BMIXSS(1,MIX)
                      B(4+MIX) = MCHAR(1)*M(3)*BMIXSS(1,MIX)+
     &                           RCHAR(1)*BMIXSS(2,MIX)     
                    ENDDO  
                    RESM(5)= MSQLT(2*K-1)
                    RESM(6)= MSQRT(2*K-1)    
                    WIDTH(5)=SQLTWD(2*K-1)
                    WIDTH(6)=SQRTWD(2*K-1)      
                  ENDIF
                  LAMCOL = 6*LAMDA2(I,J,K)**2
                  POUT1 = -9-2*I 
                  IF(J.EQ.1) POUT2 = -2
                  IF(J.GT.1) POUT2 = 1-2*J
                  IF(K.EQ.1) POUT3 = 2
                  IF(K.GT.1) POUT3 = 2*K-1
C--LLE mode
                ELSEIF(CHANEL.EQ.3) THEN
                  WIDTH(1) = SLLTWD(2*I)
                  WIDTH(3) = SLLTWD(2*J-1)
                  WIDTH(5) = SLRTWD(K)
                  M(1) = MLP(I) 
                  M(2) = ZERO 
                  M(3) = MLP(K)
                  RESM(1) = MSLLT(2*I-1)
                  RESM(3) = MSLLT(2*J)
                  RESM(5) = MSLRT(K)
                  A(3) = ZERO
                  B(3) = LCHAR(4)
                  IF(I.NE.3) THEN
                    A(1) = MCHAR(1)*M(1)
                    B(1) = LCHAR(3)
                  ELSE
C--left/right stau mixing
                    DO MIX=1,2  
                      MIXING(MIX) = LMIXSS(1,MIX)
                      A(MIX) =MCHAR(1)*M(1)*LMIXSS(1,MIX)+
     &                        RCHAR(3)*LMIXSS(2,MIX)
                      B(MIX) =MCHAR(1)*M(1)*LMIXSS(2,MIX)+
     &                        LCHAR(3)*LMIXSS(1,MIX)
                    ENDDO         
                    RESM(2) = MSLRT(3)
                    WIDTH(2)= SLRTWD(3)
                  ENDIF
                  IF(K.NE.3) THEN
                    A(5) = MCHAR(1)*M(3)
                    B(5) = RCHAR(3)
                  ELSE
C--left/right stau mixing
                    DO MIX=1,2
                      MIXING(4+MIX) = LMIXSS(2,MIX)  
                      A(4+MIX) = MCHAR(1)*M(3)*LMIXSS(2,MIX)+
     &                           LCHAR(3)*LMIXSS(1,MIX)
                      B(4+MIX) = MCHAR(1)*M(3)*LMIXSS(1,MIX)+
     &                           RCHAR(3)*LMIXSS(2,MIX)
                    ENDDO  
                    RESM(5)= MSLLT(2*K-1)
                    RESM(6)= MSLRT(K)    
                    WIDTH(5)=SLLTWD(2*K)
                    WIDTH(6)=SLRTWD(K)      
                  ENDIF
                  LAMCOL = 2*LAMDA1(I,J,K)**2
                  POUT1 = -10-2*I
                  POUT2 =  -9-2*J
                  POUT3 = 10+2*K
C--UDD mode
                ELSEIF(CHANEL.EQ.4) THEN
                  WIDTH(1) = SQRTWD(2*I)
                  WIDTH(3) = SQRTWD(2*J-1)
                  WIDTH(5) = SQRTWD(2*K-1) 
                  M(1) = MQU(2*I)
                  M(2) = MQU(2*J-1)
                  M(3) = MQU(2*K-1)
                  RESM(1) = MSQRT(2*I)
                  RESM(3) = MSQRT(2*J-1)
                  RESM(5) = MSQRT(2*K-1)
                  IF(I.NE.3) THEN
                    A(1) = MCHAR(2)*M(1)
                    B(1) = RCHAR(2)
                  ELSE
C--left/right stop mixing
                    DO MIX=1,2  
                      MIXING(MIX) = TMIXSS(2,MIX)
                      A(MIX) = MCHAR(2)*M(1)*TMIXSS(2,MIX)+
     &                         LCHAR(2)*TMIXSS(1,MIX)
                      B(MIX) = MCHAR(2)*M(1)*TMIXSS(1,MIX)+
     &                         RCHAR(2)*TMIXSS(2,MIX)
                    ENDDO  
                    RESM(1) = MSQLT(2*I)
                    RESM(2) = MSQRT(2*I)
                    WIDTH(1) = SQLTWD(2*I)
                    WIDTH(2) = SQRTWD(2*I)
                  ENDIF
                  IF(J.NE.3) THEN
                    A(3) = MCHAR(1)*M(2)
                    B(3) = RCHAR(1)
                  ELSE
C--left/right sbottom mixing
                    DO MIX=1,2  
                      MIXING(2+MIX) = BMIXSS(2,MIX)
                      A(2+MIX) = MCHAR(1)*M(2)*BMIXSS(2,MIX)+
     &                           LCHAR(1)*BMIXSS(1,MIX)
                      B(2+MIX) = MCHAR(1)*M(2)*BMIXSS(1,MIX)+
     &                           RCHAR(1)*BMIXSS(2,MIX)  
                    ENDDO     
                    RESM(3)= MSQLT(2*J-1)
                    RESM(4)= MSQRT(2*J-1)    
                    WIDTH(3)=SQLTWD(2*J-1)
                    WIDTH(4)=SQRTWD(2*J-1)  
                  ENDIF
                  IF(K.NE.3) THEN
                    A(5) = MCHAR(1)*M(3)
                    B(5) = RCHAR(1)
                  ELSE
C--left/right sbottom mixing
                    DO MIX=1,2  
                      MIXING(4+MIX) = BMIXSS(2,MIX)
                      A(4+MIX) = MCHAR(1)*M(3)*BMIXSS(2,MIX)+
     &                           LCHAR(1)*BMIXSS(1,MIX)
                      B(4+MIX) = MCHAR(1)*M(3)*BMIXSS(1,MIX)+
     &                           RCHAR(1)*BMIXSS(2,MIX)     
                    ENDDO  
                    RESM(5)= MSQLT(2*K-1)
                    RESM(6)= MSQRT(2*K-1)    
                    WIDTH(5)=SQLTWD(2*K-1)
                    WIDTH(6)=SQRTWD(2*K-1)      
                  ENDIF
                  IF(I.EQ.1) POUT1 = 1
                  IF(I.GT.1) POUT1 = 2*I
                  IF(J.EQ.1) POUT2 = 2
                  IF(J.GT.1) POUT2 = 2*J-1
                  IF(K.EQ.1) POUT3 = 2
                  IF(K.GT.1) POUT3 = 2*K-1
                  IF(J.LT.K) THEN
                    LAMCOL = 12*LAMDA3(I,J,K)**2
                  ELSE
                    LAMCOL = ZERO
                  ENDIF
                ENDIF
                IF(NMSIGN(L).LT.0) THEN
                  IF(CHANEL.EQ.4) THEN
                    DO N=1,6
                      B(N) = -B(N)
                    ENDDO
                  ELSE
                    DO N=1,2
                      A(N)   = -A(N)
                      A(N+2) = -A(N+2)
                      B(N+4) = -B(N+4)
                    ENDDO
                  ENDIF
                ENDIF
C--Decide whether to remove diagrams
                DO N=1,6
                  CRSTRM(N) = .FALSE.
                ENDDO
                IF(M(4).LT.(M(1)+M(2)+M(3))) GOTO 10
                CRSTRM(1) = (M(4).GT.(M(1)+ABS(RESM(1))).
     &                       AND.ABS(RESM(1)).GT.(M(2)+M(3))) 
                CRSTRM(2) = (M(4).GT.(M(1)+ABS(RESM(2))).
     &                       AND.ABS(RESM(2)).GT.(M(2)+M(3))) 
     &                       .OR.(ABS(RESM(2)).LT.EPS)
                CRSTRM(3) = (M(4).GT.(M(2)+ABS(RESM(3))).
     &                       AND.ABS(RESM(3)).GT.(M(1)+M(3)))
                CRSTRM(4) = (M(4).GT.(M(2)+ABS(RESM(4))).
     &                       AND.ABS(RESM(4)).GT.(M(1)+M(3)))
     &                       .OR.(ABS(RESM(4)).LT.EPS)
                CRSTRM(5) = (M(4).GT.(M(3)+ABS(RESM(5))).
     &                       AND.ABS(RESM(5)).GT.(M(1)+M(2)))
                CRSTRM(6) = (M(4).GT.(M(3)+ABS(RESM(6))).
     &                       AND.ABS(RESM(6)).GT.(M(1)+M(2)))
     &                       .OR.(ABS(RESM(6)).LT.EPS)
C--Calculation of the rate
                NEUTWD = ZERO
                IF((CRSTRM(1).AND.CRSTRM(2).AND.CRSTRM(3).AND.
     &          CRSTRM(4).AND.CRSTRM(5).AND.CRSTRM(6)).
     &          OR.LAMCOL.LT.EPS) GOTO 10
C--first the diagram squared pieces
                IF(.NOT.CRSTRM(1)) NEUTWD = NEUTWD+ MIXING(1)**2*
     &          RPINF1(M(2),M(3),M(1),M(4),WIDTH(1),RESM(1),A(1),B(1))
                IF(.NOT.CRSTRM(2)) NEUTWD = NEUTWD+MIXING(2)**2*
     &          RPINF1(M(2),M(3),M(1),M(4),WIDTH(2),RESM(2),A(2),B(2))
                IF(.NOT.CRSTRM(3)) NEUTWD = NEUTWD+MIXING(3)**2*
     &          RPINF1(M(1),M(3),M(2),M(4),WIDTH(3),RESM(3),A(3),B(3))
                IF(.NOT.CRSTRM(4)) NEUTWD = NEUTWD+ MIXING(4)**2*
     &          RPINF1(M(1),M(3),M(2),M(4),WIDTH(4),RESM(4),A(4),B(4))
                IF(.NOT.CRSTRM(5)) NEUTWD = NEUTWD+MIXING(5)**2*
     &          RPINF1(M(1),M(2),M(3),M(4),WIDTH(5),RESM(5),A(5),B(5))
                IF(.NOT.CRSTRM(6)) NEUTWD = NEUTWD+MIXING(6)**2*
     &          RPINF1(M(1),M(2),M(3),M(4),WIDTH(6),RESM(6),A(6),B(6))
C--now for the light/heavy interference due left/right mixing
                IF(.NOT.CRSTRM(1).AND..NOT.CRSTRM(2)) NEUTWD=NEUTWD+
     &          MIXING(1)*MIXING(2)*RPINF2(M(2),M(3),M(1),M(4),WIDTH(1),
     &           WIDTH(2),RESM(1),RESM(2),A(1),A(2),B(1),B(2),1)
                IF(.NOT.CRSTRM(3).AND..NOT.CRSTRM(4)) NEUTWD=NEUTWD+
     &          MIXING(3)*MIXING(4)*RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),
     &          WIDTH(4),RESM(3),RESM(4),A(3),A(4),B(3),B(4),1)
                IF(.NOT.CRSTRM(5).AND..NOT.CRSTRM(6)) NEUTWD=NEUTWD+
     &          MIXING(5)*MIXING(6)*RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),
     &          WIDTH(6),RESM(5),RESM(6),A(5),A(6),B(5),B(6),1)
C--now for the true interference terms
                IF(.NOT.CRSTRM(1)) THEN
                  IF(.NOT.CRSTRM(3)) NEUTWD=NEUTWD-MIXING(1)*MIXING(3)*
     &            RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(1),RESM(3),
     &            RESM(1),A(1),A(3),B(1),B(3),2)
                  IF(.NOT.CRSTRM(4)) NEUTWD=NEUTWD-MIXING(1)*MIXING(4)*
     &            RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(1),RESM(4),
     &            RESM(1),A(1),A(4),B(1),B(4),2)
                  IF(.NOT.CRSTRM(5)) NEUTWD=NEUTWD-MIXING(1)*MIXING(5)*
     &            RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),WIDTH(1),RESM(5),
     &            RESM(1),A(1),A(5),B(1),B(5),2)
                  IF(.NOT.CRSTRM(6)) NEUTWD=NEUTWD-MIXING(1)*MIXING(6)*
     &            RPINF2(M(1),M(2),M(3),M(4),WIDTH(6),WIDTH(1),RESM(6),
     &            RESM(1),A(1),A(6),B(1),B(6),2)
                ENDIF
                IF(.NOT.CRSTRM(2)) THEN
                  IF(.NOT.CRSTRM(3)) NEUTWD=NEUTWD-MIXING(2)*MIXING(3)*
     &            RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(2),RESM(3),
     &            RESM(2),A(2),A(3),B(2),B(3),2)
                  IF(.NOT.CRSTRM(4)) NEUTWD=NEUTWD-MIXING(2)*MIXING(4)*
     &            RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(2),RESM(4),
     &            RESM(2),A(2),A(4),B(2),B(4),2)
                  IF(.NOT.CRSTRM(5)) NEUTWD=NEUTWD-MIXING(2)*MIXING(5)*
     &            RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),WIDTH(2),RESM(5),
     &            RESM(2),A(2),A(5),B(2),B(5),2)
                  IF(.NOT.CRSTRM(6)) NEUTWD=NEUTWD-MIXING(2)*MIXING(6)*
     &            RPINF2(M(1),M(2),M(3),M(4),WIDTH(6),WIDTH(2),RESM(6),
     &            RESM(2),A(2),A(6),B(2),B(6),2)
                ENDIF
                IF(.NOT.CRSTRM(3)) THEN
                  IF(.NOT.CRSTRM(5)) NEUTWD=NEUTWD-MIXING(3)*MIXING(5)*
     &            RPINF2(M(2),M(1),M(3),M(4),WIDTH(5),WIDTH(3),RESM(5),
     &            RESM(3),A(3),A(5),B(3),B(5),2)
                  IF(.NOT.CRSTRM(6)) NEUTWD=NEUTWD-MIXING(3)*MIXING(6)*
     &            RPINF2(M(2),M(1),M(3),M(4),WIDTH(6),WIDTH(3),RESM(6),
     &            RESM(3),A(3),A(6),B(3),B(6),2)
                ENDIF
                IF(.NOT.CRSTRM(4)) THEN
                  IF(.NOT.CRSTRM(5)) NEUTWD=NEUTWD-MIXING(4)*MIXING(5)*
     &            RPINF2(M(2),M(1),M(3),M(4),WIDTH(5),WIDTH(4),RESM(5),
     &            RESM(4),A(4),A(5),B(4),B(5),2)
                  IF(.NOT.CRSTRM(6)) NEUTWD=NEUTWD-MIXING(4)*MIXING(6)*
     &            RPINF2(M(2),M(1),M(3),M(4),WIDTH(6),WIDTH(4),RESM(6),
     &            RESM(4),A(4),A(6),B(4),B(6),2)
                ENDIF
                NEUTWD = LAMCOL*NEUTWD/((2*PI)**3*32*M(4)**3)
C--Output the rate and particles to decay tables
 10             IF(NEUTWD.GT.EPS) THEN
                  CALL RPMODA(NEUTWD,20+10*L,POUT1,POUT2,POUT3)
                  CALL RPMODA(NEUTWD,20+10*L,-POUT1,-POUT2,-POUT3)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C--Gluino decay rates all the info we need is already in the arrays
C--Easier as only three gluino channels. Just calculate the total rate
C--and handle the colour structure inside HERWIG
C--Now calculate the rates
C--First the LQD rates
      DO CHANEL=1,2
        DO I=1,3
          DO J=1,3
            DO K=1,3
              DO MIX=1,2
                MIXING(2*MIX-1) = ONE
                MIXING(2*MIX)   = ZERO
                RESM(2*MIX-1)   = ZERO
                RESM(2*MIX)     = ZERO
                WIDTH(2*MIX-1)  = ZERO
                WIDTH(2*MIX)    = ZERO 
              ENDDO
C--charged lepton, ubar, down
              IF(CHANEL.EQ.1) THEN
                WIDTH(1)=SQLTWD(2*J)
                WIDTH(3)=SQRTWD(2*K-1)      
                M(1) = MLP(I)
                M(2) = MQU(2*J)
                M(3) = MQU(2*K-1)
                M(4) = AMGLSS
                RESM(1) = MSQLT(2*J)
                RESM(3) = MSQRT(2*K-1)
                IF(J.NE.3) THEN
                  A(1) =ZERO 
                  B(1) = -ONE 
                ELSE
C--left/right stop mixing
                  DO MIX=1,2
                    MIXING(MIX) = TMIXSS(1,MIX)  
                    A(MIX) = TMIXSS(2,MIX)
                    B(MIX) =-TMIXSS(1,MIX)
                  ENDDO    
C--bug fix 09/04/02 by P.R.  
                  RESM(2) = MSQRT(2*J)
                  WIDTH(2) = SQRTWD(2*J)
                ENDIF
                IF(K.NE.3) THEN
                  A(3) = ZERO
                  B(3) = ONE
                ELSE
C--left/right sbottom mixing
                  MIXING(3) = BMIXSS(2,1) 
                  MIXING(4) = BMIXSS(2,2)
                  DO MIX=1,2 
                    MIXING(2+MIX) = BMIXSS(2,MIX) 
                    A(2+MIX) = -BMIXSS(1,MIX) 
                    B(2+MIX) =  BMIXSS(2,MIX)
                  ENDDO 
                  RESM(3)= MSQLT(2*K-1)
                  RESM(4)= MSQRT(2*K-1)    
                  WIDTH(3)=SQLTWD(2*K-1)
                  WIDTH(4)=SQRTWD(2*K-1)      
                ENDIF
                LAMCOL = ALFA3*4*PI*LAMDA2(I,J,K)**2 
                POUT1 = -10-2*I
                IF(J.EQ.1) POUT2 = -1
                IF(J.GT.1) POUT2 = -2*J
                IF(K.EQ.1) POUT3 = 2
                IF(K.GT.1) POUT3 = 2*K-1
C--neutrino, dbar, down.
              ELSEIF(CHANEL.EQ.2) THEN
                WIDTH(1) = SQLTWD(2*J-1)
                WIDTH(3) = SQRTWD(2*K-1) 
                M(1) = ZERO
                M(2) = MQU(2*J-1)
                M(3) = MQU(2*K-1)
                M(4) = AMGLSS
                RESM(1) = MSQLT(2*J-1)
C--bug fix 09/04/02 by P.R. 
                RESM(3) = MSQRT(2*K-1)
                IF(J.NE.3) THEN
                  A(1) =  ZERO
                  B(1) =-ONE
                ELSE
C--left/right sbottom mixing
                  DO MIX=1,2
                    MIXING(MIX) = BMIXSS(1,MIX) 
                    A(MIX) = BMIXSS(2,MIX)
                    B(MIX) = -BMIXSS(1,MIX) 
                  ENDDO      
                  RESM(2) = MSQRT(2*J-1)
                  WIDTH(2) = SQRTWD(2*J-1)
                ENDIF
                IF(K.NE.3) THEN
                  A(3) = ZERO
                  B(3) = ONE
                ELSE
C--left/right sbottom mixing
                  DO MIX=1,2 
                    MIXING(2+MIX) = BMIXSS(2,MIX)  
                    A(2+MIX) = -BMIXSS(1,MIX) 
                    B(2+MIX) =  BMIXSS(2,MIX)
                  ENDDO  
                  RESM(3)= MSQLT(2*K-1)
                  RESM(4)= MSQRT(2*K-1)    
                  WIDTH(3)=SQLTWD(2*K-1)
                  WIDTH(4)=SQRTWD(2*K-1)      
                ENDIF
                LAMCOL = ALFA3*4*PI*LAMDA2(I,J,K)**2
                POUT1 = -9-2*I 
                IF(J.EQ.1) POUT2 = -2
                IF(J.GT.1) POUT2 = 1-2*J
                IF(K.EQ.1) POUT3 = 2
                IF(K.GT.1) POUT3 = 2*K-1
              ENDIF
C--Decide whether to remove diagrams
              DO N=1,4
                CRSTRM(N) = .FALSE.
              ENDDO
C--bug fix 09/04/02 by P.R.
              CRSTRM(1) = (M(4).GT.(M(2)+ABS(RESM(1))).
     &                     AND.ABS(RESM(1)).GT.(M(1)+M(3)))
              CRSTRM(2) = (M(4).GT.(M(2)+ABS(RESM(2))).
     &                       AND.ABS(RESM(2)).GT.(M(1)+M(3)))
     &                       .OR.(ABS(RESM(2)).LT.EPS)
              CRSTRM(3) = (M(4).GT.(M(3)+ABS(RESM(3))).
     &                       AND.ABS(RESM(3)).GT.(M(1)+M(2)))
              CRSTRM(4) = (M(4).GT.(M(3)+ABS(RESM(4))).
     &                       AND.ABS(RESM(4)).GT.(M(1)+M(2)))
     &                       .OR.(ABS(RESM(4)).LT.EPS)
C--Calculation of the rate
              GLUWD = 0
              IF((CRSTRM(1).AND.CRSTRM(2).AND.
     &            CRSTRM(3).AND.CRSTRM(4)).OR.
     &           LAMCOL.LT.EPS.OR.(M(1)+M(2)+M(3)).GT.M(4)) GOTO 20
C--First the amplitude square pieces
              IF(.NOT.CRSTRM(1)) GLUWD = GLUWD+MIXING(1)**2*
     &        RPINF1(M(1),M(3),M(2),M(4),WIDTH(1),RESM(1),A(1),B(1))
              IF(.NOT.CRSTRM(2)) GLUWD = GLUWD+ MIXING(2)**2*
     &        RPINF1(M(1),M(3),M(2),M(4),WIDTH(2),RESM(2),A(2),B(2))
              IF(.NOT.CRSTRM(3)) GLUWD = GLUWD+MIXING(3)**2*
     &        RPINF1(M(1),M(2),M(3),M(4),WIDTH(3),RESM(3),A(3),B(3))
              IF(.NOT.CRSTRM(4)) GLUWD = GLUWD+MIXING(4)**2*
     &        RPINF1(M(1),M(2),M(3),M(4),WIDTH(4),RESM(4),A(4),B(4))
C--now for the light/heavy interference due left/right mixing
              IF(.NOT.CRSTRM(1).AND..NOT.CRSTRM(2)) GLUWD=GLUWD+
     &        MIXING(1)*MIXING(2)*RPINF2(M(1),M(3),M(2),M(4),WIDTH(1),
     &        WIDTH(2),RESM(1),RESM(2),A(1),A(2),B(1),B(2),1)
              IF(.NOT.CRSTRM(3).AND..NOT.CRSTRM(4)) GLUWD=GLUWD+
     &        MIXING(3)*MIXING(4)*RPINF2(M(1),M(2),M(3),M(4),WIDTH(3),
     &        WIDTH(4),RESM(3),RESM(4),A(3),A(4),B(3),B(4),1)
C--now for the true interference terms
              IF(.NOT.CRSTRM(3)) THEN
                IF(.NOT.CRSTRM(3)) GLUWD=GLUWD-MIXING(1)*MIXING(3)*
     &          RPINF2(M(2),M(1),M(3),M(4),WIDTH(3),WIDTH(1),RESM(3),
     &          RESM(1),A(1),A(3),B(1),B(3),2)
                IF(.NOT.CRSTRM(4)) GLUWD=GLUWD-MIXING(1)*MIXING(4)*
     &          RPINF2(M(2),M(1),M(3),M(4),WIDTH(4),WIDTH(1),RESM(4),
     &          RESM(1),A(1),A(4),B(1),B(4),2)
              ENDIF
              IF(.NOT.CRSTRM(2)) THEN
                IF(.NOT.CRSTRM(3)) GLUWD=GLUWD-MIXING(2)*MIXING(3)*
     &          RPINF2(M(2),M(1),M(3),M(4),WIDTH(3),WIDTH(2),RESM(3),
     &          RESM(2),A(2),A(3),B(2),B(3),2)
                IF(.NOT.CRSTRM(4)) GLUWD=GLUWD-MIXING(2)*MIXING(4)*
     &          RPINF2(M(2),M(1),M(3),M(4),WIDTH(4),WIDTH(2),RESM(4),
     &          RESM(2),A(2),A(4),B(2),B(4),2)
              ENDIF
              GLUWD = LAMCOL*GLUWD
              GLUWD = GLUWD/((2*PI)**3*32*M(4)**3)
C--Output the rate and particles to decay tables
 20           IF(GLUWD.GT.EPS) THEN
                CALL RPMODA(GLUWD,29,POUT1,POUT2,POUT3)
                CALL RPMODA(GLUWD,29,-POUT1,-POUT2,-POUT3)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C--Now the gluino decay via UDD, there are three diagrams before LR mixing
C--i.e. gluino -> u d d
      DO I=1,3
        DO J=1,3
          DO K=1,3
            DO MIX=1,3
              MIXING(2*MIX-1) = ONE
              MIXING(2*MIX)   = ZERO
              RESM(2*MIX-1)   = ZERO
              RESM(2*MIX)     = ZERO
              WIDTH(2*MIX-1)  = ZERO
              WIDTH(2*MIX)    = ZERO 
            ENDDO
            WIDTH(1) = SQRTWD(2*I)
            WIDTH(3) = SQRTWD(2*J-1)
            WIDTH(5) = SQRTWD(2*K-1) 
            M(1) = MQU(2*I)
            M(2) = MQU(2*J-1)
            M(3) = MQU(2*K-1)
            M(4) = AMGLSS
            RESM(1) = MSQRT(2*I)
            RESM(3) = MSQRT(2*J-1)
            RESM(5) = MSQRT(2*K-1)
            IF(I.NE.3) THEN
              A(1) = ZERO
              B(1) = ONE
            ELSE
C--left/right stop mixing
              DO MIX=1,2
                MIXING(MIX) = TMIXSS(2,MIX)  
                A(MIX) = -TMIXSS(1,MIX)
                B(MIX) = TMIXSS(2,MIX)
              ENDDO  
              RESM(1) = MSQLT(2*I)
              RESM(2) = MSQRT(2*I)
              WIDTH(1) = SQLTWD(2*I)
              WIDTH(2) = SQRTWD(2*I)
            ENDIF
            IF(J.NE.3) THEN
              A(3) = ZERO
              B(3) = ONE
            ELSE
C--left/right sbottom mixing
              DO MIX=1,2  
                MIXING(2+MIX) = BMIXSS(2,MIX)
                A(2+MIX) =-BMIXSS(1,MIX)
                B(2+MIX) = BMIXSS(2,MIX)
              ENDDO     
              RESM(3)= MSQLT(2*J-1)
              RESM(4)= MSQRT(2*J-1)    
              WIDTH(3)=SQLTWD(2*J-1)
              WIDTH(4)=SQRTWD(2*J-1)  
            ENDIF
            IF(K.NE.3) THEN
              A(5) = ZERO
              B(5) = ONE 
            ELSE
C--left/right sbottom mixing
            DO MIX=1,2  
              MIXING(4+MIX)=BMIXSS(2,MIX)
              A(4+MIX) =-BMIXSS(1,MIX)
              B(4+MIX) = BMIXSS(2,MIX)
            ENDDO  
              RESM(5)= MSQLT(2*K-1)
              RESM(6)= MSQRT(2*K-1)    
              WIDTH(5)=SQLTWD(2*K-1)
              WIDTH(6)=SQRTWD(2*K-1)      
            ENDIF
            IF(I.EQ.1) POUT1 = 1
            IF(I.GT.1) POUT1 = 2*I
            IF(J.EQ.1) POUT2 = 2
            IF(J.GT.1) POUT2 = 2*J-1
            IF(K.EQ.1) POUT3 = 2
            IF(K.GT.1) POUT3 = 2*K-1
            IF(J.LT.K) THEN
              LAMCOL = 2*ALFA3*4*PI*LAMDA3(I,J,K)**2
            ELSE
              LAMCOL = 0
            ENDIF
C--Decide whether to remove diagrams
            DO N=1,6
              CRSTRM(N) = .FALSE.
            ENDDO
            GLUWD = ZERO
            IF(M(4).LT.(M(1)+M(2)+M(3))) GOTO 30
            CRSTRM(1) = (M(4).GT.(M(1)+ABS(RESM(1))).
     &                  AND.ABS(RESM(1)).GT.(M(2)+M(3))) 
            CRSTRM(2) = (M(4).GT.(M(1)+ABS(RESM(2))).
     &                  AND.ABS(RESM(2)).GT.(M(2)+M(3))) 
     &                  .OR.(ABS(RESM(2)).LT.EPS)
            CRSTRM(3) = (M(4).GT.(M(2)+ABS(RESM(3))).
     &                  AND.ABS(RESM(3)).GT.(M(1)+M(3)))
            CRSTRM(4) = (M(4).GT.(M(2)+ABS(RESM(4))).
     &                  AND.ABS(RESM(4)).GT.(M(1)+M(3)))
     &                  .OR.(ABS(RESM(4)).LT.EPS)
            CRSTRM(5) = (M(4).GT.(M(3)+ABS(RESM(5))).
     &                  AND.ABS(RESM(5)).GT.(M(1)+M(2)))
            CRSTRM(6) = (M(4).GT.(M(3)+ABS(RESM(6))).
     &                  AND.ABS(RESM(6)).GT.(M(1)+M(2)))
     &                  .OR.(ABS(RESM(6)).LT.EPS)
C--Calculation of the rate
            IF((CRSTRM(1).AND.CRSTRM(2).AND.CRSTRM(3).AND.
     &      CRSTRM(4).AND.CRSTRM(5).AND.CRSTRM(6)).
     &      OR.LAMCOL.LT.EPS) GOTO 30
C--first the diagram squared pieces
            IF(.NOT.CRSTRM(1)) GLUWD = GLUWD+ MIXING(1)**2*
     &      RPINF1(M(2),M(3),M(1),M(4),WIDTH(1),RESM(1),A(1),B(1))
            IF(.NOT.CRSTRM(2)) GLUWD = GLUWD+MIXING(2)**2*
     &      RPINF1(M(2),M(3),M(1),M(4),WIDTH(2),RESM(2),A(2),B(2))
            IF(.NOT.CRSTRM(3)) GLUWD = GLUWD+MIXING(3)**2*
     &      RPINF1(M(1),M(3),M(2),M(4),WIDTH(3),RESM(3),A(3),B(3))
            IF(.NOT.CRSTRM(4)) GLUWD = GLUWD+ MIXING(4)**2*
     &      RPINF1(M(1),M(3),M(2),M(4),WIDTH(4),RESM(4),A(4),B(4))
            IF(.NOT.CRSTRM(5)) GLUWD = GLUWD+MIXING(5)**2*
     &      RPINF1(M(1),M(2),M(3),M(4),WIDTH(5),RESM(5),A(5),B(5))
            IF(.NOT.CRSTRM(6)) GLUWD = GLUWD+MIXING(6)**2*
     &      RPINF1(M(1),M(2),M(3),M(4),WIDTH(6),RESM(6),A(6),B(6))
C--now for the light/heavy interference due left/right mixing
            IF(.NOT.CRSTRM(1).AND..NOT.CRSTRM(2)) GLUWD=GLUWD+
     &      MIXING(1)*MIXING(2)*RPINF2(M(2),M(3),M(1),M(4),WIDTH(1),
     &      WIDTH(2),RESM(1),RESM(2),A(1),A(2),B(1),B(2),1)
            IF(.NOT.CRSTRM(3).AND..NOT.CRSTRM(4)) GLUWD=GLUWD+
     &      MIXING(3)*MIXING(4)*RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),
     &      WIDTH(4),RESM(3),RESM(4),A(3),A(4),B(3),B(4),1)
            IF(.NOT.CRSTRM(5).AND..NOT.CRSTRM(6)) GLUWD=GLUWD+
     &      MIXING(5)*MIXING(6)*RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),
     &      WIDTH(6),RESM(5),RESM(6),A(5),A(6),B(5),B(6),1)
C--now for the true interference terms
            IF(.NOT.CRSTRM(1)) THEN
              IF(.NOT.CRSTRM(3)) GLUWD=GLUWD+0.5*MIXING(1)*MIXING(3)*
     &        RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(1),RESM(3),
     &        RESM(1),A(1),A(3),B(1),B(3),2)
              IF(.NOT.CRSTRM(4)) GLUWD=GLUWD+0.5*MIXING(1)*MIXING(4)*
     &        RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(1),RESM(4),
     &        RESM(1),A(1),A(4),B(1),B(4),2)
              IF(.NOT.CRSTRM(5)) GLUWD=GLUWD+0.5*MIXING(1)*MIXING(5)*
     &        RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),WIDTH(1),RESM(5),
     &        RESM(1),A(1),A(5),B(1),B(5),2)
              IF(.NOT.CRSTRM(6)) GLUWD=GLUWD+0.5*MIXING(1)*MIXING(6)*
     &        RPINF2(M(1),M(2),M(3),M(4),WIDTH(6),WIDTH(1),RESM(6),
     &        RESM(1),A(1),A(6),B(1),B(6),2)
            ENDIF
            IF(.NOT.CRSTRM(2)) THEN
              IF(.NOT.CRSTRM(3)) GLUWD=GLUWD+0.5*MIXING(2)*MIXING(3)*
     &        RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(2),RESM(3),
     &        RESM(2),A(2),A(3),B(2),B(3),2)
              IF(.NOT.CRSTRM(4)) GLUWD=GLUWD+0.5*MIXING(2)*MIXING(4)*
     &        RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(2),RESM(4),
     &        RESM(2),A(2),A(4),B(2),B(4),2)
              IF(.NOT.CRSTRM(5)) GLUWD=GLUWD+0.5*MIXING(2)*MIXING(5)*
     &        RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),WIDTH(2),RESM(5),
     &        RESM(2),A(2),A(5),B(2),B(5),2)
              IF(.NOT.CRSTRM(6)) GLUWD=GLUWD+0.5*MIXING(2)*MIXING(6)*
     &        RPINF2(M(1),M(2),M(3),M(4),WIDTH(6),WIDTH(2),RESM(6),
     &        RESM(2),A(2),A(6),B(2),B(6),2)
            ENDIF
            IF(.NOT.CRSTRM(3)) THEN
              IF(.NOT.CRSTRM(5)) GLUWD=GLUWD+0.5*MIXING(3)*MIXING(5)*
     &        RPINF2(M(2),M(1),M(3),M(4),WIDTH(5),WIDTH(3),RESM(5),
     &        RESM(3),A(3),A(5),B(3),B(5),2)
              IF(.NOT.CRSTRM(6)) GLUWD=GLUWD+0.5*MIXING(3)*MIXING(6)*
     &        RPINF2(M(2),M(1),M(3),M(4),WIDTH(6),WIDTH(3),RESM(6),
     &        RESM(3),A(3),A(6),B(3),B(6),2)
            ENDIF
            IF(.NOT.CRSTRM(4)) THEN
              IF(.NOT.CRSTRM(5)) GLUWD=GLUWD+0.5*MIXING(4)*MIXING(5)*
     &        RPINF2(M(2),M(1),M(3),M(4),WIDTH(5),WIDTH(4),RESM(5),
     &        RESM(4),A(4),A(5),B(4),B(5),2)
              IF(.NOT.CRSTRM(6)) GLUWD=GLUWD+0.5*MIXING(4)*MIXING(6)*
     &        RPINF2(M(2),M(1),M(3),M(4),WIDTH(6),WIDTH(4),RESM(6),
     &        RESM(4),A(4),A(6),B(4),B(6),2)
            ENDIF
            GLUWD = LAMCOL*GLUWD
            GLUWD = GLUWD/((2*PI)**3*32*M(4)**3)
C--Output the rate and particles to decay tables
 30         IF(GLUWD.GT.EPS) THEN
              CALL RPMODA(GLUWD,29,POUT1,POUT2,POUT3)
              CALL RPMODA(GLUWD,29,-POUT1,-POUT2,-POUT3)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C--Decay modes of the chargino 
C--we will do the one diagram via LLE and LQD modes first
      MCHAR(1) = 1./(SQRT(2.)*AMW*CBETA)
      MCHAR(2) = 1./(SQRT(2.)*AMW*SBETA)
      DO L=1,2
        DO CHANEL=1,3
          DO I=1,3
            DO J=1,3
              DO K=1,3
                MIXING(1) = ONE
                MIXING(2) = ZERO
                DO MIX=1,2
                  RESM(MIX) = ZERO
                  WIDTH(MIX) = ZERO
                ENDDO
                IF(CHANEL.EQ.1) THEN
C--LLE decay to charged lepton, neutrino, antineutrino
                  M(1) = ZERO
                  M(2) = MLP(J)
                  M(3) = ZERO
                  M(4) = ABS(CHARM(L))
                  RESM(1) = MSLRT(K)
                  WIDTH(1) = SLRTWD(K)
                  IF(K.NE.3) THEN
                    A(1) = -WMXUSS(L,2)*MLP(K)*MCHAR(1)
                    B(1) = ZERO
                  ELSE
C--left/right stau mixing
                    DO MIX=1,2  
                      MIXING(MIX) = LMIXSS(2,MIX)
                      A(MIX) = WMXUSS(L,1)*LMIXSS(1,MIX)-
     &             WMXUSS(L,2)*MLP(K)*MCHAR(1)*LMIXSS(2,MIX)
                      B(MIX) = ZERO
                    ENDDO     
                    RESM(1)= MSLLT(2*K-1)
                    RESM(2)= MSLRT(K)    
                    WIDTH(1)=SLLTWD(2*K)
                    WIDTH(2)=SLRTWD(K) 
                  ENDIF
                  LAMCOL = G**2*LAMDA1(I,J,K)**2
                  POUT1 = -9-2*I 
                  POUT2 = -10-2*J  
                  POUT3 = 9+2*K
                ELSEIF(CHANEL.EQ.2) THEN
C--LQD decay to antineutrino, dbar, up
                  M(1) = ZERO
                  M(2) = MQU(2*J-1)
                  M(3) = MQU(2*K)
                  M(4) = ABS(CHARM(L))
                  RESM(1) = MSQRT(2*K-1)
                  WIDTH(1)= SQRTWD(2*K-1)
                  IF(K.NE.3) THEN
                    A(1) = -WMXUSS(L,2)*MQU(2*K-1)*MCHAR(1)
                    B(1) = 0
                  ELSE
                    DO MIX=1,2
                      MIXING(MIX) = BMIXSS(2,MIX)
                      A(MIX) = WMXUSS(L,1)*BMIXSS(1,MIX)-
     &            WMXUSS(L,2)*BMIXSS(2,MIX)*MQU(2*K-1)*MCHAR(1)
                      B(MIX) = -MQU(2*K)*WMXVSS(L,2)*BMIXSS(1,MIX)*
     &                         MCHAR(2) 
                      IF(CMSIGN(L).LT.0) B(MIX)=-B(MIX)
                    ENDDO
                    RESM(1)= MSQLT(2*K-1)
                    RESM(2)= MSQRT(2*K-1)    
                    WIDTH(1)=SQLTWD(2*K-1)
                    WIDTH(2)=SQRTWD(2*K-1)
                  ENDIF
                  LAMCOL = 3.*G**2*LAMDA2(I,J,K)**2
                  POUT1 = -9-2*I 
                  POUT2 = 1-2*J
                  POUT3 = 2*K
                  IF(J.EQ.1) POUT2 = -2
                  IF(K.EQ.1) POUT3 = 1
                ELSEIF(CHANEL.EQ.3) THEN
C--LQD decay to charged lepton, ubar, up
                  M(1) = MLP(I)
                  M(2) = MQU(2*J)
                  M(3) = MQU(2*K)
                  M(4) = ABS(CHARM(L))
                  RESM(1) = MSQRT(2*K-1)
                  WIDTH(1)= SQRTWD(2*K-1)
                  IF(K.NE.3) THEN
                    A(1) = -WMXUSS(L,2)*MQU(2*K-1)*MCHAR(1)
                    B(1) = ZERO
                  ELSE
                    DO MIX=1,2
                      MIXING(MIX) = BMIXSS(2,MIX)
                      A(MIX) = WMXUSS(L,1)*BMIXSS(1,MIX)-
     &            WMXUSS(L,2)*BMIXSS(2,MIX)*MQU(2*K-1)*MCHAR(1)
                      B(MIX) = -MQU(2*K)*WMXVSS(L,2)*BMIXSS(1,MIX)*
     &                         MCHAR(2)
                      IF(CMSIGN(L).LT.0) B(MIX)=-B(MIX)
                    ENDDO
                    RESM(1)= MSQLT(2*K-1)
                    RESM(2)= MSQRT(2*K-1)    
                    WIDTH(1)=SQLTWD(2*K-1)
                    WIDTH(2)=SQRTWD(2*K-1) 
                  ENDIF
                  LAMCOL = 3.*G**2*LAMDA2(I,J,K)**2
                  POUT1 = -10-2*I 
                  POUT2 = -2*J
                  POUT3 = 2*K
                  IF(J.EQ.1) POUT2 = -1
                  IF(K.EQ.1) POUT3 = 1
                ENDIF  
                CHARWD = ZERO
                IF(M(4).LT.(M(1)+M(2)+M(3))) GOTO 40
C--Decide whether to remove diagrams
                DO N=1,2
                  CRSTRM(N) = .FALSE.
                ENDDO
                CRSTRM(1) = (M(4).GT.(M(3)+ABS(RESM(1))).
     &                      AND.ABS(RESM(1)).GT.(M(1)+M(2))) 
                CRSTRM(2) = (M(4).GT.(M(3)+ABS(RESM(2))).
     &                      AND.ABS(RESM(2)).GT.(M(1)+M(2))) 
     &                     .OR.(ABS(RESM(2)).LT.EPS)
C--Calculation of the rate
               IF((CRSTRM(1).AND.CRSTRM(2)).OR.LAMCOL.LT.EPS) GOTO 40
C--first the diagram squared pieces
                IF(.NOT.CRSTRM(1)) CHARWD = CHARWD+ MIXING(1)**2*
     &          RPINF1(M(1),M(2),M(3),M(4),WIDTH(1),RESM(1),A(1),B(1))
                IF(.NOT.CRSTRM(2)) CHARWD = CHARWD+MIXING(2)**2*
     &          RPINF1(M(1),M(2),M(3),M(4),WIDTH(2),RESM(2),A(2),B(2))
C--now for the light/heavy interference due left/right mixing
                IF(.NOT.CRSTRM(1).AND..NOT.CRSTRM(2)) CHARWD=CHARWD+
     &          MIXING(1)*MIXING(2)*RPINF2(M(1),M(2),M(3),M(4),WIDTH(1),
     &          WIDTH(2),RESM(1),RESM(2),A(1),A(2),B(1),B(2),1)
                CHARWD = LAMCOL*CHARWD/((2*PI)**3*32*M(4)**3)
 40             IF(CHARWD.GT.EPS) CALL RPMODA(CHARWD,29+10*L,
     &                                          POUT1,POUT2,POUT3)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C--Now for the two diagram LLE and LQD modes, these tend to have
C--higher branching ratios
      DO L=1,2
        DO CHANEL=1,4
          DO I=1,3
            DO J=1,3
              DO K=1,3
                DO N=1,6
                  A(N) = ZERO
                  B(N) = ZERO
                  RESM(N) = ZERO
                  WIDTH(N) = ZERO
                ENDDO
                DO MIX=1,2
                  MIXING(2*MIX-1) = 1.
                  MIXING(2*MIX) = ZERO
                ENDDO
C--LLE to neutrino, neutrino, charged lepton
                IF(CHANEL.EQ.1) THEN
                  M(1) = ZERO
                  M(2) = ZERO
                  M(3) = MLP(K)
                  M(4) = ABS(CHARM(L))
                  RESM(1) = MSLLT(2*I-1)
                  RESM(3) = MSLLT(2*J-1)
                  WIDTH(1)= SLLTWD(2*I)
                  WIDTH(3)= SLLTWD(2*J)
                  IF(I.NE.3) THEN
                    B(1) = WMXUSS(L,1)
                  ELSE
                    DO MIX=1,2
                      MIXING(MIX) = LMIXSS(1,MIX)
                      B(MIX) = WMXUSS(L,1)*LMIXSS(1,MIX)-
     &                WMXUSS(L,2)*LMIXSS(2,MIX)*MLP(I)*MCHAR(1)
                    ENDDO
                    RESM(2)= MSLRT(I)    
                    WIDTH(2)=SLRTWD(I) 
                  ENDIF
                  IF(J.NE.3) THEN
                    B(3) = WMXUSS(L,1)
                  ELSE
                    DO MIX=1,2
                      MIXING(2+MIX) = LMIXSS(1,MIX)
                      B(2+MIX) = WMXUSS(L,1)*LMIXSS(1,MIX)-
     &                WMXUSS(L,2)*LMIXSS(2,MIX)*MLP(J)*MCHAR(1)
                    ENDDO
                    RESM(4)= MSLRT(J)    
                    WIDTH(4)=SLRTWD(J) 
                  ENDIF
                  IF(I.GT.J) THEN
                    LAMCOL = G**2*LAMDA1(I,J,K)**2
                  ELSE
                    LAMCOL = ZERO
                  ENDIF
                  POUT1 = 9+2*I 
                  POUT2 = 9+2*J
                  POUT3 = -10-2*K
C--LLE +ve lepton, +ve lepton, -ve lepton 
                ELSEIF(CHANEL.EQ.2) THEN
                  M(1) = MLP(I)
                  M(2) = MLP(J)
                  M(3) = MLP(K)
                  M(4) = ABS(CHARM(L))
                  RESM(1) = MSLLT(2*I)
                  RESM(3) = MSLLT(2*J)
                  WIDTH(1)= SLLTWD(2*I-1)
                  WIDTH(3)= SLLTWD(2*J-1)
                  A(1) = -MLP(I)*WMXUSS(L,2)*MCHAR(1)
                  B(1) = WMXVSS(L,1)
                  A(3) = -MLP(J)*WMXUSS(L,2)*MCHAR(1)
                  B(3) = WMXVSS(L,1)
                  IF(CMSIGN(L).LT.0) THEN
                    B(1) = -B(1)
                    B(3) = -B(3)
                  ENDIF
                  IF(I.GT.J) THEN
                    LAMCOL = G**2*LAMDA1(I,J,K)**2
                  ELSE
                    LAMCOL = ZERO
                  ENDIF
                  POUT1 = -10-2*I 
                  POUT2 = -10-2*J
                  POUT3 =  10+2*K
C--LQD to charged lepton, dbar, down
                ELSEIF(CHANEL.EQ.3) THEN
                  M(1) = MLP(I)
                  M(2) = MQU(2*J-1)
                  M(3) = MQU(2*K-1)
                  M(4) = ABS(CHARM(L))
                  RESM(1) = MSLLT(2*I)
                  RESM(3) = MSQLT(2*J)
                  WIDTH(1)= SLLTWD(2*I-1)
                  WIDTH(3)= SQLTWD(2*J)
                  A(1) = -MLP(I)*WMXUSS(L,2)*MCHAR(1)
                  B(1) =  WMXVSS(L,1)
                  IF(J.NE.3) THEN
                    A(3) = -MQU(2*J-1)*WMXUSS(L,2)*MCHAR(1)
                    B(3) = WMXVSS(L,1)
                  ELSE
                    DO MIX=1,2
                      MIXING(2+MIX) = TMIXSS(1,MIX)
                      A(2+MIX) = -MQU(2*J-1)*WMXUSS(L,2)*MCHAR(1)
     &                           *TMIXSS(1,MIX)
                      B(2+MIX) = WMXVSS(L,1)*TMIXSS(1,MIX)-
     &                WMXVSS(L,2)*TMIXSS(2,MIX)*MQU(2*J)*MCHAR(2)
                    ENDDO
                    RESM(4)= MSQRT(2*J)    
                    WIDTH(4)=SQRTWD(2*J) 
                  ENDIF
                  IF(CMSIGN(L).LT.0) THEN
                    B(1) = -B(1)
                    B(3) = -B(3)
                    B(4) = -B(4)
                  ENDIF
                  LAMCOL = 3*G**2*LAMDA2(I,J,K)**2
                  POUT1 = -10-2*I 
                  POUT2 = 1-2*J
                  POUT3 = 2*K-1
                  IF(J.EQ.1) POUT2 = -2 
                  IF(K.EQ.1) POUT3 = 2
C--LQD to neutrino, up, dbar
                ELSEIF(CHANEL.EQ.4) THEN
                  M(1) = ZERO
                  M(2) = MQU(2*J)
                  M(3) = MQU(2*K-1)
                  M(4) = ABS(CHARM(L))
                  RESM(1) = MSLLT(2*I-1)
                  RESM(3) = MSQLT(2*J-1)
                  WIDTH(1)= SLLTWD(2*I)
                  WIDTH(3)= SQLTWD(2*J-1)
                  IF(I.NE.3) THEN
                    B(1) = WMXUSS(L,1)
                  ELSE
                    DO MIX=1,2
                      MIXING(MIX) = LMIXSS(1,MIX)
                      B(MIX) = WMXUSS(L,1)*LMIXSS(1,MIX)-
     &                WMXUSS(L,2)*LMIXSS(2,MIX)*MLP(I)*MCHAR(1)
                    ENDDO
                    RESM(2)= MSLRT(I)    
                    WIDTH(2)=SLRTWD(I) 
                  ENDIF
                  IF(J.NE.3) THEN
                    A(3) = -MQU(2*J)*WMXVSS(L,2)*MCHAR(2)
                    B(3) = WMXUSS(L,1)
                  ELSE
                    DO MIX=1,2
                      MIXING(2+MIX) = BMIXSS(1,MIX)
                      A(2+MIX) = -MQU(2*J)*WMXVSS(L,2)*BMIXSS(1,MIX)*
     &                         MCHAR(2)
                      B(2+MIX) = WMXUSS(L,1)*BMIXSS(1,MIX)-
     &            WMXUSS(L,2)*BMIXSS(2,MIX)*MQU(2*J-1)*MCHAR(1)
                    ENDDO
                    RESM(4)= MSQRT(2*J-1)    
                    WIDTH(4)=SQRTWD(2*J-1) 
                  ENDIF
                  IF(CMSIGN(L).LT.0) THEN
                    DO N=1,4
                      A(N) = -A(N)
                    ENDDO
                  ENDIF
                  LAMCOL = 3*G**2*LAMDA2(I,J,K)**2
                  POUT1 = 9+2*I 
                  POUT2 = 2*J
                  POUT3 = 1-2*K
                  IF(J.EQ.1) POUT2 = 1 
                  IF(K.EQ.1) POUT3 = -2
                ENDIF
C--Decide whether to remove diagrams
                CHARWD =ZERO
                IF(M(4).LT.(M(1)+M(2)+M(3))) GOTO 50
                DO N=1,4
                  CRSTRM(N) = .FALSE.
                ENDDO
                CRSTRM(1) = (M(4).GT.(M(1)+ABS(RESM(1))).
     &                       AND.ABS(RESM(1)).GT.(M(2)+M(3))) 
                CRSTRM(2) = (M(4).GT.(M(1)+ABS(RESM(2))).
     &                       AND.ABS(RESM(2)).GT.(M(2)+M(3))) 
     &                       .OR.(ABS(RESM(2)).LT.EPS)
                CRSTRM(3) = (M(4).GT.(M(2)+ABS(RESM(3))).
     &                       AND.ABS(RESM(3)).GT.(M(1)+M(3)))
                CRSTRM(4) = (M(4).GT.(M(2)+ABS(RESM(4))).
     &                       AND.ABS(RESM(4)).GT.(M(1)+M(3)))
     &                       .OR.(ABS(RESM(4)).LT.EPS)
C--Calculation of the rate
                CHARWD = ZERO
                IF((CRSTRM(1).AND.CRSTRM(2).AND.CRSTRM(3).AND.
     &          CRSTRM(4)).OR.LAMCOL.LT.EPS) GOTO 50
C--first the diagram squared pieces
                IF(.NOT.CRSTRM(1)) CHARWD = CHARWD+ MIXING(1)**2*
     &          RPINF1(M(2),M(3),M(1),M(4),WIDTH(1),RESM(1),A(1),B(1))
                IF(.NOT.CRSTRM(2)) CHARWD = CHARWD+MIXING(2)**2*
     &          RPINF1(M(2),M(3),M(1),M(4),WIDTH(2),RESM(2),A(2),B(2))
                IF(.NOT.CRSTRM(3)) CHARWD = CHARWD+MIXING(3)**2*
     &          RPINF1(M(1),M(3),M(2),M(4),WIDTH(3),RESM(3),A(3),B(3))
                IF(.NOT.CRSTRM(4)) CHARWD = CHARWD+ MIXING(4)**2*
     &          RPINF1(M(1),M(3),M(2),M(4),WIDTH(4),RESM(4),A(4),B(4))
C--now for the light/heavy interference due left/right mixing
                IF(.NOT.CRSTRM(1).AND..NOT.CRSTRM(2)) CHARWD=CHARWD+
     &          MIXING(1)*MIXING(2)*RPINF2(M(2),M(3),M(1),M(4),WIDTH(1),
     &           WIDTH(2),RESM(1),RESM(2),A(1),A(2),B(1),B(2),1)
                IF(.NOT.CRSTRM(3).AND..NOT.CRSTRM(4)) CHARWD=CHARWD+
     &          MIXING(3)*MIXING(4)*RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),
     &          WIDTH(4),RESM(3),RESM(4),A(3),A(4),B(3),B(4),1)
C--now for the true interference terms
                IF(.NOT.CRSTRM(1)) THEN
                  IF(.NOT.CRSTRM(3)) CHARWD=CHARWD+MIXING(1)*MIXING(3)*
     &            RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(1),RESM(3),
     &            RESM(1),A(1),A(3),B(1),B(3),2)
                  IF(.NOT.CRSTRM(4)) CHARWD=CHARWD+MIXING(1)*MIXING(4)*
     &            RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(1),RESM(4),
     &            RESM(1),A(1),A(4),B(1),B(4),2)
                ENDIF
                IF(.NOT.CRSTRM(2)) THEN
                  IF(.NOT.CRSTRM(3)) CHARWD=CHARWD+MIXING(2)*MIXING(3)*
     &            RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(2),RESM(3),
     &            RESM(2),A(2),A(3),B(2),B(3),2)
                  IF(.NOT.CRSTRM(4)) CHARWD=CHARWD+MIXING(2)*MIXING(4)*
     &            RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(2),RESM(4),
     &            RESM(2),A(2),A(4),B(2),B(4),2)
                ENDIF
C--final factors
                CHARWD = LAMCOL*CHARWD/((2*PI)**3*32*M(4)**3)
 50             IF(CHARWD.GT.EPS) CALL RPMODA(CHARWD,29+10*L,
     &                                          POUT1,POUT2,POUT3)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C--Finally the BV chargino
C--first the decay to two up type and one down type quark
      DO L=1,2
        DO I=1,3
          DO J=1,3
            DO K=1,3
              DO MIX=1,2
                MIXING(MIX) = LAMDA3(J,I,K)
                MIXING(2+MIX)= LAMDA3(I,J,K)
              ENDDO
              M(1) = MQU(2*I)
              M(2) = MQU(2*J)
              M(3) = MQU(2*K-1)
              M(4) = ABS(CHARM(L))
              RESM(1) = MSQRT(2*I-1)
              RESM(3) = MSQRT(2*J-1)
              WIDTH(1)= SQRTWD(2*I-1)
              WIDTH(3)= SQRTWD(2*J-1)
              RESM(2) =ZERO
              WIDTH(2)=ZERO
              RESM(4) =ZERO
              WIDTH(4)=ZERO
              IF(I.NE.3) THEN
                A(1)=-WMXUSS(L,2)*MCHAR(1)*MQU(2*I-1)
                B(1)= ZERO
                MIXING(2) = ZERO
              ELSE
                DO MIX=1,2
                  MIXING(MIX) = BMIXSS(2,MIX)*MIXING(MIX)
                  A(MIX) = WMXUSS(L,1)*BMIXSS(1,MIX)-
     &         WMXUSS(L,2)*MCHAR(1)*MQU(2*I-1)*BMIXSS(2,MIX)
                  B(MIX) = -MQU(2*I)*WMXVSS(L,2)*BMIXSS(1,MIX)*
     &                     MCHAR(2)
                ENDDO
                RESM(1) = MSQLT(2*I-1)
                RESM(2) = MSQRT(2*I-1)
                WIDTH(1)= SQLTWD(2*I-1)
                WIDTH(2)= SQRTWD(2*I-1)
              ENDIF
              IF(J.NE.3) THEN
                A(3) = -MQU(2*J-1)*WMXUSS(L,2)*MCHAR(1)
                B(3) = ZERO
                MIXING(4) = ZERO
              ELSE
                DO MIX=1,2
                  MIXING(2+MIX) = BMIXSS(2,MIX)*MIXING(MIX+2)
                  A(2+MIX) = WMXUSS(L,1)*BMIXSS(1,MIX)-
     &         WMXUSS(L,2)*MCHAR(1)*MQU(2*J-1)*BMIXSS(2,MIX)
                  B(2+MIX) = -MQU(2*J)*WMXVSS(L,2)*BMIXSS(1,MIX)*
     &                     MCHAR(2)
                ENDDO 
                RESM(3) = MSQLT(2*J-1)
                RESM(4) = MSQRT(2*J-1)
                WIDTH(3)= SQLTWD(2*J-1)
                WIDTH(4)= SQRTWD(2*J-1)
              ENDIF
              IF(CMSIGN(L).LT.0) THEN
                DO N=1,4
                  B(N) = -B(N)
                ENDDO
              ENDIF
              IF(I.LE.J) THEN
                LAMCOL = 6*G**2
              ELSE
                LAMCOL = ZERO
              ENDIF
              POUT1 = 2*I
              POUT2 = 2*J
              POUT3 = 2*K-1
              IF(I.EQ.1) POUT1 = 1
              IF(J.EQ.1) POUT2 = 1 
              IF(K.EQ.1) POUT3 = 2
C--Decide whether to remove diagrams
              DO N=1,4
                CRSTRM(N) = .FALSE.
              ENDDO
              CHARWD = ZERO
              IF(M(4).LT.(M(1)+M(2)+M(3)).OR.LAMCOL.LT.EPS) GOTO 60
              CRSTRM(1) = (M(4).GT.(M(1)+ABS(RESM(1))).
     &                     AND.ABS(RESM(1)).GT.(M(2)+M(3))) 
              CRSTRM(2) = (M(4).GT.(M(1)+ABS(RESM(2))).
     &                     AND.ABS(RESM(2)).GT.(M(2)+M(3))) 
     &                     .OR.(ABS(RESM(2)).LT.EPS)
              CRSTRM(3) = (M(4).GT.(M(2)+ABS(RESM(3))).
     &                     AND.ABS(RESM(3)).GT.(M(1)+M(3)))
              CRSTRM(4) = (M(4).GT.(M(2)+ABS(RESM(4))).
     &                     AND.ABS(RESM(4)).GT.(M(1)+M(3)))
     &                     .OR.(ABS(RESM(4)).LT.EPS)
C--Calculation of the rate
              CHARWD = ZERO
              IF((CRSTRM(1).AND.CRSTRM(2).AND.CRSTRM(3).AND.
     &        CRSTRM(4))) GOTO 60
C--first the diagram squared pieces
              IF(.NOT.CRSTRM(1)) CHARWD = CHARWD+ MIXING(1)**2*
     &        RPINF1(M(2),M(3),M(1),M(4),WIDTH(1),RESM(1),A(1),B(1))
              IF(.NOT.CRSTRM(2)) CHARWD = CHARWD+MIXING(2)**2*
     &        RPINF1(M(2),M(3),M(1),M(4),WIDTH(2),RESM(2),A(2),B(2))
              IF(.NOT.CRSTRM(3)) CHARWD = CHARWD+MIXING(3)**2*
     &        RPINF1(M(1),M(3),M(2),M(4),WIDTH(3),RESM(3),A(3),B(3))
              IF(.NOT.CRSTRM(4)) CHARWD = CHARWD+ MIXING(4)**2*
     &        RPINF1(M(1),M(3),M(2),M(4),WIDTH(4),RESM(4),A(4),B(4))
C--now for the light/heavy interference due left/right mixing
              IF(.NOT.CRSTRM(1).AND..NOT.CRSTRM(2)) CHARWD=CHARWD+
     &        MIXING(1)*MIXING(2)*RPINF2(M(2),M(3),M(1),M(4),WIDTH(1),
     &        WIDTH(2),RESM(1),RESM(2),A(1),A(2),B(1),B(2),1)
              IF(.NOT.CRSTRM(3).AND..NOT.CRSTRM(4)) CHARWD=CHARWD+
     &        MIXING(3)*MIXING(4)*RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),
     &        WIDTH(4),RESM(3),RESM(4),A(3),A(4),B(3),B(4),1)
C--now for the true interference terms
              IF(.NOT.CRSTRM(1)) THEN
                IF(.NOT.CRSTRM(3)) CHARWD=CHARWD+MIXING(1)*MIXING(3)*
     &          RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(1),RESM(3),
     &          RESM(1),A(1),A(3),B(1),B(3),2)
                IF(.NOT.CRSTRM(4)) CHARWD=CHARWD+MIXING(1)*MIXING(4)*
     &          RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(1),RESM(4),
     &          RESM(1),A(1),A(4),B(1),B(4),2)
              ENDIF
              IF(.NOT.CRSTRM(2)) THEN
                IF(.NOT.CRSTRM(3)) CHARWD=CHARWD+MIXING(2)*MIXING(3)*
     &          RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(2),RESM(3),
     &          RESM(2),A(2),A(3),B(2),B(3),2)
                IF(.NOT.CRSTRM(4)) CHARWD=CHARWD+MIXING(2)*MIXING(4)*
     &          RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(2),RESM(4),
     &          RESM(2),A(2),A(4),B(2),B(4),2)
              ENDIF
              CHARWD = LAMCOL*CHARWD/((2*PI)**3*32*M(4)**3)
C--Identical particle symmetry factor
              IF(I.EQ.J) CHARWD = CHARWD/2.0
C--Output the rate and particles to decay tables
 60           IF(CHARWD.GT.EPS) CALL RPMODA(CHARWD,29+10*L,POUT1,
     &                                                    POUT2,POUT3)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C--Now the decay to three down quarks
      DO L=1,2
        DO I=1,3
          DO J=1,3
            DO K=1,3
              DO MIX=1,2
                MIXING(MIX)   = LAMDA3(I,J,K)
                MIXING(MIX+2) = LAMDA3(J,K,I)
                MIXING(MIX+4) = LAMDA3(K,I,J)
              ENDDO   
              DO MIX=1,3
                RESM(2*MIX-1)   = ZERO
                RESM(2*MIX)     = ZERO
                WIDTH(2*MIX-1)  = ZERO
                WIDTH(2*MIX)    = ZERO 
              ENDDO
              CHARWD = ZERO
              M(1) = MQU(2*I-1)
              M(2) = MQU(2*J-1)
              M(3) = MQU(2*K-1)
              M(4) = ABS(CHARM(L))
              RESM(1) = MSQRT(2*I)
              RESM(3) = MSQRT(2*J)
              RESM(5) = MSQRT(2*K)
              WIDTH(1)= SQRTWD(2*I)
              WIDTH(3)= SQRTWD(2*J)
              WIDTH(5)= SQRTWD(2*K)
              IF(I.NE.3) THEN
                A(1) = -WMXVSS(L,2)*MQU(2*I)*MCHAR(2)
                B(1) = ZERO
                MIXING(2) = ZERO
              ELSE
                DO MIX=1,2
                  MIXING(MIX) = TMIXSS(2,MIX)*MIXING(MIX)
                  A(MIX) =WMXVSS(L,1)*TMIXSS(1,MIX)-
     &     WMXVSS(L,2)*MQU(2*I)*MCHAR(2)*TMIXSS(2,MIX)
                  B(MIX) =-MQU(2*I-1)*TMIXSS(1,MIX)*WMXUSS(L,2)*
     &                     MCHAR(1)
                ENDDO
C--bug fix 09/04/02 by P.R.
                RESM(1) = MSQLT(2*I)
                RESM(2) = MSQRT(2*I)
                WIDTH(1)= SQLTWD(2*I)
                WIDTH(2)= SQRTWD(2*I)
              ENDIF
              IF(J.NE.3) THEN
                A(3) = -WMXVSS(L,2)*MQU(2*J)*MCHAR(2)
                B(3) = ZERO 
                MIXING(4) = ZERO
              ELSE
                DO MIX=1,2
                  MIXING(MIX+2) = TMIXSS(2,MIX)*MIXING(MIX+2)
                  A(MIX+2) =WMXVSS(L,1)*TMIXSS(1,MIX)-
     &     WMXVSS(L,2)*MQU(2*J)*MCHAR(2)*TMIXSS(2,MIX)
                  B(MIX+2) =-MQU(2*I-1)*TMIXSS(1,MIX)*WMXUSS(L,2)*
     &                     MCHAR(1)
                ENDDO
C--bug fix 09/04/02 by P.R.
                RESM(3) = MSQLT(2*J)
                RESM(4) = MSQRT(2*J)
                WIDTH(3)= SQLTWD(2*J)
                WIDTH(4)= SQRTWD(2*J)
              ENDIF
              IF(K.NE.3) THEN
                A(5) = -WMXVSS(L,2)*MQU(2*K)*MCHAR(2)
                B(5) = ZERO
C--bug fix 09/04/02 by P.R.
                MIXING(6) = ZERO
              ELSE
                DO MIX=1,2
                  MIXING(MIX+4) = TMIXSS(2,MIX)*MIXING(MIX+4)
                  A(MIX+4) =WMXVSS(L,1)*TMIXSS(1,MIX)-
     &     WMXVSS(L,2)*MQU(2*K)*MCHAR(2)*TMIXSS(2,MIX)
                  B(MIX+4) =-MQU(2*K-1)*TMIXSS(1,MIX)*WMXUSS(L,2)*
     &                     MCHAR(1)
                ENDDO
C--bug fix 09/04/02 by P.R.
                RESM(5) = MSQLT(2*K)
                RESM(6) = MSQRT(2*K)
                WIDTH(5)= SQLTWD(2*K)
                WIDTH(6)= SQRTWD(2*K)
              ENDIF
              IF(CMSIGN(L).LT.0) THEN
                DO N=1,6
                  A(N) = -A(N)
                ENDDO
              ENDIF
              IF(K.LE.J.AND.J.LE.I) THEN
                LAMCOL = 6*G**2
              ELSE
                LAMCOL = ZERO
              ENDIF
              POUT1 = 1-2*I
              POUT2 = 1-2*J
              POUT3 = 1-2*K
              IF(I.EQ.1) POUT1 = -2
              IF(J.EQ.1) POUT2 = -2
              IF(K.EQ.1) POUT3 = -2
C--Decide whether to remove diagrams
              DO N=1,6
                CRSTRM(N) = .FALSE.
              ENDDO
              IF(M(4).LT.(M(1)+M(2)+M(3))) GOTO 70
              CRSTRM(1) = (M(4).GT.(M(1)+ABS(RESM(1))).
     &                     AND.ABS(RESM(1)).GT.(M(2)+M(3))) 
              CRSTRM(2) = (M(4).GT.(M(1)+ABS(RESM(2))).
     &                     AND.ABS(RESM(2)).GT.(M(2)+M(3))) 
     &                     .OR.(ABS(RESM(2)).LT.EPS)
              CRSTRM(3) = (M(4).GT.(M(2)+ABS(RESM(3))).
     &                     AND.ABS(RESM(3)).GT.(M(1)+M(3)))
              CRSTRM(4) = (M(4).GT.(M(2)+ABS(RESM(4))).
     &                     AND.ABS(RESM(4)).GT.(M(1)+M(3)))
     &                     .OR.(ABS(RESM(4)).LT.EPS)
              CRSTRM(5) = (M(4).GT.(M(3)+ABS(RESM(5))).
     &                     AND.ABS(RESM(5)).GT.(M(1)+M(2)))
              CRSTRM(6) = (M(4).GT.(M(3)+ABS(RESM(6))).
     &                     AND.ABS(RESM(6)).GT.(M(1)+M(2)))
     &                     .OR.(ABS(RESM(6)).LT.EPS)
C--Calculation of the rate
              CHARWD = ZERO
              IF((CRSTRM(1).AND.CRSTRM(2).AND.CRSTRM(3).AND.
     &        CRSTRM(4).AND.CRSTRM(5).AND.CRSTRM(6)).
     &        OR.LAMCOL.LT.EPS) GOTO 70
C--first the diagram squared pieces
              IF(.NOT.CRSTRM(1)) CHARWD = CHARWD+ MIXING(1)**2*
     &        RPINF1(M(2),M(3),M(1),M(4),WIDTH(1),RESM(1),A(1),B(1))
              IF(.NOT.CRSTRM(2)) CHARWD = CHARWD+MIXING(2)**2*
     &        RPINF1(M(2),M(3),M(1),M(4),WIDTH(2),RESM(2),A(2),B(2))
              IF(.NOT.CRSTRM(3)) CHARWD = CHARWD+MIXING(3)**2*
     &        RPINF1(M(1),M(3),M(2),M(4),WIDTH(3),RESM(3),A(3),B(3))
              IF(.NOT.CRSTRM(4)) CHARWD = CHARWD+ MIXING(4)**2*
     &        RPINF1(M(1),M(3),M(2),M(4),WIDTH(4),RESM(4),A(4),B(4))
              IF(.NOT.CRSTRM(5)) CHARWD = CHARWD+MIXING(5)**2*
     &        RPINF1(M(1),M(2),M(3),M(4),WIDTH(5),RESM(5),A(5),B(5))
              IF(.NOT.CRSTRM(6)) CHARWD = CHARWD+MIXING(6)**2*
     &        RPINF1(M(1),M(2),M(3),M(4),WIDTH(6),RESM(6),A(6),B(6))
C--now for the light/heavy interference due left/right mixing
              IF(.NOT.CRSTRM(1).AND..NOT.CRSTRM(2)) CHARWD=CHARWD+
     &        MIXING(1)*MIXING(2)*RPINF2(M(2),M(3),M(1),M(4),WIDTH(1),
     &        WIDTH(2),RESM(1),RESM(2),A(1),A(2),B(1),B(2),1)
              IF(.NOT.CRSTRM(3).AND..NOT.CRSTRM(4)) CHARWD=CHARWD+
     &        MIXING(3)*MIXING(4)*RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),
     &        WIDTH(4),RESM(3),RESM(4),A(3),A(4),B(3),B(4),1)
              IF(.NOT.CRSTRM(5).AND..NOT.CRSTRM(6)) CHARWD=CHARWD+
     &        MIXING(5)*MIXING(6)*RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),
     &        WIDTH(6),RESM(5),RESM(6),A(5),A(6),B(5),B(6),1)
C--now for the true interference terms
              IF(.NOT.CRSTRM(1)) THEN
                IF(.NOT.CRSTRM(3)) CHARWD=CHARWD-MIXING(1)*MIXING(3)*
     &          RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(1),RESM(3),
     &          RESM(1),A(1),A(3),B(1),B(3),2)
                IF(.NOT.CRSTRM(4)) CHARWD=CHARWD-MIXING(1)*MIXING(4)*
     &          RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(1),RESM(4),
     &          RESM(1),A(1),A(4),B(1),B(4),2)
                IF(.NOT.CRSTRM(5)) CHARWD=CHARWD-MIXING(1)*MIXING(5)*
     &          RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),WIDTH(1),RESM(5),
     &          RESM(1),A(1),A(5),B(1),B(5),2)
                IF(.NOT.CRSTRM(6)) CHARWD=CHARWD-MIXING(1)*MIXING(6)*
     &          RPINF2(M(1),M(2),M(3),M(4),WIDTH(6),WIDTH(1),RESM(6),
     &          RESM(1),A(1),A(6),B(1),B(6),2)
              ENDIF
              IF(.NOT.CRSTRM(2)) THEN
                IF(.NOT.CRSTRM(3)) CHARWD=CHARWD-MIXING(2)*MIXING(3)*
     &          RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(2),RESM(3),
     &          RESM(2),A(2),A(3),B(2),B(3),2)
                IF(.NOT.CRSTRM(4)) CHARWD=CHARWD-MIXING(2)*MIXING(4)*
     &          RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(2),RESM(4),
     &          RESM(2),A(2),A(4),B(2),B(4),2)
                IF(.NOT.CRSTRM(5)) CHARWD=CHARWD-MIXING(2)*MIXING(5)*
     &          RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),WIDTH(2),RESM(5),
     &          RESM(2),A(2),A(5),B(2),B(5),2)
                IF(.NOT.CRSTRM(6)) CHARWD=CHARWD-MIXING(2)*MIXING(6)*
     &          RPINF2(M(1),M(2),M(3),M(4),WIDTH(6),WIDTH(2),RESM(6),
     &          RESM(2),A(2),A(6),B(2),B(6),2)
              ENDIF
              IF(.NOT.CRSTRM(3)) THEN
                IF(.NOT.CRSTRM(5)) CHARWD=CHARWD-MIXING(3)*MIXING(5)*
     &          RPINF2(M(2),M(1),M(3),M(4),WIDTH(5),WIDTH(3),RESM(5),
     &          RESM(3),A(3),A(5),B(3),B(5),2)
                IF(.NOT.CRSTRM(6)) CHARWD=CHARWD-MIXING(3)*MIXING(6)*
     &          RPINF2(M(2),M(1),M(3),M(4),WIDTH(6),WIDTH(3),RESM(6),
     &          RESM(3),A(3),A(6),B(3),B(6),2)
              ENDIF
              IF(.NOT.CRSTRM(4)) THEN
                IF(.NOT.CRSTRM(5)) CHARWD=CHARWD-MIXING(4)*MIXING(5)*
     &          RPINF2(M(2),M(1),M(3),M(4),WIDTH(5),WIDTH(4),RESM(5),
     &          RESM(4),A(4),A(5),B(4),B(5),2)
                IF(.NOT.CRSTRM(6)) CHARWD=CHARWD-MIXING(4)*MIXING(6)*
     &          RPINF2(M(2),M(1),M(3),M(4),WIDTH(6),WIDTH(4),RESM(6),
     &          RESM(4),A(4),A(6),B(4),B(6),2)
              ENDIF
              CHARWD = LAMCOL*CHARWD/((2*PI)**3*32*M(4)**3)
C--Identical particle symmetry factor
              IF(I.EQ.J.OR.I.EQ.K.OR.J.EQ.K) CHARWD = CHARWD/2.0
C--Output the rate and particles to decay tables
 70           IF(CHARWD.GT.EPS) CALL RPMODA(CHARWD,29+10*L,POUT1,
     &                                              POUT2,POUT3)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      END
