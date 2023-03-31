CDECK  ID>, SSM1LP.
      SUBROUTINE SSM1LP(M1,M2,IALLOW)
C-----------------------------------------------------------------------
C
C          Recalculate sparticle masses including
C          radiative corrections
C          from T. Krupovnickas and H. Baer
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
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
      COMPLEX*16 SSB0,SSB1
      DOUBLE PRECISION SSA0
      REAL PI,SR2,GP,G1,G2,G3,CS2THW,CTHW,AMTLSQ,AMTRSQ,XX,
     $APD,ADMBC,BETA,SINB,COSB,HIGFRZ,
     $PITLTL,PITRTR,PITLTR,PIBLBL,PIBLBR,PIBRBR,PILLLL,PILLLR,PILRLR,
     $PIELEL,PIERER,PINENE
      INTEGER IALLOW
      REAL SIG0L,SIG0R,SIG0S,SIGPL,SIGPR,SIGPS
     $,SSIG0L(4,4),SSIG0R(4,4),SSIG0S(4,4),SSIGPL(2,2),SSIGPR(2,2)
     $,SSIGPS(2,2)
      REAL AR(4,4),NEWAR(4,4),WR(4),TEMP,WORK(4),V,VP,M1,M2
     $,MPPTRE(2,2),MPP(2,2),MPP2(2,2),ZMIX(4,4)
      REAL TMZ1SS,TMZ2SS,TMZ3SS,TMZ4SS,TMW1SS,TMW2SS
      REAL M2P,MUP,MWP,COSBEP,SINBEP,COS2BP,SIN2BP,INOSCL,ZETA,
     $ZETAS,XM,YM,THX,THY,A,B,C,D,COS2B
      INTEGER I,J,K,IERR

      PI=4.*ATAN(1.)
      SR2=SQRT(2.)
C     We will use msbar couplings at M_Z here for now
C     to give consistency between MSSM and SUGRA solutions
      G2=SQRT(4.*PI*ALFAEM/SN2THW)
      GP=G2*SQRT(SN2THW/(1.-SN2THW))
      G1=SQRT(5./3.)*GP
      G3=SQRT(4*PI*.118)
      CS2THW=1.-SN2THW
      CTHW=SQRT(CS2THW)
      V=VUQ
      VP=VDQ
      BETA=ATAN(VUQ/VDQ)
      SINB=SIN(BETA)
      COSB=COS(BETA)
      COS2B=COS(2*BETA)
      HIGFRZ=SQRT(MAX(AMZ**2,AMTLSS*AMTRSS*SIGN(1.,AMTLSS*AMTRSS)))
C
C
C     Refill MSS() for input to self energy routines
C
      MSS(1)=AMGLSS
      MSS(2)=AMULSS
      MSS(3)=AMURSS
      MSS(4)=AMDLSS
      MSS(5)=AMDRSS
      MSS(6)=AMSLSS
      MSS(7)=AMSRSS
      MSS(8)=AMCLSS
      MSS(9)=AMCRSS
      MSS(10)=AMB1SS
      MSS(11)=AMB2SS
      MSS(12)=AMT1SS
      MSS(13)=AMT2SS
      MSS(14)=AMN1SS
      MSS(15)=AMN2SS
      MSS(16)=AMN3SS
      MSS(17)=AMELSS
      MSS(18)=AMERSS
      MSS(19)=AMMLSS
      MSS(20)=AMMRSS
      MSS(21)=AML1SS
      MSS(22)=AML2SS
      MSS(23)=AMZ1SS
      MSS(24)=AMZ2SS
      MSS(25)=AMZ3SS
      MSS(26)=AMZ4SS
      MSS(27)=AMW1SS
      MSS(28)=AMW2SS
      MSS(29)=AMHL
      MSS(30)=AMHH
      MSS(31)=AMHA
      MSS(32)=AMHC
C
C     Neutralino masses
C
      AR(1,1)=0.
      AR(1,2)=-TWOM1
      AR(1,3)=-G2*V/SR2
      AR(1,4)=SQRT(3./5.)*G1*V/SR2
      AR(2,1)=-TWOM1
      AR(2,2)=0.
      AR(2,3)=G2*VP/SR2
      AR(2,4)=-SQRT(3./5.)*G1*VP/SR2
      AR(3,1)=-G2*V/SR2
      AR(3,2)=G2*VP/SR2
      AR(3,3)=M2
      AR(3,4)=0.
      AR(4,1)=SQRT(3./5.)*G1*V/SR2
      AR(4,2)=-SQRT(3./5.)*G1*VP/SR2
      AR(4,3)=0.
      AR(4,4)=M1
      XLAM=DLOG(DBLE(HIGFRZ**2))
C     Set renormalization scale to INOSCL for -ino mass computation
      INOSCL=AMZ1SS**2
      DO I=1,4
        DO J=1,4
          SSIG0L(I,J)=SIG0L(INOSCL,5-I,5-J,G1,G2,CTHW)
          SSIG0R(I,J)=SIG0R(INOSCL,5-I,5-J,G1,G2,CTHW)
          SSIG0S(I,J)=SIG0S(INOSCL,5-I,5-J,G1,G2,CTHW)
          NEWAR(I,J)=AR(I,J)
        ENDDO
      ENDDO
      DO I=1,4
        DO J=1,4
          DO K=1,4
            NEWAR(I,J)=NEWAR(I,J)-(SSIG0R(I,K)*AR(K,J)
     $+AR(I,K)*SSIG0L(K,J)
     $+SSIG0R(J,K)*AR(K,I)+AR(J,K)*SSIG0L(K,I))/2.
          ENDDO
          NEWAR(I,J)=NEWAR(I,J)-(SSIG0S(I,J)+SSIG0S(J,I))/2.
        ENDDO
      ENDDO
      CALL EISRS1(4,4,NEWAR,WR,ZMIX,IERR,WORK)
C     Sort eigenvectors and eigenvalues according to masses
      DO I=1,3
        DO J=I+1,4
          IF (ABS(WR(I)).GT.ABS(WR(J))) THEN
            TEMP=WR(J)
            WR(J)=WR(I)
            WR(I)=TEMP
            DO K=1,4
              TEMP=ZMIX(K,J)
              ZMIX(K,J)=ZMIX(K,I)
              ZMIX(K,I)=TEMP
            ENDDO
          END IF
        ENDDO
      ENDDO
      TMZ1SS=WR(1)
C     Write over with radiatively corrected -ino mixing elements 
      DO I=1,4
        DO J=1,4
          ZMIXSS(I,J)=ZMIX(I,J)
        END DO
      END DO
C
C
      INOSCL=AMZ2SS**2
      DO I=1,4
        DO J=1,4
          SSIG0L(I,J)=SIG0L(INOSCL,5-I,5-J,G1,G2,CTHW)
          SSIG0R(I,J)=SIG0R(INOSCL,5-I,5-J,G1,G2,CTHW)
          SSIG0S(I,J)=SIG0S(INOSCL,5-I,5-J,G1,G2,CTHW)
          NEWAR(I,J)=AR(I,J)
        ENDDO
      ENDDO
      DO I=1,4
        DO J=1,4
          DO K=1,4
            NEWAR(I,J)=NEWAR(I,J)-(SSIG0R(I,K)*AR(K,J)
     $+AR(I,K)*SSIG0L(K,J)
     $+SSIG0R(J,K)*AR(K,I)+AR(J,K)*SSIG0L(K,I))/2.
          ENDDO
          NEWAR(I,J)=NEWAR(I,J)-(SSIG0S(I,J)+SSIG0S(J,I))/2.
        ENDDO
      ENDDO
      CALL EISRS1(4,4,NEWAR,WR,ZMIX,IERR,WORK)
C     Sort eigenvectors and eigenvalues according to masses
      DO I=1,3
        DO J=I+1,4
          IF (ABS(WR(I)).GT.ABS(WR(J))) THEN
            TEMP=WR(J)
            WR(J)=WR(I)
            WR(I)=TEMP
            DO K=1,4
              TEMP=ZMIX(K,J)
              ZMIX(K,J)=ZMIX(K,I)
              ZMIX(K,I)=TEMP
            ENDDO
          END IF
        ENDDO
      ENDDO
      TMZ2SS=WR(2)
C
      INOSCL=AMZ3SS**2
      DO I=1,4
        DO J=1,4
          SSIG0L(I,J)=SIG0L(INOSCL,5-I,5-J,G1,G2,CTHW)
          SSIG0R(I,J)=SIG0R(INOSCL,5-I,5-J,G1,G2,CTHW)
          SSIG0S(I,J)=SIG0S(INOSCL,5-I,5-J,G1,G2,CTHW)
          NEWAR(I,J)=AR(I,J)
        ENDDO
      ENDDO
      DO I=1,4
        DO J=1,4
          DO K=1,4
            NEWAR(I,J)=NEWAR(I,J)-(SSIG0R(I,K)*AR(K,J)
     $+AR(I,K)*SSIG0L(K,J)
     $+SSIG0R(J,K)*AR(K,I)+AR(J,K)*SSIG0L(K,I))/2.
          ENDDO
          NEWAR(I,J)=NEWAR(I,J)-(SSIG0S(I,J)+SSIG0S(J,I))/2.
        ENDDO
      ENDDO
      CALL EISRS1(4,4,NEWAR,WR,ZMIX,IERR,WORK)
C     Sort eigenvectors and eigenvalues according to masses
      DO I=1,3
        DO J=I+1,4
          IF (ABS(WR(I)).GT.ABS(WR(J))) THEN
            TEMP=WR(J)
            WR(J)=WR(I)
            WR(I)=TEMP
            DO K=1,4
              TEMP=ZMIX(K,J)
              ZMIX(K,J)=ZMIX(K,I)
              ZMIX(K,I)=TEMP
            ENDDO
          END IF
        ENDDO
      ENDDO
      TMZ3SS=WR(3)
C
      INOSCL=AMZ4SS**2
      DO I=1,4
        DO J=1,4
          SSIG0L(I,J)=SIG0L(INOSCL,5-I,5-J,G1,G2,CTHW)
          SSIG0R(I,J)=SIG0R(INOSCL,5-I,5-J,G1,G2,CTHW)
          SSIG0S(I,J)=SIG0S(INOSCL,5-I,5-J,G1,G2,CTHW)
          NEWAR(I,J)=AR(I,J)
        ENDDO
      ENDDO
      DO I=1,4
        DO J=1,4
          DO K=1,4
            NEWAR(I,J)=NEWAR(I,J)-(SSIG0R(I,K)*AR(K,J)
     $+AR(I,K)*SSIG0L(K,J)
     $+SSIG0R(J,K)*AR(K,I)+AR(J,K)*SSIG0L(K,I))/2.
          ENDDO
          NEWAR(I,J)=NEWAR(I,J)-(SSIG0S(I,J)+SSIG0S(J,I))/2.
        ENDDO
      ENDDO
      CALL EISRS1(4,4,NEWAR,WR,ZMIX,IERR,WORK)
C     Sort eigenvectors and eigenvalues according to masses
      DO I=1,3
        DO J=I+1,4
          IF (ABS(WR(I)).GT.ABS(WR(J))) THEN
            TEMP=WR(J)
            WR(J)=WR(I)
            WR(I)=TEMP
            DO K=1,4
              TEMP=ZMIX(K,J)
              ZMIX(K,J)=ZMIX(K,I)
              ZMIX(K,I)=TEMP
            ENDDO
          END IF
        ENDDO
      ENDDO
      TMZ4SS=WR(4)
C
C     Charginos
C
      INOSCL=AMW1SS**2
      DO I=1,2
        DO J=1,2
          IF(I.NE.J) THEN
            SSIGPL(I,J)=-SIGPL(INOSCL,J,I,G1,G2,CTHW)
            SSIGPR(I,J)=-SIGPR(INOSCL,J,I,G1,G2,CTHW)
            SSIGPS(I,J)=-SIGPS(INOSCL,J,I,G1,G2,CTHW)
          ELSE
            SSIGPL(I,J)=SIGPL(INOSCL,I,J,G1,G2,CTHW)
            SSIGPR(I,J)=SIGPR(INOSCL,I,J,G1,G2,CTHW)
            SSIGPS(I,J)=SIGPS(INOSCL,I,J,G1,G2,CTHW)
          ENDIF
          MPP(I,J)=0.
        ENDDO
      ENDDO
      MPPTRE(1,1)=M2
      MPPTRE(1,2)=-G2*VP
      MPPTRE(2,1)=-G2*V
      MPPTRE(2,2)=TWOM1
      DO I=1,2
        DO J=1,2
          MPP(I,J)=MPPTRE(I,J)
        ENDDO
      ENDDO
      DO I=1,2
        DO J=1,2
          DO K=1,2
            MPP(I,J)=MPP(I,J)-SSIGPR(I,K)*MPPTRE(K,J)
     $-MPPTRE(I,K)*SSIGPL(K,J)
          ENDDO
          MPP(I,J)=MPP(I,J)-SSIGPS(I,J)
        ENDDO
      ENDDO
      DO I=1,2
        DO J=1,2
          MPP2(I,J)=MPP(I,1)*MPP(J,1)+MPP(I,2)*MPP(J,2)
        ENDDO
      ENDDO
C
      M2P=MPP(1,1)
      MUP=-MPP(2,2)
      MWP=SQRT((MPP(1,2)**2+MPP(2,1)**2)/2.)
      COSBEP=-MPP(1,2)/SR2/MWP
      SINBEP=-SIGN(1.,MPP(2,1))*SQRT(1.-COSBEP**2)
      COS2BP=1.-2*SINBEP**2
      SIN2BP=2*SINBEP*COSBEP
      ZETAS=(MUP**2-M2P**2)**2
     $+4*MWP**2*(MWP**2*COS2BP**2+MUP**2+M2P**2-2*MUP*M2P*SIN2BP)
      ZETA=SQRT(ZETAS)
      XM=(MUP**2-M2P**2+2*MWP**2*COS2BP-ZETA)
     $/(2*SR2*MWP*(-M2P*COSBEP+MUP*SINBEP))
      YM=(MUP**2-M2P**2-2*MWP**2*COS2BP-ZETA)
     $/(2*SR2*MWP*(-M2P*SINBEP+MUP*COSBEP))
      IF (XM.NE.0.) THEN
        GAMMAL=ATAN(1./XM)
      ELSE
        GAMMAL=PI/2.
      END IF
      IF (YM.NE.0.) THEN
        GAMMAR=ATAN(1./YM)
      ELSE
        GAMMAR=PI/2.
      END IF
      IF (GAMMAL.LT.0.) GAMMAL=GAMMAL+PI
      IF (GAMMAR.LT.0.) GAMMAR=GAMMAR+PI
      THX=SIGN(1.,XM)
      THY=SIGN(1.,YM)
      AMW2SS=THX*THY*(COS(GAMMAR)*(M2P*COS(GAMMAL)+SR2*MWP*COSBEP*
     $SIN(GAMMAL))+SIN(GAMMAR)*(SR2*MWP*SINBEP*COS(GAMMAL)-
     $MUP*SIN(GAMMAL)))
      TMW1SS=SIGN(1.,AMW1SS)*SQRT(ABS(MIN((MPP2(1,1)+MPP2(2,2)
     $-SQRT((MPP2(1,1)-MPP2(2,2))**2+4.*MPP2(1,2)*MPP2(2,1)))/2.
     $,(MPP2(1,1)+MPP2(2,2)
     $+SQRT((MPP2(1,1)-MPP2(2,2))**2+4.*MPP2(1,2)*MPP2(2,1)))/2.)))
C
      INOSCL=AMW2SS**2
      DO I=1,2
        DO J=1,2
          IF(I.NE.J) THEN
            SSIGPL(I,J)=-SIGPL(INOSCL,J,I,G1,G2,CTHW)
            SSIGPR(I,J)=-SIGPR(INOSCL,J,I,G1,G2,CTHW)
            SSIGPS(I,J)=-SIGPS(INOSCL,J,I,G1,G2,CTHW)
          ELSE
            SSIGPL(I,J)=SIGPL(INOSCL,I,J,G1,G2,CTHW)
            SSIGPR(I,J)=SIGPR(INOSCL,I,J,G1,G2,CTHW)
            SSIGPS(I,J)=SIGPS(INOSCL,I,J,G1,G2,CTHW)
          ENDIF
          MPP(I,J)=0.
        ENDDO
      ENDDO
      MPPTRE(1,1)=M2
      MPPTRE(1,2)=-G2*VP
      MPPTRE(2,1)=-G2*V
      MPPTRE(2,2)=TWOM1
      DO I=1,2
        DO J=1,2
          MPP(I,J)=MPPTRE(I,J)
        ENDDO
      ENDDO
      DO I=1,2
        DO J=1,2
          DO K=1,2
            MPP(I,J)=MPP(I,J)-SSIGPR(I,K)*MPPTRE(K,J)
     $-MPPTRE(I,K)*SSIGPL(K,J)
          ENDDO
          MPP(I,J)=MPP(I,J)-SSIGPS(I,J)
        ENDDO
      ENDDO
      DO I=1,2
        DO J=1,2
          MPP2(I,J)=MPP(I,1)*MPP(J,1)+MPP(I,2)*MPP(J,2)
        ENDDO
      ENDDO
      TMW2SS=SIGN(1.,AMW2SS)*SQRT(ABS(MAX((MPP2(1,1)+MPP2(2,2)
     $-SQRT((MPP2(1,1)-MPP2(2,2))**2+4.*MPP2(1,2)*MPP2(2,1)))/2.
     $,(MPP2(1,1)+MPP2(2,2)
     $+SQRT((MPP2(1,1)-MPP2(2,2))**2+4.*MPP2(1,2)*MPP2(2,1)))/2.)))
C
C     Do third generation squarks and sleptons
C     mass matrix is of form (A B)
C                            (C D)
C     Top squarks
      XLAM=DLOG(DBLE(HIGFRZ**2))
      AMTLSQ=SIGN(1.,AMTLSS)*AMTLSS**2
      AMTRSQ=SIGN(1.,AMTRSS)*AMTRSS**2
      XX=(AMGLSS/AMULSS)**2
      A=AMTLSQ+MTQ**2+AMZ**2*COS2B*(.5-2*SN2THW/3.)-
     $PITLTL(AMT1SS**2,G1,G2,G3,CTHW)
      D=AMTRSQ+MTQ**2+AMZ**2*COS2B*(2*SN2THW/3.)-
     $PITRTR(AMT1SS**2,G1,G2,G3,CTHW)
      B=MTQ*(-AAT-TWOM1*COSB/SINB)-PITLTR(AMT1SS**2,G1,G2,G3,CTHW)
      C=B
      APD=A+D
      ADMBC=A*D-B*C
      AMT1SS=SQRT(ABS((APD-SQRT(ABS(APD**2-4.*ADMBC)))/2.))
      THETAT=ATAN((A-AMT1SS**2)/B)
      A=AMTLSQ+MTQ**2+AMZ**2*COS2B*(.5-2*SN2THW/3.)-
     $PITLTL(AMT2SS**2,G1,G2,G3,CTHW)
      D=AMTRSQ+MTQ**2+AMZ**2*COS2B*(2*SN2THW/3.)-
     $PITRTR(AMT2SS**2,G1,G2,G3,CTHW)
      B=MTQ*(-AAT-TWOM1*COSB/SINB)-PITLTR(AMT2SS**2,G1,G2,G3,CTHW)
      C=B
      APD=A+D
      ADMBC=A*D-B*C
      AMT2SS=SQRT(ABS((APD+SQRT(ABS(APD**2-4.*ADMBC)))/2.))
C
      A=AMBLSS**2+MBQ**2+AMZ**2*COS2B*(-.5+SN2THW/3.)-
     $PIBLBL(AMB1SS**2,G1,G2,G3,CTHW)
      D=AMBRSS**2+MBQ**2+AMZ**2*COS2B*(-SN2THW/3.)-
     $PIBRBR(AMB1SS**2,G1,G2,G3,CTHW)
      B=MBQ*(-AAB-TWOM1*SINB/COSB)-PIBLBR(AMB1SS**2,G1,G2,G3,CTHW)
      C=B
      APD=A+D
      ADMBC=A*D-B*C
      AMB1SS=SQRT(ABS((APD-SQRT(ABS(APD**2-4.*ADMBC)))/2.))
      THETAB=ATAN((A-AMB1SS**2)/B)
      A=AMBLSS**2+MBQ**2+AMZ**2*COS2B*(-.5+SN2THW/3.)-
     $PIBLBL(AMB2SS**2,G1,G2,G3,CTHW)
      D=AMBRSS**2+MBQ**2+AMZ**2*COS2B*(-SN2THW/3.)-
     $PIBRBR(AMB2SS**2,G1,G2,G3,CTHW)
      B=MBQ*(-AAB-TWOM1*SINB/COSB)-PIBLBR(AMB2SS**2,G1,G2,G3,CTHW)
      C=B
      APD=A+D
      ADMBC=A*D-B*C
      AMB2SS=SQRT(ABS((APD+SQRT(ABS(APD**2-4.*ADMBC)))/2.))
C
      A=AMLLSS**2+MLQ**2+AMZ**2*COS2B*(-.5+SN2THW)-
     $PILLLL(AML1SS**2,G1,G2,G3,CTHW)
      D=AMLRSS**2+MLQ**2+AMZ**2*COS2B*(-SN2THW)-
     $PILRLR(AML1SS**2,G1,G2,G3,CTHW)
      B=MLQ*(-AAL-TWOM1*SINB/COSB)-PILLLR(AML1SS**2,G1,G2,G3,CTHW)
      C=B
      APD=A+D
      ADMBC=A*D-B*C
      AML1SS=SQRT(ABS((APD-SQRT(ABS(APD**2-4.*ADMBC)))/2.))
      THETAL=ATAN((A-AML1SS**2)/B)
      A=AMLLSS**2+MLQ**2+AMZ**2*COS2B*(-.5+SN2THW)-
     $PILLLL(AML2SS**2,G1,G2,G3,CTHW)
      D=AMLRSS**2+MLQ**2+AMZ**2*COS2B*(-SN2THW)-
     $PILRLR(AML2SS**2,G1,G2,G3,CTHW)
      B=MLQ*(-AAL-TWOM1*SINB/COSB)-PILLLR(AML2SS**2,G1,G2,G3,CTHW)
      C=B
      APD=A+D
      ADMBC=A*D-B*C
      AML2SS=SQRT(ABS((APD+SQRT(ABS(APD**2-4.*ADMBC)))/2.))
C
      XLAM=DLOG(DBLE(AMULSS**2))
      XX=(AMGLSS/AMULSS)**2
      AMULSS=SQRT(AMULSS**2*(1.+G3**2/6./PI**2
     $*(2.*SSB1(AMULSS**2,AMULSS,0.)+SSA0(AMGLSS)/AMULSS**2
     $-(1.-XX)*SSB0(AMULSS**2,AMGLSS,0.))))
C
      XLAM=DLOG(DBLE(AMURSS**2))
      XX=(AMGLSS/AMURSS)**2
      AMURSS=SQRT(AMURSS**2*(1.+G3**2/6./PI**2
     $*(2.*SSB1(AMURSS**2,AMURSS,0.)+SSA0(AMGLSS)/AMURSS**2
     $-(1.-XX)*SSB0(AMURSS**2,AMGLSS,0.))))
C
      XLAM=DLOG(DBLE(AMDLSS**2))
      XX=(AMGLSS/AMDLSS)**2
      AMDLSS=SQRT(AMDLSS**2*(1.+G3**2/6./PI**2
     $*(2.*SSB1(AMDLSS**2,AMDLSS,0.)+SSA0(AMGLSS)/AMDLSS**2
     $-(1.-XX)*SSB0(AMDLSS**2,AMGLSS,0.))))
C
      XLAM=DLOG(DBLE(AMDRSS**2))
      XX=(AMGLSS/AMDRSS)**2
      AMDRSS=SQRT(AMDRSS**2*(1.+G3**2/6./PI**2
     $*(2.*SSB1(AMDRSS**2,AMDRSS,0.)+SSA0(AMGLSS)/AMDRSS**2
     $-(1.-XX)*SSB0(AMDRSS**2,AMGLSS,0.))))
C
      XLAM=DLOG(DBLE(AMCLSS**2))
      XX=(AMGLSS/AMCLSS)**2
      AMCLSS=SQRT(AMCLSS**2*(1.+G3**2/6./PI**2
     $*(2.*SSB1(AMCLSS**2,AMCLSS,0.)+SSA0(AMGLSS)/AMCLSS**2
     $-(1.-XX)*SSB0(AMCLSS**2,AMGLSS,0.))))
C
      XLAM=DLOG(DBLE(AMCRSS**2))
      XX=(AMGLSS/AMCRSS)**2
      AMCRSS=SQRT(AMCRSS**2*(1.+G3**2/6./PI**2
     $*(2.*SSB1(AMCRSS**2,AMCRSS,0.)+SSA0(AMGLSS)/AMCRSS**2
     $-(1.-XX)*SSB0(AMCRSS**2,AMGLSS,0.))))
C
      XLAM=DLOG(DBLE(AMSLSS**2))
      XX=(AMGLSS/AMSLSS)**2
      AMSLSS=SQRT(AMSLSS**2*(1.+G3**2/6./PI**2
     $*(2.*SSB1(AMSLSS**2,AMSLSS,0.)+SSA0(AMGLSS)/AMSLSS**2
     $-(1.-XX)*SSB0(AMSLSS**2,AMGLSS,0.))))
C
      XLAM=DLOG(DBLE(AMSRSS**2))
      XX=(AMGLSS/AMSRSS)**2
      AMSRSS=SQRT(AMSRSS**2*(1.+G3**2/6./PI**2
     $*(2.*SSB1(AMSRSS**2,AMSRSS,0.)+SSA0(AMGLSS)/AMSRSS**2
     $-(1.-XX)*SSB0(AMSRSS**2,AMGLSS,0.))))
C
      XLAM=DLOG(DBLE(AMELSS**2))
      AMELSS=SQRT(AMELSS**2-PIELEL(AMELSS**2,G1,G2,G3,CTHW))     
      XLAM=DLOG(DBLE(AMERSS**2))
      AMERSS=SQRT(AMERSS**2-PIERER(AMERSS**2,G1,G2,G3,CTHW))
      XLAM=DLOG(DBLE(AMN1SS**2))
      AMN1SS=SQRT(AMN1SS**2-PINENE(AMN1SS**2,G1,G2,G3,CTHW))
      XLAM=DLOG(DBLE(AMMLSS**2))
      AMMLSS=SQRT(AMMLSS**2-PIELEL(AMMLSS**2,G1,G2,G3,CTHW))     
      XLAM=DLOG(DBLE(AMMRSS**2))
      AMMRSS=SQRT(AMMRSS**2-PIERER(AMMRSS**2,G1,G2,G3,CTHW))
      XLAM=DLOG(DBLE(AMN2SS**2))
      AMN2SS=SQRT(AMN2SS**2-PINENE(AMN2SS**2,G1,G2,G3,CTHW))
      XLAM=DLOG(DBLE(AMN3SS**2))
      AMN3SS=SQRT(AMN3SS**2-PINENE(AMN3SS**2,G1,G2,G3,CTHW))
      AMZ1SS=TMZ1SS
      AMZ2SS=TMZ2SS
      AMZ3SS=TMZ3SS
      AMZ4SS=TMZ4SS
      AMW1SS=TMW1SS
      AMW2SS=TMW2SS
C
      XLAM=DLOG(DBLE(HIGFRZ**2))
      RETURN
      END
