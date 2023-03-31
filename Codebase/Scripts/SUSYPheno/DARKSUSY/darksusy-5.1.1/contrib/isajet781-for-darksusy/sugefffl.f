CDECK  ID>, SUGEFFFL.
C-----------------------------------------------------------------
      SUBROUTINE SUGEFFFL(Q,G0,MUQ,SIG1,SIG2)
C-----------------------------------------------------------------
C
C     Compute Higgs mass shift due to 1-loop effective potential
C     remaining tadpoles added by Javier Ferrandis 5/20/03
!
!SUGEFFFL is the version of SUGEFF which has been altered to work
!with RGEFLAV. The input couplings G0 are larger in number and in
!the same format as in RGE601. Also accepts value of MUQ which is
!the Higgsino mass term. The sigmas will be used to refine the values
!of m_1 and m_2.
!
!The equations have been adapted to use the couplings present
!in the full tilde RGEs. Specifically, the combinations
!\tilde{\mu}^* f_{u,d,e} remain intact and the difference between
!\mu and \tilde{\mu} is respected.
!
!NOTE: G0 must be passed in the quark mass basis and with
!      the soft masses set equal to their decoupled values if
!      Q (=m_SUSY) is less than the decoupling scale.
!
      IMPLICIT NONE
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
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
      REAL Q,MUQ,SIG1,SIG2
      REAL DT1,DELT1S,SIG1T1,DMSDV2,FB,FT,MST2,MSB2,MSB1,SIG2B1,
     $SIG1B1,SIG2B2,SIG1B2,DB1,SIG1T2,SIG2T1,DELB1S,SIG2T2,MST1,COS2W,
     $MT,PI,COTB,TANB,MB,SIG2B,SIG1B,G,GGP,SIG2T,E,FAC,SIG1T,QS,
     $BETA,SINB,COSB,FAC4,FL,ML,SIG1L,SIG2L,MSL1,MSL2,
     $DELL1S,DL1,SIG1L1,SIG2L1,SIG1L2,SIG2L2,VUQ,VDQ
      REAL VEVQ,MWQ,MZQ
      REAL M1Q,M2Q,MPHO
      REAL COS2B
      REAL SUGFN,QSC
      REAL ZESS1,ZESS2,ZESS3,ZESS4
      REAL W1SS,W2SS
      REAL MHL,MHH,MA
      REAL SIG1C1,SIG1C2,SIG2C1,SIG2C2,SIG1C,SIG2C
      REAL SIG1HL,SIG1HH,SIG2HL,SIG2HH,SIG1HI,SIG2HI
      REAL SIGNE1,SIGNE2,SIGNEUD,SIGNEUU,SIGNEDD
      REAL G0(601)
C     Statement function
      SUGFN(QSC)=QSC**2*(LOG(QSC**2/Q**2)-1.)
C
      G=G2
      COS2W=1.-XW
      VUQ=G0(110)
      VDQ=G0(111)
      TANB=VUQ/VDQ
      COTB=1./TANB
      BETA=ATAN(TANB)
      SINB=SIN(BETA)
      COSB=COS(BETA)
      PI=4.*ATAN(1.)
      FAC=3./8./PI**2
      FAC4=FAC/3.
      E=EXP(1.)
      QS=Q**2
C-----INITIALIZE EFF. POT'L CONTRIBUTIONS
      SIG1B=0.
      SIG1B1=0.
      SIG1B2=0.
      SIG1T=0.
      SIG1T1=0.
      SIG1T2=0.
      SIG1L=0.
      SIG1L1=0.
      SIG1L2=0.
      SIGNE1=0.
      SIG1C=0.
      SIG1HI=0.
      SIG2B=0.
      SIG2B1=0.
      SIG2B2=0.
      SIG2T=0.
      SIG2T1=0.
      SIG2T2=0.
      SIG2L=0.
      SIG2L1=0.
      SIG2L2=0.
      SIGNE2=0.
      SIG2C=0.
      SIG2HI=0.
C-----CALCULATE TOP AND BOTTOM CONTRIBUTIONS; USE RUNNING MASSES--
      FL=G0(30)
      FB=G0(21)
      FT=G0(12)
      ML=FL*VDQ
      MB=FB*VDQ
      MT=FT*VUQ
      SIG1T=0.
      SIG2T=-FAC*MT**2*G0(12)**2*LOG(MT**2/E/QS)
      SIG1B=-FAC*MB**2*G0(21)**2*LOG(MB**2/E/QS)
      SIG2B=0.
      SIG1L=-FAC4*ML**2*G0(30)**2*LOG(ML**2/E/QS)
      SIG2L=0.
      GGP=(G**2+GP**2)/2.
      MST1=MSS(12)
      MST2=MSS(13)
      MSB1=MSS(10)
      MSB2=MSS(11)
      MSL1=MSS(21)
      MSL2=MSS(22)
      VEVQ=SQRT(VUQ**2+VDQ**2)
      MWQ=G/SQRT(2.)*VEVQ
      MZQ=SQRT(GGP)*VEVQ
      COS2B=SIN(BETA)
C-----CALCULATE STOP_1 CONTRIBUTION -------------------------------
      DELT1S=(.5*(G0(71)-G0(89))+(8*COS2W-5.)*GGP*
     $      (VDQ**2-VUQ**2)/12.)**2+VUQ**2
     $      *(G0(42)-G0(438)*COTB)**2
      DT1=.5*(G0(71)-G0(89))+(8*COS2W-5.)*GGP*
     $       (VDQ**2-VUQ**2)/12.
      DMSDV2=GGP/4.-(2*DT1*(8*COS2W-5.)*GGP/12.-
     $        G0(438)*(G0(42)*TANB-G0(438)))/2./SQRT(DELT1S)
      SIG1T1=FAC/2.*MST1**2*LOG(MST1**2/E/QS)*DMSDV2
      DMSDV2=-GGP/4.+FT**2-(-2*DT1*(8*COS2W-5.)*GGP/12.+
     $        G0(42)*(G0(42)-G0(438)*COTB))/2./SQRT(DELT1S)
      SIG2T1=FAC/2.*MST1**2*LOG(MST1**2/E/QS)*DMSDV2
C-----CALCULATE STOP_2 CONTRIBUTION -------------------------------
      DMSDV2=GGP/4.+(2*DT1*(8*COS2W-5.)*GGP/12.-
     $        G0(438)*(G0(42)*TANB-G0(438)))/2./SQRT(DELT1S)
      SIG1T2=FAC/2.*MST2**2*LOG(MST2**2/E/QS)*DMSDV2
      DMSDV2=-GGP/4.+FT**2+(-2*DT1*(8*COS2W-5.)*GGP/12.+
     $        G0(42)*(G0(42)-G0(438)*COTB))/2./SQRT(DELT1S)
      SIG2T2=FAC/2.*MST2**2*LOG(MST2**2/E/QS)*DMSDV2
C-----CALCULATE SBOT_1 CONTRIBUTION -------------------------------
      DELB1S=(.5*(G0(71)-G0(98))-(4*COS2W-1.)*GGP*
     $       (VDQ**2-VUQ**2)/12.)**2+VDQ**2
     $       *(G0(51)-G0(447)*TANB)**2
      DB1=.5*(G0(71)-G0(98))-(4*COS2W-1.)*GGP*
     $       (VDQ**2-VUQ**2)/12.
      DMSDV2=-GGP/4.+FB**2-(-2*DB1*(4*COS2W-1.)*GGP/12.+
     $        G0(51)*(G0(51)-G0(447)*TANB))/2./SQRT(DELB1S)
      SIG1B1=FAC/2.*MSB1**2*LOG(MSB1**2/E/QS)*DMSDV2
      DMSDV2=GGP/4.-(2*DB1*(4*COS2W-1.)*GGP/12.-
     $        G0(447)*(G0(51)*COTB-G0(447)))/2./SQRT(DELB1S)
      SIG2B1=FAC/2.*MSB1**2*LOG(MSB1**2/E/QS)*DMSDV2
C-----CALCULATE SBOT_2 CONTRIBUTION -------------------------------
      DMSDV2=-GGP/4.+FB**2+(-2*DB1*(4*COS2W-1.)*GGP/12.+
     $        G0(51)*(G0(51)-G0(447)*TANB))/2./SQRT(DELB1S)
      SIG1B2=FAC/2.*MSB2**2*LOG(MSB2**2/E/QS)*DMSDV2
      DMSDV2=GGP/4.+(2*DB1*(4*COS2W-1.)*GGP/12.-
     $        G0(447)*(G0(51)*COTB-G0(447)))/2./SQRT(DELB1S)
      SIG2B2=FAC/2.*MSB2**2*LOG(MSB2**2/E/QS)*DMSDV2
C-----CALCULATE STAU_1 CONTRIBUTION -------------------------------
      DELL1S=(.5*(G0(80)-G0(107))-(4*COS2W-3.)*GGP*
     $       (VDQ**2-VUQ**2)/4.)**2+VDQ**2
     $       *(G0(60)-G0(456)*TANB)**2
      DL1=.5*(G0(80)-G0(107))-(4*COS2W-3.)*GGP*
     $       (VDQ**2-VUQ**2)/4.
      DMSDV2=-GGP/4.+FL**2-(-2*DL1*(4*COS2W-3.)*GGP/4.+
     $        G0(60)*(G0(60)-G0(456)*TANB))/2./SQRT(DELL1S)
      SIG1L1=FAC4/2.*MSL1**2*LOG(MSL1**2/E/QS)*DMSDV2
      DMSDV2=GGP/4.-(2*DL1*(4*COS2W-3.)*GGP/4.-
     $        G0(456)*(G0(60)*COTB-G0(456)))/2./SQRT(DELL1S)
      SIG2L1=FAC4/2.*MSL1**2*LOG(MSL1**2/E/QS)*DMSDV2
C-----CALCULATE STAU_2 CONTRIBUTION -------------------------------
      DMSDV2=-GGP/4.+FL**2+(-2*DL1*(4*COS2W-3.)*GGP/4.+
     $        G0(60)*(G0(60)-G0(456)*TANB))/2./SQRT(DELL1S)
      SIG1L2=FAC4/2.*MSL2**2*LOG(MSL2**2/E/QS)*DMSDV2
      DMSDV2=GGP/4.+(2*DL1*(4*COS2W-3.)*GGP/4.-
     $        G0(456)*(G0(60)*COTB-G0(456)))/2./SQRT(DELL1S)
      SIG2L2=FAC4/2.*MSL2**2*LOG(MSL2**2/E/QS)*DMSDV2
C-----CALCULATE NEUTRALINO CONTRIBUTION -------------------------------
      ZESS1=MSS(23)
      ZESS2=MSS(24)
      ZESS3=MSS(25)
      ZESS4=MSS(26)
      if(((ZESS1-ZESS2).eq.0.).
     #     or.((ZESS1-ZESS3).eq.0.).
     #     or.((ZESS1-ZESS4).eq.0.).
     #     or.((ZESS2-ZESS3).eq.0.).
     #     or.((ZESS2-ZESS4).eq.0.).
     #     or.((ZESS3-ZESS4).eq.0.)) then
        SIGNE1=0.
        SIGNE2=0.
        goto 111
      else
        continue
      endif
      M1Q=G0(31)
      M2Q=G0(32)
      MPHO=COS2W*M1Q+XW*M2Q
      SIGNEUU=-GGP/4./2./PI**2*(
     # ZESS1**2*(ZESS1-MPHO)/(ZESS1-ZESS2)/(ZESS1-ZESS3)
     #     /(ZESS1-ZESS4)*SUGFN(ZESS1)+
     # ZESS2**2*(ZESS2-MPHO)/(ZESS2-ZESS1)/(ZESS2-ZESS3)
     #     /(ZESS2-ZESS4)*SUGFN(ZESS2)+
     # ZESS3**2*(ZESS3-MPHO)/(ZESS3-ZESS1)/(ZESS3-ZESS2)
     #    /(ZESS3-ZESS4)*SUGFN(ZESS3)+
     # ZESS4**2*(ZESS4-MPHO)/(ZESS4-ZESS1)/(ZESS4-ZESS2)
     #   /(ZESS4-ZESS3)*SUGFN(ZESS4))
      SIGNEDD=SIGNEUU
      SIGNEUD=-GGP/4./2./PI**2*MUQ*(
     # ZESS1*(ZESS1-MPHO)/(ZESS1-ZESS2)/(ZESS1-ZESS3)
     #     /(ZESS1-ZESS4)*SUGFN(ZESS1)+
     # ZESS2*(ZESS2-MPHO)/(ZESS2-ZESS1)/(ZESS2-ZESS3)
     #     /(ZESS2-ZESS4)*SUGFN(ZESS2)+
     # ZESS3*(ZESS3-MPHO)/(ZESS3-ZESS1)/(ZESS3-ZESS2)
     #    /(ZESS3-ZESS4)*SUGFN(ZESS3)+
     # ZESS4*(ZESS4-MPHO)/(ZESS4-ZESS1)/(ZESS4-ZESS2)
     #   /(ZESS4-ZESS3)*SUGFN(ZESS4))
      SIGNE1=SIGNEDD+TANB*SIGNEUD
      SIGNE2=SIGNEUU+COTB*SIGNEUD
C-----CALCULATE CHARGINO CONTRIBUTION -----------
      W1SS=MSS(23)
      W2SS=MSS(24)
      if((W1SS-W2SS).eq.0.) then
          SIG1C=0.
          SIG2C=0.
         goto 111
      else
         continue
      endif
      SIG1C1=-G**2/16./PI**2*(1.-(2.*MWQ**2*COS2B+M2Q**2+MUQ**2)/
     #   (W2SS**2-W1SS**2)-
     #   TANB*2.*M2Q*MUQ/(W2SS**2-W1SS**2))*SUGFN(W1SS)
      SIG1C2=-G**2/16./PI**2*(1.+(2.*MWQ**2*COS2B+M2Q**2+MUQ**2)/
     #   (W2SS**2-W1SS**2)+
     #   TANB*2.*M2Q*MUQ/(W2SS**2-W1SS**2))*SUGFN(W2SS)
      SIG2C1=-G**2/16./PI**2*(1.-(-2.*MWQ**2*COS2B+M2Q**2+MUQ**2)/
     #   (W2SS**2-W1SS**2)-
     #   COTB*2.*M2Q*MUQ/(W2SS**2-W1SS**2))*SUGFN(W1SS)
      SIG2C2=-G**2/16./PI**2*(1.+(-2.*MWQ**2*COS2B+M2Q**2+MUQ**2)/
     #   (W2SS**2-W1SS**2)+
     #   COTB*2.*M2Q*MUQ/(W2SS**2-W1SS**2))*SUGFN(W2SS)
      SIG1C=SIG1C1+SIG1C2
      SIG2C=SIG2C1+SIG2C2
C-----CALCULATE HIGGSES CONTRIBUTION -------------------------------
      MHL=MSS(29)
      MHH=MSS(30)
      MA=MSS(31)
      IF(MHL.LT.1.) THEN
	      SIG1HI=0.
              SIG2HI=0.
	      GOTO 111
      ELSE
	      CONTINUE
      ENDIF
      IF(MHH.LT.1.) THEN
	      SIG1HI=0.
              SIG2HI=0.
	      GOTO 111
      ELSE
	      CONTINUE
      ENDIF
      SIG1HL=GGP**2/32./PI**2*(1.-(MZQ**2+
     #    MA**2*(1.-4.*COS2B+2.*COS2B**2))/
     #   (MHH**2-MHL**2))*SUGFN(MHL)
      SIG1HH=GGP**2/32./PI**2*(1.-(MZQ**2+
     #    MA**2*(1.-4.*COS2B+2.*COS2B**2))/
     #   (MHH**2-MHL**2))*SUGFN(MHH)
      SIG2HL=GGP**2/32./PI**2*(1.-(MZQ**2+
     #   MA**2*(1.+4.*COS2B+2.*COS2B**2))/
     #   (MHH**2-MHL**2))*SUGFN(MHL)
      SIG2HH=GGP**2/32./PI**2*(1.-(MZQ**2+
     #   MA**2*(1.+4.*COS2B+2.*COS2B**2))/
     #   (MHH**2-MHL**2))*SUGFN(MHH)
      SIG1HI=SIG1HL+SIG1HH
      SIG2HI=SIG2HL+SIG2HH
111    CONTINUE
C-----ADD ALL TERMS ------------------------------------------------
      SIG1=SIG1B+SIG1B1+SIG1B2+SIG1T+SIG1T1+SIG1T2+
     $SIG1L+SIG1L1+SIG1L2+SIGNE1+SIG1C+SIG1HI
      SIG2=SIG2B+SIG2B1+SIG2B2+SIG2T+SIG2T1+SIG2T2+
     $SIG2L+SIG2L1+SIG2L2+SIGNE2+SIG2C+SIG2HI
      RETURN
      END
