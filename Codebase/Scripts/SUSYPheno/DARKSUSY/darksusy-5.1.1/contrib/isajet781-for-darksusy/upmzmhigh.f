CDECK  ID>, UPMZMHIGH.
      SUBROUTINE UPMZMHIGH
!
!Purpose: To run the RGEs up from Mz to MHIGH rotating to the current
!         basis at M_t.
!
!         Note that the SM running in *smrgedr.f contains original
!         isajet thresholds. This may not be correct for calls to 
!         *smrgedr.f once this code has derived its own thresholds,
!         but we assume here that if a threshold has not been located
!         above m_t for a particular particle, the isajet threshold is
!         a good approximation.
!
!         This first run up is to find the GUT scale values of the Yukawas
!         and mu - both of which are used to set GUT scale boundary
!         conditions. The routine requires the isajet value of mu at
!         m_SUSY as an initial guess. An alternative approach, which
!         would be independent of Isajet's running, would be to run the
!         Yukawas up and then the MSSM RGEs back down. The standard EWSB
!         conditions could then be applied to derive an initial guess
!         for the SUSY-scale mu which could then be run back up to be
!         used in the GUT scale boundary conditions.
!         If this approach was implemented, this code would be independent
!         of the Isajet RGEs and the user could be given the option (when
!         initialising Isajet) as to which RGEs they wanted to use.
!
!N.B. If RGEFLAV were indeed used as a replacement for Isajet's RGEs it
!     is essential to note the difference in conventions. We repeat here
!     that RGEFLAV is written entirely in BT (Book) notation, with the
!     exception of the MSSM RGEs and the stop decay calculation SQSIX.
!     SQSIX uses Isajet conventions and converts the output from RGEFLAV
!     appropriately. The RGEs use BT notation at all times except for the
!     MSSM RGEs which are copied from Martin and Vaughn - the output
!     is converted appropriately back to BT notation before the value for
!     the RHS of the differential equation is returned to the integrating
!     routine.
!
      IMPLICIT NONE
!
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
!
      COMMON /BSG/GISA(31),MSQISA(3),MSLISA(3),MSUISA(3),MSDISA(3),
     &            MSEISA(3),MRNISA(3),YNFRZ(3,3),MNFRZ(3,3),TNFRZ(3,3),
     &            RTISA,RBISA,RLISA
c     MSxDEC(i) - decoupling scale of i-th generation of type x sfermion
c     MRNDEC(i) - decoupling scale of i-th RH neutrino
      REAL*8 GISA,MSQISA,MSLISA,MSUISA,MSDISA,MSEISA,MRNISA,
     &       YNFRZ,MNFRZ,TNFRZ
      REAL RTISA,RBISA,RLISA
      SAVE /BSG/
!
      COMMON/WKYUK/LAMTMT,LAMBMZ,LAMTAMZ
      DOUBLE PRECISION LAMTMT,LAMBMZ,LAMTAMZ
      SAVE/WKYUK/
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
      COMMON/RGEMS/VEVMH,RGEMS,RGEMU
      DOUBLE COMPLEX VEVMH
      DOUBLE PRECISION RGEMS,RGEMU
      SAVE/RGEMS/
!
      COMMON/ATMZ/G1MZ,G2MZ,G3MZ,VSMMZ,LAMBDAMZ,LAMTMZ
      DOUBLE COMPLEX G1MZ,G2MZ,G3MZ,VSMMZ,LAMBDAMZ,LAMTMZ
      SAVE/ATMZ/
!
      COMMON/SMRGE/SMRGEMH,SMQSTEP,NU,SMDR2LP
      DOUBLE PRECISION SMRGEMH,SMQSTEP
      INTEGER NU,SMDR2LP
      SAVE/SMRGE/
!
      COMMON/LOOPS/SSQSTEP,SW2LP
      DOUBLE PRECISION SSQSTEP
      INTEGER SW2LP
      SAVE/LOOPS/
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
!
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/DEC/NEWTH,ISADEC,BELOW,NSTEPTHRESH,NLTMT,
     $           THSQ,THSU,THSD,THSL,THSE
      DOUBLE PRECISION NEWTH(20)
      INTEGER ISADEC,BELOW(20),NSTEPTHRESH(19),NLTMT
      INTEGER THSQ(3),THSU(3),THSD(3),THSL(3),THSE(3)
      SAVE/DEC/
!
      COMMON/SFMFRZ/MQSAV,MUPSAV,MDSAV,MLSAV,MESAV
      DOUBLE COMPLEX MQSAV(3,4,3),MUPSAV(3,4,3),MDSAV(3,4,3)
      DOUBLE COMPLEX MLSAV(3,4,3),MESAV(3,4,3)
      SAVE/SFMFRZ/
!
      INTEGER II,I,J,NSTEP,NSTEPSM,BELOWMS,NSTEPMT,LOOPNSTEP
      DOUBLE PRECISION TZ,TMT,TTH(20),TH,THIGH,DT,T,A1I,A2I,Q,PI
     $                ,DGSM(32),DWSM(96),DG215(215),DW215(645)
      DOUBLE PRECISION OMEGAEM,OMEGA3,A2MZMSB,A3MZMSB
      DOUBLE PRECISION SWAPVA
      DOUBLE PRECISION A1MZME,A2MZME,ALPHA3
      DOUBLE COMPLEX GSM(32),WSM(96),G215(215),W215(645)
      EXTERNAL CSMRGEDR,DSMRGEDR,CRGE215,DRGE215
!
      DO I=1,215
        G215(I)=(0.D0,0.D0)
        DG215(I)=0.D0
        IF(I.LT.33)THEN
          GSM(I)=(0.D0,0.D0)
          DGSM(I)=0.D0
        END IF
      END DO
!
      PI=4.D0*DATAN(1.D0)
!
!BELOWMS will be used when adding Yukawa finite corrections
!
      BELOWMS=1
!
!If the user wants to keep isajet's thresholds, then
!insert them into the common block
!DMAX is used to avoid negative m^2 values
!
      IF(ISADEC.EQ.1)THEN
        QTHQL(3)=DSQRT(DMAX1(1.D0,DBLE(GSS(19))))
        QTHQL(2)=DSQRT(DMAX1(1.D0,DBLE(GSS(19))))
        QTHQL(1)=DSQRT(DMAX1(1.D0,DBLE(GSS(24))))
        QTHUR(3)=DSQRT(DMAX1(1.D0,DBLE(GSS(18))))
        QTHUR(2)=DSQRT(DMAX1(1.D0,DBLE(GSS(18))))
        QTHUR(1)=DSQRT(DMAX1(1.D0,DBLE(GSS(23))))
        QTHDR(3)=DSQRT(DMAX1(1.D0,DBLE(GSS(17))))
        QTHDR(2)=DSQRT(DMAX1(1.D0,DBLE(GSS(17))))
        QTHDR(1)=DSQRT(DMAX1(1.D0,DBLE(GSS(22))))
!
        QTHLL(1)=DBLE(MSS(17))
        QTHLL(2)=DBLE(MSS(17))
        QTHLL(3)=DBLE(MSS(17))
        QTHER(1)=DBLE(MSS(17))
        QTHER(2)=DBLE(MSS(17))
        QTHER(3)=DBLE(MSS(17))
!        
        QNSH=DBLE(ABS(MU))
        QNSG=DBLE(ABS(MSS(1)))
        QNH=DBLE(MSS(30))
        QTHSB=DBLE(ABS(GSS(7)))
        QTHSW=DBLE(ABS(GSS(8)))
!
!Sort the squark and slepton eigenvectors to the same order as 
!my RGE routine
!
        DO I=1,2
          DO J=I+1,3
            IF(QTHQL(I).GT.QTHQL(J))THEN
              SWAPVA=QTHQL(J)
              QTHQL(J)=QTHQL(I)
              QTHQL(I)=SWAPVA
            END IF
            IF(QTHUR(I).GT.QTHUR(J))THEN
              SWAPVA=QTHUR(J)
              QTHUR(J)=QTHUR(I)
              QTHUR(I)=SWAPVA
            END IF
            IF(QTHDR(I).GT.QTHDR(J))THEN
              SWAPVA=QTHDR(J)
              QTHDR(J)=QTHDR(I)
              QTHDR(I)=SWAPVA
            END IF
            IF(QTHLL(I).GT.QTHLL(J))THEN
              SWAPVA=QTHLL(J)
              QTHLL(J)=QTHLL(I)
              QTHLL(I)=SWAPVA
            END IF
            IF(QTHER(I).GT.QTHER(J))THEN
              SWAPVA=QTHER(J)
              QTHER(J)=QTHER(I)
              QTHER(I)=SWAPVA
            END IF
          END DO
        END DO
!
      ELSE
!
!Right now the first guess for thresholds is the same as ISAJETs.
!DMAX is used to avoid negative m^2 values
!
        QTHQL(3)=DSQRT(DMAX1(1.D0,DBLE(GSS(19))))
        QTHQL(2)=DSQRT(DMAX1(1.D0,DBLE(GSS(19))))
        QTHQL(1)=DSQRT(DMAX1(1.D0,DBLE(GSS(24))))
        QTHUR(3)=DSQRT(DMAX1(1.D0,DBLE(GSS(18))))
        QTHUR(2)=DSQRT(DMAX1(1.D0,DBLE(GSS(18))))
        QTHUR(1)=DSQRT(DMAX1(1.D0,DBLE(GSS(23))))
        QTHDR(3)=DSQRT(DMAX1(1.D0,DBLE(GSS(17))))
        QTHDR(2)=DSQRT(DMAX1(1.D0,DBLE(GSS(17))))
        QTHDR(1)=DSQRT(DMAX1(1.D0,DBLE(GSS(22))))
!
        QTHLL(1)=DBLE(MSS(17))
        QTHLL(2)=DBLE(MSS(17))
        QTHLL(3)=DBLE(MSS(17))
        QTHER(1)=DBLE(MSS(17))
        QTHER(2)=DBLE(MSS(17))
        QTHER(3)=DBLE(MSS(17))
!        
        QNSH=DBLE(ABS(MU))
        QNSG=DBLE(ABS(MSS(1)))
        QNH=DBLE(MSS(30))
        QTHSB=DBLE(ABS(GSS(7)))
        QTHSW=DBLE(ABS(GSS(8)))
      END IF
!
!Set initial values of saved thresholds
!
      DO I=1,3
        MQSAV(I,4,I)=DCMPLX(QTHQL(I)**2,0.D0)
        MUPSAV(I,4,I)=DCMPLX(QTHUR(I)**2,0.D0)
        MDSAV(I,4,I)=DCMPLX(QTHDR(I)**2,0.D0)
        MLSAV(I,4,I)=DCMPLX(QTHLL(I)**2,0.D0)
        MESAV(I,4,I)=DCMPLX(QTHER(I)**2,0.D0)
      END DO
!
!Open a file to contain the downwards running.
!Can be reinstated (along with the appropriate
!lines below) if needed
!
!      OPEN(35,FILE='out/u1styu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!
!Sort the thresholds into size order so that the running can hit
!each one exactly
!
      CALL SORTTH
!
!Check the location of the first threshold above m_t
!
      NLTMT=0
      DO I=1,20
        IF(QTHSORT(I).LE.MT)NLTMT=I
      END DO
!
!The next line sets the upper scale so high that the couplings will
!unify before the limit is reached
!
      IF(UNI.EQ.1)MHIGH=1.D19
      NSTEP=1000
      NSTEPSM=50
!
!Use the sorted thresholds to set the number of steps between each
!threshold. Related to the number of steps in the main running.
!
      IF(NLTMT.LT.19)THEN
        IF(NLTMT.NE.0)THEN
          NSTEPTHRESH(NLTMT)=INT(ABS(DLOG(MT/QTHSORT(NLTMT+1)))
     $                                     *NSTEP/DLOG(MHIGH/MT)*10)
        END IF
        DO I=NLTMT+1,19
          NSTEPTHRESH(I)=INT(ABS(DLOG(QTHSORT(I)/QTHSORT(I+1)))
     $                                     *NSTEP/DLOG(MHIGH/MT)*10)
          IF(NSTEPTHRESH(I).EQ.0.AND.
     $                            ABS(QTHSORT(I)-QTHSORT(I+1)).GT.1D-10)
     $                                              NSTEPTHRESH(I)=10
        END DO
      END IF
!
      SMRGEMH=MHIGH
!
      TZ=LOG(MZ/MHIGH)
      TMT=LOG(MT/MHIGH)
      DO I=1,20
        TTH(I)=LOG(QTHSORT(I)/MHIGH)
      END DO
      THIGH=0.D0
!
      DT=(TMT-TZ)/FLOAT(NSTEPSM)
      NU=3 !top is decoupled at M_Z
!
!Use Arason et al. PRD46 (1992) 3945 to shift the alphas
!to account for thresholds. Include the top threshold here
!so that the top is removed from the theory at M_Z as 
!opposed to M_t. See (3.8)-(3.12) of Arason.
!
      OMEGAEM=2.D0/48.D0/PI**2*(1.D0-21.D0*DLOG(MW/MZ))
     $                         +2.D0/9.D0/PI**2*DLOG(MT/MZ)
      OMEGA3=2.D0/24.D0/PI**2*DLOG(MT/MZ)
      A1MZME=5.D0/3.D0*ALPHAEM/((1.D0-XWMSB)
     $                       *(1.D0+4.D0*PI*ALPHAEM*OMEGAEM))
      A2MZMSB=ALPHAEM/(XWMSB*(1.D0+4.D0*PI*ALPHAEM*OMEGAEM))
      A3MZMSB=ALPHASMSB/(1.D0+4.D0*PI*ALPHASMSB*OMEGA3)
!
!Now apply shifts to convert from MSBar to DRBar for alpha2
!and alpha3. There is no shift in alpha1.
!See: Martin and Vaughn, Phys.Lett.B318, 331-337 (1993)
!
      A2MZME=6.D0*PI*A2MZMSB/(6.D0*PI-A2MZMSB)
      ALPHA3=4.D0*PI*A3MZMSB/(4.D0*PI-A3MZMSB)
!
      G1MZ=DCMPLX(DBLE(SQRT(4*PI*A1MZME)))
      G2MZ=DCMPLX(DBLE(SQRT(4*PI*A2MZME)))
      G3MZ=DCMPLX(DBLE(SQRT(4*PI*ALPHA3)))
      GSM(1)=G1MZ
      GSM(2)=G2MZ
      GSM(3)=G3MZ
      GSM(4)=DCMPLX(MWEAK(1))/VSMMZ
      GSM(8)=DCMPLX(MWEAK(2))/VSMMZ
      GSM(12)=(0.D0,0.D0) !LAMTMZ is found by iteration
      GSM(13)=DCMPLX(MWEAK(3))/VSMMZ
      GSM(17)=DCMPLX(MWEAK(4))/VSMMZ
      GSM(21)=DCMPLX(LAMBMZ)
      GSM(22)=DCMPLX(MWEAK(5))/VSMMZ
      GSM(26)=DCMPLX(MWEAK(6))/VSMMZ
      GSM(30)=DCMPLX(LAMTAMZ)
      GSM(31)=(0.D0,0.D0)
      GSM(32)=VSMMZ !Set in MASS according to Pierce prescription
!
      IF(COMP.EQ.0)THEN
        DO I=1,32
          DGSM(I)=DBLE(GSM(I))
        END DO
      END IF
!
      DO II=1,NSTEPSM
        T=TZ+(TMT-TZ)*FLOAT(II-1)/FLOAT(NSTEPSM)
        SMQSTEP=SMRGEMH*EXP(T)
!
!EPS assists with locating thresholds. See the comment below.
!
        IF(II.EQ.1)EPS=ABS(SMQSTEP*(EXP(DT)-1)/(6.D0*PI))
        IF(COMP.EQ.0)THEN
          CALL DRKSTP(32,DT,T,DGSM,DSMRGEDR,DWSM)
        ELSE
          CALL CRKSTP(32,DT,T,GSM,CSMRGEDR,WSM)
        END IF
        EPS=-ABS(SMQSTEP*(EXP(DT)-1)/(6.D0*PI))
      END DO
!
      IF(COMP.EQ.0)THEN
        DO I=1,32
          GSM(I)=DCMPLX(DGSM(I))
        END DO
      END IF
!
!Now introduce the top yukawa and rotate
!
      GSM(12)=DCMPLX(LAMTMT)
      CALL ROTATESM(GSM)
!
      DO I=4,30
        IF(I.LT.7)G(I-3)=GSM(I-3)
        G(I)=(0.D0,0.D0)
        G(I+108)=GSM(I)
      END DO
      G(428)=GSM(32)
      G(429)=GSM(31)
!
      DO I=1,30
        G215(I)=G(I)
      END DO
      G215(31)=G(108)
      G215(32)=G(110)
      G215(33)=G(111)
      DO I=34,60
        G215(I)=G(I+78)
      END DO
      G215(61)=G(428)
      G215(62)=G(429)
!
!First threshold
!
      IF(NLTMT.GE.LOCMH)CALL UPMHCOND(G215)
      IF(NLTMT.NE.20)THEN
        IF(NLTMT.NE.0)THEN
          DT=(TTH(NLTMT+1)-TMT)/FLOAT(NSTEPTHRESH(NLTMT))
          LOOPNSTEP=NSTEPTHRESH(NLTMT)
        ELSE
          NSTEPMT=INT(ABS(DLOG(MT/QTHSORT(1)))*NSTEP/DLOG(MHIGH/MT)*10)
          DT=(TTH(NLTMT+1)-TMT)/FLOAT(NSTEPMT)
          LOOPNSTEP=NSTEPMT
        END IF
!
        DO II=1,LOOPNSTEP
          T=TMT+(TTH(NLTMT+1)-TMT)*FLOAT(II-1)/FLOAT(LOOPNSTEP)
          SSQSTEP=MHIGH*EXP(T) !Q for the current step
!
!Introduce the MSSM mu when we pass the SUSY scale. If m_H
!has already been passed, then we need to carry out the usual
!boundary condition replacement here
!
!Also, take care of the finite corrections to the Yukawas
!
        IF(BELOWMS.EQ.1.AND.MHIGH*EXP(T+DT).GE.RGEMS)THEN
          BELOWMS=0
          CALL ROTBACK215(G215)
          G215(215)=DCMPLX(RGEMU)
          G215(31)=G215(215)
          G215(12)=G215(12)/(1.D0-DBLE(RTISA))
          G215(21)=G215(21)/(1.D0-DBLE(RBISA))
          G215(30)=G215(30)/(1.D0-DBLE(RLISA))
          IF(RGEMS.LE.QNH)THEN
            G215(42)=G215(42)/(1.D0-DBLE(RTISA))
            G215(51)=G215(51)/(1.D0-DBLE(RBISA))
            G215(60)=G215(60)/(1.D0-DBLE(RLISA))
          END IF
          CALL ROTATE215(G215)
        END IF
!
!EPS is used in the RGE subroutine when deciding which particles are
!decoupled. If it is positive, we are calculating using the particle
!content just above the point. Negative, the particle content is that
!below the point. The sign of EPS is only tested if the difference between
!SSQSTEP and Q for the threshold is less than 1/(6*pi)*DQ. Thus the sign
!is only important when near a threshold, hence the IF(II.EQ.1). Once we
!have taken the first step, EPS can change sign, but must continue
!to be evaluated since DQ is a function of Q. 
!
          IF(II.EQ.1)EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
          IF(COMP.EQ.0)THEN
            DO I=1,215
              DG215(I)=DBLE(G215(I))
            END DO
            CALL DRKSTP(215,DT,T,DG215,DRGE215,DW215)
            DO I=1,215
              G215(I)=DCMPLX(DG215(I))
            END DO
          ELSE
            CALL CRKSTP(215,DT,T,G215,CRGE215,W215)
          END IF
!          CALL UOUTCOUP(G215,MHIGH*EXP(T))
          EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
        END DO
!
!Matching conditions at M_H
!
        IF(NLTMT+1.EQ.LOCMH)CALL UPMHCOND(G215)
!
        DO I=NLTMT+2,20
!
!Don't run between degenerate thresholds
!
          IF(NSTEPTHRESH(I-1).EQ.0)GOTO 20
!
          DT=(TTH(I)-TTH(I-1))/FLOAT(NSTEPTHRESH(I-1))
!
          DO II=1,NSTEPTHRESH(I-1)
            T=TTH(I-1)+(TTH(I)-TTH(I-1))
     $                         *FLOAT(II-1)/FLOAT(NSTEPTHRESH(I-1))
            SSQSTEP=MHIGH*EXP(T)
!
            IF(BELOWMS.EQ.1.AND.MHIGH*EXP(T+DT).GE.RGEMS)THEN
              BELOWMS=0
              CALL ROTBACK215(G215)
              G215(215)=DCMPLX(RGEMU)
              G215(31)=G215(215)
              G215(12)=G215(12)/(1.D0-DBLE(RTISA))
              G215(21)=G215(21)/(1.D0-DBLE(RBISA))
              G215(30)=G215(30)/(1.D0-DBLE(RLISA))
              IF(RGEMS.LE.QNH)THEN
                G215(42)=G215(42)/(1.D0-DBLE(RTISA))
                G215(51)=G215(51)/(1.D0-DBLE(RBISA))
                G215(60)=G215(60)/(1.D0-DBLE(RLISA))
              END IF
              CALL ROTATE215(G215)
            END IF
!
            IF(II.EQ.1)EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
            IF(COMP.EQ.0)THEN
              DO J=1,215
                DG215(J)=DBLE(G215(J))
              END DO
              CALL DRKSTP(215,DT,T,DG215,DRGE215,DW215)
              DO J=1,215
                G215(J)=DCMPLX(DG215(J))
              END DO
            ELSE
              CALL CRKSTP(215,DT,T,G215,CRGE215,W215)
            END IF
!            CALL UOUTCOUP(G215,MHIGH*EXP(T))
            EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
          END DO
!
   20     IF(I.EQ.LOCMH)CALL UPMHCOND(G215)
        END DO
      END IF
!
!Continue running to m_high. I need to be careful about the case that
!all thresholds are below m_t.
!
      IF(NLTMT.NE.20)THEN
        DT=(THIGH-TTH(20))/FLOAT(NSTEP)
      ELSE
        DT=(THIGH-TMT)/FLOAT(NSTEP)
      END IF
!
      DO II=1,NSTEP
        IF(NLTMT.NE.20)THEN
          T=TTH(20)+(THIGH-TTH(20))*FLOAT(II-1)/FLOAT(NSTEP)
        ELSE
          T=TMT+(THIGH-TMT)*FLOAT(II-1)/FLOAT(NSTEP)
        END IF
        SSQSTEP=MHIGH*EXP(T)
!
        IF(BELOWMS.EQ.1.AND.MHIGH*EXP(T+DT).GE.RGEMS)THEN
          BELOWMS=0
          CALL ROTBACK215(G215)
          G215(215)=DCMPLX(RGEMU)
          G215(31)=G215(215)
          G215(12)=G215(12)/(1.D0-DBLE(RTISA))
          G215(21)=G215(21)/(1.D0-DBLE(RBISA))
          G215(30)=G215(30)/(1.D0-DBLE(RLISA))
          IF(RGEMS.LE.QNH)THEN
            G215(42)=G215(42)/(1.D0-DBLE(RTISA))
            G215(51)=G215(51)/(1.D0-DBLE(RBISA))
            G215(60)=G215(60)/(1.D0-DBLE(RLISA))
          END IF
          CALL ROTATE215(G215)
        END IF
!
        Q=MHIGH*EXP(T)
        IF(COMP.EQ.0)THEN
          A1I=4.D0*PI/DG215(1)**2
          A2I=4.D0*PI/DG215(2)**2
        ELSE
          A1I=4.D0*PI/DBLE(G215(1)**2)
          A2I=4.D0*PI/DBLE(G215(2)**2)
        END IF
        IF(A1I.LT.A2I.AND.UNI.EQ.1)THEN
          MHIGH=Q
          GO TO 30
        END IF
        IF(II.EQ.1)EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
        IF(COMP.EQ.0)THEN
          DO I=1,215
            DG215(I)=DBLE(G215(I))
          END DO
          CALL DRKSTP(215,DT,T,DG215,DRGE215,DW215)
          DO I=1,215
            G215(I)=DCMPLX(DG215(I))
          END DO
        ELSE
          CALL CRKSTP(215,DT,T,G215,CRGE215,W215)
        END IF
!        CALL UOUTCOUP(G215,MHIGH*EXP(T))
        EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
      END DO
!
!NB: If this message is written, the programme will probably
!    not succeed. It may be better to put a stop statement here.
!
      IF(UNI.EQ.1)WRITE(*,*)'ERROR: UNIFICATION NOT FOUND'
 30   CONTINUE
!
      DO I=1,30
        G(I)=G215(I)
      END DO
      G(108)=G215(31)
      G(110)=G215(32)
      G(111)=G215(33)
      DO I=34,60
        G(I+78)=G215(I)
      END DO
      G(428)=G215(61)
      G(429)=G215(62)
!
      CLOSE(35)
      RETURN
      END
