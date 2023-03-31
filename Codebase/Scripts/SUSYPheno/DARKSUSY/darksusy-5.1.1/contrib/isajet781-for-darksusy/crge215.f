!
      SUBROUTINE CRGE215(T,G,F)
!
!Contains threshold RGEs for gauge and yukawas.
!
!CRGE215 IS FOR RUNNNG THRESHOLDS, NO TILDES, FOR THE FIRST TIME
!The sfermion thresholds are not all distinct. All left-squarks
!are at the same point, all up-right-squarks are at another point,
!all down-right-squarks at another, all left-sleptons another and
!all right-sleptons another.
!
!     G(  1) = g_1         G(  2) = g_2         G(  3) = g_3
!     G(  4) = FU(1,1)     G(  5) = FU(1,2)     G( 12) = FU(3,3)
!     G( 13) = FD(1,1)     G( 22) = FE(1,1)     G( 30) = FE(3,3)
!     G( 31) = mu          G( 32) = V_U         G( 33) = V_D
!
!     G( 34) = LU(1,1)     G( 43) = LD(1,1)     G( 52) = LE(1,1)
!
!     G( 61) = VEV_SM      G( 62) = LAM_SM
!
!     G( 63) = GTPQ(1,1)   G( 72) = GTPL(1,1)   G( 81) = GTPU(1,1)
!     G( 90) = GTPD(1,1)   G( 99) = GTPE(1,1)   G(108) = GTP_Hu
!     G(109) = GTP_Hd      G(110) = GTQ(1,1)    G(119) = GTL(1,1)
!     G(128) = GT_Hu       G(129) = GT_Hd       G(130) = GTSQ(1,1)
!     G(139) = GTSU(1,1)   G(148) = GTSD(1,1)   G(157) = FTUQ(1,1)
!     G(166) = FTDQ(1,1)   G(175) = FTEL(1,1)   G(184) = FTUU(1,1)
!     G(193) = FTDD(1,1)   G(202) = FTEE(1,1)   G(211) = sGTP_Hu
!     G(212) = cGTP_Hd     G(213) = sGT_Hu      G(214) = cGT_Hd
!
!     G(215) = mu(M)
!
!This is the BT version which receives G in book notation
!
      IMPLICIT NONE
!
      COMMON/LOOPS/SSQSTEP,SW2LP
      DOUBLE PRECISION SSQSTEP
      INTEGER SW2LP
      SAVE/LOOPS/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      DOUBLE PRECISION T
      DOUBLE COMPLEX G(215)
      DOUBLE COMPLEX F(215)
!
      DOUBLE COMPLEX FU(3,3),FD(3,3),FE(3,3)
      DOUBLE COMPLEX YU(3,3),YD(3,3),YE(3,3)
      DOUBLE COMPLEX AU(3,3),AD(3,3),AE(3,3)
      DOUBLE COMPLEX HU(3,3),HD(3,3),HE(3,3),DH(3,3,3)
      DOUBLE COMPLEX MQ(3,3),ML(3,3),MUP(3,3),MD(3,3),ME(3,3),DM(5,3,3)
      DOUBLE COMPLEX LU(3,3),LD(3,3),LE(3,3),LYU(3,3),LYD(3,3),LYE(3,3)
      DOUBLE COMPLEX GTPQ(3,3),GTPL(3,3),GTPU(3,3),GTPD(3,3),GTPE(3,3)
      DOUBLE COMPLEX GTQ(3,3),GTL(3,3),GTSQ(3,3),GTSU(3,3),GTSD(3,3)
      DOUBLE COMPLEX FTUQ(3,3),FTDQ(3,3),FTEL(3,3)
      DOUBLE COMPLEX FTUU(3,3),FTDD(3,3),FTEE(3,3)
!
      DOUBLE COMPLEX CMATMUL,CSFMUL,CTRACE,SUM
      DOUBLE COMPLEX CMODSQ,CCON,CRE
!
!These are used in the calculation of the RGEs which contain
!thresholds
!
      DOUBLE COMPLEX DUMU1(3,3),DUMU2(3,3),DUMD1(3,3),DUMD2(3,3)
      DOUBLE COMPLEX DUME1(3,3),DUME2(3,3),DUMGRKMU1(3,3),TDUMGRKMU
      DOUBLE COMPLEX DUMLUD1(3,3),DUMLUD2(3,3),DUMLUD(3,3),TDUMLUD
!
      DOUBLE COMPLEX FUS(3,3),FDS(3,3),FES(3,3)
      DOUBLE COMPLEX FTUQD(3,3),FTUUD(3,3),FTDDD(3,3),FTDQD(3,3)
      DOUBLE COMPLEX FTELD(3,3),FTEED(3,3)
      DOUBLE COMPLEX GTPUS(3,3),GTSUS(3,3),GTQS(3,3),GTLS(3,3)
      DOUBLE COMPLEX GTPQS(3,3),GTPLS(3,3),GTSQS(3,3),GTPDS(3,3)
      DOUBLE COMPLEX GTPES(3,3),GTSDS(3,3),GTPUT(3,3),GTQT(3,3)
      DOUBLE COMPLEX GTLT(3,3),GTPQT(3,3),GTPLT(3,3),GTPDT(3,3)
      DOUBLE COMPLEX GTPET(3,3),GTSUT(3,3),GTSQT(3,3),GTSDT(3,3)
!
      DOUBLE COMPLEX FUFUD(3,3),FDFDD(3,3),FEFED(3,3)
      DOUBLE COMPLEX YUDYU(3,3),YDDYD(3,3),YEDYE(3,3)
      DOUBLE COMPLEX LULUD(3,3),LDLDD(3,3),LELED(3,3)
      DOUBLE COMPLEX LYUDLYU(3,3),LYDDLYD(3,3),LYEDLYE(3,3)
      DOUBLE COMPLEX LYUDLYU2(3,3),LYDDLYD2(3,3),LYEDLYE2(3,3)
      DOUBLE COMPLEX LYUDLYU3(3,3),LYDDLYD3(3,3),LYEDLYE3(3,3)
      DOUBLE COMPLEX TLYUDLYU,TLYDDLYD,TLYEDLYE
      DOUBLE COMPLEX TLYUDLYU2,TLYDDLYD2,TLYEDLYE2
      DOUBLE COMPLEX TLYUDLYU3,TLYDDLYD3,TLYEDLYE3
      DOUBLE COMPLEX LYUDLYULYDDLYD(3,3),LYDDLYDLYUDLYU(3,3)
      DOUBLE COMPLEX TLYUDLYULYDDLYD
      DOUBLE COMPLEX FUDFU(3,3),FDDFD(3,3),FEDFE(3,3)
      DOUBLE COMPLEX YUYUD(3,3),YDYDD(3,3),YEYED(3,3)
      DOUBLE COMPLEX YUYDD(3,3),YDYUD(3,3)
      DOUBLE COMPLEX TFUDFU,TFDDFD,TFEDFE
      DOUBLE COMPLEX TYUDYU,TYDDYD,TYEDYE
      DOUBLE COMPLEX TLUDLU,TLDDLD,TLEDLE
      DOUBLE COMPLEX FUFUDFU(3,3),FUFUDFD(3,3),FDFDDFU(3,3)
      DOUBLE COMPLEX FDFDDFD(3,3),FEFEDFE(3,3)
      DOUBLE COMPLEX LULUDLU(3,3),LULUDLD(3,3),LDLDDLU(3,3)
      DOUBLE COMPLEX LDLDDLD(3,3),LELEDLE(3,3)
      DOUBLE COMPLEX SQFTUQDFTUQ(3,3),SQFTDQDFTDQ(3,3),SQGTQTGTQS(3,3)
      DOUBLE COMPLEX SQGTQTFTUQ(3,3),SQGTQTFTDQ(3,3),SQGTPQTGTPQS(3,3)
      DOUBLE COMPLEX SQGTPQTFTUQ(3,3),SQGTPQTFTDQ(3,3),SQGTSQTGTSQS(3,3)
      DOUBLE COMPLEX SUGTPUSGTPUT(3,3),SUGTSUSGTSUT(3,3)
      DOUBLE COMPLEX SUFTUUFTUUD(3,3),SUFTUUGTPUT(3,3),SDFTDDFTDDD(3,3)
      DOUBLE COMPLEX SDFTDDGTPDT(3,3),SDGTPDSGTPDT(3,3)
      DOUBLE COMPLEX SDGTSDSGTSDT(3,3),SLFTELDFTEL(3,3),SLGTLTGTLS(3,3)
      DOUBLE COMPLEX SLGTLTFTEL(3,3),SLGTPLTGTPLS(3,3),SLGTPLTFTEL(3,3)
      DOUBLE COMPLEX SEGTPESGTPET(3,3),SEFTEEFTEED(3,3),SEFTEEGTPET(3,3)
      DOUBLE COMPLEX MGTPHUSQ,MGTPHDSQ,MGTHUSQ,MGTHDSQ
      DOUBLE COMPLEX MSGTPHUSQ,MCGTPHDSQ,MSGTHUSQ,MCGTHDSQ
      DOUBLE COMPLEX Y2,H,Y4,CHI4
!
      DOUBLE COMPLEX B1U(3,3),B1D(3,3),B1E(3,3)
      DOUBLE COMPLEX BGRKMU,BETALAM1,BETALAM2

!
!The following are for two loop and soft RGEs
!
      DOUBLE COMPLEX DUM2U1(3,3),DUM2U2(3,3)
      DOUBLE COMPLEX DUM2D1(3,3),DUM2D2(3,3)
      DOUBLE COMPLEX DUM2E1(3,3),DUM2E2(3,3)
      DOUBLE COMPLEX DUM2U(3,3),DUM2D(3,3),DUM2E(3,3)
      DOUBLE COMPLEX DUM1GRKMU(3,3),DUM2GRKMU(3,3)
!
      DOUBLE COMPLEX YUYUDYUYUD(3,3),YDYDDYDYDD(3,3)
      DOUBLE COMPLEX YEYEDYEYED(3,3),YUDYUYUDYU(3,3)
      DOUBLE COMPLEX YDDYDYDDYD(3,3),YEDYEYEDYE(3,3)
      DOUBLE COMPLEX YDDYDYUDYU(3,3),YUDYUYDDYD(3,3)
      DOUBLE COMPLEX YUYDDYDYUD(3,3)
!
      DOUBLE COMPLEX TYUDYUYUDYU,TYDDYDYDDYD,TYEDYEYEDYE
      DOUBLE COMPLEX TYUDYUYDDYD,TYDDYDYUDYU
!
      DOUBLE COMPLEX BETA1U(3,3),BETA2U(3,3),BETA1D(3,3),BETA2D(3,3)
      DOUBLE COMPLEX BETA1E(3,3),BETA2E(3,3)
      DOUBLE COMPLEX BETA2GRKMU,BETA1VU,BETA1VD,BETA2VU,BETA2VD,BETA2B
      DOUBLE COMPLEX B2YMU(3,3),B2YMD(3,3),B2YME(3,3),BETAVEV1,BETAVEV2
!
      DOUBLE COMPLEX ID(3,3),MVMU,MVMUM
      DOUBLE PRECISION B1LP(3),B1LPM(3),B2LPSM(3,3),B2LPM(3,3)
      DOUBLE PRECISION CM(3,3),CSM(3,3)
      DOUBLE PRECISION PI,Q
!
      INTEGER I,J,NG,ND,NE,NNU,NU,NSQ,NSU,NSD,NSL,NSE,NSH,NH,NSW,NSG
      INTEGER THLH,THHH,THSH,THSQ(3),THSU(3),THSD(3),THSL(3)
      INTEGER THSE(3),THSB,THSW,THGL
!
      DATA ID(1,1)/(1.D0,0.D0)/,ID(1,2)/(0.D0,0.D0)/
     $    ,ID(1,3)/(0.D0,0.D0)/
      DATA ID(2,1)/(0.D0,0.D0)/,ID(2,2)/(1.D0,0.D0)/
     $    ,ID(2,3)/(0.D0,0.D0)/
      DATA ID(3,1)/(0.D0,0.D0)/,ID(3,2)/(0.D0,0.D0)/
     $    ,ID(3,3)/(1.D0,0.D0)/
      DATA B1LPM(1)/6.6D0/,B1LPM(2)/1.D0/,B1LPM(3)/-3.D0/
      DATA B2LPM(1,1)/7.96D0/,B2LPM(1,2)/5.4D0/,B2LPM(1,3)/17.6D0/
      DATA B2LPM(2,1)/1.8D0/,B2LPM(2,2)/25.D0/,B2LPM(2,3)/24.D0/
      DATA B2LPM(3,1)/2.2D0/,B2LPM(3,2)/9.D0/,B2LPM(3,3)/14.D0/
      DATA CM(1,1)/5.2D0/,CM(1,2)/2.8D0/,CM(1,3)/3.6D0/
      DATA CM(2,1)/6.D0/,CM(2,2)/6.D0/,CM(2,3)/2.D0/
      DATA CM(3,1)/4.D0/,CM(3,2)/4.D0/,CM(3,3)/0.D0/
!
!Set all F's and betas to zero
!
      DO I=1,215
        F(I)=(0.D0,0.D0)
      END DO
      DO I=1,3
        DO J=1,3
          B1U(I,J)=(0.D0,0.D0)
          B1D(I,J)=(0.D0,0.D0)
          B1E(I,J)=(0.D0,0.D0)
          BETA1U(I,J)=(0.D0,0.D0)
          BETA2U(I,J)=(0.D0,0.D0)
          BETA1D(I,J)=(0.D0,0.D0)
          BETA2D(I,J)=(0.D0,0.D0)
          BETA1E(I,J)=(0.D0,0.D0)
          BETA2E(I,J)=(0.D0,0.D0)
          B2YMU(I,J)=(0.D0,0.D0)
          B2YMD(I,J)=(0.D0,0.D0)
          B2YME(I,J)=(0.D0,0.D0)
        END DO
      END DO
      BGRKMU=(0.D0,0.D0)
      BETA2GRKMU=(0.D0,0.D0)
      BETA1VU=(0.D0,0.D0)
      BETA1VD=(0.D0,0.D0)
      BETA2VU=(0.D0,0.D0)
      BETA2VD=(0.D0,0.D0)
      BETALAM1=(0.D0,0.D0)
      BETALAM2=(0.D0,0.D0)
      BETAVEV1=(0.D0,0.D0)
      BETAVEV2=(0.D0,0.D0)
!
      PI=4.D0*ATAN(1.D0)
      Q=SSQSTEP
      IF(Q.LT.1.D0)THEN
        WRITE(*,*)'ERROR IN Q: ',Q
        STOP
      END IF
      NG=3.D0
      NU=3
      ND=3
      NE=3
      NNU=3
!
      NSQ=3
      IF((Q-QTHQL(1)).LT.-EPS.OR.
     $       (ABS(Q-QTHQL(1)).LT.ABS(EPS).AND.EPS.LT.0))NSQ=0
      NSU=3
      IF((Q-QTHUR(1)).LT.-EPS.OR.
     $       (ABS(Q-QTHUR(1)).LT.ABS(EPS).AND.EPS.LT.0))NSU=0
      NSD=3
      IF((Q-QTHDR(1)).LT.-EPS.OR.
     $       (ABS(Q-QTHDR(1)).LT.ABS(EPS).AND.EPS.LT.0))NSD=0
      NSL=3
      IF((Q-QTHLL(1)).LT.-EPS.OR.
     $       (ABS(Q-QTHLL(1)).LT.ABS(EPS).AND.EPS.LT.0))NSL=0
      NSE=3
      IF((Q-QTHER(1)).LT.-EPS.OR.
     $       (ABS(Q-QTHER(1)).LT.ABS(EPS).AND.EPS.LT.0))NSE=0
!
      IF ((Q-QNSH).GT.EPS.OR.
     $         (ABS(Q-QNSH).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSH=2
      ELSE
        NSH=0
      END IF
      IF ((Q-QNSG).GT.EPS.OR.
     $         (ABS(Q-QNSG).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSG=1
      ELSE
        NSG=0
      END IF
      IF ((Q-QNH).GT.EPS.OR.
     $         (ABS(Q-QNH).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NH=2
      ELSE
        NH=1
      END IF
      IF ((Q-QTHSB).GT.EPS.OR.
     $         (ABS(Q-QTHSB).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        THSB=1
      ELSE
        THSB=0
      END IF
      IF ((Q-QTHSW).GT.EPS.OR.
     $         (ABS(Q-QTHSW).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSW=1
        THSW=1
      ELSE
        NSW=0
        THSW=0
      END IF
      THLH=1
      THHH=NH/2 !This works so long as THHH is an integer variable
      IF(THHH.NE.1)THEN !Perform check
        IF(THHH.NE.0)WRITE(*,*)'ERROR IN THHH'
      END IF
      THSH=NSH/2
      DO I=1,3
        THSQ(I)=1
        IF((Q-QTHQL(1)).LT.-EPS.OR.
     $         (ABS(Q-QTHQL(1)).LT.ABS(EPS).AND.EPS.LT.0))THSQ(I)=0
        THSU(I)=1
        IF((Q-QTHUR(1)).LT.-EPS.OR.
     $         (ABS(Q-QTHUR(1)).LT.ABS(EPS).AND.EPS.LT.0))THSU(I)=0
        THSD(I)=1
        IF((Q-QTHDR(1)).LT.-EPS.OR.
     $         (ABS(Q-QTHDR(1)).LT.ABS(EPS).AND.EPS.LT.0))THSD(I)=0
        THSL(I)=NSL/3
        THSE(I)=NSE/3
      END DO
      THGL=NSG
!
!Convert input into 3x3 matrices.
!
      DO I=1,3
        DO J=1,3
          FU(I,J)=G(3+(I-1)*3+J)
          FD(I,J)=G(12+(I-1)*3+J)
          FE(I,J)=G(21+(I-1)*3+J)
!
!MV notation Yukawas
!
          YU(I,J)=G(3+(J-1)*3+I)
          YD(I,J)=G(12+(J-1)*3+I)
          YE(I,J)=G(21+(J-1)*3+I)
!
          LU(I,J)=G(33+(I-1)*3+J)
          LD(I,J)=G(42+(I-1)*3+J)
          LE(I,J)=G(51+(I-1)*3+J)
          GTPQ(I,J)=G(62+(I-1)*3+J)
          GTPL(I,J)=G(71+(I-1)*3+J)
          GTPU(I,J)=G(80+(I-1)*3+J)
          GTPD(I,J)=G(89+(I-1)*3+J)
          GTPE(I,J)=G(98+(I-1)*3+J)
          GTQ(I,J)=G(109+(I-1)*3+J)
          GTL(I,J)=G(118+(I-1)*3+J)
          GTSQ(I,J)=G(129+(I-1)*3+J)
          GTSU(I,J)=G(138+(I-1)*3+J)
          GTSD(I,J)=G(147+(I-1)*3+J)
          FTUQ(I,J)=G(156+(I-1)*3+J)
          FTDQ(I,J)=G(165+(I-1)*3+J)
          FTEL(I,J)=G(174+(I-1)*3+J)
          FTUU(I,J)=G(183+(I-1)*3+J)
          FTDD(I,J)=G(192+(I-1)*3+J)
          FTEE(I,J)=G(201+(I-1)*3+J)
        END DO
      END DO
!
      MVMU=G(31)
      MVMUM=G(215)
!
!SET THE TILDE TERMS TO THEIR SM COUNTERPARTS.
!
      DO I=1,3
        DO J=1,3
          GTPQ(I,J)=G(1)*SQRT(3.D0/5.D0)*ID(I,J)
          G(62+(I-1)*3+J)=GTPQ(I,J)
          GTPL(I,J)=G(1)*SQRT(3.D0/5.D0)*ID(I,J)
          G(71+(I-1)*3+J)=GTPL(I,J)
          GTPU(I,J)=G(1)*SQRT(3.D0/5.D0)*ID(I,J)
          G(80+(I-1)*3+J)=GTPU(I,J)
          GTPD(I,J)=G(1)*SQRT(3.D0/5.D0)*ID(I,J)
          G(89+(I-1)*3+J)=GTPD(I,J)
          GTPE(I,J)=G(1)*SQRT(3.D0/5.D0)*ID(I,J)
          G(98+(I-1)*3+J)=GTPE(I,J)
          GTQ(I,J)=G(2)*ID(I,J)
          G(109+(I-1)*3+J)=GTQ(I,J)
          GTL(I,J)=G(2)*ID(I,J)
          G(118+(I-1)*3+J)=GTL(I,J)
          GTSQ(I,J)=G(3)*ID(I,J)
          G(129+(I-1)*3+J)=GTSQ(I,J)
          GTSU(I,J)=G(3)*ID(I,J)
          G(138+(I-1)*3+J)=GTSU(I,J)
          GTSD(I,J)=G(3)*ID(I,J)
          G(147+(I-1)*3+J)=GTSD(I,J)
          FTUQ(I,J)=YU(J,I)
          G(156+(I-1)*3+J)=FTUQ(I,J)
          FTDQ(I,J)=YD(J,I)
          G(165+(I-1)*3+J)=FTDQ(I,J)
          FTEL(I,J)=YE(J,I)
          G(174+(I-1)*3+J)=FTEL(I,J)
          FTUU(I,J)=YU(J,I)
          G(183+(I-1)*3+J)=FTUU(I,J)
          FTDD(I,J)=YD(J,I)
          G(192+(I-1)*3+J)=FTDD(I,J)
          FTEE(I,J)=YE(J,I)
          G(201+(I-1)*3+J)=FTEE(I,J)
        END DO
      END DO
      G(108)=SQRT(3.D0/5.D0)*G(1)
      G(109)=SQRT(3.D0/5.D0)*G(1)
      G(128)=G(2)
      G(129)=G(2)
      G(211)=(0.D0,0.D0)
      G(212)=(0.D0,0.D0)
      G(213)=(0.D0,0.D0)
      G(214)=(0.D0,0.D0)
!
!The separated out contributions are from PRD 49 4882 (1194),
!Castano,Piard,Ramond
!
      B1LP(1)=2.D0/5.D0*(17.D0/12.D0*DBLE(NU)+5.D0/12.D0*DBLE(ND)
     $        +5.D0/4.D0*DBLE(NE)+1.D0/4.D0*DBLE(NNU))+1.D0/30.D0
     $        *DBLE(NSQ)+4.D0/15.D0*DBLE(NSU)+1.D0/15.D0*DBLE(NSD)
     $        +1.D0/10.D0*DBLE(NSL)+1.D0/5.D0*DBLE(NSE)
     $        +1.D0/5.D0*DBLE(NSH)+1.D0/10.D0*DBLE(NH)
      B1LP(2)=-22.D0/3.D0+1.D0/2.D0*(DBLE(NU)+DBLE(ND))+1.D0/6.D0
     $        *(DBLE(NE)+DBLE(NNU))+1.D0/2.D0*DBLE(NSQ)
     $        +1.D0/6.D0*DBLE(NSL)+1.D0/3.D0*DBLE(NSH)+1.D0/6.D0
     $        *DBLE(NH)+4.D0/3.D0*DBLE(NSW)
      B1LP(3)=-11.D0+2.D0/3.D0*(DBLE(NU)+DBLE(ND))+1.D0/3.D0
     $        *DBLE(NSQ)+1.D0/6.D0*DBLE(NSU)+1.D0/6.D0*DBLE(NSD)
     $        +2.D0*DBLE(NSG)
      IF(THHH.EQ.0)THEN
        B2LPSM(1,1)=-(-NG*19.D0/15.D0-9.D0/50.D0)
        B2LPSM(1,2)=-(-NG*3.D0/5.D0-9.D0/10.D0)
        B2LPSM(1,3)=-(-NG*44.D0/15.D0)
        B2LPSM(2,1)=-(-NG*1.D0/5.D0-3.D0/10.D0)
        B2LPSM(2,2)=-(136.D0/3.D0-NG*49.D0/3.D0-13.D0/6.D0)
        B2LPSM(2,3)=-(-NG*4.D0)
        B2LPSM(3,1)=-(-NG*11.D0/30.D0)
        B2LPSM(3,2)=-(-NG*3.D0/2.D0)
        B2LPSM(3,3)=-(102.D0-NG*76.D0/3.D0)
        CSM(1,1)=1.7D0
        CSM(1,2)=.5D0
        CSM(1,3)=1.5D0
        CSM(2,1)=1.5D0
        CSM(2,2)=1.5D0
        CSM(2,3)=.5D0
        CSM(3,1)=2.D0
        CSM(3,2)=2.D0
        CSM(3,3)=0.D0
      END IF
!
!I need many variations on the 3x3 matrices.
!
      CALL CDAGGER(FTUQ,FTUQD)
      CALL CDAGGER(FTUU,FTUUD)
      CALL CDAGGER(FTDD,FTDDD)
      CALL CDAGGER(FTDQ,FTDQD)
      CALL CDAGGER(FTEL,FTELD)
      CALL CDAGGER(FTEE,FTEED)
      DO I=1,3
        DO J=1,3
          GTPUS(I,J)=CONJG(GTPU(I,J))
          GTSUS(I,J)=CONJG(GTSU(I,J))
          GTQS(I,J)=CONJG(GTQ(I,J))
          GTLS(I,J)=CONJG(GTL(I,J))
          GTPQS(I,J)=CONJG(GTPQ(I,J))
          GTPLS(I,J)=CONJG(GTPL(I,J))
          GTSQS(I,J)=CONJG(GTSQ(I,J))
          GTPDS(I,J)=CONJG(GTPD(I,J))
          GTPES(I,J)=CONJG(GTPE(I,J))
          GTSDS(I,J)=CONJG(GTSD(I,J))
          GTPUT(I,J)=GTPU(J,I)
          GTQT(I,J)=GTQ(J,I)
          GTLT(I,J)=GTL(J,I)
          GTPQT(I,J)=GTPQ(J,I)
          GTPLT(I,J)=GTPL(J,I)
          GTPDT(I,J)=GTPD(J,I)
          GTPET(I,J)=GTPE(J,I)
          GTSUT(I,J)=GTSU(J,I)
          GTSQT(I,J)=GTSQ(J,I)
          GTSDT(I,J)=GTSD(J,I)
        END DO
      END DO
!
!Now all the matrix multiples
!
      DO I=1,3
        DO J=1,3
          YUDYU(I,J)=CMATMUL(1,YU,YU,I,J)
          YDDYD(I,J)=CMATMUL(1,YD,YD,I,J)
          YEDYE(I,J)=CMATMUL(1,YE,YE,I,J)
          YUYUD(I,J)=CMATMUL(2,YU,YU,I,J)
          YDYDD(I,J)=CMATMUL(2,YD,YD,I,J)
          YEYED(I,J)=CMATMUL(2,YE,YE,I,J)
          YUYDD(I,J)=CMATMUL(2,YU,YD,I,J)
          YDYUD(I,J)=CMATMUL(2,YD,YU,I,J)
!
          FUFUD(I,J)=CMATMUL(2,FU,FU,I,J)
          FDFDD(I,J)=CMATMUL(2,FD,FD,I,J)
          FEFED(I,J)=CMATMUL(2,FE,FE,I,J)
          FUDFU(I,J)=CMATMUL(1,FU,FU,I,J)
          FDDFD(I,J)=CMATMUL(1,FD,FD,I,J)
          FEDFE(I,J)=CMATMUL(1,FE,FE,I,J)
          LULUD(I,J)=CMATMUL(2,LU,LU,I,J)
          LDLDD(I,J)=CMATMUL(2,LD,LD,I,J)
          LELED(I,J)=CMATMUL(2,LE,LE,I,J)
        END DO
      END DO
!
      IF(THHH.EQ.0)THEN
        MSGTPHUSQ=CONJG(G(211))*G(211)
        MCGTPHDSQ=CONJG(G(212))*G(212)
        MSGTHUSQ=CONJG(G(213))*G(213)
        MCGTHDSQ=CONJG(G(214))*G(214)
      END IF
      MGTPHUSQ=CONJG(G(108))*G(108)
      MGTPHDSQ=CONJG(G(109))*G(109)
      MGTHUSQ=CONJG(G(128))*G(128)
      MGTHDSQ=CONJG(G(129))*G(129)
!
      TYUDYU=CTRACE(YUDYU)
      TYDDYD=CTRACE(YDDYD)
      TYEDYE=CTRACE(YEDYE)
!
      TFUDFU=CTRACE(FUDFU)
      TFDDFD=CTRACE(FDDFD)
      TFEDFE=CTRACE(FEDFE)
      TLUDLU=CTRACE(LULUD)
      TLDDLD=CTRACE(LDLDD)
      TLEDLE=CTRACE(LELED)
!
      DO I=1,3
        DO J=1,3
          FUFUDFU(I,J)=CMATMUL(0,FUFUD,FU,I,J)
          FUFUDFD(I,J)=CMATMUL(0,FUFUD,FD,I,J)
          FDFDDFU(I,J)=CMATMUL(0,FDFDD,FU,I,J)
          FDFDDFD(I,J)=CMATMUL(0,FDFDD,FD,I,J)
          FEFEDFE(I,J)=CMATMUL(0,FEFED,FE,I,J)
          LULUDLU(I,J)=CMATMUL(0,LULUD,LU,I,J)
          LULUDLD(I,J)=CMATMUL(0,LULUD,LD,I,J)
          LDLDDLU(I,J)=CMATMUL(0,LDLDD,LU,I,J)
          LDLDDLD(I,J)=CMATMUL(0,LDLDD,LD,I,J)
          LELEDLE(I,J)=CMATMUL(0,LELED,LE,I,J)
          SQFTUQDFTUQ(I,J)=CSFMUL(THSQ,FTUQD,FTUQ,I,J)
          SQFTDQDFTDQ(I,J)=CSFMUL(THSQ,FTDQD,FTDQ,I,J)
          SQGTQTGTQS(I,J)=CSFMUL(THSQ,GTQT,GTQS,I,J)
          SQGTQTFTUQ(I,J)=CSFMUL(THSQ,GTQT,FTUQ,I,J)
          SQGTQTFTDQ(I,J)=CSFMUL(THSQ,GTQT,FTDQ,I,J)
          SQGTPQTGTPQS(I,J)=CSFMUL(THSQ,GTPQT,GTPQS,I,J)
          SQGTPQTFTUQ(I,J)=CSFMUL(THSQ,GTPQT,FTUQ,I,J)
          SQGTPQTFTDQ(I,J)=CSFMUL(THSQ,GTPQT,FTDQ,I,J)
          SQGTSQTGTSQS(I,J)=CSFMUL(THSQ,GTSQT,GTSQS,I,J)
          SUFTUUFTUUD(I,J)=CSFMUL(THSU,FTUU,FTUUD,I,J)
          SUFTUUGTPUT(I,J)=CSFMUL(THSU,FTUU,GTPUT,I,J)
          SUGTPUSGTPUT(I,J)=CSFMUL(THSU,GTPUS,GTPUT,I,J)
          SUGTSUSGTSUT(I,J)=CSFMUL(THSU,GTSUS,GTSUT,I,J)
          SDFTDDFTDDD(I,J)=CSFMUL(THSD,FTDD,FTDDD,I,J)
          SDFTDDGTPDT(I,J)=CSFMUL(THSD,FTDD,GTPDT,I,J)
          SDGTPDSGTPDT(I,J)=CSFMUL(THSD,GTPDS,GTPDT,I,J)
          SDGTSDSGTSDT(I,J)=CSFMUL(THSD,GTSDS,GTSDT,I,J)
          SLFTELDFTEL(I,J)=CSFMUL(THSL,FTELD,FTEL,I,J)
          SLGTLTGTLS(I,J)=CSFMUL(THSL,GTLT,GTLS,I,J)
          SLGTLTFTEL(I,J)=CSFMUL(THSL,GTLT,FTEL,I,J)
          SLGTPLTGTPLS(I,J)=CSFMUL(THSL,GTPLT,GTPLS,I,J)
          SLGTPLTFTEL(I,J)=CSFMUL(THSL,GTPLT,FTEL,I,J)
          SEFTEEFTEED(I,J)=CSFMUL(THSE,FTEE,FTEED,I,J)
          SEFTEEGTPET(I,J)=CSFMUL(THSE,FTEE,GTPET,I,J)
          SEGTPESGTPET(I,J)=CSFMUL(THSE,GTPES,GTPET,I,J)
        END DO
      END DO
!
!These are the two loop terms. All in MV notation
!
      IF(SW2LP.EQ.1)THEN
!
        DO I=1,3
          DO J=1,3
            YUYUDYUYUD(I,J)=CMATMUL(0,YUYUD,YUYUD,I,J)
            YDYDDYDYDD(I,J)=CMATMUL(0,YDYDD,YDYDD,I,J)
            YEYEDYEYED(I,J)=CMATMUL(0,YEYED,YEYED,I,J)
            YUDYUYUDYU(I,J)=CMATMUL(0,YUDYU,YUDYU,I,J)
            YDDYDYDDYD(I,J)=CMATMUL(0,YDDYD,YDDYD,I,J)
            YEDYEYEDYE(I,J)=CMATMUL(0,YEDYE,YEDYE,I,J)
            YDDYDYUDYU(I,J)=CMATMUL(0,YDDYD,YUDYU,I,J)
            YUDYUYDDYD(I,J)=CMATMUL(0,YUDYU,YDDYD,I,J)
            YUYDDYDYUD(I,J)=CMATMUL(0,YUYDD,YDYUD,I,J)
          END DO
        END DO
!
        TYUDYUYUDYU=CTRACE(YUDYUYUDYU)
        TYUDYUYDDYD=CTRACE(YUDYUYDDYD)
        TYDDYDYUDYU=CTRACE(YDDYDYUDYU)
        TYDDYDYDDYD=CTRACE(YDDYDYDDYD)
        TYEDYEYEDYE=CTRACE(YEDYEYEDYE)
!
!These are SM terms for the two loop running below m_H
!I am going to use LYU for the MV notation SM Yukawa
!
        IF(THHH.EQ.0)THEN
          DO I=1,3
            DO J=1,3
              LYU(I,J)=LU(J,I)
              LYD(I,J)=LD(J,I)
              LYE(I,J)=LE(J,I)
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU(I,J)=CMATMUL(1,LYU,LYU,I,J)
              LYDDLYD(I,J)=CMATMUL(1,LYD,LYD,I,J)
              LYEDLYE(I,J)=CMATMUL(1,LYE,LYE,I,J)
            END DO
          END DO
          TLYUDLYU=CTRACE(LYUDLYU)
          TLYDDLYD=CTRACE(LYDDLYD)
          TLYEDLYE=CTRACE(LYEDLYE)
          DO I=1,3
            DO J=1,3
              LYUDLYU2(I,J)=CMATMUL(0,LYUDLYU,LYUDLYU,I,J)
              LYDDLYD2(I,J)=CMATMUL(0,LYDDLYD,LYDDLYD,I,J)
              LYEDLYE2(I,J)=CMATMUL(0,LYEDLYE,LYEDLYE,I,J)
              LYUDLYULYDDLYD(I,J)=CMATMUL(0,LYUDLYU,LYDDLYD,I,J)
              LYDDLYDLYUDLYU(I,J)=CMATMUL(0,LYDDLYD,LYUDLYU,I,J)
              DUMLUD1(I,J)=LYUDLYU(I,J)+LYDDLYD(I,J)
            END DO
          END DO
          TLYUDLYU2=CTRACE(LYUDLYU2)
          TLYDDLYD2=CTRACE(LYDDLYD2)
          TLYEDLYE2=CTRACE(LYEDLYE2)
          DO I=1,3
            DO J=1,3
              DUMLUD2(I,J)=CMATMUL(0,DUMLUD1,LYDDLYD,I,J)
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU3(I,J)=CMATMUL(0,LYUDLYU2,LYUDLYU,I,J)
              LYDDLYD3(I,J)=CMATMUL(0,LYDDLYD2,LYDDLYD,I,J)
              LYEDLYE3(I,J)=CMATMUL(0,LYEDLYE2,LYEDLYE,I,J)
              DUMLUD(I,J)=CMATMUL(0,LYUDLYU,DUMLUD2,I,J)
            END DO
          END DO
          TLYUDLYU3=CTRACE(LYUDLYU3)
          TLYDDLYD3=CTRACE(LYDDLYD3)
          TLYEDLYE3=CTRACE(LYEDLYE3)
          TLYUDLYULYDDLYD=CTRACE(LYUDLYULYDDLYD)
          TDUMLUD=CTRACE(DUMLUD)
!
          Y2=3.D0*TLYUDLYU+3.D0*TLYDDLYD+TLYEDLYE
          H=3.D0*TLYUDLYU2+3.D0*TLYDDLYD2+TLYEDLYE2
          Y4=(83.D0/40.D0*G(1)**2+27.D0/8.D0*G(2)**2
     $       +28.D0*G(3)**2)*TLYUDLYU+(-1.D0/40.D0*G(1)**2
     $       +27.D0/8.D0*G(2)**2+28.D0*G(3)**2)*TLYDDLYD
     $       +(93.D0/40.D0*G(1)**2+9.D0/8.D0*G(2)**2)*TLYEDLYE
          CHI4=9.D0/4.D0*(3.D0*TLYUDLYU2+3.D0*TLYDDLYD2+TLYEDLYE2
     $         -2.D0/3.D0*TLYUDLYULYDDLYD)
        END IF
      ELSE
        DO I=1,3
          DO J=1,3
            YUYUDYUYUD(I,J)=(0.D0,0.D0)
            YDYDDYDYDD(I,J)=(0.D0,0.D0)
            YEYEDYEYED(I,J)=(0.D0,0.D0)
            YUDYUYUDYU(I,J)=(0.D0,0.D0)
            YDDYDYDDYD(I,J)=(0.D0,0.D0)
            YEDYEYEDYE(I,J)=(0.D0,0.D0)
            YDDYDYUDYU(I,J)=(0.D0,0.D0)
            YUDYUYDDYD(I,J)=(0.D0,0.D0)
            YUYDDYDYUD(I,J)=(0.D0,0.D0)
          END DO
        END DO
        TYUDYUYUDYU=(0.D0,0.D0)
        TYUDYUYDDYD=(0.D0,0.D0)
        TYDDYDYUDYU=(0.D0,0.D0)
        TYDDYDYDDYD=(0.D0,0.D0)
        TYEDYEYEDYE=(0.D0,0.D0)
        IF(THHH.EQ.0)THEN
          DO I=1,3
            DO J=1,3
              LYU(I,J)=LU(J,I)
              LYD(I,J)=LD(J,I)
              LYE(I,J)=LE(J,I)
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU(I,J)=(0.D0,0.D0)
              LYDDLYD(I,J)=(0.D0,0.D0)
              LYEDLYE(I,J)=(0.D0,0.D0)
            END DO
          END DO
          TLYUDLYU=(0.D0,0.D0)
          TLYDDLYD=(0.D0,0.D0)
          TLYEDLYE=(0.D0,0.D0)
          DO I=1,3
            DO J=1,3
              LYUDLYU2(I,J)=(0.D0,0.D0)
              LYDDLYD2(I,J)=(0.D0,0.D0)
              LYEDLYE2(I,J)=(0.D0,0.D0)
              LYUDLYULYDDLYD(I,J)=(0.D0,0.D0)
              LYDDLYDLYUDLYU(I,J)=(0.D0,0.D0)
              DUMLUD1(I,J)=(0.D0,0.D0)
            END DO
          END DO
          TLYUDLYU2=(0.D0,0.D0)
          TLYDDLYD2=(0.D0,0.D0)
          TLYEDLYE2=(0.D0,0.D0)
          DO I=1,3
            DO J=1,3
              DUMLUD2(I,J)=(0.D0,0.D0)
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU3(I,J)=(0.D0,0.D0)
              LYDDLYD3(I,J)=(0.D0,0.D0)
              LYEDLYE3(I,J)=(0.D0,0.D0)
              DUMLUD(I,J)=(0.D0,0.D0)
            END DO
          END DO
          TLYUDLYU3=(0.D0,0.D0)
          TLYDDLYD3=(0.D0,0.D0)
          TLYEDLYE3=(0.D0,0.D0)
          TLYUDLYULYDDLYD=(0.D0,0.D0)
          TDUMLUD=(0.D0,0.D0)
!
          Y2=(0.D0,0.D0)
          H=(0.D0,0.D0)
          Y4=(0.D0,0.D0)
          CHI4=(0.D0,0.D0)
        END IF
      END IF
!
!The threshold gauge running with a change at m_H to SM
!
      DO I=1,3
        SUM=0.D0
        DO J=1,3
          IF(THHH.EQ.0)THEN
            SUM=SUM+B2LPSM(I,J)*G(J)**2
          ELSE
            SUM=SUM+B2LPM(I,J)*G(J)**2
          END IF
        END DO
        IF(THHH.EQ.0)THEN
          F(I)=G(I)**3/16.D0/PI**2*(B1LP(I)+DBLE(SW2LP)/16.D0/PI**2
     $         *(SUM-(CSM(I,1)*TLYUDLYU+CSM(I,2)*TLYDDLYD
     $         +CSM(I,3)*TLYEDLYE)))
        ELSE
          F(I)=G(I)**3/16.D0/PI**2*(B1LP(I)+DBLE(SW2LP)/16.D0/PI**2
     $         *(SUM-(CM(I,1)*TYUDYU+CM(I,2)*TYDDYD+CM(I,3)*TYEDYE)))
        END IF
      END DO
!
!Next the full Yukawas
!
      DO I=1,3
        DO J=1,3
          DUMU1(I,J)=THSH*SQFTUQDFTUQ(I,J)
     $               +4.D0/9.D0*THSB*SUGTPUSGTPUT(I,J)
     $               +4.D0/3.D0*THGL*SUGTSUSGTSUT(I,J)
          DUMU2(I,J)=2.D0*THSH*SUFTUUFTUUD(I,J)
     $               +2.D0*THSH*SDFTDDFTDDD(I,J)
     $               +3.D0*THSW*SQGTQTGTQS(I,J)
     $               +1.D0/9.D0*THSB*SQGTPQTGTPQS(I,J)
     $               +16.D0/3.D0*THGL*SQGTSQTGTSQS(I,J)
          DUMD1(I,J)=THSH*SQFTDQDFTDQ(I,J)
     $               +1.D0/9.D0*THSB*SDGTPDSGTPDT(I,J)
     $               +4.D0/3.D0*THGL*SDGTSDSGTSDT(I,J)
          DUMD2(I,J)=2.D0*THSH*SUFTUUFTUUD(I,J)
     $               +2.D0*THSH*SDFTDDFTDDD(I,J)
     $               +3.D0*THSW*SQGTQTGTQS(I,J)
     $               +1.D0/9.D0*THSB*SQGTPQTGTPQS(I,J)
     $               +16.D0/3.D0*THGL*SQGTSQTGTSQS(I,J)
          DUME1(I,J)=THSH*SLFTELDFTEL(I,J)+THSB*SEGTPESGTPET(I,J)
          DUME2(I,J)=2.D0*THSH*SEFTEEFTEED(I,J)
     $               +3.D0*THSW*SLGTLTGTLS(I,J)+THSB*SLGTPLTGTPLS(I,J)
!
          IF(SW2LP.EQ.1)THEN
            DUM2U1(I,J)=3.D0*YUYUDYUYUD(I,J)+YUYDDYDYUD(I,J)
            DUM2U2(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
            DUM2D1(I,J)=3.D0*YDYDDYDYDD(I,J)+YUYDDYDYUD(I,J)+
     $                  YEYEDYEYED(I,J)
            DUM2D2(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
            DUM2E1(I,J)=DUM2D1(I,J)
            DUM2E2(I,J)=DUM2U2(I,J)
          ELSE
            DUM2U1(I,J)=(0.D0,0.D0)
            DUM2U2(I,J)=(0.D0,0.D0)
            DUM2D1(I,J)=(0.D0,0.D0)
            DUM2D2(I,J)=(0.D0,0.D0)
            DUM2E1(I,J)=(0.D0,0.D0)
            DUM2E2(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
      DO I=1,3
        DO J=1,3
!
!Here are the two loop terms
!
          IF(SW2LP.EQ.1)THEN
            DUM2U(I,J)=-3.D0*CTRACE(DUM2U1)*ID(I,J)-YDDYD(I,J)
     $                 *CTRACE(DUM2U2)-9.D0*YUDYU(I,J)*TYUDYU
     $                 -4.D0*YUDYUYUDYU(I,J)-2.D0*YDDYDYDDYD(I,J)-2.D0
     $                 *YDDYDYUDYU(I,J)+(16.D0*G(3)**2+4.D0/5.D0
     $                 *G(1)**2)*TYUDYU*ID(I,J)+(6.D0*G(2)**2
     $                 +2.D0/5.D0*G(1)**2)*YUDYU(I,J)+2.D0/5.D0
     $                 *G(1)**2*YDDYD(I,J)+(-16.D0/9.D0*G(3)**4
     $                 +8.D0*G(3)**2*G(2)**2+136.D0/45.D0
     $                 *G(3)**2*G(1)**2+15.D0/2.D0*G(2)**4
     $                 +G(2)**2*G(1)**2+2743.D0/450.D0
     $                 *G(1)**4)*ID(I,J)
            DUM2D(I,J)=-3.D0*CTRACE(DUM2D1)*ID(I,J)-3.D0*YUDYU(I,J)
     $                 *TYUDYU-3.D0*YDDYD(I,J)*CTRACE(DUM2D2)
     $                 -4.D0*YDDYDYDDYD(I,J)-2.D0*YUDYUYUDYU(I,J)
     $                 -2.D0*YUDYUYDDYD(I,J)+(16.D0*G(3)**2
     $                 -2.D0/5.D0*G(1)**2)*TYDDYD*ID(I,J)+6.D0/5.D0
     $                 *G(1)**2*TYEDYE*ID(I,J)+4.D0/5.D0*G(1)**2
     $                 *YUDYU(I,J)+(6.D0*G(2)**2+4.D0/5.D0
     $                 *G(1)**2)*YDDYD(I,J)+(-16.D0/9.D0*G(3)**4
     $                 +8.D0*G(3)**2*G(2)**2+8.D0/9.D0
     $                 *G(3)**2*G(1)**2+15.D0/2.D0*G(2)**4
     $                 +G(2)**2*G(1)**2+287.D0/90.D0
     $                 *G(1)**4)*ID(I,J)
            DUM2E(I,J)=-3.D0*CTRACE(DUM2E1)*ID(I,J)-3.D0*YEDYE(I,J)
     $                 *CTRACE(DUM2E2)-4.D0*YEDYEYEDYE(I,J)+(16.D0
     $                 *G(3)**2-2.D0/5.D0*G(1)**2)*TYDDYD
     $                 *ID(I,J)+6.D0/5.D0*G(1)**2*TYEDYE*ID(I,J)
     $                 +6.D0*G(2)**2*YEDYE(I,J)+(15.D0/2.D0
     $                 *G(2)**4+9.D0/5.D0*G(2)**2*G(1)**2
     $                 +27.D0/2.D0*G(1)**4)*ID(I,J)
          ELSE
            DUM2U(I,J)=(0.D0,0.D0)
            DUM2D(I,J)=(0.D0,0.D0)
            DUM2E(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
!
!Now calculate the beta functions for the Yukawas
!
      DO I=1,3
        DO J=1,3
!
          B1U(I,J)=1.D0/2.D0*(3.D0*FUFUDFU(I,J)+FDFDDFU(I,J))
     $             +CMATMUL(0,FU,DUMU1,I,J)
     $             +1.D0/4.D0*CMATMUL(0,DUMU2,FU,I,J)
     $             +THSH*(-3.D0*THSW*CONJG(G(128))
     $              *SQGTQTFTUQ(I,J)+1.D0/3.D0*THSB
     $              *CONJG(G(108))*SQGTPQTFTUQ(I,J))
     $             -4.D0/3.D0*THSB*THSH*CONJG(G(108))
     $              *SUFTUUGTPUT(I,J)
     $             +FU(I,J)*3.D0*TFUDFU
     $             +1.D0/2.D0*THSH*FU(I,J)*(3.D0*THSW*MGTHUSQ
     $              +THSB*MGTPHUSQ)
     $             -FU(I,J)*(17.D0/20.D0*G(1)**2
     $              +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
!
          B1D(I,J)=1.D0/2.D0*(3.D0*FDFDDFD(I,J)+FUFUDFD(I,J))
     $             +CMATMUL(0,FD,DUMD1,I,J)
     $             +1.D0/4.D0*CMATMUL(0,DUMD2,FD,I,J)
     $             +THSH*(-3.D0*THSW*CONJG(G(129))
     $              *SQGTQTFTDQ(I,J)-1.D0/3.D0*THSB
     $              *CONJG(G(109))*SQGTPQTFTDQ(I,J))
     $             -2.D0/3.D0*THSB*THSH*CONJG(G(109))
     $              *SDFTDDGTPDT(I,J)
     $             +FD(I,J)*(3.D0*TFDDFD+TFEDFE)
     $             +1.D0/2.D0*THSH*FD(I,J)*(3.D0*THSW*MGTHDSQ
     $              +THSB*MGTPHDSQ)
     $             -FD(I,J)*(3.D0/12.D0*G(1)**2
     $              +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
!
          B1E(I,J)=3.D0/2.D0*FEFEDFE(I,J)
     $             +CMATMUL(0,FE,DUME1,I,J)
     $             +1.D0/4.D0*CMATMUL(0,DUME2,FE,I,J)
     $             +THSH*(-3.D0*THSW*CONJG(G(129))
     $              *SLGTLTFTEL(I,J)+THSB
     $              *CONJG(G(109))*SLGTPLTFTEL(I,J))
     $             -2.D0*THSB*THSH*CONJG(G(109))
     $              *SEFTEEGTPET(I,J)
     $             +FE(I,J)*(3.D0*TFDDFD+TFEDFE)
     $             +1.D0/2.D0*THSH*FE(I,J)*(3.D0*THSW*MGTHDSQ
     $              +THSB*MGTPHDSQ)
     $             -FE(I,J)*(9.D0/4.D0*G(1)**2
     $              +9.D0/4.D0*G(2)**2)
!
          IF(SW2LP.EQ.1)THEN
            B2YMU(I,J)=CMATMUL(0,YU,DUM2U,I,J)
            B2YMD(I,J)=CMATMUL(0,YD,DUM2D,I,J)
            B2YME(I,J)=CMATMUL(0,YE,DUM2E,I,J)
          END IF
        END DO
      END DO
      DO I=1,3
        DO J=1,3
!
!Convert into form readable by RKSTP. The transpose in BETA2 takes
!account of the differences in notation.
!
          F(3+(I-1)*3+J)=1.D0/16.D0/PI**2*B1U(I,J)
     $                   +1.D0/(16.D0*PI**2)**2*B2YMU(J,I)
          F(12+(I-1)*3+J)=1.D0/16.D0/PI**2*B1D(I,J)
     $                    +1.D0/(16.D0*PI**2)**2*B2YMD(J,I)
          F(21+(I-1)*3+J)=1.D0/16.D0/PI**2*B1E(I,J)
     $                    +1.D0/(16.D0*PI**2)**2*B2YME(J,I)
        END DO
      END DO
!
!The lambdas use the same dummy matrices. I will reuse the betas
!Only find the lambdas if we are below m_H
!
      IF(THHH.EQ.0)THEN
        DO I=1,3
          DO J=1,3
            B1U(I,J)=(0.D0,0.D0)
            B1D(I,J)=(0.D0,0.D0)
            B1E(I,J)=(0.D0,0.D0)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            IF(SW2LP.EQ.1)THEN
              DUM2U(I,J)=3.D0/2.D0*LYUDLYU2(I,J)
     $                  -LYUDLYULYDDLYD(I,J)
     $                  -1.D0/4.D0*LYDDLYDLYUDLYU(I,J)
     $                  +11.D0/4.D0*LYDDLYD2(I,J)
     $                  +Y2*(5.D0/4.D0*LYDDLYD(I,J)-9.D0/4.D0
     $                  *LYUDLYU(I,J))-CHI4*ID(I,J)+3.D0/2.D0
     $                  *G(62)**2*ID(I,J)-2.D0*G(62)*(3.D0
     $                  *LYUDLYU(I,J)+LYDDLYD(I,J))+(221.D0/80.D0
     $                  *G(1)**2+117.D0/16.D0*G(2)**2+20.D0
     $                  *G(3)**2)*LYUDLYU(I,J)-(17.D0/80.D0*G(1)**2
     $                  -27.D0/16.D0*G(2)**2+20.D0*G(3)**2)
     $                  *LYDDLYD(I,J)+Y4*ID(I,J)+((7.D0/150.D0
     $                  +2.D0/3.D0*NG)*G(1)**4-9.D0/20.D0*G(1)**2
     $                  *G(2)**2+19.D0/15.D0*G(1)**2*G(3)**2
     $                  -(101.D0/8.D0-2.D0*NG)*G(2)**4+9.D0*G(2)**2
     $                  *G(3)**2-(292.D0/3.D0-16.D0/3.D0*NG)
     $                  *G(3)**4)*ID(I,J)
              DUM2D(I,J)=3.D0/2.D0*LYDDLYD2(I,J)
     $                  -LYDDLYDLYUDLYU(I,J)
     $                  -1.D0/4.D0*LYUDLYULYDDLYD(I,J)
     $                  +11.D0/4.D0*LYUDLYU2(I,J)
     $                  +Y2*(5.D0/4.D0*LYUDLYU(I,J)-9.D0/4.D0
     $                  *LYDDLYD(I,J))-CHI4*ID(I,J)+3.D0/2.D0
     $                  *G(62)**2*ID(I,J)-2.D0*G(62)*(3.D0
     $                  *LYDDLYD(I,J)+LYUDLYU(I,J))+(161.D0/80.D0
     $                  *G(1)**2+117.D0/16.D0*G(2)**2+20.D0
     $                  *G(3)**2)*LYDDLYD(I,J)-(77.D0/80.D0*G(1)**2
     $                  -27.D0/16.D0*G(2)**2+20.D0*G(3)**2)
     $                  *LYUDLYU(I,J)+Y4*ID(I,J)+(-(37.D0/300.D0
     $                  -4.D0/15.D0*NG)*G(1)**4-27.D0/20.D0*G(1)**2
     $                  *G(2)**2+31.D0/15.D0*G(1)**2*G(3)**2
     $                  -(101.D0/8.D0-2.D0*NG)*G(2)**4+9.D0*G(2)**2
     $                  *G(3)**2-(292.D0/3.D0-16.D0/3.D0*NG)
     $                  *G(3)**4)*ID(I,J)
              DUM2E(I,J)=3.D0/2.D0*LYEDLYE2(I,J)
     $                  -Y2*9.D0/4.D0*LYEDLYE(I,J)-CHI4*ID(I,J)
     $                  +3.D0/2.D0*G(62)**2*ID(I,J)-6.D0*G(62)
     $                  *LYEDLYE(I,J)+(441.D0/80.D0*G(1)**2
     $                  +117.D0/16.D0*G(2)**2)*LYEDLYE(I,J)
     $                  +Y4*ID(I,J)+((21.D0/100.D0+8.D0/5.D0*NG)
     $                  *G(1)**4+27.D0/20.D0*G(1)**2*G(2)**2
     $                  -(101.D0/8.D0-2.D0*NG)*G(2)**4)*ID(I,J)
            ELSE
              DUM2U(I,J)=(0.D0,0.D0)
              DUM2D(I,J)=(0.D0,0.D0)
              DUM2E(I,J)=(0.D0,0.D0)
            END IF
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            B1U(I,J)=1.D0/2.D0*(3.D0*LULUDLU(I,J)+LDLDDLU(I,J)
     $                -4.D0*LDLDDLU(I,J))
     $               +CMATMUL(0,LU,DUMU1,I,J)
     $               +1.D0/4.D0*CMATMUL(0,DUMU2,LU,I,J)
     $               +THSH*(-3.D0*THSW*CONJG(G(213))
     $                *SQGTQTFTUQ(I,J)+1.D0/3.D0*THSB
     $                *CONJG(G(211))*SQGTPQTFTUQ(I,J))
     $               -4.D0/3.D0*THSB*THSH*CONJG(G(211))
     $                *SUFTUUGTPUT(I,J)
     $               +LU(I,J)*(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)
     $               +1.D0/2.D0*THSH*LU(I,J)*(3.D0*THSW
     $                *(MSGTHUSQ+MCGTHDSQ)+THSB
     $                *(MSGTPHUSQ+MCGTPHDSQ))
     $               -LU(I,J)*(17.D0/20.D0*G(1)**2
     $                +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
            B1D(I,J)=1.D0/2.D0*(3.D0*LDLDDLD(I,J)+LULUDLD(I,J)
     $                -4.D0*LULUDLD(I,J))
     $               +CMATMUL(0,LD,DUMD1,I,J)
     $               +1.D0/4.D0*CMATMUL(0,DUMD2,LD,I,J)
     $               +THSH*(-3.D0*THSW*CONJG(G(214))
     $                *SQGTQTFTDQ(I,J)-1.D0/3.D0*THSB
     $                *CONJG(G(212))*SQGTPQTFTDQ(I,J))
     $               -2.D0/3.D0*THSB*THSH*CONJG(G(212))
     $                *SDFTDDGTPDT(I,J)
     $               +LD(I,J)*(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)
     $               +1.D0/2.D0*THSH*LD(I,J)*(3.D0*THSW
     $                *(MSGTHUSQ+MCGTHDSQ)+THSB
     $                *(MSGTPHUSQ+MCGTPHDSQ))
     $               -LD(I,J)*(3.D0/12.D0*G(1)**2
     $                +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
            B1E(I,J)=3.D0/2.D0*LELEDLE(I,J)+CMATMUL(0,LE,DUME1,I,J)
     $               +1.D0/4.D0*CMATMUL(0,DUME2,LE,I,J)
     $               +THSH*(-3.D0*THSW*CONJG(G(214))
     $                *SLGTLTFTEL(I,J)+THSB
     $                *CONJG(G(212))*SLGTPLTFTEL(I,J))
     $               -2.D0*THSB*THSH*CONJG(G(212))
     $                *SEFTEEGTPET(I,J)
     $               +LE(I,J)*(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)
     $               +1.D0/2.D0*THSH*LE(I,J)*(3.D0*THSW
     $                *(MSGTHUSQ+MCGTHDSQ)+THSB
     $                *(MSGTPHUSQ+MCGTPHDSQ))
     $               -LE(I,J)*(9.D0/4.D0*G(1)**2
     $                +9.D0/4.D0*G(2)**2)
!
            IF(SW2LP.EQ.1)THEN
              BETA2U(I,J)=CMATMUL(0,LYU,DUM2U,I,J)
              BETA2D(I,J)=CMATMUL(0,LYD,DUM2D,I,J)
              BETA2E(I,J)=CMATMUL(0,LYE,DUM2E,I,J)
            END IF
          END DO
        END DO
        DO I=1,3
          DO J=1,3
!
!Convert into form readable by RKSTP.
!
            F(33+(I-1)*3+J)=1.D0/16.D0/PI**2*B1U(I,J)
     $              +1.D0/(16.D0*PI**2)**2*BETA2U(J,I)
            F(42+(I-1)*3+J)=1.D0/16.D0/PI**2*B1D(I,J)
     $              +1.D0/(16.D0*PI**2)**2*BETA2D(J,I)
            F(51+(I-1)*3+J)=1.D0/16.D0/PI**2*B1E(I,J)
     $              +1.D0/(16.D0*PI**2)**2*BETA2E(J,I)
          END DO
        END DO
      END IF
!
!Next I am going to work out the gaugino terms, \mu and M_{1,2,3}
!and the running of B in MV notation
!
      DO I=1,3
        DO J=1,3
          DUMGRKMU1(I,J)=3.D0*SUFTUUFTUUD(I,J)+3.D0*SDFTDDFTDDD(I,J)
     $                   +SEFTEEFTEED(I,J)+3.D0*SQFTUQDFTUQ(I,J)
     $                   +3.D0*SQFTDQDFTDQ(I,J)+SLFTELDFTEL(I,J)
!
          IF(SW2LP.EQ.1)THEN
            DUM2GRKMU(I,J)=3.D0*YUYUDYUYUD(I,J)+3.D0*YDYDDYDYDD(I,J)
     $                     +2.D0*YUYDDYDYUD(I,J)+YEYEDYEYED(I,J)
          ELSE
            DUM2GRKMU(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
      TDUMGRKMU=CTRACE(DUMGRKMU1)
!
      BGRKMU=1.D0/2.D0*G(31)*THSH*TDUMGRKMU
     $       +1.D0/4.D0*G(31)*THSH*(3.D0*THSW*MGTHUSQ+THSB*MGTPHUSQ
     $        +3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $       -G(31)*(9.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
!
      IF(SW2LP.EQ.1)THEN
        BETA2GRKMU=MVMU*(-3.D0*CTRACE(DUM2GRKMU)+(16.D0*G(3)**2
     $             +4.D0/5.D0*G(1)**2)*TYUDYU+(16.D0*G(3)**2
     $             -2.D0/5.D0*G(1)**2)*TYDDYD+6.D0/5.D0*G(1)**2
     $             *TYEDYE+15.D0/2.D0*G(2)**4+9.D0/5.D0*G(1)**2
     $             *G(2)**2+207.D0/50.D0*G(1)**4)
      END IF
!
!The RKSTP compatible derivative is:
!
        F(31)=1.D0/16.D0/PI**2*BGRKMU+1.D0/(16.D0*PI**2)**2
     $         *BETA2GRKMU
!
!V_U and V_D - from PRD 49,4882 (1994)
!
      BETA1VU=3.D0/4.D0*(1.D0/5.D0*G(1)**2+G(2)**2)-3.D0*TYUDYU
      BETA1VD=3.D0/4.D0*(1.D0/5.D0*G(1)**2+G(2)**2)-3.D0*TYDDYD
     $        -TYEDYE
      IF(SW2LP.EQ.1)THEN
        BETA2VU=3.D0/4.D0*(3.D0*TYUDYUYUDYU+3.D0*TYUDYUYDDYD)-(19.D0
     $          /10.D0*G(1)**2+9.D0/2.D0*G(2)**2+20.D0*G(3)**2)
     $          *TYUDYU-(279.D0/800.D0+1803.D0/1600.D0*3.D0)*G(1)**4
     $          -(207.D0/32.D0+357.D0/64.D0*3.D0)*G(2)**4-(27.D0/80.D0
     $          +9.D0/80.D0*3.D0)*G(1)**2*G(2)**2
        BETA2VD=3.D0/4.D0*(3.D0*TYDDYDYDDYD+3.D0*TYDDYDYUDYU
     $          +TYEDYEYEDYE)-(2.D0/5.D0*G(1)**2+9.D0/2.D0*G(2)**2
     $          +20.D0*G(3)**2)*TYDDYD-(9.D0/5.D0*G(1)**2+3.D0/2.D0
     $          *G(2)**2)*TYEDYE-(279.D0/800.D0+1803.D0/1600.D0*3.D0)
     $          *G(1)**4-(207.D0/32.D0+357.D0/64.D0*3.D0)*G(2)**4
     $          -(27.D0/80.D0+9.D0/80.D0*3.D0)*G(1)**2*G(2)**2
      END IF
!
      F(32)=G(32)*(1.D0/16.D0/PI**2*BETA1VU+1.D0/(16.D0*PI**2)**2
     $       *BETA2VU)
      F(33)=G(33)*(1.D0/16.D0/PI**2*BETA1VD+1.D0/(16.D0*PI**2)**2
     $       *BETA2VD)
!
      IF(THHH.EQ.0)THEN
!
!Finally we have the running of the Higgs Quartic Coupling and SM VEV.
!Programmed here is the MS-bar running. It therefore needs the MS-bar
!gauge and Yukawas.
!The gauge couplings and Yukawas needed to be converted to MS-bar using
!the Martin and Vaughn conversion in hep-ph/9308222.
!The following is after the conversion, so all Yukawas and Gauge
!couplings are still in the DR-bar scheme.
!
        BETALAM1=12*G(62)**2-(9.D0/5.D0*G(1)**2+9.D0*G(2)**2)
     $           *G(62)+9.D0/4.D0*(3.D0/25.D0*G(1)**4+2.D0/5.D0
     $           *G(1)**2*G(2)**2+G(2)**4)+4.D0*Y2*G(62)-4*H
        IF(SW2LP.EQ.1)THEN
          BETALAM2=-78.D0*G(62)**3+18.D0*(3.D0/5.D0*G(1)**2
     $            +3.D0*G(2)**2)*G(62)**2-((265.D0/8.D0-10*NG)
     $            *G(2)**4-117.D0/20.D0*G(1)**2*G(2)**2
     $            -9.D0/25.D0*(229.D0/24.D0+50.D0/9.D0*NG)*G(1)**4)
     $            *G(62)+(473.D0/8.D0-8.D0*NG)*G(2)**6-3.D0/5.D0
     $            *(121.D0/24.D0+8.D0/3.D0*NG)*G(1)**2*G(2)**4
     $            -9.D0/25.D0*(239.D0/24.D0+40.D0/9.D0*NG)
     $            *G(1)**4*G(2)**2-27.D0/125.D0*(59.D0/24.D0
     $            +40.D0/9.D0*NG)*G(1)**6+(-14.D0/5.D0*G(1)**2
     $            +18.D0*G(2)**2-128.D0*G(3)**2)*TLYUDLYU2
     $            +(34.D0/5.D0*G(1)**2+18.D0*G(2)**2-128.D0
     $            *G(3)**2)*TLYDDLYD2+(-42.D0/5.D0*G(1)**2
     $            +6.D0*G(2)**2)*TLYEDLYE2-3.D0/2.D0*G(2)**4
     $            *Y2+G(62)*((83.D0/10.D0*G(1)**2+27.D0/2.D0
     $            *G(2)**2+112.D0*G(3)**2)*TLYUDLYU+(-1.D0/10.D0
     $            *G(1)**2+27.D0/2.D0*G(2)**2+112.D0*G(3)**2)
     $            *TLYDDLYD+(93.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
     $            *TLYEDLYE)+3.D0/5.D0*G(1)**2*((-57.D0/10.D0
     $            *G(1)**2+21.D0*G(2)**2)*TLYUDLYU+(3.D0/2.D0
     $            *G(1)**2+9.D0*G(2)**2)*TLYDDLYD+(-15.D0/2.D0
     $            *G(1)**2+11.D0*G(2)**2)*TLYEDLYE)-24.D0
     $            *G(62)**2*Y2-G(62)*H+6.D0*G(62)
     $            *TLYUDLYULYDDLYD+20.D0*(3.D0*TLYUDLYU3+3.D0*TLYDDLYD3
     $            +TLYEDLYE3)-12.D0*TDUMLUD
        END IF
!
        F(62)=1.D0/(16.D0*PI**2)*BETALAM1+1.D0/(16.D0*PI**2)**2
     $         *BETALAM2
!
!Calculate the betas for the standard model vev.
!As with lambda this is the MS-bar running with DR-bar inputs except
!v and lambda
!
      BETAVEV1=9.D0/4.D0*(1.D0/5.D0*G(1)**2+G(2)**2)-Y2
      IF(SW2LP.EQ.1)THEN
        BETAVEV2=-3.D0/2.D0*G(62)**2-(83.D0/40.D0*G(1)**2+27.D0/8.D0
     $           *G(2)**2+28.D0*G(3)**2)*TYUDYU-(-1.D0/40.D0*G(1)**2
     $           +27.D0/8.D0*G(2)**2+28.D0*G(3)**2)*TYDDYD
     $           -(93.D0/40.D0*G(1)**2+9.D0/8.D0*G(2)**2)*TYEDYE+CHI4
     $           -27.D0/80.D0*G(1)**2*G(2)**2-(93.D0/800.D0+1.D0/2.D0
     $           *NG)*G(1)**4+(463.D0/32.D0-5.D0/2.D0*NG)*G(2)**4
      END IF
!
      F(61)=G(61)*(1.D0/(16.D0*PI**2)*BETAVEV1+1.D0/(16.D0*PI**2)**2
     $      *BETAVEV2)
!
      END IF
!
!Finally, the MSSM mu parameter
!
      DO I=1,3
        DO J=1,3
          DUMGRKMU1(I,J)=3.D0*YUYUD(I,J)+3.D0*YDYDD(I,J)+YEYED(I,J)
          IF(SW2LP.EQ.1)THEN
            DUM2GRKMU(I,J)=3.D0*YUYUDYUYUD(I,J)+3.D0*YDYDDYDYDD(I,J)
     $                     +2.D0*YUYDDYDYUD(I,J)+YEYEDYEYED(I,J)
          ELSE
            DUM2GRKMU(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
!
      BGRKMU=MVMUM*(CTRACE(DUMGRKMU1)-3.D0*G(2)**2-3.D0/5.D0
     $           *G(1)**2)
!
      IF(SW2LP.EQ.1)THEN
        BETA2GRKMU=MVMUM*(-3.D0*CTRACE(DUM2GRKMU)+(16.D0*G(3)**2+4.D0
     $             /5.D0*G(1)**2)*TYUDYU+(16.D0*G(3)**2-2.D0/5.D0
     $             *G(1)**2)*TYDDYD+6.D0/5.D0*G(1)**2*TYEDYE
     $             +15.D0/2.D0*G(2)**4+9.D0/5.D0*G(1)**2*G(2)**2
     $             +207.D0/50.D0*G(1)**4)
      END IF
!
!The RKSTP compatible derivative is.
!
        F(215)=1.D0/16.D0/PI**2*BGRKMU+1.D0/(16.D0*PI**2)**2
     $         *BETA2GRKMU
!
      RETURN
      END
