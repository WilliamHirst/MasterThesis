C-------------------------------------------------------
C#######################################################
c ==========================================================================
C   
C   ISAJET - CompHEP interface
C
      SUBROUTINE ISACHP  
C-----------------------------------------------------------------------
C          A primitive interface between ISAJET and CompHEP generated code.
C          Fills the CompHEP common blocks from ISAJET common blocks.
C-----------------------------------------------------------------------
      IMPLICIT NONE

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

      INTEGER NOUT
      PARAMETER (NOUT=33)
      INTEGER IDOUT(NOUT)
      REAL AMASS,AMPL
      REAL AMI,SUMGAM,SUMMJ,WIDMX
      REAL QSUSY,ASMB,MBMB,ASMT,MTMT,SUALFS,GG
      DOUBLE PRECISION SSMQCD
      INTEGER I,J,K,IFL1,IFL2,IFL3,JSPIN,INDEX,IITEST
C
      DATA IDOUT/
     $IDTP,ISGL,ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1,ISUPR,ISDNR,
     $ISSTR,ISCHR,ISBT2,ISTP2,ISEL,ISMUL,ISTAU1,ISNEL,ISNML,ISNTL,
     $ISER,ISMUR,ISTAU2,ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,
     $ISHL,ISHH,ISHA,ISHC/
      
      REAL FTAMZ1,FBOMZ1,FTOMZ1,VUMZ1,VDMZ1       
      

************************
c   CompHEP common blocks
************************
      Real*8 A
      COMMON/VARS/A(1800)
      CHARACTER*6 xxxxx(146)
      DATA xxxxx/
     >'EE',	'SW',  	's12',	's23',	's13',           !1
     >'MZ',	'MH3',	'Mm',	'Mt',	'Mc',		 !2 
     >'Ms',	'Mtop',	'Mb',	'wtop',	'wZ',		 !3 
     >'wW',	'sa',	'MU',	'Zp11',	'Zp12', 	 !4 
     >'Zp21',	'Zp22',	'Zm11',	'Zm12',	'Zm21', 	 !5 
     >'Zm22',	'Zn11',	'Zn12',	'Zn13',	'Zn14', 	 !6 
     >'Zn21',	'Zn22',	'Zn23',	'Zn24',	'Zn31', 	 !7 
     >'Zn32',	'Zn33',	'Zn34',	'Zn41',	'Zn42', 	 !8 
     >'Zn43',	'Zn44',	'Atau',	'Zl33',	'Zl36', 	 !9 
     >'Zl63',	'Zl66',	'Zl11',	'Zl14',	'Zl41', 	 !10 
     >'Zl44',	'Zl22',	'Zl25',	'Zl52',	'Zl55', 	 !11 
     >'Atop',	'Abot',	'Zu33',	'Zu36',	'Zu63', 	 !12 
     >'Zu66',	'Zu11',	'Zu14',	'Zu41',	'Zu44', 	 !13 
     >'Zu22',	'Zu25',	'Zu52',	'Zu55',	'Zd33', 	 !14 
     >'Zd36',	'Zd63',	'Zd66',	'Zd11',	'Zd14', 	 !15 
     >'Zd41',	'Zd44',	'Zd22',	'Zd25',	'Zd52', 	 !16 
     >'Zd55',	'Mh',	'wh',	'MHH',	'wHh',		 !17 
     >'wH3',	'MHc',	'wHc',	'MC1',	'wC1',		 !18 
     >'MC2',	'wC2',	'MNE1',	'wNE1',	'MNE2', 	 !19 
     >'wNE2',	'MNE3',	'wNE3',	'MNE4',	'wNE4', 	 !20 
     >'wSG',	'wSe1',	'wSe2',	'wSmu1','wSmu2',	 !21 
     >'MStau1',	'wStau1','MStau2','wStau2','wSne',	 !22 
     >'wSnmu',	'MSntau','wSntau','wSu1','wSu2',	 !23 
     >'wSd1',	'wSd2',	'wSc1',	'wSc2',	'wSs1', 	 !24 
     >'wSs2',	'MStop1','wStop1','MStop2','wStop2'	 !25 
     >,'MSbot1','wSbot1','MSbot2','wSbot2','tb',	 !26
     >'MSG',	'MSe1',	'MSe2',	'MSmu1','MSmu2',	 !27
     >'MSne',	'MSnmu','MSu1',	'MSu2',	'MSd1', 	 !28 
     >'MSd2',	'MSc1',	'MSc2',	'MSs1',	'MSs2', 	 !29 
     >'GG'/						 !30 
************************
c      Bases common blocks
************************
      INTEGER NPRINT
      common/printlevel/NPRINT
******************
c   Local variables
******************
      INTEGER ii
      
      
      REAL TANB,MT,PI,GPX,THX,THY,THETAX,THETAY,GAMTOT
      PI=4.*ATAN(1.)

CsB   sin of the Weinberg angle
      GPX=SQRT(.6)*GSS(1)
      IF((GSS(2)**2+GPX**2).eq.0) THEN 
        SIN2W=SN2THW
      ELSE
        SIN2W=GPX**2/(GSS(2)**2+GPX**2)
      ENDIF

C   Set constants & EW parameters
C   This *must* be done after defining SIN2W above
      CALL SETCON
      CALL SETW
      
CsB   theta_x and theta_y in the chargino mixing matrices
      CALL THETAXY(THX,THY)
      THETAX = THX
      THETAY = THY

C     CompHEP = ISAJET              ! name
      A(1) = SQRT(4.D0*PI/128.D0)   ! SQRT(4*Pi/128)
      A(2) = SQRT(SIN2W)            ! sin of the Weinberg angle
      A(3) = 0.221                  ! CKM matrix element (= 0.221 PDG'94)
      A(4) = 0.040                  ! CKM matrix element (= 0.040 PDG'94)
      A(5) = .0035                  ! CKM matrix element (= .0035 PDG'94)
     
      A(6) = AMZ                    ! Z boson mass
      A(7) = MSS(31)                ! mass of CP-odd Higgs
      A(8) = AMMU                   ! muon mass
C---------------------------------------------------C
      A(10) = AMCH                  ! charm mass
    
      A(11) = AMST                  ! strange mass

      
      A(9)  = FL2Z1*VDQ		! tau mass
      A(12) = FT2Z1*VUQ		! top mass
      A(13) = FB2Z1*VDQ		! bottom mass

c      A(9)  = MLQ             ! tau mass
c      A(12) = MTQ             ! top mass
c      A(13) = MBQ            ! bottom mas
     


      if(nprint.eq.3) then
      
 	if(nprint.ge.3) print*,'MB=',FB2Z1*VDQ
	if(nprint.ge.3) print*,'MT=',FT2Z1*VUQ
	if(nprint.ge.3) print*,'ML=',FL2Z1*VDQ
      ENDIF
        
c      print *,'MTAZMZ1,MBMZ1,MTOMZ1',
c     $      FTAMZ1*VDMZ1, FBOMZ1*VDMZ1, FTOMZ1*VUMZ1
c      print *,'FTAZMZ1,FBOZ1,FTOMZ1',FTAMZ1,FTOMZ1,FBOMZ1
c      print *,'VDMZ1,VUMZ1',VDMZ1,VUMZ1
      
c      print *,'MLQ,MBQ,MTQ',MLQ,MBQ,MTQ
       
     
ccc      FT=MTQ/VUQ
ccc      FB=MBQ/VDQ
ccc      FL=MLQ/VDQ


C---------------------------------------------------C
       A(14) = GAMTOT(IDTP)          ! top width
       A(15) = GAMZ                  ! Z boson width

       A(16) = GAMW                   ! W boson width
       A(17) = SIN(ALFAH)             ! sin(alpha) Higgs mixing angle
       A(18) = MU                     ! mu parameter

       A(19) = -SIN(GAMMAR)           ! chargino mixing angle
       A(20) =  COS(GAMMAR)           ! chargino mixing angle
       A(21) = -THETAY*COS(GAMMAR)    ! chargino mixing angle
       A(22) = -THETAY*SIN(GAMMAR)    ! chargino mixing angle

       A(23) = -SIN(GAMMAL)           ! chargino mixing angle
       A(24) = -THETAX*COS(GAMMAL)    ! chargino mixing angle
       A(25) = -COS(GAMMAL)           ! chargino mixing angle
       A(26) =  THETAX*SIN(GAMMAL)    ! chargino mixing angle

       A(27) =   ZMIXSS(4,1)          ! neutralino mixing angle 1,1
       A(28) =   ZMIXSS(4,2)          ! neutralino mixing angle 1,2
       A(29) =   ZMIXSS(4,3)	     ! neutralino mixing angle 1,3
       A(30) =   ZMIXSS(4,4)	     ! neutralino mixing angle 1,4

       A(31) =   ZMIXSS(3,1)          ! neutralino mixing angle
       A(32) =   ZMIXSS(3,2)          ! neutralino mixing angle
       A(33) =   ZMIXSS(3,3)	     ! neutralino mixing angle
       A(34) =   ZMIXSS(3,4)	     ! neutralino mixing angle

       A(35) =  -ZMIXSS(2,1)	    ! neutralino mixing angle
       A(36) =  -ZMIXSS(2,2)	    ! neutralino mixing angle
       A(37) =  -ZMIXSS(2,3)	  ! neutralino mixing angle
       A(38) =  -ZMIXSS(2,4)	  ! neutralino mixing angle

       A(39) =  -ZMIXSS(1,1)	  ! neutralino mixing angle
       A(40) =  -ZMIXSS(1,2)	  ! neutralino mixing angle
       A(41) =  -ZMIXSS(1,3)	  ! neutralino mixing angle
       A(42) =  -ZMIXSS(1,4)	  ! neutralino mixing angle

       A(43) = GSS(10)               ! Tau soft coupling

       A(44) =  COS(THETAL)          ! stau mixing angle Zl33
       A(45) =  -SIN(THETAL)          ! stau mixing angle Zl36
       A(46) =   SIN(THETAL)          ! stau mixing angle Zl63
       A(47) =  COS(THETAL)          ! stau mixing angle Zl66

       A(56) = GSS(12)               ! top tril soft coupling
       A(57) = GSS(11)               ! bottom tril  soft coupling

       A(58) =  COS(THETAT)           ! squark mixing angle
       A(59) = -SIN(THETAT)          ! squark mixing angle
       A(60) =  SIN(THETAT)          ! squark mixing angle
       A(61) =  COS(THETAT)           ! squark mixing angle

       A(70) =  COS(THETAB)           ! squark mixing angle
       A(71) =  -SIN(THETAB)           ! squark mixing angle
       A(72) =   SIN(THETAB)          ! squark mixing angle
       A(73) =  COS(THETAB)           ! squark mixing angle

       A(82) = MSS(29)               ! Mh
       A(83) = GAMTOT(ISHL)          ! width of light Higgs
       A(84) = MSS(30)               ! MHh
       A(85) = GAMTOT(ISHH)          ! width of heavy higgs
       A(86) = GAMTOT(ISHA)          ! width of CP-odd Higgs
       A(87) = MSS(32)               ! MHC
       A(88) = GAMTOT(ISHC)          ! Width of charged Higgs

       A(89) = ABS(MSS(27))          ! MC1
       A(90) = GAMTOT(ISW1)          ! width of chargino 1
       A(91) = ABS(MSS(28))          ! mass of chargino 2
       A(92) = GAMTOT(ISW2)          ! width of chargino 2
       A(93) = -MSS(23)              ! MNE1

       A(94) = GAMTOT(ISZ1)          ! width of neutralino 1
       A(94) = abs(A(93))/100.

       A(95) = -MSS(24)              ! MNe2
       A(96) = GAMTOT(ISZ2)          ! width of neutralino 2
       A(97) = -MSS(25)              ! mass of neutralino 3
       A(98) = GAMTOT(ISZ3)          ! width of neutralino 3
       A(99) = -MSS(26)              ! mass of neutralino 4
       A(100) = GAMTOT(ISZ4)         ! width of neutralino 4
       A(101) = GAMTOT(ISGL)         ! gluino width

       A(102) = GAMTOT(ISER)         ! 1st selectron width (ISAJET: e_r)
       A(103) = GAMTOT(ISEL)         ! 2nd selectron width (ISAJET: e_l)

       A(104) = GAMTOT(ISMUR)        ! 1st smuon width (ISAJET: mu_r)
       A(105) = GAMTOT(ISMUL)        ! 2nd smuon width (ISAJET: mu_l)

       A(106) = MSS(21)              ! MSTAU1
       A(107) = GAMTOT(ISTAU1)       ! 1st stau width
       A(108) = MSS(22)              ! MSTAU2
       A(109) = GAMTOT(ISTAU2)       ! 2nd stau width
       A(110) = GAMTOT(ISNEL)        ! e-sneutrino width
       A(111) = GAMTOT(ISNML)        ! mu-sneutrino width
       A(112) = MSS(16)              ! tau-sneutrino mass
       A(113) = GAMTOT(ISNTL)        ! tau-sneutrino width

       A(114) = GAMTOT(ISUPL)        ! width of u-squark 1 (ISAJET: u_l)
       A(115) = GAMTOT(ISUPR)        ! width of u-squark 2 (ISAJET: u_r)

       A(116) = GAMTOT(ISDNL)        ! width of d-squark 1 (ISAJET: d_l)
       A(117) = GAMTOT(ISDNR)        ! width of d-squark 2 (ISAJET: d_r)
       A(118) = GAMTOT(ISCHL)        ! width of c-squark 1 (ISAJET: c_l)
       A(119) = GAMTOT(ISCHR)        ! width of c-squark 1 (ISAJET: c_r)
       A(120) = GAMTOT(ISSTL)        ! width of s-squark 1 (ISAJET: s_l)
       A(121) = GAMTOT(ISSTR)        ! width of s-squark 1 (ISAJET: s_r)

       A(122) = MSS(12)              ! MST1
       A(123) = GAMTOT(ISTP1)        ! width of t-squark 1
       A(124) = MSS(13)              ! mass of t-squark 2
       A(125) = GAMTOT(ISTP2)        ! width of t-squark 2
       A(126) = MSS(10)              ! MSB1
       A(127) = GAMTOT(ISBT1)        ! width of b-squark 1
       A(128) = MSS(11)              ! mass of b-squark 2
       A(129) = GAMTOT(ISBT2)        ! width of b-squark 2
       A(130) = XTANB                ! tan beta
       A(131) = MSS(1)               ! gluino mass

       IF(MSS(17).lt.MSS(18)) THEN
         A(132) = MSS(17)              ! MSE1
         A(133) = MSS(18)              ! MSE2
         A(48) =  1
         A(49) =  0
         A(50) =  0
         A(51) =  1
       ELSE
       A(103) = GAMTOT(ISER)         ! 1st selectron width (ISAJET: e_r)
       A(104) = GAMTOT(ISEL)         ! 2nd selectron width (ISAJET: e_l)
         A(132) = MSS(18)
         A(133) = MSS(17)
         A(48) =  0                     ! slepton mixing angle
         A(49) =  1                     ! slepton mixing angle
         A(50) =  -1                     ! slepton mixing angle
         A(51) =  0                     ! slepton mixing angle
       ENDIF

       IF(MSS(19).lt.MSS(20)) THEN
         A(134) = MSS(19)              ! MSMU1
         A(135) = MSS(20)              ! MSMU2
         A(52) =  1                     ! slepton mixing angle
         A(53) =  0                     ! slepton mixing angle
         A(54) =  0                     ! slepton mixing angle
         A(55) =  1                    ! slepton mixing angle
        ELSE
       A(105) = GAMTOT(ISMUR)        ! 1st smuon width (ISAJET: mu_r)
       A(104) = GAMTOT(ISMUL)        ! 2nd smuon width (ISAJET: mu_l)
         A(134) = MSS(20)
         A(135) = MSS(19)
         A(52) =  0                     ! slepton mixing angle
         A(53) =  1                     ! slepton mixing angle
         A(54) =   -1                     ! slepton mixing angle
         A(55) =  0                    ! slepton mixing angle
       ENDIF

       A(136) = MSS(14)              ! e-sneutrino mass
       A(137) = MSS(15)              ! mu-sneutrino mass

       A(62) =  1                     ! squark mixing angle
       A(63) =  0                     ! squark mixing angle
       A(64) =  0                     ! squark mixing angle
       A(65) =  1                     ! squark mixing angle
       A(138) = MSS(2)               ! mass of u-squark 1 (ISAJET: u_l)
       A(139) = MSS(3)               ! mass of u-squark 2 (ISAJET: u_r)
       IF(MSS(2).gt.MSS(3)) then
        A(62) =  0                     ! squark mixing angle
        A(63) =  1                     ! squark mixing angle
        A(64) =  -1                     ! squark mixing angle
        A(65) =  0                     ! squark mixing angle
        A(138) = MSS(3)               ! mass of u-squark 1 (ISAJET: u_l)
        A(139) = MSS(2)               ! mass of u-squark 2 (ISAJET: u_r)
       A(115) = GAMTOT(ISUPL)        ! width of u-squark 1 (ISAJET: u_l)
       A(114) = GAMTOT(ISUPR)        ! width of u-squark 2 (ISAJET: u_r)
       ENDIF

       A(66) =  1                     ! squark mixing angle
       A(67) =  0                     ! squark mixing angle
       A(68) =  0                     ! squark mixing angle
       A(69) =  1                     ! squark mixing angle
       A(142) = MSS(8)               ! mass of c-squark 1 (ISAJET: c_l)
       A(143) = MSS(9)               ! mass of c-squark 1 (ISAJET: c_r)
       IF(MSS(8).gt.MSS(9)) then
       A(66) =  0                     ! squark mixing angle
       A(67) =  1                     ! squark mixing angle
       A(68) =  -1                     ! squark mixing angle
       A(69) =  0                     ! squark mixing angle
       A(142) = MSS(9)               ! mass of c-squark 1 (ISAJET: c_l)
       A(143) = MSS(8)               ! mass of c-squark 1 (ISAJET: c_r)
       A(119) = GAMTOT(ISCHL)        ! width of c-squark 1 (ISAJET: c_l)
       A(118) = GAMTOT(ISCHR)        ! width of c-squark 1 (ISAJET: c_r)
       ENDIF

       A(74) = 1                     ! squark mixing angle
       A(75) = 0                     ! squark mixing angle
       A(76) = 0                     ! squark mixing angle
       A(77) = 1                     ! squark mixing angle
       A(140) = MSS(4)               ! mass of d-squark 1 (ISAJET: d_l)
       A(141) = MSS(5)               ! mass of d-squark 2 (ISAJET: d_r)
       IF(MSS(4).gt.MSS(5)) THEN
       A(74) = 0                     ! squark mixing angle
       A(75) = 1                     ! squark mixing angle
       A(76) = -1                     ! squark mixing angle
       A(77) = 0                     ! squark mixing angle
       A(140) = MSS(5)               ! mass of d-squark 1 (ISAJET: d_l)
       A(141) = MSS(4)               ! mass of d-squark 2 (ISAJET: d_r)
       A(117) = GAMTOT(ISDNL)        ! width of d-squark 1 (ISAJET: d_l)
       A(116) = GAMTOT(ISDNR)        ! width of d-squark 2 (ISAJET: d_r)
       ENDIF

       A(78) = 1                     ! squark mixing angle
       A(79) = 0                     ! squark mixing angle
       A(80) = 0                     ! squark mixing angle
       A(81) = 1                     ! squark mixing angle
       A(144) = MSS(6)               ! mass of s-squark 1 (ISAJET: s_l)
       A(145) = MSS(7)               ! mass of s-squark 1 (ISAJET: s_r)
       IF(MSS(6).gt.MSS(7)) then
       A(78) = 0                     ! squark mixing angle
       A(79) = 1                     ! squark mixing angle
       A(80) = -1                     ! squark mixing angle
       A(81) = 0                     ! squark mixing angle
       A(144) = MSS(7)               ! mass of s-squark 1 (ISAJET: s_l)
       A(145) = MSS(6)               ! mass of s-squark 1 (ISAJET: s_r)
       A(121) = GAMTOT(ISSTL)        ! width of s-squark 1 (ISAJET: s_l)
       A(120) = GAMTOT(ISSTR)        ! width of s-squark 1 (ISAJET: s_r)
       ENDIF
       A(146) = GSS(3)               ! g_s(M_Z) in CompHEP, g_s(Q) in ISAJET
      IF(NPRINT.ge.2) then
      
        Do ii=1,144,3
	  Print '(3(I6,E16.3))', 
     &     ii,A(ii),ii+1,A(ii+1),ii+2,A(ii+2)
        End do
        Do ii=145,146,2   
	  Print '(2(I6,E16.3))', 
     &     ii,A(ii),ii+1,A(ii+1)
        End do
	
!      Print*, GAMTOT(IDW), GAMTOT(IDZ)
!      Stop
      ENDIF

      RETURN
      END
