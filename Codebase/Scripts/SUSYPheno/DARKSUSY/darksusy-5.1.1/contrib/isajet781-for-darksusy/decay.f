CDECK  ID>, DECAY.
      SUBROUTINE DECAY(IP)
C
C          Decay particle IP from /PARTCL/ using /DKYTAB/ branching
C          ratios and add decay products to /PARTCL/ with IORIG=IP.
C          Forced decay modes are flagged by LOOK<0.
C
C          Auxiliary routines:
C          DECPS1: generate masses for phase space
C          DECPS2: generate 2-body decays and boosts for phase space
C          DECVA:  V-A matrix elements
C          DECTAU: tau decay matrix elements with polarization
C          DECSS3: 3-body SUSY matrix element using /DKYSS3/
C          DECJET: Hadronize partons from decay.
C
C          Matrix element for Dalitz decays and W mass for TP -> W BT
C          are generated explicitly. W width is included.
C
C          Requirements for decay modes:
C          (1) For Dalitz decays, particle 1 must be GM.
C          (2) For V-A quark or lepton decays, particles 1 and 2 must
C              be from (virtual) W.
C          (3) For any decay into quarks, they must appear last.
C
C          Matrix element flags:
C          MELEM=0     phase space
C                1     Dalitz
C                2     omega/phi
C                3     V-A
C                4     top
C                5     tau -> e nu nu
C                6     tau -> pi nu
C                7     tau -> rho nu
C                8     tau -> tau (for NOTAU)
C                9     H -> W f fbar
C
      IMPLICIT NONE
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
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
      INTEGER   MXPTCL,IPACK
      PARAMETER (MXPTCL=4000,IPACK=10000)
      COMMON/PARTCL/NPTCL,PPTCL(5,MXPTCL),IORIG(MXPTCL),IDENT(MXPTCL)
     1,IDCAY(MXPTCL)
      SAVE /PARTCL/
      INTEGER   NPTCL,IORIG,IDENT,IDCAY
      REAL      PPTCL
C          LOOK must be dimensioned to the maximum value of INDEX.
      INTEGER   MXLOOK
      PARAMETER (MXLOOK=500)
      INTEGER   MXDKY
      PARAMETER (MXDKY=3000)
      COMMON/DKYTAB/LOOK(MXLOOK),CBR(MXDKY),MODE(5,MXDKY),MELEM(MXDKY)
      SAVE /DKYTAB/
      INTEGER   LOOK,MODE,MELEM
      REAL      CBR
      INTEGER   MXJSET,JPACK
      PARAMETER (MXJSET=400,JPACK=1000)
      COMMON/JETSET/NJSET,PJSET(5,MXJSET),JORIG(MXJSET),JTYPE(MXJSET),
     $JDCAY(MXJSET)
      SAVE /JETSET/
      INTEGER   NJSET,JORIG,JTYPE,JDCAY
      REAL      PJSET
      COMMON/JWORK/ZZC(MXJSET),JMATCH(MXJSET),TNEW,P1CM(4),
     1J1,J2,J3,J4,J5,E1CM,E2CM,E3CM,E4CM,E5CM
      SAVE /JWORK/
      LOGICAL TNEW
      EQUIVALENCE (J1,JJ(1)),(E1CM,EE(1))
      INTEGER   JMATCH,J1,J2,J3,J4,J5,JJ(5)
      REAL      ZZC,P1CM,E1CM,E2CM,E3CM,E4CM,E5CM,EE(5)
      COMMON/CONST/PI,SQRT2,ALFA,GF,UNITS
      SAVE /CONST/
      REAL      PI,SQRT2,ALFA,GF,UNITS
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
      COMMON/IDRUN/IDVER,IDG(2),IEVT,IEVGEN
      SAVE /IDRUN/
      INTEGER   IDVER,IDG,IEVT,IEVGEN
      INTEGER   MXFORC
      PARAMETER (MXFORC=40)
      COMMON/FORCE/NFORCE,IFORCE(MXFORC),MFORCE(5,MXFORC)
     $,LOOK2(2,MXFORC),LOOKST(MXFORC),MEFORC(MXFORC)
      SAVE /FORCE/
      INTEGER   NFORCE,IFORCE,MFORCE,LOOK2,LOOKST,MEFORC
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
C
C          Data for SUSY 3-body matrix elements. There is a double 
C          pointer structure, first to modes, and then to poles that
C          make up the matrix element for that mode:
C          MELEM=-I in /DKYTAB/ points to the mode information:
C            J1SS3(I) = start of pole list for this mode
C            J2SS3(I) = end of pole list for this mode
C            WTSS3(I) = maximum weight for this mode
C          J1SS3<J<J2SS3 points to the corresponding poles:
C            KSS3(J)    = pole type
C            AMSS3(J)   = pole mass
C            ZISS3(2,J) = initial couplings
C            ZFSS3(2,J) = final couplings
C          For gaugino -> gaugino f fbar, the pole types are
C            KSS3=1: spin-1 pole in f-fbar channel
C            KSS3=2: spin-0 pole in gaugino-f channel
C            KSS3=3: spin-0 pole in gaugino-fbar channel
C            KSS3=4: spin-0 pole in f-fbar channel
C          The two couplings are the coefficients of 1,gamma_5 or of
C          gamma_mu,gamma_mu*gamma_5. 
C
      INTEGER MXMSS3,MXPSS3
      PARAMETER (MXMSS3=1000)
      PARAMETER (MXPSS3=2000)
      COMMON/DKYSS3/NMSS3,NPSS3,
     $J1SS3(MXMSS3),J2SS3(MXMSS3),WTSS3(MXMSS3),
     $KSS3(MXPSS3),AMSS3(MXPSS3),ZISS3(2,MXPSS3),ZFSS3(2,MXPSS3)
      INTEGER NMSS3,NPSS3,KSS3,J1SS3,J2SS3
      REAL WTSS3,AMSS3
      COMPLEX ZISS3,ZFSS3
C
      REAL PGEN(5,5),BETA(3),REDUCE(5),WPROP,Z,TRY,RANF,AMASS,TWOME
      REAL PSUM(5),SUM,PREST(4,6),DOT,PCM
      REAL AMEE,REE,WTEE,SWAP,WT,A,B,C,GAMMA
      REAL SMAX,SMIN,SVAL,TANMAX,TANMIN,TANVAL
      LOGICAL WDECAY,DECVA,DECTAU,DECJET
      INTEGER IDLV1,IFL1,IFL2,IFL3,JSPIN,INDEX,IPOINT,ID1,I1,I2
      INTEGER NADD,NSTART,NEW,NADD1,J,IP,I,IDABS(5)
      INTEGER K,JETIP,IDANTI,NPASS,MEIP,MEA
      REAL DBLPCM,DECSS3,VAL
      REAL ZZSTAR
      INTEGER IW
C
      DATA REDUCE/1.,1.,2.,5.,15./
      DATA PSUM/5*0./
      DATA TWOME/1.022006E-3/
      DATA PREST/24*0./
C
C          Function definitions.
C          Use double precision for PCM on 32-bit machines
C
      PCM(A,B,C)=DBLPCM(A,B,C)
      DOT(I1,I2)=PREST(4,I1)*PREST(4,I2)-PREST(1,I1)*PREST(1,I2)
     $-PREST(2,I1)*PREST(2,I2)-PREST(3,I1)*PREST(3,I2)
C          Charged W propagator.
      WPROP(Z)=(Z-WMASS(2)**2)**2+(WMASS(2)*WGAM(2))**2
C----------------------------------------------------------------------
C          Select decay mode. Note IDENT(NPTCL+1)...IDENT(NPTCL+5)
C          are always defined even if zero.
C----------------------------------------------------------------------
      IF(IDCAY(IP).NE.0) RETURN
      IDLV1=IDENT(IP)
      CALL FLAVOR(IDLV1,IFL1,IFL2,IFL3,JSPIN,INDEX)
C          FLAVOR returns 0 for quark, but want IFL3=6 for top
      IF(IABS(IDLV1).LT.10) IFL3=IDLV1
      NPASS=0
1     CONTINUE
      NPASS=NPASS+1
      WDECAY=.FALSE.
      IF(NPASS.GT.NTRIES) GO TO 9998
      IPOINT=LOOK(INDEX)
      IF(IPOINT.EQ.0) RETURN
C          IPOINT<0 flags a forced decay.
      IF(IPOINT.LT.0) THEN
        I=1
        IF(IDENT(IP).LT.0) I=2
        IPOINT=LOOK2(I,IABS(IPOINT))
      ENDIF
C
C          Select mode.
C
      TRY=RANF()
      IPOINT=IPOINT-1
100   IPOINT=IPOINT+1
      IF(TRY.GT.CBR(IPOINT)) GO TO 100
      NADD=0
      SUM=0.
      NSTART=NPTCL+1
      IF(NPTCL+5.GT.MXPTCL) GO TO 9999
C
C          Set up masses and IDENT codes.
C
      MEIP=MELEM(IPOINT)
      DO 110 I=1,5
        NEW=NPTCL+I
        IDENT(NEW)=MODE(I,IPOINT)
        IDABS(I)=IABS(IDENT(NEW))
        IF(MODE(I,IPOINT).EQ.0) GO TO 110
        NADD=NADD+1
        IDLV1=IDENT(NEW)
        PPTCL(5,NEW)=AMASS(IDLV1)
        SUM=SUM+PPTCL(5,NEW)
110   CONTINUE
      NADD1=NADD-1
      DO 120 J=1,5
        PGEN(J,1)=PPTCL(J,IP)
120   CONTINUE
      PGEN(5,NADD)=PPTCL(5,NPTCL+NADD)
C----------------------------------------------------------------------
C          Carry out appropriate decay
C----------------------------------------------------------------------
C
C          1-body decays.
C          Determine polarization mode for 1-body tau decays
C
      IF(NADD.EQ.1) THEN
        DO 200 J=1,5
          PPTCL(J,NPTCL+1)=PPTCL(J,IP)
200     CONTINUE
        IF(MEIP.EQ.8) THEN
          IF(DECTAU(IP,NADD,MEIP,IDABS,PREST)) THEN
            IDENT(NPTCL+1)=IDTAUL
          ELSE
            IDENT(NPTCL+1)=IDTAUR
          ENDIF
        ENDIF
        GO TO 300
      ENDIF
C
C          2-body phase space decays
C
      IF(NADD.EQ.2.AND.MEIP.EQ.0) THEN
        CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
        GO TO 300
      ENDIF
C
C          N-body phase space decays
C
      IF(NADD.GT.2.AND.MEIP.EQ.0) THEN
        CALL DECPS1(IP,NADD,PGEN)
        CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
        GO TO 300
      ENDIF
C
C          Dalitz decays
C
      IF(NADD.EQ.3.AND.MEIP.EQ.1) THEN
210     AMEE=TWOME*(PPTCL(5,IP)/TWOME)**RANF()
        REE=(TWOME/AMEE)**2
        WTEE=(1.-(AMEE/PPTCL(5,IP))**2)**3*SQRT(1.-REE)*(1.+.5*REE)
        IF(WTEE.LT.RANF()) GO TO 210
        PGEN(5,2)=AMEE
        CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
        GO TO 300
      ENDIF
C
C          omega/phi decays (for reasons lost in history...)
C
      IF(NADD.EQ.3.AND.MEIP.EQ.2) THEN
220     CALL DECPS1(IP,NADD,PGEN)
        CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
        WT=(PPTCL(5,NPTCL+1)*PPTCL(5,NPTCL+2)*PPTCL(5,NPTCL+3))**2
     $  -(PPTCL(5,NPTCL+1)*DOT(2,3))**2
     $  -(PPTCL(5,NPTCL+2)*DOT(1,3))**2
     $  -(PPTCL(5,NPTCL+3)*DOT(1,2))**2
     $  +2.*DOT(1,2)*DOT(2,3)*DOT(1,3)
        IF(WT.LT.RANF()*PPTCL(5,IP)**6/108.) GO TO 220
        GO TO 300
      ENDIF
C
C          V-A decays
C
      IF(NADD.EQ.3.AND.MEIP.EQ.3) THEN
230     CALL DECPS1(IP,NADD,PGEN)
        CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
        IF(.NOT.DECVA(IP,NADD,IDABS,PREST)) GO TO 230
        GO TO 300
      ENDIF
C
C          Top decays
C          Generate mass for TP -> W BT with Breit-Wigner. 
C          W couples to 1+2 so swap 1<->3. Then m2+m3 < m < m0-m1.
C
      IF(NADD.EQ.3.AND.MEIP.EQ.4) THEN
        WDECAY=.TRUE.
        SWAP=PPTCL(5,NPTCL+1)
        PPTCL(5,NPTCL+1)=PPTCL(5,NPTCL+3)
        PPTCL(5,NPTCL+3)=SWAP
        SMAX=(PPTCL(5,IP)-PPTCL(5,NPTCL+1))**2
        SMIN=(PPTCL(5,NPTCL+2)+PPTCL(5,NPTCL+3))**2
        TANMAX=ATAN((SMAX-WMASS(2)**2)/(WMASS(2)*WGAM(2)))
        TANMIN=ATAN((SMIN-WMASS(2)**2)/(WMASS(2)*WGAM(2)))
240     TANVAL=RANF()*(TANMAX-TANMIN)+TANMIN
        SVAL=WMASS(2)**2+WMASS(2)*WGAM(2)*TAN(TANVAL)
        IF(SVAL.LT.SMIN.OR.SVAL.GT.SMAX) GO TO 240
        PGEN(5,2)=SQRT(SVAL)
        PGEN(5,3)=PPTCL(5,NPTCL+3)
        CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
        IF(.NOT.DECVA(IP,NADD,IDABS,PREST)) GO TO 240
        DO 241 K=1,5
          SWAP=PPTCL(K,NPTCL+1)
          PPTCL(K,NPTCL+1)=PPTCL(K,NPTCL+3)
          PPTCL(K,NPTCL+3)=SWAP
241     CONTINUE
        PGEN(5,3)=PPTCL(5,NPTCL+3)
        DO 242 K=1,4
          SWAP=PREST(K,1)
          PREST(K,1)=PREST(K,3)
          PREST(K,3)=SWAP
242     CONTINUE
        GO TO 300
      ENDIF
C
C          TAU decays. These are special because they take polarization
C          into account.
C
      IF(MEIP.EQ.5.OR.MEIP.EQ.6.OR.MEIP.EQ.7) THEN
250     CALL DECPS1(IP,NADD,PGEN)
        CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
        IF(.NOT.DECTAU(IP,NADD,MEIP,IDABS,PREST)) GO TO 250
        GO TO 300
      ENDIF
C
C          3-body SUSY decays
C
      IF(MEIP.LT.0.AND.NADD.EQ.3) THEN
        MEA=IABS(MEIP)
        IF(WTSS3(MEA).LE.0) THEN
          DO 260 I=1,1000
            CALL DECPS1(IP,NADD,PGEN)
            CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
            VAL=DECSS3(IP,MEA)
            WTSS3(MEA)=MAX(WTSS3(MEA),VAL)
260       CONTINUE
          IF(WTSS3(MEA).LE.0) GO TO 9998
        ENDIF
261     CALL DECPS1(IP,NADD,PGEN)
        CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
        VAL=DECSS3(IP,MEA)
        WTSS3(MEA)=MAX(WTSS3(MEA),VAL)
        IF(VAL.LT.WTSS3(MEA)*RANF()) GO TO 261
        GO TO 300
      ENDIF
C
C          H -> W f fbar decays
C          Generate f fbar mass using ZZSTAR function
C
      IF(NADD.EQ.3.AND.MEIP.EQ.9) THEN
        IF(IDENT(NPTCL+1).EQ.80) THEN
          IW=2
        ELSEIF(IDENT(NPTCL+1).EQ.-80) THEN
          IW=3
        ELSEIF(IDENT(NPTCL+1).EQ.90) THEN
          IW=4
        ELSE
          WRITE(ITLIS,*) 'ERROR IN DECAY ... BAD H -> W F FBAR'
          STOP99
        ENDIF
        PGEN(5,2)=ZZSTAR(PPTCL(5,IP),IW)
        IF(PGEN(5,2).LT.PPTCL(5,NPTCL+2)+PPTCL(5,NPTCL+3)+1.0)
     $  GO TO 1
        CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
        GO TO 300
      ENDIF
C
C          Should never fall through
C
      GO TO 9998
C----------------------------------------------------------------------
C          Swap particles and antiparticles if IDENT(IP)<0
C          Note forced modes for antiparticles are conjugated in table.
C----------------------------------------------------------------------
300   CONTINUE
      IF(IDENT(IP).LT.0.AND.IDENT(IP).NE.-20) THEN
        DO 310 I=1,NADD
          ID1=IDENT(NPTCL+I)
          IDENT(NPTCL+I)=IDANTI(ID1)
310     CONTINUE
      ENDIF
C
C          Set IORIG and IDCAY.
C
      NPTCL=NPTCL+NADD
      IDCAY(IP)=IPACK*NSTART+NPTCL
      JETIP=IABS(IORIG(IP))/IPACK
      DO 320 I=NSTART,NPTCL
        IORIG(I)=IP
        IDCAY(I)=0
320   CONTINUE
C
C          Evolve and hadronize partons. If it fails, start over.
C
      IF (.NOT.WRTLHE) THEN
      IF(IDABS(NADD).LT.10.OR.MOD(IDENT(NPTCL),100).EQ.0) THEN
        IF(.NOT.DECJET(IP,NADD,IDABS,PREST,WDECAY,BETA,GAMMA))
     $  GO TO 1
      ENDIF
      END IF
C
      RETURN
C----------------------------------------------------------------------
C          Error messages.
C----------------------------------------------------------------------
9999  CALL PRTEVT(0)
      WRITE(ITLIS,99990) NPTCL
99990 FORMAT(//5X,'ERROR IN DECAY...NPTCL > ',I6)
      RETURN
9998  CALL PRTEVT(0)
      WRITE(ITLIS,99980) IP
99980 FORMAT(//5X,'ERROR IN DECAY...NO DECAY FOUND FOR PARTICLE',I6)
      RETURN
      END
