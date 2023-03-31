CDECK  ID>, ISAEVT.
      SUBROUTINE ISAEVT(I,OK,DONE)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods :
C-
C-         Normal operation:
C-    Generate one ISAJET event and return.
C-
C-         "ISALEP" generation:
C-    Generate a TWOJET or DRELLYAN hard scattering. Then make NEVOLVE 
C-    evolutions and NHADRON fragmentations, rejecting events which 
C-    fail the desired cuts using logical functions
C-         REJJET()   tests the QCD evolution stage, e.g. by requiring
C-                    a heavy quark.
C-         REJFRG()   tests the fragmentation stage, e.g. by requiring
C-                    a high-pt lepton.
C-    These functions default to .FALSE.; i.e. they do not reject any 
C-    events. Note that one hard scattering can give more than one 
C-    event. You must choose NEVOLVE and NHADRON carefully.
C-         IEVT   = event number. This is incremented NEVOLVE * NHADRON 
C-                  times for each hard scattering; i.e. it counts the
C-                  number of potential events.
C-         IEVGEN = counter for generated events.
C-         NEVENT = maximum value of hard scatterings. Hence the limit
C-                  for IEVT is NEVENT * NEVOLVE * NHADRON.
C-    The cross section SIGF contains an extra factor of
C-         1 / (NEVOLVE * NHADRON)
C-    to produce the correct final cross section using the weight
C-         SIGF / NEVENT
C-
C-   Input:
C-    I      = number used to control printout
C-   Output:
C-    OK     = logical flag for good event.
C-    DONE   = logical flag for job completion.
C-
C-   Created   3-FEB-1988   Serban D. Protopopescu
C-   Updated  17-APR-1990 Serban D. Protopopescu (add ISALEP option) 
C-   22-JUL-1992: Move PRTEVT and GETTOT statements to end so they
C-                work for TWOJET and DRELLYAN with NOVOLVE. (FEP)
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/IDRUN/IDVER,IDG(2),IEVT,IEVGEN
      SAVE /IDRUN/
      INTEGER   IDVER,IDG,IEVT,IEVGEN
      INTEGER MXKEYS
      PARAMETER (MXKEYS=20)
      COMMON/KEYS/IKEYS,KEYON,KEYS(MXKEYS)
      COMMON/XKEYS/REAC
      SAVE /KEYS/,/XKEYS/
      LOGICAL KEYS
      LOGICAL KEYON
      CHARACTER*8 REAC
      INTEGER   IKEYS
      COMMON/NODCAY/NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOGRAV,
     $NOB,NOTAU
      LOGICAL NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOGRAV,
     $NOB,NOTAU
      SAVE /NODCAY/
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
      COMMON/JETPAR/P(3),PT(3),YJ(3),PHI(3),XJ(3),TH(3),CTH(3),STH(3)
     1 ,JETTYP(3),SHAT,THAT,UHAT,QSQ,X1,X2,PBEAM(2)
     2 ,QMW,QW,QTW,YW,XW,THW,QTMW,PHIW,SHAT1,THAT1,UHAT1,JWTYP
     3 ,ALFQSQ,CTHW,STHW,Q0W
     4 ,INITYP(2),ISIGS,PBEAMS(5)
      SAVE /JETPAR/
      INTEGER   JETTYP,JWTYP,INITYP,ISIGS
      REAL      P,PT,YJ,PHI,XJ,TH,CTH,STH,SHAT,THAT,UHAT,QSQ,X1,X2,
     +          PBEAM,QMW,QW,QTW,YW,XW,THW,QTMW,PHIW,SHAT1,THAT1,UHAT1,
     +          ALFQSQ,CTHW,STHW,Q0W,PBEAMS
      INTEGER   MXPTCL,IPACK
      PARAMETER (MXPTCL=4000,IPACK=10000)
      COMMON/PARTCL/NPTCL,PPTCL(5,MXPTCL),IORIG(MXPTCL),IDENT(MXPTCL)
     1,IDCAY(MXPTCL)
      SAVE /PARTCL/
      INTEGER   NPTCL,IORIG,IDENT,IDCAY
      REAL      PPTCL
      INTEGER   MXJSET,JPACK
      PARAMETER (MXJSET=400,JPACK=1000)
      COMMON/JETSET/NJSET,PJSET(5,MXJSET),JORIG(MXJSET),JTYPE(MXJSET),
     $JDCAY(MXJSET)
      SAVE /JETSET/
      INTEGER   NJSET,JORIG,JTYPE,JDCAY
      REAL      PJSET
      COMMON/ISLOOP/NEVOLV,NFRGMN,IEVOL,IFRG
      SAVE /ISLOOP/
      INTEGER NEVOLV,NFRGMN,IEVOL,IFRG
C
      LOGICAL REJJET,REJFRG,OK,DONE
      INTEGER NPASS,I,NLIMIT
C
      IF (WRTLHE) THEN
        NOEVOL=.TRUE.
        NOHADR=.TRUE.
      END IF
      NPASS=0
      OK=.TRUE.
      DONE=.FALSE.
      NLIMIT=NEVENT*NEVOLV*NFRGMN
C
C          Twojet or Drell-Yan events. The evolution and fragmentation
C          loops are done with GO TO statements so that we can exit
C          the loops with a good event and reenter them.
C
      IF(KEYS(1).OR.KEYS(3)) THEN
100     CONTINUE
        IF(IEVOL.EQ.1.AND.IFRG.EQ.1) THEN
          NPASS=NPASS+1
          IF(NPASS.GT.NTRIES) THEN
            WRITE(ITLIS,1001) NTRIES
1001        FORMAT(//' IT IS TAKING MORE THAN',I6,' TRIES TO MAKE',
     $      ' AN EVENT IN ISAEVT.'/
     $      ' CHECK YOUR LIMITS OR OR INCREASE NTRIES.'/
     $      ' CHECK NEVOLVE, NHADRON, AND YOUR REJJET AND REJFRG',
     $      ' FUNCTIONS IF ANY.'/
     $      ' JOB TERMINATED.')
            STOP 99
          ENDIF
          CALL RANFMT
C          Generate appropriate hard scattering
          IF(KEYS(1)) THEN
            CALL TWOJET
          ELSE
            CALL DRLLYN
          ENDIF
        ENDIF
C          QCD evolution
        IF(NOEVOL) THEN
          IEVT=IEVT+NEVOLV*NFRGMN
          GOTO 9999
        ENDIF
C          Continue if in fragmentation loop
        IF(IFRG.NE.1) GO TO 120
C          Begin multiple evolution loop
110       CONTINUE
          NJSET=0
          IEVT=IEVT+1
          CALL EVOLVE
          IEVT=IEVT-1
          IF(NJSET.LT.0) THEN
            IEVT=IEVT+NFRGMN
            GO TO 111
          ENDIF
          IF(REJJET()) THEN
            IEVT=IEVT+NFRGMN
            GO TO 111
          ENDIF
          IF(NOHADR) THEN
            IEVT=IEVT+NFRGMN
            GO TO 9999
          ENDIF
C          Begin multiple fragmentation loop
120         CONTINUE
            NPTCL=0
            CALL FRGMNT
            IEVT=IEVT+1
            IF(REJFRG()) GO TO 121
C          Finish good event
            CALL MBIAS
            IFRG=IFRG+1
            IF(IFRG.GT.NFRGMN) IFRG=1
            IF(IFRG.EQ.1) THEN
              IEVOL=IEVOL+1
              IF(IEVOL.GT.NEVOLV) IEVOL=1
            ENDIF
            GOTO 9999
C          Fragmentation failed - increment counter and loop
121         IFRG=IFRG+1
            IF(IFRG.GT.NFRGMN) THEN
              IFRG=1
            ELSE
              GO TO 120
            ENDIF
C          End of multiple fragmentation loop
C          Evolution failed - increment counter and loop
111       IEVOL=IEVOL+1
          IF(IEVOL.GT.NEVOLV) THEN
            IEVOL=1
            IFRG=1
            GO TO 100
          ELSE
            GO TO 110
          ENDIF
C
C          E+E- events
C
      ELSE IF(KEYS(2)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL ELCTRN
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(.NOT.NOHADR) CALL FRGMNT
        ENDIF
C
C          MINBIAS events
C
      ELSE IF(KEYS(4)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        NPTCL=0
        IF(.NOT.(NOEVOL.OR.NOHADR)) CALL MBIAS
C
C          SUPERSYM events
C
      ELSE IF(KEYS(5)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL TWOJET
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(NJSET.LT.0) GO TO 9999
          IF(.NOT.NOHADR) THEN
            CALL FRGMNT
            CALL MBIAS
          ENDIF
        ENDIF
C
C          WPAIR events
C
      ELSE IF(KEYS(6)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL TWOJET
        CALL WPAIR
C
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(NJSET.LT.0) GO TO 9999
          IF(.NOT.NOHADR) THEN
            CALL FRGMNT
            CALL MBIAS
          ENDIF
        ENDIF
C
C          HIGGS events
C
      ELSE IF(KEYS(7)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL DRLLYN
        CALL HIGGS
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(NJSET.LT.0) GOTO 9999
          IF(.NOT.NOHADR) THEN
            CALL FRGMNT
            CALL MBIAS
          ENDIF
        ENDIF
C
C          PHOTON events
C
      ELSEIF(KEYS(8)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL TWOJET
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(NJSET.LT.0) GOTO 9999
          IF(.NOT.NOHADR) THEN
            CALL FRGMNT
            CALL MBIAS
          ENDIF
        ENDIF
C
C          TCOLOR events, e.g. techni-rho
C
      ELSEIF(KEYS(9)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL DRLLYN
        CALL HIGGS
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(NJSET.LT.0) GOTO 9999
          IF(.NOT.NOHADR) THEN
            CALL FRGMNT
            CALL MBIAS
          ENDIF
        ENDIF
C
C          WHIGGS events
C
      ELSE IF(KEYS(10)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL TWOJET
        CALL WHIGGS
C
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(NJSET.LT.0) GO TO 9999
          IF(.NOT.NOHADR) THEN
            CALL FRGMNT
            CALL MBIAS
          ENDIF
        ENDIF
C
C          EXTRADIM events
C
      ELSE IF(KEYS(11)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL DRLLYN
C 
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(NJSET.LT.0) GO TO 9999
          IF(.NOT.NOHADR) THEN
            CALL FRGMNT
            CALL MBIAS
          ENDIF
        ENDIF
C
C          ZJJ events
C
      ELSEIF(KEYS(12)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL ZJJ
C
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(NJSET.LT.0) GO TO 9999
          IF(.NOT.NOHADR) THEN
            CALL FRGMNT
            CALL MBIAS
          ENDIF
        ENDIF
      ENDIF
C
C          Event complete
C
 9999 IEVGEN=IEVGEN+1
      IF(NJSET.LT.0) OK=.FALSE.
      IF(IEVT.GT.NLIMIT) THEN
        OK=.FALSE.
        DONE=.TRUE.
      ELSEIF(IEVT.EQ.NLIMIT) THEN
        DONE=.TRUE.
      ENDIF
      IF (WRTLHE) THEN
        CALL ISALHE
      END IF
      IF(OK) THEN
        CALL PRTEVT(I)
        CALL GETTOT(.FALSE.)
      ENDIF
      RETURN
      END
