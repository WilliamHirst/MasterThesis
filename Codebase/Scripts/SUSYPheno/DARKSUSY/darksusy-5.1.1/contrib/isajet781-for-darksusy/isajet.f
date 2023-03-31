CDECK  ID>, ISAJET.
      SUBROUTINE ISAJET(JTDKY,JTEVT,JTCOM,JTLIS)
C
C          Main subroutine for ISAJET, a Monte Carlo event generator
C          for  P P ,  AP P , and  E+ E-  interactions at high energy.
C
C          Frank E. Paige and Serban D. Protopopescu
C          Brookhaven National Laboratory
C          Upton, New York, USA
C
C          JTDKY = +/- unit number for decay table file.
C                      If it is negative, decay table is not printed.
C          JTEVT = +/- unit number for output event file.
C                      If it is negative, only stable particles are
C                      written on it.
C          JTCOM =     unit number for command file.
C          JTLIS =     unit number for listing.
C
C          Instead of calling this subroutine the user may wish to 
C          control the program himself using:
C          ISAINI      overall initialization
C          ISABEG      run initialization
C          ISAEVT      generation of one event
C          ISAEND      run termination
C          ISAWBG      initial record writing
C          ISAWEV      event record writing
C          ISAWND      end record writing
C
      IMPLICIT NONE
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/IDRUN/IDVER,IDG(2),IEVT,IEVGEN
      SAVE /IDRUN/
      INTEGER   IDVER,IDG,IEVT,IEVGEN
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
      COMMON/ISLOOP/NEVOLV,NFRGMN,IEVOL,IFRG
      SAVE /ISLOOP/
      INTEGER NEVOLV,NFRGMN,IEVOL,IFRG
C
      INTEGER JTDKY,JTEVT,JTCOM,JTLIS,IFL,ILOOP
      LOGICAL OK,DONE
      SAVE ILOOP
C
C          Initialize ISAJET
C
      CALL ISAINI(JTDKY,JTEVT,JTCOM,JTLIS)
C
C          Read instructions; terminate for STOP command or error.
C
    1 IFL=0
      CALL ISABEG(IFL)
      IF(IFL.NE.0) RETURN
C          Write begin-run record
      CALL ISAWBG
C
C          Event loop
C
      ILOOP=0
  101 CONTINUE
        ILOOP=ILOOP+1
C          Generate one event - discard if .NOT.OK
        CALL ISAEVT(ILOOP,OK,DONE)
C          Write event record
        IF(OK) CALL ISAWEV
      IF(.NOT.DONE) GO TO 101
C
C          Calculate cross section and luminosity
C
      CALL ISAEND
C          Write end-of-run record
      CALL ISAWND
      GO TO 1
C
C          Entry point for error recovery.
C          CALL RSTART will continue generation on next event.
C
      ENTRY RSTART
      IF(IEVT.EQ.0) RETURN
      IF(IEVT.GE.NEVENT*NEVOLV*NFRGMN) GO TO 1
      GO TO 101
      END
