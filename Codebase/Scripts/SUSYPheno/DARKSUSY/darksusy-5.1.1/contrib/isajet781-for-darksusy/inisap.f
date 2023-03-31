CDECK  ID>, INISAP.
      SUBROUTINE INISAP(CMSE,XREAC,BEAMS,WZ,NDCAYS,DCAYS,
     $  ETMIN,RCONE,OK)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-     initialize ISAJET for externally supplied partons
C-   Inputs  : 
C-   CMSE     = center of mass energy
C-   XREAC    = reaction
C-   BEAMS(2) = chose 'P ' or 'AP'
C-   ETMIN    = minimum ET of supplied partons
C-   RCONE    = minimum cone (R) between supplied partons
C-   WZ = option 'W' or 'Z', ' ' no W's or Z's
C-   NDCAYS= number of decay options
C-   DCAYS= list of particles W or Z are allowed to decay into
C-
C-   Controls:
C-   OK   = true if initialization is possible
C-   Created   8-OCT-1991   Serban D. Protopopescu
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER MXKEYS
      PARAMETER (MXKEYS=20)
      COMMON/KEYS/IKEYS,KEYON,KEYS(MXKEYS)
      COMMON/XKEYS/REAC
      SAVE /KEYS/,/XKEYS/
      LOGICAL KEYS
      LOGICAL KEYON
      CHARACTER*8 REAC
      INTEGER   IKEYS
      COMMON/IDRUN/IDVER,IDG(2),IEVT,IEVGEN
      SAVE /IDRUN/
      INTEGER   IDVER,IDG,IEVT,IEVGEN
      COMMON /LIMEVL/ ETTHRS,CONCUT,USELIM
      SAVE /LIMEVL/
      REAL ETTHRS,CONCUT
      LOGICAL USELIM
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
      INTEGER MXGOQ,MXGOJ
      PARAMETER (MXGOQ=85,MXGOJ=8)
      COMMON/Q1Q2/GOQ(MXGOQ,MXGOJ),GOALL(MXGOJ),GODY(4),STDDY,
     $GOWW(25,2),ALLWW(2),GOWMOD(25,MXGOJ)
      SAVE /Q1Q2/
      LOGICAL GOQ,GOALL,GODY,STDDY,GOWW,ALLWW,GOWMOD
      INTEGER MXTYPE
      PARAMETER (MXTYPE=8)
      COMMON/TYPES/LOC(100),NTYP,NJTTYP(MXTYPE),NWWTYP(2),NWMODE(3)
      COMMON/XTYPES/PARTYP(40),TITLE(10),JETYP(30,MXTYPE),WWTYP(30,2)
     $,WMODES(30,3)
      SAVE /TYPES/,/XTYPES/
      CHARACTER*8 JETYP,WWTYP,TITLE,PARTYP,WMODES
      INTEGER   LOC,NTYP,NJTTYP,NWWTYP,NWMODE
C
      REAL CMSE
      CHARACTER*8 XREAC
      CHARACTER*2 BEAMS(2)
      REAL    ETMIN,RCONE
      CHARACTER*1 WZ
      INTEGER NDCAYS
      CHARACTER*4 DCAYS(*)
      LOGICAL OK
      LOGICAL DUMY,SETTYP
      INTEGER I
C----------------------------------------------------------------------
      OK=.TRUE.
      CALL RESET
      IEVT=0
      ECM=CMSE
      SCM=ECM**2
      HALFE=ECM/2.
      ETTHRS=ETMIN
C          fudge factor 1.5 to approximate ET distributions and widths
      CONCUT=SIN(RCONE)/1.5
      IF(RCONE.GT.1.5) CONCUT=1.0
      USELIM=.TRUE.
      IKEYS=0
      DO 18 I=1,8
18    KEYS(I)=.FALSE.
      KEYON=.FALSE.
      REAC=XREAC
C
      IF(XREAC.EQ.'TWOJET  ') THEN
        KEYS(1)=.TRUE.
        IKEYS=1
C
      ELSEIF(XREAC.EQ.'DRELLYAN') THEN
        KEYS(3)=.TRUE.
        IKEYS=3
        IF(WZ.EQ.'Z') GODY(4)=.TRUE.
        IF(WZ.EQ.'W') THEN
          GODY(2)=.TRUE.
          GODY(3)=.TRUE.
        ENDIF
        NJTTYP(1)=NDCAYS
        NJTTYP(2)=0
        NJTTYP(3)=0
        DO 21 I=1,NDCAYS
          JETYP(I,1)=DCAYS(I)
   21   CONTINUE
C
      ELSEIF(XREAC.EQ.'MINBIAS ') THEN
        KEYS(4)=.TRUE.
        IKEYS=4
C
      ELSEIF(XREAC.EQ.'SUPERSYM'.OR.XREAC.EQ.'SUSY    ') THEN
        KEYS(5)=.TRUE.
        IKEYS=5
C
      ELSEIF(XREAC.EQ.'WPAIR   ') THEN
        KEYS(6)=.TRUE.
        IKEYS=6
C
      ELSEIF(XREAC.EQ.'HIGGS   ') THEN
        KEYS(7)=.TRUE.
        IKEYS=7
C
      ELSEIF(XREAC.EQ.'PHOTON  ') THEN
        KEYS(8)=.TRUE.
        IKEYS=8
      ENDIF
C
      IF(IKEYS.EQ.0) THEN
        OK=.FALSE.
        GOTO 999
      ENDIF
C
      CALL SETCON
      IDIN(1)=1120    
      IDIN(2)=-1120   
      IF (BEAMS(1).EQ.'P ') IDIN(1)=1120
      IF (BEAMS(2).EQ.'P ') IDIN(2)=1120
      IF (BEAMS(1).EQ.'AP') IDIN(1)=-1120
      IF (BEAMS(2).EQ.'AP') IDIN(2)=-1120
      DUMY=SETTYP(0)
      CALL SETW
      CALL IDGEN
      CALL SETDKY(.FALSE.)
      CALL MBSET
      CALL PRTLIM
      CALL TIMER(1)
  999 RETURN
      END
