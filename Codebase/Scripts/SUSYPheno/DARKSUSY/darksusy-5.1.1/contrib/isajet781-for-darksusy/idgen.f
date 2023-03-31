CDECK  ID>, IDGEN.
      SUBROUTINE IDGEN
C
C          Call system date and time routines (non-standard) to set up
C          run identification:
C          IDVER=100*VERSN     (integer ISAJET version number)
C          IDG(1)=YYMMDD       (integer year-month-day)
C          IDG(2)=HHMMSS       (integer hour-minute-second)
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/IDRUN/IDVER,IDG(2),IEVT,IEVGEN
      SAVE /IDRUN/
      INTEGER   IDVER,IDG,IEVT,IEVGEN
C          Default run id is zero.
      IYMD=0.
      IHMS=0.
      IDG(1)=IYMD
      IDG(2)=IHMS
      RETURN
      END
