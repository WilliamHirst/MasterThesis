CDECK  ID>, FORTOP.
      SUBROUTINE FORTOP
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-     add to force list forced decays for all heavy q particles
C-     if there was a request to force a heavy q decay
C-     Zero IFORCE after use
C-
C-   Created  15-DEC-1989   Serban D. Protopopescu
C-
C    Ver 7.30: Decay top quark rather than hadron, so no longer needed.
C----------------------------------------------------------------------
      IMPLICIT NONE
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
      INTEGER   MXFORC
      PARAMETER (MXFORC=40)
      COMMON/FORCE/NFORCE,IFORCE(MXFORC),MFORCE(5,MXFORC)
     $,LOOK2(2,MXFORC),LOOKST(MXFORC),MEFORC(MXFORC)
      SAVE /FORCE/
      INTEGER   NFORCE,IFORCE,MFORCE,LOOK2,LOOKST,MEFORC
C----------------------------------------------------------------------
      RETURN
      END
