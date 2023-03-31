C
C   N.B. THERE ARE THREE VERSIONS OF THIS ROUTINE
C   ONE FOR HDECAY 2.0, ONE FOR HDECAY3.0 AND ONE DUMMY
C
CDECK  ID>, HDCYSY
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :  Peter Richardson 
C-----------------------------------------------------------------------
      SUBROUTINE HDCYSY(M2,MEL1,MER1,MQL1,MUR1,MDR1,MEL2,
     &                     MER2,MQL2,MUR2,MDR2,MEL,MER,MSQ,MUR,MDR) 
C-----------------------------------------------------------------------
C  Dummy version of subroutine to interface HDECAY with ISAWIG for MSSM
C  Higgs Decays
C-----------------------------------------------------------------------
      IMPLICIT NONE
C--variables passed from ISAWIG
      REAL M2,MEL1,MER1,MQL1,MUR1,MDR1,MEL2,MER2,MQL2,MUR2,MDR2,
     &     MEL,MER,MSQ,MUR,MDR
      WRITE (6,10)
   10 FORMAT(/10X,'HDECAY CALLED BUT NOT LINKED')
      STOP
      END
