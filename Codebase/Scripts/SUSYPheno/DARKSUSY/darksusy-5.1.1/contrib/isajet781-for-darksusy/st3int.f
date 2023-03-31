C
C       End of functions added by Tadas
C
CDECK  ID>, ST3INT.
      SUBROUTINE ST3INT(INTMSQ)
C
C     Routine by Andrew Box 6/11/2007
CPurpose: To perform the integration and return a value for the integrated
C         matrix element. This could all be moved to the bottom of
C         STBWNEU.F instead of being a separate subroutine.
C
C         Requires DADMUL which is in CERNLIB
C
C         8/6/6
C
      IMPLICIT NONE
C
      DOUBLE PRECISION INTMSQ,A(2),B(2),WK(50000),RELERR
      INTEGER NFNEVL,IFAIL
      EXTERNAL ST3MAT
C
      A(1)=0.D0
      A(2)=0.D0
      B(1)=1.D0
      B(2)=1.D0
C
CInformation about these variables can currently (8/JUN/6) be found at
C
Chttp://wwwasdoc.web.cern.ch/wwwasdoc/shortwrupsdir/d120/top.html
C
      CALL DADMUL(ST3MAT,2,A,B,20,200000,1.D-8,WK,50000,INTMSQ,RELERR,
     .            NFNEVL,IFAIL)
C
      IF(IFAIL.GT.0)THEN
        WRITE(*,*)
        WRITE(*,*)'ERROR IN INTEGRATION OF STBWNEU'
        WRITE(*,*)'RESULT= ',INTMSQ,' TO RELATIVE ACCURACY ',RELERR
        WRITE(*,*)'NUMBER OF EVALUATIONS ',NFNEVL
        IF(IFAIL.EQ.1)THEN
          WRITE(*,*)'NUMBER OF POINTS TOO SMALL FOR SPECIFIED ACCURACY'
        ELSE IF(IFAIL.EQ.2)THEN
          WRITE(*,*)'WORKING ARRAY TOO SMALL'
        ELSE
          WRITE(*,*)'ERROR EITHER IN NUMBER OF DIMENSIONS OR MAXIMUM
     .               OR MINIMUM NUMBER OF POINTS'
        END IF
        WRITE(*,*)
      END IF
C
      RETURN
      END
