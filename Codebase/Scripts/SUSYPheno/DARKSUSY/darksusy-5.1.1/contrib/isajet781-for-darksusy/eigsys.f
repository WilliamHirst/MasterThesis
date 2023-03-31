C----------------------------------------------------------------------
      SUBROUTINE EIGSYS(NM,N,AR,WR,ZR,IERR,WORK)
C----------------------------------------------------------------------
C     Computes all eigenvalues and eigenvectors of real symmetric
C     NM x NM matrix AR.
C     This is a DOUBLE PRECISION version of the EISRS1 subroutine (F224)
C     from CERN program library.
C
C     Note: Full input matrix AR is preserved.
C     Note: Convention  ZR^T AR ZR = diagonal  is used.
C     Note: Eigenvalues are arranged in non-decreasing order.
C
C     Created: 02/23/07 by Azar Mustafayev
C
      IMPLICIT NONE
      INTEGER NM,N,IERR
      REAL*8 AR(N,N),WR(N),ZR(N,N),WORK(N)
C     WR - vector containing eigenvalues
C     ZR - matrix containing eigenvectors as columns
C     IERR - error parameter, if non-zero the computation has failed
C
      CALL TRDIAG(NM,N,AR,WR,WORK,ZR)
      CALL TQLEIG(NM,N,WR,WORK,ZR,IERR)
c      CALL IMTQLEIG(NM,N,WR,WORK,ZR,IERR) 
      
      RETURN
      END
