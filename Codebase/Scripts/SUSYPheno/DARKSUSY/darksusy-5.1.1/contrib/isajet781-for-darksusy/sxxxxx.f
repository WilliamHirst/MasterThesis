C
C
C ----------------------------------------------------------------------
C
      SUBROUTINE SXXXXX(P,NSS , SC)
C
C This subroutine computes a complex SCALAR wavefunction.               
C                                                                       
C INPUT:                                                                
C       real    P(0:3)         : four-momentum of scalar boson          
C       integer NSS  = -1 or 1 : +1 for final, -1 for initial           
C                                                                       
C OUTPUT:                                                               
C       complex SC(3)          : scalar wavefunction                   S
C
      IMPLICIT NONE
      COMPLEX*16 SC(3)
      REAL*8    P(0:3)
      INTEGER NSS
C
      SC(1) = DCMPLX( 1.0 )
      SC(2) = DCMPLX(P(0),P(3))*NSS
      SC(3) = DCMPLX(P(1),P(2))*NSS
C
      RETURN
      END
