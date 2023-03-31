C
C
C ----------------------------------------------------------------------
C
      SUBROUTINE JTIOXX(FI,FO,G , JIO)
C
C this subroutine computes an off-shell vector current from an external 
C fermion pair.  the vector boson propagator is not included in this
C routine.
C                                                                       
C input:                                                                
C       complex fi(6)          : flow-in  fermion                   |fi>
C       complex fo(6)          : flow-out fermion                   <fo|
C       real    g(2)           : coupling constants                  gvf
C                                                                       
C output:                                                               
C       complex jio(6)         : vector current          j^mu(<fo|v|fi>)
C
      IMPLICIT NONE
      COMPLEX*16 FI(6),FO(6),JIO(6)
      REAL*8    G(2)
C
      REAL*8 RXZERO, RXONE
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0 )
      COMPLEX*16 CXIMAG
      LOGICAL FIRST
      SAVE CXIMAG,FIRST
      DATA FIRST/.TRUE./
C
C          Fix compilation with g77
      IF(FIRST) THEN
        FIRST=.FALSE.
        CXIMAG=DCMPLX( RXZERO, RXONE )
      ENDIF
C
      JIO(5) = FO(5)-FI(5)
      JIO(6) = FO(6)-FI(6)
C
      IF ( G(2) .NE. RXZERO ) THEN
         JIO(1) = ( G(1)*( FO(3)*FI(1)+FO(4)*FI(2))
     &             +G(2)*( FO(1)*FI(3)+FO(2)*FI(4)) )
         JIO(2) = (-G(1)*( FO(3)*FI(2)+FO(4)*FI(1))
     &             +G(2)*( FO(1)*FI(4)+FO(2)*FI(3)) )
         JIO(3) = ( G(1)*( FO(3)*FI(2)-FO(4)*FI(1))
     &             +G(2)*(-FO(1)*FI(4)+FO(2)*FI(3)) )*CXIMAG
         JIO(4) = ( G(1)*(-FO(3)*FI(1)+FO(4)*FI(2))
     &             +G(2)*( FO(1)*FI(3)-FO(2)*FI(4)) )
C
      ELSE
         JIO(1) =  ( FO(3)*FI(1)+FO(4)*FI(2))*G(1)
         JIO(2) = -( FO(3)*FI(2)+FO(4)*FI(1))*G(1)
         JIO(3) =  ( FO(3)*FI(2)-FO(4)*FI(1))*DCMPLX(RXZERO,G(1))
         JIO(4) =  (-FO(3)*FI(1)+FO(4)*FI(2))*G(1)
      END IF
C
      RETURN
      END
