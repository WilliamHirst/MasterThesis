C
C ======================================================================
C
      SUBROUTINE IOVXXX(FI,FO,VC,G , VERTEX)
C
C this subroutine computes an amplitude of the fermion-fermion-vector   
C coupling.                                                             
C                                                                       
C input:                                                                
C       complex fi(6)          : flow-in  fermion                   |fi>
C       complex fo(6)          : flow-out fermion                   <fo|
C       complex vc(6)          : input    vector                      v 
C       real    g(2)           : coupling constants                  gvf
C                                                                       
C output:                                                               
C       complex vertex         : amplitude                     <fo|v|fi>
C
      IMPLICIT NONE
      COMPLEX*16 FI(6),FO(6),VC(6),VERTEX
      REAL*8    G(2)
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
      VERTEX =  G(1)*( (FO(3)*FI(1)+FO(4)*FI(2))*VC(1)
     &                +(FO(3)*FI(2)+FO(4)*FI(1))*VC(2)
     &                -(FO(3)*FI(2)-FO(4)*FI(1))*VC(3)*CXIMAG
     &                +(FO(3)*FI(1)-FO(4)*FI(2))*VC(4)        )
C
      IF ( G(2) .NE. RXZERO ) THEN
         VERTEX = VERTEX
     &          + G(2)*( (FO(1)*FI(3)+FO(2)*FI(4))*VC(1)
     &                  -(FO(1)*FI(4)+FO(2)*FI(3))*VC(2)
     &                  +(FO(1)*FI(4)-FO(2)*FI(3))*VC(3)*CXIMAG
     &                  -(FO(1)*FI(3)-FO(2)*FI(4))*VC(4)        )
      END IF
C
      RETURN
      END
