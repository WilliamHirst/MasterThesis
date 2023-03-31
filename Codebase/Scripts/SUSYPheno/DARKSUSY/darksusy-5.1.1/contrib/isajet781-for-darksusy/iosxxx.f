C
C ======================================================================
C
      SUBROUTINE IOSXXX(FI,FO,SC,GC , VERTEX)
C
C This subroutine computes an amplitude of the fermion-fermion-scalar   
C coupling.                                                             
C                                                                       
C INPUT:                                                                
C       complex FI(6)          : flow-in  fermion                   |FI>
C       complex FO(6)          : flow-out fermion                   <FO|
C       complex SC(3)          : input    scalar                      S 
C       complex GC(2)          : coupling constants                 GCHF
C                                                                       
C OUTPUT:                                                               
C       complex VERTEX         : amplitude                     <FO|S|FI>
C
      IMPLICIT NONE
      COMPLEX*16 FI(6),FO(6),SC(3),GC(2),VERTEX
C
      VERTEX = SC(1)*( GC(1)*(FI(1)*FO(1)+FI(2)*FO(2))
     &                +GC(2)*(FI(3)*FO(3)+FI(4)*FO(4)) )
C
      RETURN          
      END
