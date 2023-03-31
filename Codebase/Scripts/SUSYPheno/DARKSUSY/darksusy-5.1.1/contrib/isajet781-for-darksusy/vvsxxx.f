C
C
C ======================================================================
C
      SUBROUTINE VVSXXX(V1,V2,SC,G , VERTEX)
C
C this subroutine computes an amplitude of the vector-vector-scalar     
C coupling.                                                             
C                                                                       
C input:                                                                
C       complex v1(6)          : first  vector                        v1
C       complex v2(6)          : second vector                        v2
C       complex sc(3)          : input  scalar                        s 
C       real    g              : coupling constant                  gvvh
C                                                                       
C output:                                                               
C       complex vertex         : amplitude                gamma(v1,v2,s)
C
      IMPLICIT NONE
      COMPLEX*16 V1(6),V2(6),SC(3),VERTEX
      REAL*8    G
C
      VERTEX = G*SC(1)*(V1(1)*V2(1)-V1(2)*V2(2)-V1(3)*V2(3)-V1(4)*V2(4))
C
      RETURN
      END
