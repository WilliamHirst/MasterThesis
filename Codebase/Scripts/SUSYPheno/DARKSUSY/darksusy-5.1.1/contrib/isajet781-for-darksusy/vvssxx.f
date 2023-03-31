C
      SUBROUTINE VVSSXX(V1,V2,S1,S2,G , VERTEX)
C
C This subroutine computes an amplitude of the vector-vector-scalar-    
C scalar coupling.                                                      
C                                                                       
C INPUT:                                                                
C       complex V1(6)          : first  vector                        V1
C       complex V2(6)          : second vector                        V2
C       complex S1(3)          : first  scalar                        S1
C       complex S2(3)          : second scalar                        S2
C       real    G              : coupling constant                 GVVHH
C                                                                       
C OUTPUT:                                                               
C       complex VERTEX         : amplitude            Gamma(V1,V2,S1,S2)
C
      IMPLICIT NONE
      COMPLEX*16 V1(6),V2(6),S1(3),S2(3),VERTEX
      REAL*8    G
C
      VERTEX = G*S1(1)*S2(1)
     &        *(V1(1)*V2(1)-V1(2)*V2(2)-V1(3)*V2(3)-V1(4)*V2(4))
C
      RETURN
      END
