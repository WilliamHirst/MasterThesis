C ======================================================================
C
      SUBROUTINE SSSSXX(S1,S2,S3,S4,G , VERTEX)
C
C This subroutine computes an amplitude of the four-scalar coupling.    
C                                                                       
C INPUT:                                                                
C       complex S1(3)          : first  scalar                        S1
C       complex S2(3)          : second scalar                        S2
C       complex S3(3)          : third  scalar                        S3
C       complex S4(3)          : fourth scalar                        S4
C       real    G              : coupling constant                 GHHHH
C                                                                       
C OUTPUT:                                                               
C       complex VERTEX         : amplitude            Gamma(S1,S2,S3,S4)
C
      IMPLICIT NONE
      COMPLEX*16 S1(3),S2(3),S3(3),S4(3),VERTEX
      REAL*8     G
C
      VERTEX = G*S1(1)*S2(1)*S3(1)*S4(1)
C
      RETURN
      END
