C-----------------------------------------------------------------------
      SUBROUTINE MPROD5X(X,I,A,B,C,D,F)
C-----------------------------------------------------------------------
C
C     Multiplies five 3x3 matrices A, B, C, D and F and stores result in 
C     i-th submatrix of 250x3x3 matrix X.
C
C     Note: Uses unfolded DO-loop and three-stage multiplication
c           that greatly increases speed.
C           54*2=108 floating point operations
C
C     Created: 12/13/07 by Azar Mustafayev
C      
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      REAL*8 X(250,3,3),A(3,3),B(3,3),C(3,3),D(3,3),F(3,3),Y(3,3),Z(3,3)
      INTEGER I

c...compute axiliary product matrix Y      
      CALL MPROD3(Y,A,B,C)
      
c...multiply axiliary matrix Y by D      
      CALL MPROD2(Z,Y,D)

c...multiply axiliary matrix Z by remaining matrix F      
      CALL MPROD2X(X,I,Z,F)
      
      RETURN
      END      
