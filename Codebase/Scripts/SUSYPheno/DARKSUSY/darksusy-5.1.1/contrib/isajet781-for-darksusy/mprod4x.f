C-----------------------------------------------------------------------
      SUBROUTINE MPROD4X(X,I,A,B,C,D)
C-----------------------------------------------------------------------
C
C     Multiplies four 3x3 matrices A, B, C and D and stores result in 
C     i-th submatrix of 250x3x3 matrix X.
C
C     Note: Uses unfolded DO-loop and two-stage multiplication
c           that greatly increases speed.
C           54*2=108 floating point operations
C
C     Created: 12/13/07 by Azar Mustafayev
C      
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      REAL*8 X(250,3,3),A(3,3),B(3,3),C(3,3),D(3,3),Y(3,3)
      INTEGER I
      
c...compute axiliary product matrix Y      
      CALL MPROD3(Y,A,B,C)

c...multiply axiliary matrix Y by the forth matrix D      
      CALL MPROD2X(X,I,Y,D)
      
      RETURN
      END      
