C-----------------------------------------------------------------------
      SUBROUTINE MPROD2(Y,A,B)
C-----------------------------------------------------------------------
C
C     Multiplies two 3x3 matrices A and B and stores result in 
C     3x3 matrix Y.
C
C     Note: Uses unfolded DO-loop that greatly increases speed.
C           6*9=54 floating point operations
C
C     Created: 12/13/07 by Azar Mustafayev
C      
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      REAL*8 Y(3,3),A(3,3),B(3,3)
      
      Y(1,1)= A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
      Y(1,2)= A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
      Y(1,3)= A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)
      
      Y(2,1)= A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1)
      Y(2,2)= A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)
      Y(2,3)= A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
      
      Y(3,1)= A(3,1)*B(1,1)+A(3,2)*B(2,1)+A(3,3)*B(3,1)
      Y(3,2)= A(3,1)*B(1,2)+A(3,2)*B(2,2)+A(3,3)*B(3,2)
      Y(3,3)= A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3)
      
      RETURN
      END      
