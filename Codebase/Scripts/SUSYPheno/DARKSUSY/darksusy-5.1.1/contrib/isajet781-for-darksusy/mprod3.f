C-----------------------------------------------------------------------
      SUBROUTINE MPROD3(Y,A,B,C)
C-----------------------------------------------------------------------
C
C     Multiplies three 3x3 matrices A, B and C and stores result in 
C     3x3 matrix Y.
C
C     Note: Uses unfolded DO-loop and two-stage multiplication
c           that greatly increases speed.
C           54*2=108 floating point operations
C
C     Created: 12/13/07 by Azar Mustafayev
C      
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      REAL*8 Y(3,3),A(3,3),B(3,3),C(3,3),Z(3,3)

c...compute axiliary product matrix Z      
      Z(1,1)= A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
      Z(1,2)= A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
      Z(1,3)= A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)
      
      Z(2,1)= A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1)
      Z(2,2)= A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)
      Z(2,3)= A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
      
      Z(3,1)= A(3,1)*B(1,1)+A(3,2)*B(2,1)+A(3,3)*B(3,1)
      Z(3,2)= A(3,1)*B(1,2)+A(3,2)*B(2,2)+A(3,3)*B(3,2)
      Z(3,3)= A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3)

c...multiply axiliary matrix Z by the third matrix C      
      Y(1,1)= Z(1,1)*C(1,1)+Z(1,2)*C(2,1)+Z(1,3)*C(3,1)
      Y(1,2)= Z(1,1)*C(1,2)+Z(1,2)*C(2,2)+Z(1,3)*C(3,2)
      Y(1,3)= Z(1,1)*C(1,3)+Z(1,2)*C(2,3)+Z(1,3)*C(3,3)
      
      Y(2,1)= Z(2,1)*C(1,1)+Z(2,2)*C(2,1)+Z(2,3)*C(3,1)
      Y(2,2)= Z(2,1)*C(1,2)+Z(2,2)*C(2,2)+Z(2,3)*C(3,2)
      Y(2,3)= Z(2,1)*C(1,3)+Z(2,2)*C(2,3)+Z(2,3)*C(3,3)
      
      Y(3,1)= Z(3,1)*C(1,1)+Z(3,2)*C(2,1)+Z(3,3)*C(3,1)
      Y(3,2)= Z(3,1)*C(1,2)+Z(3,2)*C(2,2)+Z(3,3)*C(3,2)
      Y(3,3)= Z(3,1)*C(1,3)+Z(3,2)*C(2,3)+Z(3,3)*C(3,3)
      
      RETURN
      END      
