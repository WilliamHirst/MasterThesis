CDECK  ID>, RPINT3
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
       FUNCTION RPINT3(RHO)
C-----------------------------------------------------------------------
C      FUNCTION FOR INTEGRAND OF THE NORMAL INTERFERENCE TERM IN A 
C      3-BODY R-PARITY VIOLATING DECAY 
C-----------------------------------------------------------------------
       IMPLICIT NONE
       DOUBLE PRECISION RHO,X,RPINT3,M23MIN,M23MAX,E2STAR,E3STAR,
     &                  Y,LIMIT(2),E2MIMA,E3MIMA
       INTEGER          I
C--common block to pass the masses,etc
       COMMON /INFTRM/ M(4),MSQ(4),GAM(2),MR(2),MRSQ(2),A(2),B(2)
       DOUBLE PRECISION M,MSQ,GAM,MR,MRSQ,A,B,RPRTCH
       EXTERNAL RPRTCH
C--Calculate the change of variables
       X = MR(1)**2+GAM(1)*ABS(MR(1))*TAN(GAM(1)*ABS(MR(1))*RHO)
C--Evaulate limits on the inner integral
       E2STAR = (X-M(1)**2+M(2)**2)/(2*SQRT(X))
       E3STAR = (M(4)**2-X-M(3)**2)/(2*SQRT(X))
       E2MIMA = RPRTCH(M(2),E2STAR)
       E3MIMA = RPRTCH(M(3),E3STAR)
       M23MAX = (E2STAR+E3STAR)**2-(E2MIMA-E3MIMA)**2
       M23MIN = (E2STAR+E3STAR)**2-(E2MIMA+E3MIMA)**2
C--Do the inner integral
       Y = M23MIN
       DO I=1,2
          LIMIT(I) =Y*(X*B(1)*B(2)+A(1)*A(2)*M(1)*M(3)
     -      +A(1)*B(2)*M(1)*M(4))*
     -   (X-MRSQ(1))+(ATAN((Y-MRSQ(2))/(GAM(2)*SQRT(MRSQ(2))))*
     -     (X*A(1)*A(2)*GAM(1)*M(1)*M(3)*MR(1)*MR(2) + 
     -       X*A(2)*B(1)*GAM(1)*M(3)*M(4)*MR(1)*MR(2) - 
     -       X**2*B(1)*B(2)*GAM(2)*MRSQ(2) - 
     -       X*A(1)*A(2)*GAM(2)*M(1)*M(3)*MRSQ(2) - 
     -       X*A(1)*B(2)*GAM(2)*M(1)*M(4)*MRSQ(2) + 
     -       X*B(1)*B(2)*GAM(1)*MR(1)*MR(2)*MRSQ(2) + 
     -       A(1)*A(2)*GAM(1)*M(1)*M(3)*MR(1)*MR(2)*MRSQ(2) + 
     -       A(1)*B(2)*GAM(1)*M(1)*M(4)*MR(1)*MR(2)*MRSQ(2) + 
     -       X*B(1)*B(2)*GAM(2)*MRSQ(1)*MRSQ(2) + 
     -       A(1)*A(2)*GAM(2)*M(1)*M(3)*MRSQ(1)*MRSQ(2) + 
     -       A(1)*B(2)*GAM(2)*M(1)*M(4)*MRSQ(1)*MRSQ(2) - 
     -       A(1)*A(2)*GAM(1)*M(1)*M(3)*MR(1)*MR(2)*MSQ(1) - 
     -       A(2)*B(1)*GAM(1)*M(3)*M(4)*MR(1)*MR(2)*MSQ(1) - 
     -       A(1)*B(2)*GAM(1)*M(1)*M(4)*MR(1)*MR(2)*MSQ(2) - 
     -       A(2)*B(1)*GAM(1)*M(3)*M(4)*MR(1)*MR(2)*MSQ(2) - 
     -       A(1)*A(2)*GAM(1)*M(1)*M(3)*MR(1)*MR(2)*MSQ(3) - 
     -       A(1)*B(2)*GAM(1)*M(1)*M(4)*MR(1)*MR(2)*MSQ(3) - 
     -       B(1)*B(2)*GAM(1)*MR(1)*MR(2)*MSQ(1)*MSQ(3) - 
     - B(1)*B(2)*GAM(1)*MR(1)*MR(2)*MSQ(2)*MSQ(4)))/SQRT(MRSQ(2))+ 
     -  (LOG(Y**2 - 2*Y*MRSQ(2) + GAM(2)**2*MRSQ(2) + MRSQ(2)**2)*
     -     (X**2*A(1)*A(2)*M(1)*M(3) + X**2*A(2)*B(1)*M(3)*M(4) + 
     -       X*B(1)*B(2)*GAM(1)*GAM(2)*MR(1)*MR(2) + 
     -       A(1)*A(2)*GAM(1)*GAM(2)*M(1)*M(3)*MR(1)*MR(2) + 
     -       A(1)*B(2)*GAM(1)*GAM(2)*M(1)*M(4)*MR(1)*MR(2) - 
     - X*A(1)*A(2)*M(1)*M(3)*MRSQ(1)-X*A(2)*B(1)*M(3)*M(4)*MRSQ(1)+ 
     -       X**2*B(1)*B(2)*MRSQ(2) + X*A(1)*A(2)*M(1)*M(3)*MRSQ(2)+ 
     -    X*A(1)*B(2)*M(1)*M(4)*MRSQ(2)-X*B(1)*B(2)*MRSQ(1)*MRSQ(2)- 
     -       A(1)*A(2)*M(1)*M(3)*MRSQ(1)*MRSQ(2) - 
     -       A(1)*B(2)*M(1)*M(4)*MRSQ(1)*MRSQ(2) - 
     - X*A(1)*A(2)*M(1)*M(3)*MSQ(1) - X*A(2)*B(1)*M(3)*M(4)*MSQ(1) + 
     -       A(1)*A(2)*M(1)*M(3)*MRSQ(1)*MSQ(1) + 
     -       A(2)*B(1)*M(3)*M(4)*MRSQ(1)*MSQ(1) - 
     - X*A(1)*B(2)*M(1)*M(4)*MSQ(2) - X*A(2)*B(1)*M(3)*M(4)*MSQ(2) + 
     -       A(1)*B(2)*M(1)*M(4)*MRSQ(1)*MSQ(2) + 
     -       A(2)*B(1)*M(3)*M(4)*MRSQ(1)*MSQ(2) - 
     - X*A(1)*A(2)*M(1)*M(3)*MSQ(3) - X*A(1)*B(2)*M(1)*M(4)*MSQ(3) + 
     -     A(1)*A(2)*M(1)*M(3)*MRSQ(1)*MSQ(3) + 
     - A(1)*B(2)*M(1)*M(4)*MRSQ(1)*MSQ(3)-X*B(1)*B(2)*MSQ(1)*MSQ(3) + 
     -  B(1)*B(2)*MRSQ(1)*MSQ(1)*MSQ(3) - X*B(1)*B(2)*MSQ(2)*MSQ(4) + 
     - B(1)*B(2)*MRSQ(1)*MSQ(2)*MSQ(4)))/2.0D0
          Y=M23MAX
       ENDDO
       RPINT3 = LIMIT(2)-LIMIT(1)
       END
