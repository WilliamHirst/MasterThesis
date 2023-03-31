CDECK  ID>, RPINT1
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
       FUNCTION RPINT1(RHO)
C-----------------------------------------------------------------------
C      INTEGRAND FOR THE AMPLITUDE SQUARE PIECE OF THE THREE BODY DECAY
C-----------------------------------------------------------------------
       DOUBLE PRECISION RHO,X,RPINT1,M23MIN,M23MAX,E2STAR,E3STAR,
     &                  Y,LIMIT(2),E2MIMA,E3MIMA
       INTEGER          I
C--common block to pass the masses,etc
       COMMON /SQUTRM/ M(4),GAM,MR,A,B
       DOUBLE PRECISION M,GAM,MR,A,B,RPRTCH
       EXTERNAL RPRTCH
C--Calculate the change of variables
       X = MR**2+GAM*MR*TAN(GAM*MR*RHO)
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
          LIMIT(I) =         Y*(X - M(1)**2 - M(2)**2)*
     -  (4*A*B*M(3)*M(4) + (A**2 + B**2)*(-X + M(3)**2 + M(4)**2))
          Y=M23MAX
       ENDDO
       RPINT1 = LIMIT(2)-LIMIT(1)
       END
