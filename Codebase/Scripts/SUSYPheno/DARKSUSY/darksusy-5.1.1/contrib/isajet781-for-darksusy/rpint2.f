CDECK  ID>, RPINT2
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
       FUNCTION RPINT2(RHO)
C-----------------------------------------------------------------------
C      FUNCTION FOR THE INTEGRAND FOR THE LIGHT/HEAVY INTERFERENCE
C      TERM IN A 3-BODY R-PARITY VIOLATING DECAY
C-----------------------------------------------------------------------
       DOUBLE PRECISION RHO,X,RPINT2,M23MIN,M23MAX,E2STAR,E3STAR,
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
C--Do the inner intergral
       Y = M23MIN
       DO I=1,2
          LIMIT(I) =         (Y*(X - MSQ(1) - MSQ(2))*
     -    (2*(A(2)*B(1) + A(1)*B(2))*M(3)*M(4) + 
     -      (A(1)*A(2) + B(1)*B(2))*(-X + MSQ(3) + MSQ(4)))*
     -    (GAM(1)*GAM(2)*MR(1)*MR(2) + (X - MR(1)**2)*(X - MR(2)**2)))
     -   /(GAM(2)**2*MR(2)**2 + (X - MR(2)**2)**2)
          Y=M23MAX
       ENDDO
       RPINT2 = LIMIT(2)-LIMIT(1)
       END
